// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::VtkDiscreteFractureModule
 */
#ifndef OPM_VTK_DISCRETE_FRACTURE_MODULE_HPP
#define OPM_VTK_DISCRETE_FRACTURE_MODULE_HPP

#include <dune/common/fvector.hh>

#include <opm/material/common/Valgrind.hpp>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkdiscretefractureparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstddef>

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for quantities which make sense for all
 *        models which deal with discrete fractures in porous media.
 *
 * This module deals with the following quantities:
 * - Saturations of all fluid phases in the fracture
 * - Mobilities of all fluid phases in the fracture
 * - Relative permeabilities of all fluid phases in the fracture
 * - Porosity of the medium in the fracture
 * - Norm of the intrinsic permeability of the medium in the fracture
 */
template <class TypeTag>
class VtkDiscreteFractureModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using DiscBaseOutputModule = GetPropType<TypeTag, Properties::DiscBaseOutputModule>;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };

    using BufferType = typename ParentType::BufferType;
    using PhaseBuffer = typename ParentType::PhaseBuffer;
    using PhaseVectorBuffer = typename ParentType::PhaseVectorBuffer;
    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    explicit VtkDiscreteFractureModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        params_.read();
    }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output module.
     */
    static void registerParameters()
    {
        VtkDiscreteFractureParams::registerParameters();
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if (params_.saturationOutput_) {
            this->resizePhaseBuffer_(fractureSaturation_, BufferType::Dof);
        }
        if (params_.mobilityOutput_) {
            this->resizePhaseBuffer_(fractureMobility_, BufferType::Dof);
        }
        if (params_.relativePermeabilityOutput_) {
            this->resizePhaseBuffer_(fractureRelativePermeability_, BufferType::Dof);
        }

        if (params_.porosityOutput_) {
            this->resizeScalarBuffer_(fracturePorosity_, BufferType::Dof);
        }
        if (params_.intrinsicPermeabilityOutput_) {
            this->resizeScalarBuffer_(fractureIntrinsicPermeability_, BufferType::Dof);
        }
        if (params_.volumeFractionOutput_) {
            this->resizeScalarBuffer_(fractureVolumeFraction_, BufferType::Dof);
        }

        if (params_.velocityOutput_) {
            const std::size_t nDof = this->simulator_.model().numGridDof();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                fractureVelocity_[phaseIdx].resize(nDof);
                for (unsigned dofIdx = 0; dofIdx < nDof; ++dofIdx) {
                    fractureVelocity_[phaseIdx][dofIdx].resize(dimWorld);
                    fractureVelocity_[phaseIdx][dofIdx] = 0.0;
                }
            }
            this->resizePhaseBuffer_(fractureVelocityWeight_, BufferType::Dof);
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx) override
    {
        if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
            return;
        }

        const auto& fractureMapper = elemCtx.simulator().vanguard().fractureMapper();

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            const unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            if (!fractureMapper.isFractureVertex(I)) {
                continue;
            }

            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fractureFluidState();

            if (params_.porosityOutput_) {
                Valgrind::CheckDefined(intQuants.fracturePorosity());
                fracturePorosity_[I] = intQuants.fracturePorosity();
            }
            if (params_.intrinsicPermeabilityOutput_) {
                const auto& K = intQuants.fractureIntrinsicPermeability();
                fractureIntrinsicPermeability_[I] = K[0][0];
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (params_.saturationOutput_) {
                    Valgrind::CheckDefined(fs.saturation(phaseIdx));
                    fractureSaturation_[phaseIdx][I] = fs.saturation(phaseIdx);
                }
                if (params_.mobilityOutput_) {
                    Valgrind::CheckDefined(intQuants.fractureMobility(phaseIdx));
                    fractureMobility_[phaseIdx][I] = intQuants.fractureMobility(phaseIdx);
                }
                if (params_.relativePermeabilityOutput_) {
                    Valgrind::CheckDefined(intQuants.fractureRelativePermeability(phaseIdx));
                    fractureRelativePermeability_[phaseIdx][I] =
                        intQuants.fractureRelativePermeability(phaseIdx);
                }
                if (params_.volumeFractionOutput_) {
                    Valgrind::CheckDefined(intQuants.fractureVolume());
                    fractureVolumeFraction_[I] += intQuants.fractureVolume();
                }
            }
        }

        if (params_.velocityOutput_) {
            // calculate velocities if requested by the simulator
            for (unsigned scvfIdx = 0; scvfIdx < elemCtx.numInteriorFaces(/*timeIdx=*/0); ++scvfIdx) {
                const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, /*timeIdx=*/0);

                const unsigned i = extQuants.interiorIndex();
                const unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);

                const unsigned j = extQuants.exteriorIndex();
                const unsigned J = elemCtx.globalSpaceIndex(j, /*timeIdx=*/0);

                if (!fractureMapper.isFractureEdge(I, J)) {
                    continue;
                }

                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    Scalar weight = std::max(Scalar{1e-16},
                                             std::abs(extQuants.fractureVolumeFlux(phaseIdx)));
                    Valgrind::CheckDefined(extQuants.extrusionFactor());
                    assert(extQuants.extrusionFactor() > 0);
                    weight *= extQuants.extrusionFactor();

                    Dune::FieldVector<Scalar, dim> v(extQuants.fractureFilterVelocity(phaseIdx));
                    v *= weight;

                    for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
                        fractureVelocity_[phaseIdx][I][dimIdx] += v[dimIdx];
                        fractureVelocity_[phaseIdx][J][dimIdx] += v[dimIdx];
                    }

                    fractureVelocityWeight_[phaseIdx][I] += weight;
                    fractureVelocityWeight_[phaseIdx][J] += weight;
                }
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter) override
    {
        if (!dynamic_cast<VtkMultiWriter*>(&baseWriter)) {
            return;
        }

        if (params_.saturationOutput_) {
            this->commitPhaseBuffer_(baseWriter, "fractureSaturation_%s",
                                     fractureSaturation_, BufferType::Dof);
        }
        if (params_.mobilityOutput_) {
            this->commitPhaseBuffer_(baseWriter, "fractureMobility_%s",
                                     fractureMobility_, BufferType::Dof);
        }
        if (params_.relativePermeabilityOutput_) {
            this->commitPhaseBuffer_(baseWriter, "fractureRelativePerm_%s",
                                     fractureRelativePermeability_, BufferType::Dof);
        }

        if (params_.porosityOutput_) {
            this->commitScalarBuffer_(baseWriter, "fracturePorosity",
                                      fracturePorosity_, BufferType::Dof);
        }
        if (params_.intrinsicPermeabilityOutput_) {
            this->commitScalarBuffer_(baseWriter, "fractureIntrinsicPerm",
                                      fractureIntrinsicPermeability_, BufferType::Dof);
        }
        if (params_.volumeFractionOutput_) {
            // divide the fracture volume by the total volume of the finite volumes
            for (unsigned I = 0; I < fractureVolumeFraction_.size(); ++I) {
                fractureVolumeFraction_[I] /= this->simulator_.model().dofTotalVolume(I);
            }
            this->commitScalarBuffer_(baseWriter, "fractureVolumeFraction",
                                      fractureVolumeFraction_, BufferType::Dof);
        }

        if (params_.velocityOutput_) {
            const std::size_t nDof = this->simulator_.model().numGridDof();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // first, divide the velocity field by the
                // respective finite volume's surface area
                for (unsigned dofIdx = 0; dofIdx < nDof; ++dofIdx) {
                    fractureVelocity_[phaseIdx][dofIdx] /=
                        std::max<Scalar>(1e-20, fractureVelocityWeight_[phaseIdx][dofIdx]);
                }
                // commit the phase velocity
                char name[512];
                snprintf(name, 512, "fractureFilterVelocity_%s", FluidSystem::phaseName(phaseIdx).data());

                DiscBaseOutputModule::attachVectorDofData_(baseWriter, fractureVelocity_[phaseIdx], name);
            }
        }
    }

private:
    VtkDiscreteFractureParams params_{};
    PhaseBuffer fractureSaturation_{};
    PhaseBuffer fractureMobility_{};
    PhaseBuffer fractureRelativePermeability_{};

    ScalarBuffer fracturePorosity_{};
    ScalarBuffer fractureVolumeFraction_{};
    ScalarBuffer fractureIntrinsicPermeability_{};

    PhaseVectorBuffer fractureVelocity_{};
    PhaseBuffer fractureVelocityWeight_{};

    PhaseVectorBuffer potentialGradient_{};
    PhaseBuffer potentialWeight_{};
};

} // namespace Opm

#endif // OPM_VTK_DISCRETE_FRACTURE_MODULE_HPP
