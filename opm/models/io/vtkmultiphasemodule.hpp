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
 * \copydoc Opm::VtkMultiPhaseModule
 */
#ifndef OPM_VTK_MULTI_PHASE_MODULE_HPP
#define OPM_VTK_MULTI_PHASE_MODULE_HPP

#include <dune/common/fvector.hh>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkmultiphaseparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for quantities which make sense for all
 *        models which deal with multiple fluid phases in porous media
 *        that don't use flashy concepts like interfacial area.
 *
 * This module deals with the following quantities:
 * - Pressures of all fluid phases
 * - Densities of all fluid phases
 * - Saturations of all fluid phases
 * - Mobilities of all fluid phases
 * - Relative permeabilities of all fluid phases
 * - Viscosities of all fluid phases
 * - Average molar masses of all fluid phases
 * - Porosity of the medium
 * - Norm of the intrinsic permeability of the medium
 */
template<class TypeTag>
class VtkMultiPhaseModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using DiscBaseOutputModule = GetPropType<TypeTag, Properties::DiscBaseOutputModule>;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };

    using BufferType = typename ParentType::BufferType;
    using ScalarBuffer = typename ParentType::ScalarBuffer;
    using VectorBuffer = typename ParentType::VectorBuffer;
    using TensorBuffer = typename ParentType::TensorBuffer;
    using PhaseBuffer = typename ParentType::PhaseBuffer;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;

    using PhaseVectorBuffer = std::array<VectorBuffer, numPhases>;

public:
    explicit VtkMultiPhaseModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        params_.read();
    }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output module.
     */
    static void registerParameters()
    {
        VtkMultiPhaseParams::registerParameters();
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if (params_.extrusionFactorOutput_) {
            this->resizeScalarBuffer_(extrusionFactor_, BufferType::Dof);
        }
        if (params_.pressureOutput_) {
            this->resizePhaseBuffer_(pressure_, BufferType::Dof);
        }
        if (params_.densityOutput_) {
            this->resizePhaseBuffer_(density_, BufferType::Dof);
        }
        if (params_.saturationOutput_) {
            this->resizePhaseBuffer_(saturation_, BufferType::Dof);
        }
        if (params_.mobilityOutput_) {
            this->resizePhaseBuffer_(mobility_, BufferType::Dof);
        }
        if (params_.relativePermeabilityOutput_) {
            this->resizePhaseBuffer_(relativePermeability_, BufferType::Dof);
        }
        if (params_.viscosityOutput_) {
            this->resizePhaseBuffer_(viscosity_, BufferType::Dof);
        }
        if (params_.averageMolarMassOutput_) {
            this->resizePhaseBuffer_(averageMolarMass_, BufferType::Dof);
        }

        if (params_.porosityOutput_) {
            this->resizeScalarBuffer_(porosity_, BufferType::Dof);
        }
        if (params_.intrinsicPermeabilityOutput_) {
            this->resizeTensorBuffer_(intrinsicPermeability_, BufferType::Dof);
        }

        if (params_.velocityOutput_) {
            const std::size_t nDof = this->simulator_.model().numGridDof();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                velocity_[phaseIdx].resize(nDof);
                for (unsigned dofIdx = 0; dofIdx < nDof; ++ dofIdx) {
                    velocity_[phaseIdx][dofIdx].resize(dimWorld);
                    velocity_[phaseIdx][dofIdx] = 0.0;
                }
            }
            this->resizePhaseBuffer_(velocityWeight_, BufferType::Dof);
        }

        if (params_.potentialGradientOutput_) {
            const std::size_t nDof = this->simulator_.model().numGridDof();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                potentialGradient_[phaseIdx].resize(nDof);
                for (unsigned dofIdx = 0; dofIdx < nDof; ++ dofIdx) {
                    potentialGradient_[phaseIdx][dofIdx].resize(dimWorld);
                    potentialGradient_[phaseIdx][dofIdx] = 0.0;
                }
            }

            this->resizePhaseBuffer_(potentialWeight_, BufferType::Dof);
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities seen on
     *        an element
     */
    void processElement(const ElementContext& elemCtx) override
    {
        if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
            return;
        }

        const auto& problem = elemCtx.problem();
        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            const unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            if (params_.extrusionFactorOutput_) {
                extrusionFactor_[I] = intQuants.extrusionFactor();
            }
            if (params_.porosityOutput_) {
                porosity_[I] = getValue(intQuants.porosity());
            }

            if (params_.intrinsicPermeabilityOutput_) {
                const auto& K = problem.intrinsicPermeability(elemCtx, i, /*timeIdx=*/0);
                for (unsigned rowIdx = 0; rowIdx < K.rows; ++rowIdx) {
                    for (unsigned colIdx = 0; colIdx < K.cols; ++colIdx) {
                        intrinsicPermeability_[I][rowIdx][colIdx] = K[rowIdx][colIdx];
                    }
                }
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }
                if (params_.pressureOutput_) {
                    pressure_[phaseIdx][I] = getValue(fs.pressure(phaseIdx));
                }
                if (params_.densityOutput_) {
                    density_[phaseIdx][I] = getValue(fs.density(phaseIdx));
                }
                if (params_.saturationOutput_) {
                    saturation_[phaseIdx][I] = getValue(fs.saturation(phaseIdx));
                }
                if (params_.mobilityOutput_) {
                    mobility_[phaseIdx][I] = getValue(intQuants.mobility(phaseIdx));
                }
                if (params_.relativePermeabilityOutput_) {
                    relativePermeability_[phaseIdx][I] = getValue(intQuants.relativePermeability(phaseIdx));
                }
                if (params_.viscosityOutput_) {
                    viscosity_[phaseIdx][I] = getValue(fs.viscosity(phaseIdx));
                }
                if (params_.averageMolarMassOutput_) {
                    averageMolarMass_[phaseIdx][I] = getValue(fs.averageMolarMass(phaseIdx));
                }
            }
        }

        if (params_.potentialGradientOutput_) {
            // calculate velocities if requested
            for (unsigned faceIdx = 0; faceIdx < elemCtx.numInteriorFaces(/*timeIdx=*/0); ++faceIdx) {
                const auto& extQuants = elemCtx.extensiveQuantities(faceIdx, /*timeIdx=*/0);

                const unsigned i = extQuants.interiorIndex();
                const unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);

                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    const Scalar weight = extQuants.extrusionFactor();

                    potentialWeight_[phaseIdx][I] += weight;

                    const auto& inputPGrad = extQuants.potentialGrad(phaseIdx);
                    DimVector pGrad;
                    for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx) {
                        pGrad[dimIdx] = getValue(inputPGrad[dimIdx]) * weight;
                    }
                    potentialGradient_[phaseIdx][I] += pGrad;
                } // end for all phases
            } // end for all faces
        }

        if (params_.velocityOutput_) {
            // calculate velocities if requested
            for (unsigned faceIdx = 0; faceIdx < elemCtx.numInteriorFaces(/*timeIdx=*/0); ++faceIdx) {
                const auto& extQuants = elemCtx.extensiveQuantities(faceIdx, /*timeIdx=*/0);

                const unsigned i = extQuants.interiorIndex();
                const unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);

                const unsigned j = extQuants.exteriorIndex();
                const unsigned J = elemCtx.globalSpaceIndex(j, /*timeIdx=*/0);

                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    Scalar weight = std::max(Scalar{1e-16},
                                             std::abs(getValue(extQuants.volumeFlux(phaseIdx))));
                    Valgrind::CheckDefined(extQuants.extrusionFactor());
                    assert(extQuants.extrusionFactor() > 0);
                    weight *= extQuants.extrusionFactor();

                    const auto& inputV = extQuants.filterVelocity(phaseIdx);
                    DimVector v;
                    for (unsigned k = 0; k < dimWorld; ++k) {
                        v[k] = getValue(inputV[k]);
                    }
                    if (v.two_norm() > 1e-20) {
                        weight /= v.two_norm();
                    }
                    v *= weight;

                    velocity_[phaseIdx][I] += v;
                    velocity_[phaseIdx][J] += v;

                    velocityWeight_[phaseIdx][I] += weight;
                    velocityWeight_[phaseIdx][J] += weight;
                } // end for all phases
            } // end for all faces
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

        if (params_.extrusionFactorOutput_) {
            this->commitScalarBuffer_(baseWriter, "extrusionFactor",
                                      extrusionFactor_, BufferType::Dof);
        }
        if (params_.pressureOutput_) {
            this->commitPhaseBuffer_(baseWriter, "pressure_%s", pressure_, BufferType::Dof);
        }
        if (params_.densityOutput_) {
            this->commitPhaseBuffer_(baseWriter, "density_%s", density_, BufferType::Dof);
        }
        if (params_.saturationOutput_) {
            this->commitPhaseBuffer_(baseWriter, "saturation_%s", saturation_, BufferType::Dof);
        }
        if (params_.mobilityOutput_) {
            this->commitPhaseBuffer_(baseWriter, "mobility_%s", mobility_, BufferType::Dof);
        }
        if (params_.relativePermeabilityOutput_) {
            this->commitPhaseBuffer_(baseWriter, "relativePerm_%s",
                                     relativePermeability_, BufferType::Dof);
        }
        if (params_.viscosityOutput_) {
            this->commitPhaseBuffer_(baseWriter, "viscosity_%s", viscosity_, BufferType::Dof);
        }
        if (params_.averageMolarMassOutput_) {
            this->commitPhaseBuffer_(baseWriter, "averageMolarMass_%s",
                                     averageMolarMass_, BufferType::Dof);
        }

        if (params_.porosityOutput_) {
            this->commitScalarBuffer_(baseWriter, "porosity", porosity_, BufferType::Dof);
        }
        if (params_.intrinsicPermeabilityOutput_) {
            this->commitTensorBuffer_(baseWriter, "intrinsicPerm",
                                      intrinsicPermeability_, BufferType::Dof);
        }

        if (params_.velocityOutput_) {
            const std::size_t numDof = this->simulator_.model().numGridDof();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // first, divide the velocity field by the
                // respective finite volume's surface area
                for (unsigned i = 0; i < numDof; ++i) {
                    velocity_[phaseIdx][i] /= velocityWeight_[phaseIdx][i];
                }
                // commit the phase velocity
                char name[512];
                snprintf(name, 512, "filterVelocity_%s", FluidSystem::phaseName(phaseIdx).data());

                DiscBaseOutputModule::attachVectorDofData_(baseWriter, velocity_[phaseIdx], name);
            }
        }

        if (params_.potentialGradientOutput_) {
            const std::size_t numDof = this->simulator_.model().numGridDof();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // first, divide the velocity field by the
                // respective finite volume's surface area
                for (unsigned i = 0; i < numDof; ++i) {
                    potentialGradient_[phaseIdx][i] /= potentialWeight_[phaseIdx][i];
                }
                // commit the phase velocity
                char name[512];
                snprintf(name, 512, "gradP_%s", FluidSystem::phaseName(phaseIdx).data());

                DiscBaseOutputModule::attachVectorDofData_(baseWriter,
                                                           potentialGradient_[phaseIdx],
                                                           name);
            }
        }
    }

    /*!
     * \brief Returns true iff the module needs to access the extensive quantities of a
     * context to do its job.
     *
     * For example, this happens if velocities or gradients should be written. Always
     * returning true here does not do any harm from the correctness perspective, but it
     * slows down writing the output fields.
     */
    bool needExtensiveQuantities() const override
    {
        return params_.velocityOutput_ || params_.potentialGradientOutput_;
    }

private:
    VtkMultiPhaseParams params_{};
    ScalarBuffer extrusionFactor_{};
    PhaseBuffer pressure_{};
    PhaseBuffer density_{};
    PhaseBuffer saturation_{};
    PhaseBuffer mobility_{};
    PhaseBuffer relativePermeability_{};
    PhaseBuffer viscosity_{};
    PhaseBuffer averageMolarMass_{};

    ScalarBuffer porosity_{};
    TensorBuffer intrinsicPermeability_{};

    PhaseVectorBuffer velocity_{};
    PhaseBuffer velocityWeight_{};

    PhaseVectorBuffer potentialGradient_{};
    PhaseBuffer potentialWeight_{};
};

} // namespace Opm

#endif // OPM_VTK_MULTI_PHASE_MODULE_HPP
