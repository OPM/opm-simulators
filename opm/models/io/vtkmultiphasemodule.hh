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
#ifndef EWOMS_VTK_MULTI_PHASE_MODULE_HH
#define EWOMS_VTK_MULTI_PHASE_MODULE_HH

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>

#include <cstdio>

namespace Opm::Properties {

namespace TTag {

// create new type tag for the VTK multi-phase output
struct VtkMultiPhase {};

} // namespace TTag

// create the property tags needed for the multi phase module
template<class TypeTag, class MyTypeTag>
struct VtkWriteExtrusionFactor { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWritePressures { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteDensities { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteSaturations { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteMobilities { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteRelativePermeabilities { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteViscosities { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteAverageMolarMasses { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWritePorosity { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteIntrinsicPermeabilities { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWritePotentialGradients { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteFilterVelocities { using type = UndefinedProperty; };

// set default values for what quantities to output
template<class TypeTag>
struct VtkWriteExtrusionFactor<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWritePressures<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteDensities<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteSaturations<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteMobilities<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWriteRelativePermeabilities<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteViscosities<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWriteAverageMolarMasses<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWritePorosity<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteIntrinsicPermeabilities<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWritePotentialGradients<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWriteFilterVelocities<TypeTag, TTag::VtkMultiPhase> { static constexpr bool value = false; };

} // namespace Opm::Properties

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

    static const int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };

    using ScalarBuffer = typename ParentType::ScalarBuffer;
    using VectorBuffer = typename ParentType::VectorBuffer;
    using TensorBuffer = typename ParentType::TensorBuffer;
    using PhaseBuffer = typename ParentType::PhaseBuffer;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;

    using PhaseVectorBuffer = std::array<VectorBuffer, numPhases>;

public:
    VtkMultiPhaseModule(const Simulator& simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteExtrusionFactor,
                             "Include the extrusion factor of the degrees of freedom into the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePressures,
                             "Include the phase pressures in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteDensities,
                             "Include the phase densities in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSaturations,
                             "Include the phase saturations in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteMobilities,
                             "Include the phase mobilities in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteRelativePermeabilities,
                             "Include the phase relative permeabilities in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteViscosities,
                             "Include component phase viscosities in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteAverageMolarMasses,
                             "Include the average phase mass in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePorosity,
                             "Include the porosity in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteIntrinsicPermeabilities,
                             "Include the intrinsic permeability in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFilterVelocities,
                             "Include in the filter velocities of the phases the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePotentialGradients,
                             "Include the phase pressure potential gradients in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (extrusionFactorOutput_()) this->resizeScalarBuffer_(extrusionFactor_);
        if (pressureOutput_()) this->resizePhaseBuffer_(pressure_);
        if (densityOutput_()) this->resizePhaseBuffer_(density_);
        if (saturationOutput_()) this->resizePhaseBuffer_(saturation_);
        if (mobilityOutput_()) this->resizePhaseBuffer_(mobility_);
        if (relativePermeabilityOutput_()) this->resizePhaseBuffer_(relativePermeability_);
        if (viscosityOutput_()) this->resizePhaseBuffer_(viscosity_);
        if (averageMolarMassOutput_()) this->resizePhaseBuffer_(averageMolarMass_);

        if (porosityOutput_()) this->resizeScalarBuffer_(porosity_);
        if (intrinsicPermeabilityOutput_()) this->resizeTensorBuffer_(intrinsicPermeability_);

        if (velocityOutput_()) {
            size_t nDof = this->simulator_.model().numGridDof();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                velocity_[phaseIdx].resize(nDof);
                for (unsigned dofIdx = 0; dofIdx < nDof; ++ dofIdx) {
                    velocity_[phaseIdx][dofIdx].resize(dimWorld);
                    velocity_[phaseIdx][dofIdx] = 0.0;
                }
            }
            this->resizePhaseBuffer_(velocityWeight_);
        }

        if (potentialGradientOutput_()) {
            size_t nDof = this->simulator_.model().numGridDof();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                potentialGradient_[phaseIdx].resize(nDof);
                for (unsigned dofIdx = 0; dofIdx < nDof; ++ dofIdx) {
                    potentialGradient_[phaseIdx][dofIdx].resize(dimWorld);
                    potentialGradient_[phaseIdx][dofIdx] = 0.0;
                }
            }

            this->resizePhaseBuffer_(potentialWeight_);
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities seen on
     *        an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        const auto& problem = elemCtx.problem();
        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            if (extrusionFactorOutput_()) extrusionFactor_[I] = intQuants.extrusionFactor();
            if (porosityOutput_()) porosity_[I] = getValue(intQuants.porosity());

            if (intrinsicPermeabilityOutput_()) {
                const auto& K = problem.intrinsicPermeability(elemCtx, i, /*timeIdx=*/0);
                for (unsigned rowIdx = 0; rowIdx < K.rows; ++rowIdx)
                    for (unsigned colIdx = 0; colIdx < K.cols; ++colIdx)
                        intrinsicPermeability_[I][rowIdx][colIdx] = K[rowIdx][colIdx];
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }
                if (pressureOutput_())
                    pressure_[phaseIdx][I] = getValue(fs.pressure(phaseIdx));
                if (densityOutput_())
                    density_[phaseIdx][I] = getValue(fs.density(phaseIdx));
                if (saturationOutput_())
                    saturation_[phaseIdx][I] = getValue(fs.saturation(phaseIdx));
                if (mobilityOutput_())
                    mobility_[phaseIdx][I] = getValue(intQuants.mobility(phaseIdx));
                if (relativePermeabilityOutput_())
                    relativePermeability_[phaseIdx][I] = getValue(intQuants.relativePermeability(phaseIdx));
                if (viscosityOutput_())
                    viscosity_[phaseIdx][I] = getValue(fs.viscosity(phaseIdx));
                if (averageMolarMassOutput_())
                    averageMolarMass_[phaseIdx][I] = getValue(fs.averageMolarMass(phaseIdx));
            }
        }

        if (potentialGradientOutput_()) {
            // calculate velocities if requested
            for (unsigned faceIdx = 0; faceIdx < elemCtx.numInteriorFaces(/*timeIdx=*/0); ++ faceIdx) {
                const auto& extQuants = elemCtx.extensiveQuantities(faceIdx, /*timeIdx=*/0);

                unsigned i = extQuants.interiorIndex();
                unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);

                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    Scalar weight = extQuants.extrusionFactor();

                    potentialWeight_[phaseIdx][I] += weight;

                    const auto& inputPGrad = extQuants.potentialGrad(phaseIdx);
                    DimVector pGrad;
                    for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                        pGrad[dimIdx] = getValue(inputPGrad[dimIdx])*weight;
                    potentialGradient_[phaseIdx][I] += pGrad;
                } // end for all phases
            } // end for all faces
        }

        if (velocityOutput_()) {
            // calculate velocities if requested
            for (unsigned faceIdx = 0; faceIdx < elemCtx.numInteriorFaces(/*timeIdx=*/0); ++ faceIdx) {
                const auto& extQuants = elemCtx.extensiveQuantities(faceIdx, /*timeIdx=*/0);

                unsigned i = extQuants.interiorIndex();
                unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);

                unsigned j = extQuants.exteriorIndex();
                unsigned J = elemCtx.globalSpaceIndex(j, /*timeIdx=*/0);

                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    Scalar weight = std::max<Scalar>(1e-16,
                                                     std::abs(getValue(extQuants.volumeFlux(phaseIdx))));
                    Valgrind::CheckDefined(extQuants.extrusionFactor());
                    assert(extQuants.extrusionFactor() > 0);
                    weight *= extQuants.extrusionFactor();

                    const auto& inputV = extQuants.filterVelocity(phaseIdx);
                    DimVector v;
                    for (unsigned k = 0; k < dimWorld; ++k)
                        v[k] = getValue(inputV[k]);
                    if (v.two_norm() > 1e-20)
                        weight /= v.two_norm();
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
    void commitBuffers(BaseOutputWriter& baseWriter)
    {
        VtkMultiWriter *vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
        if (!vtkWriter)
            return;

        if (extrusionFactorOutput_())
            this->commitScalarBuffer_(baseWriter, "extrusionFactor", extrusionFactor_);
        if (pressureOutput_())
            this->commitPhaseBuffer_(baseWriter, "pressure_%s", pressure_);
        if (densityOutput_())
            this->commitPhaseBuffer_(baseWriter, "density_%s", density_);
        if (saturationOutput_())
            this->commitPhaseBuffer_(baseWriter, "saturation_%s", saturation_);
        if (mobilityOutput_())
            this->commitPhaseBuffer_(baseWriter, "mobility_%s", mobility_);
        if (relativePermeabilityOutput_())
            this->commitPhaseBuffer_(baseWriter, "relativePerm_%s", relativePermeability_);
        if (viscosityOutput_())
            this->commitPhaseBuffer_(baseWriter, "viscosity_%s", viscosity_);
        if (averageMolarMassOutput_())
            this->commitPhaseBuffer_(baseWriter, "averageMolarMass_%s", averageMolarMass_);

        if (porosityOutput_())
            this->commitScalarBuffer_(baseWriter, "porosity", porosity_);
        if (intrinsicPermeabilityOutput_())
            this->commitTensorBuffer_(baseWriter, "intrinsicPerm", intrinsicPermeability_);

        if (velocityOutput_()) {
            size_t numDof = this->simulator_.model().numGridDof();

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // first, divide the velocity field by the
                // respective finite volume's surface area
                for (unsigned i = 0; i < numDof; ++i)
                    velocity_[phaseIdx][i] /= velocityWeight_[phaseIdx][i];
                // commit the phase velocity
                char name[512];
                snprintf(name, 512, "filterVelocity_%s", FluidSystem::phaseName(phaseIdx));

                DiscBaseOutputModule::attachVectorDofData_(baseWriter, velocity_[phaseIdx], name);
            }
        }

        if (potentialGradientOutput_()) {
            size_t numDof = this->simulator_.model().numGridDof();

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // first, divide the velocity field by the
                // respective finite volume's surface area
                for (unsigned i = 0; i < numDof; ++i)
                    potentialGradient_[phaseIdx][i] /= potentialWeight_[phaseIdx][i];
                // commit the phase velocity
                char name[512];
                snprintf(name, 512, "gradP_%s", FluidSystem::phaseName(phaseIdx));

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
    virtual bool needExtensiveQuantities() const final
    {
        return velocityOutput_() || potentialGradientOutput_();
    }

private:
    static bool extrusionFactorOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteExtrusionFactor);
        return val;
    }

    static bool pressureOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWritePressures);
        return val;
    }

    static bool densityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteDensities);
        return val;
    }

    static bool saturationOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSaturations);
        return val;
    }

    static bool mobilityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteMobilities);
        return val;
    }

    static bool relativePermeabilityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteRelativePermeabilities);
        return val;
    }

    static bool viscosityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteViscosities);
        return val;
    }

    static bool averageMolarMassOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteAverageMolarMasses);
        return val;
    }

    static bool porosityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWritePorosity);
        return val;
    }

    static bool intrinsicPermeabilityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteIntrinsicPermeabilities);
        return val;
    }

    static bool velocityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFilterVelocities);
        return val;
    }

    static bool potentialGradientOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWritePotentialGradients);
        return val;
    }

    ScalarBuffer extrusionFactor_;
    PhaseBuffer pressure_;
    PhaseBuffer density_;
    PhaseBuffer saturation_;
    PhaseBuffer mobility_;
    PhaseBuffer relativePermeability_;
    PhaseBuffer viscosity_;
    PhaseBuffer averageMolarMass_;

    ScalarBuffer porosity_;
    TensorBuffer intrinsicPermeability_;

    PhaseVectorBuffer velocity_;
    PhaseBuffer velocityWeight_;

    PhaseVectorBuffer potentialGradient_;
    PhaseBuffer potentialWeight_;
};

} // namespace Opm

#endif
