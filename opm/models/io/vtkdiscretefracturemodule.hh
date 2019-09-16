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
#ifndef EWOMS_VTK_DISCRETE_FRACTURE_MODULE_HH
#define EWOMS_VTK_DISCRETE_FRACTURE_MODULE_HH

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>

#include <cstdio>

BEGIN_PROPERTIES

// create new type tag for the VTK multi-phase output
NEW_TYPE_TAG(VtkDiscreteFracture);

// create the property tags needed for the multi phase module
NEW_PROP_TAG(Vanguard);
NEW_PROP_TAG(VtkWriteFractureSaturations);
NEW_PROP_TAG(VtkWriteFractureMobilities);
NEW_PROP_TAG(VtkWriteFractureRelativePermeabilities);
NEW_PROP_TAG(VtkWriteFracturePorosity);
NEW_PROP_TAG(VtkWriteFractureIntrinsicPermeabilities);
NEW_PROP_TAG(VtkWriteFractureFilterVelocities);
NEW_PROP_TAG(VtkWriteFractureVolumeFraction);
NEW_PROP_TAG(VtkOutputFormat);
NEW_PROP_TAG(EnableVtkOutput);
NEW_PROP_TAG(DiscBaseOutputModule);

// set default values for what quantities to output
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureSaturations, true);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureMobilities, false);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureRelativePermeabilities, true);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFracturePorosity, true);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureIntrinsicPermeabilities, false);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureFilterVelocities, false);
SET_BOOL_PROP(VtkDiscreteFracture, VtkWriteFractureVolumeFraction, true);

END_PROPERTIES

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
    typedef BaseOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, Vanguard) Vanguard;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef typename GET_PROP_TYPE(TypeTag, DiscBaseOutputModule) DiscBaseOutputModule;

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    typedef Opm::VtkMultiWriter<GridView, vtkFormat> VtkMultiWriter;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::PhaseBuffer PhaseBuffer;
    typedef typename ParentType::PhaseVectorBuffer PhaseVectorBuffer;

public:
    VtkDiscreteFractureModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFractureSaturations,
                             "Include the phase saturations in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFractureMobilities,
                             "Include the phase mobilities in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFractureRelativePermeabilities,
                             "Include the phase relative permeabilities in the "
                             "VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFracturePorosity,
                             "Include the porosity in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFractureIntrinsicPermeabilities,
                             "Include the intrinsic permeability in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFractureFilterVelocities,
                             "Include in the filter velocities of the phases "
                             "the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFractureVolumeFraction,
                             "Add the fraction of the total volume which is "
                             "occupied by fractures in the VTK output");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (saturationOutput_())
            this->resizePhaseBuffer_(fractureSaturation_);
        if (mobilityOutput_())
            this->resizePhaseBuffer_(fractureMobility_);
        if (relativePermeabilityOutput_())
            this->resizePhaseBuffer_(fractureRelativePermeability_);

        if (porosityOutput_())
            this->resizeScalarBuffer_(fracturePorosity_);
        if (intrinsicPermeabilityOutput_())
            this->resizeScalarBuffer_(fractureIntrinsicPermeability_);
        if (volumeFractionOutput_())
            this->resizeScalarBuffer_(fractureVolumeFraction_);

        if (velocityOutput_()) {
            size_t nDof = this->simulator_.model().numGridDof();
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                fractureVelocity_[phaseIdx].resize(nDof);
                for (unsigned dofIdx = 0; dofIdx < nDof; ++dofIdx) {
                    fractureVelocity_[phaseIdx][dofIdx].resize(dimWorld);
                    fractureVelocity_[phaseIdx][dofIdx] = 0.0;
                }
            }
            this->resizePhaseBuffer_(fractureVelocityWeight_);
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        const auto& fractureMapper = elemCtx.simulator().vanguard().fractureMapper();

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            if (!fractureMapper.isFractureVertex(I))
                continue;

            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fractureFluidState();

            if (porosityOutput_()) {
                Opm::Valgrind::CheckDefined(intQuants.fracturePorosity());
                fracturePorosity_[I] = intQuants.fracturePorosity();
            }
            if (intrinsicPermeabilityOutput_()) {
                const auto& K = intQuants.fractureIntrinsicPermeability();
                fractureIntrinsicPermeability_[I] = K[0][0];
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (saturationOutput_()) {
                    Opm::Valgrind::CheckDefined(fs.saturation(phaseIdx));
                    fractureSaturation_[phaseIdx][I] = fs.saturation(phaseIdx);
                }
                if (mobilityOutput_()) {
                    Opm::Valgrind::CheckDefined(intQuants.fractureMobility(phaseIdx));
                    fractureMobility_[phaseIdx][I] = intQuants.fractureMobility(phaseIdx);
                }
                if (relativePermeabilityOutput_()) {
                    Opm::Valgrind::CheckDefined(intQuants.fractureRelativePermeability(phaseIdx));
                    fractureRelativePermeability_[phaseIdx][I] =
                        intQuants.fractureRelativePermeability(phaseIdx);
                }
                if (volumeFractionOutput_()) {
                    Opm::Valgrind::CheckDefined(intQuants.fractureVolume());
                    fractureVolumeFraction_[I] += intQuants.fractureVolume();
                }
            }
        }

        if (velocityOutput_()) {
            // calculate velocities if requested by the simulator
            for (unsigned scvfIdx = 0; scvfIdx < elemCtx.numInteriorFaces(/*timeIdx=*/0); ++ scvfIdx) {
                const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, /*timeIdx=*/0);

                unsigned i = extQuants.interiorIndex();
                unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);

                unsigned j = extQuants.exteriorIndex();
                unsigned J = elemCtx.globalSpaceIndex(j, /*timeIdx=*/0);

                if (!fractureMapper.isFractureEdge(I, J))
                    continue;

                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    Scalar weight =
                        std::max<Scalar>(1e-16, std::abs(extQuants.fractureVolumeFlux(phaseIdx)));
                    Opm::Valgrind::CheckDefined(extQuants.extrusionFactor());
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
    void commitBuffers(BaseOutputWriter& baseWriter)
    {
        VtkMultiWriter *vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
        if (!vtkWriter) {
            return;
        }

        if (saturationOutput_())
            this->commitPhaseBuffer_(baseWriter, "fractureSaturation_%s", fractureSaturation_);
        if (mobilityOutput_())
            this->commitPhaseBuffer_(baseWriter, "fractureMobility_%s", fractureMobility_);
        if (relativePermeabilityOutput_())
            this->commitPhaseBuffer_(baseWriter, "fractureRelativePerm_%s", fractureRelativePermeability_);

        if (porosityOutput_())
            this->commitScalarBuffer_(baseWriter, "fracturePorosity", fracturePorosity_);
        if (intrinsicPermeabilityOutput_())
            this->commitScalarBuffer_(baseWriter, "fractureIntrinsicPerm", fractureIntrinsicPermeability_);
        if (volumeFractionOutput_()) {
            // divide the fracture volume by the total volume of the finite volumes
            for (unsigned I = 0; I < fractureVolumeFraction_.size(); ++I)
                fractureVolumeFraction_[I] /= this->simulator_.model().dofTotalVolume(I);
            this->commitScalarBuffer_(baseWriter, "fractureVolumeFraction", fractureVolumeFraction_);
        }

        if (velocityOutput_()) {
            size_t nDof = this->simulator_.model().numGridDof();

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // first, divide the velocity field by the
                // respective finite volume's surface area
                for (unsigned dofIdx = 0; dofIdx < nDof; ++dofIdx)
                    fractureVelocity_[phaseIdx][dofIdx] /=
                        std::max<Scalar>(1e-20, fractureVelocityWeight_[phaseIdx][dofIdx]);
                // commit the phase velocity
                char name[512];
                snprintf(name, 512, "fractureFilterVelocity_%s", FluidSystem::phaseName(phaseIdx));

                DiscBaseOutputModule::attachVectorDofData_(baseWriter, fractureVelocity_[phaseIdx], name);
            }
        }
    }

private:
    static bool saturationOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFractureSaturations);
        return val;
    }

    static bool mobilityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFractureMobilities);
        return val;
    }

    static bool relativePermeabilityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFractureRelativePermeabilities);
        return val;
    }

    static bool porosityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFracturePorosity);
        return val;
    }

    static bool intrinsicPermeabilityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFractureIntrinsicPermeabilities);
        return val;
    }

    static bool volumeFractionOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFractureVolumeFraction);
        return val;
    }

    static bool velocityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFractureFilterVelocities);
        return val;
    }

    PhaseBuffer fractureSaturation_;
    PhaseBuffer fractureMobility_;
    PhaseBuffer fractureRelativePermeability_;

    ScalarBuffer fracturePorosity_;
    ScalarBuffer fractureVolumeFraction_;
    ScalarBuffer fractureIntrinsicPermeability_;

    PhaseVectorBuffer fractureVelocity_;
    PhaseBuffer fractureVelocityWeight_;

    PhaseVectorBuffer potentialGradient_;
    PhaseBuffer potentialWeight_;
};

} // namespace Opm

#endif
