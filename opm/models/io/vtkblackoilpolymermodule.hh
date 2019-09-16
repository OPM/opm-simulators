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
 * \copydoc Opm::VtkBlackOilPolymerModule
 */
#ifndef EWOMS_VTK_BLACK_OIL_POLYMER_MODULE_HH
#define EWOMS_VTK_BLACK_OIL_POLYMER_MODULE_HH

#include <opm/material/densead/Math.hpp>

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <dune/common/fvector.hh>

#include <cstdio>

BEGIN_PROPERTIES

// create new type tag for the VTK multi-phase output
NEW_TYPE_TAG(VtkBlackOilPolymer);

// create the property tags needed for the polymer output module
NEW_PROP_TAG(EnablePolymer);
NEW_PROP_TAG(EnableVtkOutput);
NEW_PROP_TAG(VtkWritePolymerConcentration);
NEW_PROP_TAG(VtkWritePolymerDeadPoreVolume);
NEW_PROP_TAG(VtkWritePolymerAdsorption);
NEW_PROP_TAG(VtkWritePolymerRockDensity);
NEW_PROP_TAG(VtkWritePolymerViscosityCorrection);
NEW_PROP_TAG(VtkWriteWaterViscosityCorrection);

// set default values for what quantities to output
SET_BOOL_PROP(VtkBlackOilPolymer, VtkWritePolymerConcentration, true);
SET_BOOL_PROP(VtkBlackOilPolymer, VtkWritePolymerDeadPoreVolume, true);
SET_BOOL_PROP(VtkBlackOilPolymer, VtkWritePolymerViscosityCorrection, true);
SET_BOOL_PROP(VtkBlackOilPolymer, VtkWriteWaterViscosityCorrection, true);
SET_BOOL_PROP(VtkBlackOilPolymer, VtkWritePolymerRockDensity, true);
SET_BOOL_PROP(VtkBlackOilPolymer, VtkWritePolymerAdsorption, true);

END_PROPERTIES

namespace Opm {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the black oil model's polymer related quantities.
 */
template <class TypeTag>
class VtkBlackOilPolymerModule : public BaseOutputModule<TypeTag>
{
    typedef BaseOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    typedef Opm::VtkMultiWriter<GridView, vtkFormat> VtkMultiWriter;

    enum { enablePolymer = GET_PROP_VALUE(TypeTag, EnablePolymer) };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;

public:
    VtkBlackOilPolymerModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        if (!enablePolymer)
            return;

        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePolymerConcentration,
                             "Include the concentration of the polymer component in the water phase "
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePolymerDeadPoreVolume,
                             "Include the fraction of the \"dead\" pore volume "
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePolymerRockDensity,
                             "Include the amount of already adsorbed polymer component"
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePolymerAdsorption,
                             "Include the adsorption rate of the polymer component"
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePolymerViscosityCorrection,
                             "Include the viscosity correction of the polymer component "
                             "in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteWaterViscosityCorrection,
                             "Include the viscosity correction of the water component "
                             "due to polymers in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        if (!enablePolymer)
            return;

        if (polymerConcentrationOutput_())
            this->resizeScalarBuffer_(polymerConcentration_);
        if (polymerDeadPoreVolumeOutput_())
            this->resizeScalarBuffer_(polymerDeadPoreVolume_);
        if (polymerRockDensityOutput_())
            this->resizeScalarBuffer_(polymerRockDensity_);
        if (polymerAdsorptionOutput_())
            this->resizeScalarBuffer_(polymerAdsorption_);
        if (polymerViscosityCorrectionOutput_())
            this->resizeScalarBuffer_(polymerViscosityCorrection_);
        if (waterViscosityCorrectionOutput_())
            this->resizeScalarBuffer_(waterViscosityCorrection_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        if (!enablePolymer)
            return;

        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            if (polymerConcentrationOutput_())
                polymerConcentration_[globalDofIdx] =
                    Opm::scalarValue(intQuants.polymerConcentration());

            if (polymerDeadPoreVolumeOutput_())
                polymerDeadPoreVolume_[globalDofIdx] =
                    Opm::scalarValue(intQuants.polymerDeadPoreVolume());

            if (polymerRockDensityOutput_())
                polymerRockDensity_[globalDofIdx] =
                    Opm::scalarValue(intQuants.polymerRockDensity());

            if (polymerAdsorptionOutput_())
                polymerAdsorption_[globalDofIdx] =
                    Opm::scalarValue(intQuants.polymerAdsorption());

            if (polymerViscosityCorrectionOutput_())
                polymerViscosityCorrection_[globalDofIdx] =
                    Opm::scalarValue(intQuants.polymerViscosityCorrection());

            if (waterViscosityCorrectionOutput_())
                waterViscosityCorrection_[globalDofIdx] =
                    Opm::scalarValue(intQuants.waterViscosityCorrection());
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

        if (!enablePolymer)
            return;

        if (polymerConcentrationOutput_())
            this->commitScalarBuffer_(baseWriter, "polymer concentration", polymerConcentration_);

        if (polymerDeadPoreVolumeOutput_())
            this->commitScalarBuffer_(baseWriter, "dead pore volume fraction", polymerDeadPoreVolume_);

        if (polymerRockDensityOutput_())
            this->commitScalarBuffer_(baseWriter, "polymer rock density", polymerRockDensity_);

        if (polymerAdsorptionOutput_())
            this->commitScalarBuffer_(baseWriter, "polymer adsorption", polymerAdsorption_);

        if (polymerViscosityCorrectionOutput_())
            this->commitScalarBuffer_(baseWriter, "polymer viscosity correction", polymerViscosityCorrection_);

        if (waterViscosityCorrectionOutput_())
            this->commitScalarBuffer_(baseWriter, "water viscosity correction", waterViscosityCorrection_);
    }

private:
    static bool polymerConcentrationOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWritePolymerConcentration);
        return val;
    }

    static bool polymerDeadPoreVolumeOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWritePolymerDeadPoreVolume);
        return val;
    }

    static bool polymerRockDensityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWritePolymerRockDensity);
        return val;
    }

    static bool polymerAdsorptionOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWritePolymerAdsorption);
        return val;
    }

    static bool polymerViscosityCorrectionOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWritePolymerViscosityCorrection);
        return val;
    }

    static bool waterViscosityCorrectionOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWritePolymerViscosityCorrection);
        return val;
    }

    ScalarBuffer polymerConcentration_;
    ScalarBuffer polymerDeadPoreVolume_;
    ScalarBuffer polymerRockDensity_;
    ScalarBuffer polymerAdsorption_;
    ScalarBuffer polymerViscosityCorrection_;
    ScalarBuffer waterViscosityCorrection_;
};
} // namespace Opm

#endif
