// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012-2013 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Ewoms::VcfvVtkBlackOilModule
 */
#ifndef EWOMS_VCFV_VTK_BLACK_OIL_MODULE_HH
#define EWOMS_VCFV_VTK_BLACK_OIL_MODULE_HH

#include "vcfvvtkoutputmodule.hh"

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <dune/common/fvector.hh>

#include <cstdio>

namespace Ewoms {
namespace Properties {
// create new type tag for the VTK multi-phase output
NEW_TYPE_TAG(VtkBlackOil);

// create the property tags needed for the multi phase module
NEW_PROP_TAG(VtkWriteGasFormationFactor);
NEW_PROP_TAG(VtkWriteGasFormationVolumeFactor);
NEW_PROP_TAG(VtkWriteOilFormationVolumeFactor);
NEW_PROP_TAG(VtkWriteOilSaturationPressure);

// set default values for what quantities to output
SET_BOOL_PROP(VtkBlackOil, VtkWriteGasFormationFactor, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteGasFormationVolumeFactor, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteOilFormationVolumeFactor, false);
SET_BOOL_PROP(VtkBlackOil, VtkWriteOilSaturationPressure, false);
}

/*!
 * \ingroup VcfvVtk
 *
 * \brief VTK output module for the black oil model's parameters.
 */
template<class TypeTag>
class VcfvVtkBlackOilModule : public VcfvVtkOutputModule<TypeTag>
{
    typedef VcfvVtkOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Ewoms::VtkMultiWriter<GridView> VtkMultiWriter;

    enum { oPhaseIdx = FluidSystem::oPhaseIdx };
    enum { gCompIdx = FluidSystem::gCompIdx };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;

public:
    VcfvVtkBlackOilModule(const Problem &problem)
        : ParentType(problem)
    { }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output module.
     */
    static void registerParameters()
    {
        REGISTER_PARAM(TypeTag, bool, VtkWriteGasFormationFactor, "Include the gas formation factor in the VTK output files");
        REGISTER_PARAM(TypeTag, bool, VtkWriteGasFormationVolumeFactor, "Include the gas formation volume factor in the VTK output files");
        REGISTER_PARAM(TypeTag, bool, VtkWriteOilFormationVolumeFactor, "Include the oil formation volume factor in the VTK output files");
        REGISTER_PARAM(TypeTag, bool, VtkWriteOilSaturationPressure, "Include the saturation pressure of oil in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers(VtkMultiWriter &writer)
    {
        if (gasFormationFactorOutput_()) this->resizeScalarBuffer_(gasFormationFactor_);
        if (gasFormationVolumeFactorOutput_()) this->resizeScalarBuffer_(gasFormationVolumeFactor_);
        if (oilFormationVolumeFactorOutput_()) this->resizeScalarBuffer_(oilFormationVolumeFactor_);
        if (oilSaturationPressureOutput_()) this->resizeScalarBuffer_(oilSaturationPressure_);
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const ElementContext &elemCtx)
    {
        for (int i = 0; i < elemCtx.numScv(); ++i) {
            const auto &fs = elemCtx.volVars(/*spaceIdx=*/i, /*timeIdx=*/0).fluidState();
            int I = elemCtx.globalSpaceIndex(/*spaceIdx=*/i, /*timeIdx=*/0);
            Scalar po = fs.pressure(oPhaseIdx);
            Scalar X_oG = fs.massFraction(oPhaseIdx, gCompIdx);

            if (gasFormationFactorOutput_()) gasFormationFactor_[I] = FluidSystem::gasFormationFactor(po);
            if (gasFormationVolumeFactorOutput_()) gasFormationVolumeFactor_[I] = FluidSystem::gasFormationVolumeFactor(po);
            if (oilFormationVolumeFactorOutput_()) oilFormationVolumeFactor_[I] = FluidSystem::oilFormationVolumeFactor(po);
            if (oilSaturationPressureOutput_()) oilSaturationPressure_[I] = FluidSystem::oilSaturationPressure(X_oG);
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(VtkMultiWriter &writer)
    {
        if (gasFormationFactorOutput_()) this->commitScalarBuffer_(writer, "R_s", gasFormationFactor_);
        if (gasFormationVolumeFactorOutput_()) this->commitScalarBuffer_(writer, "B_g", gasFormationVolumeFactor_);
        if (oilFormationVolumeFactorOutput_()) this->commitScalarBuffer_(writer, "B_o", oilFormationVolumeFactor_);
        if (oilSaturationPressureOutput_()) this->commitScalarBuffer_(writer, "pressure_sat,o", oilSaturationPressure_);
    }

private:
    static bool gasFormationFactorOutput_()
    { return GET_PARAM(TypeTag, bool, VtkWriteGasFormationFactor); }

    static bool gasFormationVolumeFactorOutput_()
    { return GET_PARAM(TypeTag, bool, VtkWriteGasFormationVolumeFactor); }

    static bool oilFormationVolumeFactorOutput_()
    { return GET_PARAM(TypeTag, bool, VtkWriteOilFormationVolumeFactor); }

    static bool oilSaturationPressureOutput_()
    { return GET_PARAM(TypeTag, bool, VtkWriteOilSaturationPressure); }

    ScalarBuffer gasFormationFactor_;
    ScalarBuffer gasFormationVolumeFactor_;
    ScalarBuffer oilFormationVolumeFactor_;
    ScalarBuffer oilSaturationPressure_;
};

}

#endif
