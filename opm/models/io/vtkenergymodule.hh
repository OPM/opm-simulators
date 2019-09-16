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
 * \copydoc Opm::VtkEnergyModule
 */
#ifndef EWOMS_VTK_ENERGY_MODULE_HH
#define EWOMS_VTK_ENERGY_MODULE_HH

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/material/common/MathToolbox.hpp>

BEGIN_PROPERTIES

// create new type tag for the VTK energy output
NEW_TYPE_TAG(VtkEnergy);

// create the property tags needed for the energy module
NEW_PROP_TAG(VtkWriteSolidInternalEnergy);
NEW_PROP_TAG(VtkWriteThermalConductivity);
NEW_PROP_TAG(VtkWriteInternalEnergies);
NEW_PROP_TAG(VtkWriteEnthalpies);
NEW_PROP_TAG(VtkOutputFormat);
NEW_PROP_TAG(EnableVtkOutput);

// set default values for what quantities to output
SET_BOOL_PROP(VtkEnergy, VtkWriteSolidInternalEnergy, false);
SET_BOOL_PROP(VtkEnergy, VtkWriteThermalConductivity, false);
SET_BOOL_PROP(VtkEnergy, VtkWriteInternalEnergies, false);
SET_BOOL_PROP(VtkEnergy, VtkWriteEnthalpies, false);

END_PROPERTIES

namespace Opm {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for quantities which make sense for models which
 *        assume thermal equilibrium.
 *
 * This module deals with the following quantities:
 * - Specific enthalpy of all fluid phases
 * - Specific internal energy of all fluid phases
 * - Volumetric internal energy of the solid phase
 * - Total thermal conductivity, i.e. the conductivity of the solid and all fluid phases
 *   combined
 */
template <class TypeTag>
class VtkEnergyModule : public BaseOutputModule<TypeTag>
{
    typedef BaseOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::PhaseBuffer PhaseBuffer;

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename Opm::MathToolbox<Evaluation> Toolbox;
    typedef Opm::VtkMultiWriter<GridView, vtkFormat> VtkMultiWriter;

public:
    VtkEnergyModule(const Simulator& simulator)
        : ParentType(simulator)
    {
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSolidInternalEnergy,
                             "Include the volumetric internal energy of solid"
                             "matrix in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteThermalConductivity,
                             "Include the total thermal conductivity of the"
                             "medium in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteEnthalpies,
                             "Include the specific enthalpy of the phases in "
                             "the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteInternalEnergies,
                             "Include the specific internal energy of the "
                             "phases in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (enthalpyOutput_())
            this->resizePhaseBuffer_(enthalpy_);
        if (internalEnergyOutput_())
            this->resizePhaseBuffer_(internalEnergy_);

        if (solidInternalEnergyOutput_())
            this->resizeScalarBuffer_(solidInternalEnergy_);
        if (thermalConductivityOutput_())
            this->resizeScalarBuffer_(thermalConductivity_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quanties relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            if (solidInternalEnergyOutput_())
                solidInternalEnergy_[I] = Toolbox::value(intQuants.solidInternalEnergy());
            if (thermalConductivityOutput_())
                thermalConductivity_[I] = Toolbox::value(intQuants.thermalConductivity());

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (enthalpyOutput_())
                    enthalpy_[phaseIdx][I] = Toolbox::value(fs.enthalpy(phaseIdx));
                if (internalEnergyOutput_())
                    internalEnergy_[phaseIdx][I] = Toolbox::value(fs.internalEnergy(phaseIdx));
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

        if (solidInternalEnergyOutput_())
            this->commitScalarBuffer_(baseWriter, "internalEnergySolid", solidInternalEnergy_);
        if (thermalConductivityOutput_())
            this->commitScalarBuffer_(baseWriter, "thermalConductivity", thermalConductivity_);

        if (enthalpyOutput_())
            this->commitPhaseBuffer_(baseWriter, "enthalpy_%s", enthalpy_);
        if (internalEnergyOutput_())
            this->commitPhaseBuffer_(baseWriter, "internalEnergy_%s", internalEnergy_);
    }

private:
    static bool solidInternalEnergyOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSolidInternalEnergy);
        return val;
    }

    static bool thermalConductivityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteThermalConductivity);
        return val;
    }

    static bool enthalpyOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteEnthalpies);
        return val;
    }

    static bool internalEnergyOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteInternalEnergies);
        return val;
    }

    PhaseBuffer enthalpy_;
    PhaseBuffer internalEnergy_;

    ScalarBuffer thermalConductivity_;
    ScalarBuffer solidInternalEnergy_;
};

} // namespace Opm

#endif
