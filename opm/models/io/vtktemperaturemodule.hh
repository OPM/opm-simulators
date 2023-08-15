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
 * \copydoc Opm::VtkTemperatureModule
 */
#ifndef EWOMS_VTK_TEMPERATURE_MODULE_HH
#define EWOMS_VTK_TEMPERATURE_MODULE_HH

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/material/common/MathToolbox.hpp>

namespace Opm::Properties {

namespace TTag {

// create new type tag for the VTK temperature output
struct VtkTemperature {};

} // namespace TTag

// create the property tags needed for the temperature module
template<class TypeTag, class MyTypeTag>
struct VtkWriteTemperature { using type = UndefinedProperty; };

// set default values for what quantities to output
template<class TypeTag>
struct VtkWriteTemperature<TypeTag, TTag::VtkTemperature> { static constexpr bool value = true; };

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the temperature in which assume
 *        thermal equilibrium
 */
template<class TypeTag>
class VtkTemperatureModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

    using ScalarBuffer = typename ParentType::ScalarBuffer;

    static const int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

public:
    VtkTemperatureModule(const Simulator& simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteTemperature,
                             "Include the temperature in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (temperatureOutput_()) this->resizeScalarBuffer_(temperature_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        using Toolbox = MathToolbox<Evaluation>;

        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            if (temperatureOutput_())
                temperature_[I] = Toolbox::value(fs.temperature(/*phaseIdx=*/0));
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

        if (temperatureOutput_())
            this->commitScalarBuffer_(baseWriter, "temperature", temperature_);
    }

private:
    static bool temperatureOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteTemperature);
        return val;
    }

    ScalarBuffer temperature_;
};

} // namespace Opm

#endif
