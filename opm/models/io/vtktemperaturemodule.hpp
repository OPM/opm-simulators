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
#ifndef OPM_VTK_TEMPERATURE_MODULE_HPP
#define OPM_VTK_TEMPERATURE_MODULE_HPP

#include <opm/material/common/MathToolbox.hpp>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/io/vtktemperatureparams.hpp>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

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

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

public:
    explicit VtkTemperatureModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        params_.read();
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        VtkTemperatureParams::registerParameters();
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if (params_.temperatureOutput_) {
            this->resizeScalarBuffer_(temperature_, ParentType::BufferType::Dof);
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx) override
    {
        using Toolbox = MathToolbox<Evaluation>;

        if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
            return;
        }

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            const unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            if (params_.temperatureOutput_) {
                temperature_[I] = Toolbox::value(fs.temperature(/*phaseIdx=*/0));
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

        if (params_.temperatureOutput_) {
            this->commitScalarBuffer_(baseWriter, "temperature", temperature_,
                                      ParentType::BufferType::Dof);
        }
    }

private:
    VtkTemperatureParams params_{};
    ScalarBuffer temperature_{};
};

} // namespace Opm

#endif // OPM_VTK_TEMPERATURE_MODULE_HPP
