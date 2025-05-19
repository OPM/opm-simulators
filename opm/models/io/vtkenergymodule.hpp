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
#ifndef OPM_VTK_ENERGY_MODULE_HPP
#define OPM_VTK_ENERGY_MODULE_HPP

#include <opm/material/common/MathToolbox.hpp>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkenergyparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

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
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    using BufferType = typename ParentType::BufferType;
    using ScalarBuffer = typename ParentType::ScalarBuffer;
    using PhaseBuffer = typename ParentType::PhaseBuffer;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };

    using Toolbox = typename Opm::MathToolbox<Evaluation>;
    using VtkMultiWriter = Opm::VtkMultiWriter<GridView, vtkFormat>;

public:
    explicit VtkEnergyModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        params_.read();
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        VtkEnergyParams::registerParameters();
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if (params_.enthalpyOutput_) {
            this->resizePhaseBuffer_(enthalpy_, BufferType::Dof);
        }
        if (params_.internalEnergyOutput_) {
            this->resizePhaseBuffer_(internalEnergy_, BufferType::Dof);
        }

        if (params_.solidInternalEnergyOutput_) {
            this->resizeScalarBuffer_(solidInternalEnergy_, BufferType::Dof);
        }
        if (params_.thermalConductivityOutput_) {
            this->resizeScalarBuffer_(thermalConductivity_, BufferType::Dof);
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quanties relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx) override
    {
        if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
            return;
        }

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            const unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            if (params_.solidInternalEnergyOutput_) {
                solidInternalEnergy_[I] = Toolbox::value(intQuants.solidInternalEnergy());
            }
            if (params_.thermalConductivityOutput_) {
                thermalConductivity_[I] = Toolbox::value(intQuants.thermalConductivity());
            }

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (params_.enthalpyOutput_) {
                    enthalpy_[phaseIdx][I] = Toolbox::value(fs.enthalpy(phaseIdx));
                }
                if (params_.internalEnergyOutput_) {
                    internalEnergy_[phaseIdx][I] = Toolbox::value(fs.internalEnergy(phaseIdx));
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

        if (params_.solidInternalEnergyOutput_) {
            this->commitScalarBuffer_(baseWriter, "internalEnergySolid",
                                      solidInternalEnergy_, BufferType::Dof);
        }
        if (params_.thermalConductivityOutput_) {
            this->commitScalarBuffer_(baseWriter, "thermalConductivity",
                                      thermalConductivity_, BufferType::Dof);
        }

        if (params_.enthalpyOutput_) {
            this->commitPhaseBuffer_(baseWriter, "enthalpy_%s",
                                     enthalpy_, BufferType::Dof);
        }
        if (params_.internalEnergyOutput_) {
            this->commitPhaseBuffer_(baseWriter, "internalEnergy_%s",
                                     internalEnergy_, BufferType::Dof);
        }
    }

private:
    VtkEnergyParams params_{};
    PhaseBuffer enthalpy_{};
    PhaseBuffer internalEnergy_{};

    ScalarBuffer thermalConductivity_{};
    ScalarBuffer solidInternalEnergy_{};
};

} // namespace Opm

#endif // OPM_VTK_ENERGY_MODULE_HPP
