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
 * \copydoc Opm::VtkBlackOilEnergyModule
 */
#ifndef OPM_VTK_BLACK_OIL_ENERGY_MODULE_HPP
#define OPM_VTK_BLACK_OIL_ENERGY_MODULE_HPP

#include <dune/common/fvector.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkblackoilenergyparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the black oil model's energy related quantities.
 */
template <class TypeTag>
class VtkBlackOilEnergyModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };

    using BufferType = typename ParentType::BufferType;
    using PhaseBuffer = typename ParentType::PhaseBuffer;
    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    explicit VtkBlackOilEnergyModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        if constexpr (enableEnergy) {
            params_.read();
        }
    }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        if constexpr (enableEnergy) {
            VtkBlackoilEnergyParams::registerParameters();
        }
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if constexpr (enableEnergy) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            if (params_.rockInternalEnergyOutput_) {
                this->resizeScalarBuffer_(rockInternalEnergy_, BufferType::Dof);
            }
            if (params_.totalThermalConductivityOutput_) {
                this->resizeScalarBuffer_(totalThermalConductivity_, BufferType::Dof);
            }
            if (params_.fluidInternalEnergiesOutput_) {
                this->resizePhaseBuffer_(fluidInternalEnergies_, BufferType::Dof);
            }
            if (params_.fluidEnthalpiesOutput_) {
                this->resizePhaseBuffer_(fluidEnthalpies_, BufferType::Dof);
            }
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx) override
    {
        if constexpr (enableEnergy) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
                const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                if (params_.rockInternalEnergyOutput_) {
                    rockInternalEnergy_[globalDofIdx] =
                        scalarValue(intQuants.rockInternalEnergy());
                }

                if (params_.totalThermalConductivityOutput_) {
                    totalThermalConductivity_[globalDofIdx] =
                        scalarValue(intQuants.totalThermalConductivity());
                }

                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    if (FluidSystem::phaseIsActive(phaseIdx)) {
                        if (params_.fluidInternalEnergiesOutput_) {
                            fluidInternalEnergies_[phaseIdx][globalDofIdx] =
                                scalarValue(intQuants.fluidState().internalEnergy(phaseIdx));
                        }

                        if (params_.fluidEnthalpiesOutput_) {
                            fluidEnthalpies_[phaseIdx][globalDofIdx] =
                                scalarValue(intQuants.fluidState().enthalpy(phaseIdx));
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter) override
    {
        if constexpr (enableEnergy) {
            if (!dynamic_cast<VtkMultiWriter*>(&baseWriter)) {
                return;
            }

            if (params_.rockInternalEnergyOutput_) {
                this->commitScalarBuffer_(baseWriter, "volumetric internal energy rock",
                                          rockInternalEnergy_, BufferType::Dof);
            }

            if (params_.totalThermalConductivityOutput_) {
                this->commitScalarBuffer_(baseWriter, "total thermal conductivity",
                                          totalThermalConductivity_, BufferType::Dof);
            }

            if (params_.fluidInternalEnergiesOutput_) {
                this->commitPhaseBuffer_(baseWriter, "internal energy_%s",
                                         fluidInternalEnergies_, BufferType::Dof);
            }

            if (params_.fluidEnthalpiesOutput_) {
                this->commitPhaseBuffer_(baseWriter, "enthalpy_%s",
                                         fluidEnthalpies_, BufferType::Dof);
            }
        }
    }

private:
    VtkBlackoilEnergyParams params_{};

    ScalarBuffer rockInternalEnergy_{};
    ScalarBuffer totalThermalConductivity_{};
    PhaseBuffer fluidInternalEnergies_{};
    PhaseBuffer fluidEnthalpies_{};
};

} // namespace Opm

#endif // OPM_VTK_BLACKOIL_ENERGY_MODULE_HPP
