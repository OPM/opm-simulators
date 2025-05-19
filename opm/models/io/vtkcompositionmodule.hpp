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
 * \copydoc Opm::VtkCompositionModule
 */
#ifndef OPM_VTK_COMPOSITION_MODULE_HPP
#define OPM_VTK_COMPOSITION_MODULE_HPP

#include <opm/material/common/MathToolbox.hpp>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkcompositionparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the fluid composition
 *
 * This module deals with the following quantities:
 * - Mole fraction of a component in a fluid phase
 * - Mass fraction of a component in a fluid phase
 * - Molarity (i.e. molar concentration) of a component in a fluid phase
 * - Fugacity of all components
 * - FugacityCoefficient of all components in all phases
 */
template <class TypeTag>
class VtkCompositionModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    using BufferType = typename ParentType::BufferType;
    using ComponentBuffer = typename ParentType::ComponentBuffer;
    using PhaseComponentBuffer = typename ParentType::PhaseComponentBuffer;

public:
    explicit VtkCompositionModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        params_.read();
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        VtkCompositionParams::registerParameters();
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if (params_.moleFracOutput_) {
            this->resizePhaseComponentBuffer_(moleFrac_, BufferType::Dof);
        }
        if (params_.massFracOutput_) {
            this->resizePhaseComponentBuffer_(massFrac_, BufferType::Dof);
        }
        if (params_.totalMassFracOutput_) {
            this->resizeComponentBuffer_(totalMassFrac_, BufferType::Dof);
        }
        if (params_.totalMoleFracOutput_) {
            this->resizeComponentBuffer_(totalMoleFrac_, BufferType::Dof);
        }
        if (params_.molarityOutput_) {
            this->resizePhaseComponentBuffer_(molarity_, BufferType::Dof);
        }

        if (params_.fugacityOutput_) {
            this->resizeComponentBuffer_(fugacity_, BufferType::Dof);
        }
        if (params_.fugacityCoeffOutput_) {
            this->resizePhaseComponentBuffer_(fugacityCoeff_, BufferType::Dof);
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

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    if (params_.moleFracOutput_) {
                        moleFrac_[phaseIdx][compIdx][I] = Toolbox::value(fs.moleFraction(phaseIdx, compIdx));
                    }
                    if (params_.massFracOutput_) {
                        massFrac_[phaseIdx][compIdx][I] = Toolbox::value(fs.massFraction(phaseIdx, compIdx));
                    }
                    if (params_.molarityOutput_) {
                        molarity_[phaseIdx][compIdx][I] = Toolbox::value(fs.molarity(phaseIdx, compIdx));
                    }

                    if (params_.fugacityCoeffOutput_) {
                        fugacityCoeff_[phaseIdx][compIdx][I] =
                            Toolbox::value(fs.fugacityCoefficient(phaseIdx, compIdx));
                    }
                }
            }

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                if (params_.totalMassFracOutput_) {
                    Scalar compMass = 0;
                    Scalar totalMass = 0;
                    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                        totalMass += Toolbox::value(fs.density(phaseIdx)) * Toolbox::value(fs.saturation(phaseIdx));
                        compMass +=
                            Toolbox::value(fs.density(phaseIdx)) *
                            Toolbox::value(fs.saturation(phaseIdx)) *
                            Toolbox::value(fs.massFraction(phaseIdx, compIdx));
                    }
                    totalMassFrac_[compIdx][I] = compMass / totalMass;
                }
                if (params_.totalMoleFracOutput_) {
                    Scalar compMoles = 0;
                    Scalar totalMoles = 0;
                    for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                        totalMoles +=
                            Toolbox::value(fs.molarDensity(phaseIdx)) *
                            Toolbox::value(fs.saturation(phaseIdx));
                        compMoles +=
                            Toolbox::value(fs.molarDensity(phaseIdx)) *
                            Toolbox::value(fs.saturation(phaseIdx)) *
                            Toolbox::value(fs.moleFraction(phaseIdx, compIdx));
                    }
                    totalMoleFrac_[compIdx][I] = compMoles / totalMoles;
                }
                if (params_.fugacityOutput_) {
                    fugacity_[compIdx][I] = Toolbox::value(intQuants.fluidState().fugacity(/*phaseIdx=*/0, compIdx));
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

        if (params_.moleFracOutput_) {
            this->commitPhaseComponentBuffer_(baseWriter, "moleFrac_%s^%s",
                                              moleFrac_, BufferType::Dof);
        }
        if (params_.massFracOutput_) {
            this->commitPhaseComponentBuffer_(baseWriter, "massFrac_%s^%s",
                                              massFrac_, BufferType::Dof);
        }
        if (params_.molarityOutput_) {
            this->commitPhaseComponentBuffer_(baseWriter, "molarity_%s^%s",
                                              molarity_, BufferType::Dof);
        }
        if (params_.totalMassFracOutput_) {
            this->commitComponentBuffer_(baseWriter, "totalMassFrac^%s",
                                         totalMassFrac_, BufferType::Dof);
        }
        if (params_.totalMoleFracOutput_) {
            this->commitComponentBuffer_(baseWriter, "totalMoleFrac^%s",
                                         totalMoleFrac_, BufferType::Dof);
        }

        if (params_.fugacityOutput_) {
            this->commitComponentBuffer_(baseWriter, "fugacity^%s",
                                         fugacity_, BufferType::Dof);
        }
        if (params_.fugacityCoeffOutput_) {
            this->commitPhaseComponentBuffer_(baseWriter, "fugacityCoeff_%s^%s",
                                              fugacityCoeff_, BufferType::Dof);
        }
    }

private:
    VtkCompositionParams params_{};
    PhaseComponentBuffer moleFrac_{};
    PhaseComponentBuffer massFrac_{};
    PhaseComponentBuffer molarity_{};
    ComponentBuffer totalMassFrac_{};
    ComponentBuffer totalMoleFrac_{};

    ComponentBuffer fugacity_{};
    PhaseComponentBuffer fugacityCoeff_{};
};

} // namespace Opm

#endif // OPM_VTK_COMPOSITION_MODULE_HPP
