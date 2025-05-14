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
 *
 * \copydoc Opm::VtkDiffusionModule
 */
#ifndef OPM_VTK_DIFFUSION_MODULE_HPP
#define OPM_VTK_DIFFUSION_MODULE_HPP

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkdiffusionparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hpp>

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for quantities which make sense for models which
 *        incorperate molecular diffusion.
 *
 * This module deals with the following quantities:
 * - Molecular diffusion coefficients of all components in all fluid phases
 * - Effective molecular diffusion coefficients of the porous medium of all
 *   components in all fluid phases
 */
template <class TypeTag>
class VtkDiffusionModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

    using Toolbox = MathToolbox<Evaluation>;

    using BufferType = typename ParentType::BufferType;
    using PhaseBuffer = typename ParentType::PhaseBuffer;
    using PhaseComponentBuffer = typename ParentType::PhaseComponentBuffer;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

public:
    explicit VtkDiffusionModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        params_.read();
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        VtkDiffusionParams::registerParameters();
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if (params_.tortuosityOutput_) {
            this->resizePhaseBuffer_(tortuosity_, BufferType::Dof);
        }
        if (params_.diffusionCoefficientOutput_) {
            this->resizePhaseComponentBuffer_(diffusionCoefficient_, BufferType::Dof);
        }
        if (params_.effectiveDiffusionCoefficientOutput_) {
            this->resizePhaseComponentBuffer_(effectiveDiffusionCoefficient_, BufferType::Dof);
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

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (params_.tortuosityOutput_) {
                    tortuosity_[phaseIdx][I] = Toolbox::value(intQuants.tortuosity(phaseIdx));
                }
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    if (params_.diffusionCoefficientOutput_) {
                        diffusionCoefficient_[phaseIdx][compIdx][I] =
                            Toolbox::value(intQuants.diffusionCoefficient(phaseIdx, compIdx));
                    }
                    if (params_.effectiveDiffusionCoefficientOutput_) {
                        effectiveDiffusionCoefficient_[phaseIdx][compIdx][I] =
                            Toolbox::value(intQuants.effectiveDiffusionCoefficient(phaseIdx, compIdx));
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
        if (!dynamic_cast<VtkMultiWriter*>(&baseWriter)) {
            return;
        }

        if (params_.tortuosityOutput_) {
            this->commitPhaseBuffer_(baseWriter, "tortuosity", tortuosity_, BufferType::Dof);
        }
        if (params_.diffusionCoefficientOutput_) {
            this->commitPhaseComponentBuffer_(baseWriter, "diffusionCoefficient",
                                              diffusionCoefficient_, BufferType::Dof);
        }
        if (params_.effectiveDiffusionCoefficientOutput_) {
            this->commitPhaseComponentBuffer_(baseWriter,
                                              "effectiveDiffusionCoefficient",
                                              effectiveDiffusionCoefficient_,
                                              BufferType::Dof);
        }
    }

private:
    VtkDiffusionParams params_{};
    PhaseBuffer tortuosity_{};
    PhaseComponentBuffer diffusionCoefficient_{};
    PhaseComponentBuffer effectiveDiffusionCoefficient_{};
};

} // namespace Opm

#endif // OPM_VTK_DIFFUSION_MODULE_HPP
