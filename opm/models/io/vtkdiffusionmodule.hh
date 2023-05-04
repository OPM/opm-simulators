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
#ifndef EWOMS_VTK_DIFFUSION_MODULE_HH
#define EWOMS_VTK_DIFFUSION_MODULE_HH

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>

namespace Opm::Properties {

namespace TTag {

// create new type tag for the VTK output of the quantities for molecular
// diffusion
struct VtkDiffusion {};

} // namespace TTag

// create the property tags needed for the diffusion module
template<class TypeTag, class MyTypeTag>
struct VtkWriteTortuosities { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteDiffusionCoefficients { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteEffectiveDiffusionCoefficients { using type = UndefinedProperty; };

// set default values for what quantities to output
template<class TypeTag>
struct VtkWriteTortuosities<TypeTag, TTag::VtkDiffusion> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWriteDiffusionCoefficients<TypeTag, TTag::VtkDiffusion> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWriteEffectiveDiffusionCoefficients<TypeTag, TTag::VtkDiffusion> { static constexpr bool value = false; };

} // namespace Opm::Properties

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

    using PhaseComponentBuffer = typename ParentType::PhaseComponentBuffer;
    using PhaseBuffer = typename ParentType::PhaseBuffer;

    static const int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

public:
    VtkDiffusionModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteTortuosities,
                             "Include the tortuosity for each phase in the VTK "
                             "output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteDiffusionCoefficients,
                             "Include the molecular diffusion coefficients in "
                             "the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool,
                             VtkWriteEffectiveDiffusionCoefficients,
                             "Include the effective molecular diffusion "
                             "coefficients the medium in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (tortuosityOutput_())
            this->resizePhaseBuffer_(tortuosity_);
        if (diffusionCoefficientOutput_())
            this->resizePhaseComponentBuffer_(diffusionCoefficient_);
        if (effectiveDiffusionCoefficientOutput_())
            this->resizePhaseComponentBuffer_(effectiveDiffusionCoefficient_);
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

            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (tortuosityOutput_())
                    tortuosity_[phaseIdx][I] = Toolbox::value(intQuants.tortuosity(phaseIdx));
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                    if (diffusionCoefficientOutput_())
                        diffusionCoefficient_[phaseIdx][compIdx][I] =
                            Toolbox::value(intQuants.diffusionCoefficient(phaseIdx, compIdx));
                    if (effectiveDiffusionCoefficientOutput_())
                        effectiveDiffusionCoefficient_[phaseIdx][compIdx][I] =
                            Toolbox::value(intQuants.effectiveDiffusionCoefficient(phaseIdx, compIdx));
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

        if (tortuosityOutput_())
            this->commitPhaseBuffer_(baseWriter, "tortuosity", tortuosity_);
        if (diffusionCoefficientOutput_())
            this->commitPhaseComponentBuffer_(baseWriter, "diffusionCoefficient",
                                              diffusionCoefficient_);
        if (effectiveDiffusionCoefficientOutput_())
            this->commitPhaseComponentBuffer_(baseWriter,
                                              "effectiveDiffusionCoefficient",
                                              effectiveDiffusionCoefficient_);
    }

private:
    static bool tortuosityOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteTortuosities);
        return val;
    }

    static bool diffusionCoefficientOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteDiffusionCoefficients);
        return val;
    }

    static bool effectiveDiffusionCoefficientOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteEffectiveDiffusionCoefficients);
        return val;
    }

    PhaseBuffer tortuosity_;
    PhaseComponentBuffer diffusionCoefficient_;
    PhaseComponentBuffer effectiveDiffusionCoefficient_;
};

} // namespace Opm

#endif
