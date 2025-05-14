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
 * \copydoc Opm::VtkTracerModule
 */
#ifndef OPM_VTK_TRACER_MODULE_HPP
#define OPM_VTK_TRACER_MODULE_HPP

#include <dune/common/fvector.hh>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/discretization/common/fvbaseparameters.hh>
#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <string>
#include <vector>

namespace Opm::Parameters {

// set default values for what quantities to output
struct VtkWriteTracerConcentration { static constexpr bool value = false; };

} // namespace Opm::Parameters

namespace Opm {

    /*!
     * \ingroup Vtk
     *
     * \brief VTK output module for the tracer model's parameters.
     */
    template <class TypeTag>
    class VtkTracerModule : public BaseOutputModule<TypeTag>
    {
        using ParentType = BaseOutputModule<TypeTag>;

        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
        using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

        using BufferType = typename ParentType::BufferType;
        using ScalarBuffer = typename ParentType::ScalarBuffer;

    public:
        explicit VtkTracerModule(const Simulator& simulator)
            : ParentType(simulator)
        {}

        /*!
         * \brief Register all run-time parameters for the tracer VTK output
         * module.
         */
        static void registerParameters()
        {
            Parameters::Register<Parameters::VtkWriteTracerConcentration>
                ("Include the tracer concentration in the VTK output files");
        }

        /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
        void allocBuffers() override
        {
            if (eclTracerConcentrationOutput_()) {
                const auto& tracerModel = this->simulator_.problem().tracerModel();
                eclFreeTracerConcentration_.resize(tracerModel.numTracers());
                eclSolTracerConcentration_.resize(tracerModel.numTracers());
                const auto& enableSolTracers = tracerModel.enableSolTracers();

                for (std::size_t tracerIdx = 0; tracerIdx < eclFreeTracerConcentration_.size(); ++tracerIdx) {
                    this->resizeScalarBuffer_(eclFreeTracerConcentration_[tracerIdx], BufferType::Dof);
                    if (enableSolTracers[tracerIdx]) {
                        this->resizeScalarBuffer_(eclSolTracerConcentration_[tracerIdx], BufferType::Dof);
                    }
                }
            }
        }

        /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
        void processElement(const ElementContext& elemCtx) override
        {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            if (eclTracerConcentrationOutput_()) {
                const auto& tracerModel = elemCtx.problem().tracerModel();
                const auto& enableSolTracers = tracerModel.enableSolTracers();

                for (std::size_t tracerIdx  = 0; tracerIdx < eclFreeTracerConcentration_.size(); ++tracerIdx) {
                    // free tracer
                    for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                        const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                        eclFreeTracerConcentration_[tracerIdx][globalDofIdx] =
                            tracerModel.freeTracerConcentration(tracerIdx, globalDofIdx);
                    }
                    // solution tracer (only if it exist)
                    if (enableSolTracers[tracerIdx]) {
                        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                            const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                            eclSolTracerConcentration_[tracerIdx][globalDofIdx] =
                                tracerModel.solTracerConcentration(tracerIdx, globalDofIdx);
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
            if (!dynamic_cast<VtkMultiWriter*>(&baseWriter))
                return;

            if (eclTracerConcentrationOutput_()){
                const auto& tracerModel = this->simulator_.problem().tracerModel();
                const auto& enableSolTracers = tracerModel.enableSolTracers();

                for (std::size_t tracerIdx = 0; tracerIdx < eclFreeTracerConcentration_.size(); ++tracerIdx) {
                    const std::string tmp = "freeTracerConcentration_" + tracerModel.name(tracerIdx);
                    this->commitScalarBuffer_(baseWriter, tmp.c_str(),
                                              eclFreeTracerConcentration_[tracerIdx], BufferType::Dof);
                    if (enableSolTracers[tracerIdx]) {
                        const std::string tmp2 = "solTracerConcentration_" + tracerModel.name(tracerIdx);
                        this->commitScalarBuffer_(baseWriter, tmp2.c_str(),
                                                  eclSolTracerConcentration_[tracerIdx], BufferType::Dof);
                    }
                }
            }
        }

    private:
        static bool eclTracerConcentrationOutput_()
        {
            static bool val = Parameters::Get<Parameters::VtkWriteTracerConcentration>();
            return val;
        }

        std::vector<ScalarBuffer> eclFreeTracerConcentration_;
        std::vector<ScalarBuffer> eclSolTracerConcentration_;
    };

} // namespace Opm

#endif // OPM_VTK_TRACER_MODULE_HPP
