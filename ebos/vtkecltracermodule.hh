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
 * \copydoc Opm::VtkEclTracerModule
 */
#ifndef EWOMS_VTK_ECL_TRACER_MODULE_HH
#define EWOMS_VTK_ECL_TRACER_MODULE_HH

#include <dune/common/fvector.hh>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>

#include <string>
#include <vector>

namespace Opm::Properties {

// create new type tag for the VTK tracer output
namespace TTag {
struct VtkEclTracer {};
}

// create the property tags needed for the tracer model
template<class TypeTag, class MyTypeTag>
struct VtkWriteEclTracerConcentration {
    using type = UndefinedProperty;
};

// set default values for what quantities to output
template<class TypeTag>
struct VtkWriteEclTracerConcentration<TypeTag, TTag::VtkEclTracer> {
    static constexpr bool value = false;
};

} // namespace Opm::Properties

namespace Opm {
    /*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the tracer model's parameters.
 */
    template <class TypeTag>
    class VtkEclTracerModule : public BaseOutputModule<TypeTag>
    {
        using ParentType = BaseOutputModule<TypeTag>;

        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

        using GridView = GetPropType<TypeTag, Properties::GridView>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

        static constexpr int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
        using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;


        using ScalarBuffer = typename ParentType::ScalarBuffer;

    public:
        VtkEclTracerModule(const Simulator& simulator)
            : ParentType(simulator)
        { }

        /*!
     * \brief Register all run-time parameters for the tracer VTK output
     * module.
     */
        static void registerParameters()
        {
            EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteEclTracerConcentration,
                                 "Include the tracer concentration "
                                 "in the VTK output files");
        }

        /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
        void allocBuffers()
        {
            if (eclTracerConcentrationOutput_()){
                const auto& tracerModel = this->simulator_.problem().tracerModel();
                eclTracerConcentration_.resize(tracerModel.numTracers());
                for (std::size_t tracerIdx = 0; tracerIdx < eclTracerConcentration_.size(); ++tracerIdx) {

                    this->resizeScalarBuffer_(eclTracerConcentration_[tracerIdx]);
                }
            }

        }

        /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
        void processElement(const ElementContext& elemCtx)
        {
            if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
                return;

            const auto& tracerModel = elemCtx.problem().tracerModel();

            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                if (eclTracerConcentrationOutput_()){
                    for (std::size_t tracerIdx  = 0; tracerIdx < eclTracerConcentration_.size(); ++tracerIdx) {
                        eclTracerConcentration_[tracerIdx][globalDofIdx] = tracerModel.tracerConcentration(tracerIdx, globalDofIdx);
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
            if (!vtkWriter)
                return;

            if (eclTracerConcentrationOutput_()){
                const auto& tracerModel = this->simulator_.problem().tracerModel();
                for (std::size_t tracerIdx = 0; tracerIdx < eclTracerConcentration_.size(); ++tracerIdx) {
                    const std::string tmp = "tracerConcentration_" + tracerModel.name(tracerIdx);
                    this->commitScalarBuffer_(baseWriter,tmp.c_str(), eclTracerConcentration_[tracerIdx]);
                }
            }



        }

    private:
        static bool eclTracerConcentrationOutput_()
        {
            static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteEclTracerConcentration);
            return val;
        }


        std::vector<ScalarBuffer> eclTracerConcentration_;
    };
} // namespace Opm

#endif
