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
 * \copydoc Opm::VtkPrimaryVarsModule
 */
#ifndef OPM_VTK_PRIMARY_VARS_MODULE_HPP
#define OPM_VTK_PRIMARY_VARS_MODULE_HPP

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/io/vtkprimaryvarsparams.hpp>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the fluid composition
 */
template<class TypeTag>
class VtkPrimaryVarsModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    static const int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    using ScalarBuffer = typename ParentType::ScalarBuffer;
    using EqBuffer = typename ParentType::EqBuffer;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };

public:
    VtkPrimaryVarsModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        params_.read();
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        VtkPrimaryVarsParams::registerParameters();
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (params_.primaryVarsOutput_) {
            this->resizeEqBuffer_(primaryVars_);
        }
        if (params_.processRankOutput_) {
            this->resizeScalarBuffer_(processRank_,
                                      /*bufferType=*/ParentType::ElementBuffer);
        }
        if (params_.dofIndexOutput_) {
            this->resizeScalarBuffer_(dofIndex_);
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
            return;
        }

        const auto& elementMapper = elemCtx.model().elementMapper();
        unsigned elemIdx = static_cast<unsigned>(elementMapper.index(elemCtx.element()));
        if (params_.processRankOutput_ && !processRank_.empty()) {
            processRank_[elemIdx] = static_cast<unsigned>(this->simulator_.gridView().comm().rank());
        }

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& priVars = elemCtx.primaryVars(i, /*timeIdx=*/0);

            if (params_.dofIndexOutput_) {
                dofIndex_[I] = I;
            }

            for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (params_.primaryVarsOutput_ && !primaryVars_[eqIdx].empty()) {
                    primaryVars_[eqIdx][I] = priVars[eqIdx];
                }
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter)
    {
        VtkMultiWriter* vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
        if (!vtkWriter) {
            return;
        }

        if (params_.primaryVarsOutput_) {
            this->commitPriVarsBuffer_(baseWriter, "PV_%s", primaryVars_);
        }
        if (params_.processRankOutput_) {
            this->commitScalarBuffer_(baseWriter,
                                      "process rank",
                                      processRank_,
                                      /*bufferType=*/ParentType::ElementBuffer);
        }
        if (params_.dofIndexOutput_) {
            this->commitScalarBuffer_(baseWriter, "DOF index", dofIndex_);
        }
    }

private:
    VtkPrimaryVarsParams params_{};
    EqBuffer primaryVars_{};
    ScalarBuffer processRank_{};
    ScalarBuffer dofIndex_{};
};

} // namespace Opm

#endif // OPM_VTK_PRIMARY_VARS_MODULE_HPP
