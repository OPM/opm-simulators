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
#ifndef EWOMS_VTK_PRIMARY_VARS_MODULE_HH
#define EWOMS_VTK_PRIMARY_VARS_MODULE_HH

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>

BEGIN_PROPERTIES

// create new type tag for the VTK primary variables output
NEW_PROP_TAG(EnableVtkOutput);
NEW_TYPE_TAG(VtkPrimaryVars);

// create the property tags needed for the primary variables module
NEW_PROP_TAG(VtkWritePrimaryVars);
NEW_PROP_TAG(VtkWriteProcessRank);
NEW_PROP_TAG(VtkWriteDofIndex);
NEW_PROP_TAG(VtkOutputFormat);
NEW_PROP_TAG(EnableVtkOutput);

SET_BOOL_PROP(VtkPrimaryVars, VtkWritePrimaryVars, false);
SET_BOOL_PROP(VtkPrimaryVars, VtkWriteProcessRank, false);
SET_BOOL_PROP(VtkPrimaryVars, VtkWriteDofIndex, false);

END_PROPERTIES

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the fluid composition
 */
template<class TypeTag>
class VtkPrimaryVarsModule : public BaseOutputModule<TypeTag>
{
    typedef BaseOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
    typedef Opm::VtkMultiWriter<GridView, vtkFormat> VtkMultiWriter;

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::EqBuffer EqBuffer;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

public:
    VtkPrimaryVarsModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePrimaryVars,
                             "Include the primary variables into the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteProcessRank,
                             "Include the MPI process rank into the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteDofIndex,
                             "Include the index of the degrees of freedom into the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (primaryVarsOutput_())
            this->resizeEqBuffer_(primaryVars_);
        if (processRankOutput_())
            this->resizeScalarBuffer_(processRank_,
                                      /*bufferType=*/ParentType::ElementBuffer);
        if (dofIndexOutput_())
            this->resizeScalarBuffer_(dofIndex_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        const auto& elementMapper = elemCtx.model().elementMapper();
        unsigned elemIdx = static_cast<unsigned>(elementMapper.index(elemCtx.element()));
        if (processRankOutput_() && !processRank_.empty())
            processRank_[elemIdx] = static_cast<unsigned>(this->simulator_.gridView().comm().rank());

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& priVars = elemCtx.primaryVars(i, /*timeIdx=*/0);

            if (dofIndexOutput_())
                dofIndex_[I] = I;

            for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx) {
                if (primaryVarsOutput_() && !primaryVars_[eqIdx].empty())
                    primaryVars_[eqIdx][I] = priVars[eqIdx];
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

        if (primaryVarsOutput_())
            this->commitPriVarsBuffer_(baseWriter, "PV_%s", primaryVars_);
        if (processRankOutput_())
            this->commitScalarBuffer_(baseWriter,
                                      "process rank",
                                      processRank_,
                                      /*bufferType=*/ParentType::ElementBuffer);
        if (dofIndexOutput_())
            this->commitScalarBuffer_(baseWriter, "DOF index", dofIndex_);
    }

private:
    static bool primaryVarsOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWritePrimaryVars);
        return val;
    }
    static bool processRankOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteProcessRank);
        return val;
    }
    static bool dofIndexOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteDofIndex);
        return val;
    }

    EqBuffer primaryVars_;
    ScalarBuffer processRank_;
    ScalarBuffer dofIndex_;
};

} // namespace Opm

#endif
