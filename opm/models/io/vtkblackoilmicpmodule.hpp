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
 * \copydoc Opm::VtkBlackOilMICPModule
 */
#ifndef OPM_VTK_BLACK_OIL_MICP_MODULE_HPP
#define OPM_VTK_BLACK_OIL_MICP_MODULE_HPP

#include <dune/common/fvector.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkblackoilmicpparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the MICP model's related quantities.
 */
template <class TypeTag>
class VtkBlackOilMICPModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { enableMICP = getPropValue<TypeTag, Properties::EnableMICP>() };

    using BufferType = typename ParentType::BufferType;
    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    explicit VtkBlackOilMICPModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        if constexpr (enableMICP) {
            params_.read();
        }
    }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        if constexpr (enableMICP) {
            VtkBlackoilMICPParams::registerParameters();
        }
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if constexpr (enableMICP) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            if (params_.microbialConcentrationOutput_) {
                this->resizeScalarBuffer_(microbialConcentration_, BufferType::Dof);
            }
            if (params_.oxygenConcentrationOutput_) {
                this->resizeScalarBuffer_(oxygenConcentration_, BufferType::Dof);
            }
            if (params_.ureaConcentrationOutput_) {
                this->resizeScalarBuffer_(ureaConcentration_, BufferType::Dof);
            }
            if (params_.biofilmConcentrationOutput_) {
                this->resizeScalarBuffer_(biofilmConcentration_, BufferType::Dof);
            }
            if (params_.calciteConcentrationOutput_) {
                this->resizeScalarBuffer_(calciteConcentration_, BufferType::Dof);
            }
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx) override
    {
        if constexpr (enableMICP) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
                const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                if (params_.microbialConcentrationOutput_) {
                    microbialConcentration_[globalDofIdx] =
                        scalarValue(intQuants.microbialConcentration());
                }

                if (params_.oxygenConcentrationOutput_) {
                    oxygenConcentration_[globalDofIdx] =
                        scalarValue(intQuants.oxygenConcentration());
                }

                if (params_.ureaConcentrationOutput_) {
                    ureaConcentration_[globalDofIdx] =
                        scalarValue(intQuants.ureaConcentration());
                }

                if (params_.biofilmConcentrationOutput_) {
                    biofilmConcentration_[globalDofIdx] =
                        scalarValue(intQuants.biofilmConcentration());
                }

                if (params_.calciteConcentrationOutput_) {
                    calciteConcentration_[globalDofIdx] =
                        scalarValue(intQuants.calciteConcentration());
                }
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter) override
    {
        if constexpr (enableMICP) {
            if (!dynamic_cast<VtkMultiWriter*>(&baseWriter)) {
                return;
            }

            if (params_.microbialConcentrationOutput_) {
                this->commitScalarBuffer_(baseWriter, "microbial concentration",
                                          microbialConcentration_, BufferType::Dof);
            }

            if (params_.oxygenConcentrationOutput_) {
                this->commitScalarBuffer_(baseWriter, "oxygen concentration",
                                          oxygenConcentration_, BufferType::Dof);
            }

            if (params_.ureaConcentrationOutput_) {
                this->commitScalarBuffer_(baseWriter, "urea concentration",
                                          ureaConcentration_, BufferType::Dof);
            }

            if (params_.biofilmConcentrationOutput_) {
                this->commitScalarBuffer_(baseWriter, "biofilm fraction",
                                          biofilmConcentration_, BufferType::Dof);
            }

            if (params_.calciteConcentrationOutput_) {
                this->commitScalarBuffer_(baseWriter, "calcite fraction",
                                          calciteConcentration_, BufferType::Dof);
            }
        }
    }

private:
    VtkBlackoilMICPParams params_{};
    ScalarBuffer microbialConcentration_{};
    ScalarBuffer oxygenConcentration_{};
    ScalarBuffer ureaConcentration_{};
    ScalarBuffer biofilmConcentration_{};
    ScalarBuffer calciteConcentration_{};
};

} // namespace Opm

#endif // OPM_VTK_BLACKOIL_MICP_MODULE_HPP
