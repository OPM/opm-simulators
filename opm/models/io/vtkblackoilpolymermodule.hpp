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
 * \copydoc Opm::VtkBlackOilPolymerModule
 */
#ifndef OPM_VTK_BLACK_OIL_POLYMER_MODULE_HPP
#define OPM_VTK_BLACK_OIL_POLYMER_MODULE_HPP

#include <dune/common/fvector.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkblackoilpolymerparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

namespace Opm::Properties::TTag {

// create new type tag for the VTK multi-phase output
struct VtkBlackOilPolymer {};

} // namespace Opm::Properties::TTag

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the black oil model's polymer related quantities.
 */
template <class TypeTag>
class VtkBlackOilPolymerModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>() };

    using BufferType = typename ParentType::BufferType;
    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    explicit VtkBlackOilPolymerModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        if constexpr (enablePolymer) {
            params_.read();
        }
    }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        if constexpr (enablePolymer) {
            VtkBlackoilPolymerParams::registerParameters();
        }
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if constexpr (enablePolymer) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            if (params_.polymerConcentrationOutput_) {
                this->resizeScalarBuffer_(polymerConcentration_, BufferType::Dof);
            }
            if (params_.polymerDeadPoreVolumeOutput_) {
                this->resizeScalarBuffer_(polymerDeadPoreVolume_, BufferType::Dof);
            }
            if (params_.polymerRockDensityOutput_) {
                this->resizeScalarBuffer_(polymerRockDensity_, BufferType::Dof);
            }
            if (params_.polymerAdsorptionOutput_) {
                this->resizeScalarBuffer_(polymerAdsorption_, BufferType::Dof);
            }
            if (params_.polymerViscosityCorrectionOutput_) {
                this->resizeScalarBuffer_(polymerViscosityCorrection_, BufferType::Dof);
            }
            if (params_.waterViscosityCorrectionOutput_) {
                this->resizeScalarBuffer_(waterViscosityCorrection_, BufferType::Dof);
            }
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx) override
    {
        if constexpr (enablePolymer) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
                const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                if (params_.polymerConcentrationOutput_) {
                    polymerConcentration_[globalDofIdx] =
                        scalarValue(intQuants.polymerConcentration());
                }

                if (params_.polymerDeadPoreVolumeOutput_) {
                    polymerDeadPoreVolume_[globalDofIdx] =
                        scalarValue(intQuants.polymerDeadPoreVolume());
                }

                if (params_.polymerRockDensityOutput_) {
                    polymerRockDensity_[globalDofIdx] =
                        scalarValue(intQuants.polymerRockDensity());
                }

                if (params_.polymerAdsorptionOutput_) {
                    polymerAdsorption_[globalDofIdx] =
                        scalarValue(intQuants.polymerAdsorption());
                }

                if (params_.polymerViscosityCorrectionOutput_) {
                    polymerViscosityCorrection_[globalDofIdx] =
                        scalarValue(intQuants.polymerViscosityCorrection());
                }

                if (params_.waterViscosityCorrectionOutput_) {
                    waterViscosityCorrection_[globalDofIdx] =
                        scalarValue(intQuants.waterViscosityCorrection());
                }
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter) override
    {
        if constexpr (enablePolymer) {
            if (!dynamic_cast<VtkMultiWriter*>(&baseWriter)) {
                return;
            }

            if (params_.polymerConcentrationOutput_) {
                this->commitScalarBuffer_(baseWriter, "polymer concentration",
                                          polymerConcentration_, BufferType::Dof);
            }

            if (params_.polymerDeadPoreVolumeOutput_) {
                this->commitScalarBuffer_(baseWriter, "dead pore volume fraction",
                                          polymerDeadPoreVolume_, BufferType::Dof);
            }

            if (params_.polymerRockDensityOutput_) {
                this->commitScalarBuffer_(baseWriter, "polymer rock density",
                                          polymerRockDensity_, BufferType::Dof);
            }

            if (params_.polymerAdsorptionOutput_) {
                this->commitScalarBuffer_(baseWriter, "polymer adsorption",
                                          polymerAdsorption_, BufferType::Dof);
            }

            if (params_.polymerViscosityCorrectionOutput_) {
                this->commitScalarBuffer_(baseWriter, "polymer viscosity correction",
                                          polymerViscosityCorrection_, BufferType::Dof);
            }

            if (params_.waterViscosityCorrectionOutput_) {
                this->commitScalarBuffer_(baseWriter, "water viscosity correction",
                                          waterViscosityCorrection_, BufferType::Dof);
            }
        }
    }

private:
    VtkBlackoilPolymerParams params_{};
    ScalarBuffer polymerConcentration_{};
    ScalarBuffer polymerDeadPoreVolume_{};
    ScalarBuffer polymerRockDensity_{};
    ScalarBuffer polymerAdsorption_{};
    ScalarBuffer polymerViscosityCorrection_{};
    ScalarBuffer waterViscosityCorrection_{};
};

} // namespace Opm

#endif // OPM_VTK_BLACK_OIL_POLYMER_MODULE_HPP
