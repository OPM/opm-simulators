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
 * \copydoc Opm::VtkBlackOilSolventModule
 */
#ifndef OPM_VTK_BLACK_OIL_SOLVENT_MODULE_HPP
#define OPM_VTK_BLACK_OIL_SOLVENT_MODULE_HPP

#include <dune/common/fvector.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkblackoilsolventparams.hpp>
#include <opm/models/io/vtkmultiwriter.hh>

#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the black oil model's solvent related quantities.
 */
template <class TypeTag>
class VtkBlackOilSolventModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };

    using BufferType = typename ParentType::BufferType;
    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    explicit VtkBlackOilSolventModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        if constexpr (enableSolvent) {
            params_.read();
        }
    }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        if constexpr (enableSolvent) {
            VtkBlackOilSolventParams::registerParameters();
        }
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers() override
    {
        if constexpr (enableSolvent) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            if (params_.solventSaturationOutput_) {
                this->resizeScalarBuffer_(solventSaturation_, BufferType::Dof);
            }
            if (params_.solventRswOutput_) {
                this->resizeScalarBuffer_(solventRsw_, BufferType::Dof);
            }
            if (params_.solventDensityOutput_) {
                this->resizeScalarBuffer_(solventDensity_, BufferType::Dof);
            }
            if (params_.solventViscosityOutput_) {
                this->resizeScalarBuffer_(solventViscosity_, BufferType::Dof);
            }
            if (params_.solventMobilityOutput_) {
                this->resizeScalarBuffer_(solventMobility_, BufferType::Dof);
            }
        }
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx) override
    {
        if constexpr (enableSolvent) {
            if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
                return;
            }

            using Toolbox = MathToolbox<Evaluation>;
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
                const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                if (params_.solventSaturationOutput_) {
                    solventSaturation_[globalDofIdx] =
                        Toolbox::scalarValue(intQuants.solventSaturation());
                }

                if (params_.solventRswOutput_) {
                    solventRsw_[globalDofIdx] =
                        Toolbox::scalarValue(intQuants.rsSolw());
                }

                if (params_.solventDensityOutput_) {
                    solventDensity_[globalDofIdx] =
                        Toolbox::scalarValue(intQuants.solventDensity());
                }

                if (params_.solventViscosityOutput_) {
                    solventViscosity_[globalDofIdx] =
                        Toolbox::scalarValue(intQuants.solventViscosity());
                }

                if (params_.solventMobilityOutput_) {
                    solventMobility_[globalDofIdx] =
                        Toolbox::scalarValue(intQuants.solventMobility());
                }
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter) override
    {
        if constexpr (enableSolvent) {
            if (!dynamic_cast<VtkMultiWriter*>(&baseWriter)) {
                return;
            }

            if (params_.solventSaturationOutput_) {
                this->commitScalarBuffer_(baseWriter, "saturation_solvent",
                                          solventSaturation_, BufferType::Dof);
            }

            if (params_.solventRswOutput_) {
                this->commitScalarBuffer_(baseWriter, "dissolved_solvent",
                                          solventRsw_, BufferType::Dof);
            }

            if (params_.solventDensityOutput_) {
                this->commitScalarBuffer_(baseWriter, "density_solvent",
                                          solventDensity_, BufferType::Dof);
            }

            if (params_.solventViscosityOutput_) {
                this->commitScalarBuffer_(baseWriter, "viscosity_solvent",
                                          solventViscosity_, BufferType::Dof);
            }

            if (params_.solventMobilityOutput_) {
                this->commitScalarBuffer_(baseWriter, "mobility_solvent",
                                          solventMobility_, BufferType::Dof);
            }
        }
    }

private:
    VtkBlackOilSolventParams params_{};
    ScalarBuffer solventSaturation_{};
    ScalarBuffer solventRsw_{};
    ScalarBuffer solventDensity_{};
    ScalarBuffer solventViscosity_{};
    ScalarBuffer solventMobility_{};
};

} // namespace Opm

#endif // OPM_VTK_BLACK_OIL_SOLVENT_MODULE_HPP
