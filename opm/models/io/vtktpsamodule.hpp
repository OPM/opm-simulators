// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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
#ifndef VTK_TPSA_MODULE_HPP
#define VTK_TPSA_MODULE_HPP

#include <opm/material/densead/Math.hpp>

#include <opm/models/io/baseoutputmodule.hh>
#include <opm/models/io/vtkmultiwriter.hh>
#include <opm/models/io/vtktpsaparams.hpp>
#include <opm/models/tpsa/tpsabaseproperties.hpp>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>


namespace Opm {

/*!
* \ingroup Vtk
*
* \brief VTK output module for TPSA quantities
*/
template <class TypeTag>
class VtkTpsaModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    static constexpr auto vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { enableMech = getPropValue<TypeTag, Properties::EnableMech>() };

    using BufferType = typename ParentType::BufferType;
    using ScalarBuffer = typename ParentType::ScalarBuffer;
    using VectorBuffer = typename ParentType::VectorBuffer;

public:
    /*!
    * \brief Constructor
    *
    * \param simulator Simulator object
    */
    explicit VtkTpsaModule(const Simulator& simulator)
        : ParentType(simulator)
    {
        // Read runtime parameters
        if constexpr(enableMech) {
            params_.read();
        }
    }

    /*!
    * \brief Register runtime parameters
    */
    static void registerParameters()
    {
        if constexpr(enableMech) {
            VtkTpsaParams::registerParameters();
        }
    }

    /*!
    * \brief Allocate memory for output
    */
    void allocBuffers() override
    {
        // Check if vtk output is enabled
        if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
            return;
        }

        // Allocate buffers
        if constexpr(enableMech) {
            // Displacement
            if (params_.displacementOutput_) {
                this->resizeVectorBuffer_(displacement_, BufferType::Dof);
            }

            // Rotation
            if (params_.rotationOutput_) {
                this->resizeVectorBuffer_(rotation_, BufferType::Dof);
            }

            // Solid pressure
            if (params_.solidPressureOutput_) {
                this->resizeScalarBuffer_(solidPres_, BufferType::Dof);
            }
        }
    }

    /*!
    * \brief Assign quantities to output buffers
    *
    * \param elemCtx Element context object
    */
    void processElement(const ElementContext& elemCtx) override
    {
        // Check if vtk output is enabled
        if (!Parameters::Get<Parameters::EnableVtkOutput>()) {
            return;
        }

        if constexpr(enableMech) {
            // Assign quantities from material state
            const auto& problem = elemCtx.problem();
            const auto& geoMechModel = problem.geoMechModel();
            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
                const auto& materialState = geoMechModel.materialState(globalDofIdx, /*timeIdx=*/0);

                // Loop over x-, y- and z-dir (corresponding to dirIdx = 0, 1, 2)
                for (int dirIdx = 0; dirIdx < 3; ++dirIdx) {
                    // Displacement
                    if (params_.displacementOutput_) {
                        displacement_[globalDofIdx][dirIdx] = scalarValue(materialState.displacement(dirIdx));
                    }

                    // Rotation
                    if (params_.rotationOutput_) {
                        rotation_[globalDofIdx][dirIdx] = scalarValue(materialState.rotation(dirIdx));
                    }
                }

                // Solid pressure
                if (params_.solidPressureOutput_) {
                    solidPres_[globalDofIdx] = scalarValue(materialState.solidPressure());
                }
            }
        }
    }

    /*!
    * \brief Add buffers to VTK writer
    *
    * \param baseWriter VTK output writer object
    */
    void commitBuffers(BaseOutputWriter& baseWriter) override
    {
        // Check if writer exists or mechanics output is enabled
        if (!dynamic_cast<VtkMultiWriter*>(&baseWriter)) {
            return;
        }

        if constexpr(enableMech) {
            // Displacement
            if (params_.displacementOutput_) {
                this->commitVectorBuffer_(baseWriter, "displacement", displacement_, BufferType::Dof);
            }

            // Rotation
            if (params_.rotationOutput_) {
                this->commitVectorBuffer_(baseWriter, "rotation", rotation_, BufferType::Dof);
            }

            // Solid pressure
            if (params_.solidPressureOutput_) {
                this->commitScalarBuffer_(baseWriter, "solid_pressure", solidPres_, BufferType::Dof);
            }
        }
    }

private:
    VtkTpsaParams params_{};

    VectorBuffer displacement_{};
    VectorBuffer rotation_{};
    ScalarBuffer solidPres_{};
};  // class VtkTpsaModule

}  // namespace Opm

#endif