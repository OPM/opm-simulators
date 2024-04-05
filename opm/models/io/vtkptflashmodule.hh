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
 * \copydoc Opm::VtkPTFlashModule
 */
#ifndef OPM_VTK_PTFLASH_MODULE_HH
#define OPM_VTK_PTFLASH_MODULE_HH

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/material/common/MathToolbox.hpp>

namespace Opm::Properties {

namespace TTag {

// create new type tag for the VTK PTFlash output
struct VtkPTFlash {};

} // namespace TTag

// create the property tags needed for the composition module
template<class TypeTag, class MyTypeTag>
struct VtkWriteLiquidMoleFractions { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteEquilibriumConstants { using type = UndefinedProperty; };

// set default values for what quantities to output
template<class TypeTag>
struct VtkWriteLiquidMoleFractions<TypeTag, TTag::VtkPTFlash> { static constexpr bool value = false; };
template<class TypeTag>
struct VtkWriteEquilibriumConstants<TypeTag, TTag::VtkPTFlash> { static constexpr bool value = false; };

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the PT Flash calculation
 * This module deals with the following quantities:
 * K, equilibrium ratio for all the components
 * L, liquid fraction in the two-phase system
 */
template <class TypeTag>
class VtkPTFlashModule: public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

    static const int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    using ComponentBuffer = typename ParentType::ComponentBuffer;
    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    explicit VtkPTFlashModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters()
    {
        Parameters::registerParam<TypeTag, Properties::VtkWriteLiquidMoleFractions>
            ("Include liquid mole fractions (L) in the VTK output files");
        Parameters::registerParam<TypeTag, Properties::VtkWriteEquilibriumConstants>
            ("Include equilibrium constants (K) in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (LOutput_())
            this->resizeScalarBuffer_(L_);
        if (equilConstOutput_())
            this->resizeComponentBuffer_(K_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        using Toolbox = MathToolbox<Evaluation>;

        if (!Parameters::get<TypeTag, Properties::EnableVtkOutput>())
            return;

        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            if (LOutput_())
                L_[I] = Toolbox::value(fs.L());

            for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
                if (equilConstOutput_())
                    K_[compIdx][I] = Toolbox::value(fs.K(compIdx));
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter)
    {
        auto *vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
        if (!vtkWriter) {
            return;
        }

        if (equilConstOutput_())
            this->commitComponentBuffer_(baseWriter, "K^%s", K_);
        if (LOutput_())
            this->commitScalarBuffer_(baseWriter, "L", L_);
    }

private:
    static bool LOutput_()
    {
        static bool val = Parameters::get<TypeTag, Properties::VtkWriteLiquidMoleFractions>();
        return val;
    }

    static bool equilConstOutput_()
    {
        static bool val = Parameters::get<TypeTag, Properties::VtkWriteEquilibriumConstants>();
        return val;
    }

    ComponentBuffer K_;
    ScalarBuffer L_;
};

} // namespace Opm

#endif
