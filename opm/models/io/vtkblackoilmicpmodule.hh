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
#ifndef EWOMS_VTK_BLACK_OIL_MICP_MODULE_HH
#define EWOMS_VTK_BLACK_OIL_MICP_MODULE_HH

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <dune/common/fvector.hh>

#include <opm/material/densead/Math.hpp>

#include <opm/models/discretization/common/fvbaseparameters.hh>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/blackoil/blackoilproperties.hh>

namespace Opm::Properties {

namespace TTag {

// create new type tag for the VTK multi-phase output
struct VtkBlackOilMICP {};

} // namespace TTag

// create the property tags needed for the MICP output module
template<class TypeTag, class MyTypeTag>
struct VtkWriteMicrobialConcentration { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteOxygenConcentration { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteUreaConcentration { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteBiofilmConcentration { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct VtkWriteCalciteConcentration { using type = UndefinedProperty; };

// set default values for what quantities to output
template<class TypeTag>
struct VtkWriteMicrobialConcentration<TypeTag, TTag::VtkBlackOilMICP> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteOxygenConcentration<TypeTag, TTag::VtkBlackOilMICP> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteUreaConcentration<TypeTag, TTag::VtkBlackOilMICP> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteBiofilmConcentration<TypeTag, TTag::VtkBlackOilMICP> { static constexpr bool value = true; };
template<class TypeTag>
struct VtkWriteCalciteConcentration<TypeTag, TTag::VtkBlackOilMICP> { static constexpr bool value = true; };

} // namespace Opm::Properties

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

    static const int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { enableMICP = getPropValue<TypeTag, Properties::EnableMICP>() };

    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    VtkBlackOilMICPModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        if (!enableMICP)
            return;

        Parameters::registerParam<TypeTag, Properties::VtkWriteMicrobialConcentration>
            ("Include the concentration of the microbial component in the water phase "
             "in the VTK output files");
        Parameters::registerParam<TypeTag, Properties::VtkWriteOxygenConcentration>
            ("Include the concentration of the oxygen component in the water phase "
             "in the VTK output files");
        Parameters::registerParam<TypeTag, Properties::VtkWriteUreaConcentration>
            ("Include the concentration of the urea component in the water phase "
             "in the VTK output files");
        Parameters::registerParam<TypeTag, Properties::VtkWriteBiofilmConcentration>
            ("Include the biofilm volume fraction in the VTK output files");
        Parameters::registerParam<TypeTag, Properties::VtkWriteCalciteConcentration>
            ("Include the calcite volume fraction in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (!Parameters::get<TypeTag, Parameters::EnableVtkOutput>())
            return;

        if (!enableMICP)
            return;

        if (microbialConcentrationOutput_())
            this->resizeScalarBuffer_(microbialConcentration_);
        if (oxygenConcentrationOutput_())
            this->resizeScalarBuffer_(oxygenConcentration_);
        if (ureaConcentrationOutput_())
            this->resizeScalarBuffer_(ureaConcentration_);
        if (biofilmConcentrationOutput_())
            this->resizeScalarBuffer_(biofilmConcentration_);
        if (calciteConcentrationOutput_())
            this->resizeScalarBuffer_(calciteConcentration_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!Parameters::get<TypeTag, Parameters::EnableVtkOutput>())
            return;

        if (!enableMICP)
            return;

        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            if (microbialConcentrationOutput_())
                microbialConcentration_[globalDofIdx] =
                    scalarValue(intQuants.microbialConcentration());

            if (oxygenConcentrationOutput_())
                oxygenConcentration_[globalDofIdx] =
                    scalarValue(intQuants.oxygenConcentration());

            if (ureaConcentrationOutput_())
                ureaConcentration_[globalDofIdx] =
                    10 * scalarValue(intQuants.ureaConcentration());//Multypliging by scaling factor 10 (see WellInterface_impl.hpp)

            if (biofilmConcentrationOutput_())
                biofilmConcentration_[globalDofIdx] =
                    scalarValue(intQuants.biofilmConcentration());

            if (calciteConcentrationOutput_())
                calciteConcentration_[globalDofIdx] =
                    scalarValue(intQuants.calciteConcentration());

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

        if (!enableMICP)
            return;

        if (microbialConcentrationOutput_())
            this->commitScalarBuffer_(baseWriter, "microbial concentration", microbialConcentration_);

        if (oxygenConcentrationOutput_())
            this->commitScalarBuffer_(baseWriter, "oxygen concentration", oxygenConcentration_);

        if (ureaConcentrationOutput_())
            this->commitScalarBuffer_(baseWriter, "urea concentration", ureaConcentration_);

        if (biofilmConcentrationOutput_())
            this->commitScalarBuffer_(baseWriter, "biofilm fraction", biofilmConcentration_);

        if (calciteConcentrationOutput_())
            this->commitScalarBuffer_(baseWriter, "calcite fraction", calciteConcentration_);

    }

private:
    static bool microbialConcentrationOutput_()
    {
        static bool val = Parameters::get<TypeTag, Properties::VtkWriteMicrobialConcentration>();
        return val;
    }

    static bool oxygenConcentrationOutput_()
    {
        static bool val = Parameters::get<TypeTag, Properties::VtkWriteOxygenConcentration>();
        return val;
    }

    static bool ureaConcentrationOutput_()
    {
        static bool val = Parameters::get<TypeTag, Properties::VtkWriteUreaConcentration>();
        return val;
    }

    static bool biofilmConcentrationOutput_()
    {
        static bool val = Parameters::get<TypeTag, Properties::VtkWriteBiofilmConcentration>();
        return val;
    }

    static bool calciteConcentrationOutput_()
    {
        static bool val = Parameters::get<TypeTag, Properties::VtkWriteCalciteConcentration>();
        return val;
    }

    ScalarBuffer microbialConcentration_;
    ScalarBuffer oxygenConcentration_;
    ScalarBuffer ureaConcentration_;
    ScalarBuffer biofilmConcentration_;
    ScalarBuffer calciteConcentration_;
};
} // namespace Opm

#endif
