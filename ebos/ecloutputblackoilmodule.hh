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
 * \copydoc Ewoms::EclOutputBlackOilModule
 */
#ifndef EWOMS_ECL_OUTPUT_BLACK_OIL_MODULE_HH
#define EWOMS_ECL_OUTPUT_BLACK_OIL_MODULE_HH

#include "eclwriter.hh"
#include "ecldeckunits.hh"

#include <ewoms/io/baseoutputmodule.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <opm/common/Valgrind.hpp>

#include <dune/common/fvector.hh>

#include <type_traits>

namespace Ewoms {
namespace Properties {
// create new type tag for the VTK multi-phase output
NEW_TYPE_TAG(EclOutputBlackOil);

// create the property tags needed for the multi phase module
NEW_PROP_TAG(EclOutputWriteSaturations);
NEW_PROP_TAG(EclOutputWritePressures);
NEW_PROP_TAG(EclOutputWriteGasDissolutionFactor);
NEW_PROP_TAG(EclOutputWriteOilVaporizationFactor);
NEW_PROP_TAG(EclOutputWriteGasFormationVolumeFactor);
NEW_PROP_TAG(EclOutputWriteOilFormationVolumeFactor);
NEW_PROP_TAG(EclOutputWriteOilSaturationPressure);

// set default values for what quantities to output
SET_BOOL_PROP(EclOutputBlackOil, EclOutputWriteSaturations, true);
SET_BOOL_PROP(EclOutputBlackOil, EclOutputWritePressures, true);
SET_BOOL_PROP(EclOutputBlackOil, EclOutputWriteGasDissolutionFactor, true);
SET_BOOL_PROP(EclOutputBlackOil, EclOutputWriteOilVaporizationFactor, false);
SET_BOOL_PROP(EclOutputBlackOil, EclOutputWriteGasFormationVolumeFactor, true);
SET_BOOL_PROP(EclOutputBlackOil, EclOutputWriteOilFormationVolumeFactor, true);
SET_BOOL_PROP(EclOutputBlackOil, EclOutputWriteOilSaturationPressure, false);
} // namespace Properties

// forward declaration
template <class TypeTag>
class EcfvDiscretization;

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Output module for the results black oil model writing in
 *        ECL binary format.
 */
template <class TypeTag>
class EclOutputBlackOilModule : public BaseOutputModule<TypeTag>
{
    typedef BaseOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Discretization) Discretization;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    typedef Ewoms::EclWriter<TypeTag> EclWriter;

    enum { numPhases = FluidSystem::numPhases };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;

public:
    EclOutputBlackOilModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output
     * module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclOutputWriteSaturations,
                             "Include the saturations of all fluid phases in the "
                             "ECL output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclOutputWritePressures,
                             "Include the absolute pressures of all fluid phases in the "
                             "ECL output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclOutputWriteGasDissolutionFactor,
                             "Include the gas dissolution factor of saturated oil in the "
                             "ECL output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclOutputWriteOilVaporizationFactor,
                             "Include the oil vaporization factor of saturated gas in the "
                             "ECL output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclOutputWriteGasFormationVolumeFactor,
                             "Include the gas formation volume factor in the "
                             "ECL output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclOutputWriteOilFormationVolumeFactor,
                             "Include the oil formation volume factor of saturated oil "
                             "in the ECL output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclOutputWriteOilSaturationPressure,
                             "Include the saturation pressure of oil in the "
                             "ECL output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (!std::is_same<Discretization, Ewoms::EcfvDiscretization<TypeTag> >::value)
            return;

        auto bufferType = ParentType::ElementBuffer;
        if (saturationsOutput_()) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
                this->resizeScalarBuffer_(saturation_[phaseIdx], bufferType);
        }
        if (pressuresOutput_()) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
                this->resizeScalarBuffer_(pressure_[phaseIdx], bufferType);
        }
        if (gasDissolutionFactorOutput_())
            this->resizeScalarBuffer_(gasDissolutionFactor_, bufferType);
        if (oilVaporizationFactorOutput_())
            this->resizeScalarBuffer_(oilVaporizationFactor_, bufferType);
        if (gasFormationVolumeFactorOutput_())
            this->resizeScalarBuffer_(gasFormationVolumeFactor_, bufferType);
        if (oilSaturationPressureOutput_())
            this->resizeScalarBuffer_(oilSaturationPressure_, bufferType);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quanties relevant
     *        for an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableEclOutput))
            return;

        typedef Opm::MathToolbox<Evaluation> Toolbox;

        if (!std::is_same<Discretization, Ewoms::EcfvDiscretization<TypeTag> >::value)
            return;

        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& fs = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0).fluidState();
            typedef typename std::remove_const<typename std::remove_reference<decltype(fs)>::type>::type FluidState;
            unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
            unsigned pvtRegionIdx = elemCtx.primaryVars(dofIdx, /*timeIdx=*/0).pvtRegionIndex();

            if (saturationsOutput_()) {
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx))
                        continue;

                    saturation_[phaseIdx][globalDofIdx] = Toolbox::value(fs.saturation(phaseIdx));
                    Opm::Valgrind::CheckDefined(saturation_[phaseIdx][globalDofIdx]);
                }
            }
            if (pressuresOutput_()) {
                for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx))
                        continue;

                    pressure_[phaseIdx][globalDofIdx] = Toolbox::value(fs.pressure(phaseIdx));
                    Opm::Valgrind::CheckDefined(pressure_[phaseIdx][globalDofIdx]);
                }
            }
            if (gasDissolutionFactorOutput_()) {
                Scalar SoMax = elemCtx.model().maxOilSaturation(globalDofIdx);
                gasDissolutionFactor_[globalDofIdx] =
                    FluidSystem::template saturatedDissolutionFactor<FluidState, Scalar>(fs, gasPhaseIdx, pvtRegionIdx, SoMax);
                Opm::Valgrind::CheckDefined(gasDissolutionFactor_[globalDofIdx]);
            }
            if (oilVaporizationFactorOutput_()) {
                Scalar SoMax = elemCtx.model().maxOilSaturation(globalDofIdx);
                gasDissolutionFactor_[globalDofIdx] =
                    FluidSystem::template saturatedDissolutionFactor<FluidState, Scalar>(fs, oilPhaseIdx, pvtRegionIdx, SoMax);
                Opm::Valgrind::CheckDefined(gasDissolutionFactor_[globalDofIdx]);
            }
            if (gasFormationVolumeFactorOutput_()) {
                gasFormationVolumeFactor_[globalDofIdx] =
                    1.0/FluidSystem::template inverseFormationVolumeFactor<FluidState, Scalar>(fs, gasPhaseIdx, pvtRegionIdx);
                Opm::Valgrind::CheckDefined(gasFormationVolumeFactor_[globalDofIdx]);
            }
            if (oilSaturationPressureOutput_()) {
                oilSaturationPressure_[globalDofIdx] =
                    FluidSystem::template saturationPressure<FluidState, Scalar>(fs, oilPhaseIdx, pvtRegionIdx);
                Opm::Valgrind::CheckDefined(oilSaturationPressure_[globalDofIdx]);
            }
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& writer)
    {
        if (!std::is_same<Discretization, Ewoms::EcfvDiscretization<TypeTag> >::value)
            return;

        if (!dynamic_cast<EclWriter*>(&writer))
            return; // this module only consideres ecl writers...

        typedef EclDeckUnits<TypeTag> DeckUnits;
        const DeckUnits& deckUnits = this->simulator_.problem().deckUnits();

        typename ParentType::BufferType bufferType = ParentType::ElementBuffer;
        if (pressuresOutput_()) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
                deckUnits.siToDeck(pressure_[phaseIdx], DeckUnits::pressure);

            this->commitScalarBuffer_(writer, "PRESSURE", pressure_[oilPhaseIdx], bufferType);
            this->commitScalarBuffer_(writer, "PGAS", pressure_[gasPhaseIdx], bufferType);
            this->commitScalarBuffer_(writer, "PWAT", pressure_[waterPhaseIdx], bufferType);
        }
        if (saturationsOutput_()) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
                deckUnits.siToDeck(saturation_[phaseIdx], DeckUnits::saturation);

            this->commitScalarBuffer_(writer, "SWAT", saturation_[waterPhaseIdx], bufferType);
            this->commitScalarBuffer_(writer, "SGAS", saturation_[gasPhaseIdx], bufferType);
            // the oil saturation is _NOT_ written to disk. Instead, it is calculated by
            // the visualization tool. Wondering why is probably a waste of time...
        }
        if (gasDissolutionFactorOutput_()) {
            deckUnits.siToDeck(gasDissolutionFactor_, DeckUnits::gasDissolutionFactor);
            this->commitScalarBuffer_(writer, "RS", gasDissolutionFactor_, bufferType);
        }
        if (oilVaporizationFactorOutput_()) {
            deckUnits.siToDeck(oilVaporizationFactor_, DeckUnits::oilVaporizationFactor);
            this->commitScalarBuffer_(writer, "RV", oilVaporizationFactor_, bufferType);
        }
        if (gasFormationVolumeFactorOutput_()) {
            // no unit conversion required
            this->commitScalarBuffer_(writer, "BG", gasFormationVolumeFactor_, bufferType);
        }
        if (oilSaturationPressureOutput_()) {
            deckUnits.siToDeck(oilSaturationPressure_, DeckUnits::pressure);
            this->commitScalarBuffer_(writer, "PSAT", oilSaturationPressure_);
        }
    }

private:
    static bool saturationsOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclOutputWriteSaturations); }

    static bool pressuresOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EclOutputWritePressures); }

    static bool gasDissolutionFactorOutput_()
    {
        return
            FluidSystem::phaseIsActive(oilPhaseIdx) &&
            FluidSystem::phaseIsActive(gasPhaseIdx) &&
            EWOMS_GET_PARAM(TypeTag, bool, EclOutputWriteGasDissolutionFactor);
    }

    static bool gasFormationVolumeFactorOutput_()
    {
        return
            FluidSystem::phaseIsActive(oilPhaseIdx) &&
            FluidSystem::phaseIsActive(gasPhaseIdx) &&
            EWOMS_GET_PARAM(TypeTag, bool, EclOutputWriteGasFormationVolumeFactor);
    }

    static bool oilVaporizationFactorOutput_()
    {
        return
            FluidSystem::phaseIsActive(oilPhaseIdx) &&
            FluidSystem::phaseIsActive(gasPhaseIdx) &&
            EWOMS_GET_PARAM(TypeTag, bool, EclOutputWriteOilVaporizationFactor);
    }

    static bool oilSaturationPressureOutput_()
    {
        return
            FluidSystem::phaseIsActive(oilPhaseIdx) &&
            FluidSystem::phaseIsActive(gasPhaseIdx) &&
            EWOMS_GET_PARAM(TypeTag, bool, EclOutputWriteOilSaturationPressure);
    }

    ScalarBuffer saturation_[numPhases];
    ScalarBuffer pressure_[numPhases];
    ScalarBuffer gasDissolutionFactor_;
    ScalarBuffer oilVaporizationFactor_;
    ScalarBuffer gasFormationVolumeFactor_;
    ScalarBuffer oilSaturationPressure_;
};
} // namespace Ewoms

#endif
