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
 * \copydoc Opm::EclMaterialLawManager
 */
#if ! HAVE_ECL_INPUT
#error "Eclipse input support in opm-common is required to use the ECL material manager!"
#endif

#ifndef OPM_ECL_MATERIAL_LAW_MANAGER_HPP
#define OPM_ECL_MATERIAL_LAW_MANAGER_HPP

#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsTwoPhaseLaw.hpp>
#include <opm/material/fluidmatrixinteractions/EclHysteresisTwoPhaseLaw.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsConfig.hpp>
#include <opm/material/fluidmatrixinteractions/EclHysteresisConfig.hpp>
#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/fluidstates/SimpleModularFluidState.hpp>

#if HAVE_OPM_COMMON
#include <opm/common/OpmLog/OpmLog.hpp>
#endif

#include <opm/material/common/Exceptions.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <algorithm>


namespace Opm {

/*!
 * \ingroup fluidmatrixinteractions
 *
 * \brief Provides an simple way to create and manage the material law objects
 *        for a complete ECL deck.
 */
template <class TraitsT>
class EclMaterialLawManager
{
private:
    typedef TraitsT Traits;
    typedef typename Traits::Scalar Scalar;
    enum { waterPhaseIdx = Traits::wettingPhaseIdx };
    enum { oilPhaseIdx = Traits::nonWettingPhaseIdx };
    enum { gasPhaseIdx = Traits::gasPhaseIdx };
    enum { numPhases = Traits::numPhases };

    typedef TwoPhaseMaterialTraits<Scalar, oilPhaseIdx, gasPhaseIdx> GasOilTraits;
    typedef TwoPhaseMaterialTraits<Scalar, waterPhaseIdx, oilPhaseIdx> OilWaterTraits;

    // the two-phase material law which is defined on effective (unscaled) saturations
    typedef PiecewiseLinearTwoPhaseMaterial<GasOilTraits> GasOilEffectiveTwoPhaseLaw;
    typedef PiecewiseLinearTwoPhaseMaterial<OilWaterTraits> OilWaterEffectiveTwoPhaseLaw;
    typedef typename GasOilEffectiveTwoPhaseLaw::Params GasOilEffectiveTwoPhaseParams;
    typedef typename OilWaterEffectiveTwoPhaseLaw::Params OilWaterEffectiveTwoPhaseParams;

    // the two-phase material law which is defined on absolute (scaled) saturations
    typedef EclEpsTwoPhaseLaw<GasOilEffectiveTwoPhaseLaw> GasOilEpsTwoPhaseLaw;
    typedef EclEpsTwoPhaseLaw<OilWaterEffectiveTwoPhaseLaw> OilWaterEpsTwoPhaseLaw;
    typedef typename GasOilEpsTwoPhaseLaw::Params GasOilEpsTwoPhaseParams;
    typedef typename OilWaterEpsTwoPhaseLaw::Params OilWaterEpsTwoPhaseParams;

    // the scaled two-phase material laws with hystersis
    typedef EclHysteresisTwoPhaseLaw<GasOilEpsTwoPhaseLaw> GasOilTwoPhaseLaw;
    typedef EclHysteresisTwoPhaseLaw<OilWaterEpsTwoPhaseLaw> OilWaterTwoPhaseLaw;
    typedef typename GasOilTwoPhaseLaw::Params GasOilTwoPhaseHystParams;
    typedef typename OilWaterTwoPhaseLaw::Params OilWaterTwoPhaseHystParams;

public:
    // the three-phase material law used by the simulation
    typedef EclMultiplexerMaterial<Traits, GasOilTwoPhaseLaw, OilWaterTwoPhaseLaw> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

private:
    // internal typedefs
    typedef std::vector<std::shared_ptr<GasOilEffectiveTwoPhaseParams> > GasOilEffectiveParamVector;
    typedef std::vector<std::shared_ptr<OilWaterEffectiveTwoPhaseParams> > OilWaterEffectiveParamVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPoints<Scalar> > > GasOilScalingPointsVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPoints<Scalar> > > OilWaterScalingPointsVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPointsInfo<Scalar> > > GasOilScalingInfoVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPointsInfo<Scalar> > > OilWaterScalingInfoVector;
    typedef std::vector<std::shared_ptr<GasOilTwoPhaseHystParams> > GasOilParamVector;
    typedef std::vector<std::shared_ptr<OilWaterTwoPhaseHystParams> > OilWaterParamVector;
    typedef std::vector<std::shared_ptr<MaterialLawParams> > MaterialLawParamsVector;

public:
    EclMaterialLawManager()
    {}

    void initFromDeck(const Opm::Deck& deck,
                      const Opm::EclipseState& eclState,
                      const std::vector<int>& compressedToCartesianElemIdx)
    {
        // get the number of saturation regions and the number of cells in the deck
        const size_t numSatRegions = eclState.runspec().tabdims().getNumSatTables();
        size_t numCompressedElems = compressedToCartesianElemIdx.size();

        // copy the SATNUM grid property. in some cases this is not necessary, but it
        // should not require much memory anyway...
        satnumRegionArray_.resize(numCompressedElems);
        if (eclState.get3DProperties().hasDeckIntGridProperty("SATNUM")) {
            const auto& satnumRawData = eclState.get3DProperties().getIntGridProperty("SATNUM").getData();
            for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
                unsigned cartesianElemIdx = static_cast<unsigned>(compressedToCartesianElemIdx[elemIdx]);
                satnumRegionArray_[elemIdx] = satnumRawData[cartesianElemIdx] - 1;
            }
        }
        else
            std::fill(satnumRegionArray_.begin(), satnumRegionArray_.end(), 0);

        readGlobalEpsOptions_(deck, eclState);
        readGlobalHysteresisOptions_(deck);
        readGlobalThreePhaseOptions_(deck);

        unscaledEpsInfo_.resize(numSatRegions);
        for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx)
            unscaledEpsInfo_[satRegionIdx].extractUnscaled(deck, eclState, satRegionIdx);

        initParamsForElements_(deck, eclState, compressedToCartesianElemIdx, satnumRegionArray_);
    }

    /*!
     * \brief Modify the initial condition according to the SWATINIT keyword.
     *
     * The method returns the water saturation which yields a givenn capillary
     * pressure. The reason this method is not folded directly into initFromDeck() is
     * that the capillary pressure given depends on the particuars of how the simulator
     * calculates its initial condition.
     */
    Scalar applySwatinit(unsigned elemIdx,
                         Scalar pcow,
                         Scalar Sw)
    {
        auto& elemScaledEpsInfo = *oilWaterScaledEpsInfoDrainage_[elemIdx];

        // TODO: Mixed wettability systems - see ecl kw OPTIONS switch 74

        if (pcow < 0.0)
            Sw = elemScaledEpsInfo.Swu;
        else {

            if (Sw <= elemScaledEpsInfo.Swl)
                Sw = elemScaledEpsInfo.Swl;

            // specify a fluid state which only stores the saturations
            typedef Opm::SimpleModularFluidState<Scalar,
                                                 numPhases,
                                                 /*numComponents=*/0,
                                                 /*FluidSystem=*/void, /* -> don't care */
                                                 /*storePressure=*/false,
                                                 /*storeTemperature=*/false,
                                                 /*storeComposition=*/false,
                                                 /*storeFugacity=*/false,
                                                 /*storeSaturation=*/true,
                                                 /*storeDensity=*/false,
                                                 /*storeViscosity=*/false,
                                                 /*storeEnthalpy=*/false> FluidState;
            FluidState fs;
            fs.setSaturation(waterPhaseIdx, Sw);
            fs.setSaturation(gasPhaseIdx, 0);
            fs.setSaturation(oilPhaseIdx, 0);
            Scalar pc[numPhases] = { 0 };
            MaterialLaw::capillaryPressures(pc, materialLawParams(elemIdx), fs);

            Scalar pcowAtSw = pc[oilPhaseIdx] - pc[waterPhaseIdx];
            if (pcowAtSw > 0.0) {
                elemScaledEpsInfo.maxPcow *= pcow/pcowAtSw;
                auto& elemEclEpsScalingPoints = oilWaterScaledEpsPointsDrainage(elemIdx);
                elemEclEpsScalingPoints.init(elemScaledEpsInfo, *oilWaterEclEpsConfig_, Opm::EclOilWaterSystem);
            }
        }

        return Sw;
    }

    bool enableEndPointScaling() const
    { return enableEndPointScaling_; }

    bool enableHysteresis() const
    { return hysteresisConfig_->enableHysteresis(); }

    MaterialLawParams& materialLawParams(unsigned elemIdx)
    {
        assert(0 <= elemIdx && elemIdx <  materialLawParams_.size());
        return *materialLawParams_[elemIdx];
    }

    const MaterialLawParams& materialLawParams(unsigned elemIdx) const
    {
        assert(0 <= elemIdx && elemIdx <  materialLawParams_.size());
        return *materialLawParams_[elemIdx];
    }

    /*!
     * \brief Returns a material parameter object for a given element and saturation region.
     *
     * This method changes the saturation table idx in the original material law parameter object.
     * In the context of ECL reservoir simulators, this is required to properly handle
     * wells with its own saturation table idx. In order to reset the saturation table idx
     * in the materialLawparams_ call the method with the cells satRegionIdx
     */
    const MaterialLawParams& connectionMaterialLawParams(unsigned satRegionIdx, unsigned elemIdx) const
    {
        MaterialLawParams& mlp = *materialLawParams_[elemIdx];

#if HAVE_OPM_COMMON
        if (enableHysteresis())
            OpmLog::warning("Warning: Using non-default satnum regions for conenction is not tested in combination with hysteresis");
#endif
        // Currently we don't support COMPIMP. I.e. use the same table lookup for the hysteresis curves.
        // unsigned impRegionIdx = satRegionIdx;

        // change the sat table it points to.
        switch (mlp.approach()) {
        case EclStone1Approach: {
            auto& realParams = mlp.template getRealParams<Opm::EclStone1Approach>();

            realParams.oilWaterParams().drainageParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[satRegionIdx]);
            realParams.oilWaterParams().drainageParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setUnscaledPoints(gasOilUnscaledPointsVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setEffectiveLawParams(gasOilEffectiveParamVector_[satRegionIdx]);
//            if (enableHysteresis()) {
//                realParams.oilWaterParams().imbibitionParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[impRegionIdx]);
//                realParams.oilWaterParams().imbibitionParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setUnscaledPoints(gasOilUnscaledPointsVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setEffectiveLawParams(gasOilEffectiveParamVector_[impRegionIdx]);
//            }
        }
            break;

        case EclStone2Approach: {
            auto& realParams = mlp.template getRealParams<Opm::EclStone2Approach>();
            realParams.oilWaterParams().drainageParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[satRegionIdx]);
            realParams.oilWaterParams().drainageParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setUnscaledPoints(gasOilUnscaledPointsVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setEffectiveLawParams(gasOilEffectiveParamVector_[satRegionIdx]);
//            if (enableHysteresis()) {
//                realParams.oilWaterParams().imbibitionParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[impRegionIdx]);
//                realParams.oilWaterParams().imbibitionParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setUnscaledPoints(gasOilUnscaledPointsVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setEffectiveLawParams(gasOilEffectiveParamVector_[impRegionIdx]);
//            }
        }
            break;

        case EclDefaultApproach: {
            auto& realParams = mlp.template getRealParams<Opm::EclDefaultApproach>();
            realParams.oilWaterParams().drainageParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[satRegionIdx]);
            realParams.oilWaterParams().drainageParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setUnscaledPoints(gasOilUnscaledPointsVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setEffectiveLawParams(gasOilEffectiveParamVector_[satRegionIdx]);
//            if (enableHysteresis()) {
//                realParams.oilWaterParams().imbibitionParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[impRegionIdx]);
//                realParams.oilWaterParams().imbibitionParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setUnscaledPoints(gasOilUnscaledPointsVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setEffectiveLawParams(gasOilEffectiveParamVector_[impRegionIdx]);
//            }
        }
            break;

        case EclTwoPhaseApproach: {
            auto& realParams = mlp.template getRealParams<Opm::EclTwoPhaseApproach>();
            realParams.oilWaterParams().drainageParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[satRegionIdx]);
            realParams.oilWaterParams().drainageParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setUnscaledPoints(gasOilUnscaledPointsVector_[satRegionIdx]);
            realParams.gasOilParams().drainageParams().setEffectiveLawParams(gasOilEffectiveParamVector_[satRegionIdx]);
//            if (enableHysteresis()) {
//                realParams.oilWaterParams().imbibitionParams().setUnscaledPoints(oilWaterUnscaledPointsVector_[impRegionIdx]);
//                realParams.oilWaterParams().imbibitionParams().setEffectiveLawParams(oilWaterEffectiveParamVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setUnscaledPoints(gasOilUnscaledPointsVector_[impRegionIdx]);
//                realParams.gasOilParams().imbibitionParams().setEffectiveLawParams(gasOilEffectiveParamVector_[impRegionIdx]);
//            }
        }
            break;

        default:
            throw std::logic_error("Enum value for material approach unknown!");
        }

        return mlp;
    }

    int satnumRegionIdx(unsigned elemIdx) const {
        return satnumRegionArray_[elemIdx];
    }

    std::shared_ptr<MaterialLawParams>& materialLawParamsPointerReferenceHack(unsigned elemIdx)
    {
        assert(0 <= elemIdx && elemIdx <  materialLawParams_.size());
        return materialLawParams_[elemIdx];
    }

    template <class FluidState>
    void updateHysteresis(const FluidState& fluidState, unsigned elemIdx)
    {
        if (!enableHysteresis())
            return;

        auto threePhaseParams = materialLawParams_[elemIdx];
        MaterialLaw::updateHysteresis(*threePhaseParams, fluidState);
    }

    void oilWaterHysteresisParams(Scalar& pcSwMdc,
                                  Scalar& krnSwMdc,
                                  unsigned elemIdx) const
    {
        if (!enableHysteresis()) {
            throw std::runtime_error("Cannot get hysteresis parameters if hysteresis not enabled.");
        }
        const auto& params = materialLawParams(elemIdx);
        MaterialLaw::oilWaterHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    void setOilWaterHysteresisParams(const Scalar& pcSwMdc,
                                     const Scalar& krnSwMdc,
                                     unsigned elemIdx)
    {
        if (!enableHysteresis()) {
            throw std::runtime_error("Cannot set hysteresis parameters if hysteresis not enabled.");
        }
        auto& params = materialLawParams(elemIdx);
        MaterialLaw::setOilWaterHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    void gasOilHysteresisParams(Scalar& pcSwMdc,
                                Scalar& krnSwMdc,
                                unsigned elemIdx) const
    {
        if (!enableHysteresis()) {
            throw std::runtime_error("Cannot get hysteresis parameters if hysteresis not enabled.");
        }
        const auto& params = materialLawParams(elemIdx);
        MaterialLaw::gasOilHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    void setGasOilHysteresisParams(const Scalar& pcSwMdc,
                                   const Scalar& krnSwMdc,
                                   unsigned elemIdx)
    {
        if (!enableHysteresis()) {
            throw std::runtime_error("Cannot set hysteresis parameters if hysteresis not enabled.");
        }
        auto& params = materialLawParams(elemIdx);
        MaterialLaw::setGasOilHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    EclEpsScalingPoints<Scalar>& oilWaterScaledEpsPointsDrainage(unsigned elemIdx)
    {
        auto& materialParams = *materialLawParams_[elemIdx];
        switch (materialParams.approach()) {
        case EclStone1Approach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclStone1Approach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclStone2Approach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclStone2Approach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclDefaultApproach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclDefaultApproach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclTwoPhaseApproach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclTwoPhaseApproach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }
        default:
            throw std::logic_error("Enum value for material approach unknown!");
        }
    }

    const Opm::EclEpsScalingPointsInfo<Scalar>& oilWaterScaledEpsInfoDrainage(size_t elemIdx) const
    {
        return *oilWaterScaledEpsInfoDrainage_[elemIdx];
    }

    std::shared_ptr<EclEpsScalingPointsInfo<Scalar> >& oilWaterScaledEpsInfoDrainagePointerReferenceHack(unsigned elemIdx)
    {
        return oilWaterScaledEpsInfoDrainage_[elemIdx];
    }
private:
    void readGlobalEpsOptions_(const Opm::Deck& deck, const Opm::EclipseState& eclState)
    {
        oilWaterEclEpsConfig_ = std::make_shared<Opm::EclEpsConfig>();
        oilWaterEclEpsConfig_-> initFromDeck(deck, eclState, Opm::EclOilWaterSystem);

        enableEndPointScaling_ = deck.hasKeyword("ENDSCALE");
    }

    void readGlobalHysteresisOptions_(const Opm::Deck& deck)
    {
        hysteresisConfig_ = std::make_shared<Opm::EclHysteresisConfig>();
        hysteresisConfig_->initFromDeck(deck);
    }

    void readGlobalThreePhaseOptions_(const Opm::Deck& deck)
    {
        bool gasEnabled = deck.hasKeyword("GAS");
        bool oilEnabled = deck.hasKeyword("OIL");
        bool waterEnabled = deck.hasKeyword("WATER");

        int numEnabled =
            (gasEnabled?1:0)
            + (oilEnabled?1:0)
            + (waterEnabled?1:0);

        if (numEnabled < 2)
            throw std::runtime_error("At least two fluid phases must be enabled. (Is: "+std::to_string(numEnabled)+")");

        if (numEnabled == 2) {
            threePhaseApproach_ = Opm::EclTwoPhaseApproach;
            if (!gasEnabled)
                twoPhaseApproach_ = Opm::EclTwoPhaseOilWater;
            else if (!oilEnabled)
                twoPhaseApproach_ = Opm::EclTwoPhaseGasWater;
            else if (!waterEnabled)
                twoPhaseApproach_ = Opm::EclTwoPhaseGasOil;
        }
        else {
            assert(numEnabled == 3);

            threePhaseApproach_ = Opm::EclDefaultApproach;
            if (deck.hasKeyword("STONE") || deck.hasKeyword("STONE2"))
                threePhaseApproach_ = Opm::EclStone2Approach;
            else if (deck.hasKeyword("STONE1"))
                threePhaseApproach_ = Opm::EclStone1Approach;
        }
    }

    void initParamsForElements_(const Deck& deck, const EclipseState& eclState,
                                const std::vector<int>& compressedToCartesianElemIdx,
                                const std::vector<int>& satnumRegionArray)
    {
        const size_t numSatRegions = eclState.runspec().tabdims().getNumSatTables();
        unsigned numCompressedElems = static_cast<unsigned>(compressedToCartesianElemIdx.size());

        // read the end point scaling configuration. this needs to be done only once per
        // deck.
        auto gasOilConfig = std::make_shared<Opm::EclEpsConfig>();
        auto oilWaterConfig = std::make_shared<Opm::EclEpsConfig>();
        gasOilConfig->initFromDeck(deck, eclState, Opm::EclGasOilSystem);
        oilWaterConfig->initFromDeck(deck, eclState, Opm::EclOilWaterSystem);

        // read the saturation region specific parameters from the deck
        gasOilUnscaledPointsVector_.resize(numSatRegions);
        oilWaterUnscaledPointsVector_.resize(numSatRegions);
        gasOilEffectiveParamVector_.resize(numSatRegions);
        oilWaterEffectiveParamVector_.resize(numSatRegions);
        for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
            // unscaled points for end-point scaling
            readGasOilUnscaledPoints_(gasOilUnscaledPointsVector_, gasOilConfig, deck, eclState, satRegionIdx);
            readOilWaterUnscaledPoints_(oilWaterUnscaledPointsVector_, oilWaterConfig, deck, eclState, satRegionIdx);

            // the parameters for the effective two-phase matererial laws
            readGasOilEffectiveParameters_(gasOilEffectiveParamVector_, deck, eclState, satRegionIdx);
            readOilWaterEffectiveParameters_(oilWaterEffectiveParamVector_, deck, eclState, satRegionIdx);

            // read the end point scaling info for the saturation region
            unscaledEpsInfo_[satRegionIdx].extractUnscaled(deck, eclState, satRegionIdx);

        }

        // read the scaled end point scaling parameters which are specific for each
        // element
        GasOilScalingInfoVector gasOilScaledInfoVector(numCompressedElems);
        oilWaterScaledEpsInfoDrainage_.resize(numCompressedElems);
        GasOilScalingInfoVector gasOilScaledImbInfoVector;
        OilWaterScalingInfoVector oilWaterScaledImbInfoVector;

        GasOilScalingPointsVector gasOilScaledPointsVector(numCompressedElems);
        GasOilScalingPointsVector oilWaterScaledEpsPointsDrainage(numCompressedElems);
        GasOilScalingPointsVector gasOilScaledImbPointsVector;
        OilWaterScalingPointsVector oilWaterScaledImbPointsVector;

        if (enableHysteresis()) {
            gasOilScaledImbInfoVector.resize(numCompressedElems);
            gasOilScaledImbPointsVector.resize(numCompressedElems);
            oilWaterScaledImbInfoVector.resize(numCompressedElems);
            oilWaterScaledImbPointsVector.resize(numCompressedElems);
        }

        EclEpsGridProperties epsGridProperties, epsImbGridProperties;
        epsGridProperties.initFromDeck(deck, eclState, /*imbibition=*/false);
        if (enableHysteresis())
            epsImbGridProperties.initFromDeck(deck, eclState, /*imbibition=*/true);
        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            unsigned cartElemIdx = static_cast<unsigned>(compressedToCartesianElemIdx[elemIdx]);
            readGasOilScaledPoints_(gasOilScaledInfoVector,
                                    gasOilScaledPointsVector,
                                    gasOilConfig,
                                    eclState,
                                    epsGridProperties,
                                    elemIdx,
                                    cartElemIdx);
            readOilWaterScaledPoints_(oilWaterScaledEpsInfoDrainage_,
                                      oilWaterScaledEpsPointsDrainage,
                                      oilWaterConfig,
                                      eclState,
                                      epsGridProperties,
                                      elemIdx,
                                      cartElemIdx);

            if (enableHysteresis()) {
                readGasOilScaledPoints_(gasOilScaledImbInfoVector,
                                        gasOilScaledImbPointsVector,
                                        gasOilConfig,
                                        eclState,
                                        epsImbGridProperties,
                                        elemIdx,
                                        cartElemIdx);
                readOilWaterScaledPoints_(oilWaterScaledImbInfoVector,
                                          oilWaterScaledImbPointsVector,
                                          oilWaterConfig,
                                          eclState,
                                          epsImbGridProperties,
                                          elemIdx,
                                          cartElemIdx);
            }
        }

        // create the parameter objects for the two-phase laws
        GasOilParamVector gasOilParams(numCompressedElems);
        OilWaterParamVector oilWaterParams(numCompressedElems);
        GasOilParamVector gasOilImbParams;
        OilWaterParamVector oilWaterImbParams;

        if (enableHysteresis()) {
            gasOilImbParams.resize(numCompressedElems);
            oilWaterImbParams.resize(numCompressedElems);
        }

        bool hasGas = deck.hasKeyword("GAS");
        bool hasOil = deck.hasKeyword("OIL");
        bool hasWater = deck.hasKeyword("WATER");

        const auto& imbnumData = eclState.get3DProperties().getIntGridProperty("IMBNUM").getData();
        assert(numCompressedElems == satnumRegionArray.size());
        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            unsigned satRegionIdx = static_cast<unsigned>(satnumRegionArray[elemIdx]);

            gasOilParams[elemIdx] = std::make_shared<GasOilTwoPhaseHystParams>();
            oilWaterParams[elemIdx] = std::make_shared<OilWaterTwoPhaseHystParams>();

            gasOilParams[elemIdx]->setConfig(hysteresisConfig_);
            oilWaterParams[elemIdx]->setConfig(hysteresisConfig_);

            if (hasGas && hasOil) {
                auto gasOilDrainParams = std::make_shared<GasOilEpsTwoPhaseParams>();
                gasOilDrainParams->setConfig(gasOilConfig);
                gasOilDrainParams->setUnscaledPoints(gasOilUnscaledPointsVector_[satRegionIdx]);
                gasOilDrainParams->setScaledPoints(gasOilScaledPointsVector[elemIdx]);
                gasOilDrainParams->setEffectiveLawParams(gasOilEffectiveParamVector_[satRegionIdx]);
                gasOilDrainParams->finalize();

                gasOilParams[elemIdx]->setDrainageParams(gasOilDrainParams,
                                                         *gasOilScaledInfoVector[elemIdx],
                                                         EclGasOilSystem);
            }

            if (hasOil && hasWater) {
                auto oilWaterDrainParams = std::make_shared<OilWaterEpsTwoPhaseParams>();
                oilWaterDrainParams->setConfig(oilWaterConfig);
                oilWaterDrainParams->setUnscaledPoints(oilWaterUnscaledPointsVector_[satRegionIdx]);
                oilWaterDrainParams->setScaledPoints(oilWaterScaledEpsPointsDrainage[elemIdx]);
                oilWaterDrainParams->setEffectiveLawParams(oilWaterEffectiveParamVector_[satRegionIdx]);
                oilWaterDrainParams->finalize();

                oilWaterParams[elemIdx]->setDrainageParams(oilWaterDrainParams,
                                                           *oilWaterScaledEpsInfoDrainage_[elemIdx],
                                                           EclOilWaterSystem);
            }

            if (enableHysteresis()) {
                unsigned imbRegionIdx = static_cast<unsigned>(imbnumData[elemIdx]) - 1;

                if (hasGas && hasOil) {
                    auto gasOilImbParamsHyst = std::make_shared<GasOilEpsTwoPhaseParams>();
                    gasOilImbParamsHyst->setConfig(gasOilConfig);
                    gasOilImbParamsHyst->setUnscaledPoints(gasOilUnscaledPointsVector_[imbRegionIdx]);
                    gasOilImbParamsHyst->setScaledPoints(gasOilScaledImbPointsVector[elemIdx]);
                    gasOilImbParamsHyst->setEffectiveLawParams(gasOilEffectiveParamVector_[imbRegionIdx]);
                    gasOilImbParamsHyst->finalize();

                    gasOilParams[elemIdx]->setImbibitionParams(gasOilImbParamsHyst,
                                                               *gasOilScaledImbInfoVector[elemIdx],
                                                               EclGasOilSystem);
                }

                if (hasOil && hasWater) {
                    auto oilWaterImbParamsHyst = std::make_shared<OilWaterEpsTwoPhaseParams>();
                    oilWaterImbParamsHyst->setConfig(oilWaterConfig);
                    oilWaterImbParamsHyst->setUnscaledPoints(oilWaterUnscaledPointsVector_[imbRegionIdx]);
                    oilWaterImbParamsHyst->setScaledPoints(oilWaterScaledImbPointsVector[elemIdx]);
                    oilWaterImbParamsHyst->setEffectiveLawParams(oilWaterEffectiveParamVector_[imbRegionIdx]);
                    oilWaterImbParamsHyst->finalize();

                    oilWaterParams[elemIdx]->setImbibitionParams(oilWaterImbParamsHyst,
                                                                 *gasOilScaledImbInfoVector[elemIdx],
                                                                 EclGasOilSystem);
                }
            }

            if (hasGas && hasOil)
                gasOilParams[elemIdx]->finalize();

            if (hasOil && hasWater)
                oilWaterParams[elemIdx]->finalize();
        }

        // create the parameter objects for the three-phase law
        materialLawParams_.resize(numCompressedElems);
        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            materialLawParams_[elemIdx] = std::make_shared<MaterialLawParams>();
            unsigned satRegionIdx = static_cast<unsigned>(satnumRegionArray[elemIdx]);

            initThreePhaseParams_(deck,
                                  eclState,
                                  *materialLawParams_[elemIdx],
                                  satRegionIdx,
                                  *oilWaterScaledEpsInfoDrainage_[elemIdx],
                                  oilWaterParams[elemIdx],
                                  gasOilParams[elemIdx]);

            materialLawParams_[elemIdx]->finalize();
        }
    }

    // The saturation function family.
    // If SWOF and SGOF are specified in the deck it return FamilyI
    // If SWFN, SGFN and SOF3 are specified in the deck it return FamilyII
    // If keywords are missing or mixed, an error is given.
    enum SaturationFunctionFamily {
        noFamily,
        FamilyI,
        FamilyII
    };

    SaturationFunctionFamily getSaturationFunctionFamily(const Opm::Deck& deck, const Opm::EclipseState& eclState) const
    {
        const auto& tableManager = eclState.getTableManager();
        const TableContainer& swofTables = tableManager.getSwofTables();
        const TableContainer& slgofTables= tableManager.getSlgofTables();
        const TableContainer& sgofTables = tableManager.getSgofTables();
        const TableContainer& swfnTables = tableManager.getSwfnTables();
        const TableContainer& sgfnTables = tableManager.getSgfnTables();
        const TableContainer& sof3Tables = tableManager.getSof3Tables();
        const TableContainer& sof2Tables = tableManager.getSof2Tables();


        bool hasGas = deck.hasKeyword("GAS");
        bool hasOil = deck.hasKeyword("OIL");
        bool hasWater = deck.hasKeyword("WATER");

        bool family1 = false;
        bool family2 = false;
        if (!hasGas) {
            // oil-water case
            family1 = !swofTables.empty();
            family2 = !swfnTables.empty() && !sof2Tables.empty();
        }
        else if (!hasWater) {
            // oil-gas case
            family1 = !sgofTables.empty();
            family2 = !sgfnTables.empty() && !sof2Tables.empty();
        }
        else if (!hasOil) {
            // water-gas case
            throw std::runtime_error("water-gas two-phase simulations are currently not supported");
        }
        else {
            // three-phase case
            family1 = (!sgofTables.empty() || !slgofTables.empty()) && !swofTables.empty();
            family2 = !swfnTables.empty() && !sgfnTables.empty() && !sof3Tables.empty();
        }

        if (family1 && family2) {
            throw std::invalid_argument("Saturation families should not be mixed \n"
                                        "Use either SGOF and SWOF or SGFN, SWFN and SOF3");
        }

        if (!family1 && !family2) {
            throw std::invalid_argument("Saturations function must be specified using either "
                                        "family 1 or family 2 keywords \n"
                                        "Use either SGOF and SWOF or SGFN, SWFN and SOF3" );
        }

        if (family1 && !family2)
            return SaturationFunctionFamily::FamilyI;
        else if (family2 && !family1)
            return SaturationFunctionFamily::FamilyII;
        return SaturationFunctionFamily::noFamily; // no family or two families
    }

    template <class Container>
    void readGasOilEffectiveParameters_(Container& dest,
                                        const Opm::Deck& deck,
                                        const Opm::EclipseState& eclState,
                                        unsigned satRegionIdx)
    {
        bool hasGas = deck.hasKeyword("GAS");
        bool hasOil = deck.hasKeyword("OIL");

        if (!hasGas || !hasOil)
            // we don't read anything if either the gas or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<GasOilEffectiveTwoPhaseParams>();

        auto& effParams = *dest[satRegionIdx];

        // the situation for the gas phase is complicated that all saturations are
        // shifted by the connate water saturation.
        Scalar Swco = unscaledEpsInfo_[satRegionIdx].Swl;
        const auto& tableManager = eclState.getTableManager();

        switch (getSaturationFunctionFamily(deck, eclState)) {
        case FamilyI:
        {
            const TableContainer& sgofTables = tableManager.getSgofTables();
            const TableContainer& slgofTables = tableManager.getSlgofTables();
            if (!sgofTables.empty())
                readGasOilEffectiveParametersSgof_(effParams,
                                                   Swco,
                                                   sgofTables.getTable<SgofTable>(satRegionIdx));
            else if (!slgofTables.empty())
                readGasOilEffectiveParametersSlgof_(effParams,
                                                    Swco,
                                                    slgofTables.getTable<SlgofTable>(satRegionIdx));
            break;
        }

        case FamilyII:
        {
            const SgfnTable& sgfnTable = tableManager.getSgfnTables().getTable<SgfnTable>( satRegionIdx );
            bool hasWater = deck.hasKeyword("WATER");
            if (!hasWater) {
                // oil and gas case
                const Sof2Table& sof2Table = tableManager.getSof2Tables().getTable<Sof2Table>( satRegionIdx );
                readGasOilEffectiveParametersFamily2_(effParams,
                                                      Swco,
                                                      sof2Table,
                                                      sgfnTable);
            } else {
                const Sof3Table& sof3Table = tableManager.getSof3Tables().getTable<Sof3Table>( satRegionIdx );
                readGasOilEffectiveParametersFamily2_(effParams,
                                                      Swco,
                                                      sof3Table,
                                                      sgfnTable);
            }
            break;
        }

        //default:
        case noFamily:
            throw std::domain_error("No valid saturation keyword family specified");

        }
    }

    void readGasOilEffectiveParametersSgof_(GasOilEffectiveTwoPhaseParams& effParams,
                                            Scalar Swco,
                                            const Opm::SgofTable& sgofTable)
    {
        // convert the saturations of the SGOF keyword from gas to oil saturations
        std::vector<double> SoSamples(sgofTable.numRows());
        std::vector<double> SoKroSamples(sgofTable.numRows());
        for (size_t sampleIdx = 0; sampleIdx < sgofTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = 1 - sgofTable.get("SG", sampleIdx);
            SoKroSamples[sampleIdx] = SoSamples[sampleIdx] - Swco;
        }

        effParams.setKrwSamples(SoKroSamples, sgofTable.getColumn("KROG").vectorCopy());
        effParams.setKrnSamples(SoSamples, sgofTable.getColumn("KRG").vectorCopy());
        effParams.setPcnwSamples(SoSamples, sgofTable.getColumn("PCOG").vectorCopy());
        effParams.finalize();
    }

    void readGasOilEffectiveParametersSlgof_(GasOilEffectiveTwoPhaseParams& effParams,
                                             Scalar Swco,
                                             const Opm::SlgofTable& slgofTable)
    {
        // convert the saturations of the SLGOF keyword from "liquid" to oil saturations
        std::vector<double> SoSamples(slgofTable.numRows());
        std::vector<double> SoKroSamples(slgofTable.numRows());
        for (size_t sampleIdx = 0; sampleIdx < slgofTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = slgofTable.get("SL", sampleIdx);
            SoKroSamples[sampleIdx] = slgofTable.get("SL", sampleIdx) - Swco;
        }

        effParams.setKrwSamples(SoKroSamples, slgofTable.getColumn("KROG").vectorCopy());
        effParams.setKrnSamples(SoSamples, slgofTable.getColumn("KRG").vectorCopy());
        effParams.setPcnwSamples(SoSamples, slgofTable.getColumn("PCOG").vectorCopy());
        effParams.finalize();
    }

    void readGasOilEffectiveParametersFamily2_(GasOilEffectiveTwoPhaseParams& effParams,
                                               Scalar /* Swco */,
                                               const Opm::Sof3Table& sof3Table,
                                               const Opm::SgfnTable& sgfnTable)
    {
        // convert the saturations of the SGFN keyword from gas to oil saturations
        std::vector<double> SoSamples(sgfnTable.numRows());
        std::vector<double> SoColumn = sof3Table.getColumn("SO").vectorCopy();
        for (size_t sampleIdx = 0; sampleIdx < sgfnTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = 1 - sgfnTable.get("SG", sampleIdx);
        }

        effParams.setKrwSamples(SoColumn, sof3Table.getColumn("KROG").vectorCopy());
        effParams.setKrnSamples(SoSamples, sgfnTable.getColumn("KRG").vectorCopy());
        effParams.setPcnwSamples(SoSamples, sgfnTable.getColumn("PCOG").vectorCopy());
        effParams.finalize();
    }

    void readGasOilEffectiveParametersFamily2_(GasOilEffectiveTwoPhaseParams& effParams,
                                               Scalar /* Swco */,
                                               const Opm::Sof2Table& sof2Table,
                                               const Opm::SgfnTable& sgfnTable)
    {
        // convert the saturations of the SGFN keyword from gas to oil saturations
        std::vector<double> SoSamples(sgfnTable.numRows());
        std::vector<double> SoColumn = sof2Table.getColumn("SO").vectorCopy();
        for (size_t sampleIdx = 0; sampleIdx < sgfnTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = 1 - sgfnTable.get("SG", sampleIdx);
        }

        effParams.setKrwSamples(SoColumn, sof2Table.getColumn("KRO").vectorCopy());
        effParams.setKrnSamples(SoSamples, sgfnTable.getColumn("KRG").vectorCopy());
        effParams.setPcnwSamples(SoSamples, sgfnTable.getColumn("PCOG").vectorCopy());
        effParams.finalize();
    }

    template <class Container>
    void readOilWaterEffectiveParameters_(Container& dest,
                                          const Opm::Deck& deck,
                                          const Opm::EclipseState& eclState,
                                          unsigned satRegionIdx)
    {
        bool hasWater = deck.hasKeyword("WATER");
        bool hasOil = deck.hasKeyword("OIL");

        if (!hasOil || !hasWater)
            // we don't read anything if either the water or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<OilWaterEffectiveTwoPhaseParams>();

        const auto& tableManager = eclState.getTableManager();
        auto& effParams = *dest[satRegionIdx];

        switch (getSaturationFunctionFamily(deck, eclState)) {
        case FamilyI: {
            const auto& swofTable = tableManager.getSwofTables().getTable<SwofTable>(satRegionIdx);
            std::vector<double> SwColumn = swofTable.getColumn("SW").vectorCopy();

            effParams.setKrwSamples(SwColumn, swofTable.getColumn("KRW").vectorCopy());
            effParams.setKrnSamples(SwColumn, swofTable.getColumn("KROW").vectorCopy());
            effParams.setPcnwSamples(SwColumn, swofTable.getColumn("PCOW").vectorCopy());
            effParams.finalize();
            break;
        }
        case FamilyII:
        {
            const auto& swfnTable = tableManager.getSwfnTables().getTable<SwfnTable>(satRegionIdx);
            const auto& sof3Table = tableManager.getSof3Tables().getTable<Sof3Table>(satRegionIdx);
            std::vector<double> SwColumn = swfnTable.getColumn("SW").vectorCopy();

            // convert the saturations of the SOF3 keyword from oil to water saturations
            std::vector<double> SwSamples(sof3Table.numRows());
            for (size_t sampleIdx = 0; sampleIdx < sof3Table.numRows(); ++ sampleIdx)
                SwSamples[sampleIdx] = 1 - sof3Table.get("SO", sampleIdx);

            effParams.setKrwSamples(SwColumn, swfnTable.getColumn("KRW").vectorCopy());
            effParams.setKrnSamples(SwSamples, sof3Table.getColumn("KROW").vectorCopy());
            effParams.setPcnwSamples(SwColumn, swfnTable.getColumn("PCOW").vectorCopy());
            effParams.finalize();
            break;
        }

        case noFamily:
        //default:
            throw std::domain_error("No valid saturation keyword family specified");

        }
    }


    template <class Container>
    void readGasOilUnscaledPoints_(Container& dest,
                                   std::shared_ptr<EclEpsConfig> config,
                                   const Opm::Deck& deck,
                                   const Opm::EclipseState& /* eclState */,
                                   unsigned satRegionIdx)
    {
        bool hasGas = deck.hasKeyword("GAS");
        bool hasOil = deck.hasKeyword("OIL");

        if (!hasGas || !hasOil)
            // we don't read anything if either the gas or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satRegionIdx]->init(unscaledEpsInfo_[satRegionIdx], *config, EclGasOilSystem);
    }

    template <class Container>
    void readOilWaterUnscaledPoints_(Container& dest,
                                     std::shared_ptr<EclEpsConfig> config,
                                     const Opm::Deck& deck,
                                     const Opm::EclipseState& /* eclState */,
                                     unsigned satRegionIdx)
    {
        bool hasWater = deck.hasKeyword("WATER");
        bool hasOil = deck.hasKeyword("OIL");

        if (!hasOil || !hasWater)
            // we don't read anything if either the water or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satRegionIdx]->init(unscaledEpsInfo_[satRegionIdx], *config, EclOilWaterSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readGasOilScaledPoints_(InfoContainer& destInfo,
                                 PointsContainer& destPoints,
                                 std::shared_ptr<EclEpsConfig> config,
                                 const Opm::EclipseState& eclState,
                                 const EclEpsGridProperties& epsGridProperties,
                                 unsigned elemIdx,
                                 unsigned cartElemIdx)
    {
        unsigned satRegionIdx = static_cast<unsigned>((*epsGridProperties.satnum)[cartElemIdx]) - 1; // ECL uses Fortran indices!

        destInfo[elemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satRegionIdx]);
        destInfo[elemIdx]->extractScaled(eclState, epsGridProperties, cartElemIdx);

        destPoints[elemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[elemIdx]->init(*destInfo[elemIdx], *config, EclGasOilSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readOilWaterScaledPoints_(InfoContainer& destInfo,
                                   PointsContainer& destPoints,
                                   std::shared_ptr<EclEpsConfig> config,
                                   const Opm::EclipseState& eclState,
                                   const EclEpsGridProperties& epsGridProperties,
                                   unsigned elemIdx,
                                   unsigned cartElemIdx)
    {
        unsigned satRegionIdx = static_cast<unsigned>((*epsGridProperties.satnum)[cartElemIdx]) - 1; // ECL uses Fortran indices!

        destInfo[elemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satRegionIdx]);
        destInfo[elemIdx]->extractScaled(eclState, epsGridProperties, cartElemIdx);

        destPoints[elemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[elemIdx]->init(*destInfo[elemIdx], *config, EclOilWaterSystem);
    }

    void initThreePhaseParams_(const Opm::Deck& deck,
                               const Opm::EclipseState& /* eclState */,
                               MaterialLawParams& materialParams,
                               unsigned satRegionIdx,
                               const EclEpsScalingPointsInfo<Scalar>& epsInfo,
                               std::shared_ptr<OilWaterTwoPhaseHystParams> oilWaterParams,
                               std::shared_ptr<GasOilTwoPhaseHystParams> gasOilParams)
    {
        materialParams.setApproach(threePhaseApproach_);

        switch (materialParams.approach()) {
        case EclStone1Approach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclStone1Approach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);

            if (deck.hasKeyword("STONE1EX")) {
                Scalar eta =
                    deck.getKeyword("STONE1EX").getRecord(satRegionIdx).getItem(0).getSIDouble(0);
                realParams.setEta(eta);
            }
            else
                realParams.setEta(1.0);
            realParams.finalize();
            break;
        }

        case EclStone2Approach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclStone2Approach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);
            realParams.finalize();
            break;
        }

        case EclDefaultApproach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclDefaultApproach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);
            realParams.finalize();
            break;
        }

        case EclTwoPhaseApproach: {
            auto& realParams = materialParams.template getRealParams<Opm::EclTwoPhaseApproach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setApproach(twoPhaseApproach_);
            realParams.finalize();
            break;
        }
        }
    }

    bool enableEndPointScaling_;
    std::shared_ptr<EclHysteresisConfig> hysteresisConfig_;

    std::shared_ptr<EclEpsConfig> oilWaterEclEpsConfig_;
    std::vector<Opm::EclEpsScalingPointsInfo<Scalar>> unscaledEpsInfo_;
    OilWaterScalingInfoVector oilWaterScaledEpsInfoDrainage_;

    GasOilScalingPointsVector gasOilUnscaledPointsVector_;
    OilWaterScalingPointsVector oilWaterUnscaledPointsVector_;
    GasOilEffectiveParamVector gasOilEffectiveParamVector_;
    OilWaterEffectiveParamVector oilWaterEffectiveParamVector_;

    Opm::EclMultiplexerApproach threePhaseApproach_;
    // this attribute only makes sense for twophase simulations!
    enum EclTwoPhaseApproach twoPhaseApproach_;

    std::vector<std::shared_ptr<MaterialLawParams> > materialLawParams_;

    std::vector<int> satnumRegionArray_;
};
} // namespace Opm

#endif
