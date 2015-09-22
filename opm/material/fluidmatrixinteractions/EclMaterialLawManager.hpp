// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2015 by Andreas Lauser
  Copyright (C) 2015 by IRIS AS

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
*/
/*!
 * \file
 * \copydoc Opm::EclMaterialLawManager
 */
#if ! HAVE_OPM_PARSER
#error "The opm-parser module is required to use the ECL material manager!"
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

#include <opm/material/common/Exceptions.hpp>
#include <opm/material/common/ErrorMacros.hpp>

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
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

    void initFromDeck(Opm::DeckConstPtr deck,
                      Opm::EclipseStateConstPtr eclState,
                      const std::vector<int>& compressedToCartesianElemIdx)
    {
        compressedToCartesianElemIdx_ = compressedToCartesianElemIdx;
        // get the number of saturation regions and the number of cells in the deck
        unsigned numSatRegions = static_cast<unsigned>(deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0));
        size_t numCompressedElems = compressedToCartesianElemIdx.size();

        // copy the SATNUM grid property. in some cases this is not necessary, but it
        // should not require much memory anyway...
        satnumRegionIdx_.resize(numCompressedElems);
        if (eclState->hasIntGridProperty("SATNUM")) {
            const auto& satnumRawData = eclState->getIntGridProperty("SATNUM")->getData();
            for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
                unsigned cartesianElemIdx = static_cast<unsigned>(compressedToCartesianElemIdx_[elemIdx]);
                satnumRegionIdx_[elemIdx] = satnumRawData[cartesianElemIdx] - 1;
            }
        }
        else
            std::fill(satnumRegionIdx_.begin(), satnumRegionIdx_.end(), 0);

        readGlobalEpsOptions_(deck, eclState);
        readGlobalHysteresisOptions_(deck);
        readGlobalThreePhaseOptions_(deck);

        unscaledEpsInfo_.resize(numSatRegions);
        for (unsigned satnumRegionIdx = 0; satnumRegionIdx < numSatRegions; ++satnumRegionIdx)
            unscaledEpsInfo_[satnumRegionIdx].extractUnscaled(deck, eclState, satnumRegionIdx);

        if (!hasElementSpecificParameters())
            initNonElemSpecific_(deck, eclState);
        else
            initElemSpecific_(deck, eclState);
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
        if (Sw <= elemScaledEpsInfo.Swl)
            Sw = elemScaledEpsInfo.Swl;
        else if (pcow < 0.0)
            Sw = elemScaledEpsInfo.Swu;
        else {
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
            Scalar pc[numPhases];
            MaterialLaw::capillaryPressures(pc, materialLawParams(elemIdx), fs);

            Scalar pcowAtSw = pc[oilPhaseIdx] - pc[waterPhaseIdx];
            if (pcowAtSw > 0.0) {
                elemScaledEpsInfo.maxPcow *= pcow/pcowAtSw;
                auto& elemEclEpsScalingPoints = getOilWaterScaledEpsPointsDrainage_(elemIdx);
                elemEclEpsScalingPoints.init(elemScaledEpsInfo, *oilWaterEclEpsConfig_, Opm::EclOilWaterSystem);
            }
        }

        return Sw;
    }

    bool enableEndPointScaling() const
    { return enableEndPointScaling_; }

    bool enableHysteresis() const
    { return hysteresisConfig_->enableHysteresis(); }

    bool hasElementSpecificParameters() const
    { return enableEndPointScaling() || enableHysteresis(); }

    MaterialLawParams& materialLawParams(unsigned elemIdx)
    {
        if (hasElementSpecificParameters()) {
            assert(0 <= elemIdx && elemIdx < (int) materialLawParams_.size());
            return *materialLawParams_[elemIdx];
        }
        else
            return *materialLawParams_[satnumRegionIdx_[elemIdx]];
    }

    const MaterialLawParams& materialLawParams(unsigned elemIdx) const
    {
        if (hasElementSpecificParameters()) {
            assert(0 <= elemIdx && elemIdx < (int) materialLawParams_.size());
            return *materialLawParams_[elemIdx];
        }
        else
            return *materialLawParams_[satnumRegionIdx_[elemIdx]];
    }

    template <class FluidState>
    void updateHysteresis(const FluidState& fluidState, unsigned elemIdx)
    {
        if (!enableHysteresis())
            return;

        auto threePhaseParams = materialLawParams_[elemIdx];
        MaterialLaw::updateHysteresis(*threePhaseParams, fluidState);
    }

    const Opm::EclEpsScalingPointsInfo<Scalar>& oilWaterScaledEpsInfoDrainage(unsigned elemIdx) const
    {
        if (hasElementSpecificParameters())
            return *oilWaterScaledEpsInfoDrainage_[elemIdx];

        return unscaledEpsInfo_[satnumRegionIdx_[elemIdx]];
    }

private:
    void readGlobalEpsOptions_(Opm::DeckConstPtr deck, Opm::EclipseStateConstPtr eclState)
    {
        oilWaterEclEpsConfig_ = std::make_shared<Opm::EclEpsConfig>();
        oilWaterEclEpsConfig_-> initFromDeck(deck, eclState, Opm::EclOilWaterSystem);

        enableEndPointScaling_ = deck->hasKeyword("ENDSCALE");

        if (enableEndPointScaling()) {
            // sift through the options of the ENDSCALE keyword
            Opm::DeckKeywordConstPtr endscaleKeyword = deck->getKeyword("ENDSCALE");
            Opm::DeckRecordConstPtr endscaleRecord = endscaleKeyword->getRecord(0);
            for (unsigned itemIdx = 0; itemIdx < endscaleRecord->size() && itemIdx < 2; ++ itemIdx) {
                std::string optionValue = endscaleRecord->getItem(itemIdx)->getTrimmedString(0);

                // convert the value of the option to upper case, just to be sure
                std::transform(optionValue.begin(),
                               optionValue.end(),
                               optionValue.begin(),
                               ::toupper);

                if (optionValue == "DIRECT") {
                    OPM_THROW(std::runtime_error,
                              "Directional end-point scaling (indicated by the 'DIRECT' option"
                              " of the 'ENDSCALE' keyword) is not yet supported");
                }
                if (optionValue == "IRREVERS") {
                    OPM_THROW(std::runtime_error,
                              "Irreversible end-point scaling (indicated by the 'IRREVERS' option"
                              " of the 'ENDSCALE' keyword) is not yet supported");
                }
            }
        }
    }

    void readGlobalHysteresisOptions_(Opm::DeckConstPtr deck)
    {
        hysteresisConfig_ = std::make_shared<Opm::EclHysteresisConfig>();
        hysteresisConfig_->initFromDeck(deck);
    }

    void readGlobalThreePhaseOptions_(Opm::DeckConstPtr deck)
    {
        bool gasEnabled = deck->hasKeyword("GAS");
        bool oilEnabled = deck->hasKeyword("OIL");
        bool waterEnabled = deck->hasKeyword("WATER");

        int numEnabled =
            (gasEnabled?1:0)
            + (oilEnabled?1:0)
            + (waterEnabled?1:0);

        if (numEnabled < 2)
            OPM_THROW(std::runtime_error,
                      "At least two fluid phases must be enabled. (Is: " << numEnabled << ")");

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
            if (deck->hasKeyword("STONE") || deck->hasKeyword("STONE2"))
                threePhaseApproach_ = Opm::EclStone2Approach;
            else if (deck->hasKeyword("STONE1"))
                threePhaseApproach_ = Opm::EclStone1Approach;
        }
    }

    void initNonElemSpecific_(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        unsigned numSatRegions = static_cast<unsigned>(deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0));
        unsigned numCompressedElems = static_cast<unsigned>(compressedToCartesianElemIdx_.size());

        GasOilEffectiveParamVector gasOilEffectiveParamVector(numSatRegions);
        OilWaterEffectiveParamVector oilWaterEffectiveParamVector(numSatRegions);
        GasOilParamVector gasOilParams(numSatRegions);
        OilWaterParamVector oilWaterParams(numSatRegions);
        MaterialLawParamsVector satRegionParams(numSatRegions);
        EclEpsScalingPointsInfo<Scalar> dummyInfo;
        for (unsigned satnumRegionIdx = 0; satnumRegionIdx < numSatRegions; ++satnumRegionIdx) {
            // the parameters for the effective two-phase matererial laws
            readGasOilEffectiveParameters_(gasOilEffectiveParamVector, deck, eclState, satnumRegionIdx);
            readOilWaterEffectiveParameters_(oilWaterEffectiveParamVector, deck, eclState, satnumRegionIdx);

            auto gasOilDrainParams = std::make_shared<GasOilEpsTwoPhaseParams>();
            gasOilDrainParams->setConfig(oilWaterEclEpsConfig_);
            gasOilDrainParams->setEffectiveLawParams(gasOilEffectiveParamVector[satnumRegionIdx]);
            gasOilDrainParams->finalize();

            gasOilParams[satnumRegionIdx] = std::make_shared<GasOilTwoPhaseHystParams>();
            gasOilParams[satnumRegionIdx]->setConfig(hysteresisConfig_);
            gasOilParams[satnumRegionIdx]->setDrainageParams(gasOilDrainParams, dummyInfo, Opm::EclGasOilSystem);
            gasOilParams[satnumRegionIdx]->finalize();

            auto oilWaterDrainParams = std::make_shared<OilWaterEpsTwoPhaseParams>();
            oilWaterDrainParams->setConfig(oilWaterEclEpsConfig_);
            oilWaterDrainParams->setEffectiveLawParams(oilWaterEffectiveParamVector[satnumRegionIdx]);
            oilWaterDrainParams->finalize();

            oilWaterParams[satnumRegionIdx] = std::make_shared<OilWaterTwoPhaseHystParams>();
            oilWaterParams[satnumRegionIdx]->setConfig(hysteresisConfig_);
            oilWaterParams[satnumRegionIdx]->setDrainageParams(oilWaterDrainParams, dummyInfo, Opm::EclOilWaterSystem);
            oilWaterParams[satnumRegionIdx]->finalize();

            // create the parameter objects for the three-phase law. since we don't have
            // elem specific data here, create one object per PVT region and let the
            // material law parameters for a elem point to its corresponding PVT region object.
            satRegionParams[satnumRegionIdx] = std::make_shared<MaterialLawParams>();

            // find the connate water saturation. this is pretty slow because it does a
            // lot of stuff which not needed, but for the moment it is fast enough
            // because the initialization is not performance critical
            EclEpsScalingPointsInfo<Scalar> epsInfo;
            epsInfo.extractUnscaled(deck, eclState, satnumRegionIdx);

            initThreePhaseParams_(deck,
                                  eclState,
                                  *satRegionParams[satnumRegionIdx],
                                  satnumRegionIdx,
                                  epsInfo,
                                  oilWaterParams[satnumRegionIdx],
                                  gasOilParams[satnumRegionIdx]);

            satRegionParams[satnumRegionIdx]->finalize();
        }

        materialLawParams_.resize(numCompressedElems);
        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            unsigned satnumRegionIdx = static_cast<unsigned>(satnumRegionIdx_[elemIdx]);
            materialLawParams_[elemIdx] = satRegionParams[satnumRegionIdx];
        }
    }

    void initElemSpecific_(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        unsigned numSatRegions = static_cast<unsigned>(deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0));
        unsigned numCompressedElems = static_cast<unsigned>(compressedToCartesianElemIdx_.size());

        // read the end point scaling configuration. this needs to be done only once per
        // deck.
        auto gasOilConfig = std::make_shared<Opm::EclEpsConfig>();
        auto oilWaterConfig = std::make_shared<Opm::EclEpsConfig>();
        gasOilConfig->initFromDeck(deck, eclState, Opm::EclGasOilSystem);
        oilWaterConfig->initFromDeck(deck, eclState, Opm::EclOilWaterSystem);

        // read the saturation region specific parameters from the deck
        GasOilScalingPointsVector gasOilUnscaledPointsVector(numSatRegions);
        OilWaterScalingPointsVector oilWaterUnscaledPointsVector(numSatRegions);
        GasOilEffectiveParamVector gasOilEffectiveParamVector(numSatRegions);
        OilWaterEffectiveParamVector oilWaterEffectiveParamVector(numSatRegions);
        for (unsigned satnumRegionIdx = 0; satnumRegionIdx < numSatRegions; ++satnumRegionIdx) {
            // unscaled points for end-point scaling
            readGasOilUnscaledPoints_(gasOilUnscaledPointsVector, gasOilConfig, deck, eclState, satnumRegionIdx);
            readOilWaterUnscaledPoints_(oilWaterUnscaledPointsVector, oilWaterConfig, deck, eclState, satnumRegionIdx);

            // the parameters for the effective two-phase matererial laws
            readGasOilEffectiveParameters_(gasOilEffectiveParamVector, deck, eclState, satnumRegionIdx);
            readOilWaterEffectiveParameters_(oilWaterEffectiveParamVector, deck, eclState, satnumRegionIdx);

            // read the end point scaling info for the saturation region
            unscaledEpsInfo_[satnumRegionIdx].extractUnscaled(deck, eclState, satnumRegionIdx);

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
            readGasOilScaledPoints_(gasOilScaledInfoVector,
                                    gasOilScaledPointsVector,
                                    gasOilConfig,
                                    epsGridProperties,
                                    elemIdx);
            readOilWaterScaledPoints_(oilWaterScaledEpsInfoDrainage_,
                                      oilWaterScaledEpsPointsDrainage,
                                      oilWaterConfig,
                                      epsGridProperties,
                                      elemIdx);

            if (enableHysteresis()) {
                readGasOilScaledPoints_(gasOilScaledImbInfoVector,
                                        gasOilScaledImbPointsVector,
                                        gasOilConfig,
                                        epsImbGridProperties,
                                        elemIdx);
                readOilWaterScaledPoints_(oilWaterScaledImbInfoVector,
                                          oilWaterScaledImbPointsVector,
                                          oilWaterConfig,
                                          epsImbGridProperties,
                                          elemIdx);
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

        const auto& imbnumData = eclState->getIntGridProperty("IMBNUM")->getData();
        assert(numCompressedElems == satnumRegionIdx_.size());
        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            unsigned satnumRegionIdx = static_cast<unsigned>(satnumRegionIdx_[elemIdx]);

            gasOilParams[elemIdx] = std::make_shared<GasOilTwoPhaseHystParams>();
            oilWaterParams[elemIdx] = std::make_shared<OilWaterTwoPhaseHystParams>();

            gasOilParams[elemIdx]->setConfig(hysteresisConfig_);
            oilWaterParams[elemIdx]->setConfig(hysteresisConfig_);

            auto gasOilDrainParams = std::make_shared<GasOilEpsTwoPhaseParams>();
            gasOilDrainParams->setConfig(gasOilConfig);
            gasOilDrainParams->setUnscaledPoints(gasOilUnscaledPointsVector[satnumRegionIdx]);
            gasOilDrainParams->setScaledPoints(gasOilScaledPointsVector[elemIdx]);
            gasOilDrainParams->setEffectiveLawParams(gasOilEffectiveParamVector[satnumRegionIdx]);
            gasOilDrainParams->finalize();

            auto oilWaterDrainParams = std::make_shared<OilWaterEpsTwoPhaseParams>();
            oilWaterDrainParams->setConfig(oilWaterConfig);
            oilWaterDrainParams->setUnscaledPoints(oilWaterUnscaledPointsVector[satnumRegionIdx]);
            oilWaterDrainParams->setScaledPoints(oilWaterScaledEpsPointsDrainage[elemIdx]);
            oilWaterDrainParams->setEffectiveLawParams(oilWaterEffectiveParamVector[satnumRegionIdx]);
            oilWaterDrainParams->finalize();

            gasOilParams[elemIdx]->setDrainageParams(gasOilDrainParams,
                                                         *gasOilScaledInfoVector[elemIdx],
                                                         EclGasOilSystem);
            oilWaterParams[elemIdx]->setDrainageParams(oilWaterDrainParams,
                                                           *oilWaterScaledEpsInfoDrainage_[elemIdx],
                                                           EclOilWaterSystem);

            if (enableHysteresis()) {
                unsigned imbRegionIdx = static_cast<unsigned>(imbnumData[elemIdx]) - 1;

                auto gasOilImbParamsHyst = std::make_shared<GasOilEpsTwoPhaseParams>();
                gasOilImbParamsHyst->setConfig(gasOilConfig);
                gasOilImbParamsHyst->setUnscaledPoints(gasOilUnscaledPointsVector[imbRegionIdx]);
                gasOilImbParamsHyst->setScaledPoints(gasOilScaledImbPointsVector[elemIdx]);
                gasOilImbParamsHyst->setEffectiveLawParams(gasOilEffectiveParamVector[imbRegionIdx]);
                gasOilImbParamsHyst->finalize();

                auto oilWaterImbParamsHyst = std::make_shared<OilWaterEpsTwoPhaseParams>();
                oilWaterImbParamsHyst->setConfig(oilWaterConfig);
                oilWaterImbParamsHyst->setUnscaledPoints(oilWaterUnscaledPointsVector[imbRegionIdx]);
                oilWaterImbParamsHyst->setScaledPoints(oilWaterScaledImbPointsVector[elemIdx]);
                oilWaterImbParamsHyst->setEffectiveLawParams(oilWaterEffectiveParamVector[imbRegionIdx]);
                oilWaterImbParamsHyst->finalize();

                gasOilParams[elemIdx]->setImbibitionParams(gasOilImbParamsHyst,
                                                               *gasOilScaledImbInfoVector[elemIdx],
                                                               EclGasOilSystem);
                oilWaterParams[elemIdx]->setImbibitionParams(oilWaterImbParamsHyst,
                                                                 *gasOilScaledImbInfoVector[elemIdx],
                                                                 EclGasOilSystem);
            }

            gasOilParams[elemIdx]->finalize();
            oilWaterParams[elemIdx]->finalize();
        }

        // create the parameter objects for the three-phase law
        materialLawParams_.resize(numCompressedElems);
        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            materialLawParams_[elemIdx] = std::make_shared<MaterialLawParams>();
            unsigned satnumRegionIdx = static_cast<unsigned>(satnumRegionIdx_[elemIdx]);

            initThreePhaseParams_(deck,
                                  eclState,
                                  *materialLawParams_[elemIdx],
                                  satnumRegionIdx,
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

    SaturationFunctionFamily getSaturationFunctionFamily(Opm::EclipseStateConstPtr eclState) const
    {
        const auto& tableManager = eclState->getTableManager();
        const std::vector<SwofTable>& swofTables = tableManager->getSwofTables();
        const std::vector<SlgofTable>& slgofTables = tableManager->getSlgofTables();
        const std::vector<SgofTable>& sgofTables = tableManager->getSgofTables();
        const std::vector<SwfnTable>& swfnTables = tableManager->getSwfnTables();
        const std::vector<SgfnTable>& sgfnTables = tableManager->getSgfnTables();
        const std::vector<Sof3Table>& sof3Tables = tableManager->getSof3Tables();

        bool family1 = (!sgofTables.empty() || !slgofTables.empty()) && !swofTables.empty();
        bool family2 = !swfnTables.empty() && !sgfnTables.empty() && !sof3Tables.empty();

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
                                        Opm::DeckConstPtr deck,
                                        Opm::EclipseStateConstPtr eclState,
                                        unsigned satnumRegionIdx)
    {
        dest[satnumRegionIdx] = std::make_shared<GasOilEffectiveTwoPhaseParams>();

        bool hasWater = deck->hasKeyword("WATER");
        bool hasGas = deck->hasKeyword("GAS");
        bool hasOil = deck->hasKeyword("OIL");

        auto& effParams = *dest[satnumRegionIdx];

        // the situation for the gas phase is complicated that all saturations are
        // shifted by the connate water saturation.
        Scalar Swco = unscaledEpsInfo_[satnumRegionIdx].Swl;

        // handle the twophase case
        const auto& tableManager = eclState->getTableManager();
        if (!hasWater) {
            if (!tableManager->getSgofTables().empty())
                readGasOilEffectiveParametersSgof_(effParams,
                                                   Swco,
                                                   tableManager->getSgofTables()[satnumRegionIdx]);
            else {
                assert(!tableManager->getSlgofTables().empty());
                readGasOilEffectiveParametersSlgof_(effParams,
                                                    Swco,
                                                    tableManager->getSlgofTables()[satnumRegionIdx]);
            }

            // Todo (?): support for twophase simulations using family2?
            return;
        }
        else if (!hasGas) {
            return;
        }

        // so far, only water-oil and oil-gas simulations are supported, i.e.,
        // there's no gas-water yet.
        if (!hasWater || !hasGas || !hasOil)
            throw std::domain_error("The specified phase configuration is not suppored");

        switch (getSaturationFunctionFamily(eclState)) {
        case FamilyI:
        {
            if (!tableManager->getSgofTables().empty())
                readGasOilEffectiveParametersSgof_(effParams,
                                                   Swco,
                                                   tableManager->getSgofTables()[satnumRegionIdx]);
            else if (!tableManager->getSlgofTables().empty())
                readGasOilEffectiveParametersSlgof_(effParams,
                                                    Swco,
                                                    tableManager->getSlgofTables()[satnumRegionIdx]);

            break;
        }

        case FamilyII:
        {
            readGasOilEffectiveParametersFamily2_(effParams,
                                                  Swco,
                                                  tableManager->getSof3Tables()[satnumRegionIdx],
                                                  tableManager->getSgfnTables()[satnumRegionIdx]);
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
            SoSamples[sampleIdx] = 1 - sgofTable.getSgColumn()[sampleIdx];
            SoKroSamples[sampleIdx] = SoSamples[sampleIdx] - Swco;
        }

        effParams.setKrwSamples(SoKroSamples, sgofTable.getKrogColumn());
        effParams.setKrnSamples(SoSamples, sgofTable.getKrgColumn());
        effParams.setPcnwSamples(SoSamples, sgofTable.getPcogColumn());
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
            SoSamples[sampleIdx] = slgofTable.getSlColumn()[sampleIdx];
            SoKroSamples[sampleIdx] = slgofTable.getSlColumn()[sampleIdx] - Swco;
        }

        effParams.setKrwSamples(SoKroSamples, slgofTable.getKrogColumn());
        effParams.setKrnSamples(SoSamples, slgofTable.getKrgColumn());
        effParams.setPcnwSamples(SoSamples, slgofTable.getPcogColumn());
        effParams.finalize();
    }

    void readGasOilEffectiveParametersFamily2_(GasOilEffectiveTwoPhaseParams& effParams,
                                               Scalar /* Swco */,
                                               const Opm::Sof3Table& sof3Table,
                                               const Opm::SgfnTable& sgfnTable)
    {
        // convert the saturations of the SGFN keyword from gas to oil saturations
        std::vector<double> SoSamples(sgfnTable.numRows());
        const auto &SoColumn = sof3Table.getSoColumn();
        for (size_t sampleIdx = 0; sampleIdx < sgfnTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = 1 - sgfnTable.getSgColumn()[sampleIdx];
        }

        effParams.setKrwSamples(SoColumn, sof3Table.getKrogColumn());
        effParams.setKrnSamples(SoSamples, sgfnTable.getKrgColumn());
        effParams.setPcnwSamples(SoSamples, sgfnTable.getPcogColumn());
        effParams.finalize();
    }

    template <class Container>
    void readOilWaterEffectiveParameters_(Container& dest,
                                          Opm::DeckConstPtr deck,
                                          Opm::EclipseStateConstPtr eclState,
                                          unsigned satnumRegionIdx)
    {
        dest[satnumRegionIdx] = std::make_shared<OilWaterEffectiveTwoPhaseParams>();

        bool hasWater = deck->hasKeyword("WATER");
        bool hasGas = deck->hasKeyword("GAS");
        bool hasOil = deck->hasKeyword("OIL");

        const auto tableManager = eclState->getTableManager();
        auto& effParams = *dest[satnumRegionIdx];

        // handle the twophase case
        if (!hasWater) {
            return;
        }
        else if (!hasGas) {
            const auto& swofTable = tableManager->getSwofTables()[satnumRegionIdx];
            const auto &SwColumn = swofTable.getSwColumn();

            effParams.setKrwSamples(SwColumn, swofTable.getKrwColumn());
            effParams.setKrnSamples(SwColumn, swofTable.getKrowColumn());
            effParams.setPcnwSamples(SwColumn, swofTable.getPcowColumn());
            effParams.finalize();

            // Todo (?): support for twophase simulations using family2?
            return;
        }

        // so far, only water-oil and oil-gas simulations are supported, i.e.,
        // there's no gas-water yet.
        if (!hasWater || !hasGas || !hasOil)
            throw std::domain_error("The specified phase configuration is not suppored");

        switch (getSaturationFunctionFamily(eclState)) {
        case FamilyI: {
            const auto& swofTable = tableManager->getSwofTables()[satnumRegionIdx];
            const auto &SwColumn = swofTable.getSwColumn();

            effParams.setKrwSamples(SwColumn, swofTable.getKrwColumn());
            effParams.setKrnSamples(SwColumn, swofTable.getKrowColumn());
            effParams.setPcnwSamples(SwColumn, swofTable.getPcowColumn());
            effParams.finalize();
            break;
        }
        case FamilyII:
        {
            const auto& swfnTable = tableManager->getSwfnTables()[satnumRegionIdx];
            const auto& sof3Table = tableManager->getSof3Tables()[satnumRegionIdx];
            const auto &SwColumn = swfnTable.getSwColumn();

            // convert the saturations of the SOF3 keyword from oil to water saturations
            std::vector<double> SwSamples(sof3Table.numRows());
            for (size_t sampleIdx = 0; sampleIdx < sof3Table.numRows(); ++ sampleIdx)
                SwSamples[sampleIdx] = 1 - sof3Table.getSoColumn()[sampleIdx];

            effParams.setKrwSamples(SwColumn, swfnTable.getKrwColumn());
            effParams.setKrnSamples(SwSamples, sof3Table.getKrowColumn());
            effParams.setPcnwSamples(SwColumn, swfnTable.getPcowColumn());
            effParams.finalize();
            break;
        }

        case noFamily:
        //default:
            throw std::domain_error("No valid saturation keyword family specified");

        }
    }


    template <class Container>
    void readGasOilUnscaledPoints_(Container &dest,
                                   std::shared_ptr<EclEpsConfig> config,
                                   Opm::DeckConstPtr /* deck */,
                                   Opm::EclipseStateConstPtr /* eclState */,
                                   unsigned satnumRegionIdx)
    {
        dest[satnumRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satnumRegionIdx]->init(unscaledEpsInfo_[satnumRegionIdx], *config, EclGasOilSystem);
    }

    template <class Container>
    void readOilWaterUnscaledPoints_(Container &dest,
                                     std::shared_ptr<EclEpsConfig> config,
                                     Opm::DeckConstPtr /* deck */,
                                     Opm::EclipseStateConstPtr /* eclState */,
                                     unsigned satnumRegionIdx)
    {
        dest[satnumRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satnumRegionIdx]->init(unscaledEpsInfo_[satnumRegionIdx], *config, EclOilWaterSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readGasOilScaledPoints_(InfoContainer& destInfo,
                                 PointsContainer& destPoints,
                                 std::shared_ptr<EclEpsConfig> config,
                                 const EclEpsGridProperties& epsGridProperties,
                                 unsigned elemIdx)
    {
        unsigned cartElemIdx = static_cast<unsigned>(compressedToCartesianElemIdx_[elemIdx]);
        unsigned satnumRegionIdx = (*epsGridProperties.satnum)[cartElemIdx] - 1; // ECL uses Fortran indices!

        destInfo[elemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satnumRegionIdx]);
        destInfo[elemIdx]->extractScaled(epsGridProperties, cartElemIdx);

        destPoints[elemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[elemIdx]->init(*destInfo[elemIdx], *config, EclGasOilSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readOilWaterScaledPoints_(InfoContainer& destInfo,
                                   PointsContainer& destPoints,
                                   std::shared_ptr<EclEpsConfig> config,
                                   const EclEpsGridProperties& epsGridProperties,
                                   unsigned elemIdx)
    {
        unsigned cartElemIdx = static_cast<unsigned>(compressedToCartesianElemIdx_[elemIdx]);
        unsigned satnumRegionIdx = (*epsGridProperties.satnum)[cartElemIdx] - 1; // ECL uses Fortran indices!

        destInfo[elemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satnumRegionIdx]);;
        destInfo[elemIdx]->extractScaled(epsGridProperties, cartElemIdx);

        destPoints[elemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[elemIdx]->init(*destInfo[elemIdx], *config, EclOilWaterSystem);
    }

    void initThreePhaseParams_(Opm::DeckConstPtr deck,
                               Opm::EclipseStateConstPtr /* eclState */,
                               MaterialLawParams& materialParams,
                               unsigned satnumIdx,
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
            realParams.setSowcr(epsInfo.Sowcr);

            if (deck->hasKeyword("STONE1EX")) {
                Scalar eta =
                    deck->getKeyword("STONE1EX")->getRecord(satnumIdx)->getItem(0)->getSIDouble(0);
                realParams.setSogcr(eta);
            }
            else
                realParams.setSogcr(1.0);
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

    EclEpsScalingPoints<Scalar>& getOilWaterScaledEpsPointsDrainage_(unsigned elemIdx)
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
        }
    }

    bool enableEndPointScaling_;
    std::shared_ptr<EclHysteresisConfig> hysteresisConfig_;

    std::shared_ptr<EclEpsConfig> oilWaterEclEpsConfig_;
    std::vector<Opm::EclEpsScalingPointsInfo<Scalar>> unscaledEpsInfo_;
    OilWaterScalingInfoVector oilWaterScaledEpsInfoDrainage_;

    Opm::EclMultiplexerApproach threePhaseApproach_;

    // this attribute only makes sense for twophase simulations!
    enum EclTwoPhaseApproach twoPhaseApproach_;

    std::vector<std::shared_ptr<MaterialLawParams> > materialLawParams_;

    std::vector<int> compressedToCartesianElemIdx_;
    std::vector<int> satnumRegionIdx_;
};
} // namespace Opm

#endif
