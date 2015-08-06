// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2015 by Andreas Lauser

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
        int numSatRegions = deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0);
        unsigned numCompressedElems = compressedToCartesianElemIdx.size();;

        // copy the SATNUM grid property. in some cases this is not necessary, but it
        // should not require much memory anyway...
        satnumRegionIdx_.resize(numCompressedElems);
        if (eclState->hasIntGridProperty("SATNUM")) {
            const auto& satnumRawData = eclState->getIntGridProperty("SATNUM")->getData();
            for (int elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
                int cartesianElemIdx = compressedToCartesianElemIdx_[elemIdx];
                satnumRegionIdx_[elemIdx] = satnumRawData[cartesianElemIdx] - 1;
            }
        }
        else
            std::fill(satnumRegionIdx_.begin(), satnumRegionIdx_.end(), 0);

        readGlobalEpsOptions_(deck, eclState);
        readGlobalHysteresisOptions_(deck);
        readGlobalThreePhaseOptions_(deck);

        unscaledEpsInfo_.resize(numSatRegions);
        for (int satnumRegionIdx = 0; satnumRegionIdx < numSatRegions; ++satnumRegionIdx)
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
    Scalar applySwatinit(int elemIdx,
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
                auto& elemEclEpsScalingPoints = *oilWaterScaledEpsPointsDrainage_[elemIdx];
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

    MaterialLawParams& materialLawParams(int elemIdx)
    {
        assert(0 <= elemIdx && elemIdx < (int) materialLawParams_.size());

        int paramIdx;
        if (hasElementSpecificParameters())
            paramIdx = elemIdx;
        else
            paramIdx = satnumRegionIdx_[elemIdx];

        return *materialLawParams_[paramIdx];
    }

    const MaterialLawParams& materialLawParams(int elemIdx) const
    {
        assert(0 <= elemIdx && elemIdx < materialLawParams_.size());

        int paramIdx;
        if (hasElementSpecificParameters())
            paramIdx = elemIdx;
        else
            paramIdx = satnumRegionIdx_[elemIdx];

        return *materialLawParams_[paramIdx];
    }

    template <class FluidState>
    void updateHysteresis(const FluidState& fluidState, int elemIdx)
    {
        if (!enableHysteresis())
            return;

        auto threePhaseParams = materialLawParams_[elemIdx];
        MaterialLaw::updateHysteresis(*threePhaseParams, fluidState);
    }

    const Opm::EclEpsScalingPointsInfo<Scalar>& oilWaterScaledEpsInfoDrainage(int elemIdx) const
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
        threePhaseApproach_ = Opm::EclDefaultApproach;
        if (deck->hasKeyword("STONE") || deck->hasKeyword("STONE2"))
            threePhaseApproach_ = Opm::EclStone2Approach;
        else if (deck->hasKeyword("STONE1"))
            threePhaseApproach_ = Opm::EclStone1Approach;
    }

    void initNonElemSpecific_(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        unsigned numSatRegions = deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0);
        unsigned numCompressedElems = compressedToCartesianElemIdx_.size();;

        GasOilEffectiveParamVector gasOilEffectiveParamVector(numSatRegions);
        OilWaterEffectiveParamVector oilWaterEffectiveParamVector(numSatRegions);
        GasOilParamVector gasOilParams(numSatRegions);
        OilWaterParamVector oilWaterParams(numSatRegions);
        MaterialLawParamsVector satRegionParams(numSatRegions);
        EclEpsScalingPointsInfo<Scalar> dummyInfo;
        for (int satnumRegionIdx = 0; satnumRegionIdx < numSatRegions; ++satnumRegionIdx) {
            // the parameters for the effective two-phase matererial laws
            readGasOilEffectiveParameters_(gasOilEffectiveParamVector, eclState, satnumRegionIdx);
            readOilWaterEffectiveParameters_(oilWaterEffectiveParamVector, eclState, satnumRegionIdx);

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
        for (int elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            int satnumRegionIdx = satnumRegionIdx_[elemIdx];
            materialLawParams_[elemIdx] = satRegionParams[satnumRegionIdx];
        }
    }

    void initElemSpecific_(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        unsigned numSatRegions = deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0);
        unsigned numCompressedElems = compressedToCartesianElemIdx_.size();;

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
            readGasOilEffectiveParameters_(gasOilEffectiveParamVector, eclState, satnumRegionIdx);
            readOilWaterEffectiveParameters_(oilWaterEffectiveParamVector, eclState, satnumRegionIdx);

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
        oilWaterScaledEpsPointsDrainage_.resize(numCompressedElems);
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
                                      oilWaterScaledEpsPointsDrainage_,
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
            int satnumRegionIdx = satnumRegionIdx_[elemIdx];

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
            oilWaterDrainParams->setScaledPoints(oilWaterScaledEpsPointsDrainage_[elemIdx]);
            oilWaterDrainParams->setEffectiveLawParams(oilWaterEffectiveParamVector[satnumRegionIdx]);
            oilWaterDrainParams->finalize();

            gasOilParams[elemIdx]->setDrainageParams(gasOilDrainParams,
                                                         *gasOilScaledInfoVector[elemIdx],
                                                         EclGasOilSystem);
            oilWaterParams[elemIdx]->setDrainageParams(oilWaterDrainParams,
                                                           *oilWaterScaledEpsInfoDrainage_[elemIdx],
                                                           EclOilWaterSystem);

            if (enableHysteresis()) {
                int imbRegionIdx = imbnumData[elemIdx] - 1;

                auto gasOilImbParams = std::make_shared<GasOilEpsTwoPhaseParams>();
                gasOilImbParams->setConfig(gasOilConfig);
                gasOilImbParams->setUnscaledPoints(gasOilUnscaledPointsVector[imbRegionIdx]);
                gasOilImbParams->setScaledPoints(gasOilScaledImbPointsVector[elemIdx]);
                gasOilImbParams->setEffectiveLawParams(gasOilEffectiveParamVector[imbRegionIdx]);
                gasOilImbParams->finalize();

                auto oilWaterImbParams = std::make_shared<OilWaterEpsTwoPhaseParams>();
                oilWaterImbParams->setConfig(oilWaterConfig);
                oilWaterImbParams->setUnscaledPoints(oilWaterUnscaledPointsVector[imbRegionIdx]);
                oilWaterImbParams->setScaledPoints(oilWaterScaledImbPointsVector[elemIdx]);
                oilWaterImbParams->setEffectiveLawParams(oilWaterEffectiveParamVector[imbRegionIdx]);
                oilWaterImbParams->finalize();

                gasOilParams[elemIdx]->setImbibitionParams(gasOilImbParams,
                                                               *gasOilScaledImbInfoVector[elemIdx],
                                                               EclGasOilSystem);
                oilWaterParams[elemIdx]->setImbibitionParams(oilWaterImbParams,
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
            int satnumRegionIdx = satnumRegionIdx_[elemIdx];

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

    template <class Container>
    void readGasOilEffectiveParameters_(Container& dest,
                                        Opm::EclipseStateConstPtr eclState,
                                        int satnumRegionIdx)
    {
        dest[satnumRegionIdx] = std::make_shared<GasOilEffectiveTwoPhaseParams>();

        auto& effParams = *dest[satnumRegionIdx];
        const auto& sgofTable = eclState->getSgofTables()[satnumRegionIdx];

        // convert the saturations of the SGOF keyword from gas to oil saturations
        std::vector<double> SoSamples(sgofTable.numRows());
        for (size_t sampleIdx = 0; sampleIdx < sgofTable.numRows(); ++ sampleIdx)
            SoSamples[sampleIdx] = 1 - sgofTable.getSgColumn()[sampleIdx];

        effParams.setKrwSamples(SoSamples, sgofTable.getKrogColumn());
        effParams.setKrnSamples(SoSamples, sgofTable.getKrgColumn());
        effParams.setPcnwSamples(SoSamples, sgofTable.getPcogColumn());
        effParams.finalize();

    }

    template <class Container>
    void readOilWaterEffectiveParameters_(Container& dest,
                                          Opm::EclipseStateConstPtr eclState,
                                          int satnumRegionIdx)
    {
        dest[satnumRegionIdx] = std::make_shared<OilWaterEffectiveTwoPhaseParams>();

        auto& effParams = *dest[satnumRegionIdx];
        const auto& swofTable = eclState->getSwofTables()[satnumRegionIdx];

        const auto &SwColumn = swofTable.getSwColumn();

        effParams.setKrwSamples(SwColumn, swofTable.getKrwColumn());
        effParams.setKrnSamples(SwColumn, swofTable.getKrowColumn());
        effParams.setPcnwSamples(SwColumn, swofTable.getPcowColumn());
        effParams.finalize();
    }

    template <class Container>
    void readGasOilUnscaledPoints_(Container &dest,
                                   std::shared_ptr<EclEpsConfig> config,
                                   Opm::DeckConstPtr deck,
                                   Opm::EclipseStateConstPtr eclState,
                                   int satnumRegionIdx)
    {
        dest[satnumRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satnumRegionIdx]->init(unscaledEpsInfo_[satnumRegionIdx], *config, EclGasOilSystem);
    }

    template <class Container>
    void readOilWaterUnscaledPoints_(Container &dest,
                                     std::shared_ptr<EclEpsConfig> config,
                                     Opm::DeckConstPtr deck,
                                     Opm::EclipseStateConstPtr eclState,
                                     int satnumRegionIdx)
    {
        dest[satnumRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satnumRegionIdx]->init(unscaledEpsInfo_[satnumRegionIdx], *config, EclOilWaterSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readGasOilScaledPoints_(InfoContainer& destInfo,
                                 PointsContainer& destPoints,
                                 std::shared_ptr<EclEpsConfig> config,
                                 const EclEpsGridProperties& epsGridProperties,
                                 int elemIdx)
    {
        int satnumRegionIdx = (*epsGridProperties.satnum)[elemIdx] - 1; // ECL uses Fortran indices!
        int cartElemIdx = compressedToCartesianElemIdx_[elemIdx];

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
                                   int elemIdx)
    {
        int satnumRegionIdx = (*epsGridProperties.satnum)[elemIdx] - 1; // ECL uses Fortran indices!
        int cartElemIdx = compressedToCartesianElemIdx_[elemIdx];

        destInfo[elemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satnumRegionIdx]);;
        destInfo[elemIdx]->extractScaled(epsGridProperties, cartElemIdx);

        destPoints[elemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[elemIdx]->init(*destInfo[elemIdx], *config, EclOilWaterSystem);
    }

    void initThreePhaseParams_(Opm::DeckConstPtr deck,
                               Opm::EclipseStateConstPtr eclState,
                               MaterialLawParams& materialParams,
                               int satnumIdx,
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
        }
    }

    bool enableEndPointScaling_;
    std::shared_ptr<EclHysteresisConfig> hysteresisConfig_;

    std::shared_ptr<EclEpsConfig> oilWaterEclEpsConfig_;
    std::vector<Opm::EclEpsScalingPointsInfo<Scalar>> unscaledEpsInfo_;
    OilWaterScalingInfoVector oilWaterScaledEpsInfoDrainage_;
    OilWaterScalingPointsVector oilWaterScaledEpsPointsDrainage_;

    EclMultiplexerApproach threePhaseApproach_;
    std::vector<std::shared_ptr<MaterialLawParams> > materialLawParams_;

    std::vector<int> compressedToCartesianElemIdx_;
    std::vector<int> satnumRegionIdx_;
};
} // namespace Opm

#endif
