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
                      Opm::EclipseStateConstPtr eclState)
    {
        // get the number of saturation regions and the number of cells in the deck
        int numSatRegions = deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0);
        int numCartesianElements = eclState->getEclipseGrid()->getCartesianSize();

        // copy the SATNUM grid property. in some cases this is not necessary, but it
        // should not require much memory anyway...
        if (eclState->hasIntGridProperty("SATNUM"))
            satnumData_ = eclState->getIntGridProperty("SATNUM")->getData();
        else
            satnumData_.resize(numCartesianElements, 1);

        readGlobalEpsOptions_(deck, eclState);
        readGlobalHysteresisOptions_(deck);
        readGlobalThreePhaseOptions_(deck);

        unscaledEpsInfo_.resize(numSatRegions);
        for (int satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
            unscaledEpsInfo_[satRegionIdx].extractUnscaled(deck, eclState, satRegionIdx);
        }

        if (!hasElementSpecificParameters())
            initNonCellSpecific_(deck, eclState);
        else
            initCellSpecific_(deck, eclState);
    }

    /*!
     * \brief Modify the initial condition according to the SWATINIT keyword.
     *
     * The method returns the water saturation which yields a givenn capillary
     * pressure. The reason this method is not folded directly into initFromDeck() is
     * that the capillary pressure given depends on the particuars of how the simulator
     * calculates its initial condition.
     */
    Scalar applySwatinit(int cartesianCellIdx,
                         Scalar pcow,
                         Scalar Sw)
    {
        auto& cellScaledEpsInfo = *oilWaterScaledEpsInfoDrainage_[cartesianCellIdx];

        // TODO: Mixed wettability systems - see ecl kw OPTIONS switch 74
        if (Sw <= cellScaledEpsInfo.Swl)
            Sw = cellScaledEpsInfo.Swl;
        else if (pcow < 0.0)
            Sw = cellScaledEpsInfo.Swu;
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
            MaterialLaw::capillaryPressures(pc, materialLawParams(cartesianCellIdx), fs);

            Scalar pcowAtSw = pc[oilPhaseIdx] - pc[waterPhaseIdx];
            if (pcowAtSw > 0.0) {
                cellScaledEpsInfo.maxPcow *= pcow/pcowAtSw;
                auto& cellEclEpsScalingPoints = *oilWaterScaledEpsPointsDrainage_[cartesianCellIdx];
                cellEclEpsScalingPoints.init(cellScaledEpsInfo, *oilWaterEclEpsConfig_, Opm::EclOilWaterSystem);
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

    MaterialLawParams& materialLawParams(int cartesianCellIdx)
    {
        assert(0 <= cartesianCellIdx && cartesianCellIdx < (int) materialLawParams_.size());

        int paramIdx;
        if (hasElementSpecificParameters())
            paramIdx = cartesianCellIdx;
        else
            paramIdx = satnumData_[cartesianCellIdx] - 1;

        return *materialLawParams_[paramIdx];
    }

    const MaterialLawParams& materialLawParams(int cartesianCellIdx) const
    {
        assert(0 <= cartesianCellIdx && cartesianCellIdx < materialLawParams_.size());

        int paramIdx;
        if (hasElementSpecificParameters())
            paramIdx = cartesianCellIdx;
        else
            paramIdx = satnumData_[cartesianCellIdx] - 1;

        return *materialLawParams_[paramIdx];
    }

    template <class FluidState>
    void updateHysteresis(const FluidState& fluidState, int cartesianCellIdx)
    {
        if (!enableHysteresis())
            return;

        auto threePhaseParams = materialLawParams_[cartesianCellIdx];
        MaterialLaw::updateHysteresis(*threePhaseParams, fluidState);
    }

    const Opm::EclEpsScalingPointsInfo<Scalar>& oilWaterScaledEpsInfoDrainage(int cartesianCellIdx) const
    {
        if (hasElementSpecificParameters())
            return *oilWaterScaledEpsInfoDrainage_[cartesianCellIdx];

        int satRegionIdx = satnumData_[cartesianCellIdx] - 1;
        return unscaledEpsInfo_[satRegionIdx];
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

    void initNonCellSpecific_(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        int numSatRegions = deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0);
        int numCartesianElements = eclState->getEclipseGrid()->getCartesianSize();

        GasOilEffectiveParamVector gasOilEffectiveParamVector(numSatRegions);
        OilWaterEffectiveParamVector oilWaterEffectiveParamVector(numSatRegions);
        GasOilParamVector gasOilParams(numSatRegions);
        OilWaterParamVector oilWaterParams(numSatRegions);
        MaterialLawParamsVector satRegionParams(numSatRegions);
        EclEpsScalingPointsInfo<Scalar> dummyInfo;
        for (int satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
            // the parameters for the effective two-phase matererial laws
            readGasOilEffectiveParameters_(gasOilEffectiveParamVector, eclState, satRegionIdx);
            readOilWaterEffectiveParameters_(oilWaterEffectiveParamVector, eclState, satRegionIdx);

            auto gasOilDrainParams = std::make_shared<GasOilEpsTwoPhaseParams>();
            gasOilDrainParams->setConfig(oilWaterEclEpsConfig_);
            gasOilDrainParams->setEffectiveLawParams(gasOilEffectiveParamVector[satRegionIdx]);
            gasOilDrainParams->finalize();

            gasOilParams[satRegionIdx] = std::make_shared<GasOilTwoPhaseHystParams>();
            gasOilParams[satRegionIdx]->setConfig(hysteresisConfig_);
            gasOilParams[satRegionIdx]->setDrainageParams(gasOilDrainParams, dummyInfo, Opm::EclGasOilSystem);
            gasOilParams[satRegionIdx]->finalize();

            auto oilWaterDrainParams = std::make_shared<OilWaterEpsTwoPhaseParams>();
            oilWaterDrainParams->setConfig(oilWaterEclEpsConfig_);
            oilWaterDrainParams->setEffectiveLawParams(oilWaterEffectiveParamVector[satRegionIdx]);
            oilWaterDrainParams->finalize();

            oilWaterParams[satRegionIdx] = std::make_shared<OilWaterTwoPhaseHystParams>();
            oilWaterParams[satRegionIdx]->setConfig(hysteresisConfig_);
            oilWaterParams[satRegionIdx]->setDrainageParams(oilWaterDrainParams, dummyInfo, Opm::EclOilWaterSystem);
            oilWaterParams[satRegionIdx]->finalize();

            // create the parameter objects for the three-phase law. since we don't have
            // cell specific data here, create one object per PVT region and let the
            // material law parameters for a cell point to its corresponding PVT region object.
            satRegionParams[satRegionIdx] = std::make_shared<MaterialLawParams>();

            // find the connate water saturation. this is pretty slow because it does a
            // lot of stuff which not needed, but for the moment it is fast enough
            // because the initialization is not performance critical
            EclEpsScalingPointsInfo<Scalar> epsInfo;
            epsInfo.extractUnscaled(deck, eclState, satRegionIdx);

            initThreePhaseParams_(deck,
                                  eclState,
                                  *satRegionParams[satRegionIdx],
                                  satRegionIdx,
                                  epsInfo,
                                  oilWaterParams[satRegionIdx],
                                  gasOilParams[satRegionIdx]);

            satRegionParams[satRegionIdx]->finalize();
        }

        materialLawParams_.resize(numCartesianElements);
        for (int cartElemIdx = 0; cartElemIdx < numCartesianElements; ++cartElemIdx) {
            int satRegionIdx = satnumData_[cartElemIdx] - 1;
            materialLawParams_[cartElemIdx] = satRegionParams[satRegionIdx];
        }
    }

    void initCellSpecific_(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        unsigned numSatRegions = deck->getKeyword("TABDIMS")->getRecord(0)->getItem("NTSFUN")->getInt(0);
        unsigned numCartesianElements = eclState->getEclipseGrid()->getCartesianSize();

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
        for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
            // unscaled points for end-point scaling
            readGasOilUnscaledPoints_(gasOilUnscaledPointsVector, gasOilConfig, deck, eclState, satRegionIdx);
            readOilWaterUnscaledPoints_(oilWaterUnscaledPointsVector, oilWaterConfig, deck, eclState, satRegionIdx);

            // the parameters for the effective two-phase matererial laws
            readGasOilEffectiveParameters_(gasOilEffectiveParamVector, eclState, satRegionIdx);
            readOilWaterEffectiveParameters_(oilWaterEffectiveParamVector, eclState, satRegionIdx);

            // read the end point scaling info for the saturation region
            unscaledEpsInfo_[satRegionIdx].extractUnscaled(deck, eclState, satRegionIdx);

        }

        // read the scaled end point scaling parameters which are specific for each
        // logically Cartesian cell
        GasOilScalingInfoVector gasOilScaledInfoVector(numCartesianElements);
        oilWaterScaledEpsInfoDrainage_.resize(numCartesianElements);
        GasOilScalingInfoVector gasOilScaledImbInfoVector;
        OilWaterScalingInfoVector oilWaterScaledImbInfoVector;

        GasOilScalingPointsVector gasOilScaledPointsVector(numCartesianElements);
        oilWaterScaledEpsPointsDrainage_.resize(numCartesianElements);
        GasOilScalingPointsVector gasOilScaledImbPointsVector;
        OilWaterScalingPointsVector oilWaterScaledImbPointsVector;

        if (enableHysteresis()) {
            gasOilScaledImbInfoVector.resize(numCartesianElements);
            gasOilScaledImbPointsVector.resize(numCartesianElements);
            oilWaterScaledImbInfoVector.resize(numCartesianElements);
            oilWaterScaledImbPointsVector.resize(numCartesianElements);
        }

        EclEpsGridProperties epsGridProperties, epsImbGridProperties;
        epsGridProperties.initFromDeck(deck, eclState, /*imbibition=*/false);
        if (enableHysteresis())
            epsImbGridProperties.initFromDeck(deck, eclState, /*imbibition=*/true);
        for (unsigned cartElemIdx = 0; cartElemIdx < numCartesianElements; ++cartElemIdx) {
            readGasOilScaledPoints_(gasOilScaledInfoVector,
                                    gasOilScaledPointsVector,
                                    gasOilConfig,
                                    epsGridProperties,
                                    cartElemIdx);
            readOilWaterScaledPoints_(oilWaterScaledEpsInfoDrainage_,
                                      oilWaterScaledEpsPointsDrainage_,
                                      oilWaterConfig,
                                      epsGridProperties,
                                      cartElemIdx);

            if (enableHysteresis()) {
                readGasOilScaledPoints_(gasOilScaledImbInfoVector,
                                        gasOilScaledImbPointsVector,
                                        gasOilConfig,
                                        epsImbGridProperties,
                                        cartElemIdx);
                readOilWaterScaledPoints_(oilWaterScaledImbInfoVector,
                                          oilWaterScaledImbPointsVector,
                                          oilWaterConfig,
                                          epsImbGridProperties,
                                          cartElemIdx);
            }
        }

        // create the parameter objects for the two-phase laws
        GasOilParamVector gasOilParams(numCartesianElements);
        OilWaterParamVector oilWaterParams(numCartesianElements);
        GasOilParamVector gasOilImbParams;
        OilWaterParamVector oilWaterImbParams;

        if (enableHysteresis()) {
            gasOilImbParams.resize(numCartesianElements);
            oilWaterImbParams.resize(numCartesianElements);
        }

        const auto& imbnumData = eclState->getIntGridProperty("IMBNUM")->getData();
        assert(numCartesianElements == satnumData_.size());
        for (unsigned cartElemIdx = 0; cartElemIdx < numCartesianElements; ++cartElemIdx) {
            int satRegionIdx = satnumData_[cartElemIdx] - 1;

            gasOilParams[cartElemIdx] = std::make_shared<GasOilTwoPhaseHystParams>();
            oilWaterParams[cartElemIdx] = std::make_shared<OilWaterTwoPhaseHystParams>();

            gasOilParams[cartElemIdx]->setConfig(hysteresisConfig_);
            oilWaterParams[cartElemIdx]->setConfig(hysteresisConfig_);

            auto gasOilDrainParams = std::make_shared<GasOilEpsTwoPhaseParams>();
            gasOilDrainParams->setConfig(gasOilConfig);
            gasOilDrainParams->setUnscaledPoints(gasOilUnscaledPointsVector[satRegionIdx]);
            gasOilDrainParams->setScaledPoints(gasOilScaledPointsVector[cartElemIdx]);
            gasOilDrainParams->setEffectiveLawParams(gasOilEffectiveParamVector[satRegionIdx]);
            gasOilDrainParams->finalize();

            auto oilWaterDrainParams = std::make_shared<OilWaterEpsTwoPhaseParams>();
            oilWaterDrainParams->setConfig(oilWaterConfig);
            oilWaterDrainParams->setUnscaledPoints(oilWaterUnscaledPointsVector[satRegionIdx]);
            oilWaterDrainParams->setScaledPoints(oilWaterScaledEpsPointsDrainage_[cartElemIdx]);
            oilWaterDrainParams->setEffectiveLawParams(oilWaterEffectiveParamVector[satRegionIdx]);
            oilWaterDrainParams->finalize();

            gasOilParams[cartElemIdx]->setDrainageParams(gasOilDrainParams,
                                                         *gasOilScaledInfoVector[cartElemIdx],
                                                         EclGasOilSystem);
            oilWaterParams[cartElemIdx]->setDrainageParams(oilWaterDrainParams,
                                                           *oilWaterScaledEpsInfoDrainage_[cartElemIdx],
                                                           EclOilWaterSystem);

            if (enableHysteresis()) {
                int imbRegionIdx = imbnumData[cartElemIdx] - 1;

                auto gasOilImbParams = std::make_shared<GasOilEpsTwoPhaseParams>();
                gasOilImbParams->setConfig(gasOilConfig);
                gasOilImbParams->setUnscaledPoints(gasOilUnscaledPointsVector[imbRegionIdx]);
                gasOilImbParams->setScaledPoints(gasOilScaledImbPointsVector[cartElemIdx]);
                gasOilImbParams->setEffectiveLawParams(gasOilEffectiveParamVector[imbRegionIdx]);
                gasOilImbParams->finalize();

                auto oilWaterImbParams = std::make_shared<OilWaterEpsTwoPhaseParams>();
                oilWaterImbParams->setConfig(oilWaterConfig);
                oilWaterImbParams->setUnscaledPoints(oilWaterUnscaledPointsVector[imbRegionIdx]);
                oilWaterImbParams->setScaledPoints(oilWaterScaledImbPointsVector[cartElemIdx]);
                oilWaterImbParams->setEffectiveLawParams(oilWaterEffectiveParamVector[imbRegionIdx]);
                oilWaterImbParams->finalize();

                gasOilParams[cartElemIdx]->setImbibitionParams(gasOilImbParams,
                                                               *gasOilScaledImbInfoVector[cartElemIdx],
                                                               EclGasOilSystem);
                oilWaterParams[cartElemIdx]->setImbibitionParams(oilWaterImbParams,
                                                                 *gasOilScaledImbInfoVector[cartElemIdx],
                                                                 EclGasOilSystem);
            }

            gasOilParams[cartElemIdx]->finalize();
            oilWaterParams[cartElemIdx]->finalize();
        }

        // create the parameter objects for the three-phase law
        materialLawParams_.resize(numCartesianElements);
        for (unsigned cartElemIdx = 0; cartElemIdx < numCartesianElements; ++cartElemIdx) {
            materialLawParams_[cartElemIdx] = std::make_shared<MaterialLawParams>();
            int satRegionIdx = satnumData_[cartElemIdx] - 1;

            initThreePhaseParams_(deck,
                                  eclState,
                                  *materialLawParams_[cartElemIdx],
                                  satRegionIdx,
                                  *oilWaterScaledEpsInfoDrainage_[cartElemIdx],
                                  oilWaterParams[cartElemIdx],
                                  gasOilParams[cartElemIdx]);

            materialLawParams_[cartElemIdx]->finalize();
        }
    }

    template <class Container>
    void readGasOilEffectiveParameters_(Container& dest,
                                        Opm::EclipseStateConstPtr eclState,
                                        int satRegionIdx)
    {
        dest[satRegionIdx] = std::make_shared<GasOilEffectiveTwoPhaseParams>();

        auto& effParams = *dest[satRegionIdx];
        const auto& sgofTable = eclState->getSgofTables()[satRegionIdx];

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
                                          int satRegionIdx)
    {
        dest[satRegionIdx] = std::make_shared<OilWaterEffectiveTwoPhaseParams>();

        auto& effParams = *dest[satRegionIdx];
        const auto& swofTable = eclState->getSwofTables()[satRegionIdx];

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
                                   int satRegionIdx)
    {
        dest[satRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satRegionIdx]->init(unscaledEpsInfo_[satRegionIdx], *config, EclGasOilSystem);
    }

    template <class Container>
    void readOilWaterUnscaledPoints_(Container &dest,
                                     std::shared_ptr<EclEpsConfig> config,
                                     Opm::DeckConstPtr deck,
                                     Opm::EclipseStateConstPtr eclState,
                                     int satRegionIdx)
    {
        dest[satRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satRegionIdx]->init(unscaledEpsInfo_[satRegionIdx], *config, EclOilWaterSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readGasOilScaledPoints_(InfoContainer& destInfo,
                                 PointsContainer& destPoints,
                                 std::shared_ptr<EclEpsConfig> config,
                                 const EclEpsGridProperties& epsGridProperties,
                                 int cartElemIdx)
    {
        int satRegionIdx = (*epsGridProperties.satnum)[cartElemIdx] - 1;

        destInfo[cartElemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satRegionIdx]);
        destInfo[cartElemIdx]->extractScaled(epsGridProperties, cartElemIdx);

        destPoints[cartElemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[cartElemIdx]->init(*destInfo[cartElemIdx], *config, EclGasOilSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readOilWaterScaledPoints_(InfoContainer& destInfo,
                                   PointsContainer& destPoints,
                                   std::shared_ptr<EclEpsConfig> config,
                                   const EclEpsGridProperties& epsGridProperties,
                                   int cartElemIdx)
    {
        int satRegionIdx = (*epsGridProperties.satnum)[cartElemIdx] - 1;

        destInfo[cartElemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satRegionIdx]);;
        destInfo[cartElemIdx]->extractScaled(epsGridProperties, cartElemIdx);

        destPoints[cartElemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[cartElemIdx]->init(*destInfo[cartElemIdx], *config, EclOilWaterSystem);
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

    std::vector<int> satnumData_;
};
} // namespace Opm

#endif
