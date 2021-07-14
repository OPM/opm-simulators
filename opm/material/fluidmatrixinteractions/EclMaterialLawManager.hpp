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

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableColumn.hpp>

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <vector>

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
    typedef TwoPhaseMaterialTraits<Scalar, waterPhaseIdx, gasPhaseIdx> GasWaterTraits;

    // the two-phase material law which is defined on effective (unscaled) saturations
    typedef PiecewiseLinearTwoPhaseMaterial<GasOilTraits> GasOilEffectiveTwoPhaseLaw;
    typedef PiecewiseLinearTwoPhaseMaterial<OilWaterTraits> OilWaterEffectiveTwoPhaseLaw;
    typedef PiecewiseLinearTwoPhaseMaterial<GasWaterTraits> GasWaterEffectiveTwoPhaseLaw;

    typedef typename GasOilEffectiveTwoPhaseLaw::Params GasOilEffectiveTwoPhaseParams;
    typedef typename OilWaterEffectiveTwoPhaseLaw::Params OilWaterEffectiveTwoPhaseParams;
    typedef typename GasWaterEffectiveTwoPhaseLaw::Params GasWaterEffectiveTwoPhaseParams;

    // the two-phase material law which is defined on absolute (scaled) saturations
    typedef EclEpsTwoPhaseLaw<GasOilEffectiveTwoPhaseLaw> GasOilEpsTwoPhaseLaw;
    typedef EclEpsTwoPhaseLaw<OilWaterEffectiveTwoPhaseLaw> OilWaterEpsTwoPhaseLaw;
     typedef EclEpsTwoPhaseLaw<GasWaterEffectiveTwoPhaseLaw> GasWaterEpsTwoPhaseLaw;
    typedef typename GasOilEpsTwoPhaseLaw::Params GasOilEpsTwoPhaseParams;
    typedef typename OilWaterEpsTwoPhaseLaw::Params OilWaterEpsTwoPhaseParams;
    typedef typename GasWaterEpsTwoPhaseLaw::Params GasWaterEpsTwoPhaseParams;

    // the scaled two-phase material laws with hystersis
    typedef EclHysteresisTwoPhaseLaw<GasOilEpsTwoPhaseLaw> GasOilTwoPhaseLaw;
    typedef EclHysteresisTwoPhaseLaw<OilWaterEpsTwoPhaseLaw> OilWaterTwoPhaseLaw;
    typedef EclHysteresisTwoPhaseLaw<GasWaterEpsTwoPhaseLaw> GasWaterTwoPhaseLaw;
    typedef typename GasOilTwoPhaseLaw::Params GasOilTwoPhaseHystParams;
    typedef typename OilWaterTwoPhaseLaw::Params OilWaterTwoPhaseHystParams;
    typedef typename GasWaterTwoPhaseLaw::Params GasWaterTwoPhaseHystParams;

public:
    // the three-phase material law used by the simulation
    typedef EclMultiplexerMaterial<Traits, GasOilTwoPhaseLaw, OilWaterTwoPhaseLaw, GasWaterTwoPhaseLaw> MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

private:
    // internal typedefs
    typedef std::vector<std::shared_ptr<GasOilEffectiveTwoPhaseParams> > GasOilEffectiveParamVector;
    typedef std::vector<std::shared_ptr<OilWaterEffectiveTwoPhaseParams> > OilWaterEffectiveParamVector;
    typedef std::vector<std::shared_ptr<GasWaterEffectiveTwoPhaseParams> > GasWaterEffectiveParamVector;

    typedef std::vector<std::shared_ptr<EclEpsScalingPoints<Scalar> > > GasOilScalingPointsVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPoints<Scalar> > > OilWaterScalingPointsVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPoints<Scalar> > > GasWaterScalingPointsVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPointsInfo<Scalar> > > GasOilScalingInfoVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPointsInfo<Scalar> > > OilWaterScalingInfoVector;
    typedef std::vector<std::shared_ptr<EclEpsScalingPointsInfo<Scalar> > > GasWaterScalingInfoVector;
    typedef std::vector<std::shared_ptr<GasOilTwoPhaseHystParams> > GasOilParamVector;
    typedef std::vector<std::shared_ptr<OilWaterTwoPhaseHystParams> > OilWaterParamVector;
    typedef std::vector<std::shared_ptr<GasWaterTwoPhaseHystParams> > GasWaterParamVector;
    typedef std::vector<std::shared_ptr<MaterialLawParams> > MaterialLawParamsVector;

public:
    EclMaterialLawManager()
    {}

    void initFromState(const EclipseState& eclState)
    {
        // get the number of saturation regions and the number of cells in the deck
        const auto&  runspec       = eclState.runspec();
        const size_t numSatRegions = runspec.tabdims().getNumSatTables();

        const auto& ph = runspec.phases();
        this->hasGas = ph.active(Phase::GAS);
        this->hasOil = ph.active(Phase::OIL);
        this->hasWater = ph.active(Phase::WATER);

        readGlobalEpsOptions_(eclState);
        readGlobalHysteresisOptions_(eclState);
        readGlobalThreePhaseOptions_(runspec);

        // Read the end point scaling configuration (once per run).
        gasOilConfig = std::make_shared<EclEpsConfig>();
        oilWaterConfig = std::make_shared<EclEpsConfig>();
        gasWaterConfig = std::make_shared<EclEpsConfig>();
        gasOilConfig->initFromState(eclState, EclGasOilSystem);
        oilWaterConfig->initFromState(eclState, EclOilWaterSystem);
        gasWaterConfig->initFromState(eclState, EclGasWaterSystem);


        const auto& tables = eclState.getTableManager();

        {
            const auto& stone1exTables = tables.getStone1exTable();

            if (! stone1exTables.empty()) {
                stoneEtas.clear();
                stoneEtas.reserve(numSatRegions);

                for (const auto& table : stone1exTables) {
                    stoneEtas.push_back(table.eta);
                }
            }
        }

        this->unscaledEpsInfo_.resize(numSatRegions);

        if (this->hasGas + this->hasOil + this->hasWater == 1) {
            // Single-phase simulation.  Special case.  Nothing to do here.
            return;
        }

        // Multiphase simulation.  Common case.
        const auto tolcrit = runspec.saturationFunctionControls()
            .minimumRelpermMobilityThreshold();

        const auto rtep  = satfunc::getRawTableEndpoints(tables, ph, tolcrit);
        const auto rfunc = satfunc::getRawFunctionValues(tables, ph, rtep);

        for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
            this->unscaledEpsInfo_[satRegionIdx]
                .extractUnscaled(rtep, rfunc, satRegionIdx);
        }
    }

    void initParamsForElements(const EclipseState& eclState, size_t numCompressedElems)
    {
        // get the number of saturation regions
        const size_t numSatRegions = eclState.runspec().tabdims().getNumSatTables();

        // setup the saturation region specific parameters
        gasOilUnscaledPointsVector_.resize(numSatRegions);
        oilWaterUnscaledPointsVector_.resize(numSatRegions);
        gasWaterUnscaledPointsVector_.resize(numSatRegions);

        gasOilEffectiveParamVector_.resize(numSatRegions);
        oilWaterEffectiveParamVector_.resize(numSatRegions);
        gasWaterEffectiveParamVector_.resize(numSatRegions);
        for (unsigned satRegionIdx = 0; satRegionIdx < numSatRegions; ++satRegionIdx) {
            // unscaled points for end-point scaling
            readGasOilUnscaledPoints_(gasOilUnscaledPointsVector_, gasOilConfig, eclState, satRegionIdx);
            readOilWaterUnscaledPoints_(oilWaterUnscaledPointsVector_, oilWaterConfig, eclState, satRegionIdx);
            readGasWaterUnscaledPoints_(gasWaterUnscaledPointsVector_, gasWaterConfig, eclState, satRegionIdx);

            // the parameters for the effective two-phase matererial laws
            readGasOilEffectiveParameters_(gasOilEffectiveParamVector_, eclState, satRegionIdx);
            readOilWaterEffectiveParameters_(oilWaterEffectiveParamVector_, eclState, satRegionIdx);
            readGasWaterEffectiveParameters_(gasWaterEffectiveParamVector_, eclState, satRegionIdx);
        }

        // copy the SATNUM grid property. in some cases this is not necessary, but it
        // should not require much memory anyway...
        satnumRegionArray_.resize(numCompressedElems);
        if (eclState.fieldProps().has_int("SATNUM")) {
            const auto& satnumRawData = eclState.fieldProps().get_int("SATNUM");
            for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
                satnumRegionArray_[elemIdx] = satnumRawData[elemIdx] - 1;
            }
        }
        else
            std::fill(satnumRegionArray_.begin(), satnumRegionArray_.end(), 0);

        // create the information for the imbibition region (IMBNUM). By default this is
        // the same as the saturation region (SATNUM)
        imbnumRegionArray_ = satnumRegionArray_;
        if (eclState.fieldProps().has_int("IMBNUM")) {
            const auto& imbnumRawData = eclState.fieldProps().get_int("IMBNUM");
            for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
                imbnumRegionArray_[elemIdx] = imbnumRawData[elemIdx] - 1;
            }
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
        
        GasWaterScalingInfoVector gasWaterScaledInfoVector(numCompressedElems);
        GasWaterScalingPointsVector gasWaterScaledPointsVector(numCompressedElems);
        GasWaterScalingInfoVector gasWaterScaledImbInfoVector;
        GasWaterScalingPointsVector gasWaterScaledImbPointsVector;
     
        if (enableHysteresis()) {
            gasOilScaledImbInfoVector.resize(numCompressedElems);
            gasOilScaledImbPointsVector.resize(numCompressedElems);
            oilWaterScaledImbInfoVector.resize(numCompressedElems);
            oilWaterScaledImbPointsVector.resize(numCompressedElems);
            gasWaterScaledImbInfoVector.resize(numCompressedElems);
            gasWaterScaledImbPointsVector.resize(numCompressedElems);
        }

        EclEpsGridProperties epsGridProperties(eclState, false);

        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            readGasOilScaledPoints_(gasOilScaledInfoVector,
                                    gasOilScaledPointsVector,
                                    gasOilConfig,
                                    eclState,
                                    epsGridProperties,
                                    elemIdx);

            readOilWaterScaledPoints_(oilWaterScaledEpsInfoDrainage_,
                                      oilWaterScaledEpsPointsDrainage,
                                      oilWaterConfig,
                                      eclState,
                                      epsGridProperties,
                                      elemIdx);
                                      
            readGasWaterScaledPoints_(gasWaterScaledInfoVector,
                                      gasWaterScaledPointsVector,
                                      gasWaterConfig,
                                      eclState,
                                      epsGridProperties,
                                      elemIdx);

        }

        if (enableHysteresis()) {
            EclEpsGridProperties epsImbGridProperties(eclState, true);
            for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
                readGasOilScaledPoints_(gasOilScaledImbInfoVector,
                                        gasOilScaledImbPointsVector,
                                        gasOilConfig,
                                        eclState,
                                        epsImbGridProperties,
                                        elemIdx);

                readOilWaterScaledPoints_(oilWaterScaledImbInfoVector,
                                          oilWaterScaledImbPointsVector,
                                          oilWaterConfig,
                                          eclState,
                                          epsImbGridProperties,
                                          elemIdx);
                
                readGasWaterScaledPoints_(gasWaterScaledImbInfoVector,
                                          gasWaterScaledImbPointsVector,
                                          gasWaterConfig,
                                          eclState,
                                          epsImbGridProperties,
                                          elemIdx);
            }
        }

        // create the parameter objects for the two-phase laws
        GasOilParamVector gasOilParams(numCompressedElems);
        OilWaterParamVector oilWaterParams(numCompressedElems);
        GasWaterParamVector gasWaterParams(numCompressedElems);

        GasOilParamVector gasOilImbParams;
        OilWaterParamVector oilWaterImbParams;
        GasWaterParamVector gasWaterImbParams;

        if (enableHysteresis()) {
            gasOilImbParams.resize(numCompressedElems);
            oilWaterImbParams.resize(numCompressedElems);
            gasWaterImbParams.resize(numCompressedElems);
        }

        assert(numCompressedElems == satnumRegionArray_.size());
        assert(!enableHysteresis() || numCompressedElems == imbnumRegionArray_.size());
        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            unsigned satRegionIdx = static_cast<unsigned>(satnumRegionArray_[elemIdx]);

            gasOilParams[elemIdx] = std::make_shared<GasOilTwoPhaseHystParams>();
            oilWaterParams[elemIdx] = std::make_shared<OilWaterTwoPhaseHystParams>();
            gasWaterParams[elemIdx] = std::make_shared<GasWaterTwoPhaseHystParams>();
            gasOilParams[elemIdx]->setConfig(hysteresisConfig_);
            oilWaterParams[elemIdx]->setConfig(hysteresisConfig_);
            gasWaterParams[elemIdx]->setConfig(hysteresisConfig_);

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

            if (hasGas && hasWater && !hasOil) {
                auto gasWaterDrainParams = std::make_shared<GasWaterEpsTwoPhaseParams>();
                gasWaterDrainParams->setConfig(gasWaterConfig);
                gasWaterDrainParams->setUnscaledPoints(gasWaterUnscaledPointsVector_[satRegionIdx]);
                gasWaterDrainParams->setScaledPoints(gasWaterScaledPointsVector[elemIdx]);
                gasWaterDrainParams->setEffectiveLawParams(gasWaterEffectiveParamVector_[satRegionIdx]);
                gasWaterDrainParams->finalize();

                gasWaterParams[elemIdx]->setDrainageParams(gasWaterDrainParams,
                                                           *gasWaterScaledInfoVector[elemIdx],
                                                           EclGasWaterSystem);
            }


            if (enableHysteresis()) {
                unsigned imbRegionIdx = imbnumRegionArray_[elemIdx];

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

                if (hasGas && hasWater && !hasOil) {
                    auto gasWaterImbParamsHyst = std::make_shared<GasWaterEpsTwoPhaseParams>();
                    gasWaterImbParamsHyst->setConfig(gasWaterConfig);
                    gasWaterImbParamsHyst->setUnscaledPoints(gasWaterUnscaledPointsVector_[imbRegionIdx]);
                    gasWaterImbParamsHyst->setScaledPoints(gasWaterScaledImbPointsVector[elemIdx]);
                    gasWaterImbParamsHyst->setEffectiveLawParams(gasWaterEffectiveParamVector_[imbRegionIdx]);
                    gasWaterImbParamsHyst->finalize();

                    gasWaterParams[elemIdx]->setImbibitionParams(gasWaterImbParamsHyst,
                                                                 *gasWaterScaledImbInfoVector[elemIdx],
                                                                 EclGasWaterSystem);
                }
            }

            if (hasGas && hasOil)
                gasOilParams[elemIdx]->finalize();

            if (hasOil && hasWater)
                oilWaterParams[elemIdx]->finalize();

            if (hasGas && hasWater && !hasOil)
                gasWaterParams[elemIdx]->finalize();
        }

        // create the parameter objects for the three-phase law
        materialLawParams_.resize(numCompressedElems);
        for (unsigned elemIdx = 0; elemIdx < numCompressedElems; ++elemIdx) {
            materialLawParams_[elemIdx] = std::make_shared<MaterialLawParams>();
            unsigned satRegionIdx = static_cast<unsigned>(satnumRegionArray_[elemIdx]);

            initThreePhaseParams_(eclState,
                                  *materialLawParams_[elemIdx],
                                  satRegionIdx,
                                  *oilWaterScaledEpsInfoDrainage_[elemIdx],
                                  oilWaterParams[elemIdx],
                                  gasOilParams[elemIdx],
                                  gasWaterParams[elemIdx]);

            materialLawParams_[elemIdx]->finalize();
        }
    }


    /*!
     * \brief Modify the initial condition according to the SWATINIT keyword.
     *
     * The method returns the water saturation which yields a givenn capillary
     * pressure. The reason this method is not folded directly into initFromState() is
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
            typedef SimpleModularFluidState<Scalar,
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
            const Scalar pcowAtSwThreshold = 1.0; //Pascal
            // avoid divison by very small number
            if (std::abs(pcowAtSw) > pcowAtSwThreshold) {
                elemScaledEpsInfo.maxPcow *= pcow/pcowAtSw;
                auto& elemEclEpsScalingPoints = oilWaterScaledEpsPointsDrainage(elemIdx);
                elemEclEpsScalingPoints.init(elemScaledEpsInfo, *oilWaterEclEpsConfig_, EclOilWaterSystem);
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
        assert(elemIdx <  materialLawParams_.size());
        return *materialLawParams_[elemIdx];
    }

    const MaterialLawParams& materialLawParams(unsigned elemIdx) const
    {
        assert(elemIdx <  materialLawParams_.size());
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
        case EclMultiplexerApproach::EclStone1Approach: {
            auto& realParams = mlp.template getRealParams<EclMultiplexerApproach::EclStone1Approach>();

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

        case EclMultiplexerApproach::EclStone2Approach: {
            auto& realParams = mlp.template getRealParams<EclMultiplexerApproach::EclStone2Approach>();
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

        case EclMultiplexerApproach::EclDefaultApproach: {
            auto& realParams = mlp.template getRealParams<EclMultiplexerApproach::EclDefaultApproach>();
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

        case EclMultiplexerApproach::EclTwoPhaseApproach: {
            auto& realParams = mlp.template getRealParams<EclMultiplexerApproach::EclTwoPhaseApproach>();
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

    int satnumRegionIdx(unsigned elemIdx) const
    { return satnumRegionArray_[elemIdx]; }

    int imbnumRegionIdx(unsigned elemIdx) const
    { return imbnumRegionArray_[elemIdx]; }

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
        if (!enableHysteresis())
            throw std::runtime_error("Cannot get hysteresis parameters if hysteresis not enabled.");

        const auto& params = materialLawParams(elemIdx);
        MaterialLaw::oilWaterHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    void setOilWaterHysteresisParams(const Scalar& pcSwMdc,
                                     const Scalar& krnSwMdc,
                                     unsigned elemIdx)
    {
        if (!enableHysteresis())
            throw std::runtime_error("Cannot set hysteresis parameters if hysteresis not enabled.");

        auto& params = materialLawParams(elemIdx);
        MaterialLaw::setOilWaterHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    void gasOilHysteresisParams(Scalar& pcSwMdc,
                                Scalar& krnSwMdc,
                                unsigned elemIdx) const
    {
        if (!enableHysteresis())
            throw std::runtime_error("Cannot get hysteresis parameters if hysteresis not enabled.");

        const auto& params = materialLawParams(elemIdx);
        MaterialLaw::gasOilHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    void setGasOilHysteresisParams(const Scalar& pcSwMdc,
                                   const Scalar& krnSwMdc,
                                   unsigned elemIdx)
    {
        if (!enableHysteresis())
            throw std::runtime_error("Cannot set hysteresis parameters if hysteresis not enabled.");

        auto& params = materialLawParams(elemIdx);
        MaterialLaw::setGasOilHysteresisParams(pcSwMdc, krnSwMdc, params);
    }

    EclEpsScalingPoints<Scalar>& oilWaterScaledEpsPointsDrainage(unsigned elemIdx)
    {
        auto& materialParams = *materialLawParams_[elemIdx];
        switch (materialParams.approach()) {
        case EclMultiplexerApproach::EclStone1Approach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclStone1Approach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclMultiplexerApproach::EclStone2Approach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclStone2Approach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclMultiplexerApproach::EclDefaultApproach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclDefaultApproach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }

        case EclMultiplexerApproach::EclTwoPhaseApproach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclTwoPhaseApproach>();
            return realParams.oilWaterParams().drainageParams().scaledPoints();
        }
        default:
            throw std::logic_error("Enum value for material approach unknown!");
        }
    }

    const EclEpsScalingPointsInfo<Scalar>& oilWaterScaledEpsInfoDrainage(size_t elemIdx) const
    { return *oilWaterScaledEpsInfoDrainage_[elemIdx]; }

    std::shared_ptr<EclEpsScalingPointsInfo<Scalar> >& oilWaterScaledEpsInfoDrainagePointerReferenceHack(unsigned elemIdx)
    { return oilWaterScaledEpsInfoDrainage_[elemIdx]; }

private:
    void readGlobalEpsOptions_(const EclipseState& eclState)
    {
        oilWaterEclEpsConfig_ = std::make_shared<EclEpsConfig>();
        oilWaterEclEpsConfig_->initFromState(eclState, EclOilWaterSystem);

        enableEndPointScaling_ = eclState.getTableManager().hasTables("ENKRVD");
    }

    void readGlobalHysteresisOptions_(const EclipseState& state)
    {
        hysteresisConfig_ = std::make_shared<EclHysteresisConfig>();
        hysteresisConfig_->initFromState(state.runspec());
    }

    void readGlobalThreePhaseOptions_(const Runspec& runspec)
    {
        bool gasEnabled = runspec.phases().active(Phase::GAS);
        bool oilEnabled = runspec.phases().active(Phase::OIL);
        bool waterEnabled = runspec.phases().active(Phase::WATER);

        int numEnabled =
            (gasEnabled?1:0)
            + (oilEnabled?1:0)
            + (waterEnabled?1:0);

        if (numEnabled == 0) {
            throw std::runtime_error("At least one fluid phase must be enabled. (Is: "+std::to_string(numEnabled)+")");
        } else if (numEnabled == 1) {
            threePhaseApproach_ = EclMultiplexerApproach::EclOnePhaseApproach;
        } else if ( numEnabled == 2) {
            threePhaseApproach_ = EclMultiplexerApproach::EclTwoPhaseApproach;
            if (!gasEnabled)
                twoPhaseApproach_ = EclTwoPhaseApproach::EclTwoPhaseOilWater;
            else if (!oilEnabled)
                twoPhaseApproach_ = EclTwoPhaseApproach::EclTwoPhaseGasWater;
            else if (!waterEnabled)
                twoPhaseApproach_ = EclTwoPhaseApproach::EclTwoPhaseGasOil;
        }
        else {
            assert(numEnabled == 3);

            threePhaseApproach_ = EclMultiplexerApproach::EclDefaultApproach;
            const auto& satctrls = runspec.saturationFunctionControls();
            if (satctrls.krModel() == SatFuncControls::ThreePhaseOilKrModel::Stone2)
                threePhaseApproach_ = EclMultiplexerApproach::EclStone2Approach;
            else if (satctrls.krModel() == SatFuncControls::ThreePhaseOilKrModel::Stone1)
                threePhaseApproach_ = EclMultiplexerApproach::EclStone1Approach;
        }
    }

    template <class Container>
    void readGasOilEffectiveParameters_(Container& dest,
                                        const EclipseState& eclState,
                                        unsigned satRegionIdx)
    {
        if (!hasGas || !hasOil)
            // we don't read anything if either the gas or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<GasOilEffectiveTwoPhaseParams>();

        auto& effParams = *dest[satRegionIdx];

        // the situation for the gas phase is complicated that all saturations are
        // shifted by the connate water saturation.
        const Scalar Swco = unscaledEpsInfo_[satRegionIdx].Swl;
        const auto tolcrit = eclState.runspec().saturationFunctionControls()
            .minimumRelpermMobilityThreshold();

        const auto& tableManager = eclState.getTableManager();

        switch (eclState.runspec().saturationFunctionControls().family()) {
        case SatFuncControls::KeywordFamily::Family_I:
        {
            const TableContainer& sgofTables = tableManager.getSgofTables();
            const TableContainer& slgofTables = tableManager.getSlgofTables();
            if (!sgofTables.empty())
                readGasOilEffectiveParametersSgof_(effParams, Swco, tolcrit,
                                                   sgofTables.getTable<SgofTable>(satRegionIdx));
            else if (!slgofTables.empty())
                readGasOilEffectiveParametersSlgof_(effParams, Swco, tolcrit,
                                                    slgofTables.getTable<SlgofTable>(satRegionIdx));
            break;
        }

        case SatFuncControls::KeywordFamily::Family_II:
        {
            const SgfnTable& sgfnTable = tableManager.getSgfnTables().getTable<SgfnTable>( satRegionIdx );
            if (!hasWater) {
                // oil and gas case
                const Sof2Table& sof2Table = tableManager.getSof2Tables().getTable<Sof2Table>( satRegionIdx );
                readGasOilEffectiveParametersFamily2_(effParams, Swco, tolcrit, sof2Table, sgfnTable);
            }
            else {
                const Sof3Table& sof3Table = tableManager.getSof3Tables().getTable<Sof3Table>( satRegionIdx );
                readGasOilEffectiveParametersFamily2_(effParams, Swco, tolcrit, sof3Table, sgfnTable);
            }
            break;
        }

        case SatFuncControls::KeywordFamily::Undefined:
            throw std::domain_error("No valid saturation keyword family specified");
        }
    }

    void readGasOilEffectiveParametersSgof_(GasOilEffectiveTwoPhaseParams& effParams,
                                            const Scalar Swco,
                                            const double tolcrit,
                                            const SgofTable& sgofTable)
    {
        // convert the saturations of the SGOF keyword from gas to oil saturations
        std::vector<double> SoSamples(sgofTable.numRows());
        for (size_t sampleIdx = 0; sampleIdx < sgofTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = (1.0 - Swco) - sgofTable.get("SG", sampleIdx);
        }

        effParams.setKrwSamples(SoSamples, normalizeKrValues_(tolcrit, sgofTable.getColumn("KROG")));
        effParams.setKrnSamples(SoSamples, normalizeKrValues_(tolcrit, sgofTable.getColumn("KRG")));
        effParams.setPcnwSamples(SoSamples, sgofTable.getColumn("PCOG").vectorCopy());
        effParams.finalize();
    }

    void readGasOilEffectiveParametersSlgof_(GasOilEffectiveTwoPhaseParams& effParams,
                                             const Scalar Swco,
                                             const double tolcrit,
                                             const SlgofTable& slgofTable)
    {
        // convert the saturations of the SLGOF keyword from "liquid" to oil saturations
        std::vector<double> SoSamples(slgofTable.numRows());
        for (size_t sampleIdx = 0; sampleIdx < slgofTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = slgofTable.get("SL", sampleIdx) - Swco;
        }

        effParams.setKrwSamples(SoSamples, normalizeKrValues_(tolcrit, slgofTable.getColumn("KROG")));
        effParams.setKrnSamples(SoSamples, normalizeKrValues_(tolcrit, slgofTable.getColumn("KRG")));
        effParams.setPcnwSamples(SoSamples, slgofTable.getColumn("PCOG").vectorCopy());
        effParams.finalize();
    }

    void readGasOilEffectiveParametersFamily2_(GasOilEffectiveTwoPhaseParams& effParams,
                                               const Scalar Swco,
                                               const double tolcrit,
                                               const Sof3Table& sof3Table,
                                               const SgfnTable& sgfnTable)
    {
        // convert the saturations of the SGFN keyword from gas to oil saturations
        std::vector<double> SoSamples(sgfnTable.numRows());
        std::vector<double> SoColumn = sof3Table.getColumn("SO").vectorCopy();
        for (size_t sampleIdx = 0; sampleIdx < sgfnTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = (1.0 - Swco) - sgfnTable.get("SG", sampleIdx);
        }

        effParams.setKrwSamples(SoColumn, normalizeKrValues_(tolcrit, sof3Table.getColumn("KROG")));
        effParams.setKrnSamples(SoSamples, normalizeKrValues_(tolcrit, sgfnTable.getColumn("KRG")));
        effParams.setPcnwSamples(SoSamples, sgfnTable.getColumn("PCOG").vectorCopy());
        effParams.finalize();
    }

    void readGasOilEffectiveParametersFamily2_(GasOilEffectiveTwoPhaseParams& effParams,
                                               const Scalar Swco,
                                               const double tolcrit,
                                               const Sof2Table& sof2Table,
                                               const SgfnTable& sgfnTable)
    {
        // convert the saturations of the SGFN keyword from gas to oil saturations
        std::vector<double> SoSamples(sgfnTable.numRows());
        std::vector<double> SoColumn = sof2Table.getColumn("SO").vectorCopy();
        for (size_t sampleIdx = 0; sampleIdx < sgfnTable.numRows(); ++ sampleIdx) {
            SoSamples[sampleIdx] = (1.0 - Swco) - sgfnTable.get("SG", sampleIdx);
        }

        effParams.setKrwSamples(SoColumn, normalizeKrValues_(tolcrit, sof2Table.getColumn("KRO")));
        effParams.setKrnSamples(SoSamples, normalizeKrValues_(tolcrit, sgfnTable.getColumn("KRG")));
        effParams.setPcnwSamples(SoSamples, sgfnTable.getColumn("PCOG").vectorCopy());
        effParams.finalize();
    }

    template <class Container>
    void readOilWaterEffectiveParameters_(Container& dest,
                                          const EclipseState& eclState,
                                          unsigned satRegionIdx)
    {
        if (!hasOil || !hasWater)
            // we don't read anything if either the water or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<OilWaterEffectiveTwoPhaseParams>();

        const auto tolcrit = eclState.runspec().saturationFunctionControls()
            .minimumRelpermMobilityThreshold();

        const auto& tableManager = eclState.getTableManager();
        auto& effParams = *dest[satRegionIdx];

        switch (eclState.runspec().saturationFunctionControls().family()) {
        case SatFuncControls::KeywordFamily::Family_I:
        {
            const auto& swofTable = tableManager.getSwofTables().getTable<SwofTable>(satRegionIdx);
            const std::vector<double> SwColumn = swofTable.getColumn("SW").vectorCopy();

            effParams.setKrwSamples(SwColumn, normalizeKrValues_(tolcrit, swofTable.getColumn("KRW")));
            effParams.setKrnSamples(SwColumn, normalizeKrValues_(tolcrit, swofTable.getColumn("KROW")));
            effParams.setPcnwSamples(SwColumn, swofTable.getColumn("PCOW").vectorCopy());
            effParams.finalize();
            break;
        }

        case SatFuncControls::KeywordFamily::Family_II:
        {
            const auto& swfnTable = tableManager.getSwfnTables().getTable<SwfnTable>(satRegionIdx);
            const auto& sof3Table = tableManager.getSof3Tables().getTable<Sof3Table>(satRegionIdx);
            const std::vector<double> SwColumn = swfnTable.getColumn("SW").vectorCopy();

            // convert the saturations of the SOF3 keyword from oil to water saturations
            std::vector<double> SwSamples(sof3Table.numRows());
            for (size_t sampleIdx = 0; sampleIdx < sof3Table.numRows(); ++ sampleIdx)
                SwSamples[sampleIdx] = 1 - sof3Table.get("SO", sampleIdx);

            effParams.setKrwSamples(SwColumn, normalizeKrValues_(tolcrit, swfnTable.getColumn("KRW")));
            effParams.setKrnSamples(SwSamples, normalizeKrValues_(tolcrit, sof3Table.getColumn("KROW")));
            effParams.setPcnwSamples(SwColumn, swfnTable.getColumn("PCOW").vectorCopy());
            effParams.finalize();
            break;
        }

        case SatFuncControls::KeywordFamily::Undefined:
            throw std::domain_error("No valid saturation keyword family specified");
        }
    }

    template <class Container>
    void readGasWaterEffectiveParameters_(Container& dest,
                                        const EclipseState& eclState,
                                        unsigned satRegionIdx)
    {
        if (!hasGas || !hasWater || hasOil)
            // we don't read anything if either the gas or the water phase is not active or if oil is present
            return;

        dest[satRegionIdx] = std::make_shared<GasWaterEffectiveTwoPhaseParams>();

        auto& effParams = *dest[satRegionIdx];

        const auto tolcrit = eclState.runspec().saturationFunctionControls()
            .minimumRelpermMobilityThreshold();

        const auto& tableManager = eclState.getTableManager();

        switch (eclState.runspec().saturationFunctionControls().family()) {
        case SatFuncControls::KeywordFamily::Family_I:
        {
            throw std::domain_error("Saturation keyword family I is not applicable for a gas-water system");
        }

        case SatFuncControls::KeywordFamily::Family_II:
        {
            //Todo: allow also for Sgwfn table input as alternative to Sgfn and Swfn table input
            const SgfnTable& sgfnTable = tableManager.getSgfnTables().getTable<SgfnTable>( satRegionIdx );
            const SwfnTable& swfnTable = tableManager.getSwfnTables().getTable<SwfnTable>( satRegionIdx );

            std::vector<double> SwColumn = swfnTable.getColumn("SW").vectorCopy();
            
            effParams.setKrwSamples(SwColumn, normalizeKrValues_(tolcrit, swfnTable.getColumn("KRW")));
            std::vector<double> SwSamples(sgfnTable.numRows());
            for (size_t sampleIdx = 0; sampleIdx < sgfnTable.numRows(); ++ sampleIdx)
                SwSamples[sampleIdx] = 1 - sgfnTable.get("SG", sampleIdx);
            effParams.setKrnSamples(SwSamples, normalizeKrValues_(tolcrit, sgfnTable.getColumn("KRG")));
            //Capillary pressure is read from SWFN. 
            //For gas-water system the capillary pressure column values are set to 0 in SGFN
            effParams.setPcnwSamples(SwColumn, swfnTable.getColumn("PCOW").vectorCopy());
            effParams.finalize();
                       
            break;
        }

        case SatFuncControls::KeywordFamily::Undefined:
            throw std::domain_error("No valid saturation keyword family specified");
        }
    }

    template <class Container>
    void readGasOilUnscaledPoints_(Container& dest,
                                   std::shared_ptr<EclEpsConfig> config,
                                   const EclipseState& /* eclState */,
                                   unsigned satRegionIdx)
    {
        if (!hasGas || !hasOil)
            // we don't read anything if either the gas or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satRegionIdx]->init(unscaledEpsInfo_[satRegionIdx], *config, EclGasOilSystem);
    }

    template <class Container>
    void readOilWaterUnscaledPoints_(Container& dest,
                                     std::shared_ptr<EclEpsConfig> config,
                                     const EclipseState& /* eclState */,
                                     unsigned satRegionIdx)
    {
        if (!hasOil || !hasWater)
            // we don't read anything if either the water or the oil phase is not active
            return;

        dest[satRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satRegionIdx]->init(unscaledEpsInfo_[satRegionIdx], *config, EclOilWaterSystem);
    }

    template <class Container>
    void readGasWaterUnscaledPoints_(Container& dest,
                                     std::shared_ptr<EclEpsConfig> config,
                                     const EclipseState& /* eclState */,
                                     unsigned satRegionIdx)
    {
        if (hasOil)
            // we don't read anything if oil phase is active
            return;

        dest[satRegionIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        dest[satRegionIdx]->init(unscaledEpsInfo_[satRegionIdx], *config, EclGasWaterSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readGasOilScaledPoints_(InfoContainer& destInfo,
                                 PointsContainer& destPoints,
                                 std::shared_ptr<EclEpsConfig> config,
                                 const EclipseState& eclState,
                                 const EclEpsGridProperties& epsGridProperties,
                                 unsigned elemIdx)
    {
        unsigned satRegionIdx = epsGridProperties.satRegion( elemIdx );

        destInfo[elemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satRegionIdx]);
        destInfo[elemIdx]->extractScaled(eclState, epsGridProperties, elemIdx);

        destPoints[elemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[elemIdx]->init(*destInfo[elemIdx], *config, EclGasOilSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readOilWaterScaledPoints_(InfoContainer& destInfo,
                                   PointsContainer& destPoints,
                                   std::shared_ptr<EclEpsConfig> config,
                                   const EclipseState& eclState,
                                   const EclEpsGridProperties& epsGridProperties,
                                   unsigned elemIdx)
    {
        unsigned satRegionIdx = epsGridProperties.satRegion( elemIdx );

        destInfo[elemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satRegionIdx]);
        destInfo[elemIdx]->extractScaled(eclState, epsGridProperties, elemIdx);

        destPoints[elemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[elemIdx]->init(*destInfo[elemIdx], *config, EclOilWaterSystem);
    }

    template <class InfoContainer, class PointsContainer>
    void readGasWaterScaledPoints_(InfoContainer& destInfo,
                                   PointsContainer& destPoints,
                                   std::shared_ptr<EclEpsConfig> config,
                                   const EclipseState& eclState,
                                   const EclEpsGridProperties& epsGridProperties,
                                   unsigned elemIdx)
    {
        unsigned satRegionIdx = epsGridProperties.satRegion( elemIdx );

        destInfo[elemIdx] = std::make_shared<EclEpsScalingPointsInfo<Scalar> >(unscaledEpsInfo_[satRegionIdx]);
        destInfo[elemIdx]->extractScaled(eclState, epsGridProperties, elemIdx);

        destPoints[elemIdx] = std::make_shared<EclEpsScalingPoints<Scalar> >();
        destPoints[elemIdx]->init(*destInfo[elemIdx], *config, EclGasWaterSystem);
    }

    void initThreePhaseParams_(const EclipseState& /* eclState */,
                               MaterialLawParams& materialParams,
                               unsigned satRegionIdx,
                               const EclEpsScalingPointsInfo<Scalar>& epsInfo,
                               std::shared_ptr<OilWaterTwoPhaseHystParams> oilWaterParams,
                               std::shared_ptr<GasOilTwoPhaseHystParams> gasOilParams,
                               std::shared_ptr<GasWaterTwoPhaseHystParams> gasWaterParams)
    {
        materialParams.setApproach(threePhaseApproach_);

        switch (materialParams.approach()) {
        case EclMultiplexerApproach::EclStone1Approach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclStone1Approach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);

            if (!stoneEtas.empty()) {
                realParams.setEta(stoneEtas[satRegionIdx]);
            }
            else
                realParams.setEta(1.0);
            realParams.finalize();
            break;
        }

        case EclMultiplexerApproach::EclStone2Approach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclStone2Approach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);
            realParams.finalize();
            break;
        }

        case EclMultiplexerApproach::EclDefaultApproach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclDefaultApproach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setSwl(epsInfo.Swl);
            realParams.finalize();
            break;
        }

        case EclMultiplexerApproach::EclTwoPhaseApproach: {
            auto& realParams = materialParams.template getRealParams<EclMultiplexerApproach::EclTwoPhaseApproach>();
            realParams.setGasOilParams(gasOilParams);
            realParams.setOilWaterParams(oilWaterParams);
            realParams.setGasWaterParams(gasWaterParams);
            realParams.setApproach(twoPhaseApproach_);
            realParams.finalize();
            break;
        }

        case EclMultiplexerApproach::EclOnePhaseApproach: {
            // Nothing to do, no parameters.
            break;
        }
        }
    }

    // Relative permeability values not strictly greater than 'tolcrit' treated as zero.
    std::vector<double> normalizeKrValues_(const double tolcrit,
                                           const TableColumn& krValues) const
    {
        auto kr = krValues.vectorCopy();
        std::transform(kr.begin(), kr.end(), kr.begin(),
            [tolcrit](const double kri)
        {
            return (kri > tolcrit) ? kri : 0.0;
        });

        return kr;
    }

    bool enableEndPointScaling_;
    std::shared_ptr<EclHysteresisConfig> hysteresisConfig_;

    std::shared_ptr<EclEpsConfig> oilWaterEclEpsConfig_;
    std::vector<EclEpsScalingPointsInfo<Scalar>> unscaledEpsInfo_;
    OilWaterScalingInfoVector oilWaterScaledEpsInfoDrainage_;

    std::shared_ptr<EclEpsConfig> gasWaterEclEpsConfig_;
    GasWaterScalingInfoVector gasWaterScaledEpsInfoDrainage_;

    GasOilScalingPointsVector gasOilUnscaledPointsVector_;
    OilWaterScalingPointsVector oilWaterUnscaledPointsVector_;
    GasWaterScalingPointsVector gasWaterUnscaledPointsVector_;

    GasOilEffectiveParamVector gasOilEffectiveParamVector_;
    OilWaterEffectiveParamVector oilWaterEffectiveParamVector_;
    GasWaterEffectiveParamVector gasWaterEffectiveParamVector_;

    EclMultiplexerApproach threePhaseApproach_ = EclMultiplexerApproach::EclDefaultApproach;
    // this attribute only makes sense for twophase simulations!
    enum EclTwoPhaseApproach twoPhaseApproach_ = EclTwoPhaseApproach::EclTwoPhaseGasOil;

    std::vector<std::shared_ptr<MaterialLawParams> > materialLawParams_;

    std::vector<int> satnumRegionArray_;
    std::vector<int> imbnumRegionArray_;
    std::vector<Scalar> stoneEtas;

    bool hasGas;
    bool hasOil;
    bool hasWater;

    std::shared_ptr<EclEpsConfig> gasOilConfig;
    std::shared_ptr<EclEpsConfig> oilWaterConfig;
    std::shared_ptr<EclEpsConfig> gasWaterConfig;
};
} // namespace Opm

#endif
