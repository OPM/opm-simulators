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
 * \copydoc Opm::EclOutputBlackOilModule
 */
#ifndef EWOMS_ECL_GENERIC_OUTPUT_BLACK_OIL_MODULE_HH
#define EWOMS_ECL_GENERIC_OUTPUT_BLACK_OIL_MODULE_HH

#include <array>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>

#include <opm/output/data/Wells.hpp>
#include <opm/output/eclipse/Inplace.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>

namespace Opm {

namespace data { struct Solution; }
class EclipseState;

template<class FluidSystem, class Scalar>
class EclGenericOutputBlackoilModule {
public:
    using Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>;

    // write cumulative production and injection reports to output
    void outputCumLog(size_t reportStepNum,
                      const bool substep,
                      bool forceDisableCumOutput);

    // write production report to output
    void outputProdLog(size_t reportStepNum,
                       const bool substep,
                       bool forceDisableProdOutput);

    // write injection report to output
    void outputInjLog(size_t reportStepNum,
                      const bool substep,
                      bool forceDisableInjOutput);

    // write Fluid In Place to output log
    Inplace outputFipLog(std::map<std::string, double>& miscSummaryData,
                         std::map<std::string, std::vector<double>>& regionData,
                         const bool substep,
                         const Comm& comm);


    void outputErrorLog(const Comm& comm) const;

    void addRftDataToWells(data::Wells& wellDatas,
                           size_t reportStepNum);

    /*!
     * \brief Move all buffers to data::Solution.
     */
    void assignToSolution(data::Solution& sol);

    void setRestart(const data::Solution& sol,
                    unsigned elemIdx,
                    unsigned globalDofIndex);

    Scalar getSolventSaturation(unsigned elemIdx) const
    {
        if (sSol_.size() > elemIdx)
            return sSol_[elemIdx];

        return 0;
    }

    Scalar getPolymerConcentration(unsigned elemIdx) const
    {
        if (cPolymer_.size() > elemIdx)
            return cPolymer_[elemIdx];

        return 0;
    }

    Scalar getPolymerMW(unsigned elemIdx) const
    {
        if (mwPolymer_.size() > elemIdx)
            return mwPolymer_[elemIdx];

        return 0;
    }

    Scalar getFoamConcentration(unsigned elemIdx) const
    {
        if (cFoam_.size() > elemIdx)
            return cFoam_[elemIdx];

        return 0;
    }

    Scalar getSaltConcentration(unsigned elemIdx) const
    {
        if (cSalt_.size() > elemIdx)
            return cSalt_[elemIdx];

        return 0;
    }

    const std::map<std::size_t, double>& getWBPData() const
    {
        return this->wbpData_;
    }

    const std::map<std::pair<std::string, int>, double>& getBlockData()
    {
        return blockData_;
    }

    const Inplace& initialInplace() const
    {
        return this->initialInplace_.value();
    }

protected:
    using ScalarBuffer = std::vector<Scalar>;
    using StringBuffer = std::vector<std::string>;
    enum { numPhases = FluidSystem::numPhases };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    EclGenericOutputBlackoilModule(const EclipseState& eclState,
                                   const Schedule& schedule,
                                   const SummaryConfig& summaryConfig,
                                   const SummaryState& summaryState,
                                   bool enableEnergy,
                                   bool enableSolvent,
                                   bool enablePolymer,
                                   bool enablePolymerMolarWeight,
                                   bool enableFoam,
                                   bool enableBrine,
                                   bool enableExtbo);

    struct WellProdDataType
    {
        enum WPId
        {
            WellLocationi = 0, //WLi
            WellLocationj = 1, //WLj
            OilRate = 2, //OR
            WaterRate = 3, //WR
            GasRate = 4, //GR
            FluidResVol = 5, //FRV
            WaterCut = 6, //WC
            GasOilRatio = 7, //GOR
            WatGasRatio = 8, //WGR
            BHP = 9, //BHP
            THP = 10, //THP
            SteadyStatePI = 11, //SteadyStatePI
            WellName = 0, //WName
            CTRLMode = 1, //CTRL
        };

        static constexpr int numWPValues = 12;
        static constexpr int numWPNames = 2;
    };

    struct WellInjDataType
    {
        enum WIId
        {
            WellLocationi = 0, //WLi
            WellLocationj = 1, //WLj
            OilRate = 2, //OR
            WaterRate = 3, //WR
            GasRate = 4, //GR
            FluidResVol = 5, //FRV
            BHP = 6, //BHP
            THP = 7, //THP
            SteadyStateII = 8, //SteadyStateII
            WellName = 0, //WName
            CTRLModeOil = 1, //CTRLo
            CTRLModeWat = 2, //CTRLw
            CTRLModeGas = 3, //CTRLg
        };
        static constexpr int numWIValues = 9;
        static constexpr int numWINames = 4;
    };

    struct WellCumDataType
    {
        enum WCId
        {
            WellLocationi = 0, //WLi
            WellLocationj = 1, //WLj
            OilProd = 2, //OP
            WaterProd = 3, //WP
            GasProd = 4, //GP
            FluidResVolProd = 5, //FRVP
            OilInj = 6, //OI
            WaterInj = 7, //WI
            GasInj = 8, //GI
            FluidResVolInj = 9, //FRVI
            WellName = 0, //WName
            WellType = 1, //WType
            WellCTRL = 2, //WCTRL
        };
        static constexpr int numWCValues = 10;
        static constexpr int numWCNames = 3;
    };

    void doAllocBuffers(unsigned bufferSize,
                        unsigned reportStepNum,
                        const bool substep,
                        const bool log,
                        const bool isRestart,
                        const bool vapparsActive,
                        const bool enableHysteresis,
                        unsigned numTracers);

    void fipUnitConvert_(std::unordered_map<Inplace::Phase, Scalar>& fip) const;

    void pressureUnitConvert_(Scalar& pav) const;

    void outputRegionFluidInPlace_(std::unordered_map<Inplace::Phase, Scalar> oip,
                                   std::unordered_map<Inplace::Phase, Scalar> cip,
                                   const Scalar& pav, const int reg = 0) const;
    void outputProductionReport_(const ScalarBuffer& wellProd,
                                 const StringBuffer& wellProdNames,
                                 const bool forceDisableProdOutput);
    void outputInjectionReport_(const ScalarBuffer& wellInj,
                                const StringBuffer& wellInjNames,
                                const bool forceDisableInjOutput);
    void outputCumulativeReport_(const ScalarBuffer& wellCum,
                                 const StringBuffer& wellCumNames,
                                 const bool forceDisableCumOutput);

    void outputFipLogImpl(const Inplace& inplace) const;

    void makeRegionSum(Inplace& inplace,
                       const std::string& region_name,
                       const Comm& comm);

    Inplace accumulateRegionSums(const Comm& comm);

    void updateSummaryRegionValues(const Inplace& inplace,
                                   std::map<std::string, double>& miscSummaryData,
                                   std::map<std::string, std::vector<double>>& regionData) const;

    static bool isOutputCreationDirective_(const std::string& keyword);

    static Scalar pressureAverage_(const Scalar& pressurePvHydrocarbon,
                                   const Scalar& pvHydrocarbon,
                                   const Scalar& pressurePv,
                                   const Scalar& pv,
                                   bool hydrocarbon);

    static ScalarBuffer pressureAverage_(const ScalarBuffer& pressurePvHydrocarbon,
                                         const ScalarBuffer& pvHydrocarbon,
                                         const ScalarBuffer& pressurePv,
                                         const ScalarBuffer& pv,
                                         bool hydrocarbon);
    // Sum Fip values over regions.
    static ScalarBuffer regionSum(const ScalarBuffer& property,
                                  const std::vector<int>& regionId,
                                  size_t maxNumberOfRegions,
                                  const Comm& comm);

    static int regionMax(const std::vector<int>& region,
                         const Comm& comm);

    static void update(Inplace& inplace,
                       const std::string& region_name,
                       Inplace::Phase phase,
                       std::size_t ntFip,
                       const std::vector<double>& values);

    static Scalar sum(const ScalarBuffer& v);

    virtual bool isDefunctParallelWell(std::string wname) const = 0;

    const EclipseState& eclState_;
    const Schedule& schedule_;
    const SummaryConfig& summaryConfig_;
    const SummaryState& summaryState_;
    bool enableEnergy_;

    bool enableSolvent_;
    bool enablePolymer_;
    bool enablePolymerMolarWeight_;
    bool enableFoam_;
    bool enableBrine_;
    bool enableExtbo_;

    bool forceDisableFipOutput_;
    bool outputFipRestart_;
    bool computeFip_;

    std::unordered_map<Inplace::Phase, ScalarBuffer> fip_;
    std::unordered_map<std::string, std::vector<int>> regions_;
    std::unordered_map<Inplace::Phase, std::vector<SummaryConfigNode>> regionNodes_;

    std::vector<SummaryConfigNode> RPRNodes_;
    std::vector<SummaryConfigNode> RPRPNodes_;

    std::vector<int> failedCellsPb_;
    std::vector<int> failedCellsPd_;

    ScalarBuffer gasFormationVolumeFactor_;
    ScalarBuffer hydrocarbonPoreVolume_;
    ScalarBuffer pressureTimesPoreVolume_;
    ScalarBuffer pressureTimesHydrocarbonVolume_;
    ScalarBuffer oilPressure_;
    ScalarBuffer temperature_;
    ScalarBuffer rs_;
    ScalarBuffer rv_;
    ScalarBuffer overburdenPressure_;
    ScalarBuffer oilSaturationPressure_;
    ScalarBuffer sSol_;
    ScalarBuffer cPolymer_;
    ScalarBuffer mwPolymer_;
    ScalarBuffer cFoam_;
    ScalarBuffer cSalt_;
    ScalarBuffer extboX_;
    ScalarBuffer extboY_;
    ScalarBuffer extboZ_;
    ScalarBuffer mFracOil_;
    ScalarBuffer mFracGas_;
    ScalarBuffer mFracCo2_;
    ScalarBuffer soMax_;
    ScalarBuffer pcSwMdcOw_;
    ScalarBuffer krnSwMdcOw_;
    ScalarBuffer pcSwMdcGo_;
    ScalarBuffer krnSwMdcGo_;
    ScalarBuffer ppcw_;
    ScalarBuffer gasDissolutionFactor_;
    ScalarBuffer oilVaporizationFactor_;
    ScalarBuffer bubblePointPressure_;
    ScalarBuffer dewPointPressure_;
    ScalarBuffer rockCompPorvMultiplier_;
    ScalarBuffer swMax_;
    ScalarBuffer minimumOilPressure_;
    ScalarBuffer saturatedOilFormationVolumeFactor_;
    ScalarBuffer rockCompTransMultiplier_;

    std::array<ScalarBuffer, numPhases> saturation_;
    std::array<ScalarBuffer, numPhases> invB_;
    std::array<ScalarBuffer, numPhases> density_;
    std::array<ScalarBuffer, numPhases> viscosity_;
    std::array<ScalarBuffer, numPhases> relativePermeability_;

    std::vector<ScalarBuffer> tracerConcentrations_;

    std::map<size_t, Scalar> oilConnectionPressures_;
    std::map<size_t, Scalar> waterConnectionSaturations_;
    std::map<size_t, Scalar> gasConnectionSaturations_;
    std::map<std::pair<std::string, int>, double> blockData_;
    std::map<std::size_t , double> wbpData_;

    std::optional<Inplace> initialInplace_;
};

} // namespace Opm

#endif
