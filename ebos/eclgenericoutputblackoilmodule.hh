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
#include <functional>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>

#include <opm/output/data/Wells.hpp>
#include <opm/output/eclipse/Inplace.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <ebos/eclinterregflows.hh>

namespace Opm {

namespace data { class Solution; }
class EclipseState;
class Schedule;
class SummaryConfig;
class SummaryConfigNode;
class SummaryState;

template<class FluidSystem, class Scalar>
class EclGenericOutputBlackoilModule {
public:
     Scalar* getPRESSURE_ptr(void) {
        return (this->fluidPressure_.data()) ;
    };
    
    int  getPRESSURE_size( void ) {
        return (this->fluidPressure_.size()) ;
    };
    
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
                         const Parallel::Communication& comm);

    // write Reservoir Volumes to output log
    Inplace outputFipresvLog(std::map<std::string, double>& miscSummaryData,
                         std::map<std::string, std::vector<double>>& regionData,
                         const bool substep,
                         const Parallel::Communication& comm);



    void outputErrorLog(const Parallel::Communication& comm) const;

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

    Scalar getSaltSaturation(unsigned elemIdx) const
    {
        if (pSalt_.size() > elemIdx)
            return pSalt_[elemIdx];

        return 0;
    }

    Scalar getPermFactor(unsigned elemIdx) const
    {
        if (permFact_.size() > elemIdx)
            return permFact_[elemIdx];

        return 0;
    }

    Scalar getMicrobialConcentration(unsigned elemIdx) const
    {
        if (cMicrobes_.size() > elemIdx)
            return cMicrobes_[elemIdx];

        return 0;
    }

    Scalar getOxygenConcentration(unsigned elemIdx) const
    {
        if (cOxygen_.size() > elemIdx)
            return cOxygen_[elemIdx];

        return 0;
    }

    Scalar getUreaConcentration(unsigned elemIdx) const
    {
        if (cUrea_.size() > elemIdx)
            return cUrea_[elemIdx];

        return 0;
    }

    Scalar getBiofilmConcentration(unsigned elemIdx) const
    {
        if (cBiofilm_.size() > elemIdx)
            return cBiofilm_[elemIdx];

        return 0;
    }

    Scalar getCalciteConcentration(unsigned elemIdx) const
    {
        if (cCalcite_.size() > elemIdx)
            return cCalcite_[elemIdx];

        return 0;
    }

    const std::array<std::pair<std::string, std::pair<std::vector<int>, std::vector<double>>>, 3>& getFlowsn() const
    {
        return this->flowsn_;
    }

    bool hasFlowsn() const
    {
        return enableFlowsn_;
    }

    bool hasFlows() const
    {
        return enableFlows_;
    }

    bool anyFlows() const
    {
        return anyFlows_;
    }

    const std::array<std::pair<std::string, std::pair<std::vector<int>, std::vector<double>>>, 3>& getFloresn() const
    {
        return this->floresn_;
    }

    bool hasFloresn() const
    {
        return enableFloresn_;
    }

    bool hasFlores() const
    {
        return enableFlores_;
    }

    bool anyFlores() const
    {
        return anyFlores_;
    }

    bool needInterfaceFluxes([[maybe_unused]] const bool isSubStep) const
    {
        return this->interRegionFlows_.wantInterRegflowSummary();
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

    bool localDataValid() const{
        return local_data_valid_;
    }

    void invalidateLocalData(){
        local_data_valid_ = false;
    }

    void validateLocalData(){
        local_data_valid_ = true;
    }

    // Virtual destructor for safer inheritance.
    virtual ~EclGenericOutputBlackoilModule();

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(initialInplace_);
    }

protected:
    using ScalarBuffer = std::vector<Scalar>;
    using StringBuffer = std::vector<std::string>;
    enum { numPhases = FluidSystem::numPhases };
    enum { numComponents = FluidSystem::numComponents };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };

    EclGenericOutputBlackoilModule(const EclipseState& eclState,
                                   const Schedule& schedule,
                                   const SummaryConfig& summaryConfig,
                                   const SummaryState& summaryState,
                                   bool enableEnergy,
                                   bool enableTemperature,
                                   bool enableSolvent,
                                   bool enablePolymer,
                                   bool enableFoam,
                                   bool enableBrine,
                                   bool enableSaltPrecipitation,
                                   bool enableExtbo,
                                   bool enableMICP);

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
                        unsigned numTracers,
                        unsigned numOutputNnc);

    void fipUnitConvert_(std::unordered_map<Inplace::Phase, Scalar>& fip) const;

    void pressureUnitConvert_(Scalar& pav) const;

    void outputRegionFluidInPlace_(std::unordered_map<Inplace::Phase, Scalar> oip,
                                   std::unordered_map<Inplace::Phase, Scalar> cip,
                                   const Scalar& pav, const int reg = 0) const;
    void outputResvFluidInPlace_(std::unordered_map<Inplace::Phase, Scalar> cipr,
                                 const int reg = 0) const;
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

    void outputFipresvLogImpl(const Inplace& inplace) const;

    void makeRegionSum(Inplace& inplace,
                       const std::string& region_name,
                       const Parallel::Communication& comm) const;

    Inplace accumulateRegionSums(const Parallel::Communication& comm);

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
                                  const std::size_t maxNumberOfRegions,
                                  const Parallel::Communication& comm);

    static int regionMax(const std::vector<int>& region,
                         const Parallel::Communication& comm);

    static void update(Inplace& inplace,
                       const std::string& region_name,
                       const Inplace::Phase phase,
                       const std::size_t ntFip,
                       const ScalarBuffer& values);

    static Scalar sum(const ScalarBuffer& v);

    void setupBlockData(std::function<bool(int)> isCartIdxOnThisRank);

    virtual bool isDefunctParallelWell(std::string wname) const = 0;

    const EclipseState& eclState_;
    const Schedule& schedule_;
    const SummaryConfig& summaryConfig_;
    const SummaryState& summaryState_;

    EclInterRegFlowMap interRegionFlows_;

    bool enableEnergy_;
    bool enableTemperature_;

    bool enableSolvent_;
    bool enablePolymer_;
    bool enableFoam_;
    bool enableBrine_;
    bool enableSaltPrecipitation_;
    bool enableExtbo_;
    bool enableMICP_;

    bool forceDisableFipOutput_;
    bool forceDisableFipresvOutput_;
    bool outputFipRestart_;
    bool computeFip_;

    bool anyFlows_;
    bool anyFlores_;
    bool enableFlows_;
    bool enableFlores_;
    bool enableFlowsn_;
    bool enableFloresn_;

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
    ScalarBuffer dynamicPoreVolume_;
    ScalarBuffer fluidPressure_;
    ScalarBuffer temperature_;
    ScalarBuffer rs_;
    ScalarBuffer rsw_;
    ScalarBuffer rv_;
    ScalarBuffer rvw_;
    ScalarBuffer overburdenPressure_;
    ScalarBuffer oilSaturationPressure_;
    ScalarBuffer drsdtcon_;
    ScalarBuffer sSol_;
    ScalarBuffer cPolymer_;
    ScalarBuffer cFoam_;
    ScalarBuffer cSalt_;
    ScalarBuffer pSalt_;
    ScalarBuffer permFact_;
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
    ScalarBuffer cMicrobes_;
    ScalarBuffer cOxygen_;
    ScalarBuffer cUrea_;
    ScalarBuffer cBiofilm_;
    ScalarBuffer cCalcite_;
    ScalarBuffer pcow_;
    ScalarBuffer pcog_;

    std::array<ScalarBuffer, numPhases> saturation_;
    std::array<ScalarBuffer, numPhases> invB_;
    std::array<ScalarBuffer, numPhases> density_;
    std::array<ScalarBuffer, numPhases> viscosity_;
    std::array<ScalarBuffer, numPhases> relativePermeability_;

    std::vector<ScalarBuffer> tracerConcentrations_;

    std::array<ScalarBuffer, numPhases> flowsi_;
    std::array<ScalarBuffer, numPhases> flowsj_;
    std::array<ScalarBuffer, numPhases> flowsk_;
    std::array<ScalarBuffer, numPhases> floresi_;
    std::array<ScalarBuffer, numPhases> floresj_;
    std::array<ScalarBuffer, numPhases> floresk_;

    std::array<std::pair<std::string, std::pair<std::vector<int>, ScalarBuffer>>, 3> floresn_;
    std::array<std::pair<std::string, std::pair<std::vector<int>, ScalarBuffer>>, 3> flowsn_;

    std::map<size_t, Scalar> oilConnectionPressures_;
    std::map<size_t, Scalar> waterConnectionSaturations_;
    std::map<size_t, Scalar> gasConnectionSaturations_;
    std::map<std::pair<std::string, int>, double> blockData_;
    std::map<std::size_t , double> wbpData_;

    std::optional<Inplace> initialInplace_;
    bool local_data_valid_;
};

} // namespace Opm

#endif
