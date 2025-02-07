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
 * \copydoc Opm::OutputBlackOilModule
 */
#ifndef OPM_GENERIC_OUTPUT_BLACK_OIL_MODULE_HPP
#define OPM_GENERIC_OUTPUT_BLACK_OIL_MODULE_HPP

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/output/data/Wells.hpp>
#include <opm/output/eclipse/Inplace.hpp>

#include <opm/simulators/flow/ExtboContainer.hpp>
#include <opm/simulators/flow/FIPContainer.hpp>
#include <opm/simulators/flow/FlowsData.hpp>
#include <opm/simulators/flow/InterRegFlows.hpp>
#include <opm/simulators/flow/LogOutputHelper.hpp>
#include <opm/simulators/flow/MechContainer.hpp>
#include <opm/simulators/flow/MICPContainer.hpp>
#include <opm/simulators/flow/RegionPhasePVAverage.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm::Parameters {

struct ForceDisableFluidInPlaceOutput { static constexpr bool value = false; };
struct ForceDisableResvFluidInPlaceOutput { static constexpr bool value = false; };

} // namespace Opm::Parameters

namespace Opm {

namespace data { class Solution; }
class EclipseState;
class Schedule;
class SummaryConfig;
class SummaryConfigNode;
class SummaryState;

template<class FluidSystem>
class GenericOutputBlackoilModule {
public:
    using Scalar = typename FluidSystem::Scalar;

    // Virtual destructor for safer inheritance.
    virtual ~GenericOutputBlackoilModule();

    Scalar* getPRESSURE_ptr()
    {
        return this->fluidPressure_.data();
    }

    int getPRESSURE_size()
    {
        return this->fluidPressure_.size();
    }

    /*!
     * \brief Register all run-time parameters for the Vtk output module.
     */
    static void registerParameters();

    void outputTimeStamp(const std::string& lbl,
                         double elapsed,
                         int rstep,
                         boost::posix_time::ptime currentDate);

    /// Clear internal arrays for parallel accumulation of per-region phase
    /// density averages.
    void prepareDensityAccumulation();

    /// Run cross-rank parallel accumulation of per-region phase density
    /// running sums (average values).
    void accumulateDensityParallel();

    // write cumulative production and injection reports to output
    void outputCumLog(std::size_t reportStepNum);

    // write production report to output
    void outputProdLog(std::size_t reportStepNum);

    // write injection report to output
    void outputInjLog(std::size_t reportStepNum);

    // calculate Initial Fluid In Place
    Inplace calc_initial_inplace(const Parallel::Communication& comm);

    // calculate Fluid In Place
    Inplace calc_inplace(std::map<std::string, double>& miscSummaryData,
                         std::map<std::string, std::vector<double>>& regionData,
                         const Parallel::Communication& comm);

    void outputFipAndResvLog(const Inplace& inplace,
                         const std::size_t reportStepNum,
                         double elapsed,
                         boost::posix_time::ptime currentDate,
                         const bool substep,
                         const Parallel::Communication& comm);

    void outputErrorLog(const Parallel::Communication& comm) const;

    void accumulateRftDataParallel(const Parallel::Communication& comm);

    void addRftDataToWells(data::Wells& wellDatas,
                           std::size_t reportStepNum);

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

    Scalar getSolventRsw(unsigned elemIdx) const
    {
        if (rswSol_.size() > elemIdx)
            return rswSol_[elemIdx];

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

    const MICPContainer<Scalar>& getMICP() const
    { return this->micpC_; }

    const std::array<FlowsData<double>, 3>& getFlowsn() const
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

    bool hasBlockFlows() const
    {
        return blockFlows_;
    }

    bool anyFlows() const
    {
        return anyFlows_;
    }

    const std::array<FlowsData<double>, 3>& getFloresn() const
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

    void setCnvData(const std::vector<std::vector<int>>& data)
    {
        cnvData_ = data;
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(initialInplace_);
    }

    //! \brief Assign fields that are in global numbering to the solution.
    //! \detail This is used to add fields that for some reason cannot be collected
    //!         using the regular collect mechanism to the solution. In particular this
    //!         is used with RPTRST CONV output.
    void assignGlobalFieldsToSolution(data::Solution& sol);

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
    using Dir = FaceDir::DirEnum;

    GenericOutputBlackoilModule(const EclipseState& eclState,
                                const Schedule& schedule,
                                const SummaryConfig& summaryConfig,
                                const SummaryState& summaryState,
                                const std::string& moduleVersionName,
                                bool enableEnergy,
                                bool enableTemperature,
                                bool enableMech,
                                bool enableSolvent,
                                bool enablePolymer,
                                bool enableFoam,
                                bool enableBrine,
                                bool enableSaltPrecipitation,
                                bool enableExtbo,
                                bool enableMICP,
                                bool isCompositional = false);

    void doAllocBuffers(unsigned bufferSize,
                        unsigned reportStepNum,
                        const bool substep,
                        const bool log,
                        const bool isRestart,
                        const bool vapparsActive = false,
                        const bool enablePCHysteresis = false,
                        const bool enableNonWettingHysteresis =false,
                        const bool enableWettingHysteresis = false,
                        unsigned numTracers = 0,
                        const std::vector<bool>& enableSolTracers = {},
                        unsigned numOutputNnc = 0);

    void makeRegionSum(Inplace& inplace,
                       const std::string& region_name,
                       const Parallel::Communication& comm) const;

    Inplace accumulateRegionSums(const Parallel::Communication& comm);

    void updateSummaryRegionValues(const Inplace& inplace,
                                   std::map<std::string, double>& miscSummaryData,
                                   std::map<std::string, std::vector<double>>& regionData) const;

    static bool isOutputCreationDirective_(const std::string& keyword);

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

    void gatherAndUpdateRftMap(std::map<std::size_t, Scalar>& local_map, const Parallel::Communication& comm);

    const EclipseState& eclState_;
    const Schedule& schedule_;
    const SummaryState& summaryState_;

    SummaryConfig summaryConfig_;

    InterRegFlowMap interRegionFlows_;
    LogOutputHelper<Scalar> logOutput_;

    bool enableEnergy_{false};
    bool enableTemperature_{false};
    bool enableMech_{false};

    bool enableSolvent_{false};
    bool enablePolymer_{false};
    bool enableFoam_{false};
    bool enableBrine_{false};
    bool enableSaltPrecipitation_{false};
    bool enableExtbo_{false};
    bool enableMICP_{false};
    bool isCompositional_{false};

    bool forceDisableFipOutput_{false};
    bool forceDisableFipresvOutput_{false};
    bool computeFip_{false};

    bool anyFlows_{false};
    bool anyFlores_{false};
    bool blockFlows_{false};
    bool enableFlows_{false};
    bool enableFlores_{false};
    bool enableFlowsn_{false};
    bool enableFloresn_{false};

    FIPContainer<FluidSystem> fipC_;
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
    ScalarBuffer rPorV_;
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
    ScalarBuffer rswSol_;
    ScalarBuffer cPolymer_;
    ScalarBuffer cFoam_;
    ScalarBuffer cSalt_;
    ScalarBuffer pSalt_;
    ScalarBuffer permFact_;
    ExtboContainer<Scalar> extboC_;
    ScalarBuffer soMax_;
    ScalarBuffer swMax_;
    ScalarBuffer sgmax_;
    ScalarBuffer shmax_;
    ScalarBuffer somin_;
    ScalarBuffer swmin_;
    ScalarBuffer ppcw_;
    ScalarBuffer gasDissolutionFactor_;
    ScalarBuffer oilVaporizationFactor_;
    ScalarBuffer gasDissolutionFactorInWater_;
    ScalarBuffer waterVaporizationFactor_;
    ScalarBuffer bubblePointPressure_;
    ScalarBuffer dewPointPressure_;
    ScalarBuffer rockCompPorvMultiplier_;
    ScalarBuffer minimumOilPressure_;
    ScalarBuffer saturatedOilFormationVolumeFactor_;
    ScalarBuffer rockCompTransMultiplier_;
    MICPContainer<Scalar> micpC_;
    ScalarBuffer pcgw_;
    ScalarBuffer pcow_;
    ScalarBuffer pcog_;

    // buffers for mechanical output
    MechContainer<Scalar> mech_;

    std::array<ScalarBuffer, numPhases> saturation_;
    std::array<ScalarBuffer, numPhases> invB_;
    std::array<ScalarBuffer, numPhases> density_;
    std::array<ScalarBuffer, numPhases> viscosity_;
    std::array<ScalarBuffer, numPhases> relativePermeability_;

    // total mole fractions for each component
    std::array<ScalarBuffer, numComponents> moleFractions_;
    // mole fractions for each component in each phase
    std::array<std::array<ScalarBuffer, numComponents>, numPhases> phaseMoleFractions_;
    std::vector<ScalarBuffer> freeTracerConcentrations_;
    std::vector<ScalarBuffer> solTracerConcentrations_;

    std::array<ScalarBuffer, numPhases> residual_;

    std::array<std::array<ScalarBuffer, numPhases>, 6> flows_;
    std::array<std::array<ScalarBuffer, numPhases>, 6> flores_;

    std::array<FlowsData<double>, 3> floresn_;
    std::array<FlowsData<double>, 3> flowsn_;

    std::map<std::size_t, Scalar> oilConnectionPressures_;
    std::map<std::size_t, Scalar> waterConnectionSaturations_;
    std::map<std::size_t, Scalar> gasConnectionSaturations_;
    std::map<std::pair<std::string, int>, double> blockData_;

    std::vector<std::vector<int>> cnvData_; //!< Data for CNV_xxx arrays

    std::optional<Inplace> initialInplace_;
    bool local_data_valid_{false};

    std::optional<RegionPhasePoreVolAverage> regionAvgDensity_;
};

} // namespace Opm

#endif // OPM_GENERIC_OUTPUT_BLACK_OIL_MODULE_HPP
