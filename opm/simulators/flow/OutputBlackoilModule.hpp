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
#ifndef OPM_OUTPUT_BLACK_OIL_MODULE_HPP
#define OPM_OUTPUT_BLACK_OIL_MODULE_HPP

#include <dune/common/fvector.hh>

#include <opm/simulators/utils/moduleVersion.hpp>

#include <opm/common/Exceptions.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/Visitor.hpp>

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/utils/parametersystem.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/output/data/Cells.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/Inplace.hpp>

#include <opm/simulators/flow/CollectDataOnIORank.hpp>
#include <opm/simulators/flow/FlowBaseVanguard.hpp>
#include <opm/simulators/flow/GenericOutputBlackoilModule.hpp>
#include <opm/simulators/flow/OutputExtractor.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm {

// forward declaration
template <class TypeTag>
class EcfvDiscretization;

/*!
 * \ingroup BlackOilSimulator
 *
 * \brief Output module for the results black oil model writing in
 *        ECL binary format.
 */
template <class TypeTag>
class OutputBlackOilModule : public GenericOutputBlackoilModule<GetPropType<TypeTag, Properties::FluidSystem>>
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Discretization = GetPropType<TypeTag, Properties::Discretization>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = typename IntensiveQuantities::FluidState;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    using BaseType = GenericOutputBlackoilModule<FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Dir = FaceDir::DirEnum;
    using BlockExtractor = detail::BlockExtractor<TypeTag>;
    using Extractor = detail::Extractor<TypeTag>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using CollectDataOnIORankType = CollectDataOnIORank<Grid,EquilGrid,GridView>;

    static constexpr int conti0EqIdx = Indices::conti0EqIdx;
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static constexpr int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;
    static constexpr int gasCompIdx = FluidSystem::gasCompIdx;
    static constexpr int oilCompIdx = FluidSystem::oilCompIdx;
    static constexpr int waterCompIdx = FluidSystem::waterCompIdx;
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableMICP = getPropValue<TypeTag, Properties::EnableMICP>() };
    enum { enableVapwat = getPropValue<TypeTag, Properties::EnableVapwat>() };
    enum { enableDisgasInWater = getPropValue<TypeTag, Properties::EnableDisgasInWater>() };
    enum { enableDissolvedGas = Indices::compositionSwitchIdx >= 0 };

    template<int idx, class VectorType>
    static Scalar value_or_zero(const VectorType& v)
    {
        if constexpr (idx == -1) {
            return 0.0;
        } else {
            return v.empty() ? 0.0 : v[idx];
        }
    }

public:
    OutputBlackOilModule(const Simulator& simulator,
                         const SummaryConfig& smryCfg,
                         const CollectDataOnIORankType& collectOnIORank)
        : BaseType(simulator.vanguard().eclState(),
                   simulator.vanguard().schedule(),
                   smryCfg,
                   simulator.vanguard().summaryState(),
                   moduleVersionName(),
                   [this](const int idx)
                   { return simulator_.problem().eclWriter().collectOnIORank().localIdxToGlobalIdx(idx); },
                   simulator.vanguard().grid().comm(),
                   getPropValue<TypeTag, Properties::EnableEnergy>(),
                   getPropValue<TypeTag, Properties::EnableTemperature>(),
                   getPropValue<TypeTag, Properties::EnableMech>(),
                   getPropValue<TypeTag, Properties::EnableSolvent>(),
                   getPropValue<TypeTag, Properties::EnablePolymer>(),
                   getPropValue<TypeTag, Properties::EnableFoam>(),
                   getPropValue<TypeTag, Properties::EnableBrine>(),
                   getPropValue<TypeTag, Properties::EnableSaltPrecipitation>(),
                   getPropValue<TypeTag, Properties::EnableExtbo>(),
                   getPropValue<TypeTag, Properties::EnableMICP>())
        , simulator_(simulator)
        , collectOnIORank_(collectOnIORank)
    {
        for (auto& region_pair : this->regions_) {
            this->createLocalRegion_(region_pair.second);
        }

        auto isCartIdxOnThisRank = [&collectOnIORank](const int idx) {
            return collectOnIORank.isCartIdxOnThisRank(idx);
        };

        this->setupBlockData(isCartIdxOnThisRank);

        if (! Parameters::Get<Parameters::OwnerCellsFirst>()) {
            const std::string msg = "The output code does not support --owner-cells-first=false.";
            if (collectOnIORank.isIORank()) {
                OpmLog::error(msg);
            }
            OPM_THROW_NOLOG(std::runtime_error, msg);
        }

        if (smryCfg.match("[FB]PP[OGW]") || smryCfg.match("RPP[OGW]*")) {
            auto rset = this->eclState_.fieldProps().fip_regions();
            rset.push_back("PVTNUM");

            // Note: We explicitly use decltype(auto) here because the
            // default scheme (-> auto) will deduce an undesirable type.  We
            // need the "reference to vector" semantics in this instance.
            this->regionAvgDensity_
                .emplace(this->simulator_.gridView().comm(),
                         FluidSystem::numPhases, rset,
                         [fp = std::cref(this->eclState_.fieldProps())]
                         (const std::string& rsetName) -> decltype(auto)
                         { return fp.get().get_int(rsetName); });
        }
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to ECL output files
     */
    void
    allocBuffers(const unsigned bufferSize,
                 const unsigned reportStepNum,
                 const bool     substep,
                 const bool     log,
                 const bool     isRestart)
    {
        if (! std::is_same<Discretization, EcfvDiscretization<TypeTag>>::value) {
            return;
        }

        const auto& problem = this->simulator_.problem();

        this->doAllocBuffers(bufferSize,
                             reportStepNum,
                             substep,
                             log,
                             isRestart,
                             &problem.materialLawManager()->hysteresisConfig(),
                             problem.eclWriter().getOutputNnc().size());
    }

    //! \brief Setup list of active element-level data extractors
    void setupExtractors(const bool isSubStep,
                         const int  reportStepNum)
    {
        this->setupElementExtractors_();
        this->setupBlockExtractors_(isSubStep, reportStepNum);
    }

    //! \brief Clear list of active element-level data extractors
    void clearExtractors()
    {
        this->extractors_.clear();
        this->blockExtractors_.clear();
        this->extraBlockExtractors_.clear();
    }

    /*!
     * \brief Modify the internal buffers according to the intensive
     *        quanties relevant for an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        OPM_TIMEBLOCK_LOCAL(processElement);
        if (!std::is_same<Discretization, EcfvDiscretization<TypeTag>>::value) {
            return;
        }

        if (this->extractors_.empty()) {
            assert(0);
        }

        const auto& matLawManager = simulator_.problem().materialLawManager();

        typename Extractor::HysteresisParams hysterParams;
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            const typename Extractor::Context ectx{
                elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0),
                elemCtx.primaryVars(dofIdx, /*timeIdx=*/0).pvtRegionIndex(),
                elemCtx.simulator().episodeIndex(),
                fs,
                intQuants,
                hysterParams
            };

            if (matLawManager->enableHysteresis()) {
                if (FluidSystem::phaseIsActive(oilPhaseIdx) && FluidSystem::phaseIsActive(waterPhaseIdx)) {
                    matLawManager->oilWaterHysteresisParams(hysterParams.somax,
                                                            hysterParams.swmax,
                                                            hysterParams.swmin,
                                                            ectx.globalDofIdx);
                }
                if (FluidSystem::phaseIsActive(oilPhaseIdx) && FluidSystem::phaseIsActive(gasPhaseIdx)) {
                    matLawManager->gasOilHysteresisParams(hysterParams.sgmax,
                                                          hysterParams.shmax,
                                                          hysterParams.somin,
                                                          ectx.globalDofIdx);
                }
            }

            Extractor::process(ectx, extractors_);
        }
    }

    void processElementBlockData(const ElementContext& elemCtx)
    {
        OPM_TIMEBLOCK_LOCAL(processElementBlockData);
        if (!std::is_same<Discretization, EcfvDiscretization<TypeTag>>::value) {
            return;
        }

        if (this->blockExtractors_.empty() && this->extraBlockExtractors_.empty()) {
            return;
        }

        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            // Adding block data
            const auto globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
            const auto cartesianIdx = elemCtx.simulator().vanguard().cartesianIndex(globalDofIdx);

            const auto be_it = this->blockExtractors_.find(cartesianIdx);
            const auto bee_it = this->extraBlockExtractors_.find(cartesianIdx);
            if (be_it == this->blockExtractors_.end() &&
                bee_it == this->extraBlockExtractors_.end())
            {
                continue;
            }

            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            const typename BlockExtractor::Context ectx{
                globalDofIdx,
                dofIdx,
                fs,
                intQuants,
                elemCtx,
            };

            if (be_it != this->blockExtractors_.end()) {
                BlockExtractor::process(be_it->second, ectx);
            }
            if (bee_it != this->extraBlockExtractors_.end()) {
                BlockExtractor::process(bee_it->second, ectx);
            }
        }
    }

    void outputFipAndResvLog(const Inplace& inplace,
                             const std::size_t reportStepNum,
                             double elapsed,
                             boost::posix_time::ptime currentDate,
                             const bool substep,
                             const Parallel::Communication& comm)
    {

        if (comm.rank() != 0) {
            return;
        }

        // For report step 0 we use the RPTSOL config, else derive from RPTSCHED
        std::unique_ptr<FIPConfig> fipSched;
        if (reportStepNum > 0) {
            const auto& rpt = this->schedule_[reportStepNum-1].rpt_config.get();
            fipSched = std::make_unique<FIPConfig>(rpt);
        }
        const FIPConfig& fipc = reportStepNum == 0 ? this->eclState_.getEclipseConfig().fip()
                                                   : *fipSched;

        if (!substep && !this->forceDisableFipOutput_ && fipc.output(FIPConfig::OutputField::FIELD)) {

            this->logOutput_.timeStamp("BALANCE", elapsed, reportStepNum, currentDate);

            const auto& initial_inplace = this->initialInplace().value();
            this->logOutput_.fip(inplace, initial_inplace, "");

            if (fipc.output(FIPConfig::OutputField::FIPNUM)) {
                this->logOutput_.fip(inplace, initial_inplace, "FIPNUM");

                if (fipc.output(FIPConfig::OutputField::RESV))
                    this->logOutput_.fipResv(inplace, "FIPNUM");
            }

            if (fipc.output(FIPConfig::OutputField::FIP)) {
                for (const auto& reg : this->regions_) {
                    if (reg.first != "FIPNUM") {
                        std::ostringstream ss;
                        ss << "BAL" << reg.first.substr(3);
                        this->logOutput_.timeStamp(ss.str(), elapsed, reportStepNum, currentDate);
                        this->logOutput_.fip(inplace, initial_inplace, reg.first);

                        if (fipc.output(FIPConfig::OutputField::RESV))
                            this->logOutput_.fipResv(inplace, reg.first);
                    }
                }
            }
        }
    }

    void outputFipAndResvLogToCSV(const std::size_t reportStepNum,
                                  const bool substep,
                                  const Parallel::Communication& comm)
    {
        if (comm.rank() != 0) {
            return;
        }

        if ((reportStepNum == 0) && (!substep) &&
            (this->schedule_.initialReportConfiguration().has_value()) &&
            (this->schedule_.initialReportConfiguration()->contains("CSVFIP"))) {

            std::ostringstream csv_stream;

            this->logOutput_.csv_header(csv_stream);

            const auto& initial_inplace = this->initialInplace().value();

            this->logOutput_.fip_csv(csv_stream, initial_inplace, "FIPNUM");

            for (const auto& reg : this->regions_) {
                if (reg.first != "FIPNUM") {
                    this->logOutput_.fip_csv(csv_stream, initial_inplace, reg.first);
                }
            }

            const IOConfig& io =  this->eclState_.getIOConfig();
            auto csv_fname = io.getOutputDir() + "/" + io.getBaseName() + ".CSV";

            std::ofstream outputFile(csv_fname);

            outputFile <<  csv_stream.str();

            outputFile.close();
        }
    }

    /*!
     * \brief Capture connection fluxes, particularly to account for inter-region flows.
     *
     * \tparam ActiveIndex Callable type, typically a lambda, that enables
     *    retrieving the active index, on the local MPI rank, of a
     *    particular cell/element.  Must support a function call operator of
     *    the form
     \code
        int operator()(const Element& elem) const
     \endcode
     *
     * \tparam CartesianIndex Callable type, typically a lambda, that
     *    enables retrieving the globally unique Cartesian index of a
     *    particular cell/element given its active index on the local MPI
     *    rank.  Must support a function call operator of the form
     \code
        int operator()(const int activeIndex) const
     \endcode
     *
     * \param[in] elemCtx Primary lookup structure for per-cell/element
     *    dynamic information.
     *
     * \param[in] activeIndex Mapping from cell/elements to linear indices
     *    on local MPI rank.
     *
     * \param[in] cartesianIndex Mapping from active index on local MPI rank
     *    to globally unique Cartesian cell/element index.
     */
    template <class ActiveIndex, class CartesianIndex>
    void processFluxes(const ElementContext& elemCtx,
                       ActiveIndex&&         activeIndex,
                       CartesianIndex&&      cartesianIndex)
    {
        OPM_TIMEBLOCK_LOCAL(processFluxes);
        const auto identifyCell = [&activeIndex, &cartesianIndex](const Element& elem)
            -> InterRegFlowMap::Cell
        {
            const auto cellIndex = activeIndex(elem);

            return {
                static_cast<int>(cellIndex),
                cartesianIndex(cellIndex),
                elem.partitionType() == Dune::InteriorEntity
            };
        };

        const auto timeIdx = 0u;
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto numInteriorFaces = elemCtx.numInteriorFaces(timeIdx);

        for (auto scvfIdx = 0 * numInteriorFaces; scvfIdx < numInteriorFaces; ++scvfIdx) {
            const auto& face = stencil.interiorFace(scvfIdx);
            const auto left  = identifyCell(stencil.element(face.interiorIndex()));
            const auto right = identifyCell(stencil.element(face.exteriorIndex()));

            const auto rates = this->
                getComponentSurfaceRates(elemCtx, face.area(), scvfIdx, timeIdx);

            this->interRegionFlows_.addConnection(left, right, rates);
        }
    }

    /*!
     * \brief Prepare for capturing connection fluxes, particularly to
     *    account for inter-region flows.
     */
    void initializeFluxData()
    {
        // Inter-region flow rates.  Note: ".clear()" prepares to accumulate
        // contributions per bulk connection between FIP regions.
        this->interRegionFlows_.clear();
    }

    /*!
     * \brief Finalize capturing connection fluxes.
     */
    void finalizeFluxData()
    {
        this->interRegionFlows_.compress();
    }

    /*!
     * \brief Get read-only access to collection of inter-region flows.
     */
    const InterRegFlowMap& getInterRegFlows() const
    {
        return this->interRegionFlows_;
    }

    template <class FluidState>
    void assignToFluidState(FluidState& fs, unsigned elemIdx) const
    {
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (this->saturation_[phaseIdx].empty())
                continue;

            fs.setSaturation(phaseIdx, this->saturation_[phaseIdx][elemIdx]);
        }

        if (!this->fluidPressure_.empty()) {
            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            std::array<Scalar, numPhases> pc = {0};
            const MaterialLawParams& matParams = simulator_.problem().materialLawParams(elemIdx);
            MaterialLaw::capillaryPressures(pc, matParams, fs);
            Valgrind::CheckDefined(this->fluidPressure_[elemIdx]);
            Valgrind::CheckDefined(pc);
            const auto& pressure = this->fluidPressure_[elemIdx];
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                if (Indices::oilEnabled)
                    fs.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[oilPhaseIdx]));
                else if (Indices::gasEnabled)
                    fs.setPressure(phaseIdx, pressure + (pc[phaseIdx] - pc[gasPhaseIdx]));
                else if (Indices::waterEnabled)
                    //single (water) phase
                    fs.setPressure(phaseIdx, pressure);
            }
        }

        if (!this->temperature_.empty())
            fs.setTemperature(this->temperature_[elemIdx]);
        if constexpr (enableDissolvedGas) {
            if (!this->rs_.empty())
                fs.setRs(this->rs_[elemIdx]);
            if (!this->rv_.empty())
                fs.setRv(this->rv_[elemIdx]);
        }
        if constexpr (enableDisgasInWater) {
            if (!this->rsw_.empty())
                fs.setRsw(this->rsw_[elemIdx]);
        }
        if constexpr (enableVapwat) {
            if (!this->rvw_.empty())
                fs.setRvw(this->rvw_[elemIdx]);
        }
    }

    void initHysteresisParams(Simulator& simulator, unsigned elemIdx) const
    {
        if (!this->soMax_.empty())
            simulator.problem().setMaxOilSaturation(elemIdx, this->soMax_[elemIdx]);

        if (simulator.problem().materialLawManager()->enableHysteresis()) {
            auto matLawManager = simulator.problem().materialLawManager();

            if (FluidSystem::phaseIsActive(oilPhaseIdx)
                && FluidSystem::phaseIsActive(waterPhaseIdx)) {
                    Scalar somax = 2.0;
                    Scalar swmax = -2.0;
                    Scalar swmin = 2.0;

                if (matLawManager->enableNonWettingHysteresis()) {
                    if (!this->soMax_.empty()) {
                        somax = this->soMax_[elemIdx];
                    }
                }
                if (matLawManager->enableWettingHysteresis()) {
                    if (!this->swMax_.empty()) {
                        swmax = this->swMax_[elemIdx];
                    }
                }
                if (matLawManager->enablePCHysteresis()) {
                    if (!this->swmin_.empty()) {
                        swmin = this->swmin_[elemIdx];
                    }
                }
                matLawManager->setOilWaterHysteresisParams(
                        somax, swmax, swmin, elemIdx);
            }
            if (FluidSystem::phaseIsActive(oilPhaseIdx)
                && FluidSystem::phaseIsActive(gasPhaseIdx)) {
                    Scalar sgmax = 2.0;
                    Scalar shmax = -2.0;
                    Scalar somin = 2.0;

                if (matLawManager->enableNonWettingHysteresis()) {
                    if (!this->sgmax_.empty()) {
                        sgmax = this->sgmax_[elemIdx];
                    }
                }
                if (matLawManager->enableWettingHysteresis()) {
                    if (!this->shmax_.empty()) {
                        shmax = this->shmax_[elemIdx];
                    }
                }
                if (matLawManager->enablePCHysteresis()) {
                    if (!this->somin_.empty()) {
                        somin = this->somin_[elemIdx];
                    }
                }
                matLawManager->setGasOilHysteresisParams(
                        sgmax, shmax, somin, elemIdx);
            }

        }

        if (simulator_.vanguard().eclState().fieldProps().has_double("SWATINIT")) {
            simulator.problem().materialLawManager()
                ->applyRestartSwatInit(elemIdx, this->ppcw_[elemIdx]);
        }
    }

    void updateFluidInPlace(const ElementContext& elemCtx)
    {
        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            updateFluidInPlace_(elemCtx, dofIdx);
        }
    }

    void updateFluidInPlace(const unsigned             globalDofIdx,
                            const IntensiveQuantities& intQuants,
                            const double               totVolume)
    {
        this->updateFluidInPlace_(globalDofIdx, intQuants, totVolume);
    }

private:
    template <typename T>
    using RemoveCVR = std::remove_cv_t<std::remove_reference_t<T>>;

    template <typename, class = void>
    struct HasGeoMech : public std::false_type {};

    template <typename Problem>
    struct HasGeoMech<
        Problem, std::void_t<decltype(std::declval<Problem>().geoMechModel())>
    > : public std::true_type {};

    bool isDefunctParallelWell(const std::string& wname) const override
    {
        if (simulator_.gridView().comm().size() == 1)
            return false;
        const auto& parallelWells = simulator_.vanguard().parallelWells();
        std::pair<std::string, bool> value {wname, true};
        auto candidate = std::lower_bound(parallelWells.begin(), parallelWells.end(), value);
        return candidate == parallelWells.end() || *candidate != value;
    }

    bool isOwnedByCurrentRank(const std::string& wname) const override
    {
        return this->simulator_.problem().wellModel().isOwner(wname);
    }

    void updateFluidInPlace_(const ElementContext& elemCtx, const unsigned dofIdx)
    {
        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
        const unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        const auto totVolume = elemCtx.simulator().model().dofTotalVolume(globalDofIdx);

        this->updateFluidInPlace_(globalDofIdx, intQuants, totVolume);
    }

    void updateFluidInPlace_(const unsigned             globalDofIdx,
                             const IntensiveQuantities& intQuants,
                             const double               totVolume)
    {
        OPM_TIMEBLOCK_LOCAL(updateFluidInPlace);

        this->updateTotalVolumesAndPressures_(globalDofIdx, intQuants, totVolume);

        if (this->computeFip_) {
            this->updatePhaseInplaceVolumes_(globalDofIdx, intQuants, totVolume);
        }
    }

    void createLocalRegion_(std::vector<int>& region)
    {
        // For CpGrid with LGRs, where level zero grid has been distributed,
        // resize region is needed, since in this case the total amount of
        // element - per process - in level zero grid and leaf grid do not
        // coincide, in general.
        region.resize(simulator_.gridView().size(0));
        std::size_t elemIdx = 0;
        for (const auto& elem : elements(simulator_.gridView())) {
            if (elem.partitionType() != Dune::InteriorEntity) {
                region[elemIdx] = 0;
            }

            ++elemIdx;
        }
    }

    template <typename FluidState>
    void aggregateAverageDensityContributions_(const FluidState&  fs,
                                               const unsigned int globalDofIdx,
                                               const double       porv)
    {
        auto pvCellValue = RegionPhasePoreVolAverage::CellValue{};
        pvCellValue.porv = porv;

        for (auto phaseIdx = 0*FluidSystem::numPhases;
             phaseIdx < FluidSystem::numPhases; ++phaseIdx)
        {
            if (! FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            pvCellValue.value = getValue(fs.density(phaseIdx));
            pvCellValue.sat   = getValue(fs.saturation(phaseIdx));

            this->regionAvgDensity_
                ->addCell(globalDofIdx,
                          RegionPhasePoreVolAverage::Phase { phaseIdx },
                          pvCellValue);
        }
    }

    /*!
     * \brief Compute surface level component flow rates across a single
     *   intersection.
     *
     * \param[in] elemCtx Primary lookup structure for per-cell/element
     *    dynamic information.
     *
     * \param[in] scvfIdx Linear index of current interior bulk connection.
     *
     * \param[in] timeIdx Historical time-point at which to evaluate dynamic
     *    quantities (e.g., reciprocal FVF or dissolved gas concentration).
     *    Zero for the current time.
     *
     * \return Surface level component flow rates.
     */
    data::InterRegFlowMap::FlowRates
    getComponentSurfaceRates(const ElementContext& elemCtx,
                             const Scalar          faceArea,
                             const std::size_t     scvfIdx,
                             const std::size_t     timeIdx) const
    {
        using Component = data::InterRegFlowMap::Component;

        auto rates = data::InterRegFlowMap::FlowRates {};

        const auto& extQuant = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        const auto alpha = getValue(extQuant.extrusionFactor()) * faceArea;

        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            const auto& up = elemCtx
                .intensiveQuantities(extQuant.upstreamIndex(oilPhaseIdx), timeIdx);

            const auto pvtReg = up.pvtRegionIndex();

            const auto bO = getValue(getInvB_<FluidSystem, FluidState, Scalar>
                                     (up.fluidState(), oilPhaseIdx, pvtReg));

            const auto qO = alpha * bO * getValue(extQuant.volumeFlux(oilPhaseIdx));

            rates[Component::Oil] += qO;

            if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                const auto Rs = getValue(
                    BlackOil::getRs_<FluidSystem, FluidState, Scalar>
                    (up.fluidState(), pvtReg));

                rates[Component::Gas]    += qO * Rs;
                rates[Component::Disgas] += qO * Rs;
            }
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            const auto& up = elemCtx
                .intensiveQuantities(extQuant.upstreamIndex(gasPhaseIdx), timeIdx);

            const auto pvtReg = up.pvtRegionIndex();

            const auto bG = getValue(getInvB_<FluidSystem, FluidState, Scalar>
                                     (up.fluidState(), gasPhaseIdx, pvtReg));

            const auto qG = alpha * bG * getValue(extQuant.volumeFlux(gasPhaseIdx));

            rates[Component::Gas] += qG;

            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                const auto Rv = getValue(
                    BlackOil::getRv_<FluidSystem, FluidState, Scalar>
                    (up.fluidState(), pvtReg));

                rates[Component::Oil]    += qG * Rv;
                rates[Component::Vapoil] += qG * Rv;
            }
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
            const auto& up = elemCtx
                .intensiveQuantities(extQuant.upstreamIndex(waterPhaseIdx), timeIdx);

            const auto pvtReg = up.pvtRegionIndex();

            const auto bW = getValue(getInvB_<FluidSystem, FluidState, Scalar>
                                     (up.fluidState(), waterPhaseIdx, pvtReg));

            rates[Component::Water] +=
                alpha * bW * getValue(extQuant.volumeFlux(waterPhaseIdx));
        }

        return rates;
    }

    template <typename FluidState>
    Scalar hydroCarbonFraction(const FluidState& fs) const
    {
        if (this->eclState_.runspec().co2Storage()) {
            // CO2 storage: Hydrocarbon volume is full pore-volume.
            return 1.0;
        }

        // Common case.  Hydrocarbon volume is fraction occupied by actual
        // hydrocarbons.
        auto hydrocarbon = Scalar {0};
        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            hydrocarbon += getValue(fs.saturation(oilPhaseIdx));
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            hydrocarbon += getValue(fs.saturation(gasPhaseIdx));
        }

        return hydrocarbon;
    }

    void updateTotalVolumesAndPressures_(const unsigned             globalDofIdx,
                                         const IntensiveQuantities& intQuants,
                                         const double               totVolume)
    {
        const auto& fs = intQuants.fluidState();

        const double pv = totVolume * intQuants.porosity().value();
        const auto hydrocarbon = this->hydroCarbonFraction(fs);

        if (! this->hydrocarbonPoreVolume_.empty()) {
            this->fipC_.assignPoreVolume(globalDofIdx,
                                         totVolume * intQuants.referencePorosity());

            this->dynamicPoreVolume_[globalDofIdx] = pv;
            this->hydrocarbonPoreVolume_[globalDofIdx] = pv * hydrocarbon;
        }

        if (!this->pressureTimesHydrocarbonVolume_.empty() &&
            !this->pressureTimesPoreVolume_.empty())
        {
            assert(this->hydrocarbonPoreVolume_.size() == this->pressureTimesHydrocarbonVolume_.size());
            assert(this->fipC_.get(Inplace::Phase::PoreVolume).size() == this->pressureTimesPoreVolume_.size());

            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                this->pressureTimesPoreVolume_[globalDofIdx] =
                    getValue(fs.pressure(oilPhaseIdx)) * pv;

                this->pressureTimesHydrocarbonVolume_[globalDofIdx] =
                    this->pressureTimesPoreVolume_[globalDofIdx] * hydrocarbon;
            }
            else if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                this->pressureTimesPoreVolume_[globalDofIdx] =
                    getValue(fs.pressure(gasPhaseIdx)) * pv;

                this->pressureTimesHydrocarbonVolume_[globalDofIdx] =
                    this->pressureTimesPoreVolume_[globalDofIdx] * hydrocarbon;
            }
            else if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
                this->pressureTimesPoreVolume_[globalDofIdx] =
                    getValue(fs.pressure(waterPhaseIdx)) * pv;
            }
        }
    }

    void updatePhaseInplaceVolumes_(const unsigned             globalDofIdx,
                                    const IntensiveQuantities& intQuants,
                                    const double               totVolume)
    {
        std::array<Scalar, FluidSystem::numPhases> fip {};
        std::array<Scalar, FluidSystem::numPhases> fipr{}; // at reservoir condition

        const auto& fs = intQuants.fluidState();
        const auto  pv = totVolume * intQuants.porosity().value();

        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const auto b = getValue(fs.invB(phaseIdx));
            const auto s = getValue(fs.saturation(phaseIdx));

            fipr[phaseIdx] = s * pv;
            fip [phaseIdx] = b * fipr[phaseIdx];
        }

        this->fipC_.assignVolumesSurface(globalDofIdx, fip);
        this->fipC_.assignVolumesReservoir(globalDofIdx,
                                           fs.saltConcentration().value(),
                                           fipr);

        if (FluidSystem::phaseIsActive(oilPhaseIdx) &&
            FluidSystem::phaseIsActive(gasPhaseIdx))
        {
            this->updateOilGasDistribution(globalDofIdx, fs, fip);
        }

        if (FluidSystem::phaseIsActive(waterPhaseIdx) &&
            FluidSystem::phaseIsActive(gasPhaseIdx))
        {
            this->updateGasWaterDistribution(globalDofIdx, fs, fip);
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx) &&
            this->fipC_.hasCo2InGas())
        {
            this->updateCO2InGas(globalDofIdx, pv, intQuants);
        }

        if (this->fipC_.hasCo2InWater() &&
            (FluidSystem::phaseIsActive(waterPhaseIdx) ||
             FluidSystem::phaseIsActive(oilPhaseIdx)))
        {
            this->updateCO2InWater(globalDofIdx, pv, fs);
        }

        if constexpr(enableMICP) {
            const auto surfVolWat = pv * getValue(fs.invB(waterPhaseIdx));
            if (this->fipC_.hasMicrobialMass()) {
                this->updateMicrobialMass(globalDofIdx, intQuants, surfVolWat);
            }
            if (this->fipC_.hasOxygenMass()) {
                this->updateOxygenMass(globalDofIdx, intQuants, surfVolWat);
            }
            if (this->fipC_.hasUreaMass()) {
                this->updateUreaMass(globalDofIdx, intQuants, surfVolWat);
            }
            if (this->fipC_.hasBiofilmMass()) {
                this->updateBiofilmMass(globalDofIdx, intQuants, totVolume);
            }
            if (this->fipC_.hasCalciteMass()) {
                this->updateCalciteMass(globalDofIdx, intQuants, totVolume);
            }
        }

        if (this->fipC_.hasWaterMass() && FluidSystem::phaseIsActive(waterPhaseIdx))
        {
            this->updateWaterMass(globalDofIdx, fs, fip);
        }
    }

    template <typename FluidState, typename FIPArray>
    void updateOilGasDistribution(const unsigned    globalDofIdx,
                                  const FluidState& fs,
                                  const FIPArray&   fip)
    {
        // Gas dissolved in oil and vaporized oil
        const auto gasInPlaceLiquid = getValue(fs.Rs()) * fip[oilPhaseIdx];
        const auto oilInPlaceGas    = getValue(fs.Rv()) * fip[gasPhaseIdx];

        this->fipC_.assignOilGasDistribution(globalDofIdx, gasInPlaceLiquid, oilInPlaceGas);
    }

    template <typename FluidState, typename FIPArray>
    void updateGasWaterDistribution(const unsigned    globalDofIdx,
                                    const FluidState& fs,
                                    const FIPArray&   fip)
    {
        // Gas dissolved in water and vaporized water
        const auto gasInPlaceWater = getValue(fs.Rsw()) * fip[waterPhaseIdx];
        const auto waterInPlaceGas = getValue(fs.Rvw()) * fip[gasPhaseIdx];

        this->fipC_.assignGasWater(globalDofIdx, fip, gasInPlaceWater, waterInPlaceGas);
    }

    template <typename IntensiveQuantities>
    void updateCO2InGas(const unsigned    globalDofIdx,
                        const double      pv,
                        const IntensiveQuantities& intQuants)
    {
        const auto& scaledDrainageInfo = this->simulator_.problem().materialLawManager()
            ->oilWaterScaledEpsInfoDrainage(globalDofIdx);

        const auto& fs = intQuants.fluidState();
        Scalar sgcr = scaledDrainageInfo.Sgcr;
        if (this->simulator_.problem().materialLawManager()->enableHysteresis()) {
            const auto& matParams = simulator_.problem().materialLawParams(globalDofIdx);
            sgcr = MaterialLaw::trappedGasSaturation(matParams, /*maximumTrapping*/false);
        }

        Scalar trappedGasSaturation = scaledDrainageInfo.Sgcr;
        if (this->fipC_.has(Inplace::Phase::CO2MassInGasPhaseMaximumTrapped) ||
            this->fipC_.has(Inplace::Phase::CO2MassInGasPhaseMaximumUnTrapped))
        {
            if (this->simulator_.problem().materialLawManager()->enableHysteresis()) {
                const auto& matParams = simulator_.problem().materialLawParams(globalDofIdx);
                // Get the maximum trapped gas saturation
                trappedGasSaturation = MaterialLaw::trappedGasSaturation(matParams, /*maximumTrapping*/true);
            }
        }

        const Scalar sg = getValue(fs.saturation(gasPhaseIdx));
        Scalar strandedGasSaturation = scaledDrainageInfo.Sgcr;
        if (this->fipC_.has(Inplace::Phase::CO2MassInGasPhaseEffectiveTrapped) ||
            this->fipC_.has(Inplace::Phase::CO2MassInGasPhaseEffectiveUnTrapped))
        {
            if (this->simulator_.problem().materialLawManager()->enableHysteresis()) {
                const auto& matParams = simulator_.problem().materialLawParams(globalDofIdx);
                const double krg = getValue(intQuants.relativePermeability(gasPhaseIdx));
                strandedGasSaturation = MaterialLaw::strandedGasSaturation(matParams, sg, krg);
            }
        }

        const typename FIPContainer<FluidSystem>::Co2InGasInput v{
            pv,
            sg,
            sgcr,
            getValue(fs.density(gasPhaseIdx)),
            FluidSystem::phaseIsActive(waterPhaseIdx)
              ? FluidSystem::convertRvwToXgW(getValue(fs.Rvw()), fs.pvtRegionIndex())
              : FluidSystem::convertRvToXgO(getValue(fs.Rv()), fs.pvtRegionIndex()),
            FluidSystem::molarMass(gasCompIdx, fs.pvtRegionIndex()),
            trappedGasSaturation,
            strandedGasSaturation,
        };

        this->fipC_.assignCo2InGas(globalDofIdx, v);
    }

    template <typename FluidState>
    void updateCO2InWater(const unsigned    globalDofIdx,
                          const double      pv,
                          const FluidState& fs)
    {
        const auto co2InWater = FluidSystem::phaseIsActive(oilPhaseIdx)
            ? this->co2InWaterFromOil(fs, pv)
            : this->co2InWaterFromWater(fs, pv);

        const Scalar mM = FluidSystem::molarMass(gasCompIdx, fs.pvtRegionIndex());

        this->fipC_.assignCo2InWater(globalDofIdx, co2InWater, mM);
    }

    template <typename FluidState>
    Scalar co2InWaterFromWater(const FluidState& fs, const double pv) const
    {
        const double rhow = getValue(fs.density(waterPhaseIdx));
        const double sw   = getValue(fs.saturation(waterPhaseIdx));
        const double xwG  = FluidSystem::convertRswToXwG(getValue(fs.Rsw()), fs.pvtRegionIndex());

        const Scalar mM = FluidSystem::molarMass(gasCompIdx, fs.pvtRegionIndex());

        return xwG * pv * rhow * sw / mM;
    }

    template <typename FluidState>
    Scalar co2InWaterFromOil(const FluidState& fs, const double pv) const
    {
        const double rhoo = getValue(fs.density(oilPhaseIdx));
        const double so   = getValue(fs.saturation(oilPhaseIdx));
        const double xoG  = FluidSystem::convertRsToXoG(getValue(fs.Rs()), fs.pvtRegionIndex());

        const Scalar mM = FluidSystem::molarMass(gasCompIdx, fs.pvtRegionIndex());

        return xoG * pv * rhoo * so / mM;
    }

    template <typename FluidState, typename FIPArray>
    void updateWaterMass(const unsigned    globalDofIdx,
                         const FluidState& fs,
                         const FIPArray&   fip
                         )
    {
        const Scalar rhoW = FluidSystem::referenceDensity(waterPhaseIdx, fs.pvtRegionIndex());

        this->fipC_.assignWaterMass(globalDofIdx, fip, rhoW);
    }   

    template <typename IntensiveQuantities>
    void updateMicrobialMass(const unsigned             globalDofIdx,
                             const IntensiveQuantities& intQuants,
                             const double               surfVolWat)
    {
        const Scalar mass = surfVolWat * intQuants.microbialConcentration().value();

        this->fipC_.assignMicrobialMass(globalDofIdx, mass);
    }

    template <typename IntensiveQuantities>
    void updateOxygenMass(const unsigned             globalDofIdx,
                          const IntensiveQuantities& intQuants,
                          const double               surfVolWat)
    {
        const Scalar mass = surfVolWat * intQuants.oxygenConcentration().value();

        this->fipC_.assignOxygenMass(globalDofIdx, mass);
    }

    template <typename IntensiveQuantities>
    void updateUreaMass(const unsigned             globalDofIdx,
                        const IntensiveQuantities& intQuants,
                        const double               surfVolWat)
    {
        const Scalar mass = surfVolWat * intQuants.ureaConcentration().value();

        this->fipC_.assignUreaMass(globalDofIdx, mass);
    }

    template <typename IntensiveQuantities>
    void updateBiofilmMass(const unsigned             globalDofIdx,
                           const IntensiveQuantities& intQuants,
                           const double               totVolume)
    {
        const Scalar mass = totVolume * intQuants.biofilmMass().value();

        this->fipC_.assignBiofilmMass(globalDofIdx, mass);
    }

    template <typename IntensiveQuantities>
    void updateCalciteMass(const unsigned             globalDofIdx,
                           const IntensiveQuantities& intQuants,
                           const double               totVolume)
    {
        const Scalar mass = totVolume * intQuants.calciteMass().value();

        this->fipC_.assignCalciteMass(globalDofIdx, mass);
    }

    //! \brief Setup extractors for element-level data.
    void setupElementExtractors_()
    {
        using Entry = typename Extractor::Entry;
        using Context = typename Extractor::Context;
        using ScalarEntry = typename Extractor::ScalarEntry;
        using PhaseEntry = typename Extractor::PhaseEntry;

        const bool hasResidual = simulator_.model().linearizer().residual().size() > 0;
        const auto& hysteresisConfig = simulator_.problem().materialLawManager()->hysteresisConfig();

        auto extractors = std::array{
            Entry{PhaseEntry{&this->saturation_,
                             [](const unsigned phase, const Context& ectx)
                             { return getValue(ectx.fs.saturation(phase)); }
                  }
            },
            Entry{PhaseEntry{&this->invB_,
                            [](const unsigned phase, const Context& ectx)
                            { return getValue(ectx.fs.invB(phase)); }
                  }
            },
            Entry{PhaseEntry{&this->density_,
                             [](const unsigned phase, const Context& ectx)
                             { return getValue(ectx.fs.density(phase)); }
                  }
            },
            Entry{PhaseEntry{&this->relativePermeability_,
                            [](const unsigned phase, const Context& ectx)
                            { return getValue(ectx.intQuants.relativePermeability(phase)); }
                  }
            },
            Entry{PhaseEntry{&this->viscosity_,
                             [this](const unsigned phaseIdx, const Context& ectx)
                             {
                                if (this->extboC_.allocated() && phaseIdx == oilPhaseIdx) {
                                    return getValue(ectx.intQuants.oilViscosity());
                                }
                                else if (this->extboC_.allocated() && phaseIdx == gasPhaseIdx) {
                                    return getValue(ectx.intQuants.gasViscosity());
                                }
                                else {
                                    return getValue(ectx.fs.viscosity(phaseIdx));
                                }
                             }
                  }
            },
            Entry{PhaseEntry{&this->residual_,
                             [&modelResid = this->simulator_.model().linearizer().residual()]
                             (const unsigned phaseIdx, const Context& ectx)
                             {
                                const unsigned sIdx = FluidSystem::solventComponentIndex(phaseIdx);
                                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(sIdx);
                                return modelResid[ectx.globalDofIdx][activeCompIdx];
                             }
                  },
                  hasResidual
            },
            Entry{ScalarEntry{&this->rockCompPorvMultiplier_,
                             [&problem = this->simulator_.problem()](const Context& ectx)
                             {
                                return problem.template
                                    rockCompPoroMultiplier<Scalar>(ectx.intQuants,
                                                                   ectx.globalDofIdx);
                             }
                  }
            },
            Entry{ScalarEntry{&this->rockCompTransMultiplier_,
                             [&problem = this->simulator_.problem()](const Context& ectx)
                             {
                                return problem.
                                    template rockCompTransMultiplier<Scalar>(ectx.intQuants,
                                                                             ectx.globalDofIdx);
                  }}
            },
            Entry{ScalarEntry{&this->minimumOilPressure_,
                              [&problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  return std::min(getValue(ectx.fs.pressure(oilPhaseIdx)),
                                                           problem.minOilPressure(ectx.globalDofIdx));
                              }
                  }
            },
            Entry{ScalarEntry{&this->bubblePointPressure_,
                             [&failedCells = this->failedCellsPb_,
                              &vanguard = this->simulator_.vanguard()](const Context& ectx)
                             {
                                try {
                                    return getValue(
                                              FluidSystem::bubblePointPressure(ectx.fs,
                                                                               ectx.intQuants.pvtRegionIndex())
                                           );
                                } catch (const NumericalProblem&) {
                                    const auto cartesianIdx = vanguard.cartesianIndex(ectx.globalDofIdx);
                                    failedCells.push_back(cartesianIdx);
                                    return Scalar{0};
                                }
                             }
                  }
            },
            Entry{ScalarEntry{&this->dewPointPressure_,
                              [&failedCells = this->failedCellsPd_,
                               &vanguard = this->simulator_.vanguard()](const Context& ectx)
                              {
                                  try {
                                      return getValue(
                                          FluidSystem::dewPointPressure(ectx.fs,
                                                                        ectx.intQuants.pvtRegionIndex())
                                      );
                                  } catch (const NumericalProblem&) {
                                      const auto cartesianIdx =  vanguard.cartesianIndex(ectx.globalDofIdx);
                                      failedCells.push_back(cartesianIdx);
                                      return Scalar{0};
                                  }
                              }
                  }
            },
            Entry{ScalarEntry{&this->overburdenPressure_,
                             [&problem = simulator_.problem()](const Context& ectx)
                             { return problem.overburdenPressure(ectx.globalDofIdx); }
                  }
            },
            Entry{ScalarEntry{&this->temperature_,
                              [](const Context& ectx)
                              { return getValue(ectx.fs.temperature(oilPhaseIdx)); }
                  }
            },
            Entry{ScalarEntry{&this->sSol_,
                              [](const Context& ectx)
                              { return getValue(ectx.intQuants.solventSaturation()); }
                  }
            },
            Entry{ScalarEntry{&this->rswSol_,
                              [](const Context& ectx)
                              { return getValue(ectx.intQuants.rsSolw()); }
                  }
            },
            Entry{ScalarEntry{&this->cPolymer_,
                              [](const Context& ectx)
                              { return getValue(ectx.intQuants.polymerConcentration()); }
                  }
            },
            Entry{ScalarEntry{&this->cFoam_,
                              [](const Context& ectx)
                              { return getValue(ectx.intQuants.foamConcentration()); }
                  }
            },
            Entry{ScalarEntry{&this->cSalt_,
                              [](const Context& ectx)
                              { return getValue(ectx.fs.saltConcentration()); }
                  }
            },
            Entry{ScalarEntry{&this->pSalt_,
                              [](const Context& ectx)
                              { return getValue(ectx.fs.saltSaturation()); }
                  }
            },
            Entry{ScalarEntry{&this->permFact_,
                              [](const Context& ectx)
                              { return getValue(ectx.intQuants.permFactor()); }
                  }
            },
            Entry{ScalarEntry{&this->rPorV_,
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  const auto totVolume = model.dofTotalVolume(ectx.globalDofIdx);
                                  return totVolume * getValue(ectx.intQuants.porosity());
                              }
                  }
            },
            Entry{ScalarEntry{&this->rs_,
                              [](const Context& ectx)
                              { return getValue(ectx.fs.Rs()); }
                  }
            },
            Entry{ScalarEntry{&this->rv_,
                              [](const Context& ectx)
                              { return getValue(ectx.fs.Rv()); }
                  }
            },
            Entry{ScalarEntry{&this->rsw_,
                              [](const Context& ectx)
                              { return getValue(ectx.fs.Rsw()); }
                  }
            },
            Entry{ScalarEntry{&this->rvw_,
                              [](const Context& ectx)
                              { return getValue(ectx.fs.Rvw()); }
                  }
            },
            Entry{ScalarEntry{&this->ppcw_,
                              [&matLawManager = *this->simulator_.problem().materialLawManager()]
                              (const Context& ectx)
                              {
                                  return matLawManager.
                                      oilWaterScaledEpsInfoDrainage(ectx.globalDofIdx).maxPcow;
                              }
                  }
            },
            Entry{ScalarEntry{&this->drsdtcon_,
                              [&problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  return problem.drsdtcon(ectx.globalDofIdx,
                                                          ectx.episodeIndex);
                              }
                  }
            },
            Entry{ScalarEntry{&this->pcgw_,
                              [](const Context& ectx)
                              {
                                  return getValue(ectx.fs.pressure(gasPhaseIdx)) -
                                         getValue(ectx.fs.pressure(waterPhaseIdx));
                              }
                  }
            },
            Entry{ScalarEntry{&this->pcow_,
                              [](const Context& ectx)
                              {
                                  return getValue(ectx.fs.pressure(oilPhaseIdx)) -
                                         getValue(ectx.fs.pressure(waterPhaseIdx));
                              }
                  }
            },
            Entry{ScalarEntry{&this->pcog_,
                              [](const Context& ectx)
                              {
                                  return getValue(ectx.fs.pressure(gasPhaseIdx)) -
                                         getValue(ectx.fs.pressure(oilPhaseIdx));
                              }
                  }
            },
            Entry{ScalarEntry{&this->fluidPressure_,
                              [](const Context& ectx)
                              {
                                  if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                                      // Output oil pressure as default
                                      return getValue(ectx.fs.pressure(oilPhaseIdx));
                                  }
                                  else if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                                      // Output gas if oil is not present
                                      return getValue(ectx.fs.pressure(gasPhaseIdx));
                                  }
                                  else {
                                      // Output water if neither oil nor gas is present
                                      return getValue(ectx.fs.pressure(waterPhaseIdx));
                                  }
                              }
                  }
            },
            Entry{ScalarEntry{&this->gasDissolutionFactor_,
                              [&problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const Scalar SoMax = problem.maxOilSaturation(ectx.globalDofIdx);
                                  return FluidSystem::template
                                      saturatedDissolutionFactor<FluidState, Scalar>(ectx.fs,
                                                                                     oilPhaseIdx,
                                                                                     ectx.pvtRegionIdx,
                                                                                     SoMax);
                              }
                  }
            },
            Entry{ScalarEntry{&this->oilVaporizationFactor_,
                              [&problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const Scalar SoMax = problem.maxOilSaturation(ectx.globalDofIdx);
                                  return FluidSystem::template
                                      saturatedDissolutionFactor<FluidState, Scalar>(ectx.fs,
                                                                                     gasPhaseIdx,
                                                                                     ectx.pvtRegionIdx,
                                                                                     SoMax);
                              }
                  }
            },
            Entry{ScalarEntry{&this->gasDissolutionFactorInWater_,
                              [&problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const Scalar SwMax = problem.maxWaterSaturation(ectx.globalDofIdx);
                                  return FluidSystem::template
                                      saturatedDissolutionFactor<FluidState, Scalar>(ectx.fs,
                                                                                     waterPhaseIdx,
                                                                                     ectx.pvtRegionIdx,
                                                                                     SwMax);
                              }
                  }
            },
            Entry{ScalarEntry{&this->waterVaporizationFactor_,
                              [](const Context& ectx)
                              {
                                  return FluidSystem::template
                                      saturatedVaporizationFactor<FluidState, Scalar>(ectx.fs,
                                                                                      gasPhaseIdx,
                                                                                      ectx.pvtRegionIdx);
                              }
                  }
            },
            Entry{ScalarEntry{&this->gasFormationVolumeFactor_,
                              [](const Context& ectx)
                              {
                                  return 1.0 / FluidSystem::template
                                                   inverseFormationVolumeFactor<FluidState, Scalar>(ectx.fs,
                                                                                                    gasPhaseIdx,
                                                                                                    ectx.pvtRegionIdx);
                              }
                  }
            },
            Entry{ScalarEntry{&this->saturatedOilFormationVolumeFactor_,
                              [](const Context& ectx)
                              {
                                  return 1.0 / FluidSystem::template
                                             saturatedInverseFormationVolumeFactor<FluidState, Scalar>(ectx.fs,
                                                                                                       oilPhaseIdx,
                                                                                                       ectx.pvtRegionIdx);
                              }
                  }
            },
            Entry{ScalarEntry{&this->oilSaturationPressure_,
                              [](const Context& ectx)
                              {
                                  return FluidSystem::template
                                      saturationPressure<FluidState, Scalar>(ectx.fs,
                                                                             oilPhaseIdx,
                                                                             ectx.pvtRegionIdx);
                              }
                  }
            },
            Entry{ScalarEntry{&this->soMax_,
                             [&problem = this->simulator_.problem()](const Context& ectx)
                             {
                                 return std::max(getValue(ectx.fs.saturation(oilPhaseIdx)),
                                                 problem.maxOilSaturation(ectx.globalDofIdx));
                             }
                  },
                  !hysteresisConfig.enableHysteresis()
            },
            Entry{ScalarEntry{&this->swMax_,
                             [&problem = this->simulator_.problem()](const Context& ectx)
                             {
                                 return std::max(getValue(ectx.fs.saturation(waterPhaseIdx)),
                                                 problem.maxWaterSaturation(ectx.globalDofIdx));
                             }
                  },
                  !hysteresisConfig.enableHysteresis()
            },
            Entry{ScalarEntry{&this->soMax_,
                              [](const Context& ectx)
                              { return ectx.hParams.somax; }
                  },
                  hysteresisConfig.enableHysteresis() &&
                  hysteresisConfig.enableNonWettingHysteresis() &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(waterPhaseIdx)
            },
            Entry{ScalarEntry{&this->swMax_,
                              [](const Context& ectx)
                              { return ectx.hParams.swmax; }
                  },
                  hysteresisConfig.enableHysteresis() &&
                  hysteresisConfig.enableWettingHysteresis() &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(waterPhaseIdx)
            },
            Entry{ScalarEntry{&this->swmin_,
                             [](const Context& ectx)
                             { return ectx.hParams.swmin; }
                  },
                  hysteresisConfig.enableHysteresis() &&
                  hysteresisConfig.enablePCHysteresis() &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(waterPhaseIdx)
            },
            Entry{ScalarEntry{&this->sgmax_,
                              [](const Context& ectx)
                              { return ectx.hParams.sgmax; }
                  },
                  hysteresisConfig.enableHysteresis() &&
                  hysteresisConfig.enableNonWettingHysteresis() &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(gasPhaseIdx)
            },
            Entry{ScalarEntry{&this->shmax_,
                              [](const Context& ectx)
                              { return ectx.hParams.shmax; }
                  },
                  hysteresisConfig.enableHysteresis() &&
                  hysteresisConfig.enableWettingHysteresis() &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(gasPhaseIdx)
            },
            Entry{ScalarEntry{&this->somin_,
                              [](const Context& ectx)
                              { return ectx.hParams.somin; }
                  },
                  hysteresisConfig.enableHysteresis() &&
                  hysteresisConfig.enablePCHysteresis() &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(gasPhaseIdx)
            },
            Entry{[&model = this->simulator_.model(), this](const Context& ectx)
                  {
                      // Note: We intentionally exclude effects of rock
                      // compressibility by using referencePorosity() here.
                      const auto porv = ectx.intQuants.referencePorosity()
                          * model.dofTotalVolume(ectx.globalDofIdx);

                      this->aggregateAverageDensityContributions_(ectx.fs, ectx.globalDofIdx,
                                                                  static_cast<double>(porv));
                  }, this->regionAvgDensity_.has_value()
            },
            Entry{[&extboC = this->extboC_](const Context& ectx)
                  {
                      extboC.assignVolumes(ectx.globalDofIdx,
                                           ectx.intQuants.xVolume().value(),
                                           ectx.intQuants.yVolume().value());
                      extboC.assignZFraction(ectx.globalDofIdx,
                                             ectx.intQuants.zFraction().value());

                      const Scalar stdVolOil = getValue(ectx.fs.saturation(oilPhaseIdx)) *
                                               getValue(ectx.fs.invB(oilPhaseIdx)) +
                                               getValue(ectx.fs.saturation(gasPhaseIdx)) *
                                               getValue(ectx.fs.invB(gasPhaseIdx)) *
                                               getValue(ectx.fs.Rv());
                      const Scalar stdVolGas = getValue(ectx.fs.saturation(gasPhaseIdx)) *
                                               getValue(ectx.fs.invB(gasPhaseIdx)) *
                                               (1.0 - ectx.intQuants.yVolume().value()) +
                                               getValue(ectx.fs.saturation(oilPhaseIdx)) *
                                               getValue(ectx.fs.invB(oilPhaseIdx)) *
                                               getValue(ectx.fs.Rs()) *
                                               (1.0 - ectx.intQuants.xVolume().value());
                      const Scalar stdVolCo2 = getValue(ectx.fs.saturation(gasPhaseIdx)) *
                                               getValue(ectx.fs.invB(gasPhaseIdx)) *
                                               ectx.intQuants.yVolume().value() +
                                               getValue(ectx.fs.saturation(oilPhaseIdx)) *
                                               getValue(ectx.fs.invB(oilPhaseIdx)) *
                                               getValue(ectx.fs.Rs()) *
                                               ectx.intQuants.xVolume().value();
                      const Scalar rhoO = FluidSystem::referenceDensity(gasPhaseIdx, ectx.pvtRegionIdx);
                      const Scalar rhoG = FluidSystem::referenceDensity(gasPhaseIdx, ectx.pvtRegionIdx);
                      const Scalar rhoCO2 = ectx.intQuants.zRefDensity();
                      const Scalar stdMassTotal = 1.0e-10 + stdVolOil * rhoO + stdVolGas * rhoG + stdVolCo2 * rhoCO2;
                      extboC.assignMassFractions(ectx.globalDofIdx,
                                                 stdVolGas * rhoG / stdMassTotal,
                                                 stdVolOil * rhoO / stdMassTotal,
                                                 stdVolCo2 * rhoCO2 / stdMassTotal);
                    }, this->extboC_.allocated()
            },
            Entry{[&micpC = this->micpC_](const Context& ectx)
                  {
                      micpC.assign(ectx.globalDofIdx,
                                   ectx.intQuants.microbialConcentration().value(),
                                   ectx.intQuants.oxygenConcentration().value(),
                                   ectx.intQuants.ureaConcentration().value(),
                                   ectx.intQuants.biofilmConcentration().value(),
                                   ectx.intQuants.calciteConcentration().value());
                  }, this->micpC_.allocated()
            },
            Entry{[&rftC = this->rftC_,
                   &vanguard = this->simulator_.vanguard()](const Context& ectx)
                  {
                      const auto cartesianIdx = vanguard.cartesianIndex(ectx.globalDofIdx);
                      rftC.assign(cartesianIdx,
                                  [&fs = ectx.fs]() { return getValue(fs.pressure(oilPhaseIdx)); },
                                  [&fs = ectx.fs]() { return getValue(fs.saturation(waterPhaseIdx)); },
                                  [&fs = ectx.fs]() { return getValue(fs.saturation(gasPhaseIdx)); });
                   }
            },
            Entry{[&tC = this->tracerC_,
                   &tM = this->simulator_.problem().tracerModel()](const Context& ectx)
                  {
                      tC.assignFreeConcentrations(ectx.globalDofIdx,
                                                  [gIdx = ectx.globalDofIdx, &tM](const unsigned tracerIdx)
                                                      { return tM.freeTracerConcentration(tracerIdx, gIdx); });
                      tC.assignSolConcentrations(ectx.globalDofIdx,
                                                 [gIdx = ectx.globalDofIdx, &tM](const unsigned tracerIdx)
                                                 { return tM.solTracerConcentration(tracerIdx, gIdx); });
                    }
            },
            Entry{[&flowsInf = this->simulator_.problem().model().linearizer().getFlowsInfo(),
                   &flowsC = this->flowsC_](const Context& ectx)
                  {
                      constexpr auto gas_idx = Indices::gasEnabled ?
                          conti0EqIdx + Indices::canonicalToActiveComponentIndex(gasCompIdx) : -1;
                      constexpr auto oil_idx = Indices::oilEnabled ?
                          conti0EqIdx + Indices::canonicalToActiveComponentIndex(oilCompIdx) : -1;
                      constexpr auto water_idx = Indices::waterEnabled ?
                          conti0EqIdx + Indices::canonicalToActiveComponentIndex(waterCompIdx) : -1;
                      const auto& flowsInfos = flowsInf[ectx.globalDofIdx];
                      for (const auto& flowsInfo : flowsInfos) {
                          flowsC.assignFlows(ectx.globalDofIdx,
                                             flowsInfo.faceId,
                                             flowsInfo.nncId,
                                             value_or_zero<gas_idx>(flowsInfo.flow),
                                             value_or_zero<oil_idx>(flowsInfo.flow),
                                             value_or_zero<water_idx>(flowsInfo.flow));
                        }
                 }, !this->simulator_.problem().model().linearizer().getFlowsInfo().empty()
            },
            Entry{[&floresInf = this->simulator_.problem().model().linearizer().getFloresInfo(),
                   &flowsC = this->flowsC_](const Context& ectx)
                  {
                      constexpr auto gas_idx = Indices::gasEnabled ?
                          conti0EqIdx + Indices::canonicalToActiveComponentIndex(gasCompIdx) : -1;
                      constexpr auto oil_idx = Indices::oilEnabled ?
                          conti0EqIdx + Indices::canonicalToActiveComponentIndex(oilCompIdx) : -1;
                      constexpr auto water_idx = Indices::waterEnabled ?
                          conti0EqIdx + Indices::canonicalToActiveComponentIndex(waterCompIdx) : -1;
                      const auto& floresInfos = floresInf[ectx.globalDofIdx];
                      for (const auto& floresInfo : floresInfos) {
                          flowsC.assignFlores(ectx.globalDofIdx,
                                              floresInfo.faceId,
                                              floresInfo.nncId,
                                              value_or_zero<gas_idx>(floresInfo.flow),
                                              value_or_zero<oil_idx>(floresInfo.flow),
                                              value_or_zero<water_idx>(floresInfo.flow));
                      }
                 }, !this->simulator_.problem().model().linearizer().getFloresInfo().empty()
            },
            // hack to make the intial output of rs and rv Ecl compatible.
            // For cells with swat == 1 Ecl outputs; rs = rsSat and rv=rvSat, in all but the initial step
            // where it outputs rs and rv values calculated by the initialization. To be compatible we overwrite
            // rs and rv with the values computed in the initially.
            // Volume factors, densities and viscosities need to be recalculated with the updated rs and rv values.
            Entry{ScalarEntry{&this->rv_,
                             [&problem = this->simulator_.problem()](const Context& ectx)
                             { return problem.initialFluidState(ectx.globalDofIdx).Rv(); }
                  },
                  simulator_.episodeIndex() < 0 &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(gasPhaseIdx)
            },
            Entry{ScalarEntry{&this->rs_,
                              [&problem = this->simulator_.problem()](const Context& ectx)
                              { return problem.initialFluidState(ectx.globalDofIdx).Rs(); }
                  },
                  simulator_.episodeIndex() < 0 &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(gasPhaseIdx)
            },
            Entry{ScalarEntry{&this->rsw_,
                              [&problem = this->simulator_.problem()](const Context& ectx)
                              { return problem.initialFluidState(ectx.globalDofIdx).Rsw(); }
                  },
                  simulator_.episodeIndex() < 0 &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(gasPhaseIdx)
            },
            Entry{ScalarEntry{&this->rvw_,
                              [&problem = this->simulator_.problem()](const Context& ectx)
                              { return problem.initialFluidState(ectx.globalDofIdx).Rvw(); }
                  },
                  simulator_.episodeIndex() < 0 &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(gasPhaseIdx)
            },
            // re-compute the volume factors, viscosities and densities if asked for
            Entry{PhaseEntry{&this->density_,
                            [&problem = this->simulator_.problem()](const unsigned phase,
                                                                    const Context& ectx)
                            {
                                const auto& fsInitial = problem.initialFluidState(ectx.globalDofIdx);
                                return FluidSystem::density(fsInitial,
                                                            phase,
                                                            ectx.intQuants.pvtRegionIndex());
                            }
                  },
                  simulator_.episodeIndex() < 0 &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(gasPhaseIdx)
            },
            Entry{PhaseEntry{&this->invB_,
                            [&problem = this->simulator_.problem()](const unsigned phase,
                                                                    const Context& ectx)
                            {
                                const auto& fsInitial = problem.initialFluidState(ectx.globalDofIdx);
                                return FluidSystem::inverseFormationVolumeFactor(fsInitial,
                                                                                 phase,
                                                                                 ectx.intQuants.pvtRegionIndex());
                            }
                  },
                  simulator_.episodeIndex() < 0 &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(gasPhaseIdx)
            },
            Entry{PhaseEntry{&this->viscosity_,
                            [&problem = this->simulator_.problem()](const unsigned phase,
                                                                    const Context& ectx)
                            {
                                const auto& fsInitial = problem.initialFluidState(ectx.globalDofIdx);
                                return FluidSystem::viscosity(fsInitial,
                                                              phase,
                                                              ectx.intQuants.pvtRegionIndex());
                            }
                  },
                  simulator_.episodeIndex() < 0 &&
                  FluidSystem::phaseIsActive(oilPhaseIdx) &&
                  FluidSystem::phaseIsActive(gasPhaseIdx)
            },
        };

        // Setup active extractors
        this->extractors_ = Extractor::removeInactive(extractors);

        if constexpr (getPropValue<TypeTag, Properties::EnableMech>()) {
            if (this->mech_.allocated()) {
                this->extractors_.push_back(
                    Entry{[&mech = this->mech_,
                           &model = simulator_.problem().geoMechModel()](const Context& ectx)
                           {
                              mech.assignDelStress(ectx.globalDofIdx,
                                                   model.delstress(ectx.globalDofIdx));

                              mech.assignDisplacement(ectx.globalDofIdx,
                                                      model.disp(ectx.globalDofIdx, /*include_fracture*/true));

                              // is the tresagii stress which make rock fracture
                              mech.assignFracStress(ectx.globalDofIdx,
                                                    model.fractureStress(ectx.globalDofIdx));

                              mech.assignLinStress(ectx.globalDofIdx,
                                                   model.linstress(ectx.globalDofIdx));

                              mech.assignPotentialForces(ectx.globalDofIdx,
                                                         model.mechPotentialForce(ectx.globalDofIdx),
                                                         model.mechPotentialPressForce(ectx.globalDofIdx),
                                                         model.mechPotentialTempForce(ectx.globalDofIdx));

                              mech.assignStrain(ectx.globalDofIdx,
                                                model.strain(ectx.globalDofIdx, /*include_fracture*/true));

                              // Total stress is not stored but calculated result is Voigt notation
                              mech.assignStress(ectx.globalDofIdx,
                                                model.stress(ectx.globalDofIdx, /*include_fracture*/true));
                           }
                    }
                );
            }
        }
    }

    //! \brief Setup extractor execution map for block data.
    void setupBlockExtractors_(const bool isSubStep,
                               const int  reportStepNum)
    {
        using Entry = typename BlockExtractor::Entry;
        using Context = typename BlockExtractor::Context;
        using PhaseEntry = typename BlockExtractor::PhaseEntry;
        using ScalarEntry = typename BlockExtractor::ScalarEntry;

        using namespace std::string_view_literals;

        const auto pressure_handler =
            Entry{ScalarEntry{std::vector{"BPR"sv, "BPRESSUR"sv},
                              [](const Context& ectx)
                              {
                                  if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                                      return getValue(ectx.fs.pressure(oilPhaseIdx));
                                  }
                                  else if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                                      return getValue(ectx.fs.pressure(gasPhaseIdx));
                                  }
                                  else { //if (FluidSystem::phaseIsActive(waterPhaseIdx))
                                      return getValue(ectx.fs.pressure(waterPhaseIdx));
                                  }
                              }
                  }
            };

        const auto handlers = std::array{
            pressure_handler,
            Entry{PhaseEntry{std::array{
                                std::array{"BWSAT"sv, "BOSAT"sv, "BGSAT"sv},
                                std::array{"BSWAT"sv, "BSOIL"sv, "BSGAS"sv}
                             },
                             [](const unsigned phaseIdx, const Context& ectx)
                             {
                                 return getValue(ectx.fs.saturation(phaseIdx));
                             }
                  }
            },
            Entry{ScalarEntry{"BNSAT",
                              [](const Context& ectx)
                              {
                                  return ectx.intQuants.solventSaturation().value();
                              }
                  }
            },
            Entry{ScalarEntry{std::vector{"BTCNFHEA"sv, "BTEMP"sv},
                              [](const Context& ectx)
                              {
                                  if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                                      return getValue(ectx.fs.temperature(oilPhaseIdx));
                                  }
                                  else if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                                      return getValue(ectx.fs.temperature(gasPhaseIdx));
                                  }
                                  else { //if (FluidSystem::phaseIsActive(waterPhaseIdx))
                                      return getValue(ectx.fs.temperature(waterPhaseIdx));
                                  }
                              }
                  }
            },
            Entry{PhaseEntry{std::array{
                                std::array{"BWKR"sv, "BOKR"sv, "BGKR"sv},
                                std::array{"BKRW"sv, "BKRO"sv, "BKRG"sv}
                             },
                             [](const unsigned phaseIdx, const Context& ectx)
                             {
                                 return getValue(ectx.intQuants.relativePermeability(phaseIdx));
                             }
                  }
            },
            Entry{ScalarEntry{"BKROG",
                              [&problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const auto& materialParams =
                                      problem.materialLawParams(ectx.elemCtx,
                                                                ectx.dofIdx,
                                                                /* timeIdx = */ 0);
                                  return getValue(MaterialLaw::template
                                                    relpermOilInOilGasSystem<Evaluation>(materialParams,
                                                                                         ectx.fs));
                              }
                 }
            },
            Entry{ScalarEntry{"BKROW",
                             [&problem = this->simulator_.problem()](const Context& ectx)
                             {
                                 const auto& materialParams = problem.materialLawParams(ectx.elemCtx,
                                                                                          ectx.dofIdx,
                                                                                          /* timeIdx = */ 0);
                                 return getValue(MaterialLaw::template
                                                     relpermOilInOilWaterSystem<Evaluation>(materialParams,
                                                                                            ectx.fs));
                            }
                 }
            },
            Entry{ScalarEntry{"BWPC",
                              [](const Context& ectx)
                              {
                                  return getValue(ectx.fs.pressure(oilPhaseIdx)) -
                                         getValue(ectx.fs.pressure(waterPhaseIdx));
                              }
                  }
            },
            Entry{ScalarEntry{"BGPC",
                              [](const Context& ectx)
                              {
                                  return getValue(ectx.fs.pressure(gasPhaseIdx)) -
                                         getValue(ectx.fs.pressure(oilPhaseIdx));
                              }
                  }
            },
            Entry{ScalarEntry{"BWPR",
                              [](const Context& ectx)
                              {
                                  return getValue(ectx.fs.pressure(waterPhaseIdx));
                              }
                  }
            },
            Entry{ScalarEntry{"BGPR",
                              [](const Context& ectx)
                              {
                                  return getValue(ectx.fs.pressure(gasPhaseIdx));
                              }
                  }
            },
            Entry{PhaseEntry{std::array{
                                std::array{"BVWAT"sv, "BVOIL"sv, "BVGAS"sv},
                                std::array{"BWVIS"sv, "BOVIS"sv, "BGVIS"sv}
                             },
                             [](const unsigned phaseIdx, const Context& ectx)
                             {
                                 return getValue(ectx.fs.viscosity(phaseIdx));
                             }
                  }
            },
            Entry{PhaseEntry{std::array{
                                std::array{"BWDEN"sv, "BODEN"sv, "BGDEN"sv},
                                std::array{"BDENW"sv, "BDENO"sv, "BDENG"sv}
                             },
                             [](const unsigned phaseIdx, const Context& ectx)
                             {
                                 return getValue(ectx.fs.density(phaseIdx));
                             }
                  }
            },
            Entry{ScalarEntry{"BFLOWI",
                              [&flowsC = this->flowsC_](const Context& ectx)
                              {
                                  return flowsC.getFlow(ectx.globalDofIdx, Dir::XPlus, waterCompIdx);
                              }
                  }
            },
            Entry{ScalarEntry{"BFLOWJ",
                              [&flowsC = this->flowsC_](const Context& ectx)
                              {
                                  return flowsC.getFlow(ectx.globalDofIdx, Dir::YPlus, waterCompIdx);
                              }
                  }
            },
            Entry{ScalarEntry{"BFLOWK",
                              [&flowsC = this->flowsC_](const Context& ectx)
                              {
                                  return flowsC.getFlow(ectx.globalDofIdx, Dir::ZPlus, waterCompIdx);
                              }
                  }
            },
            Entry{ScalarEntry{"BRPV",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                   return getValue(ectx.intQuants.porosity()) *
                                          model.dofTotalVolume(ectx.globalDofIdx);
                              }
                  }
            },
            Entry{PhaseEntry{std::array{"BWPV"sv, "BOPV"sv, "BGPV"sv},
                             [&model = this->simulator_.model()](const unsigned phaseIdx,
                                                                 const Context& ectx)
                             {
                                 return getValue(ectx.fs.saturation(phaseIdx)) *
                                        getValue(ectx.intQuants.porosity()) *
                                        model.dofTotalVolume(ectx.globalDofIdx);
                             }
                  }
            },
            Entry{ScalarEntry{"BRS",
                              [](const Context& ectx)
                              {
                                  return getValue(ectx.fs.Rs());
                              }
                  }
            },
            Entry{ScalarEntry{"BRV",
                              [](const Context& ectx)
                              {
                                  return getValue(ectx.fs.Rv());
                              }
                  }
            },
            Entry{ScalarEntry{"BOIP",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  return (getValue(ectx.fs.invB(oilPhaseIdx)) *
                                          getValue(ectx.fs.saturation(oilPhaseIdx)) +
                                          getValue(ectx.fs.Rv()) *
                                          getValue(ectx.fs.invB(gasPhaseIdx)) *
                                          getValue(ectx.fs.saturation(gasPhaseIdx))) *
                                          model.dofTotalVolume(ectx.globalDofIdx) *
                                          getValue(ectx.intQuants.porosity());
                              }
                  }
            },
            Entry{ScalarEntry{"BGIP",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  Scalar result = getValue(ectx.fs.invB(gasPhaseIdx)) *
                                                  getValue(ectx.fs.saturation(gasPhaseIdx));

                                  if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                                       result += getValue(ectx.fs.Rs()) *
                                                 getValue(ectx.fs.invB(oilPhaseIdx)) *
                                                 getValue(ectx.fs.saturation(oilPhaseIdx));
                                  }
                                  else {
                                      result += getValue(ectx.fs.Rsw()) *
                                                getValue(ectx.fs.invB(waterPhaseIdx)) *
                                                getValue(ectx.fs.saturation(waterPhaseIdx));
                                  }

                                  return result *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity());
                              }
                  }
            },
            Entry{ScalarEntry{"BWIP",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  return getValue(ectx.fs.invB(waterPhaseIdx)) *
                                         getValue(ectx.fs.saturation(waterPhaseIdx)) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity());
                              }
                  }
            },
            Entry{ScalarEntry{"BOIPL",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  return getValue(ectx.fs.invB(oilPhaseIdx)) *
                                         getValue(ectx.fs.saturation(oilPhaseIdx)) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity());
                              }
                  }
            },
            Entry{ScalarEntry{"BGIPL",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  Scalar result;
                                  if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                                      result  = getValue(ectx.fs.Rs()) *
                                                getValue(ectx.fs.invB(oilPhaseIdx)) *
                                                getValue(ectx.fs.saturation(oilPhaseIdx));
                                  }
                                  else {
                                      result = getValue(ectx.fs.Rsw()) *
                                               getValue(ectx.fs.invB(waterPhaseIdx)) *
                                               getValue(ectx.fs.saturation(waterPhaseIdx));
                                  }
                                  return result *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity());
                              }
                  }
            },
            Entry{ScalarEntry{"BGIPG",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  return getValue(ectx.fs.invB(gasPhaseIdx)) *
                                         getValue(ectx.fs.saturation(gasPhaseIdx)) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity());
                              }
                  }
            },
            Entry{ScalarEntry{"BOIPG",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  return getValue(ectx.fs.Rv()) *
                                         getValue(ectx.fs.invB(gasPhaseIdx)) *
                                         getValue(ectx.fs.saturation(gasPhaseIdx)) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity());
                              }
                  }
            },
            Entry{PhaseEntry{std::array{"BPPW"sv, "BPPO"sv, "BPPG"sv},
                             [&simConfig = this->eclState_.getSimulationConfig(),
                              &grav = this->simulator_.problem().gravity(),
                              &regionAvgDensity = this->regionAvgDensity_,
                              &problem = this->simulator_.problem(),
                              &regions = this->regions_](const unsigned phaseIdx, const Context& ectx)
                             {
                                auto phase = RegionPhasePoreVolAverage::Phase{};
                                phase.ix = phaseIdx;

                                // Note different region handling here.  FIPNUM is
                                // one-based, but we need zero-based lookup in
                                // DatumDepth.  On the other hand, pvtRegionIndex is
                                // zero-based but we need one-based lookup in
                                // RegionPhasePoreVolAverage.

                                // Subtract one to convert FIPNUM to region index.
                                const auto datum = simConfig.datumDepths()(regions["FIPNUM"][ectx.dofIdx] - 1);

                                // Add one to convert region index to region ID.
                                const auto region = RegionPhasePoreVolAverage::Region {
                                    ectx.elemCtx.primaryVars(ectx.dofIdx, /*timeIdx=*/0).pvtRegionIndex() + 1
                                };

                                const auto density = regionAvgDensity->value("PVTNUM", phase, region);

                                const auto press = getValue(ectx.fs.pressure(phase.ix));
                                const auto dz = problem.dofCenterDepth(ectx.globalDofIdx) - datum;
                                return press - density*dz*grav[GridView::dimensionworld - 1];
                            }
                  }
            },
            Entry{ScalarEntry{"BAMIP",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  const Scalar rhoW = FluidSystem::referenceDensity(waterPhaseIdx,
                                                                   ectx.intQuants.pvtRegionIndex());
                                  return getValue(ectx.fs.invB(waterPhaseIdx)) *
                                         getValue(ectx.fs.saturation(waterPhaseIdx)) *
                                         rhoW *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity());
                              }
                  }
            },
            Entry{ScalarEntry{"BMMIP",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  return getValue(ectx.intQuants.microbialConcentration()) *
                                         getValue(ectx.intQuants.porosity()) *
                                         model.dofTotalVolume(ectx.globalDofIdx);
                              }
                  }
            },
            Entry{ScalarEntry{"BMOIP",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  return getValue(ectx.intQuants.oxygenConcentration()) *
                                         getValue(ectx.intQuants.porosity()) *
                                         model.dofTotalVolume(ectx.globalDofIdx);
                              }
                  }
            },
            Entry{ScalarEntry{"BMUIP",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  return getValue(ectx.intQuants.ureaConcentration()) *
                                         getValue(ectx.intQuants.porosity()) *
                                         model.dofTotalVolume(ectx.globalDofIdx) * 1;
                              }
                  }
            },
            Entry{ScalarEntry{"BMBIP",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  return model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.biofilmMass());
                              }
                  }
            },
            Entry{ScalarEntry{"BMCIP",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  return model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.calciteMass());
                              }
                  }
            },
            Entry{ScalarEntry{"BGMIP",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  Scalar result = getValue(ectx.fs.invB(gasPhaseIdx)) *
                                                  getValue(ectx.fs.saturation(gasPhaseIdx));

                                  if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                                       result += getValue(ectx.fs.Rs()) *
                                                 getValue(ectx.fs.invB(oilPhaseIdx)) *
                                                 getValue(ectx.fs.saturation(oilPhaseIdx));
                                  }
                                  else {
                                      result += getValue(ectx.fs.Rsw()) *
                                                getValue(ectx.fs.invB(waterPhaseIdx)) *
                                                getValue(ectx.fs.saturation(waterPhaseIdx));
                                  }
                                  const Scalar rhoG = FluidSystem::referenceDensity(gasPhaseIdx,
                                                                                    ectx.intQuants.pvtRegionIndex());
                                  return result *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity()) *
                                         rhoG;
                              }
                  }
            },
            Entry{ScalarEntry{"BGMGP",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  const Scalar rhoG = FluidSystem::referenceDensity(gasPhaseIdx,
                                                                                    ectx.intQuants.pvtRegionIndex());
                                  return getValue(ectx.fs.invB(gasPhaseIdx)) *
                                         getValue(ectx.fs.saturation(gasPhaseIdx)) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity()) *
                                         rhoG;
                              }
                  }
            },
            Entry{ScalarEntry{"BGMDS",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  Scalar result;
                                  if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                                      result  = getValue(ectx.fs.Rs()) *
                                                getValue(ectx.fs.invB(oilPhaseIdx)) *
                                                getValue(ectx.fs.saturation(oilPhaseIdx));
                                  }
                                  else {
                                      result = getValue(ectx.fs.Rsw()) *
                                               getValue(ectx.fs.invB(waterPhaseIdx)) *
                                               getValue(ectx.fs.saturation(waterPhaseIdx));
                                  }
                                  const Scalar rhoG = FluidSystem::referenceDensity(gasPhaseIdx,
                                                                                    ectx.intQuants.pvtRegionIndex());
                                  return result *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity()) *
                                         rhoG;
                              }
                  }
            },
            Entry{ScalarEntry{"BGMST",
                              [&model = this->simulator_.model(),
                               &problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const auto& scaledDrainageInfo = problem.materialLawManager()
                                                                   ->oilWaterScaledEpsInfoDrainage(ectx.dofIdx);
                                  const Scalar sg = getValue(ectx.fs.saturation(gasPhaseIdx));
                                  Scalar strandedGas = scaledDrainageInfo.Sgcr;
                                  if (problem.materialLawManager()->enableHysteresis()) {
                                      const auto& matParams = problem.materialLawParams(ectx.dofIdx);
                                      const Scalar krg = getValue(ectx.intQuants.relativePermeability(gasPhaseIdx));
                                      strandedGas = MaterialLaw::strandedGasSaturation(matParams, sg, krg);
                                  }
                                  const Scalar xgW = FluidSystem::phaseIsActive(waterPhaseIdx) ?
                                    FluidSystem::convertRvwToXgW(getValue(ectx.fs.Rvw()), ectx.intQuants.pvtRegionIndex())
                                  : FluidSystem::convertRvToXgO(getValue(ectx.fs.Rv()), ectx.intQuants.pvtRegionIndex());
                                  return (1.0 - xgW) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity()) * 
                                         getValue(ectx.fs.density(gasPhaseIdx)) *
                                         std::min(strandedGas, sg);
                              }
                  }
            },
            Entry{ScalarEntry{"BGMUS",
                              [&model = this->simulator_.model(),
                               &problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const auto& scaledDrainageInfo = problem.materialLawManager()
                                                                   ->oilWaterScaledEpsInfoDrainage(ectx.dofIdx);
                                  const Scalar sg = getValue(ectx.fs.saturation(gasPhaseIdx));
                                  Scalar strandedGas = scaledDrainageInfo.Sgcr;
                                  if (problem.materialLawManager()->enableHysteresis()) {
                                      const auto& matParams = problem.materialLawParams(ectx.dofIdx);
                                      const Scalar krg = getValue(ectx.intQuants.relativePermeability(gasPhaseIdx));
                                      strandedGas = MaterialLaw::strandedGasSaturation(matParams, sg, krg);
                                  }
                                  const Scalar xgW = FluidSystem::phaseIsActive(waterPhaseIdx) ?
                                    FluidSystem::convertRvwToXgW(getValue(ectx.fs.Rvw()), ectx.intQuants.pvtRegionIndex())
                                  : FluidSystem::convertRvToXgO(getValue(ectx.fs.Rv()), ectx.intQuants.pvtRegionIndex());
                                  return (1.0 - xgW) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity()) * 
                                         getValue(ectx.fs.density(gasPhaseIdx)) *
                                         std::max(Scalar{0.0}, sg - strandedGas);
                              }
                  }
            },
            Entry{ScalarEntry{"BGMTR",
                              [&model = this->simulator_.model(),
                               &problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const auto& scaledDrainageInfo = problem.materialLawManager()
                                                                   ->oilWaterScaledEpsInfoDrainage(ectx.dofIdx);
                                  Scalar trappedGas = scaledDrainageInfo.Sgcr;
                                  if (problem.materialLawManager()->enableHysteresis()) {
                                      const auto& matParams = problem.materialLawParams(ectx.dofIdx);
                                      trappedGas = MaterialLaw::trappedGasSaturation(matParams, /*maxTrapping*/true);
                                  }
                                  const Scalar xgW = FluidSystem::phaseIsActive(waterPhaseIdx) ?
                                    FluidSystem::convertRvwToXgW(getValue(ectx.fs.Rvw()), ectx.intQuants.pvtRegionIndex())
                                  : FluidSystem::convertRvToXgO(getValue(ectx.fs.Rv()), ectx.intQuants.pvtRegionIndex());
                                  return (1.0 - xgW) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity()) * 
                                         getValue(ectx.fs.density(gasPhaseIdx)) *
                                         std::min(trappedGas, getValue(ectx.fs.saturation(gasPhaseIdx)));
                              }
                  }
            },
            Entry{ScalarEntry{"BGMMO",
                              [&model = this->simulator_.model(),
                               &problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const auto& scaledDrainageInfo = problem.materialLawManager()
                                                                   ->oilWaterScaledEpsInfoDrainage(ectx.dofIdx);
                                  Scalar trappedGas = scaledDrainageInfo.Sgcr;
                                  if (problem.materialLawManager()->enableHysteresis()) {
                                      const auto& matParams = problem.materialLawParams(ectx.dofIdx);
                                      trappedGas = MaterialLaw::trappedGasSaturation(matParams, /*maxTrapping*/true);
                                  }
                                  const Scalar xgW = FluidSystem::phaseIsActive(waterPhaseIdx) ?
                                    FluidSystem::convertRvwToXgW(getValue(ectx.fs.Rvw()), ectx.intQuants.pvtRegionIndex())
                                  : FluidSystem::convertRvToXgO(getValue(ectx.fs.Rv()), ectx.intQuants.pvtRegionIndex());
                                  return (1.0 - xgW) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity()) * 
                                         getValue(ectx.fs.density(gasPhaseIdx)) *
                                         std::max(Scalar{0.0}, getValue(ectx.fs.saturation(gasPhaseIdx)) - trappedGas);
                              }
                  }
            },
            Entry{ScalarEntry{"BGKTR",
                              [&model = this->simulator_.model(),
                               &problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const auto& scaledDrainageInfo = problem.materialLawManager()
                                                                   ->oilWaterScaledEpsInfoDrainage(ectx.dofIdx);
                                  const Scalar sg = getValue(ectx.fs.saturation(gasPhaseIdx));                                
                                  Scalar sgcr = scaledDrainageInfo.Sgcr;
                                  if (problem.materialLawManager()->enableHysteresis()) {
                                      const auto& matParams = problem.materialLawParams(ectx.dofIdx);
                                      sgcr = MaterialLaw::trappedGasSaturation(matParams, /*maxTrapping*/false);
                                  }
                                  if (sg > sgcr) {
                                      return 0.0;
                                  }
                                  else {
                                      const Scalar xgW = FluidSystem::phaseIsActive(waterPhaseIdx) ?
                                        FluidSystem::convertRvwToXgW(getValue(ectx.fs.Rvw()), ectx.intQuants.pvtRegionIndex())
                                      : FluidSystem::convertRvToXgO(getValue(ectx.fs.Rv()), ectx.intQuants.pvtRegionIndex());
                                      return (1.0 - xgW) *
                                             model.dofTotalVolume(ectx.globalDofIdx) *
                                             getValue(ectx.intQuants.porosity()) * 
                                             getValue(ectx.fs.density(gasPhaseIdx)) *
                                             getValue(ectx.fs.saturation(gasPhaseIdx));
                                  }
                              }
                  }
            },
            Entry{ScalarEntry{"BGKMO",
                              [&model = this->simulator_.model(),
                               &problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const auto& scaledDrainageInfo = problem.materialLawManager()
                                                                   ->oilWaterScaledEpsInfoDrainage(ectx.dofIdx);
                                  const Scalar sg = getValue(ectx.fs.saturation(gasPhaseIdx));                                
                                  Scalar sgcr = scaledDrainageInfo.Sgcr;
                                  if (problem.materialLawManager()->enableHysteresis()) {
                                      const auto& matParams = problem.materialLawParams(ectx.dofIdx);
                                      sgcr = MaterialLaw::trappedGasSaturation(matParams, /*maxTrapping*/false);
                                  }
                                  if (sgcr >= sg) {
                                      return 0.0;
                                  }
                                  else {
                                      const Scalar xgW = FluidSystem::phaseIsActive(waterPhaseIdx) ?
                                        FluidSystem::convertRvwToXgW(getValue(ectx.fs.Rvw()), ectx.intQuants.pvtRegionIndex())
                                      : FluidSystem::convertRvToXgO(getValue(ectx.fs.Rv()), ectx.intQuants.pvtRegionIndex());
                                      return (1.0 - xgW) *
                                             model.dofTotalVolume(ectx.globalDofIdx) *
                                             getValue(ectx.intQuants.porosity()) * 
                                             getValue(ectx.fs.density(gasPhaseIdx)) *
                                             getValue(ectx.fs.saturation(gasPhaseIdx));
                                  }
                              }
                  }
            },
            Entry{ScalarEntry{"BGCDI",
                              [&model = this->simulator_.model(),
                               &problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const auto& scaledDrainageInfo = problem.materialLawManager()
                                                                   ->oilWaterScaledEpsInfoDrainage(ectx.dofIdx);
                                  Scalar sgcr = scaledDrainageInfo.Sgcr;
                                  if (problem.materialLawManager()->enableHysteresis()) {
                                      const auto& matParams = problem.materialLawParams(ectx.dofIdx);
                                      sgcr = MaterialLaw::trappedGasSaturation(matParams, /*maxTrapping*/false);
                                  }
                                  const Scalar xgW = FluidSystem::phaseIsActive(waterPhaseIdx) ?
                                    FluidSystem::convertRvwToXgW(getValue(ectx.fs.Rvw()), ectx.intQuants.pvtRegionIndex())
                                  : FluidSystem::convertRvToXgO(getValue(ectx.fs.Rv()), ectx.intQuants.pvtRegionIndex());
                                  return (1.0 - xgW) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity()) * 
                                         getValue(ectx.fs.density(gasPhaseIdx)) *
                                         std::min(sgcr, getValue(ectx.fs.saturation(gasPhaseIdx))) /
                                         FluidSystem::molarMass(gasCompIdx, ectx.intQuants.pvtRegionIndex());
                              }
                  }
            },
            Entry{ScalarEntry{"BGCDM",
                              [&model = this->simulator_.model(),
                               &problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const auto& scaledDrainageInfo = problem.materialLawManager()
                                                                   ->oilWaterScaledEpsInfoDrainage(ectx.dofIdx);
                                  Scalar sgcr = scaledDrainageInfo.Sgcr;
                                  if (problem.materialLawManager()->enableHysteresis()) {
                                      const auto& matParams = problem.materialLawParams(ectx.dofIdx);
                                      sgcr = MaterialLaw::trappedGasSaturation(matParams, /*maxTrapping*/false);
                                  }
                                  const Scalar xgW = FluidSystem::phaseIsActive(waterPhaseIdx) ?
                                    FluidSystem::convertRvwToXgW(getValue(ectx.fs.Rvw()), ectx.intQuants.pvtRegionIndex())
                                  : FluidSystem::convertRvToXgO(getValue(ectx.fs.Rv()), ectx.intQuants.pvtRegionIndex());
                                  return (1.0 - xgW) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity()) * 
                                         getValue(ectx.fs.density(gasPhaseIdx)) *
                                         std::max(Scalar{0.0}, getValue(ectx.fs.saturation(gasPhaseIdx)) - sgcr) /
                                         FluidSystem::molarMass(gasCompIdx, ectx.intQuants.pvtRegionIndex());
                              }
                  }
            },
            Entry{ScalarEntry{"BGKDI",
                              [&model = this->simulator_.model(),
                               &problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const auto& scaledDrainageInfo = problem.materialLawManager()
                                                                   ->oilWaterScaledEpsInfoDrainage(ectx.dofIdx);
                                  const Scalar sg = getValue(ectx.fs.saturation(gasPhaseIdx));                                
                                  Scalar sgcr = scaledDrainageInfo.Sgcr;
                                  if (problem.materialLawManager()->enableHysteresis()) {
                                      const auto& matParams = problem.materialLawParams(ectx.dofIdx);
                                      sgcr = MaterialLaw::trappedGasSaturation(matParams, /*maxTrapping*/false);
                                  }
                                  if (sg > sgcr) {
                                      return 0.0;
                                  }
                                  else {
                                      const Scalar xgW = FluidSystem::phaseIsActive(waterPhaseIdx) ?
                                        FluidSystem::convertRvwToXgW(getValue(ectx.fs.Rvw()), ectx.intQuants.pvtRegionIndex())
                                      : FluidSystem::convertRvToXgO(getValue(ectx.fs.Rv()), ectx.intQuants.pvtRegionIndex());
                                      return (1.0 - xgW) *
                                             model.dofTotalVolume(ectx.globalDofIdx) *
                                             getValue(ectx.intQuants.porosity()) * 
                                             getValue(ectx.fs.density(gasPhaseIdx)) *
                                             getValue(ectx.fs.saturation(gasPhaseIdx)) /
                                             FluidSystem::molarMass(gasCompIdx, ectx.intQuants.pvtRegionIndex());
                                  }
                              }
                  }
            },
            Entry{ScalarEntry{"BGKDM",
                              [&model = this->simulator_.model(),
                               &problem = this->simulator_.problem()](const Context& ectx)
                              {
                                  const auto& scaledDrainageInfo = problem.materialLawManager()
                                                                   ->oilWaterScaledEpsInfoDrainage(ectx.dofIdx);
                                  const Scalar sg = getValue(ectx.fs.saturation(gasPhaseIdx));                                
                                  Scalar sgcr = scaledDrainageInfo.Sgcr;
                                  if (problem.materialLawManager()->enableHysteresis()) {
                                      const auto& matParams = problem.materialLawParams(ectx.dofIdx);
                                      sgcr = MaterialLaw::trappedGasSaturation(matParams, /*maxTrapping*/false);
                                  }
                                  if (sgcr >= sg) {
                                      return 0.0;
                                  }
                                  else {
                                      const Scalar xgW = FluidSystem::phaseIsActive(waterPhaseIdx) ?
                                        FluidSystem::convertRvwToXgW(getValue(ectx.fs.Rvw()), ectx.intQuants.pvtRegionIndex())
                                      : FluidSystem::convertRvToXgO(getValue(ectx.fs.Rv()), ectx.intQuants.pvtRegionIndex());
                                      return (1.0 - xgW) *
                                             model.dofTotalVolume(ectx.globalDofIdx) *
                                             getValue(ectx.intQuants.porosity()) * 
                                             getValue(ectx.fs.density(gasPhaseIdx)) *
                                             getValue(ectx.fs.saturation(gasPhaseIdx)) /
                                             FluidSystem::molarMass(gasCompIdx, ectx.intQuants.pvtRegionIndex());
                                  }
                              }
                  }
            },
            Entry{ScalarEntry{"BWCD",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  Scalar result;
                                  if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                                      result  = getValue(ectx.fs.Rs()) *
                                                getValue(ectx.fs.invB(oilPhaseIdx)) *
                                                getValue(ectx.fs.saturation(oilPhaseIdx));
                                  }
                                  else {
                                      result = getValue(ectx.fs.Rsw()) *
                                               getValue(ectx.fs.invB(waterPhaseIdx)) *
                                               getValue(ectx.fs.saturation(waterPhaseIdx));
                                  }
                                  const Scalar rhoG = FluidSystem::referenceDensity(gasPhaseIdx,
                                                                                    ectx.intQuants.pvtRegionIndex());
                                  return result *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity()) *
                                         rhoG /
                                         FluidSystem::molarMass(gasCompIdx, ectx.intQuants.pvtRegionIndex());
                              }
                  }
            },
            Entry{ScalarEntry{"BWIPG",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  Scalar result = 0.0;
                                  if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
                                      result  = getValue(ectx.fs.Rvw()) *
                                                getValue(ectx.fs.invB(gasPhaseIdx)) *
                                                getValue(ectx.fs.saturation(gasPhaseIdx));
                                  }
                                  return result *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity());
                              }
                  }
            },
            Entry{ScalarEntry{"BWIPL",
                              [&model = this->simulator_.model()](const Context& ectx)
                              {
                                  return getValue(ectx.fs.invB(waterPhaseIdx)) *
                                         getValue(ectx.fs.saturation(waterPhaseIdx)) *
                                         model.dofTotalVolume(ectx.globalDofIdx) *
                                         getValue(ectx.intQuants.porosity());
                              }
                  }
            },
        };

        this->blockExtractors_ = BlockExtractor::setupExecMap(this->blockData_, handlers);

        this->extraBlockData_.clear();
        if (reportStepNum > 0 && !isSubStep) {
            // check we need extra block pressures for RPTSCHED
            const auto& rpt = this->schedule_[reportStepNum - 1].rpt_config.get();
            if (rpt.contains("WELLS") && rpt.at("WELLS") > 1) {
                this->setupExtraBlockData(reportStepNum,
                                          [&c = this->collectOnIORank_](const int idx)
                                          { return c.isCartIdxOnThisRank(idx); });

                const auto extraHandlers = std::array{
                    pressure_handler,
                };

                this->extraBlockExtractors_ = BlockExtractor::setupExecMap(this->extraBlockData_, extraHandlers);
            }
        }
    }

    const Simulator& simulator_;
    const CollectDataOnIORankType& collectOnIORank_;
    std::vector<typename Extractor::Entry> extractors_;
    typename BlockExtractor::ExecMap blockExtractors_;
    typename BlockExtractor::ExecMap extraBlockExtractors_;
};

} // namespace Opm

#endif // OPM_OUTPUT_BLACK_OIL_MODULE_HPP
