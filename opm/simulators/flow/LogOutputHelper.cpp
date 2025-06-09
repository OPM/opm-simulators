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

#include <config.h>

#include <opm/simulators/flow/LogOutputHelper.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/Well/NameOrder.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/output/eclipse/report/WellSpecification.hpp>

#include <opm/simulators/utils/moduleVersion.hpp>
#include <opm/simulators/utils/PressureAverage.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <string_view>
#include <vector>

#include <fmt/format.h>

namespace {

template <typename IJKString>
void logUniqueFailedCells(const std::string& messageTag,
                          std::string_view   prefix,
                          const std::size_t  maxNumCellsFaillog,
                          const std::vector<int>& cells,
                          IJKString&&        ijkString)
{
    if (cells.empty()) {
        return;
    }

    std::vector<int> sorted(cells);
    std::sort(sorted.begin(), sorted.end());
    auto u = std::unique(sorted.begin(), sorted.end());

    const auto numFailed = static_cast<std::size_t>
        (std::distance(sorted.begin(), u));

    std::ostringstream errlog;
    errlog << prefix << " failed for " << numFailed << " cell"
           << ((numFailed != std::size_t{1}) ? "s" : "")
           << " [" << ijkString(cells[0]);

    const auto maxElems = std::min(maxNumCellsFaillog, numFailed);
    for (auto i = 1 + 0*maxElems; i < maxElems; ++i) {
        errlog << ", " << ijkString(cells[i]);
    }

    if (numFailed > maxNumCellsFaillog) {
        errlog << ", ...";
    }

    errlog << ']';

    Opm::OpmLog::warning(messageTag, errlog.str());
}

std::string injectorCModeToString(const Opm::WellInjectorCMode cmode)
{
    auto cmodeStr = std::string{};

    // CMODE_UNDEFINED => ""
    if (cmode != Opm::WellInjectorCMode::CMODE_UNDEFINED) {
        cmodeStr = Opm::WellInjectorCMode2String(cmode);
    }

    return cmodeStr;
}

std::string producerCModeToString(const Opm::WellProducerCMode cmode)
{
    auto cmodeStr = std::string{};

    // NONE && CMODE_UNDEFINED => ""
    if ((cmode != Opm::WellProducerCMode::NONE) &&
        (cmode != Opm::WellProducerCMode::CMODE_UNDEFINED))
    {
        cmodeStr = Opm::WellProducerCMode2String(cmode);
    }

    return cmodeStr;
}

template<class Array>
std::string formatBorder(const Array& widths)
{
    std::string ss;
    ss = " ";
    std::for_each(widths.begin(), widths.end(),
                  [&ss](const auto w)
                  { ss += fmt::format(":{:->{}}", "", w); });
    ss += ':';

    return ss;
}

template<std::size_t size>
std::string formatTextRow(const std::array<int, size>& widths,
                          const std::array<std::string_view, size>& entries)
{
    std::string ss = " ";
    std::for_each(widths.begin(), widths.end(),
                  [&entries, &ss, i = 0](const auto w) mutable
                  { ss += fmt::format(":{:^{}}", entries[i++], w); });
    ss += ":\n";

    return ss;
}

double conn_value_or_zero(const Opm::SummaryState& st,
                          const std::string& wname,
                          const std::string& name,
                          const std::size_t gindex)
{
    return st.has_conn_var(wname, name, gindex)
        ? st.get_conn_var(wname, name, gindex)
        : 0.0;
}

} // Namespace anonymous

namespace Opm {

template<class Scalar>
LogOutputHelper<Scalar>::ConnData::ConnData(const Connection& conn)
    : I(conn.getI() + 1)
    , J(conn.getJ() + 1)
    , K(conn.getK() + 1)
{}

template<class Scalar>
LogOutputHelper<Scalar>::LogOutputHelper(const EclipseState& eclState,
                                         const Schedule& schedule,
                                         const SummaryState& summaryState,
                                         const std::string& moduleVersionName)
    : eclState_(eclState)
    , schedule_(schedule)
    , summaryState_(summaryState)
    , flowVersionName_(moduleVersionName)
{}

template<class Scalar>
void LogOutputHelper<Scalar>::
cumulative(const std::size_t reportStepNum,
           const bool        withConns) const
{
    this->beginCumulativeReport_();

    using Ix = typename WellCumDataType::WCId;

    const auto& st = summaryState_;
    for (const auto& gname : this->schedule_.groupNames()) {
        std::vector<Scalar> values(WellCumDataType::numWCValues);
        std::vector<std::string> names(WellCumDataType::numWCNames);

        names[Ix::WellName] = gname;

        const auto isField = gname == "FIELD";

        values[Ix::OilProd]         = isField ? st.get("FOPT", 0.0) : st.get_group_var(gname, "GOPT", 0.0);
        values[Ix::WaterProd]       = isField ? st.get("FWPT", 0.0) : st.get_group_var(gname, "GWPT", 0.0);
        values[Ix::GasProd]         = isField ? st.get("FGPT", 0.0) : st.get_group_var(gname, "GGPT", 0.0);
        values[Ix::FluidResVolProd] = isField ? st.get("FVPT", 0.0) : st.get_group_var(gname, "GVPT", 0.0);
        values[Ix::OilInj]          = isField ? st.get("FOIT", 0.0) : st.get_group_var(gname, "GOIT", 0.0);
        values[Ix::WaterInj]        = isField ? st.get("FWIT", 0.0) : st.get_group_var(gname, "GWIT", 0.0);
        values[Ix::GasInj]          = isField ? st.get("FGIT", 0.0) : st.get_group_var(gname, "GGIT", 0.0);
        values[Ix::FluidResVolInj]  = isField ? st.get("FVIT", 0.0) : st.get_group_var(gname, "GVIT", 0.0);

        this->outputCumulativeReportRecord_(values, names, {});
    }

    for (const auto& wname : this->schedule_[reportStepNum].well_order()) {
        const auto& well = this->schedule_[reportStepNum].wells.get(wname);

        std::vector<Scalar> values(WellCumDataType::numWCValues);
        std::vector<std::string> names(WellCumDataType::numWCNames);

        names[Ix::WellName] = wname;

        if (well.isInjector()) {
            const auto& controls = well.injectionControls(st);
            const auto ctlMode = controls.cmode;

            names[Ix::WellType] = "INJ";

            if (ctlMode != Well::InjectorCMode::RATE) {
                names[Ix::WellCTRL] = injectorCModeToString(ctlMode);
            }
            else {
                const auto injType = controls.injector_type;
                if      (injType == InjectorType::OIL)   { names[Ix::WellCTRL] = "ORAT"; }
                else if (injType == InjectorType::WATER) { names[Ix::WellCTRL] = "WRAT"; }
                else if (injType == InjectorType::GAS)   { names[Ix::WellCTRL] = "GRAT"; }
            }
        }
        else if (well.isProducer()) {
            const auto cmode = well.productionControls(st).cmode;

            names[Ix::WellType] = "PROD";
            names[Ix::WellCTRL] = producerCModeToString(cmode);
        }

        values[Ix::WellLocationi] = well.getHeadI() + 1;
        values[Ix::WellLocationj] = well.getHeadJ() + 1;

        values[Ix::OilProd]         = st.get_well_var(wname, "WOPT", 0.0);
        values[Ix::WaterProd]       = st.get_well_var(wname, "WWPT", 0.0);
        values[Ix::GasProd]         = st.get_well_var(wname, "WGPT", 0.0);
        values[Ix::FluidResVolProd] = st.get_well_var(wname, "WVPT", 0.0);
        values[Ix::OilInj]          = st.get_well_var(wname, "WOIT", 0.0);
        values[Ix::WaterInj]        = st.get_well_var(wname, "WWIT", 0.0);
        values[Ix::GasInj]          = st.get_well_var(wname, "WGIT", 0.0);
        values[Ix::FluidResVolInj]  = st.get_well_var(wname, "WVIT", 0.0);

        std::vector<ConnData> connData;
        if (withConns) {
            for (const auto& connection : well.getConnections()) {
                ConnData& conn = connData.emplace_back(connection);
                conn.data.resize(WellCumDataType::numWCValues);
                const auto gindex = connection.global_index() + 1;
                conn.data[Ix::OilProd]         = conn_value_or_zero(st, wname, "COPT", gindex);
                conn.data[Ix::WaterProd]       = conn_value_or_zero(st, wname, "CWPT", gindex);
                conn.data[Ix::GasProd]         = conn_value_or_zero(st, wname, "CGPT", gindex);
                conn.data[Ix::FluidResVolProd] = conn_value_or_zero(st, wname, "CVPT", gindex);
                conn.data[Ix::OilInj]          = conn_value_or_zero(st, wname, "COIT", gindex);
                conn.data[Ix::WaterInj]        = conn_value_or_zero(st, wname, "CWIT", gindex);
                conn.data[Ix::GasInj]          = conn_value_or_zero(st, wname, "CGIT", gindex);
                conn.data[Ix::FluidResVolInj]  = conn_value_or_zero(st, wname, "CVIT", gindex);
            }
        }

        this->outputCumulativeReportRecord_(values, names, connData);
    }

    this->endCumulativeReport_();
}

template<class Scalar>
void LogOutputHelper<Scalar>::
error(const std::vector<int>& failedCellsPbub,
      const std::vector<int>& failedCellsPdew) const
{
    auto ijkString = [this](const std::size_t globalIndex)
    {
        const auto ijk = this->eclState_.gridDims().getIJK(globalIndex);

        return fmt::format("({},{},{})", ijk[0] + 1, ijk[1] + 1, ijk[2] + 1);
    };

    constexpr auto maxNumCellsFaillog = static_cast<std::size_t>(20);

    logUniqueFailedCells("Bubble point numerical problem",
                         "Finding the bubble point pressure",
                         maxNumCellsFaillog,
                         failedCellsPbub,
                         ijkString);

    logUniqueFailedCells("Dew point numerical problem",
                         "Finding the dew point pressure",
                         maxNumCellsFaillog,
                         failedCellsPdew,
                         ijkString);
}

template<class Scalar>
void LogOutputHelper<Scalar>::
fip(const Inplace& inplace,
    const Inplace& initialInplace,
    const std::string& name) const
{
    auto iget = [&name](const Inplace& ip,
                        Inplace::Phase phase,
                        std::size_t idx)
    {
        if (name.empty()) {
            return ip.get(phase);
        }

        return ip.get(name, phase, idx);
    };

    for (std::size_t reg = 1; reg <= (name.empty() ? 1 : inplace.max_region(name)); ++reg) {
        std::unordered_map<Inplace::Phase, Scalar> initial_values;
        std::unordered_map<Inplace::Phase, Scalar> current_values;

        for (const auto& phase : Inplace::phases()) {
            initial_values[phase] = iget(initialInplace, phase, reg);
            current_values[phase] = iget(inplace, phase, reg);
        }

        current_values[Inplace::Phase::DynamicPoreVolume] =
            iget(inplace, Inplace::Phase::DynamicPoreVolume, reg);

        this->fipUnitConvert_(initial_values);
        this->fipUnitConvert_(current_values);

        Scalar regHydroCarbonPoreVolumeAveragedPressure =
            detail::pressureAverage(iget(inplace, Inplace::Phase::PressureHydroCarbonPV, reg),
                                    iget(inplace, Inplace::Phase::HydroCarbonPV, reg),
                                    iget(inplace, Inplace::Phase::PressurePV, reg),
                                    iget(inplace, Inplace::Phase::DynamicPoreVolume, reg),
                                    true);
        this->pressureUnitConvert_(regHydroCarbonPoreVolumeAveragedPressure);
        this->outputRegionFluidInPlace_(std::move(initial_values),
                                        std::move(current_values),
                                        regHydroCarbonPoreVolumeAveragedPressure,
                                        name,
                                        name.empty() ? 0 : reg);
    }
}

template<class Scalar>
void LogOutputHelper<Scalar>::
fipResv(const Inplace& inplace, const std::string& name) const
{
    {
        std::unordered_map<Inplace::Phase, Scalar> current_values;

        for (const auto& phase : Inplace::phases()) {
            current_values[phase] = inplace.get(phase);
        }
        Scalar field_dyn_pv = 0.0;

        for (auto nreg = inplace.max_region(name), reg = 0*nreg + 1; reg <= nreg; ++reg) {
            field_dyn_pv = field_dyn_pv + inplace.get(name, Inplace::Phase::DynamicPoreVolume, reg);
        }
        current_values[Inplace::Phase::DynamicPoreVolume] = field_dyn_pv;

        this->fipUnitConvert_(current_values);
        this->outputResvFluidInPlace_(current_values, 0);
    }

    for (auto nreg = inplace.max_region(), reg = 0*nreg + 1; reg <= nreg; ++reg) {
        std::unordered_map<Inplace::Phase, Scalar> current_values;

        for (const auto& phase : Inplace::phases()) {
            if (reg <= inplace.max_region(name)) {
                current_values[phase] = inplace.get(name, phase, reg);
            }
            else {
                current_values[phase] = 0.0;
            }
        }

        if (reg <= inplace.max_region(name)) {
            current_values[Inplace::Phase::DynamicPoreVolume] =
                inplace.get(name, Inplace::Phase::DynamicPoreVolume, reg);
        }
        else {
            current_values[Inplace::Phase::DynamicPoreVolume] = 0.0;
        }

        this->fipUnitConvert_(current_values);
        this->outputResvFluidInPlace_(current_values, reg);
    }

    OpmLog::note(fmt::format(" {:=^91}", ""));
}


template<class Scalar>
void LogOutputHelper<Scalar>::
csv_header(std::ostringstream& ss) const
{
    ss << "RegName,RegNum"
       << ",OilLiquid,OilVapour,Oil,GasFree,GasDissolved,Gas,Water"
       << ",TotalPoreVolume,PoreVolumeOil,PoreVolumeWater,PoreVolumeGas,PoreVolumeHC"
       << "\n";
}

template<class Scalar>
void LogOutputHelper<Scalar>::
fip_csv(std::ostringstream& ss, const Inplace& initial_inplace, const std::string& name) const
{
    for (size_t n = 0; n < initial_inplace.max_region(name); n++) {
        ss << name << "," << n+1
           << fmt::format(",{:4.2f}", initial_inplace.get(name, Inplace::Phase::OilInLiquidPhase, n+1))
           << fmt::format(",{:4.2f}", initial_inplace.get(name, Inplace::Phase::OilInGasPhase, n+1))
           << fmt::format(",{:4.2f}", initial_inplace.get(name, Inplace::Phase::OIL, n+1))
           << fmt::format(",{:4.2f}", initial_inplace.get(name, Inplace::Phase::GasInGasPhase, n+1))
           << fmt::format(",{:4.2f}", initial_inplace.get(name, Inplace::Phase::GasInLiquidPhase, n+1))
           << fmt::format(",{:4.2f}", initial_inplace.get(name, Inplace::Phase::GAS, n+1))
           << fmt::format(",{:4.2f}", initial_inplace.get(name, Inplace::Phase::WATER, n+1))
           << fmt::format(",{:4.2f}", initial_inplace.get(name, Inplace::Phase::DynamicPoreVolume, n+1))
           << fmt::format(",{:4.2f}", initial_inplace.get(name, Inplace::Phase::OilResVolume, n+1))
           << fmt::format(",{:4.2f}", initial_inplace.get(name, Inplace::Phase::WaterResVolume, n+1))
           << fmt::format(",{:4.2f}", initial_inplace.get(name, Inplace::Phase::GasResVolume, n+1));

        auto hvpv = initial_inplace.get(name, Inplace::Phase::OilResVolume, n+1) +
                    initial_inplace.get(name, Inplace::Phase::GasResVolume, n+1);

        ss << fmt::format(",{:4.2f}", hvpv) << "\n";
    }
}

template<class Scalar>
void LogOutputHelper<Scalar>::
timeStamp(const std::string& lbl, double elapsed, int rstep, boost::posix_time::ptime currentDate) const
{
    std::ostringstream ss;
    boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d %b %Y");
    ss.imbue(std::locale(std::locale::classic(), facet));

    ss  << "\n                              **************************************************************************\n"
        << "  " << std::left << std::setw(9) << lbl << "AT" << std::right << std::setw(10)
        << (double)unit::convert::to(elapsed, unit::day) << "  DAYS" << " *" << std::setw(30) << eclState_.getTitle() << "                                          *\n"
        << "  REPORT " << std::setw(4) << rstep << "    " << currentDate
        << "  *                                             Flow  version " << std::setw(11) << flowVersionName_ << "  *\n"
        << "                              **************************************************************************\n";

    OpmLog::note(ss.str());
}

template<class Scalar>
void LogOutputHelper<Scalar>::
injection(const std::size_t reportStepNum,
          const std::map<std::pair<std::string,int>, double>& block_pressures) const
{
    this->beginInjectionReport_();

    using Ix = typename WellInjDataType::WIId;

    const auto& st = this->summaryState_;

    for (const auto& gname : this->schedule_.groupNames()) {
        std::vector<Scalar> values(WellInjDataType::numWIValues);
        std::vector<std::string> names(WellInjDataType::numWINames);

        names[Ix::WellName] = gname;

        const auto isField = gname == "FIELD";

        values[Ix::OilRate]     = isField ? st.get("FOIR", 0.0) : st.get_group_var(gname, "GOIR", 0.0);
        values[Ix::WaterRate]   = isField ? st.get("FWIR", 0.0) : st.get_group_var(gname, "GWIR", 0.0);
        values[Ix::GasRate]     = isField ? st.get("FGIR", 0.0) : st.get_group_var(gname, "GGIR", 0.0);
        values[Ix::FluidResVol] = isField ? st.get("FVIR", 0.0) : st.get_group_var(gname, "GVIR", 0.0);

        this->outputInjectionReportRecord_(values, names, {});
    }

    for (const auto& wname : this->schedule_[reportStepNum].well_order()) {
        const auto& well = this->schedule_[reportStepNum].wells.get(wname);

        // Ignore production wells for the injection report.
        if (well.isProducer()) {
            continue;
        }

        std::vector<Scalar> values(WellInjDataType::numWIValues);
        std::vector<std::string> names(WellInjDataType::numWINames);

        names[Ix::WellName] = wname;

        const auto& controls = well.injectionControls(st);
        const auto ctlMode = controls.cmode;
        const auto injType = controls.injector_type;

        const auto isRate = ctlMode == Well::InjectorCMode::RATE;

        if (injType == InjectorType::OIL) {
            names[Ix::CTRLModeOil] = isRate
                ? "ORAT"
                : injectorCModeToString(ctlMode);
        }
        else if (injType == InjectorType::WATER) {
            names[Ix::CTRLModeWat] = isRate
                ? "WRAT"
                : injectorCModeToString(ctlMode);
        }
        else if (injType == InjectorType::GAS) {
            names[Ix::CTRLModeGas] = isRate
                ? "GRAT"
                : injectorCModeToString(ctlMode);
        }

        values[Ix::WellLocationi] = well.getHeadI() + 1;
        values[Ix::WellLocationj] = well.getHeadJ() + 1;

        values[Ix::OilRate]     = st.get_well_var(wname, "WOIR", 0.0);
        values[Ix::WaterRate]   = st.get_well_var(wname, "WWIR", 0.0);
        values[Ix::GasRate]     = st.get_well_var(wname, "WGIR", 0.0);
        values[Ix::FluidResVol] = st.get_well_var(wname, "WVIR", 0.0);
        values[Ix::BHP]         = st.get_well_var(wname, "WBHP", 0.0);
        values[Ix::THP]         = st.get_well_var(wname, "WTHP", 0.0);

        // values[Ix::SteadyStateII] = 0;

        std::vector<ConnData> connData;
        if (!block_pressures.empty()) {
            const auto& units = this->eclState_.getUnits();
            for (const auto& connection : well.getConnections()) {
                const auto gindex = connection.global_index() + 1;
                const auto bpr_it = block_pressures.find({"BPR", gindex});
                const auto bpr =
                    bpr_it == block_pressures.end()
                         ? 0.0
                         : units.from_si(UnitSystem::measure::pressure, bpr_it->second);
                ConnData& conn = connData.emplace_back(connection);
                conn.data.resize(WellInjDataType::numWIValues);
                conn.data[Ix::OilRate]     = conn_value_or_zero(st, wname, "COIR", gindex);
                conn.data[Ix::WaterRate]   = conn_value_or_zero(st, wname, "CWIR", gindex);
                conn.data[Ix::GasRate]     = conn_value_or_zero(st, wname, "CGIR", gindex);
                conn.data[Ix::FluidResVol] = conn_value_or_zero(st, wname, "CVIR", gindex);
                conn.data[Ix::CPR]         = conn_value_or_zero(st, wname, "CPR", gindex);
                conn.data[Ix::BPR]         = bpr;
            }
        }

        this->outputInjectionReportRecord_(values, names, connData);
    }

    this->endInjectionReport_();
}

template<class Scalar>
void LogOutputHelper<Scalar>::
msw(const std::size_t reportStepNum) const
{
    const auto& sched = this->schedule_[reportStepNum];

    auto msWells = std::vector<std::reference_wrapper<const Well>>{};
    for (const auto& wname : sched.well_order()) {
        if (const auto& well = sched.wells.get(wname); well.isMultiSegment()) {
            msWells.push_back(std::cref(well));
        }
    }

    if (msWells.empty()) {
        return;
    }

    this->beginMSWReport_();

    std::for_each(msWells.begin(), msWells.end(),
                  [this](const Well& well)
                  { this->outputMSWReportRecord_(well); });

    this->endMSWReport_();
}


template<class Scalar>
void LogOutputHelper<Scalar>::
production(const std::size_t reportStepNum,
           const std::map<std::pair<std::string,int>,double>& block_pressures) const
{
    this->beginProductionReport_();

    using Ix = typename WellProdDataType::WPId;

    const auto& st = this->summaryState_;

    for (const auto& gname : this->schedule_.groupNames()) {
        std::vector<Scalar> values(WellProdDataType::numWPValues);
        std::vector<std::string> names(WellProdDataType::numWPNames);

        names[0] = gname;

        const auto isField = gname == "FIELD";

        values[Ix::OilRate]     = isField ? st.get("FOPR", 0.0) : st.get_group_var(gname, "GOPR", 0.0);
        values[Ix::WaterRate]   = isField ? st.get("FWPR", 0.0) : st.get_group_var(gname, "GWPR", 0.0);
        values[Ix::GasRate]     = isField ? st.get("FGPR", 0.0) : st.get_group_var(gname, "GGPR", 0.0);
        values[Ix::FluidResVol] = isField ? st.get("FVPR", 0.0) : st.get_group_var(gname, "GVPR", 0.0);
        values[Ix::WaterCut]    = isField ? st.get("FWCT", 0.0) : st.get_group_var(gname, "GWCT", 0.0);
        values[Ix::GasOilRatio] = isField ? st.get("FGOR", 0.0) : st.get_group_var(gname, "GGOR", 0.0);

        if (values[Ix::WaterRate] == 0.0) {
            values[Ix::WatGasRatio] = 0;
        } else if (values[Ix::GasRate] != 0.0) {
            values[Ix::WatGasRatio] = values[Ix::WaterRate] / values[Ix::GasRate];
        } else {
            values[Ix::WatGasRatio] = 0.0;
        }

        this->outputProductionReportRecord_(values, names, {});
    }

    for (const auto& wname : this->schedule_[reportStepNum].well_order()) {
        const auto& well = this->schedule_[reportStepNum].wells.get(wname);

        // Ignore injection wells for the production report.
        if (well.isInjector()) {
            continue;
        }

        std::vector<Scalar> values(WellProdDataType::numWPValues);
        std::vector<std::string> names(WellProdDataType::numWPNames);

        names[Ix::WellName] = wname;
        names[Ix::CTRLMode] = producerCModeToString(well.productionControls(st).cmode);

        values[Ix::WellLocationi] = well.getHeadI() + 1;
        values[Ix::WellLocationj] = well.getHeadJ() + 1;

        values[Ix::OilRate]     = st.get_well_var(wname, "WOPR", 0.0);
        values[Ix::WaterRate]   = st.get_well_var(wname, "WWPR", 0.0);
        values[Ix::GasRate]     = st.get_well_var(wname, "WGPR", 0.0);
        values[Ix::FluidResVol] = st.get_well_var(wname, "WVPR", 0.0);
        values[Ix::WaterCut]    = st.get_well_var(wname, "WWCT", 0.0);
        values[Ix::GasOilRatio] = st.get_well_var(wname, "WGOR", 0.0);
        values[Ix::BHP]         = st.get_well_var(wname, "WBHP", 0.0);
        values[Ix::THP]         = st.get_well_var(wname, "WTHP", 0.0);

        // values[Ix::SteadyStatePI] = 0;

        auto calc_vat_rate = [](auto& v)
        {
            if (v[Ix::WaterRate] == 0.0) {
                v[Ix::WatGasRatio] = 0.0;
            } else if (v[Ix::GasRate] != 0.0) {
                v[Ix::WatGasRatio] = v[Ix::WaterRate] / v[Ix::GasRate];
            } else {
                v[Ix::WatGasRatio] = 0.0;
            }
        };

        calc_vat_rate(values);

        std::vector<ConnData> connData;
        if (!block_pressures.empty()) {
            const auto& units = this->eclState_.getUnits();
            for (const auto& connection : well.getConnections()) {
                const auto gindex = connection.global_index() + 1;
                const auto bpr_it = block_pressures.find({"BPR", gindex});
                const auto bpr =
                    bpr_it == block_pressures.end()
                         ? 0.0
                         : units.from_si(UnitSystem::measure::pressure, bpr_it->second);
                ConnData& conn = connData.emplace_back(connection);
                conn.data.resize(WellProdDataType::numWPValues);
                conn.data[Ix::OilRate]     = conn_value_or_zero(st, wname, "COPR", gindex);
                conn.data[Ix::WaterRate]   = conn_value_or_zero(st, wname, "CWPR", gindex);
                conn.data[Ix::GasRate]     = conn_value_or_zero(st, wname, "CGPR", gindex);
                conn.data[Ix::FluidResVol] = conn_value_or_zero(st, wname, "CVPR", gindex);
                conn.data[Ix::WaterCut]    = conn_value_or_zero(st, wname, "CWCT", gindex);
                conn.data[Ix::GasOilRatio] = conn_value_or_zero(st, wname, "CGOR", gindex);
                conn.data[Ix::CPR]         = conn_value_or_zero(st, wname, "CPR", gindex);
                conn.data[Ix::BPR]         = bpr;

                calc_vat_rate(conn.data);
            }
        }

        this->outputProductionReportRecord_(values, names, connData);
    }

    this->endProductionReport_();
}

template<class Scalar>
void LogOutputHelper<Scalar>::
wellSpecification(const std::vector<std::string>& changedWells,
                  const std::size_t               reportStepNum) const
{
    std::ostringstream ss;

    auto blockDepth = PrtFile::Reports::BlockDepthCallback {
        [&grid = this->eclState_.getInputGrid()]
        (const std::size_t globalCellIndex)
        { return grid.getCellDepth(globalCellIndex); }
    };

    PrtFile::Reports::wellSpecification(changedWells, reportStepNum,
                                        this->schedule_, blockDepth, ss);

    OpmLog::note(ss.str());
}

template <typename Scalar>
void LogOutputHelper<Scalar>::beginCumulativeReport_() const
{
    const auto& units = this->eclState_.getUnits();
    const auto widths = std::array{8, 11, 8, 4, 11, 11, 11, 11, 11, 11, 11, 11};

    using namespace std::string_view_literals;

    auto oil_unit =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "MSCM"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "MSTB"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "MSCC"sv;
            default: return "unsupp"sv;
            }
        };
    auto wat_unit =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "MSCM"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "MSTB"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "MSCC"sv;
            default: return "unsupp"sv;
            }
        };
    auto gas_unit =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "MMSCM"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "MMSCF"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "MMSCC"sv;
            default: return "unsupp"sv;
            }
        };
    auto res_unit =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "MRCM"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "MRB"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "MRCC"sv;
            default: return "unsupp"sv;
            }
        };

    const auto unitStrings = std::array{
                                ""sv,
                                ""sv,
                                ""sv,
                                ""sv,
                                oil_unit(),
                                wat_unit(),
                                gas_unit(),
                                res_unit(),
                                oil_unit(),
                                wat_unit(),
                                gas_unit(),
                                res_unit(),
                             };

    std::ostringstream ss;
    ss << fmt::format("\n{:=^132}\n", " CUMULATIVE PRODUCTION/INJECTION TOTALS ")
       << formatTextRow(widths,
                        std::array{
                            "WELL"sv,
                            "LOCATION"sv,
                            "WELL"sv,
                            "CTRL"sv,
                            "OIL"sv,
                            "WATER"sv,
                            "GAS"sv,
                            "Prod"sv,
                            "OIL"sv,
                            "WATER"sv,
                            "GAS"sv,
                            "INJ"sv,
                        })
       << formatTextRow(widths,
                        std::array{
                            "NAME"sv,
                            "(I,J,K)"sv,
                            "TYPE"sv,
                            "MODE"sv,
                            "PROD"sv,
                            "PROD"sv,
                            "PROD"sv,
                            "RES.VOL."sv,
                            "INJ"sv,
                            "INJ"sv,
                            "INJ"sv,
                            "RES.VOL."sv,
                        })
       << formatTextRow(widths, unitStrings)
       << fmt::format("{:=>132}", "");

    OpmLog::note(ss.str());
}

template <typename Scalar>
void LogOutputHelper<Scalar>::endCumulativeReport_() const
{
    const auto widths = std::array{8, 11, 8, 4, 11, 11, 11, 11, 11, 11, 11, 11};
    OpmLog::note(formatBorder(widths));
}

template<class Scalar>
void LogOutputHelper<Scalar>::
outputCumulativeReportRecord_(const std::vector<Scalar>& wellCum,
                              const std::vector<std::string>& wellCumNames,
                              const std::vector<ConnData>& connData) const
{
    std::ostringstream ss;

    ss << fmt::format(":{:<8}:", wellCumNames[WellCumDataType::WellName]);

    if (wellCum[WellCumDataType::WellLocationi] < 1) {
        ss << fmt::format("{:11}:", "");
    } else {
        ss << fmt::format("{:>3},{:>3} {:>3}:",
                          wellCum[WellCumDataType::WellLocationi],
                          wellCum[WellCumDataType::WellLocationj],
                          "");
    }

    auto scaledValue = [](const auto& wc, const auto quantity)
    {
        // Unit M*
        return wc[quantity] / 1000.0;
    };

    auto scaledGasValue = [](const auto& wc, const auto quantity)
    {
        // Unit MM*
        return wc[quantity] / (1000.0 * 1000.0);
    };

    ss << fmt::format("{:>8}:{:>4}:{:>11.1f}:{:>11.1f}:{:>11.1f}:{:>11.1f}:"
                      "{:>11.1f}:{:>11.1f}:{:>11.1f}:{:>11.1f}:",
                      wellCumNames[WellCumDataType::WCId::WellType],
                      wellCumNames[WellCumDataType::WCId::WellCTRL],
                      scaledValue(wellCum, WellCumDataType::WCId::OilProd),
                      scaledValue(wellCum, WellCumDataType::WCId::WaterProd),
                      scaledGasValue(wellCum, WellCumDataType::WCId::GasProd),
                      scaledValue(wellCum, WellCumDataType::WCId::FluidResVolProd),
                      scaledValue(wellCum, WellCumDataType::WCId::OilInj),
                      scaledValue(wellCum, WellCumDataType::WCId::WaterInj),
                      scaledGasValue(wellCum, WellCumDataType::WCId::GasInj),
                      scaledValue(wellCum, WellCumDataType::WCId::FluidResVolInj));

    for (const auto& conn : connData) {
        ss << fmt::format("\n:  BLOCK :{0:>3},{1:>3},{2:>3}:{3:>8}:{3:>4}:",
                          conn.I, conn.J, conn.K, "")
           << fmt::format("{:>11.1f}:{:>11.1f}:{:>11.1f}:{:>11.1f}:"
                          "{:>11.1f}:{:>11.1f}:{:>11.1f}:{:>11.1f}:",
                          scaledValue(conn.data, WellCumDataType::WCId::OilProd),
                          scaledValue(conn.data, WellCumDataType::WCId::WaterProd),
                          scaledGasValue(conn.data, WellCumDataType::WCId::GasProd),
                          scaledValue(conn.data, WellCumDataType::WCId::FluidResVolProd),
                          scaledValue(conn.data, WellCumDataType::WCId::OilInj),
                          scaledValue(conn.data, WellCumDataType::WCId::WaterInj),
                          scaledGasValue(conn.data, WellCumDataType::WCId::GasInj),
                          scaledValue(conn.data, WellCumDataType::WCId::FluidResVolInj));
    }

    OpmLog::note(ss.str());
}

template <typename Scalar>
void LogOutputHelper<Scalar>::beginInjectionReport_() const
{
    const auto& units = this->eclState_.getUnits();
    const auto widths = std::array{8, 11, 6, 6, 6, 11, 11, 11, 11, 8, 8};

    using namespace std::string_view_literals;

    auto oil_rate =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "SCM/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "STB/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "SCC/HR"sv;
            default: return "unsupp"sv;
            }
        };
    auto wat_rate =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "SCM/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "STB/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "SCC/HR"sv;
            default: return "unsupp"sv;
            }
        };
    auto gas_rate =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "SCM/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "MSCF/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "SCC/HR"sv;
            default: return "unsupp"sv;
            }
        };
    auto res_vol =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "RCM/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "RB/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "RCC/HR"sv;
            default: return "unsupp"sv;
            }
        };

    const auto unitStrings = std::array{
                                ""sv,
                                ""sv,
                                "OIL"sv,
                                "WAT"sv,
                                "GAS"sv,
                                oil_rate(),
                                wat_rate(),
                                gas_rate(),
                                res_vol(),
                                std::string_view{units.name(UnitSystem::measure::pressure)},
                                std::string_view{units.name(UnitSystem::measure::pressure)},
                             };

    std::ostringstream ss;
    ss << fmt::format("\n{:=^109}\n", " INJECTION REPORT ")
       << formatTextRow(widths,
                        std::array{
                            "WELL"sv,
                            "LOCATION"sv,
                            "CTRL"sv,
                            "CTRL"sv,
                            "CTRL"sv,
                            "OIL"sv,
                            "WATER"sv,
                            "GAS"sv,
                            "FLUID"sv,
                            "BHP OR"sv,
                            "THP OR"sv,
                        })
       << formatTextRow(widths,
                        std::array{
                            "NAME"sv,
                            "(I,J,K)"sv,
                            "MODE"sv,
                            "MODE"sv,
                            "MODE"sv,
                            "RATE"sv,
                            "RATE"sv,
                            "RATE"sv,
                            "RES.VOL."sv,
                            "CON.PR."sv,
                            "BLK.PR."sv,
                        })
       << formatTextRow(widths, unitStrings)
       << fmt::format("{:=>109}", "");

    OpmLog::note(ss.str());
}

template <typename Scalar>
void LogOutputHelper<Scalar>::endInjectionReport_() const
{
    const auto widths = std::array{8, 11, 6, 6, 6, 11, 11, 11, 11, 8, 8};
    OpmLog::note(formatBorder(widths));
}

template<class Scalar>
void LogOutputHelper<Scalar>::
outputInjectionReportRecord_(const std::vector<Scalar>& wellInj,
                             const std::vector<std::string>& wellInjNames,
                             const std::vector<ConnData>& connData) const
{
    const auto isWellRecord =
        wellInj[WellProdDataType::WellLocationi] >= 1;

    std::ostringstream ss;
    ss << fmt::format(":{:<8}:", wellInjNames[WellInjDataType::WellName]);

    if (! isWellRecord) {
        ss << fmt::format("{:11}:", "");
    } else {
        ss << fmt::format("{:>3},{:>3} {:>3}:",
                          wellInj[WellInjDataType::WellLocationi],
                          wellInj[WellInjDataType::WellLocationj],
                          "");
    }

    ss << fmt::format("{:>6}:{:>6}:{:>6}:{:>11.1f}:{:>11.1f}:{:>11.1f}:{:>11.1f}:",
                      wellInjNames[WellInjDataType::CTRLModeOil],
                      wellInjNames[WellInjDataType::CTRLModeWat],
                      wellInjNames[WellInjDataType::CTRLModeGas],
                      wellInj[WellInjDataType::OilRate],
                      wellInj[WellInjDataType::WaterRate],
                      wellInj[WellInjDataType::GasRate],
                      wellInj[WellInjDataType::FluidResVol]);

    if (! isWellRecord) {
        ss << fmt::format("{0:8}:{0:8}:", "");
    } else {
        ss << fmt::format("{:>8.1f}:{:>8.1f}:",
                          wellInj[WellInjDataType::BHP],
                          wellInj[WellInjDataType::THP]);
    }

    for (const auto& conn : connData) {
        ss << fmt::format("\n:  BLOCK :{0:>3},{1:>3},{2:>3}:"
                          "{3:>6}:{3:>6}:{3:>6}:",
                          conn.I, conn.J, conn.K, "")
           << fmt::format("{:>11.1f}:{:>11.1f}:{:>11.1f}:{:>11.1f}:{:>8.1f}:{:>8.1f}:",
                          conn.data[WellInjDataType::OilRate],
                          conn.data[WellInjDataType::WaterRate],
                          conn.data[WellInjDataType::GasRate],
                          conn.data[WellInjDataType::FluidResVol],
                          conn.data[WellInjDataType::CPR],
                          conn.data[WellInjDataType::BPR]);
    }

    OpmLog::note(ss.str());
}

template <typename Scalar>
void LogOutputHelper<Scalar>::beginMSWReport_() const
{
    const auto& units = this->eclState_.getUnits();
    const auto widths  = std::array{10, 5, 5, 11, 11, 11, 9, 16, 9, 8, 8, 8};
    const auto widths2 = std::array{10, 5, 5, 11, 11, 11, 9, 16, 9, 26};

    using namespace std::string_view_literals;

    std::ostringstream ss;
    ss << fmt::format("\n{:=^124}\n", " MULTI-SEGMENT WELL REPORT ")
       << formatTextRow(widths2,
                        std::array{
                            "WELL"sv,
                            "BRN"sv,
                            "SEG"sv,
                            "OIL"sv,
                            "WATER"sv,
                            "GAS"sv,
                            "MIXTURE"sv,
                            "HOLDUP FRACTION"sv,
                            "PRESSURE"sv,
                            "PRESSURE HEAD LOSSES"sv,
                        })
       << formatTextRow(widths,
                        std::array{
                            "NAME"sv,
                            "NO."sv,
                            "NO."sv,
                            "FLOW"sv,
                            "FLOW"sv,
                            "FLOW"sv,
                            "VELOCITY"sv,
                            "OIL  WAT  GAS"sv,
                            ""sv,
                            "H-STATIC"sv,
                            "FRICTION"sv,
                            "ACCELRTN"sv,
                        })
       << formatTextRow(widths,
                        std::array{
                            ""sv,
                            ""sv,
                            ""sv,
                            std::string_view{units.name(UnitSystem::measure::liquid_surface_rate)},
                            std::string_view{units.name(UnitSystem::measure::liquid_surface_rate)},
                            std::string_view{units.name(UnitSystem::measure::gas_surface_rate)},
                            std::string_view{units.name(UnitSystem::measure::pipeflow_velocity)},
                            ""sv,
                            std::string_view{units.name(UnitSystem::measure::pressure)},
                            std::string_view{units.name(UnitSystem::measure::pressure_drop)},
                            std::string_view{units.name(UnitSystem::measure::pressure_drop)},
                            std::string_view{units.name(UnitSystem::measure::pressure_drop)},
                        })
       << fmt::format("{:=>124}", "");

    OpmLog::note(ss.str());
}

template <typename Scalar>
void LogOutputHelper<Scalar>::endMSWReport_() const
{
    const auto widths = std::array{10, 5, 5, 11, 11, 11, 9, 16, 9, 8, 8, 8};
    OpmLog::note(formatBorder(widths));
}

template<class Scalar>
void LogOutputHelper<Scalar>::
outputMSWReportRecord_(const Well& well) const
{
    const auto& st = this->summaryState_;

    using namespace std::string_view_literals;

    auto clampToZero = [](const auto& number)
    {
        return std::fabs(number) > 1e-6 ? number : 0.0;
    };

    auto get_phase_values =
        [&st,
         &clampToZero,
         &wname = well.name()](const int segmentNumber,
                               const std::array<std::string_view, 3>& fields)
        {
            std::array<Scalar, 3> result;
            for (std::size_t p = 0; p < 3; ++p) {
                const auto value =
                    st.has_segment_var(wname, std::string{fields[p]}, segmentNumber)
                        ? st.get_segment_var(wname, std::string{fields[p]}, segmentNumber)
                        : 0.0;
                result[p] = clampToZero(value);
            }
            return result;
        };

    std::ostringstream ss;
    ss << fmt::format(": {:<9}:",  well.name());

    for (int i = 1; i <= well.maxBranchID(); ++i) {
        if (i != 1) {
            ss << fmt::format(": {:<9}:", "");
        }
        ss << fmt::format("{:^5}:", i);
        const auto& segments = well.getSegments().branchSegments(i);
        bool first = true;

        for (const auto& segment : segments) {
            if (!first) {
                ss << fmt::format("\n:{0:>10}:{0:>5}:", "");
            }
            const auto rates = get_phase_values(segment.segmentNumber(),
                                                std::array{"SOFR"sv, "SWFR"sv, "SGFR"sv});
            const auto holdups = get_phase_values(segment.segmentNumber(),
                                                  std::array{"SOHF"sv, "SWHF"sv, "SGHF"sv});
            const auto velocities = get_phase_values(segment.segmentNumber(),
                                                  std::array{"SOFV"sv, "SWFV"sv, "SGFV"sv});
            const auto press_drop = get_phase_values(segment.segmentNumber(),
                                                     std::array{"SPRDH"sv, "SPRDF"sv, "SPRDA"sv});

            const auto mixture_vel = std::inner_product(holdups.begin(), holdups.end(),
                                                        velocities.begin(), 0.0);

            ss << fmt::format("{:^5}:{:>11.1f}:{:>11.1f}:{:>11.1f}:{:>9.3f}:"
                              " {:>3.2f} {:>3.2f} {:>3.2f} :"
                              "{:>9.1f}:{:>8.3f}:{:>8.3f}:{:>8.3f}:",
                              segment.segmentNumber(),
                              rates[0],
                              rates[1],
                              rates[2],
                              mixture_vel,
                              holdups[0],
                              holdups[1],
                              holdups[2],
                              clampToZero(st.has_segment_var(well.name(), "SPR",
                                                             segment.segmentNumber())
                                   ? st.get_segment_var(well.name(), "SPR",
                                                        segment.segmentNumber())
                                   : 0.0
                              ),
                              press_drop[0],
                              press_drop[1],
                              press_drop[2]);
            first = false;
        }
    }
    OpmLog::note(ss.str());
}

template <typename Scalar>
void LogOutputHelper<Scalar>::beginProductionReport_() const
{
    const auto& units = this->eclState_.getUnits();
    const auto widths = std::array{8, 11, 4, 11, 11, 11, 11, 11, 10, 12, 8, 8};

    using namespace std::string_view_literals;

    auto oil_rate =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "SCM/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "STB/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "SCC/HR"sv;
            default: return "unsupp"sv;
            }
        };
    auto wat_rate =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "SCM/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "STB/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "SCC/HR"sv;
            default: return "unsupp"sv;
            }
        };
    auto gas_rate =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "RCM/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "MSCF/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "RCC"sv;
            default: return "unsupp"sv;
            }
        };
    auto res_vol =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "SCM/SCM"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "RB/DAY"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "SCC/SCC"sv;
            default: return "unsupp"sv;
            }
        };
    auto water_cut =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "SCM/SCM"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return ""sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "SCC/SCC"sv;
            default: return "unsupp"sv;
            }
        };
    auto gas_oil_rate =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "SCM/SCM"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "MSCF/STB"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "SCC/SCC"sv;
            default: return "unsupp"sv;
            }
        };
    auto wat_gas_rate =
        [utype = units.getType()]()
        {
            switch (utype) {
            case UnitSystem::UnitType::UNIT_TYPE_METRIC: return "SCM/SCM"sv;
            case UnitSystem::UnitType::UNIT_TYPE_FIELD: return "STB/MSCF"sv;
            case UnitSystem::UnitType::UNIT_TYPE_LAB: return "SCC/SCC"sv;
            default: return "unsupp"sv;
            }
        };

    const auto unitStrings = std::array{
                                ""sv,
                                ""sv,
                                ""sv,
                                oil_rate(),
                                wat_rate(),
                                gas_rate(),
                                res_vol(),
                                water_cut(),
                                gas_oil_rate(),
                                wat_gas_rate(),
                                std::string_view{units.name(UnitSystem::measure::pressure)},
                                std::string_view{units.name(UnitSystem::measure::pressure)},
                             };

    std::ostringstream ss;
    ss << fmt::format("\n{:=^129}\n", " PRODUCTION REPORT ")
       << formatTextRow(widths,
                        std::array{
                            "WELL"sv,
                            "LOCATION"sv,
                            "CTRL"sv,
                            "OIL"sv,
                            "WATER"sv,
                            "GAS"sv,
                            "FLUID"sv,
                            "WATER"sv,
                            "GAS/OIL"sv,
                            "WAT/GAS"sv,
                            "BHP OR"sv,
                            "THP OR"sv,
                        })
       << formatTextRow(widths,
                        std::array{
                            "NAME"sv,
                            "(I,J,K)"sv,
                            "MODE"sv,
                            "RATE"sv,
                            "RATE"sv,
                            "RATE"sv,
                            "RES.VOL."sv,
                            "CUT"sv,
                            "RATIO"sv,
                            "RATIO"sv,
                            "CON.PR."sv,
                            "BLK.PR."sv,
                        })
       << formatTextRow(widths, unitStrings)
       << fmt::format("{:=>129}", "");

    OpmLog::note(ss.str());
}

template <typename Scalar>
void LogOutputHelper<Scalar>::endProductionReport_() const
{
    const auto widths = std::array{8, 11, 4, 11, 11, 11, 11, 11, 10, 12, 8, 8};
    OpmLog::note(formatBorder(widths));
}

template<class Scalar>
void LogOutputHelper<Scalar>::
outputProductionReportRecord_(const std::vector<Scalar>& wellProd,
                              const std::vector<std::string>& wellProdNames,
                              const std::vector<ConnData>& connData) const
{
    const auto isWellRecord =
        wellProd[WellProdDataType::WellLocationi] >= 1;

    std::ostringstream ss;

    ss << fmt::format(":{:<8}:", wellProdNames[WellProdDataType::WellName]);

    if (! isWellRecord) {
        ss << fmt::format("{:11}:", "");
    } else {
        ss << fmt::format("{:>3},{:>3} {:>3}:",
                          wellProd[WellProdDataType::WellLocationi] ,
                          wellProd[WellProdDataType::WellLocationj],
                          "");
    }

    ss << fmt::format("{:>4}:{:>11.1f}:{:>11.1f}:{:>11.1f}:{:>11.1f}:{:>11.3f}:{:>10.2f}:{:>12.4f}:",
                      wellProdNames[WellProdDataType::CTRLMode],
                      wellProd[WellProdDataType::OilRate],
                      wellProd[WellProdDataType::WaterRate],
                      wellProd[WellProdDataType::GasRate],
                      wellProd[WellProdDataType::FluidResVol],
                      wellProd[WellProdDataType::WaterCut],
                      wellProd[WellProdDataType::GasOilRatio],
                      wellProd[WellProdDataType::WatGasRatio]);

    if (! isWellRecord) {
        ss << fmt::format("{0:8}:{0:8}:", "");
    } else {
        ss << fmt::format("{:>8.1f}:{:>8.1f}:",
                          wellProd[WellProdDataType::BHP],
                          wellProd[WellProdDataType::THP]);
    }

    for (const auto& conn : connData) {
        ss << fmt::format("\n:  BLOCK :{0:>3},{1:>3},{2:>3}:{3:>4}:",
                          conn.I, conn.J, conn.K, "")
           << fmt::format("{:>11.1f}:{:>11.1f}:{:>11.1f}:{:>11.1f}:"
                          "{:>11.3f}:{:>10.2f}:{:>12.4f}:{:>8.1f}:{:>8.1f}:",
                          conn.data[WellProdDataType::OilRate],
                          conn.data[WellProdDataType::WaterRate],
                          conn.data[WellProdDataType::GasRate],
                          conn.data[WellProdDataType::FluidResVol],
                          conn.data[WellProdDataType::WaterCut],
                          conn.data[WellProdDataType::GasOilRatio],
                          conn.data[WellProdDataType::WatGasRatio],
                          conn.data[WellProdDataType::CPR],
                          conn.data[WellProdDataType::BPR]);
    }

    OpmLog::note(ss.str());
}

template<class Scalar>
void LogOutputHelper<Scalar>::
outputRegionFluidInPlace_(std::unordered_map<Inplace::Phase, Scalar> oip,
                          std::unordered_map<Inplace::Phase, Scalar> cip,
                          const Scalar pav,
                          const std::string& name,
                          const int reg) const
{
    // don't output FIPNUM report if the region has no porv.
    if (! (cip[Inplace::Phase::PoreVolume] > Scalar{0})) {
        return;
    }

    const UnitSystem& units = eclState_.getUnits();
    const auto widths = std::array{25, 43, 16, 43};

    auto topEntry = [&units](const char* text,
                             const Scalar value,
                             const UnitSystem::measure m)
    {
        return fmt::format(" {0: >52}:{0: >8}{1:<5}={2:>14.0f}  {3:<18}:\n",
                           "", text, value, units.name(m));
    };

    auto topLabel = [](const char* text)
    {
        return fmt::format(" {: >52}: {:^47}:\n", "", text);
    };

    auto formatEntry = [](auto& fip, std::string_view label)
    {
        return fmt::format(" :{:<25}:{:^15.0f}{:^14.0f}{:^14.0f}:{:^16.0f}:"
                           "{:^15.0f}{:^14.0f}{:^14.0f}:\n",
                           label,
                           fip[Inplace::Phase::OilInLiquidPhase],
                           fip[Inplace::Phase::OilInGasPhase],
                           fip[Inplace::Phase::OIL],
                           fip[Inplace::Phase::WATER],
                           fip[Inplace::Phase::GasInGasPhase],
                           fip[Inplace::Phase::GasInLiquidPhase],
                           fip[Inplace::Phase::GAS]);
    };

    std::ostringstream ss;
    ss << fmt::format("\n {0: >52}{0:=>50}\n", "");
    if (reg == 0) {
        ss << fmt::format(" {0: >52}:{1:^48}:\n", "", "FIELD TOTALS");
    }
    else {
        ss << fmt::format(" {0: >52}:{1:>20} REPORT REGION {2:>{3}}{0: ^11}:\n",
                          "", name, reg, 8 - name.size());
    }
    ss << topEntry("PAV", pav, UnitSystem::measure::pressure)
       << topEntry("PORV", cip[Inplace::Phase::PoreVolume], UnitSystem::measure::volume);
    if (!reg) {
        ss << topLabel("Pressure is weighted by hydrocarbon pore volume")
           << topLabel("Pore volumes are taken at reference conditions");
    }

    ss << fmt::format(" {0: >26}:{0:->15} OIL {1:>4} {0:->18}:"
                      "-- WAT{0: >2} {1:>4} --:{0:->14}  GAS{0: >3} {2:>4} {0:-^15}:\n",
                      "",
                      units.name(UnitSystem::measure::liquid_surface_volume),
                      units.name(UnitSystem::measure::gas_surface_volume))
       << fmt::format(" {0: >26}:{1:^15}{2:^14}{3:^14}:{3:^16}:{4:^15}{5:^14}{3:^14}:\n",
                      "",
                      "LIQUID",
                      "VAPOUR",
                      "TOTAL",
                      "FREE",
                      "DISSOLVED")
       << formatBorder(widths) << '\n'
       << formatEntry(cip, "CURRENTLY IN PLACE")
       << formatBorder(widths) << '\n'
       << formatEntry(oip, "ORIGINALLY IN PLACE");

    if (reg != 0) {
        ss << formatBorder(widths) << '\n';
    }
    ss << fmt::format(" {:=^132}\n\n", "");

    OpmLog::note(ss.str());
}

template<class Scalar>
void LogOutputHelper<Scalar>::
outputResvFluidInPlace_(std::unordered_map<Inplace::Phase, Scalar> cipr,
                        const int reg) const
{
    const UnitSystem& units = eclState_.getUnits();
    using namespace std::string_view_literals;

    std::ostringstream ss;
    if (reg == 0) {
        const auto widths = std::array{9, 15, 15, 15, 15, 15};
        ss << fmt::format("\n {0: >52}{0:=>35}", "")
           << fmt::format("\n {0: >52}:  RESERVOIR VOLUMES {1:^13}:\n",
                          "",  units.name(UnitSystem::measure::volume))
           << formatBorder(widths) << '\n'
           << formatTextRow(widths,
                            std::array{
                                "REGION"sv,
                                "TOTAL PORE"sv,
                                "PORE VOLUME"sv,
                                "PORE VOLUME"sv,
                                "PORE VOLUME"sv,
                                "PORE VOLUME"sv,
                            })
           << formatTextRow(widths,
                            std::array{
                                ""sv,
                                "VOLUME"sv,
                                "CONTAINING"sv,
                                "CONTAINING"sv,
                                "CONTAINING"sv,
                                "CONTAINING"sv,
                            })
           << formatTextRow(widths,
                            std::array{
                                ""sv,
                                ""sv,
                                "OIL"sv,
                                "WATER"sv,
                                "GAS"sv,
                                "HYDRO-CARBON"sv,
                            })
           << formatBorder(widths) << '\n'
           << fmt::format(" :{:<9}:", "FIELD");
    } else {
        ss << fmt::format(" :{:<9}:", reg);
    }

    ss << fmt::format("{0:>15.0f}:{1:>15.0f}:{2:>15.0f}:{3:>15.0f}:{4:>15.0f}:",
                      cipr[Inplace::Phase::DynamicPoreVolume],
                      cipr[Inplace::Phase::OilResVolume],
                      cipr[Inplace::Phase::WaterResVolume],
                      cipr[Inplace::Phase::GasResVolume],
                      cipr[Inplace::Phase::OilResVolume] + cipr[Inplace::Phase::GasResVolume]);

    OpmLog::note(ss.str());
}

template<class Scalar>
void LogOutputHelper<Scalar>::
fipUnitConvert_(std::unordered_map<Inplace::Phase, Scalar>& fip) const
{
    const UnitSystem& units = eclState_.getUnits();
    using M = UnitSystem::measure;
    const auto unit_map = std::unordered_map<Inplace::Phase, M> {
        {Inplace::Phase::WATER,                     M::liquid_surface_volume},
        {Inplace::Phase::OIL,                       M::liquid_surface_volume},
        {Inplace::Phase::OilInLiquidPhase,          M::liquid_surface_volume},
        {Inplace::Phase::OilInGasPhase,             M::liquid_surface_volume},
        {Inplace::Phase::GAS,                       M::gas_surface_volume},
        {Inplace::Phase::GasInLiquidPhase,          M::gas_surface_volume},
        {Inplace::Phase::GasInGasPhase,             M::gas_surface_volume},
        {Inplace::Phase::PoreVolume,                M::volume},
        {Inplace::Phase::DynamicPoreVolume,         M::volume},
        {Inplace::Phase::WaterResVolume,            M::volume},
        {Inplace::Phase::OilResVolume,              M::volume},
        {Inplace::Phase::GasResVolume,              M::volume},
        {Inplace::Phase::SALT,                      M::mass},
        {Inplace::Phase::CO2InWaterPhase,           M::moles},
        {Inplace::Phase::CO2InGasPhaseInMob,        M::moles},
        {Inplace::Phase::CO2InGasPhaseMob,          M::moles},
        {Inplace::Phase::CO2InGasPhaseInMobKrg,     M::moles},
        {Inplace::Phase::CO2InGasPhaseMobKrg,       M::moles},
        {Inplace::Phase::WaterInWaterPhase,         M::liquid_surface_volume},
        {Inplace::Phase::WaterInGasPhase,           M::liquid_surface_volume},
        {Inplace::Phase::WaterMass,                 M::mass},
        {Inplace::Phase::CO2Mass,                   M::mass},
        {Inplace::Phase::CO2MassInWaterPhase,       M::mass},
        {Inplace::Phase::CO2MassInGasPhase,         M::mass},
        {Inplace::Phase::CO2MassInGasPhaseInMob,    M::mass},
        {Inplace::Phase::CO2MassInGasPhaseMob,      M::mass},
        {Inplace::Phase::CO2MassInGasPhaseInMobKrg, M::mass},
        {Inplace::Phase::CO2MassInGasPhaseMobKrg,   M::mass},
        {Inplace::Phase::CO2MassInGasPhaseEffectiveTrapped, M::mass},
        {Inplace::Phase::CO2MassInGasPhaseEffectiveUnTrapped, M::mass},
        {Inplace::Phase::CO2MassInGasPhaseMaximumTrapped,  M::mass},
        {Inplace::Phase::CO2MassInGasPhaseMaximumUnTrapped,  M::mass},

    };

    for (auto& [phase, value] : fip) {
        auto unitPos = unit_map.find(phase);
        if (unitPos != unit_map.end()) {
            value = units.from_si(unitPos->second, value);
        }
    }
}

template<class Scalar>
void LogOutputHelper<Scalar>::
pressureUnitConvert_(Scalar& pav) const
{
    pav = eclState_.getUnits()
        .from_si(UnitSystem::measure::pressure, pav);
}

template class LogOutputHelper<double>;

#if FLOW_INSTANTIATE_FLOAT
template class LogOutputHelper<float>;
#endif


} // namespace Opm
