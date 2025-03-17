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
#include <opm/simulators/utils/moduleVersion.hpp>
#include <opm/simulators/flow/LogOutputHelper.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEnums.hpp>

#include <opm/simulators/utils/PressureAverage.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
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
    std::string ss;
    std::for_each(widths.begin(), widths.end(),
                  [&entries, &ss, i = 0](const auto w) mutable
                  { ss += fmt::format(":{:^{}}", entries[i++], w); });
    ss += ":\n";

    return ss;
}

} // Namespace anonymous

namespace Opm {

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
cumulative(const std::size_t reportStepNum) const
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

        this->outputCumulativeReportRecord_(values, names);
    }

    for (const auto& wname : this->schedule_.wellNames(reportStepNum)) {
        const auto& well = this->schedule_.getWell(wname, reportStepNum);

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

        this->outputCumulativeReportRecord_(values, names);
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

        for (const auto& phase : Inplace::phases()) 
            current_values[phase] = inplace.get(phase);
            
        Scalar field_dyn_pv = 0.0;


        for (auto nreg = inplace.max_region(name), reg = 0*nreg + 1; reg <= nreg; ++reg)       
            field_dyn_pv = field_dyn_pv + inplace.get(name, Inplace::Phase::DynamicPoreVolume, reg);
        
        current_values[Inplace::Phase::DynamicPoreVolume] = field_dyn_pv;
        
        this->fipUnitConvert_(current_values);
        this->outputResvFluidInPlace_(current_values, 0);
    }
    
    for (auto nreg = inplace.max_region(), reg = 0*nreg + 1; reg <= nreg; ++reg) {        
        std::unordered_map<Inplace::Phase, Scalar> current_values;
              
        for (const auto& phase : Inplace::phases()) {
            if (reg <= inplace.max_region(name))
                current_values[phase] = inplace.get(name, phase, reg);
            else    
                current_values[phase] = 0.0;
        }
        
        if (reg <= inplace.max_region(name))
            current_values[Inplace::Phase::DynamicPoreVolume] =
                inplace.get(name, Inplace::Phase::DynamicPoreVolume, reg);
        else
            current_values[Inplace::Phase::DynamicPoreVolume] = 0.0;
 
        this->fipUnitConvert_(current_values);
        this->outputResvFluidInPlace_(current_values, reg);
    }
    
    std::ostringstream ss;
    ss << " ===========================================================================================";
    OpmLog::note(ss.str());
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
injection(const std::size_t reportStepNum) const
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

        this->outputInjectionReportRecord_(values, names);
    }

    for (const auto& wname : this->schedule_.wellNames(reportStepNum)) {
        const auto& well = this->schedule_.getWell(wname, reportStepNum);

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

        this->outputInjectionReportRecord_(values, names);
    }

    this->endInjectionReport_();
}

template<class Scalar>
void LogOutputHelper<Scalar>::
production(const std::size_t reportStepNum) const
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
        } else {
            values[Ix::WatGasRatio] = values[Ix::WaterRate] / values[Ix::GasRate];
        }

        if (std::isnan(values[Ix::WatGasRatio])) {
            values[Ix::WatGasRatio] = 0.0;
        }

        this->outputProductionReportRecord_(values, names);
    }

    for (const auto& wname : this->schedule_.wellNames(reportStepNum)) {
        const auto& well = this->schedule_.getWell(wname, reportStepNum);

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

        if (values[Ix::WaterRate] == 0.0) {
            values[Ix::WatGasRatio] = 0.0;
        } else {
            values[Ix::WatGasRatio] = values[Ix::WaterRate] / values[Ix::GasRate];
        }

        if (std::isnan(values[Ix::WatGasRatio])) {
            values[Ix::WatGasRatio] = 0.0;
        }

        this->outputProductionReportRecord_(values, names);
    }

    this->endProductionReport_();
}

template <typename Scalar>
void LogOutputHelper<Scalar>::beginCumulativeReport_() const
{
    const auto unitType = this->eclState_.getUnits().getType();

    std::ostringstream ss;

    ss << "\n=================================================== CUMULATIVE PRODUCTION/INJECTION REPORT =========================================\n"
       << ":  WELL  :  LOCATION :  WELL  :CTRL:    OIL    :   WATER   :    GAS    :   Prod    :    OIL    :   WATER   :    GAS    :   INJ     :\n"
       << ":  NAME  :  (I,J,K)  :  TYPE  :MODE:    PROD   :   PROD    :    PROD   :  RES.VOL. :    INJ    :   INJ     :    INJ    :  RES.VOL. :\n";

    if (unitType == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
        ss << ":        :           :        :    :    MSCM   :   MSCM    :    MMSCM  :   MRCM    :    MSCM   :   MSCM    :    MMSCM  :   MRCM    :\n";
    } else if (unitType == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
        ss << ":        :           :        :    :    MSTB   :   MSTB    :    MMSCF  :   MRB     :    MSTB   :   MSTB    :    MMSCF  :   MRB     :\n";
    } else if (unitType == UnitSystem::UnitType::UNIT_TYPE_LAB) {
        ss << ":        :           :        :    :     MSCC  :   MSCC    :    MMSCC  :   MRCC    :    MSCC   :   MSCC    :    MMSCC  :   MRCC    :\n";
    }

    ss << "====================================================================================================================================";

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
                              const std::vector<std::string>& wellCumNames) const
{
    std::ostringstream ss;

    ss << std::right << std::fixed << std::setprecision(0) << ':'
       << std::setw(8) << wellCumNames[WellCumDataType::WellName] << ':';

    if (wellCum[WellCumDataType::WellLocationi] < 1) {
        ss << std::setw(11) <<  "" << ':';
    } else {
        ss << std::setw( 5) << wellCum[WellCumDataType::WellLocationi] << ','
           << std::setw( 5) << wellCum[WellCumDataType::WellLocationj] << ':';
    }

    auto scaledValue = [&wellCum](const typename WellCumDataType::WCId quantity)
    {
        // Unit M*
        return wellCum[quantity] / 1000.0;
    };

    auto scaledGasValue = [&wellCum](const typename WellCumDataType::WCId quantity)
    {
        // Unit MM*
        return wellCum[quantity] / (1000.0 * 1000.0);
    };

    ss << std::setw( 8) << wellCumNames[WellCumDataType::WCId::WellType]       << ':'
       << std::setw( 4) << wellCumNames[WellCumDataType::WCId::WellCTRL]       << ':' << std::setprecision(1)
       << std::setw(11) << scaledValue(WellCumDataType::WCId::OilProd)         << ':'
       << std::setw(11) << scaledValue(WellCumDataType::WCId::WaterProd)       << ':'
       << std::setw(11) << scaledGasValue(WellCumDataType::WCId::GasProd)      << ':'
       << std::setw(11) << scaledValue(WellCumDataType::WCId::FluidResVolProd) << ':'
       << std::setw(11) << scaledValue(WellCumDataType::WCId::OilInj)          << ':'
       << std::setw(11) << scaledValue(WellCumDataType::WCId::WaterInj)        << ':'
       << std::setw(11) << scaledGasValue(WellCumDataType::WCId::GasInj)       << ':'
       << std::setw(11) << scaledValue(WellCumDataType::WCId::FluidResVolInj)  << ':';

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
                             const std::vector<std::string>& wellInjNames) const
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
    ss << fmt::format("\n{:=^{}}\n", " PRODUCTION REPORT ", 129)
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
       << fmt::format("{:=>{}}", "", 129);

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
                              const std::vector<std::string>& wellProdNames) const
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
    std::ostringstream ss;

    ss << '\n';

    if (reg == 0) {
        ss << "                                                     ==================================================\n"
           << "                                                     :               FIELD TOTALS                     :\n";
    }
    else {
        ss << "                                                     ==================================================\n"
           << "                                                     :        " << name << " REPORT REGION  "
           << std::setw(8 - name.size()) << reg << "                :\n";
    }
    if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
        ss << "                                                     :         PAV = " << std::setw(14) << pav << " BARSA             :\n"
           << std::fixed << std::setprecision(0)
           << "                                                     :         PORV= " << std::setw(14) << cip[Inplace::Phase::PoreVolume] << "  RM3              :\n";
        if (!reg) {
            ss << "                                                     : Pressure is weighted by hydrocarbon pore volume:\n"
               << "                                                     : Porv volumes are taken at reference conditions :\n";
        }
        ss << "                           :--------------- OIL    SM3 ----------------:-- WAT    SM3 --:--------------- GAS    SM3  ---------------:\n";
    } else if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
        ss << std::fixed << std::setprecision(0) 
           << "                                                     :         PAV =" << std::setw(14) << pav << "  PSIA              :\n"
           << std::fixed << std::setprecision(0)
           << "                                                     :         PORV=" << std::setw(14) << cip[Inplace::Phase::PoreVolume] << "   RB               :\n";
        if (!reg) {
            ss << "                                                     : Pressure is weighted by hydrocarbon pore volume:\n"
               << "                                                     : Pore volumes are taken at reference conditions :\n";
        }
        ss << "                           :--------------- OIL    STB ----------------:-- WAT    STB --:--------------- GAS   MSCF ----------------:\n";
    }
    ss << "                           :      LIQUID        VAPOUR         TOTAL   :      TOTAL     :       FREE      DISSOLVED         TOTAL   :" << "\n"
       << " :-------------------------:-------------------------------------------:----------------:-------------------------------------------:" << "\n"
       << " :CURRENTLY IN PLACE       :" << std::setw(14) << cip[Inplace::Phase::OilInLiquidPhase]
       << std::setw(14) << cip[Inplace::Phase::OilInGasPhase]
       << std::setw(15) << cip[Inplace::Phase::OIL] << ":"
       << std::setw(14) << cip[Inplace::Phase::WATER] << "  :"
       << std::setw(14) << (cip[Inplace::Phase::GasInGasPhase])
       << std::setw(14) << cip[Inplace::Phase::GasInLiquidPhase]
       << std::setw(15) << cip[Inplace::Phase::GAS] << ":\n"
       << " :-------------------------:-------------------------------------------:----------------:-------------------------------------------:\n"
       << " :ORIGINALLY IN PLACE      :" << std::setw(14) << oip[Inplace::Phase::OilInLiquidPhase]
       << std::setw(14) << oip[Inplace::Phase::OilInGasPhase]
       << std::setw(15) << oip[Inplace::Phase::OIL] << ":"
       << std::setw(14) << oip[Inplace::Phase::WATER] << "  :"
       << std::setw(14) << oip[Inplace::Phase::GasInGasPhase]
       << std::setw(14) << oip[Inplace::Phase::GasInLiquidPhase]
       << std::setw(15) << oip[Inplace::Phase::GAS] << ":\n";
    
    if (reg == 0){
       ss << " ====================================================================================================================================\n\n";
    
    } else {
       ss << " :-------------------------:-------------------------------------------:----------------:-------------------------------------------:\n";
       ss << " ====================================================================================================================================\n\n";
    }
       
    OpmLog::note(ss.str());
}

template<class Scalar>
void LogOutputHelper<Scalar>::
outputResvFluidInPlace_(std::unordered_map<Inplace::Phase, Scalar> cipr,
                        const int reg) const
{
    const UnitSystem& units = eclState_.getUnits();
    std::ostringstream ss;

    if (reg == 0) {
        ss << "\n                                                     ===================================\n";
        if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
            ss << "                                                     :  RESERVOIR VOLUMES      RM3     :\n";
        } else if (units.getType() == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
            ss << "                                                     :  RESERVOIR VOLUMES      RB      :\n";
        }
        ss << " :---------:---------------:---------------:---------------:---------------:---------------:\n"
           << " : REGION  :  TOTAL PORE   :  PORE VOLUME  :  PORE VOLUME  : PORE VOLUME   :  PORE VOLUME  :\n"
           << " :         :   VOLUME      :  CONTAINING   :  CONTAINING   : CONTAINING    :  CONTAINING   :\n"
           << " :         :               :     OIL       :    WATER      :    GAS        :  HYDRO-CARBON :\n"
           << " :---------:---------------:---------------:---------------:---------------:---------------\n";

        ss << std::right << std::fixed << std::setprecision(0) << " :"
           << std::setw (8) <<  "FIELD" << " :";

    } else {
        ss << std::right << std::fixed << std::setprecision(0) << " :"
           << std::setw (8) <<  reg << " :";
    }    
        
    ss << std::setw(15) << cipr[Inplace::Phase::DynamicPoreVolume] << ":"
       << std::setw(15) << cipr[Inplace::Phase::OilResVolume] << ":"
       << std::setw(15) << cipr[Inplace::Phase::WaterResVolume] << ":"
       << std::setw(15) << cipr[Inplace::Phase::GasResVolume] << ":"
       << std::setw(15) << cipr[Inplace::Phase::OilResVolume] +
                           cipr[Inplace::Phase::GasResVolume] << ":";
    
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
