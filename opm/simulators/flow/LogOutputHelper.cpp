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

#include <opm/simulators/utils/PressureAverage.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

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

} // Namespace anonymous

namespace Opm {

template<class Scalar>
LogOutputHelper<Scalar>::LogOutputHelper(const EclipseState& eclState,
                                         const Schedule& schedule,
                                         const SummaryState& summaryState)
    : eclState_(eclState)
    , schedule_(schedule)
    , summaryState_(summaryState)
{ 
    flowVersionName_ = moduleVersionName();
}

template<class Scalar>
void LogOutputHelper<Scalar>::
cumulative(const std::size_t reportStepNum) const
{
    this->beginCumulativeReport_();

    std::vector<Scalar> tmp_values(WellCumDataType::numWCValues, 0.0);
    std::vector<std::string> tmp_names(WellCumDataType::numWCNames, "");

    const auto& st = summaryState_;
    for (const auto& gname : schedule_.groupNames()) {
        auto gName = static_cast<std::string>(gname);
        auto get = [&st, &gName](const std::string& vector)
        {
            const auto key = vector + ':' + gName;
            return st.has(key) ? st.get(key) : 0.0;
        };

        tmp_names[0] = gname;

        if (tmp_names[0] == "FIELD") {
            tmp_values[2] = st.get("FOPT", 0.0); // WellCumDataType::OilProd
            tmp_values[3] = st.get("FWPT", 0.0); // WellCumDataType::WaterProd
            tmp_values[4] = st.get("FGPT", 0.0); // WellCumDataType::GasProd
            tmp_values[5] = st.get("FVPT", 0.0); // WellCumDataType::FluidResVolProd
            tmp_values[6] = st.get("FOIT", 0.0); // WellCumDataType::OilInj
            tmp_values[7] = st.get("FWIT", 0.0); // WellCumDataType::WaterInj
            tmp_values[8] = st.get("FGIT", 0.0); // WellCumDataType::GasInj
            tmp_values[9] = st.get("FVIT", 0.0); // WellCumDataType::FluidResVolInj
        } else {
            tmp_values[2] = get("GOPT"); // WellCumDataType::OilProd
            tmp_values[3] = get("GWPT"); // WellCumDataType::WaterProd
            tmp_values[4] = get("GGPT"); // WellCumDataType::GasProd
            tmp_values[5] = get("GVPT"); // WellCumDataType::FluidResVolProd
            tmp_values[6] = get("GOIT"); // WellCumDataType::OilInj
            tmp_values[7] = get("GWIT"); // WellCumDataType::WaterInj
            tmp_values[8] = get("GGIT"); // WellCumDataType::GasInj
            tmp_values[9] = get("GVIT"); // WellCumDataType::FluidResVolInj
        }

        this->outputCumulativeReportRecord_(tmp_values, tmp_names);
    }

    for (const auto& wname : schedule_.wellNames(reportStepNum)) {
        const auto& well = schedule_.getWell(wname, reportStepNum);
        tmp_names[0] = wname; // WellCumDataType::WellName
        auto wName = static_cast<std::string>(wname);
        auto get = [&st, &wName](const std::string& vector)
        {
            const auto key = vector + ':' + wName;
            return st.has(key) ? st.get(key) : 0.0;
        };

        if (well.isInjector()) {
            const auto& controls = well.injectionControls(st);
            const auto ctlMode = controls.cmode;
            const auto injType = controls.injector_type;
            using CMode = ::Opm::Well::InjectorCMode;
            using WType = ::Opm::InjectorType;

            auto ftype = [](const auto wtype) -> std::string
            {
                switch (wtype) {
                case WType::OIL:   return "Oil";
                case WType::WATER: return "Wat";
                case WType::GAS:   return "Gas";
                case WType::MULTI: return "Multi";
                default: return "";
                }
            };

            auto fctl = [](const auto wmctl) -> std::string
            {
                switch (wmctl) {
                case CMode::RATE: return "RATE";
                case CMode::RESV: return "RESV";
                case CMode::THP:  return "THP";
                case CMode::BHP:  return "BHP";
                case CMode::GRUP: return "GRUP";
                default: return "";
                }
            };

            tmp_names[1] = "INJ"; // WellCumDataType::WellType
            const auto flowctl = fctl(ctlMode);
            if (flowctl == "RATE") { // WellCumDataType::WellCTRL
                const auto flowtype = ftype(injType);
                if (flowtype == "Oil") {
                    tmp_names[2] = "ORAT";
                } else if (flowtype == "Wat") {
                    tmp_names[2] = "WRAT";
                } else if (flowtype == "Gas") {
                    tmp_names[2] = "GRAT";
                }
            } else {
                tmp_names[2] = flowctl;
            }
        } else if (well.isProducer()) {
            const auto& controls = well.productionControls(st);
            using CMode = ::Opm::Well::ProducerCMode;

            auto fctl = [](const auto wmctl) -> std::string
            {
                switch (wmctl) {
                case CMode::ORAT: return "ORAT";
                case CMode::WRAT: return "WRAT";
                case CMode::GRAT: return "GRAT";
                case CMode::LRAT: return "LRAT";
                case CMode::RESV: return "RESV";
                case CMode::THP:  return "THP";
                case CMode::BHP:  return "BHP";
                case CMode::CRAT: return "CRAT";
                case CMode::GRUP: return "GRUP";
                default: return "none";
                }
            };
            tmp_names[1] = "PROD"; // WellProdDataType::CTRLMode
            tmp_names[2] = fctl(controls.cmode); // WellProdDataType::CTRLMode
        }

        tmp_values[0] = well.getHeadI() + 1; // WellCumDataType::wellLocationi
        tmp_values[1] = well.getHeadJ() + 1; // WellCumDataType::wellLocationj
        tmp_values[2] = get("WOPT"); // WellCumDataType::OilProd
        tmp_values[3] = get("WWPT"); // WellCumDataType::WaterProd
        tmp_values[4] = get("WGPT"); // WellCumDataType::GasProd
        tmp_values[5] = get("WVPT"); // WellCumDataType::FluidResVolProd
        tmp_values[6] = get("WOIT"); // WellCumDataType::OilInj
        tmp_values[7] = get("WWIT"); // WellCumDataType::WaterInj
        tmp_values[8] = get("WGIT"); // WellCumDataType::GasInj
        tmp_values[9] = get("WVIT"); // WellCumDataType::FluidResVolInj

        this->outputCumulativeReportRecord_(tmp_values, tmp_names);
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

    std::vector<Scalar> tmp_values(WellInjDataType::numWIValues, 0.0);
    std::vector<std::string> tmp_names(WellInjDataType::numWINames, "");

    const auto& st = summaryState_;
    for (const auto& gname : schedule_.groupNames()) {
        auto gName = static_cast<std::string>(gname);
        auto get = [&st, &gName](const std::string& vector)
        {
            const auto key = vector + ':' + gName;
            return st.has(key) ? st.get(key) : 0.0;
        };

        tmp_names[0] = gname;
        if (tmp_names[0] == "FIELD") {
            tmp_values[2] = st.get("FOIR", 0.0); // WellInjDataType::OilRate
            tmp_values[3] = st.get("FWIR", 0.0); // WellInjDataType::WaterRate
            tmp_values[4] = st.get("FGIR", 0.0); // WellInjDataType::GasRate
            tmp_values[5] = st.get("FVIR", 0.0); // WellInjDataType::FluidResVol
        } else {
            tmp_values[2] = get("GOIR"); // WellInjDataType::OilRate
            tmp_values[3] = get("GWIR"); // WellInjDataType::WaterRate
            tmp_values[4] = get("GGIR"); // WellInjDataType::GasRate
            tmp_values[5] = get("GVIR"); // WellInjDataType::FluidResVol
        }

        this->outputInjectionReportRecord_(tmp_values, tmp_names);
    }

    for (const auto& wname : schedule_.wellNames(reportStepNum)) {
        const auto& well = schedule_.getWell(wname, reportStepNum);

        // Ignore Producer wells
        if (well.isProducer()) {
            continue;
        }

        tmp_names[0] = wname; // WellInjDataType::WellName
        auto wName = static_cast<std::string>(wname);
        auto get = [&st, &wName](const std::string& vector)
        {
            const auto key = vector + ':' + wName;
            return st.has(key) ? st.get(key) : 0.0;
        };

        const auto& controls = well.injectionControls(st);
        const auto ctlMode = controls.cmode;
        const auto injType = controls.injector_type;
        using CMode = Well::InjectorCMode;
        using WType = InjectorType;

        auto ftype = [](const auto wtype) -> std::string
        {
            switch (wtype) {
            case WType::OIL:   return "Oil";
            case WType::WATER: return "Wat";
            case WType::GAS:   return "Gas";
            case WType::MULTI: return "Multi";
            default: return "";
            }
        };

        auto fctl = [](const auto wmctl) -> std::string
        {
            switch (wmctl) {
            case CMode::RATE: return "RATE";
            case CMode::RESV: return "RESV";
            case CMode::THP:  return "THP";
            case CMode::BHP:  return "BHP";
            case CMode::GRUP: return "GRUP";
            default: return "";
            }
        };

        const auto flowtype = ftype(injType);
        const auto flowctl = fctl(ctlMode);
        if (flowtype == "Oil") { // WellInjDataType::CTRLModeOil
            if (flowctl == "RATE") {
                tmp_names[1] = "ORAT";
            } else {
                tmp_names[1] =  flowctl;
            }
        }
        else if (flowtype == "Wat") { // WellInjDataType::CTRLModeWat
            if (flowctl == "RATE") {
                tmp_names[3] = "WRAT";
            } else {
                tmp_names[2] =  flowctl;
            }
        }
        else if (flowtype == "Gas") // WellInjDataType::CTRLModeGas
        {
            if (flowctl == "RATE") {
                tmp_names[3] = "GRAT";
            } else {
                tmp_names[3] =  flowctl;
            }
        }

        tmp_values[0] = well.getHeadI() + 1; // WellInjDataType::wellLocationi
        tmp_values[1] = well.getHeadJ() + 1; // WellInjDataType::wellLocationj
        tmp_values[2] = get("WOIR"); // WellInjDataType::OilRate
        tmp_values[3] = get("WWIR"); // WellInjDataType::WaterRate
        tmp_values[4] = get("WGIR"); // WellInjDataType::GasRate
        tmp_values[5] = get("WVIR");// WellInjDataType::FluidResVol
        tmp_values[6] = get("WBHP"); // WellInjDataType::BHP
        tmp_values[7] = get("WTHP"); // WellInjDataType::THP
        //tmp_values[8] = 0; // WellInjDataType::SteadyStateII

        this->outputInjectionReportRecord_(tmp_values, tmp_names);
    }

    this->endInjectionReport_();
}

template<class Scalar>
void LogOutputHelper<Scalar>::
production(const std::size_t reportStepNum) const
{
    this->beginProductionReport_();

    std::vector<Scalar> tmp_values(WellProdDataType::numWPValues, 0.0);
    std::vector<std::string> tmp_names(WellProdDataType::numWPNames, "");

    const auto& st = summaryState_;
    for (const auto& gname : schedule_.groupNames()) {
        auto gName = static_cast<std::string>(gname);
        auto get = [&st, &gName](const std::string& vector)
        {
            const auto key = vector + ':' + gName;
            return st.has(key) ? st.get(key) : 0.0;
        };

        tmp_names[0] = gname;
        if (tmp_names[0] == "FIELD") {
            tmp_values[2] = st.get("FOPR", 0.0); // WellProdDataType::OilRate
            tmp_values[3] = st.get("FWPR", 0.0); // WellProdDataType::WaterRate
            tmp_values[4] = st.get("FGPR", 0.0); // WellProdDataType::GasRate
            tmp_values[5] = st.get("FVPR", 0.0); // WellProdDataType::FluidResVol
            tmp_values[6] = st.get("FWCT", 0.0); // WellProdDataType::WaterCut
            tmp_values[7] = st.get("FGOR", 0.0); // WellProdDataType::GasOilRatio
        } else {
            tmp_values[2] = get("GOPR"); // WellProdDataType::OilRate
            tmp_values[3] = get("GWPR"); // WellProdDataType::WaterRate
            tmp_values[4] = get("GGPR"); // WellProdDataType::GasRate
            tmp_values[5] = get("GVPR"); // WellProdDataType::FluidResVol
            tmp_values[6] = get("GWCT"); // WellProdDataType::WaterCut
            tmp_values[7] = get("GGOR"); // WellProdDataType::GasOilRatio
        }

        tmp_values[8] = tmp_values[3] / tmp_values[4]; // WellProdDataType::WaterGasRatio
        if (std::isnan(tmp_values[8])) {
            tmp_values[8] = 0.0;
        }

        this->outputProductionReportRecord_(tmp_values, tmp_names);
    }

    for (const auto& wname : schedule_.wellNames(reportStepNum)) {
        const auto& well = schedule_.getWell(wname, reportStepNum);

        // Ignore injector wells
        if (well.isInjector()) {
            continue;
        }

        tmp_names[0] = wname; // WellProdDataType::WellName


        auto wName = static_cast<std::string>(wname);
        auto get = [&st, &wName](const std::string& vector)
        {
            const auto key = vector + ':' + wName;
            return st.has(key) ? st.get(key) : 0.0;
        };

        const auto& controls = well.productionControls(st);
        using CMode = Well::ProducerCMode;

        auto fctl = [](const auto wmctl) -> std::string
        {
            switch (wmctl) {
            case CMode::ORAT: return "ORAT";
            case CMode::WRAT: return "WRAT";
            case CMode::GRAT: return "GRAT";
            case CMode::LRAT: return "LRAT";
            case CMode::RESV: return "RESV";
            case CMode::THP:  return "THP";
            case CMode::BHP:  return "BHP";
            case CMode::CRAT: return "CRate";
            case CMode::GRUP: return "GRUP";
            default: return "none";
            }
        };

        tmp_names[1] = fctl(controls.cmode); // WellProdDataType::CTRLMode

        tmp_values[0] = well.getHeadI() + 1; // WellProdDataType::WellLocationi
        tmp_values[1] = well.getHeadJ() + 1; // WellProdDataType::WellLocationj
        tmp_values[2] = get("WOPR"); // WellProdDataType::OilRate
        tmp_values[3] = get("WWPR"); // WellProdDataType::WaterRate
        tmp_values[4] = get("WGPR"); // WellProdDataType::GasRate
        tmp_values[5] = get("WVPR"); // WellProdDataType::FluidResVol
        tmp_values[6] = get("WWCT"); // WellProdDataType::WaterCut
        tmp_values[7] = get("WGOR"); // WellProdDataType::GasOilRatio
        tmp_values[9] = get("WBHP"); // WellProdDataType::BHP
        tmp_values[10] = get("WTHP"); // WellProdDataType::THP
        //tmp_values[11] = 0; //WellProdDataType::SteadyStatePI //

        tmp_values[8] = tmp_values[3] / tmp_values[4]; // WellProdDataType::WaterGasRatio
        if (std::isnan(tmp_values[8])) {
            tmp_values[8] = 0.0;
        }

        this->outputProductionReportRecord_(tmp_values, tmp_names);
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
    const auto ss = std::string { ":--------:-----------:--------:----:-----------:-----------:-----------:-----------:-----------:-----------:-----------:-----------:" };

    OpmLog::note(ss);
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
    const auto unitType = this->eclState_.getUnits().getType();

    std::ostringstream ss;
    ss << "\n=================================================== INJECTION REPORT ========================================\n"//===================== \n"
       << ":  WELL  :  LOCATION : CTRL : CTRL : CTRL :    OIL    :   WATER   :    GAS    :   FLUID   : BHP OR : THP OR :\n"// STEADY-ST II       :\n"
       << ":  NAME  :  (I,J,K)  : MODE : MODE : MODE :    RATE   :   RATE    :    RATE   :  RES.VOL. : CON.PR.: BLK.PR.:\n";// OR POTENTIAL       :\n";

    if (unitType == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
        ss << ":        :           : OIL  : WAT  : GAS  :  SCM/DAY  :  SCM/DAY  :  SCM/DAY  :  RCM/DAY  :  BARSA :  BARSA :\n";//                    :\n";
    } else if (unitType == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
        ss << ":        :           : OIL  : WAT  : GAS  :  STB/DAY  :  STB/DAY  :  MSCF/DAY :  RB/DAY   :  PSIA  :  PSIA  :\n";//                    :\n";
    } else if (unitType == UnitSystem::UnitType::UNIT_TYPE_LAB) {
        ss << ":        :           : OIL  : WAT  : GAS  :   SCC/HR  :  SCC/HR   :  SCC/HR   :  RCC/HR   :  ATMA  :  ATMA  :\n";//                    :\n";
    }

    ss << "=============================================================================================================";//=====================";

    OpmLog::note(ss.str());
}

template <typename Scalar>
void LogOutputHelper<Scalar>::endInjectionReport_() const
{
    const auto ss = std::string { ":--------:-----------:------:------:------:-----------:-----------:-----------:-----------:--------:--------:" };

    OpmLog::note(ss);
}

template<class Scalar>
void LogOutputHelper<Scalar>::
outputInjectionReportRecord_(const std::vector<Scalar>& wellInj,
                             const std::vector<std::string>& wellInjNames) const
{
    const auto isWellRecord =
        wellInj[WellProdDataType::WellLocationi] >= 1;

    std::ostringstream ss;

    ss << std::right << std::fixed << std::setprecision(0) << ':'
       << std::setw(8) << wellInjNames[WellInjDataType::WellName] << ':';

    if (! isWellRecord) {
        ss << std::setw(11) << "" << ':';
    } else {
        ss << std::setw( 5) << wellInj[WellInjDataType::WellLocationi] << ','
           << std::setw( 5) << wellInj[WellInjDataType::WellLocationj] << ':';
    }

    ss << std::setw( 6) << wellInjNames[WellInjDataType::CTRLModeOil] << ':'
       << std::setw( 6) << wellInjNames[WellInjDataType::CTRLModeWat] << ':'
       << std::setw( 6) << wellInjNames[WellInjDataType::CTRLModeGas] << ':' << std::setprecision(1)
       << std::setw(11) << wellInj[WellInjDataType::OilRate]          << ':'
       << std::setw(11) << wellInj[WellInjDataType::WaterRate]        << ':'
       << std::setw(11) << wellInj[WellInjDataType::GasRate]          << ':'
       << std::setw(11) << wellInj[WellInjDataType::FluidResVol]      << ':';

    if (! isWellRecord) {
        ss << std::setw(8) << "" << ':' << std::setw(8) << "" << ':'; //wellInj[WellInjDataType::SteadyStateII] << std::setw(10) << "\n"
    } else {
        ss << std::setw(8) << wellInj[WellInjDataType::BHP] << ':'
           << std::setw(8) << wellInj[WellInjDataType::THP] << ':'; //wellInj[WellInjDataType::SteadyStateII] << std::setw(10) << "\n"
    }

    OpmLog::note(ss.str());
}

template <typename Scalar>
void LogOutputHelper<Scalar>::beginProductionReport_() const
{
    const auto unitType = this->eclState_.getUnits().getType();

    std::ostringstream ss;
    ss << "\n======================================================= PRODUCTION REPORT =======================================================\n"//=================== \n"
       << ":  WELL  :  LOCATION :CTRL:    OIL    :   WATER   :    GAS    :   FLUID   :   WATER   : GAS/OIL  :  WAT/GAS   : BHP OR : THP OR :\n"// STEADY-ST PI       :\n"
       << ":  NAME  :  (I,J,K)  :MODE:    RATE   :   RATE    :    RATE   :  RES.VOL. :    CUT    :  RATIO   :   RATIO    : CON.PR.: BLK.PR.:\n";// OR POTN OF PREF. PH:\n";

    if (unitType == UnitSystem::UnitType::UNIT_TYPE_METRIC) {
        ss << ":        :           :    :  SCM/DAY  :  SCM/DAY  :  SCM/DAY  :  RCM/DAY  :  SCM/SCM  :  SCM/SCM :  SCM/SCM   :  BARSA :  BARSA :\n";//                    :\n";
    } else if (unitType == UnitSystem::UnitType::UNIT_TYPE_FIELD) {
        ss << ":        :           :    :  STB/DAY  :  STB/DAY  :  MSCF/DAY :  RB/DAY   :           : MSCF/STB :  STB/MSCF  :  PSIA  :  PSIA  :\n";//                    :\n";
    } else if (unitType == UnitSystem::UnitType::UNIT_TYPE_LAB) {
        ss << ":        :           :    :  SCC/HR   :  SCC/HR   :  SCC/HR   :    RCC    :  SCC/SCC  :  SCC/SCC :  SCC/SCC   :  ATMA  :  ATMA  :\n";//                    :\n";
    }

    ss << "=================================================================================================================================";//===================";

    OpmLog::note(ss.str());
}

template <typename Scalar>
void LogOutputHelper<Scalar>::endProductionReport_() const
{
    std::ostringstream ss;

    ss << ':' << std::setfill ('-') << std::setw (9) << ':'
       << std::setfill ('-') << std::setw (12) << ':'
       << std::setfill ('-') << std::setw ( 5) << ':'
       << std::setfill ('-') << std::setw (12) << ':'
       << std::setfill ('-') << std::setw (12) << ':'
       << std::setfill ('-') << std::setw (12) << ':'
       << std::setfill ('-') << std::setw (12) << ':'
       << std::setfill ('-') << std::setw (12) << ':'
       << std::setfill ('-') << std::setw (11) << ':'
       << std::setfill ('-') << std::setw (13) << ':'
       << std::setfill ('-') << std::setw ( 9) << ':'
       << std::setfill ('-') << std::setw ( 9) << ':';

    OpmLog::note(ss.str());
}

template<class Scalar>
void LogOutputHelper<Scalar>::
outputProductionReportRecord_(const std::vector<Scalar>& wellProd,
                              const std::vector<std::string>& wellProdNames) const
{
    const auto isWellRecord =
        wellProd[WellProdDataType::WellLocationi] >= 1;

    std::ostringstream ss;

    ss << std::right << std::fixed << ':'
       << std::setw(8) << wellProdNames[WellProdDataType::WellName] << ':';

    if (! isWellRecord) {
        ss << std::setprecision(0) << std::setw(11) << "" << ':';
    } else {
        ss << std::setprecision(0)
           << std::setw(5) << wellProd[WellProdDataType::WellLocationi] << ','
           << std::setw(5) << wellProd[WellProdDataType::WellLocationj] << ':';
    }

    ss << std::setw( 4) << wellProdNames[WellProdDataType::CTRLMode] << ':' << std::setprecision(1)
       << std::setw(11) << wellProd[WellProdDataType::OilRate]       << ':'
       << std::setw(11) << wellProd[WellProdDataType::WaterRate]     << ':'
       << std::setw(11) << wellProd[WellProdDataType::GasRate]       << ':'
       << std::setw(11) << wellProd[WellProdDataType::FluidResVol]   << ':' << std::setprecision(3)
       << std::setw(11) << wellProd[WellProdDataType::WaterCut]      << ':' << std::setprecision(2)
       << std::setw(10) << wellProd[WellProdDataType::GasOilRatio]   << ':' << std::setprecision(4)
       << std::setw(12) << wellProd[WellProdDataType::WatGasRatio]   << ':' << std::setprecision(1);

    if (! isWellRecord) {
        ss << std::setw(8) << "" << ':' << std::setw(8) << "" << ':';
    } else {
        ss << std::setw(8) << wellProd[WellProdDataType::BHP] << ':'
           << std::setw(8) << wellProd[WellProdDataType::THP] << ':';
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
        {Inplace::Phase::WATER,             M::liquid_surface_volume},
        {Inplace::Phase::OIL,               M::liquid_surface_volume},
        {Inplace::Phase::OilInLiquidPhase,  M::liquid_surface_volume},
        {Inplace::Phase::OilInGasPhase,     M::liquid_surface_volume},
        {Inplace::Phase::GAS,               M::gas_surface_volume},
        {Inplace::Phase::GasInLiquidPhase,  M::gas_surface_volume},
        {Inplace::Phase::GasInGasPhase,     M::gas_surface_volume},
        {Inplace::Phase::PoreVolume,        M::volume},
        {Inplace::Phase::DynamicPoreVolume, M::volume},
        {Inplace::Phase::WaterResVolume,    M::volume},
        {Inplace::Phase::OilResVolume,      M::volume},
        {Inplace::Phase::GasResVolume,      M::volume},
        {Inplace::Phase::SALT,              M::mass},
        {Inplace::Phase::CO2InWaterPhase,   M::moles},
        {Inplace::Phase::CO2InGasPhaseInMob,M::moles},
        {Inplace::Phase::CO2InGasPhaseMob,  M::moles},
        {Inplace::Phase::WaterInWaterPhase, M::liquid_surface_volume},
        {Inplace::Phase::WaterInGasPhase,   M::liquid_surface_volume},
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

} // namespace Opm
