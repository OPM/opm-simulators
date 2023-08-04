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
 */
#ifndef LOG_OUTPUT_HELPER_HPP
#define LOG_OUTPUT_HELPER_HPP

#include <cstddef>
#include <functional>
#include <string>

namespace Opm {

class EclipseState;
class Schedule;
class SummaryState;

template<class Scalar>
class LogOutputHelper {
public:
    LogOutputHelper(const EclipseState& eclState,
                    const Schedule& schedule,
                    const SummaryState& st);

    //! \brief Write cumulative production and injection reports to output.
    void cumulative(const std::size_t reportStepNum,
                    std::function<bool(const std::string&)> isDefunct) const;

    //! \brief Write production report to output.
    void production(const std::size_t reportStepNum,
                    std::function<bool(const std::string&)> isDefunct) const;

private:
    void outputCumulativeReport_(const std::vector<Scalar>& wellCum,
                                 const std::vector<std::string>& wellCumNames) const;

    void outputProductionReport_(const std::vector<Scalar>& wellProd,
                                 const std::vector<std::string>& wellProdNames) const;

    struct WellCumDataType
    {
        enum WCId
        {
            WellLocationi = 0, // WLi
            WellLocationj = 1, // WLj
            OilProd = 2, // OP
            WaterProd = 3, // WP
            GasProd = 4, // GP
            FluidResVolProd = 5, // FRVP
            OilInj = 6, // OI
            WaterInj = 7, // WI
            GasInj = 8, // GI
            FluidResVolInj = 9, // FRVI
            WellName = 0, // WName
            WellType = 1, // WType
            WellCTRL = 2, // WCTRL
        };
        static constexpr int numWCValues = 10;
        static constexpr int numWCNames = 3;
    };

    struct WellProdDataType
    {
        enum WPId
        {
            WellLocationi = 0, // WLi
            WellLocationj = 1, // WLj
            OilRate = 2, // OR
            WaterRate = 3, // WR
            GasRate = 4, // GR
            FluidResVol = 5, // FRV
            WaterCut = 6, // WC
            GasOilRatio = 7, // GOR
            WatGasRatio = 8, // WGR
            BHP = 9, // BHP
            THP = 10, // THP
            SteadyStatePI = 11, // SteadyStatePI
            WellName = 0, // WName
            CTRLMode = 1, // CTRL
        };

        static constexpr int numWPValues = 12;
        static constexpr int numWPNames = 2;
    };

    const EclipseState& eclState_;
    const Schedule& schedule_;
    const SummaryState& summaryState_;
};

} // namespace Opm

#endif // LOG_OUTPUT_HELPER_HPP
