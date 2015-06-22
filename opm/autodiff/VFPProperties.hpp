/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_AUTODIFF_VFPPROPERTIES_HPP_
#define OPM_AUTODIFF_VFPPROPERTIES_HPP_

#include <opm/core/wells.h>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <boost/multi_array.hpp>

namespace Opm {

/**
 * Class which linearly interpolates BHP as a function of rate type, tubing head pressure,
 * water fraction, gas fraction, and artificial lift.
 */
class VFPProperties {
public:
    typedef boost::multi_array<double, 5> array_type;
    typedef boost::array<array_type::index, 5> extents;
    typedef AutoDiffBlock<double> ADB;

    ///Rate type
    enum FLO_TYPE {
        FLO_OIL, //< Oil rate
        FLO_LIQ, //< Liquid rate
        FLO_GAS, //< Gas rate
        //FLO_WG
        //FLO_TM
        FLO_INVALID
    };

    ///Water fraction variable
    enum WFR_TYPE {
        WFR_WOR, //< Water-oil ratio
        WFR_WCT, //< Water cut
        WFR_WGR, //< Water-gas ratio
        WFR_INVALID
    };

    ///Gas fraction variable
    enum GFR_TYPE {
        GFR_GOR, //< Gas-oil ratio
        GFR_GLR, //< Gas-liquid ratio
        GFR_OGR, //< Oil-gas ratio
        GFR_INVALID
    };

    ///Artificial lift quantity
    enum ALQ_TYPE {
        ALQ_GRAT, //< Lift as injection rate
        ALQ_IGLR, //< Injection gas-liquid ratio
        ALQ_TGLR, //< Total gas-liquid ratio
        ALQ_PUMP, //< Pump rating
        ALQ_COMP, //< Compressor power
        ALQ_BEAN, //< Choke diameter
        ALQ_UNDEF, //< Undefined
        ALQ_INVALID
    };

    /**
     * Constructor
     * @param table_num VFP table number
     * @param datum_depth Reference depth for BHP
     * @param flo_type Specifies what flo_data represents
     * @param wfr_type Specifies what wfr_data represents
     * @param gfr_type Specifies what gfr_data represents
     * @param alq_type Specifies what alq_data represents
     * @param flo_data Axis for flo_type
     * @param thp_data Axis for thp_type
     * @param wfr_data Axis for wfr_type
     * @param gfr_data Axis for gfr_type
     * @param alq_data Axis for alq_type
     * @param data BHP to be interpolated. Given as a 5D array so that
     *        BHP = data[thp][wfr][gfr][alq][flo] for the indices thp, wfr, etc.
     */
    VFPProperties(int table_num,
        double datum_depth,
        FLO_TYPE flo_type,
        WFR_TYPE wfr_type,
        GFR_TYPE gfr_type,
        ALQ_TYPE alq_type,
        const std::vector<double>& flo_data,
        const std::vector<double>& thp_data,
        const std::vector<double>& wfr_data,
        const std::vector<double>& gfr_data,
        const std::vector<double>& alq_data,
        array_type data
        ) :
            table_num_(table_num),
            datum_depth_(datum_depth),
            flo_type_(flo_type),
            wfr_type_(wfr_type),
            gfr_type_(gfr_type),
            alq_type_(alq_type),
            flo_data_(flo_data),
            thp_data_(thp_data),
            wfr_data_(wfr_data),
            gfr_data_(gfr_data),
            alq_data_(alq_data),
            data_(data) {

    }

    /**
     * Constructor which parses a deck keyword and retrieves the relevant parts for a
     * VFP table.
     */
    VFPProperties(DeckKeywordConstPtr table);

    /**
     * Linear interpolation of bhp as function of the input parameters.
     * @param wells Wells structure with information about wells in qs
     * @param qs Flow quantities
     * @param thp Tubing head pressure
     * @param alq Artificial lift or other parameter
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    ADB bhp(const Wells& wells, const ADB& qs, const ADB& thp, const ADB& alq);

    /**
     * Linear interpolation of bhp as a function of the input parameters
     * @param flo Production rate of oil, gas or liquid
     * @param thp Tubing head pressure
     * @param wfr Water-oil ratio, water cut, or water-gas ratio
     * @param gfr Gas-oil ratio, gas-liquid ratio, or oil-gas ratio
     * @param alq Artificial lift or other parameter
     *
     * @return The bottom hole pressure, interpolated/extrapolated linearly using
     * the above parameters from the values in the input table.
     */
    double bhp(const double& flo, const double& thp, const double& wfr, const double& gfr, const double& alq);

    ADB bhp(const ADB& flo, const ADB& thp, const ADB& wfr, const ADB& gfr, const ADB& alq);

private:
    /**
     * Helper struct for linear interpolation
     */
    struct InterpData {
        InterpData() : ind_{0, 0}, factor_(0.0) {}
        int ind_[2]; //[First element greater than or equal to value, Last element smaller than or equal to value]
        double factor_; //Interpolation factor
    };

    /**
     * Helper function to find indices etc. for linear interpolation
     */
    static InterpData find_interp_data(const double& value, const std::vector<double>& values);

    /**
     * Helper function which interpolates data using the indices etc. given in the inputs.
     */
    double interpolate(const InterpData& flo_i, const InterpData& thp_i,
            const InterpData& wfr_i, const InterpData& gfr_i, const InterpData& alq_i);

    //"Header" variables
    int table_num_;
    double datum_depth_;
    FLO_TYPE flo_type_;
    WFR_TYPE wfr_type_;
    GFR_TYPE gfr_type_;
    ALQ_TYPE alq_type_;

    //The actual table axes
    std::vector<double> flo_data_;
    std::vector<double> thp_data_;
    std::vector<double> wfr_data_;
    std::vector<double> gfr_data_;
    std::vector<double> alq_data_;

    //The data itself
    array_type data_;
};

}

#endif /* OPM_AUTODIFF_VFPPROPERTIES_HPP_ */
