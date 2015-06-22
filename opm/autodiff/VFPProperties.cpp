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

#include "config.h"


#include <opm/autodiff/VFPProperties.hpp>

#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <algorithm>

namespace Opm {

VFPProperties::VFPProperties(DeckKeywordConstPtr table) {
    auto iter = table->begin();

    auto header = (*iter++);
    table_num_ = header->getItem("TABLE")->getInt(0);
    datum_depth_ = header->getItem("DATUM_DEPTH")->getRawDouble(0);

    //Rate type
    try {
        std::string flo_string = header->getItem("RATE_TYPE")->getString(0);
        if (flo_string == "OIL") {
            flo_type_ = FLO_OIL;
        }
        else if (flo_string == "LIQ") {
            flo_type_ = FLO_LIQ;
        }
        else if (flo_string == "GAS") {
            flo_type_ = FLO_GAS;
        }
        else {
            flo_type_ = FLO_INVALID;
        }
    }
    catch (std::invalid_argument& e) {
        //TODO: log here
        flo_type_ = FLO_INVALID;
    }

    //Water fraction
    try {
        std::string wfr_string = header->getItem("WFR")->getString(0);
        if (wfr_string == "WOR") {
            wfr_type_ = WFR_WOR;
        }
        else if (wfr_string == "WCT") {
            wfr_type_ = WFR_WCT;
        }
        else if (wfr_string == "WGR") {
            wfr_type_ = WFR_WGR;
        }
        else {
            wfr_type_ = WFR_INVALID;
        }
    }
    catch (std::invalid_argument& e) {
        //TODO: log here
        wfr_type_ = WFR_INVALID;
    }

    //Gas fraction
    try {
        std::string gfr_string = header->getItem("GFR")->getString(0);
        if (gfr_string == "GOR") {
            gfr_type_ = GFR_GOR;
        }
        else if (gfr_string == "GLR") {
            gfr_type_ = GFR_GLR;
        }
        else if (gfr_string == "OGR") {
            gfr_type_ = GFR_OGR;
        }
        else {
            gfr_type_ = GFR_INVALID;
        }
    }
    catch (std::invalid_argument& e) {
        //TODO: log here
        gfr_type_ = GFR_INVALID;
    }

    //Artificial lift
    try {
        std::string alq_string = header->getItem("ALQ")->getString(0);
        if (alq_string == "GRAT") {
            alq_type_ = ALQ_GRAT;
        }
        else if (alq_string == "IGLR") {
            alq_type_ = ALQ_IGLR;
        }
        else if (alq_string == "TGLR") {
            alq_type_ = ALQ_TGLR;
        }
        else if (alq_string == "PUMP") {
            alq_type_ = ALQ_PUMP;
        }
        else if (alq_string == "COMP") {
            alq_type_ = ALQ_COMP;
        }
        else if (alq_string == "BEAN") {
            alq_type_ = ALQ_BEAN;
        }
        else if (alq_string == "UNDEF") {
            alq_type_ = ALQ_UNDEF;
        }
        else {
            alq_type_ = ALQ_INVALID;
        }
    }
    catch (std::invalid_argument& e) {
        //TODO: log here
        //FIXME:  header->getItem("ALQ")->getString(0); appears to fail, ignoring in current implementation
        //alq_type_ = ALQ_INVALID;
        alq_type_ = ALQ_UNDEF;
    }

    //Get actual rate / flow values
    flo_data_ = (*iter++)->getItem("FLOW_VALUES")->getRawDoubleData();

    //Get actual tubing head pressure values
    thp_data_ = (*iter++)->getItem("THP_VALUES")->getRawDoubleData();

    //Get actual water fraction values
    wfr_data_ = (*iter++)->getItem("WFR_VALUES")->getRawDoubleData();

    //Get actual gas fraction values
    gfr_data_ = (*iter++)->getItem("GFR_VALUES")->getRawDoubleData();

    //Get actual gas fraction values
    alq_data_ = (*iter++)->getItem("ALQ_VALUES")->getRawDoubleData();

    //Finally, read the actual table itself.
    size_t nt = thp_data_.size();
    size_t nw = wfr_data_.size();
    size_t ng = gfr_data_.size();
    size_t na = alq_data_.size();
    size_t nf = flo_data_.size();
    extents shape;
    shape[0] = nt;
    shape[1] = nw;
    shape[2] = ng;
    shape[3] = na;
    shape[4] = nf;
    data_.resize(shape);

    for (; iter!=table->end(); ++iter) {
        //Get indices (subtract 1 to get 0-based index)
        int t = (*iter)->getItem("THP_INDEX")->getInt(0) - 1;
        int w = (*iter)->getItem("WFR_INDEX")->getInt(0) - 1;
        int g = (*iter)->getItem("GFR_INDEX")->getInt(0) - 1;
        int a = (*iter)->getItem("ALQ_INDEX")->getInt(0) - 1;

        //Rest of values (bottom hole pressure or tubing head temperature) have index of flo value
        const std::vector<double>& bhp_tht = (*iter)->getItem("VALUES")->getRawDoubleData();
        std::copy(bhp_tht.begin(), bhp_tht.end(), &data_[t][w][g][a][0]);
    }

    check();
}


VFPProperties::VFPProperties(int table_num,
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
    check();
}


void VFPProperties::check() {
    //Table number
    assert(table_num_ > 0);

    //Misc types
    assert(flo_type_ >= FLO_OIL && flo_type_ < FLO_INVALID);
    assert(wfr_type_ >= WFR_WOR && flo_type_ < WFR_INVALID);
    assert(gfr_type_ >= GFR_GOR && flo_type_ < GFR_INVALID);
    assert(alq_type_ >= ALQ_GRAT && flo_type_ < ALQ_INVALID);

    //Data axis size
    assert(flo_data_.size() > 0);
    assert(thp_data_.size() > 0);
    assert(wfr_data_.size() > 0);
    assert(gfr_data_.size() > 0);
    assert(alq_data_.size() > 0);

    //Data axis sorted?
    assert(is_sorted(flo_data_.begin(), flo_data_.end()));
    assert(is_sorted(thp_data_.begin(), thp_data_.end()));
    assert(is_sorted(wfr_data_.begin(), wfr_data_.end()));
    assert(is_sorted(gfr_data_.begin(), gfr_data_.end()));
    assert(is_sorted(alq_data_.begin(), alq_data_.end()));

    //Check data size matches axes
    assert(data_.num_dimensions() == 5);
    assert(data_.shape()[0] == thp_data_.size());
    assert(data_.shape()[1] == wfr_data_.size());
    assert(data_.shape()[2] == gfr_data_.size());
    assert(data_.shape()[3] == alq_data_.size());
    assert(data_.shape()[4] == flo_data_.size());

    //Finally, check that all data is within reasonable ranges, defined to be up-to 1.0e10...
    for (int t=0; t<data_.shape()[0]; ++t) {
        for (int w=0; w<data_.shape()[1]; ++w) {
            for (int g=0; g<data_.shape()[2]; ++g) {
                for (int a=0; a<data_.shape()[3]; ++a) {
                    for (int f=0; f<data_.shape()[4]; ++f) {
                        if (data_[t][w][g][a][f] > 1.0e10) {
                            //TODO: Replace with proper log message
                            std::cerr << "Too large value encountered in VFPPROD in ["
                                    << t << "," << w << "," << g << "," << a << "," << f << "]="
                                    << data_[t][w][g][a][f] << std::endl;
                        }
                    }
                }
            }
        }
    }
}



double VFPProperties::bhp(const double& flo, const double& thp, const double& wfr, const double& gfr, const double& alq) {
    //First, find the values to interpolate between
    auto flo_i = find_interp_data(flo, flo_data_);
    auto thp_i = find_interp_data(thp, thp_data_);
    auto wfr_i = find_interp_data(wfr, wfr_data_);
    auto gfr_i = find_interp_data(gfr, gfr_data_);
    auto alq_i = find_interp_data(alq, alq_data_);

    //Then perform the interpolation itself
    return interpolate(flo_i, thp_i, wfr_i, gfr_i, alq_i);
}

VFPProperties::ADB VFPProperties::bhp(const ADB& flo, const ADB& thp, const ADB& wfr, const ADB& gfr, const ADB& alq) {
    const ADB::V& f_v = flo.value();
    const ADB::V& t_v = thp.value();
    const ADB::V& w_v = wfr.value();
    const ADB::V& g_v = gfr.value();
    const ADB::V& a_v = alq.value();

    const int nw = f_v.size();

    //Compute the BHP for each well independently
    ADB::V bhp_vals;
    bhp_vals.resize(nw);
    for (int i=0; i<nw; ++i) {
        bhp_vals[i] = bhp(f_v[i], t_v[i], w_v[i], g_v[i], a_v[i]);
    }
    //Create an ADB constant value.
    return ADB::constant(bhp_vals);
}



VFPProperties::ADB VFPProperties::bhp(const Wells& wells, const ADB& qs, const ADB& thp, const ADB& alq) {
    const int np = wells.number_of_phases;
    const int nw = wells.number_of_wells;

    //Short-hands for water / oil / gas phases
    const ADB& w = subset(qs, Span(nw, 1, BlackoilPhases::Aqua*nw));
    const ADB& o = subset(qs, Span(nw, 1, BlackoilPhases::Liquid*nw));
    const ADB& g = subset(qs, Span(nw, 1, BlackoilPhases::Vapour*nw));

    ADB flo = getFlo(w, o, g);
    ADB wfr = getWFR(w, o, g);
    ADB gfr = getGFR(w, o, g);

    //TODO: What is this actually supposed to be used for?
    switch(alq_type_) {
        case ALQ_GRAT: // Lift as injection rate
            break;
        case ALQ_IGLR: // Injection gas-liquid ratio
            break;
        case ALQ_TGLR: // Total gas-liquid ratio
            break;
        case ALQ_PUMP: // Pump rating
            break;
        case ALQ_COMP: // Compressor power
            break;
        case ALQ_BEAN: // Choke diameter
            break;
        case ALQ_UNDEF: // Undefined
            break;
        case ALQ_INVALID: //Intentional fall-through
        default:
            //TODO: Log/throw
            std::cerr << "ERROR, ALQ_INVALID" << std::endl;
    }

    return bhp(flo, thp, wfr, gfr, alq);
}


VFPProperties::ADB VFPProperties::getFlo(const ADB& aqua, const ADB& liquid, const ADB& vapour) {
    switch (flo_type_) {
        case FLO_OIL:
            //Oil = liquid phase
            //TODO assert("qs has oil phase")
            return liquid;
        case FLO_LIQ:
            //Liquid = aqua + liquid phases
            return aqua + liquid;
        case FLO_GAS:
            //Gas = vapor phase
            return vapour;
        case FLO_INVALID: //Intentional fall-through
        default:
            //TODO: Log/throw
            std::cerr << "ERROR, FLO_INVALID" << std::endl;
            return ADB::null();
    }
}

VFPProperties::ADB VFPProperties::getWFR(const ADB& aqua, const ADB& liquid, const ADB& vapour) {
    switch(wfr_type_) {
        case WFR_WOR:
            //Water-oil ratio = water / oil
            return aqua / liquid;
        case WFR_WCT:
            //Water cut = water / (water + oil + gas)
            return aqua / (aqua + liquid + vapour);
        case WFR_WGR:
            //Water-gas ratio = water / gas
            return aqua / vapour;
        case WFR_INVALID: //Intentional fall-through
        default:
            //TODO: Log/throw
            std::cerr << "ERROR, WFR_INVALID" << std::endl;
            return ADB::null();
    }
}


VFPProperties::ADB VFPProperties::getGFR(const ADB& aqua, const ADB& liquid, const ADB& vapour) {
    switch(gfr_type_) {
        case GFR_GOR:
            // Gas-oil ratio = gas / oil
            return vapour / liquid;
        case GFR_GLR:
            // Gas-liquid ratio = gas / (oil + water)
            return vapour / (liquid + aqua);
        case GFR_OGR:
            // Oil-gas ratio = oil / gas
            return liquid / vapour;
        case GFR_INVALID: //Intentional fall-through
        default:
            //TODO: Log/throw
            std::cerr << "ERROR, GFR_INVALID" << std::endl;
            return ADB::null();
    }
}


VFPProperties::InterpData VFPProperties::find_interp_data(const double& value, const std::vector<double>& values) {
    InterpData retval;

    //First element greater than or equal to value
    //Start with the second element, so that floor_iter does not go out of range
    //Don't access out-of-range, therefore values.end()-1
    auto ceil_iter = std::lower_bound(values.begin()+1, values.end()-1, value);

    //Find last element smaller than range
    auto floor_iter = ceil_iter-1;

    //Find the indices
    const int a = floor_iter - values.begin();
    const int b = ceil_iter - values.begin();
    const int max_size = std::max(static_cast<int>(values.size()) - 1, 0);

    //Clamp indices to range of vector
    retval.ind_[0] = a;
    retval.ind_[1] = std::min(b, max_size);

    //Find interpolation ratio
    double dist = (*ceil_iter - *floor_iter);
    assert(dist >= 0.0);
    if (dist > 0.0) {
        //Possible source for floating point error here if value and floor are large,
        //but very close to each other
        retval.factor_ = (value-*floor_iter) / dist;
    }
    else {
        retval.factor_ = 1.0;
    }

    return retval;
}


#ifdef __GNUC__
#pragma GCC push_options
#pragma GCC optimize ("unroll-loops")
#endif

double VFPProperties::interpolate(const InterpData& flo_i, const InterpData& thp_i,
        const InterpData& wfr_i, const InterpData& gfr_i, const InterpData& alq_i) {
    double nn[2][2][2][2][2];

    //Pick out nearest neighbors (nn) to our evaluation point
    //This is not really required, but performance-wise it may pay off, since the 32-elements
    //we copy to (nn) will fit better in cache than the full original table for the
    //interpolation below.
    //The following ladder of for loops will presumably be unrolled by a reasonable compiler.
    for (int t=0; t<=1; ++t) {
        for (int w=0; w<=1; ++w) {
            for (int g=0; g<=1; ++g) {
                for (int a=0; a<=1; ++a) {
                    for (int f=0; f<=1; ++f) {
                        //Shorthands for indexing
                        const int ti = thp_i.ind_[t];
                        const int wi = wfr_i.ind_[w];
                        const int gi = gfr_i.ind_[g];
                        const int ai = alq_i.ind_[a];
                        const int fi = flo_i.ind_[f];

                        //Copy element
                        nn[t][w][g][a][f] = data_[ti][wi][gi][ai][fi];
                    }
                }
            }
        }
    }

    //Remove dimensions one by one
    // Example: going from 3D to 2D to 1D, we start by interpolating along
    // the z axis first, leaving a 2D problem. Then interpolating along the y
    // axis, leaving a 1D, problem, etc.
    double tf = flo_i.factor_;
    for (int t=0; t<=1; ++t) {
        for (int w=0; w<=1; ++w) {
            for (int g=0; g<=1; ++g) {
                for (int a=0; a<=1; ++a) {
                    nn[t][w][g][a][0] = (1.0-tf)*nn[t][w][g][a][0] + tf*nn[t][w][g][a][1];
                }
            }
        }
    }

    tf = alq_i.factor_;
    for (int t=0; t<=1; ++t) {
        for (int w=0; w<=1; ++w) {
            for (int g=0; g<=1; ++g) {
                nn[t][w][g][0][0] = (1.0-tf)*nn[t][w][g][0][0] + tf*nn[t][w][g][1][0];
            }
        }
    }

    tf = gfr_i.factor_;
    for (int t=0; t<=1; ++t) {
        for (int w=0; w<=1; ++w) {
            nn[t][w][0][0][0] = (1.0-tf)*nn[t][w][0][0][0] + tf*nn[t][w][1][0][0];
        }
    }

    tf = wfr_i.factor_;
    for (int t=0; t<=1; ++t) {
        nn[t][0][0][0][0] = (1.0-tf)*nn[t][0][0][0][0] + tf*nn[t][1][0][0][0];
    }

    tf = thp_i.factor_;
    return (1.0-tf)*nn[0][0][0][0][0] + tf*nn[1][0][0][0][0];
}

#ifdef __GNUC__
#pragma GCC pop_options //unroll loops
#endif


} //Namespace
