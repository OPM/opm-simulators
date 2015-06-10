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

#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <boost/multi_array.hpp>

namespace Opm {

    class VFPProperties {
	public:
    	VFPProperties(DeckKeywordConstPtr table) {
            auto iter = table->begin();

            auto header = (*iter++);
    		table_num_ = header->getItem("TABLE")->getInt(0);
    		datum_depth_ = header->getItem("DATUM_DEPTH")->getRawDouble(0);

    		//Rate type
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

    		//Water fraction
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

    		//Gas fraction
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

    		//Artificial lift
    		/*
    		 * ALQ not implemented properly in parser?

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
    		*/

    		//Get actual rate / flow values
            const std::vector<double>& flo = (*iter++)->getItem("FLOW_VALUES")->getRawDoubleData();
            flo_data_.resize(flo.size());
            std::copy(flo.begin(), flo.end(), flo_data_.begin());

    		//Get actual tubing head pressure values
            const std::vector<double>& thp = (*iter++)->getItem("THP_VALUES")->getRawDoubleData();
            thp_data_.resize(thp.size());
            std::copy(thp.begin(), thp.end(), thp_data_.begin());

    		//Get actual water fraction values
            const std::vector<double>& wfr = (*iter++)->getItem("WFR_VALUES")->getRawDoubleData();
            wfr_data_.resize(wfr.size());
            std::copy(wfr.begin(), wfr.end(), wfr_data_.begin());

    		//Get actual gas fraction values
            const std::vector<double>& gfr = (*iter++)->getItem("GFR_VALUES")->getRawDoubleData();
            gfr_data_.resize(gfr.size());
            std::copy(gfr.begin(), gfr.end(), gfr_data_.begin());

    		//Get actual gas fraction values
            const std::vector<double>& alq = (*iter++)->getItem("ALQ_VALUES")->getRawDoubleData();
            alq_data_.resize(alq.size());
            std::copy(alq.begin(), alq.end(), alq_data_.begin());

            //Finally, read the actual table itself.
            unsigned int nt = thp_data_.size();
            unsigned int nw = wfr_data_.size();
            unsigned int ng = gfr_data_.size();
            unsigned int na = alq_data_.size();
            unsigned int nf = flo_data_.size();
            extents shape;
            shape[0] = nt;
            shape[1] = nw;
            shape[2] = ng;
            shape[3] = na;
            shape[4] = nf;
            data_.resize(shape);

			for (; iter!=table->end(); ++iter) {
				//Get indices (subtract 1 to get 0-based index)
            	unsigned int t = (*iter)->getItem("THP_INDEX")->getInt(0) - 1;
            	unsigned int w = (*iter)->getItem("WFR_INDEX")->getInt(0) - 1;
            	unsigned int g = (*iter)->getItem("GFR_INDEX")->getInt(0) - 1;
            	unsigned int a = (*iter)->getItem("ALQ_INDEX")->getInt(0) - 1;

            	//Rest of values (bottom hole pressure or tubing head temperature) have index of flo value
            	const std::vector<double>& bhp_tht = (*iter)->getItem("VALUES")->getRawDoubleData();
            	std::copy(bhp_tht.begin(), bhp_tht.end(), &data_[t][w][g][a][0]);

            	//Check for large values
            	for (unsigned int i = 0; i<bhp_tht.size(); ++i) {
            		if (bhp_tht[i] > 1.0e10) {
            			//TODO: Replace with proper log message
            			std::cerr << "Too large value encountered in VFPPROD in ["
            					<< t << "," << w << "," << g << "," << a << "]="
								<< bhp_tht[i] << std::endl;
            		}
            	}
            }
    	}

    	struct InterpData {
    		InterpData() : ind_({0}), factor_(0.0) {}
    		unsigned int ind_[2]; //[First element greater than or equal to value, Last element smaller than or equal to value]
    		double factor_; //Interpolation factor
    	};

    	InterpData find_interp_data(double value, const std::vector<double>& values) {
    		InterpData retval;

    		//First element greater than or equal to value
    		//Don't access out-of-range, therefore values.end()-1
    		auto ceil_iter = std::lower_bound(values.begin(), values.end()-1, value);

    		//Find last element smaller than or equal to range
    		auto floor_iter = ceil_iter;
    		if (*floor_iter == value) {
    			// floor_iter == ceil_iter == value
    		}
    		else if (floor_iter > values.begin()) {
    			// floor_iter <= value <= ceil_iter
    			--floor_iter;
    		}

    		//Now set these in the retval struct
    		retval.ind_[0] = floor_iter - values.begin();
    		retval.ind_[1] = ceil_iter - values.begin();

    		//Find interpolation ratio
    		double dist = (*ceil_iter - *floor_iter);
    		if (dist > 0) {
    			retval.factor_ = (value-*floor_iter) / dist;
    		}
    		else {
    			retval.factor_ = 1.0;
    		}

    		return retval;
    	}

    	double interpolate(const InterpData& flo_i, const InterpData& thp_i,
    			const InterpData& wfr_i, const InterpData& gfr_i, const InterpData& alq_i) {
    		extents shape({{2, 2, 2, 2, 2}});
        	array_type nn(shape);

    		//Pick out nearest neighbors (nn) to our evaluation point
        	//The following ladder of for loops will presumably be unrolled by a reasonable compiler.
        	//This is not really required, but performance-wise it may pay off, since the 32-elements
        	//we copy to (nn) will fit better in cache than the full original table for the
        	//interpolation below.
        	for (unsigned int t=0; t<=1; ++t) {
        		for (unsigned int w=0; w<=1; ++w) {
        			for (unsigned int g=0; g<=1; ++g) {
        				for (unsigned int a=0; a<=1; ++a) {
        					for (unsigned int f=0; f<=1; ++f) {
        						//Shorthands for indexing
        						unsigned int ti = thp_i.ind_[t];
        						unsigned int wi = wfr_i.ind_[w];
        						unsigned int gi = gfr_i.ind_[g];
        						unsigned int ai = alq_i.ind_[a];
        						unsigned int fi = flo_i.ind_[f];

        						//Copy element
        			        	nn[t][w][g][a][f] = data_[ti][wi][gi][ai][fi];
        					}
        				}
        			}
        		}
        	}

    		//Remove dimensions iteratively
    		// Example: going from 3D to 2D to 1D, we start by interpolating along
    		// the z axis first, leaving a 2D problem. Then interpolating along the y
    		// axis, leaving a 1D, problem, etc.
        	double tf = flo_i.factor_;
        	for (unsigned int t=0; t<=1; ++t) {
        		for (unsigned int w=0; w<=1; ++w) {
        			for (unsigned int g=0; g<=1; ++g) {
        				for (unsigned int a=0; a<=1; ++a) {
    			        	nn[t][w][g][a][0] = (1.0-tf)*nn[t][w][g][a][0] + tf*nn[t][w][g][a][1];
        				}
        			}
        		}
        	}

        	tf = alq_i.factor_;
        	for (unsigned int t=0; t<=1; ++t) {
        		for (unsigned int w=0; w<=1; ++w) {
        			for (unsigned int g=0; g<=1; ++g) {
						nn[t][w][g][0][0] = (1.0-tf)*nn[t][w][g][0][0] + tf*nn[t][w][g][1][0];
        			}
        		}
        	}

        	tf = gfr_i.factor_;
        	for (unsigned int t=0; t<=1; ++t) {
        		for (unsigned int w=0; w<=1; ++w) {
					nn[t][w][0][0][0] = (1.0-tf)*nn[t][w][0][0][0] + tf*nn[t][w][1][0][0];
        		}
        	}

        	tf = wfr_i.factor_;
        	for (unsigned int t=0; t<=1; ++t) {
				nn[t][0][0][0][0] = (1.0-tf)*nn[t][0][0][0][0] + tf*nn[t][1][0][0][0];
        	}

        	tf = thp_i.factor_;
        	return (1.0-tf)*nn[0][0][0][0][0] + tf*nn[1][0][0][0][0];
    	}

    	/**
    	 * Linear interpolation of bhp as a function of the input parameters
    	 * @param flo Production rate of oil, gas or liquid
    	 * @param thp Tubing head pressure
    	 * @param wfr Water-oil ratio, water cut, or water-gas ratio
    	 * @param gfr Gas-oil ratio, gas-liquid ratio, or oil-gas ratio
    	 * @param alq Artificial lift or other parameter
    	 *
    	 * @return The bottom hole pressure, interpolated linearly using
    	 * the above parameters from the values in the input table.
    	 */
    	double bhp(double flo, double thp, float wfr, float gfr, float alq) {
    		//First, find the floor value of the inputs
    		auto flo_i = find_interp_data(flo, flo_data_);
    		auto thp_i = find_interp_data(thp, thp_data_);
    		auto wfr_i = find_interp_data(wfr, wfr_data_);
    		auto gfr_i = find_interp_data(gfr, gfr_data_);
    		auto alq_i = find_interp_data(alq, alq_data_);

    		return interpolate(flo_i, thp_i, wfr_i, gfr_i, alq_i);
    	}

    	///Rate type
    	enum FLO_TYPE {
    		FLO_OIL,
			FLO_LIQ,
			FLO_GAS,
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

	private:
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
    	typedef boost::multi_array<double, 5> array_type;
    	typedef boost::array<array_type::index, 5> extents;
    	array_type data_;
    };

}

#endif /* OPM_AUTODIFF_VFPPROPERTIES_HPP_ */
