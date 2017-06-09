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

#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>

namespace Opm
{


    BlackoilModelParameters::BlackoilModelParameters()
    {
        // set default values
        reset();
    }




    BlackoilModelParameters::BlackoilModelParameters( const ParameterGroup& param )
    {
        // set default values
        reset();

        // overload with given parameters
        dp_max_rel_  = param.getDefault("dp_max_rel", dp_max_rel_);
        ds_max_      = param.getDefault("ds_max", ds_max_);
        dr_max_rel_  = param.getDefault("dr_max_rel", dr_max_rel_);
        dbhp_max_rel_=  param.getDefault("dbhp_max_rel", dbhp_max_rel_);
        dwell_fraction_max_ = param.getDefault("dwell_fraction_max", dwell_fraction_max_);
        max_residual_allowed_ = param.getDefault("max_residual_allowed", max_residual_allowed_);
        tolerance_mb_    = param.getDefault("tolerance_mb", tolerance_mb_);
        tolerance_cnv_   = param.getDefault("tolerance_cnv", tolerance_cnv_);
        tolerance_wells_ = param.getDefault("tolerance_wells", tolerance_wells_ );
        tolerance_well_control_ = param.getDefault("tolerance_well_control", tolerance_well_control_);
        maxSinglePrecisionTimeStep_ = unit::convert::from(
                param.getDefault("max_single_precision_days", unit::convert::to( maxSinglePrecisionTimeStep_, unit::day) ), unit::day );
        max_iter_ = param.getDefault("max_iter",10);
        solve_welleq_initially_ = param.getDefault("solve_welleq_initially",solve_welleq_initially_);
        update_equations_scaling_ = param.getDefault("update_equations_scaling", update_equations_scaling_);
        use_update_stabilization_ = param.getDefault("use_update_stabilization", use_update_stabilization_);
        deck_file_name_ = param.template get<std::string>("deck_filename");
    }




    void BlackoilModelParameters::reset()
    {
        // default values for the solver parameters
        dp_max_rel_      = 1.0;
        ds_max_          = 0.2;
        dr_max_rel_      = 1.0e9;
        dbhp_max_rel_    = 1.0;
        dwell_fraction_max_ = 0.2;
        max_residual_allowed_ = 1e7;
        tolerance_mb_    = 1.0e-5;
        tolerance_cnv_   = 1.0e-2;
        tolerance_wells_ = 1.0e-4;
        tolerance_well_control_ = 1.0e-7;
        maxSinglePrecisionTimeStep_ = unit::convert::from( 20.0, unit::day );
        solve_welleq_initially_ = true;
        update_equations_scaling_ = false;
        use_update_stabilization_ = true;
    }


} // namespace Opm
