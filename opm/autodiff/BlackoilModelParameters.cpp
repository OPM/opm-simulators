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

namespace Opm
{


    BlackoilModelParameters::BlackoilModelParameters()
    {
        // set default values
        reset();
    }




    BlackoilModelParameters::BlackoilModelParameters( const parameter::ParameterGroup& param )
    {
        // set default values
        reset();

        // overload with given parameters
        dp_max_rel_  = param.getDefault("dp_max_rel", dp_max_rel_);
        ds_max_      = param.getDefault("ds_max", ds_max_);
        dr_max_rel_  = param.getDefault("dr_max_rel", dr_max_rel_);
        max_residual_allowed_ = param.getDefault("max_residual_allowed", max_residual_allowed_);
        tolerance_mb_    = param.getDefault("tolerance_mb", tolerance_mb_);
        tolerance_cnv_   = param.getDefault("tolerance_cnv", tolerance_cnv_);
        tolerance_wells_ = param.getDefault("tolerance_wells", tolerance_wells_ );
    }




    void BlackoilModelParameters::reset()
    {
        // default values for the solver parameters
        dp_max_rel_      = 1.0e9;
        ds_max_          = 0.2;
        dr_max_rel_      = 1.0e9;
        max_residual_allowed_ = 1e7;
        tolerance_mb_    = 1.0e-5;
        tolerance_cnv_   = 1.0e-2;
        tolerance_wells_ = 5.0e-1;
    }


} // namespace Opm
