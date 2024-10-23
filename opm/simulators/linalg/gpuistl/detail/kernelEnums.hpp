/*
  Copyright 2024 SINTEF AS

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

#ifndef OPM_GPUISTL_PRECONDITIONER_STORAGE_OPTION_HPP
#define OPM_GPUISTL_PRECONDITIONER_STORAGE_OPTION_HPP

/*
    This file is here to organize a growing amount of different mixed precision options for the preconditioners.
*/

namespace Opm::gpuistl {
    enum MixedPrecisionScheme {
        DEFAULT,
        STORE_ENTIRE_FACTORIZATION_AS_FLOAT,
        STORE_ONLY_FACTORIZED_DIAGONAL_AS_DOUBLE
    };

    inline bool isValidMixedPrecisionScheme(int scheme) {
        switch (scheme) {
            case DEFAULT:
            case STORE_ENTIRE_FACTORIZATION_AS_FLOAT:
            case STORE_ONLY_FACTORIZED_DIAGONAL_AS_DOUBLE:
                return true;
            default:
                return false;
        }
    }
}

#endif // OPM_GPUISTL_PRECONDITIONER_STORAGE_OPTION_HPP
