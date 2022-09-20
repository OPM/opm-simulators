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
 *
 * \brief Contains the parameters required to extend the black-oil model by brine.
 */
#ifndef EWOMS_BLACK_OIL_BRINE_PARAMS_HH
#define EWOMS_BLACK_OIL_BRINE_PARAMS_HH

#include <opm/material/common/Tabulated1DFunction.hpp>

#include <vector>

namespace Opm {

//! \brief Struct holding the parameters for the BlackoilBrineModule class.
template<class Scalar>
struct BlackOilBrineParams {
    using TabulatedFunction = Tabulated1DFunction<Scalar>;

    std::vector<TabulatedFunction> bdensityTable_;
    std::vector<TabulatedFunction> permfactTable_;
    std::vector<Scalar> saltsolTable_;
    std::vector<Scalar> referencePressure_;
};

} // namespace Opm

#endif
