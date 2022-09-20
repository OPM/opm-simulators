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
 * \brief Contains the parameters required to extend the black-oil model by solvent component.
 *  For details, refer:
 *  [*] T.H. Sandve, O. Sævareid and I. Aavatsmark: “Improved Extended Blackoil Formulation
 *  for CO2 EOR Simulations.” in ECMOR XVII – The 17th European Conference on the
 *  Mathematics of Oil Recovery,  September 2020.
 */
#ifndef EWOMS_BLACK_OIL_EXTBO_PARAMS_HH
#define EWOMS_BLACK_OIL_EXTBO_PARAMS_HH

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>

#include <vector>

namespace Opm {

//! \brief Struct holding the parameters for the BlackoilExtboModule class.
template<class Scalar>
struct BlackOilExtboParams {
    using TabulatedFunction = Tabulated1DFunction<Scalar>;
    using Tabulated2DFunction = UniformXTabulated2DFunction<Scalar>;

    std::vector<Tabulated2DFunction> X_;
    std::vector<Tabulated2DFunction> Y_;
    std::vector<Tabulated2DFunction> PBUB_RS_;
    std::vector<Tabulated2DFunction> PBUB_RV_;
    std::vector<Tabulated2DFunction> VISCO_;
    std::vector<Tabulated2DFunction> VISCG_;
    std::vector<Tabulated2DFunction> BO_;
    std::vector<Tabulated2DFunction> BG_;
    std::vector<Tabulated2DFunction> RS_;
    std::vector<Tabulated2DFunction> RV_;

    std::vector<Scalar> zReferenceDensity_;

    std::vector<Scalar> zLim_;
    std::vector<TabulatedFunction> oilCmp_;
    std::vector<TabulatedFunction> gasCmp_;
};

} // namespace Opm

#endif
