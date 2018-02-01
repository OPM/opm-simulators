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
 * \brief This test ensures that the API of the black-oil fluid state conforms to the
 *        fluid state specification
 */
#include "config.h"

#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/material/checkFluidSystem.hpp>

#include <dune/common/parallel/mpihelper.hh>

int main()
{
    {
        typedef double Scalar;
        typedef double Evaluation;
        typedef typename Opm::FluidSystems::BlackOil<Scalar> FluidSystem;
        typedef Opm::BlackOilFluidState<Scalar, FluidSystem> FluidState;

        FluidState fs;
        checkFluidState<Evaluation>(fs);
    }

    {
        typedef float Scalar;
        typedef float Evaluation;
        typedef typename Opm::FluidSystems::BlackOil<Scalar> FluidSystem;
        typedef Opm::BlackOilFluidState<Scalar, FluidSystem> FluidState;

        FluidState fs;
        checkFluidState<Evaluation>(fs);
    }

    {
        typedef float Scalar;
        typedef Opm::DenseAd::Evaluation<Scalar, 2> Evaluation;
        typedef typename Opm::FluidSystems::BlackOil<Scalar> FluidSystem;
        typedef Opm::BlackOilFluidState<Evaluation, FluidSystem> FluidState;

        FluidState fs;
        checkFluidState<Evaluation>(fs);
    }

    return 0;
}
