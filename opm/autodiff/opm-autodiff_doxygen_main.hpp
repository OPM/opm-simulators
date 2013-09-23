/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_OPM-AUTODIFF_DOXYGEN_MAIN_HEADER_INCLUDED
#define OPM_OPM-AUTODIFF_DOXYGEN_MAIN_HEADER_INCLUDED


/** \mainpage Documentation for the opm-autodiff library.

<h3>Automatic differentiation</h3>

This library implements automatic differentiation for vector data with
multiple blocks of sparse jacobians. This is contained in the class
Opm::AutoDiffBlock. Also available is Opm::AutoDiff, a much simpler
single-value single-derivative AD class.

There are also some helper classes and functions that are intended to
aid in the development of solvers and simulators with AD, these
include Opm::HelperOps, Opm::UpwindSelector, Opm::subset,
Opm::superset, Opm::Selector, Opm::collapseJacs, Opm::vertcat,
Opm::Span and Opm::sign.

<h3>Solvers and simulators</h3>

There are some solvers and simulators in opm-autodiff. They should all
be considered experimental prototypes at this point. Notable simulator
prototypes include
- examples/sim_fibo_ad.cpp, a fully implicit black-oil simulator.
- examples/sim_2p_incomp_ad.cpp, a sequential incompressible 2-phase simulator.

*/

#endif // OPM_OPM-AUTODIFF_DOXYGEN_MAIN_HEADER_INCLUDED
