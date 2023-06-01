/*
  Copyright 2022-2023 SINTEF AS

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
#ifndef OPM_CUISTL_PRECONDITIONERHOLDER_HPP
#define OPM_CUISTL_PRECONDITIONERHOLDER_HPP
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>


namespace Opm::cuistl
{
//! \brief Common interface for adapters that hold preconditioners.
//!
//! \note The goal is that this class will be made useless after further
//! restructuring of the solver interfaces (FlexibleSolver and friends), but
//! for this is needed as an intermediate layer. See specifically SolverAdapter.hpp
//! for how this is used.
template <class X, class Y>
class PreconditionerHolder
{
public:
    /**
     * @brief getUnderlyingPreconditioner gets the underlying preconditioner (preconditioner being held)
     */
    virtual std::shared_ptr<Dune::PreconditionerWithUpdate<X, Y>> getUnderlyingPreconditioner() = 0;
};
} // end namespace Opm::cuistl

#endif
