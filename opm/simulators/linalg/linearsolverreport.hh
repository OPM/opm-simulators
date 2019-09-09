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
 * \copydoc Opm::Linear::SolverReport
 */
#ifndef EWOMS_LINEAR_SOLVER_REPORT_HH
#define EWOMS_LINEAR_SOLVER_REPORT_HH

#include "convergencecriterion.hh"

#include <ewoms/common/timer.hh>
#include <ewoms/common/timerguard.hh>

namespace Opm {
namespace Linear {

/*!
 * \brief Collects summary information about the execution of the linear solver.
 */
class SolverReport
{
public:
    SolverReport()
    { reset(); }

    void reset()
    {
        timer_.halt();
        iterations_ = 0;
        converged_ = 0;
    }

    const Opm::Timer& timer() const
    { return timer_; }

    Opm::Timer& timer()
    { return timer_; }

    unsigned iterations() const
    { return iterations_; }

    void increment()
    { ++iterations_; }

    SolverReport& operator++()
    { ++iterations_; return *this; }

    bool converged() const
    { return converged_; }

    void setConverged(bool value)
    { converged_ = value; }

private:
    Opm::Timer timer_;
    unsigned iterations_;
    bool converged_;
};

}} // end namespace Linear, Ewoms

#endif
