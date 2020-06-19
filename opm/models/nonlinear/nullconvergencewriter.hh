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
 * \copydoc Opm::NullConvergenceWriter
 */
#ifndef EWOMS_NULL_CONVERGENCE_WRITER_HH
#define EWOMS_NULL_CONVERGENCE_WRITER_HH

#include <opm/models/utils/propertysystem.hh>

#include <opm/material/common/Unused.hpp>

namespace Opm::Properties {

template <class TypeTag, class MyTypeTag>
struct NewtonMethod;

} // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup Newton
 *
 * \brief A convergence writer for the Newton method which does nothing
 */
template <class TypeTag>
class NullConvergenceWriter
{
    using NewtonMethod = GetPropType<TypeTag, Properties::NewtonMethod>;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;

public:
    NullConvergenceWriter(NewtonMethod& method  OPM_UNUSED)
    {}

    /*!
     * \brief Called by the Newton method before the actual algorithm
     *        is started for any given timestep.
     */
    void beginTimeStep()
    {}

    /*!
     * \brief Called by the Newton method before an iteration of the
     *        Newton algorithm is started.
     */
    void beginIteration()
    {}

    /*!
     * \brief Write the Newton update to disk.
     *
     * Called after the linear solution is found for an iteration.
     *
     * \param uLastIter The solution vector of the previous iteration.
     * \param deltaU The negative difference between the solution
     *        vectors of the previous and the current iteration.
     */
    void writeFields(const SolutionVector& uLastIter  OPM_UNUSED,
                     const GlobalEqVector& deltaU  OPM_UNUSED)
    {}

    /*!
     * \brief Called by the Newton method after an iteration of the
     *        Newton algorithm has been completed.
     */
    void endIteration()
    {}

    /*!
     * \brief Called by the Newton method after Newton algorithm
     *        has been completed for any given timestep.
     *
     * This method is called regardless of whether the Newton method
     * converged or not.
     */
    void endTimeStep()
    {}
};

} // namespace Opm

#endif
