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
 * \copydoc Opm::ThreadManager
 */
#ifndef OPM_THREAD_MANAGER_HPP
#define OPM_THREAD_MANAGER_HPP

namespace Opm {

/*!
 * \brief Simplifies multi-threaded capabilities.
 */
class ThreadManager
{
public:
    enum {
#if defined(_OPENMP) || DOXYGEN
        //! Specify whether OpenMP is really available or not
        isFake = false
#else
        isFake = true
#endif
    };

    /*!
     * \brief Register all run-time parameters of the thread manager.
     */
    static void registerParameters();

    /*!
     * \brief Initialize number of threads used thread manager.
     *
     * \param queryCommandLineParameter if set to true we will query ThreadsPerProcess
     *        and if set (disregard the environment variable OPM_NUM_THREADS).
     *        If false we will assume that the number of OpenMP threads is already set
     *        outside of this function (e.g. by OPM_NUM_THREADS or in the simulator by
     *        the ThreadsPerProcess parameter).
     */
    static void init(bool queryCommandLineParameter = true);

    /*!
     * \brief Return the maximum number of threads of the current process.
     */
    static unsigned maxThreads()
    { return static_cast<unsigned>(numThreads_); }

    /*!
     * \brief Return the index of the current OpenMP thread
     */
    static unsigned threadId();

private:
    static int numThreads_;
};

} // namespace Opm

#endif // OPM_THREAD_MANAGER_HPP
