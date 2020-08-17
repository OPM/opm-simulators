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
 * \brief Provides convenience routines to bring up the simulation at runtime.
 */
#ifndef EWOMS_STARTEBOS_HH
#define EWOMS_STARTEBOS_HH

#include "ebos.hh"
#include <opm/models/utils/start.hh>

#if HAVE_MPI
#include <mpi.h>
#endif

#if HAVE_ECL_INPUT
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/EclipsePRTLog.hpp>
#include <opm/common/OpmLog/LogUtil.hpp>
#endif

#include <opm/common/utility/String.hpp>

#include <opm/simulators/flow/Main.hpp>
//! \cond SKIP_THIS

namespace Opm {

//! \endcond

/*!
 * \ingroup Common
 *
 * \brief Wrapper around the main function that set up the OPM
 *        logging (.PRT, .DBG) for ebos.
 *
 * \tparam TypeTag  The type tag of the problem which needs to be solved
 *
 * \param argc The number of command line arguments
 * \param argv The array of the command line arguments
 */
template <class TypeTag>
static inline int startEbos(int argc, char **argv)
{

    int myRank = 0;
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argc, argv);
        myRank = Dune::Fem::MPIManager::rank();
#else
        myRank = Dune::MPIHelper::instance(argc, argv).rank();
#endif

        int paramStatus = setupParameters_<TypeTag>(argc, const_cast<const char**>(argv), /*registerParams=*/true);
        if (paramStatus == 1)
            return 1;
        if (paramStatus == 2)
            return 0;

    bool outputCout = false;
    if (myRank == 0)
        outputCout = EWOMS_GET_PARAM(TypeTag, bool, EnableTerminalOutput);

    std::string deckFilename = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);
    setupLogging(myRank,
                 deckFilename,
                 EWOMS_GET_PARAM(TypeTag, std::string, OutputDir),
                 EWOMS_GET_PARAM(TypeTag, std::string, OutputMode),
                 outputCout, "STDOUT_LOGGER");

    // Call the main function. Parameters are already registered
    // They should not be registered again
    return Opm::start<TypeTag>(argc, argv, /*registerParams=*/false);

}

} // namespace Opm

#endif
