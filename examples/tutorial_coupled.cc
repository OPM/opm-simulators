// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Main file of the tutorial for a fully coupled twophase box model.
 */
#include "config.h" /*@\label{tutorial-coupled:include-begin}@*/
#include "tutorialproblem_coupled.hh"  /*@\label{tutorial-coupled:include-problem-header}@*/

#include <dune/common/mpihelper.hh>
#include <iostream> /*@\label{tutorial-coupled:include-end}@*/

void usage(const char *progname)
{
    std::cout << "usage: " << progname << " [--restart restartTime] tEnd dt\n";
    exit(1);
};

int main(int argc, char** argv)
{
    try {
        typedef TTAG(TutorialProblemCoupled) TypeTag; /*@\label{tutorial-coupled:set-type-tag}@*/
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid; /*@\label{tutorial-coupled:retrieve-types-begin}@*/
        typedef GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem; /*@\label{tutorial-coupled:retrieve-types-end}@*/

        // Initialize MPI
        Dune::MPIHelper::instance(argc, argv); /*@\label{tutorial-coupled:init-mpi}@*/

        // parse the command line arguments
        if (argc < 3) /*@\label{tutorial-coupled:parse-args-begin}@*/
            usage(argv[0]);

        // parse restart time if restart is requested
        int argPos = 1;
        bool restart = false;
        double restartTime = 0;
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;

            std::istringstream(argv[argPos++]) >> restartTime;
        }

        // read the initial time step and the end time
        if (argc - argPos != 2)
            usage(argv[0]);

        double tEnd, dt; /*@\label{tutorial-coupled:parse-tEn-and-dt}@*/
        std::istringstream(argv[argPos++]) >> tEnd;
        std::istringstream(argv[argPos++]) >> dt; /*@\label{tutorial-coupled:parse-args-end}@*/

        // create the grid
        Grid *gridPtr = GET_PROP(TypeTag, PTAG(Grid))::create(); /*@\label{tutorial-coupled:create-grid}@*/

        // create time manager responsible for global simulation control
        TimeManager timeManager;

        // instantiate the problem on the leaf grid
        Problem problem(timeManager, gridPtr->leafView()); /*@\label{tutorial-coupled:instantiate-problem}@*/
        timeManager.init(problem, 0, dt, tEnd, !restart);
        // load some previously saved state from disk
        if (restart)
            problem.restart(restartTime); /*@\label{tutorial-coupled:begin-restart}@*/
        // run the simulation
        timeManager.run(); /*@\label{tutorial-coupled:execute}@*/
        return 0;
    }
    catch (Dune::Exception &e) { /*@\label{tutorial-coupled:catch-dune-exceptions}@*/
        // Catch exceptions thrown somewhere in DUNE
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {  /*@\label{tutorial-coupled:catch-other-exceptions}@*/
        // Catch exceptions thrown elsewhere
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }

    return 3;
}
