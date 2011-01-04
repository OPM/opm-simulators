// $Id$
/*****************************************************************************
 *   Copyright (C) 20010 by Markus Wolff                                     *
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
 * \brief tutorial for the sequential two-phase model
 */
#include "config.h" /*@\label{tutorial-decoupled:include-begin}@*/

#include "tutorialproblem_decoupled.hh" /*@\label{tutorial-decoupled:include-problem-header}@*/

#include <dune/grid/common/gridinfo.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp> /*@\label{tutorial-decoupled:include-end}@*/


////////////////////////////////////////////
// function to check the input parameters
////////////////////////////////////////////
void usage(const char *progname)
{
    std::cout << boost::format("usage: %s [--restart restartTime] tEnd\n")%progname;
    exit(1);
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    try {
        typedef TTAG(TutorialProblemDecoupled) TypeTag; /*@\label{tutorial-decoupled:set-type-tag}@*/
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;    /*@\label{tutorial-decoupled:retrieve-types-begin}@*/
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition; /*@\label{tutorial-decoupled:retrieve-types-end}@*/

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);  /*@\label{tutorial-decoupled:init-mpi}@*/

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc < 2)   /*@\label{tutorial-decoupled:parse-args-begin}@*/
            usage(argv[0]);

        // deal with the restart stuff
        int argPos = 1;
        bool restart = false;
        double restartTime = 0;
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;

            std::istringstream(argv[argPos++]) >> restartTime;
        }
        // output in case of wrong numbers of input parameters
        if (argc - argPos != 1) {
            usage(argv[0]);
        }

        // read the initial time step and the end time
        double tEnd, dt;
        std::istringstream(argv[argPos++]) >> tEnd;
        dt = tEnd;  /*@\label{tutorial-decoupled:parse-args-end}@*/

        // create the grid
        Grid *gridPtr = GET_PROP(TypeTag, PTAG(Grid))::create(); /*@\label{tutorial-decoupled:create-grid}@*/


        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////

        Problem problem(gridPtr->leafView()); /*@\label{tutorial-decoupled:instantiate-problem}@*/

        // load restart file if necessary
        if (restart)    /*@\label{tutorial-decoupled:mainRestart}@*/
            problem.deserialize(restartTime);

        // define simulation parameters
        problem.timeManager().init(problem, 0, dt, tEnd, !restart); /*@\label{tutorial-decoupled:initTimeManager}@*/
        // run the simulation
        problem.timeManager().run();    /*@\label{tutorial-decoupled:execute}@*/
        return 0;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }

    return 3;
}
