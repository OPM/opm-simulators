// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 20010 by Markus Wolff                                     *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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

#include <iostream> /*@\label{tutorial-decoupled:include-end}@*/

////////////////////////////////////////////
// function to check the input parameters
////////////////////////////////////////////
void usage(const char *progname)
{
    std::cout << "usage: "<<progname<<" [--restart restartTime] tEnd\n";
    exit(1);
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    try {
        typedef TTAG(TutorialProblemDecoupled) TypeTag; /*@\label{tutorial-decoupled:set-type-tag}@*/
        typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;    /*@\label{tutorial-decoupled:retrieve-types-begin}@*/
        typedef GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
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
        double startTime = 0.;
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;
            // use restart time as start time
            std::istringstream(argv[argPos++]) >> startTime;
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
        Grid *gridPtr = GET_PROP(TypeTag, Grid)::create(); /*@\label{tutorial-decoupled:create-grid}@*/

        // create time manager responsible for global simulation control
        TimeManager timeManager;

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        Problem problem(timeManager, gridPtr->leafView()); /*@\label{tutorial-decoupled:instantiate-problem}@*/

        // define simulation parameters
        timeManager.init(problem, startTime, dt, tEnd, restart); /*@\label{tutorial-decoupled:initTimeManager}@*/

        // run the simulation
        timeManager.run();    /*@\label{tutorial-decoupled:execute}@*/
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
