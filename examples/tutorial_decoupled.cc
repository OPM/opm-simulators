// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> /*@\label{tutorial-decoupled:include-begin}@*/
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/fractionalflow/variableclass2p.hh"
#include "dumux/material/fluids/water.hh"
#include "dumux/material/fluids/oil.hh"
#include "tutorial_soilproperties_decoupled.hh"
#include "dumux/material/twophaserelations.hh"
#include "tutorialproblem_decoupled.hh"
#include "dumux/diffusion/fv/fvtotalvelocity2p.hh"
#include "dumux/transport/fv/fvsaturationwetting2p.hh"
#include "dumux/fractionalflow/impes/impes.hh"
#include "dumux/timedisc/timeloop.hh" /*@\label{tutorial-decoupled:include-end}@*/


int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=2; /*@\label{tutorial-decoupled:dim}@*/

        // create a grid object
        typedef double Scalar; /*@\label{tutorial-decoupled:grid-begin}@*/
        typedef Dune::SGrid<dim,dim> Grid;
        typedef Grid::LevelGridView GridView;
        typedef Dune::FieldVector<Grid::ctype,dim> FieldVector;
        Dune::FieldVector<int,dim> N(10); N[0] = 30;
        FieldVector L(0);
        FieldVector H(300); H[0] = 600;
        Grid grid(N,L,H);
        GridView gridView(grid.levelView(0));/*@\label{tutorial-decoupled:grid-end}@*/


        // define fluid and solid properties and constitutive relationships
        Dune::Water wettingfluid; /*@\label{tutorial-decoupled:water}@*/
        Dune::Oil nonwettingfluid; /*@\label{tutorial-decoupled:oil}@*/
        Dune::TutorialSoil<Grid, Scalar> soil; /*@\label{tutorial-decoupled:soil}@*/
        Dune::TwoPhaseRelations<Grid, Scalar> materialLaw(soil, wettingfluid, nonwettingfluid);/*@\label{tutorial-decoupled:twophaserelations}@*/

        // create object containing the variables
        typedef Dune::VariableClass<GridView, Scalar> VariableClass;
        VariableClass variables(gridView);

        // create object including the problem definition
        typedef Dune::TutorialProblemDecoupled<GridView, Scalar, VariableClass> Problem;
        Problem problem(variables, wettingfluid, nonwettingfluid, soil, materialLaw,L, H); /*@\label{tutorial-decoupled:problem}@*/

        // create object including the discretisation of the pressure equation
        typedef Dune::FVTotalVelocity2P<GridView, Scalar, VariableClass, Problem> Diffusion;
        Diffusion diffusion(gridView, problem, "pw"); /*@\label{tutorial-decoupled:diffusion}@*/

        // create object including the space discretisation of the saturation equation
        typedef Dune::FVSaturationWetting2P<GridView, Scalar, VariableClass, Problem> Transport;
        Transport transport(gridView, problem, "vt"); /*@\label{tutorial-decoupled:transport}@*/

        // some parameters used in the IMPES-object
        int iterFlag = 2;
        int nIter = 30;
        double maxDefect = 1e-5;

        // create object including the IMPES (IMplicit Pressure Explicit Saturation) algorithm
        typedef Dune::IMPES<GridView, Diffusion, Transport, VariableClass> IMPES;
        IMPES impes(diffusion, transport, iterFlag, nIter, maxDefect); /*@\label{tutorial-decoupled:impes}@*/

        // some parameters needed for the TimeLoop-object
        double tStart = 0; // start simulation at t = tStart
        double tEnd = 1e8; // stop simulation at t = tEnd
        const char* fileName = "tutorial_decoupled"; // name of the output files
        int modulo = 1; // define time step interval in which output files are generated
        double cFLFactor = 0.9; // security factor for the Courant-Friedrichs-Lewy-Criterion

        // create TimeLoop-object
        Dune::TimeLoop<Grid, IMPES> timeloop(tStart, tEnd, fileName, modulo, cFLFactor); /*@\label{tutorial-decoupled:timeloop}@*/

        Dune::Timer timer;
        timer.reset();

        // start simulation
        timeloop.execute(impes); /*@\label{tutorial-decoupled:execute}@*/

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
    }
}
