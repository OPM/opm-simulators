/*
  Copyright 2014 Dr. Markus Blatt - HPC-Simulation-Software & Services

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

#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE OPM-ParallelIterativeSolverTest
#include <boost/test/unit_test.hpp>

// MPI header
#if HAVE_MPI
#include <dune/common/version.hh>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#else
#error "This file needs to compiled with MPI support!"
#endif
#include "DuneIstlTestHelpers.hpp"
#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <memory>
#include <cstdlib>
#include <string>

void run_test(const Opm::ParameterGroup& param)
{
    int N=100;
    int start, end, istart, iend;
    std::tie(start,istart,iend,end) = computeRegions(N);
    Opm::ParallelISTLInformation comm(MPI_COMM_WORLD);
    auto mat = create1DLaplacian(*comm.indexSet(), N, start, end, istart, iend);
    std::vector<double> x(end-start), b(end-start);
    createRandomVectors(comm, end-start, x, b, *mat);
    std::vector<double> exact(x);
    std::fill(x.begin(), x.end(), 0.0);
    Opm::LinearSolverFactory ls(param);
    boost::any anyComm(comm);
    ls.solve(b.size(), mat->data.size(), &(mat->rowStart[0]),
             &(mat->colIndex[0]), &(mat->data[0]), &(b[0]),
             &(x[0]), anyComm);
}

#ifdef HAVE_DUNE_ISTL
BOOST_AUTO_TEST_CASE(CGAMGTest)
{
    Opm::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("istl"));
    param.insertParameter(std::string("linsolver_type"), std::string("1"));
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    run_test(param);
}

BOOST_AUTO_TEST_CASE(CGILUTest)
{
    Opm::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("istl"));
    param.insertParameter(std::string("linsolver_type"), std::string("0"));
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    run_test(param);
}

BOOST_AUTO_TEST_CASE(BiCGILUTest)
{
    Opm::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("istl"));
    param.insertParameter(std::string("linsolver_type"), std::string("2"));
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    run_test(param);
}
#endif

