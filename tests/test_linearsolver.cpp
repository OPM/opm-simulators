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

#define BOOST_TEST_MODULE OPM-IterativeSolverTest
#include <boost/test/unit_test.hpp>

#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <dune/common/version.hh>
#include <memory>
#include <cstdlib>
#include <string>

struct MyMatrix
{
    MyMatrix(int rows, int nnz)
        : data(nnz, 0.0), rowStart(rows+1, -1),
          colIndex(nnz, -1)
    {}
    MyMatrix()
        : data(), rowStart(), colIndex()
    {}

    std::vector<double> data;
    std::vector<int> rowStart;
    std::vector<int> colIndex;
};

std::shared_ptr<MyMatrix> createLaplacian(int N)
{
    MyMatrix* mm=new MyMatrix(N*N, N*N*5);
    int nnz=0;
    mm->rowStart[0]=0;
    for(int row=0; row<N*N; row++)
    {

        int x=row%N;
        int y=row/N;
        if(y>0)
        {
            mm->colIndex[nnz]=row-N;
            mm->data[nnz++]=-1;
        }
        if(x>0)
        {
            mm->colIndex[nnz]=row-1;
            mm->data[nnz++]=-1;
        }
        mm->colIndex[nnz]=row;
        mm->data[nnz++]=4;
        if(x<N-1)
        {
            mm->colIndex[nnz]=row+1;
            mm->data[nnz++]=-1;
        }
        if(y<N-1)
        {
            mm->colIndex[nnz]=row+N;
            mm->data[nnz++]=-1;
        }
        mm->rowStart[row+1]=nnz;
    }
    mm->data.resize(nnz);
    mm->colIndex.resize(nnz);
    return std::shared_ptr<MyMatrix>(mm);
}

void createRandomVectors(int NN, std::vector<double>& x, std::vector<double>& b,
                         const MyMatrix& mat)
{
    x.resize(NN);
    for(auto entry=x.begin(), end =x.end(); entry!=end; ++entry)
        *entry=((double) (rand()%100))/10.0;

    b.resize(NN);
    std::fill(b.begin(), b.end(), 0.0);

    // Construct the right hand side as b=A*x
    for(std::size_t row=0; row<mat.rowStart.size()-1; ++row)
    {
        for(int i=mat.rowStart[row], end=mat.rowStart[row+1]; i!=end; ++i)
        {
            b[row]+= mat.data[i]*x[mat.colIndex[i]];
        }
    }
}

void run_test(const Opm::parameter::ParameterGroup& param)
{
    int N=4;
    auto mat = createLaplacian(N);
    std::vector<double> x(N*N), b(N*N);
    createRandomVectors(100*100, x, b, *mat);
    std::vector<double> exact(x);
    std::fill(x.begin(), x.end(), 0.0);
    Opm::LinearSolverFactory ls(param);
    ls.solve(N*N, mat->data.size(), &(mat->rowStart[0]),
             &(mat->colIndex[0]), &(mat->data[0]), &(b[0]),
             &(x[0]));
}


BOOST_AUTO_TEST_CASE(DefaultTest)
{
    Opm::parameter::ParameterGroup param;
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    param.insertParameter(std::string("linsolver_verbosity"), std::string("2"));
    run_test(param);
}

#ifdef HAVE_DUNE_ISTL
BOOST_AUTO_TEST_CASE(CGAMGTest)
{
    Opm::parameter::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("istl"));
    param.insertParameter(std::string("linsolver_type"), std::string("1"));
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    param.insertParameter(std::string("linsolver_verbosity"), std::string("2"));
    run_test(param);
}

BOOST_AUTO_TEST_CASE(CGILUTest)
{
    Opm::parameter::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("istl"));
    param.insertParameter(std::string("linsolver_type"), std::string("0"));
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    param.insertParameter(std::string("linsolver_verbosity"), std::string("2"));
    run_test(param);
}

BOOST_AUTO_TEST_CASE(BiCGILUTest)
{
    Opm::parameter::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("istl"));
    param.insertParameter(std::string("linsolver_type"), std::string("2"));
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    param.insertParameter(std::string("linsolver_verbosity"), std::string("2"));
    run_test(param);
}

#if defined(HAS_DUNE_FAST_AMG) || DUNE_VERSION_NEWER(DUNE_ISTL, 2, 3)
BOOST_AUTO_TEST_CASE(FastAMGTest)
{
    Opm::parameter::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("istl"));
    param.insertParameter(std::string("linsolver_type"), std::string("3"));
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    param.insertParameter(std::string("linsolver_verbosity"), std::string("2"));
    run_test(param);
}

BOOST_AUTO_TEST_CASE(KAMGTest)
{
    Opm::parameter::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("istl"));
    param.insertParameter(std::string("linsolver_type"), std::string("4"));
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    run_test(param);
}
#endif
#endif

#if HAVE_PETSC
BOOST_AUTO_TEST_CASE(PETScTest)
{
    Opm::parameter::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("petsc"));
    param.insertParameter(std::string("ksp_type"), std::string("cg"));
    param.insertParameter(std::string("pc_type"), std::string("jacobi"));
    param.insertParameter(std::string("ksp_rtol"), std::string("1e-10"));
    param.insertParameter(std::string("ksp_view"), std::string("0"));
    run_test(param);
}
#endif
