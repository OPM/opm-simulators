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
#include <mpi.h>
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/communicator.hh>
#include <dune/common/parallel/remoteindices.hh>
#include <dune/common/mpicollectivecommunication.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#else
#error "This file needs to compiled with MPI support!"
#endif
#include <opm/core/linalg/LinearSolverFactory.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <memory>
#include <cstdlib>
#include <string>

struct MPIFixture {
    MPIFixture()
    {    
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        MPI_Init(&m_argc, &m_argv);
    }
    ~MPIFixture()
    {
        MPI_Finalize();
    }
};

BOOST_GLOBAL_FIXTURE(MPIFixture);

struct MyMatrix
{
    MyMatrix(std::size_t rows, std::size_t nnz)
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

typedef int LocalId;
typedef int GlobalId;
typedef Dune::OwnerOverlapCopyCommunication<GlobalId,LocalId> Communication;
typedef Dune::OwnerOverlapCopyAttributeSet GridAttributes;
typedef GridAttributes::AttributeSet GridFlag;
typedef Dune::ParallelLocalIndex<GridFlag> LocalIndex;

/// \brief Sets up a paralle Laplacian.
///
/// The process stores the unknowns with indices in the range [start, end).
/// As we use an overlapping domain decomposition, the process owns the indices
/// in the range [istart, iend]. If we would only used the indices in this range then
/// they form a partitioning of the whole index set.
/// \tparam I The type of the parallel index set (for convenience)
/// \param indexset The parallel index set for marking owner and copy region.
/// \param N The global number of unknowns of the system.
/// \param start The first index stored on this process
/// \param end One past the last index stored on this process
/// \param istart The first index that the process owns.
/// \param iend One past the last index the process owns.
template<class I>
std::shared_ptr<MyMatrix> create1DLaplacian(I& indexset, int N, int start, int end,
                                            int istart, int iend)
{
    indexset.beginResize();
    MyMatrix* mm=new MyMatrix(end-start, (end-start)*3);
    int nnz=0;
    mm->rowStart[0]=0;
    assert(start==0||start<istart);
    assert(end==N||iend<end);
    
    for(int row=start, localRow=0; row<end; row++, localRow++)
    {
        if(row<istart || row>=iend)
        {
            // We are in the overlap region of the grid
            // therefore we setup the system such that
            // right hand side will equal the left hand side
            // of the linear system.
            mm->colIndex[nnz]=localRow;
            mm->data[nnz++]=1.0;
            indexset.add(row, LocalIndex(localRow, GridAttributes::copy, true));
                         mm->rowStart[localRow+1]=nnz;
            continue;
        }
        
        double dval=0;
        if(row>0)
        {
            mm->colIndex[nnz]=localRow-1;
            mm->data[nnz++]=-1;
            dval+=1;
        }
        mm->colIndex[nnz]=localRow;
        mm->data[nnz++]=2;//dval+(row<N-1);
        if(row<N-1)
        {
            mm->colIndex[nnz]=localRow+1;
            mm->data[nnz++]=-1;
            dval+=1;
        }
        mm->rowStart[localRow+1]=nnz;
        indexset.add(row, LocalIndex(localRow, GridAttributes::owner, true));
    }
    mm->data.resize(nnz);
    mm->colIndex.resize(nnz);
    indexset.endResize();
    return std::shared_ptr<MyMatrix>(mm);
}

template<class O>
void createRandomVectors(O& pinfo, int NN, std::vector<double>& x, std::vector<double>& b,
                         const MyMatrix& mat)
{
    x.resize(NN);
    for(auto entry=x.begin(), end =x.end(); entry!=end; ++entry)
        *entry=((double) (rand()%100))/10.0;

    pinfo.copyOwnerToAll(x,x);
    
    b.resize(NN);

    // Construct the right hand side as b=A*x
    std::fill(b.begin(), b.end(), 0.0);
    for(std::size_t row=0; row<mat.rowStart.size()-1; ++row)
    {
        for(int i=mat.rowStart[row], end=mat.rowStart[row+1]; i!=end; ++i)
        {
            b[row]+= mat.data[i]*x[mat.colIndex[i]];
        }
    }
    pinfo.copyOwnerToAll(b,b);
}

void run_test(const Opm::parameter::ParameterGroup& param)
{
    int N=100;
    int procs, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);
    int n = N/procs; // number of unknowns per process
    int bigger = N%procs; // number of process with n+1 unknows


    int start, end, istart, iend;
    // Compute owner region
    if(rank<bigger) {
        start = rank*(n+1);
        end   = start+(n+1);
    }else{
        start = bigger*(n+1) + (rank-bigger) * n;
        end   = start+n;
    }
    // Compute owner region
    if(rank<bigger) {
        istart = rank*(n+1);
        iend   = start+(n+1);
    }else{
        istart = bigger*(n+1) + (rank-bigger) * n;
        iend   = start+n;
    }

    // Compute overlap region
    if(istart>0)
        start = istart - 1;
    else
        start = istart;
    
    if(iend<N)
        end = iend + 1;
    else
        end = iend;

    Opm::ParallelISTLInformation comm(MPI_COMM_WORLD);
    auto mat = create1DLaplacian(*comm.indexSet(), N, start, end, istart, iend);
    std::vector<double> x(end-start), b(end-start);
    createRandomVectors(comm, end-start, x, b, *mat);
    std::vector<double> exact(x);
    std::fill(x.begin(), x.end(), 0.0);
    Opm::LinearSolverFactory ls(param);
    boost::any anyComm(comm);
    ls.solve(N, mat->data.size(), &(mat->rowStart[0]),
             &(mat->colIndex[0]), &(mat->data[0]), &(b[0]),
             &(x[0]), anyComm);
}

#ifdef HAVE_DUNE_ISTL
BOOST_AUTO_TEST_CASE(CGAMGTest)
{
    Opm::parameter::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("istl"));
    param.insertParameter(std::string("linsolver_type"), std::string("1"));
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    run_test(param);
}

BOOST_AUTO_TEST_CASE(CGILUTest)
{
    Opm::parameter::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("istl"));
    param.insertParameter(std::string("linsolver_type"), std::string("0"));
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    run_test(param);
}

BOOST_AUTO_TEST_CASE(BiCGILUTest)
{
    Opm::parameter::ParameterGroup param;
    param.insertParameter(std::string("linsolver"), std::string("istl"));
    param.insertParameter(std::string("linsolver_type"), std::string("2"));
    param.insertParameter(std::string("linsolver_max_iterations"), std::string("200"));
    run_test(param);
}
#endif

