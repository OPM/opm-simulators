/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 Statoil AS

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

#include <boost/test/unit_test.hpp>

#include <opm/autodiff/SendReceiveCommunicator.hpp>
//#include "2dslicedinterface.hpp"
// MPI header
#if HAVE_MPI
#include <dune/common/parallel/interface.hh>
#include <mpi.h>

Dune::Interface::InformationMap setup2DSlicedInterface(int size, int rank,
                                                       int N=5)
{
        // We assume a partitioning on regular structured 2D grid distributed
        // in 1D (along the x-axis) with N cells along the y-axis.
        // Each processes partition is N cells wide.
        // each process additionally stores one overlap (also called
        // ghost or halo layer). Each process (except for the last) needs to
        // send one layer of cells
        // next to the right hand side boundary of its partition to the next
        // high process.
        typedef Dune::Interface::InformationMap InterfaceMap;
        InterfaceMap interface;

        int x_start = rank *N, x_end = (rank+1) * N;
        int y_global_width = N;
        int x_overlap_start = x_start, x_overlap_end = x_end;
        if ( rank > 0 ){
            --x_overlap_start;
        }
        if ( rank < size - 1){
            ++x_overlap_end;
        }
        int x_overlap_width = x_overlap_end - x_overlap_start;
        auto& right_interface = interface[rank+1].first;
        auto& left_interface = interface[rank-1].second;
        if ( rank < size - 1)
        {
            // Setup the send interface.
            right_interface.reserve(y_global_width);
            // send the indices left to the right partition boundary.
            for(int i=0; i < y_global_width; i++)
            {
                right_interface.add((i + 1) * x_overlap_width - 2);
            }
        }
        if( rank > 0)
        {
            // Setup the receive interface.
            left_interface.reserve(y_global_width);
            // receive in the left overlap cells.
            for(int i=0; i < y_global_width; i++)
            {
                left_interface.add(i * x_overlap_width);
            }
        }
        return interface;
}

class MPIError {
public:
    /** @brief Constructor. */
    MPIError(std::string s, int e) : errorstring(s), errorcode(e){}
    /** @brief The error string. */
    std::string errorstring;
    /** @brief The mpi error code. */
    int errorcode;
};

void MPI_err_handler(MPI_Comm *, int *err_code, ...){
    char *err_string=new char[MPI_MAX_ERROR_STRING];
    int err_length;
    MPI_Error_string(*err_code, err_string, &err_length);
    std::string s(err_string, err_length);
    std::cerr << "An MPI Error ocurred:"<<std::endl<<s<<std::endl;
    delete[] err_string;
    throw MPIError(s, *err_code);
}

struct GraphConnectionHandle
{
    typedef std::size_t DataType;

    std::size_t x_overlap_start_;
    std::size_t x_overlap_width_;
    std::size_t x_global_width_;
    std::size_t y_global_width_;
    int rank;

    GraphConnectionHandle(int size, int rank, int N=5)
        : y_global_width_(N)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int x_start = rank *N, x_end = (rank+1) * N;
        x_overlap_start_ = x_start;
        int x_overlap_end = x_end;
        if ( rank > 0 ) --x_overlap_start_;
        if ( rank < size - 1) ++x_overlap_end;
        x_overlap_width_ = x_overlap_end - x_overlap_start_;
        x_global_width_ = size * N;
    }

    bool fixedsize() const
    {
        return false;
    }
    std::size_t size(std::size_t i) const
    {
        std::size_t size = 1;
        std::size_t local_x = i % x_overlap_width_;
        std::size_t local_y = i / x_overlap_width_;
        if( local_x > 0 )
            ++size;
        if( local_x < x_overlap_width_ -1 )
            ++size;
        if( local_y > 0 )
            ++size;
        if( local_y < y_global_width_ - 1 )
            ++size;
        return size;
    }
    template<class B>
    void gather(B& buf, std::size_t i)
    {
        std::vector<std::size_t> neighbours=computeGlobalNeighbours(i);
        for( auto& neighbour : neighbours)
        {
            buf.write(neighbour);
        }
    }

    std::vector<size_t> computeGlobalNeighbours(std::size_t i)
    {
        std::vector<std::size_t> neighbours;
        neighbours.reserve(5);
        std::size_t local_x = i % x_overlap_width_;
        std::size_t local_y = i / x_overlap_width_;
        std::size_t global  = local_y * x_global_width_ +
            x_overlap_start_ + local_x;
        // Note that neighbours will be sorted!
        if( local_y > 0 )
            neighbours.push_back(global-x_global_width_);
        if( local_x > 0 )
            neighbours.push_back(global - 1);
        neighbours.push_back(global );
        if( local_x < x_overlap_width_ -1 )
            neighbours.push_back(global + 1);
        if( local_y < y_global_width_ - 1 )
            neighbours.push_back(global+x_global_width_);
        return neighbours;
    }

    template<class B>
    void scatter(B& buf, std::size_t i, std::size_t n)
    {
        std::size_t local_x = i % x_overlap_width_;
        std::size_t local_y = i / x_overlap_width_;
        std::size_t global  = local_y * x_global_width_ +
            x_overlap_start_ + local_x;

        if ( n > 0)
        {
            std::vector<std::size_t> message(n);

            for( auto& item : message )
            {
                buf.read(item);
            }

            std::vector<std::size_t> neighbours = computeGlobalNeighbours(i);
            // Test that each neighor is in the received values.
            BOOST_CHECK(neighbours.size() <= message.size());
            // Check that we receive each global neighbour.
            // Note that message and neighbours is sorted.
            auto start= message.begin();
            for( auto& neighbour : neighbours)
            {
                auto found = std::lower_bound(start, message.end(), neighbour);
                BOOST_CHECK_MESSAGE( found != message.end(),
                                     "Did not receive expected connection "
                                     << " to global index "<<neighbour<<
                                     " for global index "<<global);
                if( found != message.end() )
                {
                    start=found;
                }
            }
        }
    }
};


void sendRecvTest()
{
    // Mimic communication pattern for sending matrix
    // rows to the next higher process number.
    int rank=-1, size=0;
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
    MPI_Comm_size(MPI_COMM_WORLD, & size);
    if ( size == 1)
    {
        std::cout<<"Cannot run test with just one process"<<std::endl;
    }else
    {
        // We assume a partitioning on regular structured 2D grid distributed
        // in 1D (along the x-axis) with 3 cells along the y-axis.
        // Each processes partition is 5 cells wide.
        // each process additionally stores one overlap (also called
        // ghost or halo layer). Each process (except for the last) needs to
        // send one layer of cells
        // next to the right hand side boundary of its partition to the next
        // high process.
        typedef Dune::Interface::InformationMap InterfaceMap;
        InterfaceMap interface = setup2DSlicedInterface(size, rank);
        Opm::SendReceiveCommunicator comm(MPI_COMM_WORLD, interface);
        GraphConnectionHandle handle(size, rank);
        comm.receiveData(handle);
        comm.sendData(handle);
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

#else
#include<iostream>

void sendRecvTest()
{
    std::cout<<"Cannot run sequential test."<<std::endl;
    std::cout<<"Please compile with MPI and run with multiple processes."<<std::endl;
}
#endif

bool init_function() {
    boost::unit_test::framework::master_test_suite().
        add( BOOST_TEST_CASE(&sendRecvTest) );
    return true;
}
int main(int argc, char** argv)
{
#if HAVE_MPI
    MPI_Init(&argc, &argv);
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif
    ::boost::unit_test::unit_test_main( &init_function, argc, argv );
#if HAVE_MPI
    MPI_Finalize();
#else
#warning "This file needs to compiled with MPI support to make any sense!"
#endif
}
