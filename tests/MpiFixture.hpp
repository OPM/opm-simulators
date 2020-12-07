/*
  Copyright 2020 OPM-OP AS
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services.

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
#ifndef OPM_SIMULATOR_TEST_MPIFIXTURE_HPP
#define OPM_SIMULATOR_TEST_MPIFIXTURE_HPP

#include <dune/common/parallel/mpihelper.hh>
#include <boost/test/unit_test.hpp>

class MPIError {
public:
  /** @brief Constructor. */
  MPIError(std::string s, int e) : errorstring(s), errorcode(e){}
  /** @brief The error string. */
  std::string errorstring;
  /** @brief The mpi error code. */
  int errorcode;
};

#ifdef HAVE_MPI
void MPI_err_handler(MPI_Comm *, int *err_code, ...){
  char *err_string=new char[MPI_MAX_ERROR_STRING];
  int err_length;
  MPI_Error_string(*err_code, err_string, &err_length);
  std::string s(err_string, err_length);
  std::cerr << "An MPI Error ocurred:"<<std::endl<<s<<std::endl;
  delete[] err_string;
  throw MPIError(s, *err_code);
}
#endif

struct MPIFixture
{
    MPIFixture()
    {
#if HAVE_MPI
        int m_argc = boost::unit_test::framework::master_test_suite().argc;
        char** m_argv = boost::unit_test::framework::master_test_suite().argv;
        helper = &Dune::MPIHelper::instance(m_argc, m_argv);
#ifdef MPI_2
        MPI_Comm_create_errhandler(MPI_err_handler, &handler);
        MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#else
        MPI_Errhandler_create(MPI_err_handler, &handler);
        MPI_Errhandler_set(MPI_COMM_WORLD, handler);
#endif
#endif
    }
    ~MPIFixture()
    {
#if HAVE_MPI
        MPI_Finalize();
#endif
    }
    Dune::MPIHelper* helper;
#if HAVE_MPI
    MPI_Errhandler handler;
#endif
};
#endif
