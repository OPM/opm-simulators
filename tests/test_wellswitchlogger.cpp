#include <config.h>
#include <dune/common/version.hh>

#define BOOST_TEST_MODULE DistributedCpGridTests
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/simulators/WellSwitchingLogger.hpp>

#if HAVE_MPI
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
#endif

bool
init_unit_test_func()
{
    return true;
}
BOOST_AUTO_TEST_CASE(wellswitchlog)
{
    auto cc = Dune::MPIHelper::getCollectiveCommunication();

    Opm::wellhelpers::WellSwitchingLogger logger(cc);
    std::ostringstream name;
    name <<"Well on rank "<<cc.rank()<<std::flush;

    logger.wellSwitched(name.str(), BHP, THP);

}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
#if HAVE_MPI
    // register a throwing error handler to allow for
    // debugging with "catch throw" in gdb
    MPI_Errhandler handler;
    MPI_Comm_create_errhandler(MPI_err_handler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#endif
    boost::unit_test::unit_test_main(&init_unit_test_func,
                                     argc, argv);
}
