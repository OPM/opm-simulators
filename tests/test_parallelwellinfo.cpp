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
#include<config.h>
#include<opm/simulators/wells/ParallelWellInfo.hpp>
#include<vector>
#include<string>
#include<tuple>
#include<ostream>

#define BOOST_TEST_MODULE ParallelWellInfo
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

BOOST_GLOBAL_FIXTURE(MPIFixture);

// Needed for BOOST_CHECK_EQUAL_COLLECTIONS
namespace std
{
std::ostream& operator<<(std::ostream& os, const std::pair<std::string, bool>& p)
{
    return os << "{" << p.first << " "<< p.second << "}";
}
}
namespace Opm
{
std::ostream& operator<<(std::ostream& os, const Opm::ParallelWellInfo& w)
{
    return os << "{" << w.name() << " "<< w.hasLocalCells() << " "<<
        w.isOwner() << "}";
}
}

BOOST_AUTO_TEST_CASE(ParallelWellComparison)
{
    int argc = 0;
    char** argv = nullptr;
    const auto& helper = Dune::MPIHelper::instance(argc, argv);
    std::vector<std::pair<std::string,bool>> pairs;
    if (helper.rank() == 0)
        pairs = {{"Test1", true},{"Test2", true}, {"Test1", false} };
    else
        pairs = {{"Test1", false},{"Test2", true}, {"Test1", true} };

    std::vector<Opm::ParallelWellInfo> well_info;
    well_info.assign(pairs.begin(), pairs.end());

    BOOST_CHECK_EQUAL_COLLECTIONS(pairs.begin(), pairs.end(),
                                  well_info.begin(), well_info.end());

    BOOST_CHECK_EQUAL_COLLECTIONS(well_info.begin(), well_info.end(),
                                  pairs.begin(), pairs.end());

    BOOST_CHECK(well_info[0] < pairs[1]);
    BOOST_CHECK(pairs[0] != well_info[1]);
    BOOST_CHECK(pairs[0] < well_info[1]);
    BOOST_CHECK(well_info[0] == pairs[0]);

    BOOST_CHECK(well_info[0] != well_info[1]);

    Opm::ParallelWellInfo well0, well1;

    BOOST_CHECK(well0 == well1);
#if HAVE_MPI
    BOOST_CHECK(well0.communication()==helper.getLocalCommunicator());
#endif
    Opm::ParallelWellInfo well2("Test", false);
    std::pair<std::string, bool> pwell={"Test", true};
    BOOST_CHECK(well2 < pwell);
    Opm::ParallelWellInfo well3("Test", true);
    BOOST_CHECK(! (well3 < pwell));
    pwell.second = false;
    BOOST_CHECK(! (well3 < pwell));

    if (helper.rank() == 0)
        BOOST_CHECK(well_info[0].communication().size()==1);

#if HAVE_MPI
    Opm::ParallelWellInfo::Communication comm{MPI_COMM_WORLD};

    BOOST_CHECK(well_info[1].communication().size() == comm.size());

    if (helper.rank() > 0)
    {
        BOOST_CHECK(well_info[2].communication().size() == comm.size()-1);
    }
#endif

}

BOOST_AUTO_TEST_CASE(CommunicateAboveBelowSelf)
{
    auto comm = Dune::MPIHelper::getLocalCommunicator();
    Opm::CommunicateAboveBelow commAboveBelow{ comm };
    for(std::size_t count=0; count < 2; ++count)
    {
        std::vector<int> eclIndex = {0, 1, 2, 3, 7 , 8, 10, 11};
        std::vector<double> current(eclIndex.size());
        std::transform(eclIndex.begin(), eclIndex.end(), current.begin(),
                       [](double v){ return 1+10.0*v;});
        commAboveBelow.beginReset();
        for (std::size_t i = 0; i < current.size(); ++i)
        {
            if (i==0)
                commAboveBelow.pushBackEclIndex(-1, eclIndex[i]);
            else
                commAboveBelow.pushBackEclIndex(eclIndex[i-1], eclIndex[i]);
        }
        commAboveBelow.endReset();
        auto above = commAboveBelow.communicateAbove(-10.0, current.data(), current.size());
        BOOST_CHECK(above[0]==-10.0);
        BOOST_CHECK(above.size() == current.size());
        auto a = above.begin()+1;
        std::for_each(current.begin(), current.begin() + (current.size()-1),
                      [&a](double v){ BOOST_CHECK(*(a++) == v);});
        auto below = commAboveBelow.communicateBelow(-10.0, current.data(), current.size());
        BOOST_CHECK(below.back() == -10.0);
        BOOST_CHECK(below.size() == current.size());
        auto b = below.begin();
        std::for_each(current.begin()+1, current.end(),
                      [&b](double v){ BOOST_CHECK(*(b++) == v);});
    }
}


BOOST_AUTO_TEST_CASE(CommunicateAboveBelowSelf1)
{
    auto comm = Dune::MPIHelper::getLocalCommunicator();
    Opm::CommunicateAboveBelow commAboveBelow{ comm };
    for(std::size_t count=0; count < 2; ++count)
    {
        std::vector<int> eclIndex = {0};
        std::vector<double> current(eclIndex.size());
        std::transform(eclIndex.begin(), eclIndex.end(), current.begin(),
                       [](double v){ return 1+10.0*v;});
        commAboveBelow.beginReset();
        for (std::size_t i = 0; i < current.size(); ++i)
        {
            if (i==0)
                commAboveBelow.pushBackEclIndex(-1, eclIndex[i]);
            else
                commAboveBelow.pushBackEclIndex(eclIndex[i-1], eclIndex[i]);
        }
        commAboveBelow.endReset();
        auto above = commAboveBelow.communicateAbove(-10.0, current.data(), current.size());
        BOOST_CHECK(above[0]==-10.0);
        BOOST_CHECK(above.size() == current.size());
        auto a = above.begin()+1;
        std::for_each(current.begin(), current.begin() + (current.size()-1),
                      [&a](double v){ BOOST_CHECK(*(a++) == v);});
        auto below = commAboveBelow.communicateBelow(-10.0, current.data(), current.size());
        BOOST_CHECK(below.back() == -10.0);
        BOOST_CHECK(below.size() == current.size());
        auto b = below.begin();
        std::for_each(current.begin()+1, current.end(),
                      [&b](double v){ BOOST_CHECK(*(b++) == v);});
    }
}

using MPIComm = typename Dune::MPIHelper::MPICommunicator;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
using Communication = Dune::Communication<MPIComm>;
#else
using Communication = Dune::CollectiveCommunication<MPIComm>;
#endif
std::vector<int> createGlobalEclIndex(const Communication& comm)
{
    std::vector<int> globalEclIndex = {0, 1, 2, 3, 7 , 8, 10, 11};
    auto oldSize = globalEclIndex.size();
    std::size_t globalSize = 3 * comm.size();
    auto lastIndex = globalEclIndex.back();
    globalEclIndex.resize(globalSize);
    if ( globalSize > oldSize)
    {
        ++lastIndex;
        for(auto entry = globalEclIndex.begin() + oldSize;
            entry != globalEclIndex.end(); ++entry, ++lastIndex)
        {
            *entry = lastIndex;
        }
    }
    return globalEclIndex;
}

template<class C>
std::vector<double> populateCommAbove(C& commAboveBelow,
                                      const Communication& comm,
                                      const std::vector<int>& globalEclIndex,
                                      const std::vector<double> globalCurrent)
{
    std::vector<double> current(3);

    commAboveBelow.beginReset();
    for (std::size_t i = 0; i < current.size(); ++i)
    {
        auto gi = comm.rank() + comm.size() * i;

        if (gi==0)
        {
            commAboveBelow.pushBackEclIndex(-1, globalEclIndex[gi]);
        }
        else
        {
            commAboveBelow.pushBackEclIndex(globalEclIndex[gi-1], globalEclIndex[gi]);
        }
        current[i] = globalCurrent[gi];
    }
    commAboveBelow.endReset();
    return current;
}

BOOST_AUTO_TEST_CASE(CommunicateAboveBelowParallel)
{
    auto comm = Communication(Dune::MPIHelper::getCommunicator());

    Opm::CommunicateAboveBelow commAboveBelow{ comm };
    for(std::size_t count=0; count < 2; ++count)
    {
        auto globalEclIndex = createGlobalEclIndex(comm);
        std::vector<double> globalCurrent(globalEclIndex.size());
        std::transform(globalEclIndex.begin(), globalEclIndex.end(), globalCurrent.begin(),
                       [](double v){ return 1+10.0*v;});

        auto current = populateCommAbove(commAboveBelow, comm, globalEclIndex, globalCurrent);
        auto above = commAboveBelow.communicateAbove(-10.0, current.data(), current.size());
        if (comm.rank() == 0)
            BOOST_CHECK(above[0]==-10.0);

        BOOST_CHECK(above.size() == current.size());

        for (std::size_t i = 0; i < current.size(); ++i)
        {
            auto gi = comm.rank() + comm.size() * i;
            if (gi > 0)
            {
                BOOST_CHECK(above[i]==globalCurrent[gi-1]);
            }
        }
        auto below = commAboveBelow.communicateBelow(-10.0, current.data(), current.size());
        if (comm.rank() == comm.size() - 1)
            BOOST_CHECK(below.back() == -10.0);

        BOOST_CHECK(below.size() == current.size());

        for (std::size_t i = 0; i < current.size(); ++i)
        {
            auto gi = comm.rank() + comm.size() * i;
            if (gi < globalCurrent.size() - 1)
            {
                BOOST_CHECK(below[i]==globalCurrent[gi+1]);
            }
        }
    }
}
