/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#define BOOST_TEST_MODULE TestParallelRestart
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/ParallelRestart.hpp>
#include <opm/output/eclipse/RestartValue.hpp>


namespace {


#if HAVE_MPI
Opm::data::Solution getSolution()
{
    Opm::data::Solution sol1;
    sol1.insert("testdata", Opm::UnitSystem::measure::length,
                {1.0, 2.0, 3.0}, Opm::data::TargetType::RESTART_SOLUTION);
    return sol1;
}


Opm::data::Rates getRates()
{
    Opm::data::Rates rat1;
    rat1.set(Opm::data::Rates::opt::wat, 1.0);
    rat1.set(Opm::data::Rates::opt::oil, 2.0);
    rat1.set(Opm::data::Rates::opt::gas, 3.0);
    rat1.set(Opm::data::Rates::opt::polymer, 4.0);
    rat1.set(Opm::data::Rates::opt::solvent, 5.0);
    rat1.set(Opm::data::Rates::opt::energy, 6.0);
    rat1.set(Opm::data::Rates::opt::dissolved_gas, 7.0);
    rat1.set(Opm::data::Rates::opt::vaporized_oil, 8.0);
    rat1.set(Opm::data::Rates::opt::reservoir_water, 9.0);
    rat1.set(Opm::data::Rates::opt::reservoir_oil, 10.0);
    rat1.set(Opm::data::Rates::opt::reservoir_gas, 11.0);
    rat1.set(Opm::data::Rates::opt::productivity_index_water, 12.0);
    rat1.set(Opm::data::Rates::opt::productivity_index_oil, 13.0);
    rat1.set(Opm::data::Rates::opt::productivity_index_gas, 14.0);
    rat1.set(Opm::data::Rates::opt::well_potential_water, 15.0);
    rat1.set(Opm::data::Rates::opt::well_potential_oil, 16.0);
    rat1.set(Opm::data::Rates::opt::well_potential_gas, 17.0);

    return rat1;
}


Opm::data::Connection getConnection()
{
    Opm::data::Connection con1;
    con1.rates = getRates();
    con1.index = 1;
    con1.pressure = 2.0;
    con1.reservoir_rate = 3.0;
    con1.cell_pressure = 4.0;
    con1.cell_saturation_water = 5.0;
    con1.cell_saturation_gas = 6.0;
    con1.effective_Kh = 7.0;
    return con1;
}


Opm::data::Segment getSegment()
{
    Opm::data::Segment seg1;
    seg1.rates = getRates();
    seg1.segNumber = 1;
    seg1.pressure = 2.0;
    return seg1;
}


Opm::data::Well getWell()
{
    Opm::data::Well well1;
    well1.rates = getRates();
    well1.bhp = 1.0;
    well1.thp = 2.0;
    well1.temperature = 3.0;
    well1.control = 4;
    well1.connections.push_back(getConnection());
    well1.segments.insert({0, getSegment()});
    return well1;
}
#endif


}


template<class T>
std::tuple<T,int,int> PackUnpack(const T& in)
{
    auto comm = Dune::MPIHelper::getCollectiveCommunication();
    std::size_t packSize = Opm::Mpi::packSize(in, comm);
    std::vector<char> buffer(packSize);
    int pos1 = 0;
    Opm::Mpi::pack(in, buffer, pos1, comm);
    int pos2 = 0;
    T out;
    Opm::Mpi::unpack(out, buffer, pos2, comm);

    return {out, pos1, pos2};
}


BOOST_AUTO_TEST_CASE(Solution)
{
#if HAVE_MPI
    Opm::data::Solution sol1 = getSolution();
    auto sol2 = PackUnpack(sol1);
    BOOST_CHECK(std::get<1>(sol2) == std::get<2>(sol2));
    BOOST_CHECK(sol1 == std::get<0>(sol2));
#endif
}


BOOST_AUTO_TEST_CASE(Rates)
{
#if HAVE_MPI
    Opm::data::Rates rat1 = getRates();
    auto rat2 = PackUnpack(rat1);
    BOOST_CHECK(std::get<1>(rat2) == std::get<2>(rat2));
    BOOST_CHECK(rat1 == std::get<0>(rat2));
#endif
}


BOOST_AUTO_TEST_CASE(Connection)
{
#if HAVE_MPI
    Opm::data::Connection con1 = getConnection();
    auto con2 = PackUnpack(con1);
    BOOST_CHECK(std::get<1>(con2) == std::get<2>(con2));
    BOOST_CHECK(con1 == std::get<0>(con2));
#endif
}


BOOST_AUTO_TEST_CASE(Segment)
{
#if HAVE_MPI
    Opm::data::Segment seg1 = getSegment();
    auto seg2 = PackUnpack(seg1);
    BOOST_CHECK(std::get<1>(seg2) == std::get<2>(seg2));
    BOOST_CHECK(seg1 == std::get<0>(seg2));
#endif
}


BOOST_AUTO_TEST_CASE(Well)
{
#if HAVE_MPI
    Opm::data::Well well1 = getWell();
    auto well2 = PackUnpack(well1);
    BOOST_CHECK(std::get<1>(well2) == std::get<2>(well2));
    BOOST_CHECK(well1 == std::get<0>(well2));
#endif
}


BOOST_AUTO_TEST_CASE(WellRates)
{
#if HAVE_MPI
    Opm::data::WellRates wells1;
    wells1.insert({"test_well", getWell()});
    auto wells2 = PackUnpack(wells1);
    BOOST_CHECK(std::get<1>(wells2) == std::get<2>(wells2));
    BOOST_CHECK(wells1 == std::get<0>(wells2));
#endif
}


BOOST_AUTO_TEST_CASE(CellData)
{
#if HAVE_MPI
    Opm::data::CellData data1;
    data1.dim = Opm::UnitSystem::measure::length;
    data1.data = {1.0, 2.0, 3.0};
    data1.target = Opm::data::TargetType::RESTART_SOLUTION;
    auto data2 = PackUnpack(data1);
    BOOST_CHECK(std::get<1>(data2) == std::get<2>(data2));
    BOOST_CHECK(data1 == std::get<0>(data2));
#endif
}


BOOST_AUTO_TEST_CASE(RestartKey)
{
#if HAVE_MPI
    Opm::RestartKey key1("key", Opm::UnitSystem::measure::length, true);
    auto key2 = PackUnpack(key1);
    BOOST_CHECK(std::get<1>(key2) == std::get<2>(key2));
    BOOST_CHECK(key1 == std::get<0>(key2));
#endif
}


BOOST_AUTO_TEST_CASE(RestartValue)
{
#if HAVE_MPI
    Opm::data::WellRates wells1;
    wells1.insert({"test_well", getWell()});
    Opm::RestartValue val1(getSolution(), wells1);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}


bool init_unit_test_func()
{
    return true;
}


int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
