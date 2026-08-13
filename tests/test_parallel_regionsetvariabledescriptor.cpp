/*
  Copyright 2026 Equinor AS

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

#define BOOST_TEST_MODULE TestParallelRegionSetVariableDescriptor

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/simulators/utils/ParallelRegionsetVariableDescriptor.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <array>
#include <cstddef>
#include <deque>
#include <iostream>
#include <forward_list>
#include <string>
#include <string_view>
#include <vector>

namespace {

    struct MPIError
    {
        MPIError(std::string_view errstr, const int ec)
            : errorstring { errstr }
            , errorcode   { ec }
        {}

        std::string errorstring;
        int errorcode;
    };

    void MPI_err_handler(MPI_Comm*, int* err_code, ...)
    {
        std::array<char, MPI_MAX_ERROR_STRING> err_string_vec{'\0'};
        auto err_length = 0;

        MPI_Error_string(*err_code, err_string_vec.data(), &err_length);

        auto err_string = std::string_view {
            err_string_vec.data(),
            static_cast<std::string_view::size_type>(err_length)
        };

        std::cerr << "An MPI Error occurred:\n  -> " << err_string << '\n';

        throw MPIError { err_string, *err_code };
    }

    // Register a throwing error handler to allow for debugging with
    //
    //   catch throw
    //
    // in GDB.
    void register_error_handler()
    {
        MPI_Errhandler handler{};

        MPI_Comm_create_errhandler(MPI_err_handler, &handler);
        MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
    }

    class NProc_Is
    {
    public:
        explicit NProc_Is(const int expectNP)
            : expectNP_ { expectNP }
        {}

        boost::test_tools::assertion_result
        operator()(boost::unit_test::test_unit_id) const
        {
            auto comm = Opm::Parallel::Communication {
                Dune::MPIHelper::getCommunicator()
            };

            if (comm.size() == this->expectNP_) {
                return true;
            }

            boost::test_tools::assertion_result response(false);
            response.message() << "Number of MPI processes ("
                               << comm.size()
                               << ") differs from expected "
                               << this->expectNP_;

            return response;
        }

    private:
        int expectNP_{};
    };

    bool init_unit_test_func()
    {
        return true;
    }

} // Anonymous namespace

BOOST_AUTO_TEST_SUITE(Single_Regset)

BOOST_AUTO_TEST_SUITE(NProc_2, * boost::unit_test::precondition(NProc_Is{2}))

BOOST_AUTO_TEST_CASE(Empty)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 2);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();
    d.addRegionSet(/* maxRegionID = */ 0); // FIELD or similar.
    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 2);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* maxRegionID = */ 1); // FIELD or similar.
    }
    else {
        d.addRegionSet(/* maxRegionID = */ 3);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{4});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Common_Max_Declared)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 2);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* declaredMaxRegionID = */ 3,
                   std::array { 1, 1, 1, 2, 2, 0, 1, });

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{4});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Common_Max_Value)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 2);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* declaredMaxRegionID = */ 1,
                   std::array { 1, 1, 1, 2, 2, 5, 1, });

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{6});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Distinct_Max_Declared)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 2);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 1, 1, 1, });
    }
    else {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 2, 0, 1, });
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{4});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Distinct_Max_Value)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 2);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 1, 1, 5, });
    }
    else {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 2, 0, 1, });
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{6});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_2

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_3, * boost::unit_test::precondition(NProc_Is{3}))

BOOST_AUTO_TEST_CASE(Empty)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();
    d.addRegionSet(/* maxRegionID = */ 0); // FIELD or similar.
    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* maxRegionID = */ 1);
    }
    else if (comm.rank() == 1) {
        d.addRegionSet(/* maxRegionID = */ 3);
    }
    else {
        d.addRegionSet(/* maxRegionID = */ 5);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{6});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Common_Max_Declared)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* declaredMaxRegionID = */ 5,
                   std::array { 1, 1, 1, 2, 2, 0, 1, });

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{6});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Common_Max_Value)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* declaredMaxRegionID = */ 1,
                   std::array { 1, 1, 1, 2, 2, 5, 1, });

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{6});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Distinct_Max_Declared)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 1, 1, 1, });
    }
    else if (comm.rank() == 1) {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 2, 0, 1, });
    }
    else {
        d.addRegionSet(/* declaredMaxRegionID = */ 4,
                       std::array { 3, 2, 1, });
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{5});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Distinct_Max_Value)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 1, 1, 5, });
    }
    else if (comm.rank() == 1) {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 2, 0, 1, });
    }
    else {
        d.addRegionSet(/* declaredMaxRegionID = */ 4,
                       std::array { 3, 2, 1, });
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{6});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_3

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_4, * boost::unit_test::precondition(NProc_Is{4}))

BOOST_AUTO_TEST_CASE(Empty)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();
    d.addRegionSet(/* maxRegionID = */ 0); // FIELD or similar.
    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* maxRegionID = */ 1);
    }
    else if (comm.rank() == 1) {
        d.addRegionSet(/* maxRegionID = */ 3);
    }
    else if (comm.rank() == 2) {
        d.addRegionSet(/* maxRegionID = */ 5);
    }
    else {
        d.addRegionSet(/* maxRegionID = */ 7);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{8});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Common_Max_Declared)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* declaredMaxRegionID = */ 7,
                   std::array { 1, 1, 1, 2, 2, 0, 1, });

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{8});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Common_Max_Value)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* declaredMaxRegionID = */ 1,
                   std::array { 1, 1, 1, 2, 2, 5, 1, });

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{6});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Distinct_Max_Declared)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 1, 1, 1, });
    }
    else if (comm.rank() == 1) {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 2, 0, 1, });
    }
    else if (comm.rank() == 2) {
        d.addRegionSet(/* declaredMaxRegionID = */ 4,
                       std::array { 3, 2, 1, });
    }
    else {
        d.addRegionSet(/* declaredMaxRegionID = */ 5,
                       std::array { 4, 4, 2, });
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{6});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Distinct_Max_Value)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 1, 1, 5, });
    }
    else if (comm.rank() == 1) {
        d.addRegionSet(/* declaredMaxRegionID = */ 3,
                       std::array { 2, 0, 1, });
    }
    else if (comm.rank() == 2) {
        d.addRegionSet(/* declaredMaxRegionID = */ 4,
                       std::array { 3, 2, 1, });
    }
    else {
        d.addRegionSet(/* declaredMaxRegionID = */ 4,
                       std::array { 4, 4, 8, });
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{1});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{9});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_4

BOOST_AUTO_TEST_SUITE_END()     // Single_Regset

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Two_Regsets)

BOOST_AUTO_TEST_SUITE(NProc_2, * boost::unit_test::precondition(NProc_Is{2}))

BOOST_AUTO_TEST_CASE(Single_Region)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto d = Opm::ParallelRegionsetVariableDescriptor{ comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{1});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto d = Opm::ParallelRegionsetVariableDescriptor{ comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* maxRegionID = */ 5); // Supports regions 0..5
    }
    else {
        d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
    }

    d.addRegionSet(/* maxRegionID = */ 3); // Supports regions 0..3

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{10});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{6});
}

BOOST_AUTO_TEST_CASE(Single_Region_Scan_Regions)
{
    const auto reg_1 = std::vector { 0, 0, 0, 0, 0, };
    const auto reg_2 = std::deque { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, };

    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto d = Opm::ParallelRegionsetVariableDescriptor{ comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_1);
    d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_2);

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{1});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_I)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto d = Opm::ParallelRegionsetVariableDescriptor{ comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 1, 1, 3, };
        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{10});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{6});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_II)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto d = Opm::ParallelRegionsetVariableDescriptor{ comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 2, 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 2, 1, 1, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{10});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{4});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_III)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto d = Opm::ParallelRegionsetVariableDescriptor{ comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 2, 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 2, 1, 1, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{8});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{4});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_IV)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    auto d = Opm::ParallelRegionsetVariableDescriptor{ comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 17, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 29, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 2, 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 2, 1, 1, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 17, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 29, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{48});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{18});
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_2

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_3, * boost::unit_test::precondition(NProc_Is{3}))

BOOST_AUTO_TEST_CASE(Empty)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();
    d.addRegionSet(/* maxRegionID = */ 0); // FIELD or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // FIELD or similar.
    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{1});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* maxRegionID = */ 5); // Supports regions 0..5
    }
    else if (comm.rank() == 1) {
        d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
    }
    else {
        d.addRegionSet(/* maxRegionID = */ 3); // Supports regions 0..3
    }

    d.addRegionSet(/* maxRegionID = */ 3); // Supports regions 0..3

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{10});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{6});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_I)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_2);
    }
    else if (comm.rank() == 1) {
        const auto reg_1 = std::vector { 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 1, 1, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 3, 3, };
        const auto reg_2 = std::array { 3, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{10});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{6});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_II)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_2);
    }
    else if (comm.rank() == 1) {
        const auto reg_1 = std::vector { 2, 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 2, 1, 1, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 3, };
        const auto reg_2 = std::array { 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{10});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{4});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_III)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_2);
    }
    else if (comm.rank() == 1) {
        const auto reg_1 = std::vector { 2, 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 2, 1, 1, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 3, };
        const auto reg_2 = std::array { 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{8});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{4});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_IV)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 17, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 29, reg_2);
    }
    else if (comm.rank() == 1) {
        const auto reg_1 = std::vector { 2, 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 2, 1, 1, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 17, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 29, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 3, };
        const auto reg_2 = std::array { 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 17, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 29, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{48});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{18});
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_3

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_4, * boost::unit_test::precondition(NProc_Is{4}))

BOOST_AUTO_TEST_CASE(Empty)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();
    d.addRegionSet(/* maxRegionID = */ 0); // FIELD or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // FIELD or similar.
    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{1});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* maxRegionID = */ 5); // Supports regions 0..5
    }
    else if (comm.rank() == 1) {
        d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
    }
    else if (comm.rank() == 2) {
        d.addRegionSet(/* maxRegionID = */ 3); // Supports regions 0..3
    }
    else {
        d.addRegionSet(/* maxRegionID = */ 7); // Supports regions 0..7
    }

    d.addRegionSet(/* maxRegionID = */ 3); // Supports regions 0..3

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{12});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{8});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_I)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_2);
    }
    else if (comm.rank() == 1) {
        const auto reg_1 = std::vector { 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 1, 1, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_2);
    }
    else if (comm.rank() == 2) {
        const auto reg_1 = std::vector { 3, };
        const auto reg_2 = std::array { 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 4, };
        const auto reg_2 = std::array { 4, };

        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{11});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{6});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_II)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_2);
    }
    else if (comm.rank() == 1) {
        const auto reg_1 = std::vector { 2, 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 2, 1, 1, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_2);
    }
    else if (comm.rank() == 2) {
        const auto reg_1 = std::vector { 3, };
        const auto reg_2 = std::array { 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 4, };
        const auto reg_2 = std::array { 4, };

        d.addRegionSet(/* declaredMaxRegionID = */ 3, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 5, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{11});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{5});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_III)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_2);
    }
    else if (comm.rank() == 1) {
        const auto reg_1 = std::vector { 2, 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 2, 1, 1, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_2);
    }
    else if (comm.rank() == 2) {
        const auto reg_1 = std::vector { 3, };
        const auto reg_2 = std::array { 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 4, };
        const auto reg_2 = std::array { 4, };

        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 1, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{10});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{5});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions_Scan_Regions_IV)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        const auto reg_1 = std::vector { 1, 1, 2, };
        const auto reg_2 = std::array { 1, 1, 2, };

        d.addRegionSet(/* declaredMaxRegionID = */ 17, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 29, reg_2);
    }
    else if (comm.rank() == 1) {
        const auto reg_1 = std::vector { 2, 2, 1, 1, 3, };
        const auto reg_2 = std::array { 2, 2, 1, 1, 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 17, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 29, reg_2);
    }
    else if (comm.rank() == 2) {
        const auto reg_1 = std::vector { 3, };
        const auto reg_2 = std::array { 3, };

        d.addRegionSet(/* declaredMaxRegionID = */ 17, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 29, reg_2);
    }
    else {
        const auto reg_1 = std::vector { 4, };
        const auto reg_2 = std::array { 4, };

        d.addRegionSet(/* declaredMaxRegionID = */ 17, reg_1);
        d.addRegionSet(/* declaredMaxRegionID = */ 29, reg_2);
    }

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{2});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{48});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{18});
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_4

BOOST_AUTO_TEST_SUITE_END()     // Two_Regsets

// ===========================================================================

BOOST_AUTO_TEST_SUITE(Multiple_Regsets)

BOOST_AUTO_TEST_SUITE(NProc_2, * boost::unit_test::precondition(NProc_Is{2}))

BOOST_AUTO_TEST_CASE(Single_Region)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 2);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{5});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{5});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{1});
    BOOST_CHECK_EQUAL(d.startIndex(2), std::size_t{2});
    BOOST_CHECK_EQUAL(d.startIndex(3), std::size_t{3});
    BOOST_CHECK_EQUAL(d.startIndex(4), std::size_t{4});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 2);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* maxRegionID = */ 4); // Supports regions 0..4
        d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
    }
    else {
        d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
        d.addRegionSet(/* maxRegionID = */ 3); // Supports regions 0..3
    }

    d.addRegionSet(/* maxRegionID = */ 2); // Supports regions 0..2
    d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
    d.addRegionSet(/* maxRegionID = */ 0); // Supports regions 0..0

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{5});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{15});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{ 0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{ 5});
    BOOST_CHECK_EQUAL(d.startIndex(2), std::size_t{ 9});
    BOOST_CHECK_EQUAL(d.startIndex(3), std::size_t{12});
    BOOST_CHECK_EQUAL(d.startIndex(4), std::size_t{14});
}

BOOST_AUTO_TEST_CASE(Scan_Regions)
{
    const auto reg_1 = std::vector<int> {};
    const auto reg_2 = std::vector { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, };
    const auto reg_5 = std::deque { 11, 22, 33, 17, 29, };
    const auto reg_6 = std::array { 0, 1, 0, 2, 3, 0, 1, };

    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 2);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* declaredMaxRegionID = */  5, reg_1); // max ID =  5
    d.addRegionSet(/* declaredMaxRegionID = */ 42, reg_2); // max ID = 42

    if (comm.rank() == 0) {
        const auto reg_3 = std::vector { 5, 9, 26, };
        const auto reg_4 = std::forward_list { 0, 0, 0 };

        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_3); // max ID = 26
        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_4); // max ID =  0
    }
    else {
        const auto reg_3 = std::vector { 3, 14, 1, };
        const auto reg_4 = std::forward_list { 0, 0, 0 };

        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_3); // max ID = 26
        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_4); // max ID =  0
    }

    d.addRegionSet(/* declaredMaxRegionID = */ 11, reg_5); // max ID = 33
    d.addRegionSet(/* declaredMaxRegionID = */  5, reg_6); // max ID =  5

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{6});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{117});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{  0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{  6});
    BOOST_CHECK_EQUAL(d.startIndex(2), std::size_t{ 49});
    BOOST_CHECK_EQUAL(d.startIndex(3), std::size_t{ 76});
    BOOST_CHECK_EQUAL(d.startIndex(4), std::size_t{ 77});
    BOOST_CHECK_EQUAL(d.startIndex(5), std::size_t{111});
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_2

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_3, * boost::unit_test::precondition(NProc_Is{3}))

BOOST_AUTO_TEST_CASE(Single_Region)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{5});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{5});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{1});
    BOOST_CHECK_EQUAL(d.startIndex(2), std::size_t{2});
    BOOST_CHECK_EQUAL(d.startIndex(3), std::size_t{3});
    BOOST_CHECK_EQUAL(d.startIndex(4), std::size_t{4});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* maxRegionID = */ 5); // Supports regions 0..5
        d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
    }
    else if (comm.rank() == 1) {
        d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
        d.addRegionSet(/* maxRegionID = */ 3); // Supports regions 0..3
    }
    else {
        d.addRegionSet(/* maxRegionID = */ 3); // Supports regions 0..3
        d.addRegionSet(/* maxRegionID = */ 7); // Supports regions 0..7
    }

    d.addRegionSet(/* maxRegionID = */ 2); // Supports regions 0..2
    d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
    d.addRegionSet(/* maxRegionID = */ 0); // Supports regions 0..0

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{5});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{20});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{ 0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{ 6});
    BOOST_CHECK_EQUAL(d.startIndex(2), std::size_t{14});
    BOOST_CHECK_EQUAL(d.startIndex(3), std::size_t{17});
    BOOST_CHECK_EQUAL(d.startIndex(4), std::size_t{19});
}

BOOST_AUTO_TEST_CASE(Scan_Regions)
{
    const auto reg_1 = std::vector<int> {};
    const auto reg_2 = std::vector { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, };
    const auto reg_5 = std::deque { 11, 22, 33, 17, 29, };
    const auto reg_6 = std::array { 0, 1, 0, 2, 3, 0, 1, };

    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 3);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* declaredMaxRegionID = */  5, reg_1); // max ID =  5
    d.addRegionSet(/* declaredMaxRegionID = */ 42, reg_2); // max ID = 42

    if (comm.rank() == 0) {
        const auto reg_3 = std::vector { 5, 9, 26, };
        const auto reg_4 = std::forward_list { 0, 0, 0 };

        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_3); // max ID = 26
        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_4); // max ID =  0
    }
    else if (comm.rank() == 1) {
        const auto reg_3 = std::vector { 3, 14, 1, };
        const auto reg_4 = std::forward_list { 0, 0, 0 };

        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_3); // max ID = 26
        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_4); // max ID =  0
    }
    else {
        const auto reg_3 = std::vector { 5, 3, 5, };
        const auto reg_4 = std::forward_list<int> {};

        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_3); // max ID = 5
        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_4); // max ID = 0
    }

    d.addRegionSet(/* declaredMaxRegionID = */ 11, reg_5); // max ID = 33
    d.addRegionSet(/* declaredMaxRegionID = */  5, reg_6); // max ID =  5

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{6});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{117});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{  0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{  6});
    BOOST_CHECK_EQUAL(d.startIndex(2), std::size_t{ 49});
    BOOST_CHECK_EQUAL(d.startIndex(3), std::size_t{ 76});
    BOOST_CHECK_EQUAL(d.startIndex(4), std::size_t{ 77});
    BOOST_CHECK_EQUAL(d.startIndex(5), std::size_t{111});
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_3

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(NProc_4, * boost::unit_test::precondition(NProc_Is{4}))

BOOST_AUTO_TEST_CASE(Single_Region)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.
    d.addRegionSet(/* maxRegionID = */ 0); // "FIELD" or similar.

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{5});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{5});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{1});
    BOOST_CHECK_EQUAL(d.startIndex(2), std::size_t{2});
    BOOST_CHECK_EQUAL(d.startIndex(3), std::size_t{3});
    BOOST_CHECK_EQUAL(d.startIndex(4), std::size_t{4});
}

BOOST_AUTO_TEST_CASE(Multiple_Regions)
{
    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    if (comm.rank() == 0) {
        d.addRegionSet(/* maxRegionID = */ 5); // Supports regions 0..5
        d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
    }
    else if (comm.rank() == 1) {
        d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
        d.addRegionSet(/* maxRegionID = */ 3); // Supports regions 0..3
    }
    else if (comm.rank() == 2) {
        d.addRegionSet(/* maxRegionID = */ 3); // Supports regions 0..3
        d.addRegionSet(/* maxRegionID = */ 7); // Supports regions 0..7
    }
    else {
        d.addRegionSet(/* maxRegionID = */ 7); // Supports regions 0..7
        d.addRegionSet(/* maxRegionID = */ 9); // Supports regions 0..9
    }

    d.addRegionSet(/* maxRegionID = */ 2); // Supports regions 0..2
    d.addRegionSet(/* maxRegionID = */ 1); // Supports regions 0..1
    d.addRegionSet(/* maxRegionID = */ 0); // Supports regions 0..0

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{5});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{24});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{ 0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{ 8});
    BOOST_CHECK_EQUAL(d.startIndex(2), std::size_t{18});
    BOOST_CHECK_EQUAL(d.startIndex(3), std::size_t{21});
    BOOST_CHECK_EQUAL(d.startIndex(4), std::size_t{23});
}

BOOST_AUTO_TEST_CASE(Scan_Regions)
{
    const auto reg_1 = std::vector<int> {};
    const auto reg_2 = std::vector { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, };
    const auto reg_5 = std::deque { 11, 22, 33, 17, 29, };
    const auto reg_6 = std::array { 0, 1, 0, 2, 3, 0, 1, };

    auto comm = Opm::Parallel::Communication {
        Dune::MPIHelper::getCommunicator()
    };

    BOOST_REQUIRE_EQUAL(comm.size(), 4);

    auto d = Opm::ParallelRegionsetVariableDescriptor { comm };

    d.prepareDescriptorSet();

    d.addRegionSet(/* declaredMaxRegionID = */  5, reg_1); // max ID =  5
    d.addRegionSet(/* declaredMaxRegionID = */ 42, reg_2); // max ID = 42

    if (comm.rank() == 0) {
        const auto reg_3 = std::vector { 5, 9, 26, };
        const auto reg_4 = std::forward_list { 0, 0, 0 };

        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_3); // max ID = 26
        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_4); // max ID = 0
    }
    else if (comm.rank() == 1) {
        const auto reg_3 = std::vector { 3, 14, 1, };
        const auto reg_4 = std::forward_list { 0, 0, 0 };

        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_3); // max ID = 14
        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_4); // max ID = -1
    }
    else if (comm.rank() == 2) {
        const auto reg_3 = std::vector<int> {};
        const auto reg_4 = std::forward_list<int> {};

        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_3); // max ID = -1
        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_4); // max ID = -1
    }
    else if (comm.rank() == 3) {
        const auto reg_3 = std::vector<int> {};
        const auto reg_4 = std::forward_list<int> {};

        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_3); // max ID = -1
        d.addRegionSet(/* declaredMaxRegionID = */ 0, reg_4); // max ID = -1
    }

    d.addRegionSet(/* declaredMaxRegionID = */ 11, reg_5); // max ID = 33
    d.addRegionSet(/* declaredMaxRegionID = */  5, reg_6); // max ID =  5

    d.finaliseDescriptorSet();

    BOOST_CHECK_EQUAL(d.numRegionSets(), std::size_t{6});
    BOOST_CHECK_EQUAL(d.numVariableSlots(), std::size_t{117});
    BOOST_CHECK_EQUAL(d.startIndex(0), std::size_t{  0});
    BOOST_CHECK_EQUAL(d.startIndex(1), std::size_t{  6});
    BOOST_CHECK_EQUAL(d.startIndex(2), std::size_t{ 49});
    BOOST_CHECK_EQUAL(d.startIndex(3), std::size_t{ 76});
    BOOST_CHECK_EQUAL(d.startIndex(4), std::size_t{ 77});
    BOOST_CHECK_EQUAL(d.startIndex(5), std::size_t{111});
}

BOOST_AUTO_TEST_SUITE_END()     // NProc_4

BOOST_AUTO_TEST_SUITE_END()     // Multiple_Regsets

// ===========================================================================

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    register_error_handler();

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
