/*
  Copyright 2025 SINTEF Digital, Mathematics and Cybernetics.

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

#define BOOST_TEST_MODULE TestPropertyTree

#ifndef HAVE_MPI
// Suppress GCC diagnostics of the form
//
//   warning: "HAVE_MPI" is not defined, evaluates to 0
//
// when compiling with "-Wundef".
#define HAVE_MPI 0
#endif  // HAVE_MPI

#include <boost/test/unit_test.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>

#include <opm/common/utility/FileSystem.hpp>

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

BOOST_AUTO_TEST_SUITE(Put_And_Get)

BOOST_AUTO_TEST_CASE(Top_Node_Only)
{
    auto t = Opm::PropertyTree{};

    t.put("a", 1234);
    t.put("b", 123.4);
    t.put("c", std::string { "hello" });
    t.put("d", 12.34f);
    t.put("e", true);

    {
        const auto a = t.get<int>("a");
        BOOST_CHECK_EQUAL(a, 1234);

        const auto aa = t.get<int>("aa", 42);
        BOOST_CHECK_EQUAL(aa, 42);
    }

    {
        const auto b = t.get<double>("b");
        BOOST_CHECK_CLOSE(b, 123.4, 1.0e-8);

        const auto bb = t.get("bb", 2.71828);
        BOOST_CHECK_CLOSE(bb, 2.71828, 1.0e-8);
    }

    {
        const auto c = t.get<std::string>("c");
        BOOST_CHECK_EQUAL(c, "hello");

        const auto cc = t.get("cc", std::string { "world" });
        BOOST_CHECK_EQUAL(cc, "world");
    }

    {
        const auto d = t.get<float>("d");
        BOOST_CHECK_CLOSE(d, 12.34f, 1.0e-6f);

        const auto dd = t.get("dd", -1.618f);
        BOOST_CHECK_CLOSE(dd, -1.618f, 1.0e-6f);
    }

    {
        const auto e = t.get<bool>("e");
        BOOST_CHECK_EQUAL(e, true);

        const auto ee = t.get<bool>("ee", false);
        BOOST_CHECK_EQUAL(ee, false);
    }
}

BOOST_AUTO_TEST_CASE(Size_T)
{
    auto t = Opm::PropertyTree{};

    t.put("s.ramanujan", std::size_t{1729});
    BOOST_CHECK_EQUAL(t.get<std::size_t>("s.ramanujan"), std::size_t{1729});

    t.put("m", static_cast<std::size_t>(-1));
    BOOST_CHECK_EQUAL(t.get<std::size_t>("m"), std::numeric_limits<std::size_t>::max());

    // Conversion: int -> std::size_t
    t.put("n", -1);
    BOOST_CHECK_EQUAL(t.get<std::size_t>("n"), std::numeric_limits<std::size_t>::max());
}

BOOST_AUTO_TEST_CASE(Missing_Keys)
{
    auto t = Opm::PropertyTree{};

    BOOST_CHECK_THROW(t.get<int>("a"), std::exception);
}

BOOST_AUTO_TEST_CASE(Hierarchy)
{
    auto t = Opm::PropertyTree{};

    t.put("a.b.c", 123);

    {
        const auto c = t.get<int>("a.b.c");
        BOOST_CHECK_EQUAL(c, 123);
    }

    {
        const auto a = t.get_child("a");
        const auto c = a.get<int>("b.c");
        BOOST_CHECK_EQUAL(c, 123);
    }

    {
        const auto a = t.get_child("a");
        const auto b = a.get_child("b");
        const auto c = b.get<int>("c");
        BOOST_CHECK_EQUAL(c, 123);
        BOOST_CHECK_CLOSE(a.get<double>("b.d", 3.1415), 3.1415, 1.0e-8);
    }

    {
        const auto f = t.get_child_optional("d.e.f");
        BOOST_CHECK_MESSAGE(! f.has_value(), R"(Node "f" must not exist)");
    }
}

BOOST_AUTO_TEST_SUITE_END() // Put_And_Get

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Load_From_File)

namespace {

    class TempFile
    {
    public:
        TempFile()
            : fname_ { std::filesystem::temp_directory_path() /
                       Opm::unique_path("wrk-%%%%") }
        {}

        ~TempFile()
        {
            std::filesystem::remove_all(this->fname_);
        }

        void append(std::string_view s)
        {
            std::ofstream { this->fname_, std::ios::app } << s;
        }

        std::string name() const
        {
            return this->fname_.generic_string();
        }

    private:
        std::filesystem::path fname_;
    };

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(Top_Node_Only)
{
    auto f = TempFile{};
    f.append(R"({
 "a" : 1234,
 "b" : 123.4,
 "c" : "hello",
 "d" : 12.34
}
)");

    const auto t = Opm::PropertyTree { f.name() };
    {
        const auto a = t.get<int>("a");
        BOOST_CHECK_EQUAL(a, 1234);

        const auto aa = t.get<int>("aa", 42);
        BOOST_CHECK_EQUAL(aa, 42);
    }

    {
        const auto b = t.get<double>("b");
        BOOST_CHECK_CLOSE(b, 123.4, 1.0e-8);

        const auto bb = t.get("bb", 2.71828);
        BOOST_CHECK_CLOSE(bb, 2.71828, 1.0e-8);
    }

    {
        const auto c = t.get<std::string>("c");
        BOOST_CHECK_EQUAL(c, "hello");

        const auto cc = t.get("cc", std::string { "world" });
        BOOST_CHECK_EQUAL(cc, "world");
    }

    {
        const auto d = t.get<float>("d");
        BOOST_CHECK_CLOSE(d, 12.34f, 1.0e-6f);

        const auto dd = t.get("dd", -1.618f);
        BOOST_CHECK_CLOSE(dd, -1.618f, 1.0e-6f);
    }
}

BOOST_AUTO_TEST_CASE(Hierarchy)
{
    auto f = TempFile{};
    f.append(R"({
 "a" : { "b" : { "c" : 123 } }
}
)");

    const auto t = Opm::PropertyTree { f.name() };
    {
        const auto c = t.get<int>("a.b.c");
        BOOST_CHECK_EQUAL(c, 123);
    }

    {
        const auto a = t.get_child("a");
        const auto c = a.get<int>("b.c");
        BOOST_CHECK_EQUAL(c, 123);
    }

    {
        const auto a = t.get_child("a");
        const auto b = a.get_child("b");
        const auto c = b.get<int>("c");
        BOOST_CHECK_EQUAL(c, 123);
        BOOST_CHECK_CLOSE(a.get<double>("b.d", 3.1415), 3.1415, 1.0e-8);
    }

    {
        const auto d_e_f = t.get_child_optional("d.e.f");
        BOOST_CHECK_MESSAGE(! d_e_f.has_value(), R"(Node "f" must not exist)");
    }
}

BOOST_AUTO_TEST_CASE(Vector)
{
    auto f = TempFile{};
    f.append(R"({
 "a" : [ 1, 2, 3, 4 ],
 "b" : [ 11.22, 33.44 ]
})");

    const auto t = Opm::PropertyTree { f.name() };

    {
        const auto a = t.get_child_items_as_vector<int>("a");
        BOOST_REQUIRE_MESSAGE(a.has_value(), R"(Node "a" must exist)");

        const auto expect = std::vector { 1, 2, 3, 4, };
        BOOST_CHECK_EQUAL_COLLECTIONS(a    ->begin(), a    ->end(),
                                      expect.begin(), expect.end());
    }

    {
        const auto aa = t.get_child_items_as_vector<int>("aa");
        BOOST_CHECK_MESSAGE(!aa.has_value(), R"(Node "aa" must NOT exist)");
    }

    {
        const auto b = t.get_child_items_as_vector<double>("b");
        BOOST_REQUIRE_MESSAGE(b.has_value(), R"(Node "b" must exist)");

        const auto expect = std::vector { 11.22, 33.44, };
        BOOST_REQUIRE_EQUAL(b->size(), expect.size());

        for (auto i = 0*b->size(); i < b->size(); ++i) {
            BOOST_TEST_MESSAGE("Element " << i);
            BOOST_CHECK_CLOSE((*b)[i], expect[i], 1.0e-8);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END() // Load_From_File
