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
#include <sstream>
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

BOOST_AUTO_TEST_CASE(Matches_Json_String)
{
    // A tree loaded from a file and one built from the same text in memory
    // must be indistinguishable -- this is the point of fromJsonString().
    const auto json = std::string { R"({
 "a" : 1234,
 "b" : 123.4,
 "c" : "hello",
 "d" : { "e" : { "f" : 42 } },
 "g" : [ 1, 2, 3 ]
}
)" };

    auto f = TempFile{};
    f.append(json);

    const auto fromFile = Opm::PropertyTree { f.name() };
    const auto fromString = Opm::PropertyTree::fromJsonString(json);

    BOOST_CHECK_EQUAL(fromString.get<int>("a"), fromFile.get<int>("a"));
    BOOST_CHECK_CLOSE(fromString.get<double>("b"), fromFile.get<double>("b"), 1.0e-8);
    BOOST_CHECK_EQUAL(fromString.get<std::string>("c"), fromFile.get<std::string>("c"));
    BOOST_CHECK_EQUAL(fromString.get<int>("d.e.f"), fromFile.get<int>("d.e.f"));

    const auto gString = fromString.get_child_items_as_vector<int>("g");
    const auto gFile = fromFile.get_child_items_as_vector<int>("g");
    BOOST_REQUIRE_MESSAGE(gString.has_value(), R"(Node "g" must exist)");
    BOOST_REQUIRE_MESSAGE(gFile.has_value(), R"(Node "g" must exist)");
    BOOST_CHECK_EQUAL_COLLECTIONS(gString->begin(), gString->end(),
                                  gFile  ->begin(), gFile  ->end());
}

BOOST_AUTO_TEST_SUITE_END() // Load_From_File

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(From_Json_String)

BOOST_AUTO_TEST_CASE(Top_Node_Only)
{
    const auto t = Opm::PropertyTree::fromJsonString(R"({
 "a" : 1234,
 "b" : 123.4,
 "c" : "hello",
 "d" : 12.34
}
)");

    BOOST_CHECK_EQUAL(t.get<int>("a"), 1234);
    BOOST_CHECK_EQUAL(t.get<int>("aa", 42), 42);

    BOOST_CHECK_CLOSE(t.get<double>("b"), 123.4, 1.0e-8);
    BOOST_CHECK_CLOSE(t.get("bb", 2.71828), 2.71828, 1.0e-8);

    BOOST_CHECK_EQUAL(t.get<std::string>("c"), "hello");
    BOOST_CHECK_EQUAL(t.get("cc", std::string { "world" }), "world");

    BOOST_CHECK_CLOSE(t.get<float>("d"), 12.34f, 1.0e-6f);
}

BOOST_AUTO_TEST_CASE(Hierarchy)
{
    const auto t = Opm::PropertyTree::fromJsonString(R"({
 "a" : { "b" : { "c" : 123 } }
}
)");

    BOOST_CHECK_EQUAL(t.get<int>("a.b.c"), 123);

    const auto a = t.get_child("a");
    BOOST_CHECK_EQUAL(a.get<int>("b.c"), 123);

    const auto d_e_f = t.get_child_optional("d.e.f");
    BOOST_CHECK_MESSAGE(! d_e_f.has_value(), R"(Node "f" must not exist)");
}

BOOST_AUTO_TEST_CASE(Vector)
{
    const auto t = Opm::PropertyTree::fromJsonString(R"({
 "a" : [ 1, 2, 3, 4 ],
 "b" : [ 11.22, 33.44 ]
})");

    const auto a = t.get_child_items_as_vector<int>("a");
    BOOST_REQUIRE_MESSAGE(a.has_value(), R"(Node "a" must exist)");

    const auto expect = std::vector { 1, 2, 3, 4, };
    BOOST_CHECK_EQUAL_COLLECTIONS(a     ->begin(), a     ->end(),
                                  expect.begin(), expect.end());

    const auto aa = t.get_child_items_as_vector<int>("aa");
    BOOST_CHECK_MESSAGE(! aa.has_value(), R"(Node "aa" must NOT exist)");
}

BOOST_AUTO_TEST_CASE(Solver_Config)
{
    // The motivating use: a linear-solver configuration held inline rather
    // than in a file on disk, as opm-flowgeomechanics builds them.
    const auto t = Opm::PropertyTree::fromJsonString(R"({
 "tol" : 0.01,
 "maxiter" : 200,
 "solver" : "cg",
 "preconditioner" : { "type" : "ParOverILU0", "relaxation" : 0.9 }
})");

    BOOST_CHECK_CLOSE(t.get<double>("tol"), 0.01, 1.0e-8);
    BOOST_CHECK_EQUAL(t.get<int>("maxiter"), 200);
    BOOST_CHECK_EQUAL(t.get<std::string>("solver"), "cg");
    BOOST_CHECK_EQUAL(t.get<std::string>("preconditioner.type"), "ParOverILU0");
    BOOST_CHECK_CLOSE(t.get<double>("preconditioner.relaxation"), 0.9, 1.0e-8);
}

BOOST_AUTO_TEST_CASE(Malformed_Json_Throws)
{
    // Parse errors must surface, not yield a silently empty tree.
    BOOST_CHECK_THROW(Opm::PropertyTree::fromJsonString(R"({ "a" : )"),
                      std::exception);
}

BOOST_AUTO_TEST_SUITE_END() // From_Json_String

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Toml)

namespace {

    class TempTomlFile
    {
    public:
        explicit TempTomlFile(std::string_view content, std::string_view ext = ".toml")
            : fname_ { std::filesystem::temp_directory_path() /
                       Opm::unique_path("wrk-%%%%") }
        {
            fname_ += std::string { ext };
            std::ofstream { this->fname_ } << content;
        }

        ~TempTomlFile()
        {
            std::filesystem::remove_all(this->fname_);
        }

        std::string name() const
        {
            return this->fname_.generic_string();
        }

    private:
        std::filesystem::path fname_;
    };

    std::string asJson(const Opm::PropertyTree& t)
    {
        auto os = std::ostringstream{};
        t.write_json(os, false);
        return os.str();
    }

    std::string asToml(const Opm::PropertyTree& t,
                       const Opm::PropertyTree::TomlCommentProvider& comments = {})
    {
        auto os = std::ostringstream{};
        t.write_toml(os, comments);
        return os.str();
    }

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(Value_Types)
{
    const auto t = Opm::PropertyTree::fromTomlString(R"(
a = 1234
b = 123.4
c = "hello"
d = true
e = 1e-5
f = 'literal \string'
g = [1, 2, 3]
h = { x = 1, y = "z" }
i.j.k = 42

[table]
x = 2.5

[table.sub]
y = -7
)");

    BOOST_CHECK_EQUAL(t.get<int>("a"), 1234);
    BOOST_CHECK_CLOSE(t.get<double>("b"), 123.4, 1.0e-8);
    BOOST_CHECK_EQUAL(t.get<std::string>("c"), "hello");
    BOOST_CHECK_EQUAL(t.get<bool>("d"), true);
    BOOST_CHECK_EQUAL(t.get<std::string>("d"), "true");
    BOOST_CHECK_CLOSE(t.get<double>("e"), 1.0e-5, 1.0e-8);
    BOOST_CHECK_EQUAL(t.get<std::string>("f"), R"(literal \string)");
    BOOST_CHECK_EQUAL(t.get<int>("h.x"), 1);
    BOOST_CHECK_EQUAL(t.get<std::string>("h.y"), "z");
    BOOST_CHECK_EQUAL(t.get<int>("i.j.k"), 42);
    BOOST_CHECK_CLOSE(t.get<double>("table.x"), 2.5, 1.0e-8);
    BOOST_CHECK_EQUAL(t.get<int>("table.sub.y"), -7);

    const auto g = t.get_child_items_as_vector<int>("g");
    BOOST_REQUIRE_MESSAGE(g.has_value(), R"(Node "g" must exist)");
    const auto expect = std::vector { 1, 2, 3, };
    BOOST_CHECK_EQUAL_COLLECTIONS(g->begin(), g->end(), expect.begin(), expect.end());
}

BOOST_AUTO_TEST_CASE(Key_Order_Is_Preserved)
{
    const auto t = Opm::PropertyTree::fromTomlString(R"(
zeta = 1
alpha = 2
[features.SGAS]
[features.PRESSURE]
[features.AAA]
)");

    const auto top = t.get_child_keys();
    BOOST_CHECK_EQUAL(top.size(), std::size_t{3});
    BOOST_CHECK_EQUAL(top[0], "zeta");
    BOOST_CHECK_EQUAL(top[1], "alpha");
    BOOST_CHECK_EQUAL(top[2], "features");

    const auto features = t.get_child("features").get_child_keys();
    const auto expect = std::vector<std::string> { "SGAS", "PRESSURE", "AAA" };
    BOOST_CHECK_EQUAL_COLLECTIONS(features.begin(), features.end(),
                                  expect.begin(), expect.end());
}

BOOST_AUTO_TEST_CASE(Matches_Json)
{
    // A TOML file and a JSON file with the same content must give identical
    // trees, so that either may be used for the same configuration.
    const auto toml = TempTomlFile { R"(
tol = 0.005
maxiter = 20
solver = "bicgstab"
apply = [1, 2]

[preconditioner]
type = "cprw"
add_wells = true
relaxation = 1.0
)" };

    const auto json = TempTomlFile { R"({
 "tol" : 0.005,
 "maxiter" : 20,
 "solver" : "bicgstab",
 "apply" : [1, 2],
 "preconditioner" : {
   "type" : "cprw",
   "add_wells" : true,
   "relaxation" : 1
 }
})", ".json" };

    const auto fromToml = Opm::PropertyTree { toml.name() };
    const auto fromJson = Opm::PropertyTree { json.name() };

    BOOST_CHECK_EQUAL(asJson(fromToml), asJson(fromJson));
}

BOOST_AUTO_TEST_CASE(Malformed_Toml_Throws)
{
    BOOST_CHECK_THROW(Opm::PropertyTree::fromTomlString("a = 1\nb = = 2\n"),
                      std::runtime_error);

    BOOST_CHECK_THROW(Opm::PropertyTree::fromTomlString("a = 1\na = 2\n"),
                      std::runtime_error);

    const auto f = TempTomlFile { "a = [1, 2\n" };
    try {
        const auto t = Opm::PropertyTree { f.name() };
        BOOST_ERROR("Parsing malformed TOML file must throw");
    }
    catch (const std::runtime_error& e) {
        const auto msg = std::string { e.what() };
        BOOST_CHECK_MESSAGE(msg.find(f.name()) != std::string::npos,
                            "Error message must hold the file name: " << msg);
    }
}

BOOST_AUTO_TEST_CASE(File_Names)
{
    BOOST_CHECK(Opm::PropertyTree::isTomlFileName("a.toml"));
    BOOST_CHECK(Opm::PropertyTree::isTomlFileName("dir/a.b.toml"));
    BOOST_CHECK(!Opm::PropertyTree::isTomlFileName(".toml"));
    BOOST_CHECK(!Opm::PropertyTree::isTomlFileName("a.json"));
    BOOST_CHECK(!Opm::PropertyTree::isTomlFileName("toml"));

    BOOST_CHECK(Opm::PropertyTree::isConfigFileName("a.toml"));
    BOOST_CHECK(Opm::PropertyTree::isConfigFileName("a.json"));
    BOOST_CHECK(!Opm::PropertyTree::isConfigFileName("cprw"));
    BOOST_CHECK(!Opm::PropertyTree::isConfigFileName("json"));
    BOOST_CHECK(!Opm::PropertyTree::isConfigFileName(""));
}

BOOST_AUTO_TEST_CASE(Write_Round_Trip)
{
    auto t = Opm::PropertyTree{};
    t.put("int", 42);
    t.put("negative", -3);
    t.put("double", 0.333333333333);
    t.put("small", 1.0e-5);
    t.put("whole", 2.0);
    t.put("flag", true);
    t.put("text", std::string { "cprw" });
    t.put("tricky", std::string { "a \"quoted\" \\ string\nwith newline" });
    t.put("leading_zero", std::string { "007" });
    t.put("leading_dot", std::string { ".35" });
    t.put("empty", std::string {});
    t.put("a.b.c", 1);
    t.put("a.b.d", std::string { "x" });
    t.put("a.e", 3.5);
    t.put("dotted key", std::string { "needs quotes" });

    const auto toml = asToml(t);
    BOOST_TEST_MESSAGE("TOML:\n" << toml);

    const auto u = Opm::PropertyTree::fromTomlString(toml);

    BOOST_CHECK_EQUAL(u.get<int>("int"), 42);
    BOOST_CHECK_EQUAL(u.get<int>("negative"), -3);
    BOOST_CHECK_CLOSE(u.get<double>("double"), 0.333333333333, 1.0e-12);
    BOOST_CHECK_CLOSE(u.get<double>("small"), 1.0e-5, 1.0e-12);
    BOOST_CHECK_CLOSE(u.get<double>("whole"), 2.0, 1.0e-12);
    BOOST_CHECK_EQUAL(u.get<bool>("flag"), true);
    BOOST_CHECK_EQUAL(u.get<std::string>("text"), "cprw");
    BOOST_CHECK_EQUAL(u.get<std::string>("tricky"), "a \"quoted\" \\ string\nwith newline");
    BOOST_CHECK_EQUAL(u.get<std::string>("leading_zero"), "007");
    BOOST_CHECK_EQUAL(u.get<std::string>("leading_dot"), ".35");
    BOOST_CHECK_EQUAL(u.get<std::string>("empty"), "");
    BOOST_CHECK_EQUAL(u.get<int>("a.b.c"), 1);
    BOOST_CHECK_EQUAL(u.get<std::string>("a.b.d"), "x");
    BOOST_CHECK_CLOSE(u.get<double>("a.e"), 3.5, 1.0e-12);
    BOOST_CHECK_EQUAL(u.get_child_keys().back(), "a");

    // Numbers are written as TOML numbers.  The tree stores 2.0 as "2", so
    // it is written as an integer.
    BOOST_CHECK(toml.find("int = 42\n") != std::string::npos);
    BOOST_CHECK(toml.find("whole = 2\n") != std::string::npos);
    BOOST_CHECK(toml.find("double = 0.333333333333\n") != std::string::npos);
    BOOST_CHECK(toml.find("flag = true\n") != std::string::npos);
    BOOST_CHECK(toml.find("\"dotted key\" = \"needs quotes\"\n") != std::string::npos);
    BOOST_CHECK(toml.find("[a.b]\n") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(Write_Arrays)
{
    const auto t = Opm::PropertyTree::fromJsonString(R"({
 "v" : [ 1, 2.5, "x" ],
 "t" : [ { "a" : 1 }, { "a" : 2 } ],
 "sub" : { "w" : [ 3, 4 ] }
})");

    const auto toml = asToml(t);
    BOOST_TEST_MESSAGE("TOML:\n" << toml);

    const auto u = Opm::PropertyTree::fromTomlString(toml);
    BOOST_CHECK_EQUAL(asJson(u), asJson(t));
}

BOOST_AUTO_TEST_CASE(Write_Comments)
{
    auto t = Opm::PropertyTree{};
    t.put("tol", 0.01);
    t.put("preconditioner.type", std::string { "ilu0" });

    const auto toml = asToml(t, [](const std::string& key) -> std::string
    {
        if (key == "tol") {
            return "Tolerance.\nSecond line.";
        }
        if (key == "preconditioner") {
            return "The preconditioner.";
        }
        return {};
    });
    BOOST_TEST_MESSAGE("TOML:\n" << toml);

    BOOST_CHECK(toml.find("# Tolerance.\n# Second line.\ntol = 0.01\n") != std::string::npos);
    BOOST_CHECK(toml.find("# The preconditioner.\n[preconditioner]\ntype = \"ilu0\"\n") != std::string::npos);

    const auto u = Opm::PropertyTree::fromTomlString(toml);
    BOOST_CHECK_EQUAL(asJson(u), asJson(t));
}

BOOST_AUTO_TEST_CASE(Nesting_Limit_After_Table_Headers)
{
    // The parser limits the nesting depth of arrays to 30 to protect
    // against stack overflow.  Table headers must not affect this limit.
    const auto document = [](const int numHeaders, const int depth)
    {
        auto toml = std::string{};
        for (auto i = 0; i < numHeaders; ++i) {
            toml += "[table" + std::to_string(i) + "]\n";
        }
        return toml + "x = " + std::string(depth, '[') + std::string(depth, ']') + '\n';
    };

    BOOST_CHECK_NO_THROW(Opm::PropertyTree::fromTomlString(document(0, 30)));
    BOOST_CHECK_NO_THROW(Opm::PropertyTree::fromTomlString(document(1, 30)));
    BOOST_CHECK_NO_THROW(Opm::PropertyTree::fromTomlString(document(5, 30)));

    BOOST_CHECK_THROW(Opm::PropertyTree::fromTomlString(document(0, 31)), std::runtime_error);
    BOOST_CHECK_THROW(Opm::PropertyTree::fromTomlString(document(1, 31)), std::runtime_error);
    BOOST_CHECK_THROW(Opm::PropertyTree::fromTomlString(document(200, 200)), std::runtime_error);

    // The same holds for inline tables.
    auto inlineTables = std::string { "[a]\n[b]\nx = " };
    for (auto i = 0; i < 31; ++i) {
        inlineTables += "{ y = ";
    }
    inlineTables += "1" + std::string(31, '}') + '\n';
    BOOST_CHECK_THROW(Opm::PropertyTree::fromTomlString(inlineTables), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END() // Toml
