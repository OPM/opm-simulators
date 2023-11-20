#include <config.h>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/graph.hh>

#include <opm/simulators/linalg/GraphColoring.hpp>

#define BOOST_TEST_MODULE GraphColoringTest
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

///! \brief check that all indices are represented in the new ordering.
void checkAllIndices(const std::vector<std::size_t>& ordering)
{
    std::vector<int> counters(ordering.size(), 0);
    for(auto index: ordering)
    {
        ++counters[index];
    }

    for(auto count: counters)
    {
        BOOST_CHECK(count==1);
    }
}

BOOST_AUTO_TEST_CASE(TestWelschPowell)
{
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>;
    using Graph = Dune::Amg::MatrixGraph<Matrix>;
    int N = 10;
    Matrix matrix(N*N, N*N, 5, 0.4, Matrix::implicit);
    for( int j = 0; j < N; j++)
    {
        for(int i = 0; i < N; i++)
        {
            auto index = j*10+i;
            matrix.entry(index,index) = 1;

            if ( i > 0 )
            {
                matrix.entry(index,index-1) = 1;
            }
            if ( i  < N - 1)
            {
                matrix.entry(index,index+1) = 1;
            }

            if ( j > 0 )
            {
                matrix.entry(index,index-N) = 1;
            }
            if ( j  < N - 1)
            {
                matrix.entry(index,index+N) = 1;
            }

        }
    }
    matrix.compress();

    Graph graph(matrix);
    auto colorsTuple = Opm::colorVerticesWelshPowell(graph);
    const auto& colors = std::get<0>(colorsTuple);
    const auto& verticesPerColor = std::get<2>(colorsTuple);
    auto noColors = std::get<1>(colorsTuple);
    auto firstCornerColor = colors[0];
    BOOST_CHECK(noColors == 2);

    // Check for checkerboard coloring

    for( int j = 0, index = 0; j < N; j++)
    {
        auto expectedColor = firstCornerColor;

        for(int i = 0; i < N; i++)
        {
            BOOST_CHECK(colors[index]==expectedColor);
            index++;
            expectedColor = (expectedColor + 1) % 2;
        }
        firstCornerColor=(firstCornerColor + 1) % 2;
    }
    auto newOrder = Opm::reorderVerticesPreserving(colors, noColors, verticesPerColor,
                                                   graph);
    std::vector<std::size_t> colorIndex(noColors, 0);
    std::partial_sum(verticesPerColor.begin(),
                    verticesPerColor.begin()+verticesPerColor.size()-1,
                    colorIndex.begin()+1);

    for (auto vertex : graph)
    {
        BOOST_CHECK(colorIndex[colors[vertex]]++ == newOrder[vertex]);
    }
    checkAllIndices(newOrder);
    newOrder = Opm::reorderVerticesSpheres(colors, noColors, verticesPerColor,
                                           graph, 0);
    checkAllIndices(newOrder);
}

// The following tests verify the graph coloring in the context of revealing which rows
// can be operated on at the same time in the DILU preconditioner
BOOST_AUTO_TEST_CASE(TestColoredDiluParallelisms3x3Matrix)
{
     /*
    Matrix on this form:
    |x  |
    | xx|
    | xx|
    We only expect a DILU dependency from the second to the third row,
    hence row 1 and 2 should have color 0, row 3 should have color 1
    */
    const int N = 3;
    const int bz = 3;
    const int nonZeroes = 5;

    // creating some shorthand typenames
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;

    Matrix testMatrix(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = testMatrix.createbegin(); row != testMatrix.createend(); ++row) {
        if (row.index() == 0) {
            row.insert(row.index());
        }
        else if (row.index() == 1) {
            row.insert(row.index());
            row.insert(row.index() + 1);
        }
        else if (row.index() == 2) {
            row.insert(row.index() - 1);
            row.insert(row.index());
        }
    }

    testMatrix[0][0][0][0] = 3.0;
    testMatrix[1][1][0][0] = 1.0;
    testMatrix[1][2][0][0] = 1.0;
    testMatrix[2][1][0][0] = 1.0;
    testMatrix[2][2][0][0] = 1.0;

    Opm::SparseTable<std::size_t> coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::SYMMETRIC);

    std::vector<std::vector<std::size_t>> correctColor = {{0, 1}, {2}};

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }

    coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::UPPER);
    correctColor = {{0, 2}, {1}};

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }
    coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::LOWER);
    correctColor = {{0, 1}, {2}};

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestColoredDiluParallelisms5x5Simple)
{
     /*
    Test matrix:
    |xxx  |
    |xx   |
    |x xx |
    |   x |
    |   xx|
    */
    const int N = 5;
    const int bz = 3;
    const int nonZeroes = 11;

    // creating some shorthand typenames
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;

    Matrix testMatrix(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = testMatrix.createbegin(); row != testMatrix.createend(); ++row) {
        if (row.index() == 0) {
            row.insert(row.index());
            row.insert(row.index()+1);
            row.insert(row.index()+2);
        }
        else if (row.index() == 1) {
            row.insert(row.index());
            row.insert(row.index() - 1);
        }
        else if (row.index() == 2) {
            row.insert(row.index() - 2);
            row.insert(row.index());
            row.insert(row.index() + 1);
        }
        else if (row.index() == 3) {
            row.insert(row.index());
        }
        else if (row.index() == 4) {
            row.insert(row.index() - 1);
            row.insert(row.index());
        }
    }

    testMatrix[0][0][0][0] = 1.0;
    testMatrix[0][1][0][0] = 1.0;
    testMatrix[0][2][0][0] = 1.0;
    testMatrix[1][0][0][0] = 1.0;
    testMatrix[1][1][0][0] = 1.0;
    testMatrix[2][0][0][0] = 1.0;
    testMatrix[2][2][0][0] = 1.0;
    testMatrix[2][3][0][0] = 1.0;
    testMatrix[3][3][0][0] = 1.0;
    testMatrix[4][3][0][0] = 1.0;
    testMatrix[4][4][0][0] = 1.0;

    Opm::SparseTable<std::size_t> coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::SYMMETRIC);

    std::vector<std::vector<std::size_t>> correctColor = {{0, 3, 4}, {1, 2}};

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }

    coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::UPPER);
    correctColor = {{1, 3, 4}, {2}, {0}};

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }

    coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::LOWER);
    correctColor = {{0, 3}, {1, 2, 4}};

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestColoredDiluParallelisms5x5Tridiagonal)
{
     /*
    Test matrix:
    |xx   |
    |xxx  |
    | xxx |
    |  xxx|
    |   xx|
    The tridiagonal structure will force a strictly serial computation stage
    */
    const int N = 5;
    const int bz = 3;
    const int nonZeroes = 13;

    // creating some shorthand typenames
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;

    Matrix testMatrix(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = testMatrix.createbegin(); row != testMatrix.createend(); ++row) {
        if (row.index() == 0) {
            row.insert(row.index());
            row.insert(row.index()+1);
        }
        else if (row.index() > 0 && row.index() < 4) {
            row.insert(row.index() - 1);
            row.insert(row.index());
            row.insert(row.index() + 1);
        }
        else if (row.index() == 4) {
            row.insert(row.index() - 1);
            row.insert(row.index());
        }
    }

    testMatrix[0][0][0][0] = 1.0;
    testMatrix[0][1][0][0] = 1.0;
    testMatrix[1][0][0][0] = 1.0;
    testMatrix[1][1][0][0] = 1.0;
    testMatrix[1][2][0][0] = 1.0;
    testMatrix[2][1][0][0] = 1.0;
    testMatrix[2][2][0][0] = 1.0;
    testMatrix[2][3][0][0] = 1.0;
    testMatrix[3][2][0][0] = 1.0;
    testMatrix[3][3][0][0] = 1.0;
    testMatrix[3][4][0][0] = 1.0;
    testMatrix[4][3][0][0] = 1.0;
    testMatrix[4][4][0][0] = 1.0;

    Opm::SparseTable<std::size_t> coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::SYMMETRIC);

    std::vector<std::vector<std::size_t>> correctColor = {{0}, {1}, {2}, {3}, {4}};

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }

    coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::LOWER);

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }


    coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::UPPER);
    correctColor = {{4}, {3}, {2}, {1}, {0}};

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestColoredDiluParallelisms5x5Complex)
{
     /*
    Test matrix:
    |xxx x|
    |xx x |
    |x x x|
    | x x |
    |x x x|
    */
    const int N = 5;
    const int bz = 3;
    const int nonZeroes = 15;

    // creating some shorthand typenames
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;

    Matrix testMatrix(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = testMatrix.createbegin(); row != testMatrix.createend(); ++row) {
        if (row.index() == 0) {
            row.insert(row.index());
            row.insert(row.index()+1);
            row.insert(row.index()+2);
            row.insert(row.index()+4);
        }
        else if (row.index() == 1) {
            row.insert(row.index() - 1);
            row.insert(row.index());
            row.insert(row.index() + 2);
        }
        else if (row.index() == 2) {
            row.insert(row.index() - 2);
            row.insert(row.index());
            row.insert(row.index() + 2);
        }
        else if (row.index() == 3) {
            row.insert(row.index() - 2);
            row.insert(row.index());
        }
        else if (row.index() == 4) {
            row.insert(row.index() - 4);
            row.insert(row.index() - 2);
            row.insert(row.index());
        }
    }

    testMatrix[0][0][0][0] = 1.0;
    testMatrix[0][1][0][0] = 1.0;
    testMatrix[0][2][0][0] = 1.0;
    testMatrix[0][4][0][0] = 1.0;
    testMatrix[1][0][0][0] = 1.0;
    testMatrix[1][1][0][0] = 1.0;
    testMatrix[1][3][0][0] = 1.0;
    testMatrix[2][0][0][0] = 1.0;
    testMatrix[2][2][0][0] = 1.0;
    testMatrix[2][4][0][0] = 1.0;
    testMatrix[3][1][0][0] = 1.0;
    testMatrix[3][3][0][0] = 1.0;
    testMatrix[4][0][0][0] = 1.0;
    testMatrix[4][2][0][0] = 1.0;
    testMatrix[4][4][0][0] = 1.0;

    Opm::SparseTable<std::size_t> coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::SYMMETRIC);

    std::vector<std::vector<std::size_t>> correctColor = {{0}, {1, 2}, {3, 4}};

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }

    coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::LOWER);

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }

    coloring = Opm::getMatrixRowColoring(testMatrix, Opm::ColoringType::UPPER);
    correctColor = {{3, 4}, {1, 2}, {0}};
    coloring.print(std::cout);

    for (std::size_t i = 0; i < correctColor.size(); ++i){
        for (std::size_t j = 0; j < correctColor[i].size(); ++j){
            BOOST_CHECK(coloring[i][j] == correctColor[i][j]);
        }
    }
}