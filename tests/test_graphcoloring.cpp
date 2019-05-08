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
