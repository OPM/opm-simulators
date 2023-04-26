/*
  The content of this file is based on the file dune/istl/matrixmarket.hh in
  the Dune module dune-istl.

  The license of this file is therefore the same as that of Dune, see
  https://www.dune-project.org/about/license/
*/

#ifndef OPM_MATRIXMARKETSPECIALIZATIONS_HEADER_INCLUDED
#define OPM_MATRIXMARKETSPECIALIZATIONS_HEADER_INCLUDED

#include <dune/istl/matrixmarket.hh>

namespace Opm
{
template<typename T, int i, int j>
class MatrixBlock;
}

namespace Dune
{

namespace MatrixMarketImpl
{

    template <typename T, int i, int j, typename A>
    struct mm_header_printer<BCRSMatrix<Opm::MatrixBlock<T,i,j>, A>>
    {
        static void print(std::ostream& os)
        {
            os << "%%MatrixMarket matrix coordinate ";
            os << mm_numeric_type<T>::str() << " general" << std::endl;
        }
    };

    template <typename T, int i, int j, typename A>
    struct mm_block_structure_header<BCRSMatrix<Opm::MatrixBlock<T,i,j>, A>>
    {
        using M = BCRSMatrix<Opm::MatrixBlock<T,i,j>, A>;
        static void print(std::ostream& os, const M&)
        {
            os << "% ISTL_STRUCT blocked ";
            os << i << " " << j << std::endl;
        }
    };
} // namespace MatrixMarketImpl

namespace MatrixMarketImpl
{
    template <typename T, int i, int j, typename A>
    struct mm_multipliers<BCRSMatrix<Opm::MatrixBlock<T,i,j>, A>>
    {
        enum {
            rows = i,
            cols = j
        };
    };
} // namespace MatrixMarketImpl

} // namespace Dune

#endif
