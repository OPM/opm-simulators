/*
  Copyright 2015, 2022 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 Statoil AS

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
#include <opm/simulators/linalg/MILU.hpp>
#include <opm/common/TimingMacros.hpp>
#include <dune/common/version.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/ilu.hh>

#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/matrixblock.hh>

#include <array>

template <typename T>
T Opm::detail::identityFunctor(const T& t){return t;};

template <typename T>
T Opm::detail::oneFunctor(const T&){return 1.0;};

template <typename T>
T Opm::detail::signFunctor(const T& t){if (t< 0){return -1;} else{return 1;}};

template <typename T>
T Opm::detail::isPositiveFunctor(const T& t){if (t<0){return 0;} else{return 1;}};

template <typename T>
T Opm::detail::absFunctor(const T& t){return std::abs(t);};

namespace Opm
{

MILU_VARIANT convertString2Milu(const std::string& milu)
{
    if( 0 == milu.compare("MILU_1") )
    {
        return MILU_VARIANT::MILU_1;
    }
    if ( 0 == milu.compare("MILU_2") )
    {
        return MILU_VARIANT::MILU_2;
    }
    if ( 0 == milu.compare("MILU_3") )
    {
        return MILU_VARIANT::MILU_3;
    }
    return MILU_VARIANT::ILU;
}

namespace detail
{

template<class M>
void milu0_decomposition(M& A, FieldFunct<M> absFunctor, FieldFunct<M> signFunctor,
                         std::vector<typename M::block_type>* diagonal)
{
    if( diagonal )
    {
        diagonal->reserve(A.N());
    }

    for ( auto irow = A.begin(), iend = A.end(); irow != iend; ++irow)
    {
        auto a_i_end = irow->end();
        auto a_ik    = irow->begin();

        std::array<typename M::field_type, M::block_type::rows> sum_dropped{};

        // Eliminate entries in lower triangular matrix
        // and store factors for L
        for ( ; a_ik.index() < irow.index(); ++a_ik )
        {
            auto k = a_ik.index();
            auto a_kk = A[k].find(k);
            // L_ik = A_kk^-1 * A_ik
            a_ik->rightmultiply(*a_kk);

            // modify the rest of the row, everything right of a_ik
            // a_i* -=a_ik * a_k*
            auto a_k_end = A[k].end();
            auto a_kj = a_kk, a_ij = a_ik;
            ++a_kj; ++a_ij;

            while ( a_kj != a_k_end)
            {
                auto modifier = *a_kj;
                modifier.leftmultiply(*a_ik);

                while( a_ij != a_i_end && a_ij.index() < a_kj.index())
                {
                    ++a_ij;
                }

                if ( a_ij != a_i_end && a_ij.index() == a_kj.index() )
                {
                    // Value is not dropped
                    *a_ij -= modifier;
                    ++a_ij; ++a_kj;
                }
                else
                {
                    auto entry = sum_dropped.begin();
                    for( const auto& row: modifier )
                    {
                        for( const auto& colEntry: row )
                        {
                            *entry += absFunctor(-colEntry);
                        }
                        ++entry;
                    }
                    ++a_kj;
                }
            }
        }

        if ( a_ik.index() != irow.index() )
            OPM_THROW(std::logic_error,
                      "Matrix is missing diagonal for row " + std::to_string(irow.index()));

        int index = 0;
        for(const auto& entry: sum_dropped)
        {
            auto& bdiag = (*a_ik)[index][index];
            bdiag += signFunctor(bdiag) * entry;
            ++index;
        }

        if ( diagonal )
        {
            diagonal->push_back(*a_ik);
        }
        a_ik->invert();   // compute inverse of diagonal block
    }
}

template<class M>
void milun_decomposition(const M& A, int n, MILU_VARIANT milu, M& ILU,
                         Reorderer& ordering, Reorderer& inverseOrdering)
{
    using Map = std::map<std::size_t, int>;

    auto iluRow = ILU.createbegin();

    for(std::size_t i = 0, iend = A.N(); i < iend; ++i)
    {
        auto& orow = A[inverseOrdering[i]];

        Map rowPattern;
        for ( auto col = orow.begin(), cend = orow.end(); col != cend; ++col)
        {
            rowPattern[ordering[col.index()]] = 0;
        }

        for(auto ik = rowPattern.begin(); ik->first < i; ++ik)
        {
            if ( ik->second < n )
            {
                auto& rowk = ILU[ik->first];

                for ( auto kj = rowk.find(ik->first), endk = rowk.end();
                      kj != endk; ++kj)
                {
                    // Assume double and block_type FieldMatrix
                    // first element is misused to store generation number
                    int generation = (*kj)[0][0];
                    if(generation < n)
                    {
                        auto ij = rowPattern.find(kj.index());
                        if ( ij == rowPattern.end() )
                        {
                            rowPattern[ordering[kj.index()]] = generation + 1;
                        }
                    }
                }
            }
        }
        // create the row
        for (const auto& entry : rowPattern)
        {
            iluRow.insert(entry.first);
        }
        ++iluRow;

        // write generation to newly created row.
        auto generationPair = rowPattern.begin();
        for ( auto col = ILU[i].begin(), cend = ILU[i].end(); col != cend;
              ++col, ++generationPair)
        {
            assert(col.index() == generationPair->first);
            (*col)[0][0] = generationPair->second;
        }
    }

    // copy Entries from A
    for(auto iter=A.begin(), iend = A.end(); iter != iend; ++iter)
    {
        auto& newRow = ILU[ordering[iter.index()]];
        // reset stored generation
        for ( auto& col: newRow)
        {
            col = 0;
        }
        // copy row.
        for(auto col = iter->begin(), cend = iter->end(); col != cend; ++col)
        {
            newRow[ordering[col.index()]] = *col;
        }
    }
    // call decomposition on pattern
    switch ( milu )
    {
    case MILU_VARIANT::MILU_1:
        detail::milu0_decomposition ( ILU);
        break;
    case MILU_VARIANT::MILU_2:
        detail::milu0_decomposition ( ILU, identityFunctor<typename M::field_type>,
                                      signFunctor<typename M::field_type>);

        break;
    case MILU_VARIANT::MILU_3:
        detail::milu0_decomposition ( ILU, absFunctor<typename M::field_type>,
                                      signFunctor<typename M::field_type>);
        break;
    case MILU_VARIANT::MILU_4:
        detail::milu0_decomposition ( ILU, identityFunctor<typename M::field_type> ,
                                      isPositiveFunctor<typename M::field_type>);
        break;
    default:
#if DUNE_VERSION_LT(DUNE_GRID, 2, 8)
        bilu0_decomposition( ILU );
#else
        Dune::ILU::blockILU0Decomposition( ILU );
#endif
        break;
    }
}


template double Opm::detail::identityFunctor(const double&);
template double Opm::detail::oneFunctor(const double&);
template double Opm::detail::signFunctor(const double&);
template double Opm::detail::isPositiveFunctor(const double&);
template double Opm::detail::absFunctor(const double&);

#define INSTANCE(...)                                                   \
    template void milu0_decomposition<__VA_ARGS__>                      \
    (__VA_ARGS__&,std::function<double(const double&)>, std::function<double(const double&)>, \
     std::vector<typename __VA_ARGS__::block_type>*);

#define INSTANCE_ILUN(...)                                              \
    template void milun_decomposition(const __VA_ARGS__&, int, MILU_VARIANT, \
                                      __VA_ARGS__&,Reorderer&,Reorderer&);

#define INSTANCE_FULL(...)                      \
    INSTANCE(__VA_ARGS__)                       \
    INSTANCE_ILUN(__VA_ARGS__)

#define INSTANCE_BLOCK(Dim)                                             \
    INSTANCE_FULL(Dune::BCRSMatrix<MatrixBlock<double,Dim,Dim>>)

#define INSTANCE_FM(Dim)                                                \
    INSTANCE_FULL(Dune::BCRSMatrix<Dune::FieldMatrix<double,Dim,Dim>>)

INSTANCE_FM(1)
INSTANCE_FM(2)
INSTANCE_FM(3)
INSTANCE_FM(4)
INSTANCE_FM(5)
INSTANCE_FM(6)

INSTANCE_BLOCK(1)
INSTANCE_BLOCK(2)
INSTANCE_BLOCK(3)
INSTANCE_BLOCK(4)
INSTANCE_BLOCK(5)
INSTANCE_BLOCK(6)

} // end namespace detail

} // end namespace Opm
