/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_AUTODIFFBLOCKARMA_HEADER_INCLUDED
#define OPM_AUTODIFFBLOCKARMA_HEADER_INCLUDED


// #include "AutoDiff.hpp"
// #include <Eigen/Eigen>
// #include <Eigen/Sparse>
#include <armadillo>
#include <vector>
#include <cassert>

namespace AutoDiff
{

    template <typename Scalar>
    class ForwardBlock
    {
    public:
        /// Underlying types for scalar vectors and jacobians.
        typedef arma::Col<Scalar> V;
        typedef arma::SpMat<Scalar> M;

        /// Named constructor pattern used here.
        static ForwardBlock constant(const V& val, const std::vector<int>& blocksizes)
        {
            std::vector<M> jac;
            const int num_elem = val.size();
            const int num_blocks = blocksizes.size();
            // For constants, all jacobian blocks are zero.
            jac.resize(num_blocks);
            for (int i = 0; i < num_blocks; ++i) {
                jac[i] = M(num_elem, blocksizes[i]);
            }
            return ForwardBlock(val, jac);
        }

        static ForwardBlock variable(const int index, const V& val, const std::vector<int>& blocksizes)
        {
            std::vector<M> jac;
            const int num_elem = val.size();
            const int num_blocks = blocksizes.size();
            // First, set all jacobian blocks to zero...
            jac.resize(num_blocks);
            for (int i = 0; i < num_blocks; ++i) {
                jac[i] = M(num_elem, blocksizes[i]);
            }
            // ... then set the one corrresponding to this variable to identity.
            assert(blocksizes[index] == num_elem);
            jac[index].eye(num_elem, num_elem);
            return ForwardBlock(val, jac);
        }

        static ForwardBlock function(const V& val, const std::vector<M>& jac)
        {
            return ForwardBlock(val, jac);
        }

        /// Operator +
        ForwardBlock operator+(const ForwardBlock& rhs)
        {
            std::vector<M> jac = jac_;
            assert(numBlocks() == rhs.numBlocks());
            int num_blocks = numBlocks();
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac[block].n_rows == rhs.jac_[block].n_rows);
                assert(jac[block].n_cols == rhs.jac_[block].n_cols);
                jac[block] += rhs.jac_[block];
            }
            return function(val_ + rhs.val_, jac);
        }

        /// Operator -
        ForwardBlock operator-(const ForwardBlock& rhs)
        {
            std::vector<M> jac = jac_;
            assert(numBlocks() == rhs.numBlocks());
            int num_blocks = numBlocks();
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac[block].n_rows == rhs.jac_[block].n_rows);
                assert(jac[block].n_cols == rhs.jac_[block].n_cols);
                jac[block] -= rhs.jac_[block];
            }
            return function(val_ - rhs.val_, jac);
        }

        /// Operator *
        ForwardBlock operator*(const ForwardBlock& rhs)
        {
            int num_blocks = numBlocks();
            std::vector<M> jac(num_blocks);
            assert(numBlocks() == rhs.numBlocks());
            typedef Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D;
            D D1 = val_.matrix().asDiagonal();
            D D2 = rhs.val_.matrix().asDiagonal();
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac_[block].n_rows == rhs.jac_[block].n_rows);
                assert(jac_[block].n_cols == rhs.jac_[block].n_cols);
                jac[block] = D2*jac_[block] + D1*rhs.jac_[block];
            }
            return function(val_ * rhs.val_, jac);
        }

        /// Operator /
        ForwardBlock operator/(const ForwardBlock& rhs)
        {
            int num_blocks = numBlocks();
            std::vector<M> jac(num_blocks);
            assert(numBlocks() == rhs.numBlocks());
            typedef Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D;
            D D1 = val_.matrix().asDiagonal();
            D D2 = rhs.val_.matrix().asDiagonal();
            D D3 = std::pow(rhs.val_, -2).matrix().asDiagonal();
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac_[block].n_rows == rhs.jac_[block].n_rows);
                assert(jac_[block].n_cols == rhs.jac_[block].n_cols);
                jac[block] = D3 * (D2*jac_[block] - D1*rhs.jac_[block]);
            }
            return function(val_ / rhs.val_, jac);
        }
        /// I/O.
        template <class Ostream>
        Ostream&
        print(Ostream& os) const
        {
            int num_blocks = jac_.size();
            os << "Value =\n" << val_ << "\n\nJacobian =\n";
            for (int i = 0; i < num_blocks; ++i) {
                os << "Sub Jacobian #" << i << '\n' << jac_[i] << "\n";
            }
            return os;
        }

        /// Number of variables or Jacobian blocks.
        int numBlocks() const
        {
            return jac_.size();
        }

        /// Function value
        const V& value() const
        {
            return val_;
        }

        /// Function derivatives
        const std::vector<M>& derivative() const
        {
            return jac_;
        }

    private:
        ForwardBlock(const V& val,
                     const std::vector<M>& jac)
            : val_(val), jac_(jac)
        {
#ifndef NDEBUG
            const int num_elem = val_.size();
            const int num_blocks = jac_.size();
            for (int block = 0; block < num_blocks; ++block) {
                assert(num_elem == jac_[block].n_rows);
            }
#endif
        }

        V val_;
        std::vector<M> jac_;
    };


    template <class Ostream, typename Scalar>
    Ostream&
    operator<<(Ostream& os, const ForwardBlock<Scalar>& fw)
    {
        return fw.print(os);
    }

    /// Multiply with sparse matrix from the left.
    template <typename Scalar>
    ForwardBlock<Scalar> operator*(const typename ForwardBlock<Scalar>::M& lhs,
                                   const ForwardBlock<Scalar>& rhs)
    {
        int num_blocks = rhs.numBlocks();
        std::vector<typename ForwardBlock<Scalar>::M> jac(num_blocks);
        assert(lhs.n_cols == rhs.value().n_rows);
        for (int block = 0; block < num_blocks; ++block) {
            jac[block] = lhs*rhs.derivative()[block];
        }
        typename ForwardBlock<Scalar>::V val = lhs*rhs.value().matrix();
        return ForwardBlock<Scalar>::function(val, jac);
    }


} // namespace Autodiff



#endif // OPM_AUTODIFFBLOCKARMA_HEADER_INCLUDED
