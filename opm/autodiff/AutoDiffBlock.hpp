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

#ifndef OPM_AUTODIFFBLOCK_HEADER_INCLUDED
#define OPM_AUTODIFFBLOCK_HEADER_INCLUDED

#include "disable_warning_pragmas.h"

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include "reenable_warning_pragmas.h"

#include <vector>
#include <cassert>

namespace Opm
{

    /// A class for forward-mode automatic differentiation with vector
    /// values and sparse jacobian matrices.
    ///
    /// The class contains a (column) vector of values and multiple
    /// sparse matrices representing its partial derivatives. Each
    /// such matrix has a number of rows equal to the number of rows
    /// in the value vector, and a number of columns equal to the
    /// number of discrete variables we want to compute the
    /// derivatives with respect to. The reason to have multiple such
    /// jacobians instead of just one is to allow simpler grouping of
    /// variables, making it easier to implement various
    /// preconditioning schemes. Only basic arithmetic operators are
    /// implemented for this class, reflecting our needs so far.
    ///
    /// The class is built on the Eigen library, using an Eigen array
    /// type to contain the values and Eigen sparse matrices for the
    /// jacobians. The overloaded operators are intended to behave in
    /// a similar way to Eigen arrays, meaning for example that the *
    /// operator is elementwise multiplication. The only exception is
    /// multiplication with a sparse matrix from the left, which is
    /// treated as an Eigen matrix operation.
    ///
    /// There are no public constructors, instead we use the Named
    /// Constructor pattern. In general, one needs to know which
    /// variables one wants to compute the derivatives with respect to
    /// before constructing an AutoDiffBlock. Some of the constructors
    /// require you to pass a block pattern. This should be a vector
    /// containing the number of columns you want for each jacobian
    /// matrix.
    ///
    /// For example: you want the derivatives with respect to three
    /// different variables p, r and s. Assuming that there are 10
    /// elements in p, and 20 in each of r and s, the block pattern is
    /// { 10, 20, 20 }. When creating the variables p, r and s in your
    /// program you have two options:
    ///     - Use the variable() constructor three times, passing the
    ///       index (0 for p, 1 for r and 2 for s), initial value of
    ///       each variable and the block pattern.
    ///     - Use the variables() constructor passing only the initial
    ///       values of each variable. The block pattern will be
    ///       inferred from the size of the initial value vectors.
    ///       This is usually the simplest option if you have multiple
    ///       variables. Note that this constructor returns a vector
    ///       of all three variables, so you need to use index access
    ///       (operator[]) to get the individual variables (that is p,
    ///       r and s).
    ///
    /// After this, the r variable for example will have a size() of
    /// 20 and three jacobian matrices. The first is a 20 by 10 zero
    /// matrix, the second is a 20 by 20 identity matrix, and the
    /// third is a 20 by 20 zero matrix.
    template <typename Scalar>
    class AutoDiffBlock
    {
    public:
        /// Underlying type for values.
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> V;
        /// Underlying type for jacobians.
        typedef Eigen::SparseMatrix<Scalar> M;

        /// Construct an empty AutoDiffBlock.
        static AutoDiffBlock null()
        {
            V val;
            std::vector<M> jac;
            return AutoDiffBlock(val, jac);
        }

        /// Create an AutoDiffBlock representing a constant.
        /// \param[in] val         values 
        static AutoDiffBlock constant(const V& val)
        {
            return AutoDiffBlock(val);
        }

        /// Create an AutoDiffBlock representing a constant.
        /// This variant requires specifying the block sizes used
        /// for the Jacobians even though the Jacobian matrices
        /// themselves will be zero.
        /// \param[in] val         values 
        /// \param[in] blocksizes  block pattern
        static AutoDiffBlock constant(const V& val, const std::vector<int>& blocksizes)
        {
            std::vector<M> jac;
            const int num_elem = val.size();
            const int num_blocks = blocksizes.size();
            // For constants, all jacobian blocks are zero.
            jac.resize(num_blocks);
            for (int i = 0; i < num_blocks; ++i) {
                jac[i] = M(num_elem, blocksizes[i]);
            }
            return AutoDiffBlock(val, jac);
        }

        /// Create an AutoDiffBlock representing a single variable block.
        /// \param[in] index       index of the variable you are constructing
        /// \param[in] val         values 
        /// \param[in] blocksizes  block pattern
        /// The resulting object will have size() equal to block_pattern[index].
        /// Its jacobians will all be zero, except for derivative()[index], which
        /// will be an identity matrix.
        static AutoDiffBlock variable(const int index, const V& val, const std::vector<int>& blocksizes)
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
            jac[index].reserve(Eigen::VectorXi::Constant(val.size(), 1));
            for (typename M::Index row = 0; row < val.size(); ++row) {
                jac[index].insert(row, row) = Scalar(1.0);
            }
            return AutoDiffBlock(val, jac);
        }

        /// Create an AutoDiffBlock by directly specifying values and jacobians.
        /// \param[in] val         values 
        /// \param[in] jac         vector of jacobians
        static AutoDiffBlock function(const V& val, const std::vector<M>& jac)
        {
            return AutoDiffBlock(val, jac);
        }

        /// Construct a set of primary variables, each initialized to
        /// a given vector.
        static std::vector<AutoDiffBlock> variables(const std::vector<V>& initial_values)
        {
            const int num_vars = initial_values.size();
            std::vector<int> bpat;
            for (int v = 0; v < num_vars; ++v) {
                bpat.push_back(initial_values[v].size());
            }
            std::vector<AutoDiffBlock> vars;
            for (int v = 0; v < num_vars; ++v) {
                vars.emplace_back(variable(v, initial_values[v], bpat));
            }
            return vars;
        }

        /// Elementwise operator +=
        AutoDiffBlock& operator+=(const AutoDiffBlock& rhs)
        {
            if (jac_.empty()) {
                jac_ = rhs.jac_;
            } else if (!rhs.jac_.empty()) {
                assert (numBlocks()    == rhs.numBlocks());
                assert (value().size() == rhs.value().size());

                const int num_blocks = numBlocks();
                for (int block = 0; block < num_blocks; ++block) {
                    assert(jac_[block].rows() == rhs.jac_[block].rows());
                    assert(jac_[block].cols() == rhs.jac_[block].cols());
                    jac_[block] += rhs.jac_[block];
                }
            }

            val_ += rhs.val_;

            return *this;
        }

        /// Elementwise operator -=
        AutoDiffBlock& operator-=(const AutoDiffBlock& rhs)
        {
            if (jac_.empty()) {
                const int num_blocks = rhs.numBlocks();
                jac_.resize(num_blocks);
                for (int block = 0; block < num_blocks; ++block) {
                    jac_[block] = -rhs.jac_[block];
                }
            } else if (!rhs.jac_.empty()) {
                assert (numBlocks()    == rhs.numBlocks());
                assert (value().size() == rhs.value().size());

                const int num_blocks = numBlocks();
                for (int block = 0; block < num_blocks; ++block) {
                    assert(jac_[block].rows() == rhs.jac_[block].rows());
                    assert(jac_[block].cols() == rhs.jac_[block].cols());
                    jac_[block] -= rhs.jac_[block];
                }
            }

            val_ -= rhs.val_;

            return *this;
        }

        /// Elementwise operator +
        AutoDiffBlock operator+(const AutoDiffBlock& rhs) const
        {
            if (jac_.empty() && rhs.jac_.empty()) {
                return constant(val_ + rhs.val_);
            }
            if (jac_.empty()) {
                return val_ + rhs;
            }
            if (rhs.jac_.empty()) {
                return *this + rhs.val_;
            }
            std::vector<M> jac = jac_;
            assert(numBlocks() == rhs.numBlocks());
            int num_blocks = numBlocks();
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac[block].rows() == rhs.jac_[block].rows());
                assert(jac[block].cols() == rhs.jac_[block].cols());
                jac[block] += rhs.jac_[block];
            }
            return function(val_ + rhs.val_, jac);
        }

        /// Elementwise operator -
        AutoDiffBlock operator-(const AutoDiffBlock& rhs) const
        {
            if (jac_.empty() && rhs.jac_.empty()) {
                return constant(val_ - rhs.val_);
            }
            if (jac_.empty()) {
                return val_ - rhs;
            }
            if (rhs.jac_.empty()) {
                return *this - rhs.val_;
            }
            std::vector<M> jac = jac_;
            assert(numBlocks() == rhs.numBlocks());
            int num_blocks = numBlocks();
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac[block].rows() == rhs.jac_[block].rows());
                assert(jac[block].cols() == rhs.jac_[block].cols());
                jac[block] -= rhs.jac_[block];
            }
            return function(val_ - rhs.val_, jac);
        }

        /// Elementwise operator *
        AutoDiffBlock operator*(const AutoDiffBlock& rhs) const
        {
            if (jac_.empty() && rhs.jac_.empty()) {
                return constant(val_ * rhs.val_);
            }
            if (jac_.empty()) {
                return val_ * rhs;
            }
            if (rhs.jac_.empty()) {
                return *this * rhs.val_;
            }
            int num_blocks = numBlocks();
            std::vector<M> jac(num_blocks);
            assert(numBlocks() == rhs.numBlocks());
            typedef Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D;
            D D1 = val_.matrix().asDiagonal();
            D D2 = rhs.val_.matrix().asDiagonal();
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac_[block].rows() == rhs.jac_[block].rows());
                assert(jac_[block].cols() == rhs.jac_[block].cols());
                jac[block] = D2*jac_[block] + D1*rhs.jac_[block];
            }
            return function(val_ * rhs.val_, jac);
        }

        /// Elementwise operator /
        AutoDiffBlock operator/(const AutoDiffBlock& rhs) const
        {
            if (jac_.empty() && rhs.jac_.empty()) {
                return constant(val_ / rhs.val_);
            }
            if (jac_.empty()) {
                return val_ / rhs;
            }
            if (rhs.jac_.empty()) {
                return *this / rhs.val_;
            }
            int num_blocks = numBlocks();
            std::vector<M> jac(num_blocks);
            assert(numBlocks() == rhs.numBlocks());
            typedef Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D;
            D D1 = val_.matrix().asDiagonal();
            D D2 = rhs.val_.matrix().asDiagonal();
            D D3 = (1.0/(rhs.val_*rhs.val_)).matrix().asDiagonal();
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac_[block].rows() == rhs.jac_[block].rows());
                assert(jac_[block].cols() == rhs.jac_[block].cols());
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

        /// Number of elements
        int size() const
        {
            return val_.size();
        }

        /// Number of Jacobian blocks.
        int numBlocks() const
        {
            return jac_.size();
        }

        /// Sizes (number of columns) of Jacobian blocks.
        std::vector<int> blockPattern() const
        {
            const int nb = numBlocks();
            std::vector<int> bp(nb);
            for (int block = 0; block < nb; ++block) {
                bp[block] = jac_[block].cols();
            }
            return bp;
        }

        /// Function value.
        const V& value() const
        {
            return val_;
        }

        /// Function derivatives.
        const std::vector<M>& derivative() const
        {
            return jac_;
        }

    private:
        AutoDiffBlock(const V& val)
        : val_(val)
        {
        }

        AutoDiffBlock(const V& val,
                      const std::vector<M>& jac)
            : val_(val), jac_(jac)
        {
#ifndef NDEBUG
            const int num_elem = val_.size();
            const int num_blocks = jac_.size();
            for (int block = 0; block < num_blocks; ++block) {
                assert(num_elem == jac_[block].rows());
            }
#endif
        }

        V val_;
        std::vector<M> jac_;
    };


    // ---------  Free functions and operators for AutoDiffBlock  ---------

    /// Stream output.
    template <class Ostream, typename Scalar>
    Ostream&
    operator<<(Ostream& os, const AutoDiffBlock<Scalar>& fw)
    {
        return fw.print(os);
    }


    /// Multiply with sparse matrix from the left.
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator*(const typename AutoDiffBlock<Scalar>::M& lhs,
                                    const AutoDiffBlock<Scalar>& rhs)
    {
        int num_blocks = rhs.numBlocks();
        std::vector<typename AutoDiffBlock<Scalar>::M> jac(num_blocks);
        assert(lhs.cols() == rhs.value().rows());
        for (int block = 0; block < num_blocks; ++block) {
            jac[block] = lhs*rhs.derivative()[block];
        }
        typename AutoDiffBlock<Scalar>::V val = lhs*rhs.value().matrix();
        return AutoDiffBlock<Scalar>::function(val, jac);
    }


    /// Elementwise multiplication with constant on the left.
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator*(const typename AutoDiffBlock<Scalar>::V& lhs,
                                    const AutoDiffBlock<Scalar>& rhs)
    {
        return AutoDiffBlock<Scalar>::constant(lhs, rhs.blockPattern()) * rhs;
    }


    /// Elementwise multiplication with constant on the right.
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator*(const AutoDiffBlock<Scalar>& lhs,
                                    const typename AutoDiffBlock<Scalar>::V& rhs)
    {
        return rhs * lhs; // Commutative operation.
    }


    /// Elementwise addition with constant on the left.
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator+(const typename AutoDiffBlock<Scalar>::V& lhs,
                                    const AutoDiffBlock<Scalar>& rhs)
    {
        return AutoDiffBlock<Scalar>::constant(lhs, rhs.blockPattern()) + rhs;
    }


    /// Elementwise addition with constant on the right.
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator+(const AutoDiffBlock<Scalar>& lhs,
                                    const typename AutoDiffBlock<Scalar>::V& rhs)
    {
        return rhs + lhs; // Commutative operation.
    }


    /// Elementwise subtraction with constant on the left.
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator-(const typename AutoDiffBlock<Scalar>::V& lhs,
                                    const AutoDiffBlock<Scalar>& rhs)
    {
        return AutoDiffBlock<Scalar>::constant(lhs, rhs.blockPattern()) - rhs;
    }


    /// Elementwise subtraction with constant on the right.
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator-(const AutoDiffBlock<Scalar>& lhs,
                                    const typename AutoDiffBlock<Scalar>::V& rhs)
    {
        return lhs - AutoDiffBlock<Scalar>::constant(rhs, lhs.blockPattern());
    }


    /// Elementwise division with constant on the left.
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator/(const typename AutoDiffBlock<Scalar>::V& lhs,
                                    const AutoDiffBlock<Scalar>& rhs)
    {
        return AutoDiffBlock<Scalar>::constant(lhs, rhs.blockPattern()) / rhs;
    }


    /// Elementwise division with constant on the right.
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator/(const AutoDiffBlock<Scalar>& lhs,
                                    const typename AutoDiffBlock<Scalar>::V& rhs)
    {
        return lhs / AutoDiffBlock<Scalar>::constant(rhs, lhs.blockPattern());
    }


    /**
     * @brief Operator for multiplication with a scalar on the right-hand side
     *
     * @param lhs The left-hand side AD forward block
     * @param rhs The scalar to multiply with
     * @return The product
     */
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator*(const AutoDiffBlock<Scalar>& lhs,
                                    const Scalar& rhs)
    {
        std::vector< typename AutoDiffBlock<Scalar>::M > jac;
        jac.reserve( lhs.numBlocks() );
        for (int block=0; block<lhs.numBlocks(); block++) {
            jac.emplace_back( lhs.derivative()[block] * rhs );
        }
        return AutoDiffBlock<Scalar>::function( lhs.value() * rhs, jac );
    }


    /**
     * @brief Operator for multiplication with a scalar on the left-hand side
     *
     * @param lhs The scalar to multiply with
     * @param rhs The right-hand side AD forward block
     * @return The product
     */
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator*(const Scalar& lhs,
                                    const AutoDiffBlock<Scalar>& rhs)
    {
        return rhs * lhs; // Commutative operation.
    }


} // namespace Opm



#endif // OPM_AUTODIFFBLOCK_HEADER_INCLUDED
