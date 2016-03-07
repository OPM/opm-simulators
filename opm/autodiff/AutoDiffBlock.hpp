/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2016 IRIS AS

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

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <opm/autodiff/fastSparseOperations.hpp>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/autodiff/AutoDiffMatrix.hpp>


#include <utility>
#include <vector>
#include <cassert>
#include <iostream>

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
        typedef AutoDiffMatrix M;

        /// Construct an empty AutoDiffBlock.
        static AutoDiffBlock null()
        {
            return AutoDiffBlock(V(), {});
        }

        /// Create an AutoDiffBlock representing a constant.
        /// \param[in] val         values
        static AutoDiffBlock constant(V&& val)
        {
            return AutoDiffBlock(std::move(val));
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
            V val_copy(val);
            return AutoDiffBlock(std::move(val_copy), std::move(jac));
        }

        /// Create an AutoDiffBlock representing a single variable block.
        /// \param[in] index       index of the variable you are constructing
        /// \param[in] val         values
        /// \param[in] blocksizes  block pattern
        /// The resulting object will have size() equal to block_pattern[index].
        /// Its jacobians will all be zero, except for derivative()[index], which
        /// will be an identity matrix.
        static AutoDiffBlock variable(const int index, V&& val, const std::vector<int>& blocksizes)
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
            jac[index] = M::createIdentity(val.size());
            return AutoDiffBlock(std::move(val), std::move(jac));
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
            V value = val;
            return variable(index, std::move(value), blocksizes);
        }

        /// Create an AutoDiffBlock by directly specifying values and jacobians.
        /// This version of function() moves its arguments and is therefore
        /// quite efficient, but leaves the argument variables empty (but valid).
        /// \param[in] val         values
        /// \param[in] jac         vector of jacobians
        static AutoDiffBlock function(V&& val, std::vector<M>&& jac)
        {
            return AutoDiffBlock(std::move(val), std::move(jac));
        }

        /// Create an AutoDiffBlock by directly specifying values and jacobians.
        /// This version of function() copies its arguments and is therefore
        /// less efficient than the other (moving) overload.
        /// \param[in] val         values
        /// \param[in] jac         vector of jacobians
        static AutoDiffBlock function(const V& val, const std::vector<M>& jac)
        {
            V val_copy(val);
            std::vector<M> jac_copy(jac);
            return AutoDiffBlock(std::move(val_copy), std::move(jac_copy));
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
#pragma omp parallel for schedule(static)
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
#pragma omp parallel for schedule(static)
                for (int block = 0; block < num_blocks; ++block) {
                    jac_[block] = rhs.jac_[block] * (-1.0);
                }
            } else if (!rhs.jac_.empty()) {
                assert (numBlocks()    == rhs.numBlocks());
                assert (value().size() == rhs.value().size());

                const int num_blocks = numBlocks();
#pragma omp parallel for schedule(static)
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
#pragma omp parallel for schedule(static)
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac[block].rows() == rhs.jac_[block].rows());
                assert(jac[block].cols() == rhs.jac_[block].cols());
                jac[block] += rhs.jac_[block];
            }
            return function(val_ + rhs.val_, std::move(jac));
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
#pragma omp parallel for schedule(static)
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac[block].rows() == rhs.jac_[block].rows());
                assert(jac[block].cols() == rhs.jac_[block].cols());
                jac[block] -= rhs.jac_[block];
            }
            return function(val_ - rhs.val_, std::move(jac));
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
            M D1(val_.matrix().asDiagonal());
            M D2(rhs.val_.matrix().asDiagonal());
#pragma omp parallel for schedule(dynamic)
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac_[block].rows() == rhs.jac_[block].rows());
                assert(jac_[block].cols() == rhs.jac_[block].cols());
                if( jac_[block].nonZeros() == 0 && rhs.jac_[block].nonZeros() == 0 ) {
                    jac[block] = M( D2.rows(), jac_[block].cols() );
                }
                else if( jac_[block].nonZeros() == 0 )
                    jac[block] = D1*rhs.jac_[block];
                else if ( rhs.jac_[block].nonZeros() == 0 ) {
                    jac[block] = D2*jac_[block];
                }
                else {
                    jac[block]  = D2*jac_[block];
                    jac[block] += D1*rhs.jac_[block];
                }
            }
            return function(val_ * rhs.val_, std::move(jac));
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
            M D1(val_.matrix().asDiagonal());
            M D2(rhs.val_.matrix().asDiagonal());
            M D3((1.0/(rhs.val_*rhs.val_)).matrix().asDiagonal());
#pragma omp parallel for schedule(dynamic)
            for (int block = 0; block < num_blocks; ++block) {
                assert(jac_[block].rows() == rhs.jac_[block].rows());
                assert(jac_[block].cols() == rhs.jac_[block].cols());
                if( jac_[block].nonZeros() == 0 && rhs.jac_[block].nonZeros() == 0 ) {
                    jac[block] = M( D3.rows(), jac_[block].cols() );
                }
                else if( jac_[block].nonZeros() == 0 ) {
                    jac[block] =  D3 * ( D1*rhs.jac_[block]) * (-1.0);
                }
                else if ( rhs.jac_[block].nonZeros() == 0 )
                {
                    jac[block] = D3 * (D2*jac_[block]);
                }
                else {
                    jac[block] = D3 * (D2*jac_[block] + (D1*rhs.jac_[block]*(-1.0)));
                }
            }
            return function(val_ / rhs.val_, std::move(jac));
        }

        /// I/O.
        template <class Ostream>
        Ostream&
        print(Ostream& os) const
        {
            int num_blocks = jac_.size();
            os << "Value =\n" << val_ << "\n\nJacobian =\n";
            for (int i = 0; i < num_blocks; ++i) {
                Eigen::SparseMatrix<double> m;
                jac_[i].toSparse(m);
                os << "Sub Jacobian #" << i << '\n' << m << "\n";
            }
            return os;
        }

        /// Efficient swap function.
        void swap(AutoDiffBlock& other)
        {
            val_.swap(other.val_);
            jac_.swap(other.jac_);
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

        AutoDiffBlock(V&& val)
            : val_(std::move(val))
        {
        }

        AutoDiffBlock(V&& val, std::vector<M>&& jac)
            : val_(std::move(val)), jac_(std::move(jac))
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


    /// Multiply with AutoDiffMatrix from the left.
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator*(const typename AutoDiffBlock<Scalar>::M& lhs,
                                    const AutoDiffBlock<Scalar>& rhs)
    {
        int num_blocks = rhs.numBlocks();
        std::vector<typename AutoDiffBlock<Scalar>::M> jac(num_blocks);
        assert(lhs.cols() == rhs.value().rows());
#pragma omp parallel for schedule(dynamic)
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(lhs, rhs.derivative()[block], jac[block]);
        }
        typename AutoDiffBlock<Scalar>::V val = lhs*rhs.value().matrix();
        return AutoDiffBlock<Scalar>::function(std::move(val), std::move(jac));
    }


    /// Multiply with Eigen sparse matrix from the left.
    template <typename Scalar>
    AutoDiffBlock<Scalar> operator*(const Eigen::SparseMatrix<Scalar>& lhs,
                                    const AutoDiffBlock<Scalar>& rhs)
    {
        int num_blocks = rhs.numBlocks();
        std::vector<typename AutoDiffBlock<Scalar>::M> jac(num_blocks);
        assert(lhs.cols() == rhs.value().rows());
        for (int block = 0; block < num_blocks; ++block) {
            fastSparseProduct(lhs, rhs.derivative()[block], jac[block]);
        }
        typename AutoDiffBlock<Scalar>::V val = lhs*rhs.value().matrix();
        return AutoDiffBlock<Scalar>::function(std::move(val), std::move(jac));
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
        return AutoDiffBlock<Scalar>::function( lhs.value() * rhs, std::move(jac) );
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

    /**
     * @brief Computes the value of base raised to the power of exp elementwise
     *
     * @param base The AD forward block
     * @param exp  array of exponents
     * @return The value of base raised to the power of exp elementwise
     */
    template <typename Scalar>
    AutoDiffBlock<Scalar> pow(const AutoDiffBlock<Scalar>& base,
                              const typename AutoDiffBlock<Scalar>::V& exp)
    {
        const int num_elem = base.value().size();
        typename AutoDiffBlock<Scalar>::V val (num_elem);
        typename AutoDiffBlock<Scalar>::V derivative = exp;
        assert(exp.size() == num_elem);
        for (int i = 0; i < num_elem; ++i) {
            val[i] = std::pow(base.value()[i], exp[i]);
            derivative[i] *= std::pow(base.value()[i], exp[i] - 1.0);
        }
        const typename AutoDiffBlock<Scalar>::M derivative_diag(derivative.matrix().asDiagonal());

        std::vector< typename AutoDiffBlock<Scalar>::M > jac (base.numBlocks());
        for (int block = 0; block < base.numBlocks(); block++) {
             fastSparseProduct(derivative_diag, base.derivative()[block], jac[block]);
        }

        return AutoDiffBlock<Scalar>::function( std::move(val), std::move(jac) );
    }


    /**
     * @brief Computes the value of base raised to the power of exp
     *
     * @param base The AD forward block
     * @param exp  exponent
     * @return The value of base raised to the power of exp
     */
    template <typename Scalar>
    AutoDiffBlock<Scalar> pow(const AutoDiffBlock<Scalar>& base,
                              const double exp)
    {
        const typename AutoDiffBlock<Scalar>::V val = base.value().pow(exp);
        const typename AutoDiffBlock<Scalar>::V derivative = exp * base.value().pow(exp - 1.0);
        const typename AutoDiffBlock<Scalar>::M derivative_diag(derivative.matrix().asDiagonal());

        std::vector< typename AutoDiffBlock<Scalar>::M > jac (base.numBlocks());
        for (int block = 0; block < base.numBlocks(); block++) {
             fastSparseProduct(derivative_diag, base.derivative()[block], jac[block]);
        }

        return AutoDiffBlock<Scalar>::function( std::move(val), std::move(jac) );
    }

    /**
     * @brief Computes the value of base raised to the power of exp
     *
     * @param base The base AD forward block
     * @param exp  The exponent AD forward block
     * @return The value of base raised to the power of exp
     */    template <typename Scalar>
    AutoDiffBlock<Scalar> pow(const AutoDiffBlock<Scalar>& base,
                                    const AutoDiffBlock<Scalar>& exp)
    {
        const int num_elem = base.value().size();
        assert(exp.value().size() == num_elem);
        typename AutoDiffBlock<Scalar>::V val (num_elem);
        for (int i = 0; i < num_elem; ++i) {
            val[i] = std::pow(base.value()[i], exp.value()[i]);
        }

        // (f^g)' = f^g * ln(f) * g' + g * f^(g-1) * f'
        typename AutoDiffBlock<Scalar>::V der1 = val;
        for (int i = 0; i < num_elem; ++i) {
            der1[i] *= std::log(base.value()[i]);
        }
        std::vector< typename AutoDiffBlock<Scalar>::M > jac1 (base.numBlocks());
        const typename AutoDiffBlock<Scalar>::M der1_diag(der1.matrix().asDiagonal());
        for (int block = 0; block < base.numBlocks(); block++) {
             fastSparseProduct(der1_diag, exp.derivative()[block], jac1[block]);
        }
        typename AutoDiffBlock<Scalar>::V der2 = exp.value();
        for (int i = 0; i < num_elem; ++i) {
            der2[i] *= std::pow(base.value()[i], exp.value()[i] - 1.0);
        }
        std::vector< typename AutoDiffBlock<Scalar>::M > jac2 (base.numBlocks());
        const typename AutoDiffBlock<Scalar>::M der2_diag(der2.matrix().asDiagonal());
        for (int block = 0; block < base.numBlocks(); block++) {
             fastSparseProduct(der2_diag, base.derivative()[block], jac2[block]);
        }
        std::vector< typename AutoDiffBlock<Scalar>::M > jac (base.numBlocks());
        for (int block = 0; block < base.numBlocks(); block++) {
             jac[block] = jac1[block] + jac2[block];
        }

        return AutoDiffBlock<Scalar>::function(std::move(val), std::move(jac));
    }



} // namespace Opm



#endif // OPM_AUTODIFFBLOCK_HEADER_INCLUDED
