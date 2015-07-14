/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Statoil AS.

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

#ifndef OPM_AUTODIFFDENSEBLOCK_HEADER_INCLUDED
#define OPM_AUTODIFFDENSEBLOCK_HEADER_INCLUDED

#include <opm/core/utility/ErrorMacros.hpp>

#include <opm/core/utility/platform_dependent/disable_warnings.h>

#include <Eigen/Eigen>

#include <opm/core/utility/platform_dependent/reenable_warnings.h>

#include <utility>
#include <vector>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <array>

namespace Opm
{

    /// A class for forward-mode automatic differentiation with array
    /// values and a small, dense set of derivatives.
    ///
    /// The class contains a (column) array of values and a rectangular
    /// dense array representing its partial derivatives. Each
    /// such array has a number of rows equal to the number of rows
    /// in the value array, and a number of columns equal to the
    /// number of discrete variables we want to compute the
    /// derivatives with respect to. Only basic arithmetic operators.
    ///
    /// The class is built on the Eigen library, using Eigen array
    /// types. The overloaded operators are intended to behave in
    /// a similar way to Eigen arrays, meaning for example that the *
    /// operator is elementwise multiplication.
    ///
    /// There are no public constructors, instead we use the Named
    /// Constructor pattern. In general, one needs to know which
    /// variables one wants to compute the derivatives with respect to
    /// before constructing an AutoDiffDenseBlock.
    ///
    /// For example: you want the derivatives with respect to three
    /// different variables p, r and s, each of which has 20 elements.
    /// When creating the variables p, r and s in your program you
    /// have two options:
    ///     - Use the variable() constructor three times, passing the
    ///       index (0 for p, 1 for r and 2 for s) and initial values
    ///       for each variable.
    ///     - Use the variables() constructor passing only the initial
    ///       values of each variable. This is usually the simplest
    ///       option if you have multiple variables. Note that this
    ///       constructor returns a vector of all three variables, so
    ///       you need to use index access (operator[]) to get the
    ///       individual variables (that is p, r and s).
    ///
    /// After this, the r variable for example will have a size() of
    /// 20, its value() will have 20 values equal to the initial
    /// values passed in during construction, and it derivative() will
    /// be a 20x3 array with zeros in columns 0 and 2, and ones in
    /// column 1.
    template <typename Scalar, int NumDerivs>
    class AutoDiffDenseBlock
    {
    public:
        /// Underlying type for values.
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> Value;
        /// Underlying type for derivatives.
        typedef Eigen::Array<Scalar, Eigen::Dynamic, NumDerivs> Derivative;

        /// Construct an empty AutoDiffDenseBlock.
        static AutoDiffDenseBlock null()
        {
            return AutoDiffDenseBlock(Value());
        }

        /// Create an AutoDiffDenseBlock representing a constant.
        /// \param[in] val         values
        static AutoDiffDenseBlock constant(Value&& val)
        {
            return AutoDiffDenseBlock(std::move(val));
        }

        /// Create an AutoDiffDenseBlock representing a constant.
        /// \param[in] val         values
        static AutoDiffDenseBlock constant(const Value& val)
        {
            return AutoDiffDenseBlock(val);
        }

        /// Create an AutoDiffDenseBlock representing a single variable block.
        /// \param[in] index       index of the variable you are constructing
        /// \param[in] val         values
        /// The resulting object will have size() val.size().
        /// Its derivative() will all be zero, except for column index, which
        /// will be all ones.
        static AutoDiffDenseBlock variable(const int index, Value&& val)
        {
            if (index >= NumDerivs) {
                OPM_THROW(std::runtime_error, "Cannot create variable with index "
                          << index << " since NumDerivs = " << NumDerivs);
            }
            Derivative jac = Derivative::Zero(val.rows(), 3);
            jac.col(index) = Value::Ones(val.rows());
            AutoDiffDenseBlock retval(std::move(val), std::move(jac));
            retval.zero_jac_.fill(true);
            retval.zero_jac_[index] = false;
            return std::move(retval);
        }

        /// Create an AutoDiffDenseBlock representing a single variable block.
        /// \param[in] index       index of the variable you are constructing
        /// \param[in] val         values
        /// The resulting object will have size() val.size().
        /// Its derivative() will all be zero, except for column index, which
        /// will be all ones.
        static AutoDiffDenseBlock variable(const int index, const Value& val)
        {
            Value val_copy = val;
            return variable(index, std::move(val_copy));
        }

        /// Create an AutoDiffDenseBlock by directly specifying values and jacobians.
        /// This version of function() moves its arguments and is therefore
        /// quite efficient, but leaves the argument variables empty (but valid).
        /// \param[in] val         values
        /// \param[in] jac         derivatives
        static AutoDiffDenseBlock function(Value&& val, Derivative&& jac)
        {
            return AutoDiffDenseBlock(std::move(val), std::move(jac));
        }

        /// Create an AutoDiffDenseBlock by directly specifying values and jacobians.
        /// This version of function() copies its arguments and is therefore
        /// less efficient than the other (moving) overload.
        /// \param[in] val         values
        /// \param[in] jac         derivatives
        static AutoDiffDenseBlock function(const Value& val, const Derivative& jac)
        {
            Value val_copy(val);
            Derivative jac_copy(jac);
            return AutoDiffDenseBlock(std::move(val_copy), std::move(jac_copy));
        }

        /// Construct a set of primary variables, each initialized to
        /// a given vector.
        static std::vector<AutoDiffDenseBlock> variables(const std::vector<Value>& initial_values)
        {
            const int num_vars = initial_values.size();
            std::vector<AutoDiffDenseBlock> vars;
            for (int v = 0; v < num_vars; ++v) {
                vars.emplace_back(variable(v, initial_values[v]));
            }
            return vars;
        }

        /// Elementwise operator +=
        AutoDiffDenseBlock& operator+=(const AutoDiffDenseBlock& rhs)
        {
            val_ += rhs.val_;
            jac_ += rhs.jac_;
            zeroDerivIfBothZero(rhs);
            return *this;
        }

        /// Elementwise operator -=
        AutoDiffDenseBlock& operator-=(const AutoDiffDenseBlock& rhs)
        {
            val_ -= rhs.val_;
            jac_ -= rhs.jac_;
            zeroDerivIfBothZero(rhs);
            return *this;
        }

        /// Elementwise operator +
        AutoDiffDenseBlock operator+(const AutoDiffDenseBlock& rhs) const
        {
            return function(val_ + rhs.val_, jac_ + rhs.jac_);
        }

        /// Elementwise operator -
        AutoDiffDenseBlock operator-(const AutoDiffDenseBlock& rhs) const
        {
            return function(val_ - rhs.val_, jac_ - rhs.jac_);
        }

        /// Elementwise operator *
        AutoDiffDenseBlock operator*(const AutoDiffDenseBlock& rhs) const
        {
            Derivative jac = jac_.colwise() * rhs.val_  +  rhs.jac_.colwise() * val_;
            return function(val_ * rhs.val_, std::move(jac));
        }

        /// Elementwise operator /
        AutoDiffDenseBlock operator/(const AutoDiffDenseBlock& rhs) const
        {
            Derivative jac = (jac_.colwise() * rhs.val_ - rhs.jac_.colwise() * val_).colwise() / (rhs.val_*rhs.val_);
            return function(val_ / rhs.val_, std::move(jac));
        }

        /// I/O.
        template <class Ostream>
        Ostream&
        print(Ostream& os) const
        {
            os << "Value =\n" << val_ << "\n\nJacobian =\n" << jac_ << '\n';
            return os;
        }

        /// Efficient swap function.
        void swap(AutoDiffDenseBlock& other)
        {
            val_.swap(other.val_);
            jac_.swap(other.jac_);
            zero_jac_.swap(other.zero_jac_);
        }

        /// Number of elements
        int size() const
        {
            return val_.size();
        }

        /// Function value.
        const Value& value() const
        {
            return val_;
        }

        /// Function derivatives.
        const Derivative& derivative() const
        {
            return jac_;
        }

        bool hasZeroDerivative(int variable) const
        {
            return zero_jac_[variable];
        }

        bool isConstant() const
        {
            for (int i = 0; i < NumDerivs; ++i) {
                if (!zero_jac_[i]) {
                    return false;
                }
            }
            return true;
        }

    private:
        AutoDiffDenseBlock(const Value& val)
        : val_(val), jac_(Derivative::Zero(val.rows(), NumDerivs))
        {
            zero_jac_.fill(true);
        }

        AutoDiffDenseBlock(Value&& val)
        : val_(std::move(val)), jac_(Derivative::Zero(val.rows(), NumDerivs))
        {
            zero_jac_.fill(true);
        }

        AutoDiffDenseBlock(Value&& val, Derivative&& jac)
            : val_(std::move(val)), jac_(std::move(jac))
        {
            zero_jac_.fill(false);
            assert(val_.rows() == jac_.rows());
        }

        void zeroDerivIfBothZero(const AutoDiffDenseBlock& other)
        {
            for (int i = 0; i < NumDerivs; ++i) {
                zero_jac_[i] &= other.zero_jac_[i];
            }
        }

        Value val_;
        Derivative jac_;
        std::array<bool, NumDerivs> zero_jac_;
    };


    // ---------  Free functions and operators for AutoDiffDenseBlock  ---------

    /// Stream output.
    template <class Ostream, typename Scalar, int NumDerivs>
    Ostream&
    operator<<(Ostream& os, const AutoDiffDenseBlock<Scalar, NumDerivs>& fw)
    {
        return fw.print(os);
    }


    /// Elementwise multiplication with constant on the left.
    template <typename Scalar, int NumDerivs>
    AutoDiffDenseBlock<Scalar, NumDerivs> operator*(const typename AutoDiffDenseBlock<Scalar, NumDerivs>::Value& lhs,
                                                    const AutoDiffDenseBlock<Scalar, NumDerivs>& rhs)
    {
        return AutoDiffDenseBlock<Scalar, NumDerivs>::constant(lhs) * rhs;
    }


    /// Elementwise multiplication with constant on the right.
    template <typename Scalar, int NumDerivs>
    AutoDiffDenseBlock<Scalar, NumDerivs> operator*(const AutoDiffDenseBlock<Scalar, NumDerivs>& lhs,
                                                    const typename AutoDiffDenseBlock<Scalar, NumDerivs>::Value& rhs)
    {
        return rhs * lhs; // Commutative operation.
    }


    /// Elementwise addition with constant on the left.
    template <typename Scalar, int NumDerivs>
    AutoDiffDenseBlock<Scalar, NumDerivs> operator+(const typename AutoDiffDenseBlock<Scalar, NumDerivs>::Value& lhs,
                                                    const AutoDiffDenseBlock<Scalar, NumDerivs>& rhs)
    {
        return AutoDiffDenseBlock<Scalar, NumDerivs>::constant(lhs) + rhs;
    }


    /// Elementwise addition with constant on the right.
    template <typename Scalar, int NumDerivs>
    AutoDiffDenseBlock<Scalar, NumDerivs> operator+(const AutoDiffDenseBlock<Scalar, NumDerivs>& lhs,
                                                    const typename AutoDiffDenseBlock<Scalar, NumDerivs>::Value& rhs)
    {
        return rhs + lhs; // Commutative operation.
    }


    /// Elementwise subtraction with constant on the left.
    template <typename Scalar, int NumDerivs>
    AutoDiffDenseBlock<Scalar, NumDerivs> operator-(const typename AutoDiffDenseBlock<Scalar, NumDerivs>::Value& lhs,
                                                    const AutoDiffDenseBlock<Scalar, NumDerivs>& rhs)
    {
        return AutoDiffDenseBlock<Scalar, NumDerivs>::constant(lhs, rhs.blockPattern()) - rhs;
    }


    /// Elementwise subtraction with constant on the right.
    template <typename Scalar, int NumDerivs>
    AutoDiffDenseBlock<Scalar, NumDerivs> operator-(const AutoDiffDenseBlock<Scalar, NumDerivs>& lhs,
                                                    const typename AutoDiffDenseBlock<Scalar, NumDerivs>::Value& rhs)
    {
        return lhs - AutoDiffDenseBlock<Scalar, NumDerivs>::constant(rhs);
    }


    /// Elementwise division with constant on the left.
    template <typename Scalar, int NumDerivs>
    AutoDiffDenseBlock<Scalar, NumDerivs> operator/(const typename AutoDiffDenseBlock<Scalar, NumDerivs>::Value& lhs,
                                                    const AutoDiffDenseBlock<Scalar, NumDerivs>& rhs)
    {
        return AutoDiffDenseBlock<Scalar, NumDerivs>::constant(lhs) / rhs;
    }


    /// Elementwise division with constant on the right.
    template <typename Scalar, int NumDerivs>
    AutoDiffDenseBlock<Scalar, NumDerivs> operator/(const AutoDiffDenseBlock<Scalar, NumDerivs>& lhs,
                                                    const typename AutoDiffDenseBlock<Scalar, NumDerivs>::Value& rhs)
    {
        return lhs / AutoDiffDenseBlock<Scalar, NumDerivs>::constant(rhs);
    }


    /**
     * @brief Operator for multiplication with a scalar on the right-hand side
     *
     * @param lhs The left-hand side AD forward block
     * @param rhs The scalar to multiply with
     * @return The product
     */
    template <typename Scalar, int NumDerivs>
    AutoDiffDenseBlock<Scalar, NumDerivs> operator*(const AutoDiffDenseBlock<Scalar, NumDerivs>& lhs,
                                                    const Scalar& rhs)
    {
        return AutoDiffDenseBlock<Scalar, NumDerivs>::function( lhs.value() * rhs, lhs.derivative() * rhs );
    }


    /**
     * @brief Operator for multiplication with a scalar on the left-hand side
     *
     * @param lhs The scalar to multiply with
     * @param rhs The right-hand side AD forward block
     * @return The product
     */
    template <typename Scalar, int NumDerivs>
    AutoDiffDenseBlock<Scalar, NumDerivs> operator*(const Scalar& lhs,
                                                    const AutoDiffDenseBlock<Scalar, NumDerivs>& rhs)
    {
        return rhs * lhs; // Commutative operation.
    }


} // namespace Opm



#endif // OPM_AUTODIFFDENSEBLOCK_HEADER_INCLUDED
