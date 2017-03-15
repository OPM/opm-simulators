// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief Representation of an evaluation of a function and its derivatives w.r.t. a set
 *        of variables in the localized OPM automatic differentiation (AD) framework.
 */

#ifndef OPM_LOCAL_AD_EVALUATION_HPP
#define OPM_LOCAL_AD_EVALUATION_HPP

#include "Math.hpp"

#include <opm/common/Valgrind.hpp>

#include <dune/common/version.hh>

#include <array>
#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <algorithm>

namespace Opm {
namespace DenseAd {
template <class ValueT, int numVars>
class Evaluation;

/*!
 * \brief Provides operations on Evaluations, so that these can be specialized without
 *        having to copy-and-paste the whole Evaluation class.
 *
 * \internal
 */
template <class ValueT, int numVars>
struct EvaluationOps
{
private:
    typedef Evaluation<ValueT, numVars> Eval;

    static constexpr int length_ = Eval::length_;
    static constexpr int valuepos_ = Eval::valuepos_;
    static constexpr int dstart_ = Eval::dstart_;
    static constexpr int dend_ = Eval::dend_ ;

public:
    typedef std::array<ValueT, length_> DataVector;

    static inline void clearDerivatives(Eval& a)
    {
        for (int i = dstart_; i < dend_; ++i)
            a.data_[i] = 0.0;
    }

    static inline void assign(Eval& a, const Eval& b)
    {
        a.data_ = b.data_;
    }

    static inline void assignNegative(Eval& a, const Eval& b)
    {
        for (int idx = 0; idx < length_; ++idx)
            a.data_[idx] = - b.data_[idx];
    }

    static inline Eval& addEq(Eval& a, const Eval& b)
    {
        for (int i = 0; i < length_; ++i)
            a.data_[i] += b.data_[i];

        return a;
    }

    static inline Eval& subEq(Eval& a, const Eval& b)
    {
        for (int i = 0; i < length_; ++i)
            a.data_[i] -= b.data_[i];

        return a;
    }

    static inline Eval& mulEq(Eval& a, const Eval& b)
    {
        // while the values are multiplied, the derivatives follow the product rule,
        // i.e., (u*v)' = (v'u + u'v).
        const ValueT u = a.value();
        const ValueT v = b.value();

        a.data_[valuepos_] *= v ;
        for (int idx = dstart_; idx < dend_; ++idx)
            a.data_[idx] = a.data_[idx] * v + b.data_[idx] * u;

        return a;
    }

    template <class RhsValueType>
    static inline Eval& scalarMulEq(Eval& a, const RhsValueType& other)
    {
        // values and derivatives are multiplied
        for (int idx = 0 ; idx < length_; ++ idx)
            a.data_[idx] *= other;

        return a;
    }

    static inline Eval& divEq(Eval& a, const Eval& b)
    {
        // values are divided, derivatives follow the rule for division, i.e., (u/v)' = (v'u - u'v)/v^2.
        const ValueT v_vv = 1.0 / b.value();
        const ValueT u_vv = a.value() * v_vv * v_vv;

        a.data_[valuepos_] *= v_vv;
        for (int idx = dstart_; idx < dend_; ++idx)
            a.data_[idx] = a.data_[idx] * v_vv - b.data_[idx] * u_vv;

        return a;
    }

    template <class RhsValueType>
    static inline Eval& scalarDivEq(Eval& a, const RhsValueType& other)
    {
        // values and derivatives are divided
        for (int idx = 0 ; idx <= length_; ++ idx)
            a.data_[idx] /= other;

        return a;
    }

    // a/b with 'a' being a scalar and 'b' an Evaluation
    template <class RhsValueType>
    static inline Eval divide(const RhsValueType& a, const Eval& b )
    {
        Eval result;

        result.setValue( a/b.value() );

        const ValueT df_dg = - result.value()/b.value();
        for (int idx = dstart_; idx < dend_; ++idx)
            result.data_[idx] = df_dg*b.data_[idx];

        return result;
    }
};


/*!
 * \brief Represents a function evaluation and its derivatives w.r.t. a fixed set of
 *        variables.
 */
template <class ValueT, int numVars>
class Evaluation
{
    typedef EvaluationOps<ValueT, numVars> Ops;

    friend Ops;

public:
    //! field type
    typedef ValueT ValueType;

    //! number of derivatives
    static constexpr int size = numVars;

protected:
    //! length of internal data vector
    static constexpr int length_ = numVars + 1 ;

    //! position index for value
    static constexpr int valuepos_ = 0;
    //! start index for derivatives
    static constexpr int dstart_   = 1;
    //! end+1 index for derivatives
    static constexpr int dend_     = length_ ;
public:

    //! default constructor
    Evaluation() : data_()
    {}

    //! copy other function evaluation
    Evaluation(const Evaluation& other)
    {
        Ops::assign(*this, other);
    }

    // create an evaluation which represents a constant function
    //
    // i.e., f(x) = c. this implies an evaluation with the given value and all
    // derivatives being zero.
    template <class RhsValueType>
    Evaluation(const RhsValueType& c)
    {
        setValue( c );
        clearDerivatives();
        Valgrind::CheckDefined( data_ );
    }

    // create an evaluation which represents a constant function
    //
    // i.e., f(x) = c. this implies an evaluation with the given value and all
    // derivatives being zero.
    template <class RhsValueType>
    Evaluation(const RhsValueType& c, int varPos)
    {
        setValue( c );
        clearDerivatives();
        // The variable position must be in represented by the given variable descriptor
        assert(0 <= varPos && varPos < numVars);

        data_[varPos + dstart_] = 1.0;
        Valgrind::CheckDefined(data_);
    }

    // set all derivatives to zero
    void clearDerivatives()
    {
        Ops::clearDerivatives(*this);
    }

    // create a function evaluation for a "naked" depending variable (i.e., f(x) = x)
    template <class RhsValueType>
    static Evaluation createVariable(const RhsValueType& value, int varPos)
    {
        // copy function value and set all derivatives to 0, except for the variable
        // which is represented by the value (which is set to 1.0)
        return Evaluation( value, varPos );
    }

    // "evaluate" a constant function (i.e. a function that does not depend on the set of
    // relevant variables, f(x) = c).
    template <class RhsValueType>
    static Evaluation createConstant(const RhsValueType& value)
    {
        return Evaluation( value );
    }

    // print the value and the derivatives of the function evaluation
    void print(std::ostream& os = std::cout) const
    {
        // print value
        os << "v: " << value() << " / d:";
        // print derivatives
        for (int varIdx = 0; varIdx < numVars; ++varIdx)
            os << " " << derivative(varIdx);
    }

    // copy all derivatives from other
    void copyDerivatives(const Evaluation& other)
    {
        for (int varIdx = dstart_; varIdx < dend_; ++varIdx)
            data_[ varIdx ] = other.data_[ varIdx ];
    }


    // add value and derivatives from other to this values and derivatives
    Evaluation& operator+=(const Evaluation& other)
    { return Ops::addEq(*this, other); }

    // add value from other to this values
    template <class RhsValueType>
    Evaluation& operator+=(const RhsValueType& other)
    {
        // value is added, derivatives stay the same
        data_[valuepos_] += other;
        return *this;
    }

    // subtract other's value and derivatives from this values
    Evaluation& operator-=(const Evaluation& other)
    { return Ops::subEq(*this, other); }

    // subtract other's value from this values
    template <class RhsValueType>
    Evaluation& operator-=(const RhsValueType& other)
    {
        // for constants, values are subtracted, derivatives stay the same
        data_[ valuepos_ ] -= other;

        return *this;
    }

    // multiply values and apply chain rule to derivatives: (u*v)' = (v'u + u'v)
    Evaluation& operator*=(const Evaluation& other)
    { return Ops::mulEq(*this, other); }

    // m(u*v)' = (v'u + u'v)
    template <class RhsValueType>
    Evaluation& operator*=(const RhsValueType& other)
    { return Ops::scalarMulEq(*this, other); }

    // m(u*v)' = (v'u + u'v)
    Evaluation& operator/=(const Evaluation& other)
    { return Ops::divEq(*this, other); }

    // multiply value and derivatives by value of other
    template <class RhsValueType>
    Evaluation& operator/=(const RhsValueType& other)
    { return Ops::scalarDivEq(*this, other); }

    // add two evaluation objects
    Evaluation operator+(const Evaluation& other) const
    {
        Evaluation result(*this);
        result += other;
        return result;
    }

    // add constant to this object
    template <class RhsValueType>
    Evaluation operator+(const RhsValueType& other) const
    {
        Evaluation result(*this);
        result += other;
        return result;
    }

    // subtract two evaluation objects
    Evaluation operator-(const Evaluation& other) const
    {
        Evaluation result(*this);
        return (result -= other);
    }

    // subtract constant from evaluation object
    template <class RhsValueType>
    Evaluation operator-(const RhsValueType& other) const
    {
        Evaluation result(*this);
        return (result -= other);
    }

    // negation (unary minus) operator
    Evaluation operator-() const
    {
        Evaluation result;
        // set value and derivatives to negative
        Ops::assignNegative(result, *this);

        return result;
    }

    Evaluation operator*(const Evaluation& other) const
    {
        Evaluation result(*this);
        result *= other;
        return result;
    }

    template <class RhsValueType>
    Evaluation operator*(const RhsValueType& other) const
    {
        Evaluation result(*this);
        result *= other;
        return result;
    }

    Evaluation operator/(const Evaluation& other) const
    {
        Evaluation result(*this);
        result /= other;
        return result;
    }

    template <class RhsValueType>
    Evaluation operator/(const RhsValueType& other) const
    {
        Evaluation result(*this);
        result /= other;
        return result;
    }

    template <class RhsValueType>
    Evaluation& operator=(const RhsValueType& other)
    {
        setValue( other );
        clearDerivatives();
        return *this;
    }

    // copy assignment from evaluation
    Evaluation& operator=(const Evaluation& other)
    {
        Ops::assign(*this, other);
        return *this;
    }

    template <class RhsValueType>
    bool operator==(const RhsValueType& other) const
    { return value() == other; }

    bool operator==(const Evaluation& other) const
    {
        for (int idx = 0; idx < length_; ++idx)
            if (data_[idx] != other.data_[idx])
                return false;

        return true;
    }

    bool operator!=(const Evaluation& other) const
    { return !operator==(other); }

    template <class RhsValueType>
    bool operator>(RhsValueType other) const
    { return value() > other; }

    bool operator>(const Evaluation& other) const
    { return value() > other.value(); }

    template <class RhsValueType>
    bool operator<(RhsValueType other) const
    { return value() < other; }

    bool operator<(const Evaluation& other) const
    { return value() < other.value(); }

    template <class RhsValueType>
    bool operator>=(RhsValueType other) const
    { return value() >= other; }

    bool operator>=(const Evaluation& other) const
    { return value() >= other.value(); }

    template <class RhsValueType>
    bool operator<=(RhsValueType other) const
    { return value() <= other; }

    bool operator<=(const Evaluation& other) const
    { return value() <= other.value(); }

    // return value of variable
    const ValueType& value() const
    { return data_[valuepos_]; }

    // set value of variable
    template <class RhsValueType>
    void setValue(const RhsValueType& val)
    { data_[valuepos_] = val; }

    // return varIdx'th derivative
    const ValueType& derivative(int varIdx) const
    {
        assert(varIdx < numVars);
        return data_[varIdx + dstart_];
    }

    // set derivative at position varIdx
    void setDerivative(int varIdx, const ValueType& derVal)
    {
        assert(varIdx < numVars);
        data_[varIdx + dstart_] = derVal;
    }

protected:
    typename Ops::DataVector data_;
};

template <class RhsValueType, class ValueType, int numVars>
bool operator<(const RhsValueType& a, const Evaluation<ValueType, numVars>& b)
{ return b > a; }

template <class RhsValueType, class ValueType, int numVars>
bool operator>(const RhsValueType& a, const Evaluation<ValueType, numVars>& b)
{ return b < a; }

template <class RhsValueType, class ValueType, int numVars>
bool operator<=(const RhsValueType& a, const Evaluation<ValueType, numVars>& b)
{ return b >= a; }

template <class RhsValueType, class ValueType, int numVars>
bool operator>=(const RhsValueType& a, const Evaluation<ValueType, numVars>& b)
{ return b <= a; }

template <class RhsValueType, class ValueType, int numVars>
bool operator!=(const RhsValueType& a, const Evaluation<ValueType, numVars>& b)
{ return a != b.value(); }

template <class RhsValueType, class ValueType, int numVars>
Evaluation<ValueType, numVars> operator+(const RhsValueType& a, const Evaluation<ValueType, numVars>& b)
{
    Evaluation<ValueType, numVars> result(b);
    result += a;
    return result;
}

template <class RhsValueType, class ValueType, int numVars>
Evaluation<ValueType, numVars> operator-(const RhsValueType& a, const Evaluation<ValueType, numVars>& b)
{
    Evaluation<ValueType, numVars> result(a);
    result -= b;
    return result;
}

template <class RhsValueType, class ValueType, int numVars>
Evaluation<ValueType, numVars> operator/(const RhsValueType& a, const Evaluation<ValueType, numVars>& b)
{
    typedef EvaluationOps<ValueType, numVars> Ops;

    return Ops::divide(a, b);
}

template <class RhsValueType, class ValueType, int numVars>
Evaluation<ValueType, numVars> operator*(const RhsValueType& a, const Evaluation<ValueType, numVars>& b)
{
    Evaluation<ValueType, numVars> result(b);
    result *= a;
    return result;
}

template <class ValueType, int numVars>
std::ostream& operator<<(std::ostream& os, const Evaluation<ValueType, numVars>& eval)
{
    os << eval.value();
    return os;
}

} // namespace DenseAd
} // namespace Opm

// In Dune 2.3, the Evaluation.hpp header must be included before the fmatrix.hh
// header. Dune 2.4+ does not suffer from this because of some c++-foo.
//
// for those who are wondering: in C++ function templates cannot be partially
// specialized, and function argument overloads must be known _before_ they are used. The
// latter is what we do for the 'Dune::fvmeta::absreal()' function.
//
// consider the following test program:
//
// double foo(double i)
// { return i; }
//
// void bar()
// { std::cout << foo(0) << "\n"; }
//
// int foo(int i)
// { return i + 1; }
//
// void foobar()
// { std::cout << foo(0) << "\n"; }
//
// this will print '0' for bar() and '1' for foobar()...
#if !(DUNE_VERSION_NEWER(DUNE_COMMON, 2,4))

namespace Opm {
namespace DenseAd {
template <class ValueType, int numVars>
Evaluation<ValueType, numVars> abs(const Evaluation<ValueType, numVars>&);
}}

namespace std {
template <class ValueType, int numVars>
const Opm::DenseAd::Evaluation<ValueType, numVars> abs(const Opm::DenseAd::Evaluation<ValueType, numVars>& x)
{ return Opm::DenseAd::abs(x); }

} // namespace std

#if defined DUNE_DENSEMATRIX_HH
#warning \
 "Due to some C++ peculiarity regarding function overloads, the 'Evaluation.hpp'" \
 "header file must be included before Dune's 'densematrix.hh' for Dune < 2.4. " \
 "(If Evaluations are to be used in conjunction with a dense matrix.)"
#endif

#endif

// this makes the Dune matrix/vector classes happy...
#include <dune/common/ftraits.hh>

namespace Dune {
template <class ValueType, int numVars>
struct FieldTraits<Opm::DenseAd::Evaluation<ValueType, numVars> >
{
public:
    typedef Opm::DenseAd::Evaluation<ValueType, numVars> field_type;
    // setting real_type to field_type here potentially leads to slightly worse
    // performance, but at least it makes things compile.
    typedef field_type real_type;
};

} // namespace Dune

#include <opm/material/densead/EvaluationSpecializations.hpp>

#endif
