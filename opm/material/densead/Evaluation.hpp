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

#include <opm/material/common/Valgrind.hpp>

#include <dune/common/version.hh>

#include <array>
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>

namespace Opm {
namespace DenseAd {
/*!
 * \brief Represents a function evaluation and its derivatives w.r.t. a fixed set of
 *        variables.
 */
template <class ValueT, int numVars>
class Evaluation
{
public:
    typedef ValueT ValueType;

    enum { size = numVars };

    Evaluation()
    {}

    // copy other function evaluation
    Evaluation(const Evaluation& other)
    {
        // copy evaluated function value and derivatives
        value = other.value;
        std::copy(other.derivatives.begin(), other.derivatives.end(), derivatives.begin());
    }

    // create an evaluation which represents a constant function
    //
    // i.e., f(x) = c. this implies an evaluation with the given value and all
    // derivatives being zero.
    template <class RhsValueType>
    Evaluation(const RhsValueType& c)
    {
        value = c;
        std::fill(derivatives.begin(), derivatives.end(), 0.0);
    }

    // create a function evaluation for a "naked" depending variable (i.e., f(x) = x)
    template <class RhsValueType>
    static Evaluation createVariable(const RhsValueType& value, unsigned varPos)
    {
        // The variable position must be in represented by the given variable descriptor
        assert(0 <= varPos && varPos < size);

        Evaluation result;

        // copy function value and set all derivatives to 0, except for the variable
        // which is represented by the value (which is set to 1.0)
        result.value = value;
        std::fill(result.derivatives.begin(), result.derivatives.end(), 0.0);
        result.derivatives[varPos] = 1.0;

        return result;
    }

    // "evaluate" a constant function (i.e. a function that does not depend on the set of
    // relevant variables, f(x) = c).
    template <class RhsValueType>
    static Evaluation createConstant(const RhsValueType& value)
    {
        Evaluation result;
        result.value = value;
        std::fill(result.derivatives.begin(), result.derivatives.end(), 0.0);
        Valgrind::CheckDefined(result.value);
        Valgrind::CheckDefined(result.derivatives);
        return result;
    }

    // print the value and the derivatives of the function evaluation
    void print(std::ostream& os = std::cout) const
    {
        os << "v: " << value << " / d:";
        for (int varIdx = 0; varIdx < derivatives.size(); ++varIdx)
            os << " " << derivatives[varIdx];
    }


    Evaluation& operator+=(const Evaluation& other)
    {
        // value and derivatives are added
        this->value += other.value;
        for (unsigned varIdx = 0; varIdx < size; ++varIdx)
            this->derivatives[varIdx] += other.derivatives[varIdx];

        return *this;
    }

    template <class RhsValueType>
    Evaluation& operator+=(const RhsValueType& other)
    {
        // value is added, derivatives stay the same
        this->value += other;

        return *this;
    }

    Evaluation& operator-=(const Evaluation& other)
    {
        // value and derivatives are subtracted
        this->value -= other.value;
        for (unsigned varIdx = 0; varIdx < size; ++varIdx)
            this->derivatives[varIdx] -= other.derivatives[varIdx];

        return *this;
    }

    template <class RhsValueType>
    Evaluation& operator-=(const RhsValueType& other)
    {
        // for constants, values are subtracted, derivatives stay the same
        this->value -= other;

        return *this;
    }

    Evaluation& operator*=(const Evaluation& other)
    {
        // while the values are multiplied, the derivatives follow the product rule,
        // i.e., (u*v)' = (v'u + u'v).
        const ValueType& u = this->value;
        const ValueType& v = other.value;
        for (unsigned varIdx = 0; varIdx < size; ++varIdx) {
            const ValueType& uPrime = this->derivatives[varIdx];
            const ValueType& vPrime = other.derivatives[varIdx];

            this->derivatives[varIdx] = (v*uPrime + u*vPrime);
        }
        this->value *= v;

        return *this;
    }

    template <class RhsValueType>
    Evaluation& operator*=(RhsValueType other)
    {
        // values and derivatives are multiplied
        this->value *= other;
        for (unsigned varIdx = 0; varIdx < size; ++varIdx)
            this->derivatives[varIdx] *= other;

        return *this;
    }

    Evaluation& operator/=(const Evaluation& other)
    {
        // values are divided, derivatives follow the rule for division, i.e., (u/v)' = (v'u -
        // u'v)/v^2.
        const ValueType& u = this->value;
        const ValueType& v = other.value;
        for (unsigned varIdx = 0; varIdx < size; ++varIdx) {
            const ValueType& uPrime = this->derivatives[varIdx];
            const ValueType& vPrime = other.derivatives[varIdx];

            this->derivatives[varIdx] = (v*uPrime - u*vPrime)/(v*v);
        }
        this->value /= v;

        return *this;
    }

    template <class RhsValueType>
    Evaluation& operator/=(const RhsValueType& other)
    {
        // values and derivatives are divided
        auto tmp = 1.0/other;
        this->value *= tmp;
        for (unsigned varIdx = 0; varIdx < size; ++varIdx)
            this->derivatives[varIdx] *= tmp;

        return *this;
    }

    Evaluation operator+(const Evaluation& other) const
    {
        Evaluation result(*this);
        result += other;
        return result;
    }

    template <class RhsValueType>
    Evaluation operator+(const RhsValueType& other) const
    {
        Evaluation result(*this);
        result += other;
        return result;
    }

    Evaluation operator-(const Evaluation& other) const
    {
        Evaluation result(*this);
        result -= other;
        return result;
    }

    template <class RhsValueType>
    Evaluation operator-(const RhsValueType& other) const
    {
        Evaluation result(*this);
        result -= other;
        return result;
    }

    // negation (unary minus) operator
    Evaluation operator-() const
    {
        Evaluation result;
        result.value = -this->value;
        for (unsigned varIdx = 0; varIdx < size; ++varIdx)
            result.derivatives[varIdx] = - this->derivatives[varIdx];

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
        this->value = other;
        std::fill(this->derivatives.begin(), this->derivatives.end(), 0.0);
        return *this;
    }

    // copy assignment from evaluation
    Evaluation& operator=(const Evaluation& other)
    {
        this->value = other.value;
        std::copy(other.derivatives.begin(), other.derivatives.end(), this->derivatives.begin());
        return *this;
    }

    template <class RhsValueType>
    bool operator==(const RhsValueType& other) const
    { return this->value == other; }

    bool operator==(const Evaluation& other) const
    {
        if (this->value != other.value)
            return false;

        for (unsigned varIdx = 0; varIdx < size; ++varIdx)
            if (this->derivatives[varIdx] != other.derivatives[varIdx])
                return false;

        return true;
    }

    bool operator!=(const Evaluation& other) const
    { return !operator==(other); }

    template <class RhsValueType>
    bool operator>(RhsValueType other) const
    { return this->value > other; }

    bool operator>(const Evaluation& other) const
    { return this->value > other.value; }

    template <class RhsValueType>
    bool operator<(RhsValueType other) const
    { return this->value < other; }

    bool operator<(const Evaluation& other) const
    { return this->value < other.value; }

    template <class RhsValueType>
    bool operator>=(RhsValueType other) const
    { return this->value >= other; }

    bool operator>=(const Evaluation& other) const
    { return this->value >= other.value; }

    template <class RhsValueType>
    bool operator<=(RhsValueType other) const
    { return this->value <= other; }

    bool operator<=(const Evaluation& other) const
    { return this->value <= other.value; }

    // maybe this should be made 'private'...
    ValueType value;
    std::array<ValueType, size> derivatives;
};

template <class RhsValueType, class ValueType, int numVars>
bool operator<(const RhsValueType& a, const Evaluation<ValueType, numVars> &b)
{ return b > a; }

template <class RhsValueType, class ValueType, int numVars>
bool operator>(const RhsValueType& a, const Evaluation<ValueType, numVars> &b)
{ return b < a; }

template <class RhsValueType, class ValueType, int numVars>
bool operator<=(const RhsValueType& a, const Evaluation<ValueType, numVars> &b)
{ return b >= a; }

template <class RhsValueType, class ValueType, int numVars>
bool operator>=(const RhsValueType& a, const Evaluation<ValueType, numVars> &b)
{ return b <= a; }

template <class RhsValueType, class ValueType, int numVars>
bool operator!=(const RhsValueType& a, const Evaluation<ValueType, numVars> &b)
{ return a != b.value; }

template <class RhsValueType, class ValueType, int numVars>
Evaluation<ValueType, numVars> operator+(const RhsValueType& a, const Evaluation<ValueType, numVars> &b)
{
    Evaluation<ValueType, numVars> result(b);

    result += a;

    return result;
}

template <class RhsValueType, class ValueType, int numVars>
Evaluation<ValueType, numVars> operator-(const RhsValueType& a, const Evaluation<ValueType, numVars> &b)
{
    Evaluation<ValueType, numVars> result;

    result.value = a - b.value;
    for (unsigned varIdx = 0; varIdx < numVars; ++varIdx)
        result.derivatives[varIdx] = - b.derivatives[varIdx];

    return result;
}

template <class RhsValueType, class ValueType, int numVars>
Evaluation<ValueType, numVars> operator/(const RhsValueType& a, const Evaluation<ValueType, numVars> &b)
{
    Evaluation<ValueType, numVars> result;

    result.value = a/b.value;

    // outer derivative
    const ValueType& df_dg = - a/(b.value*b.value);
    for (unsigned varIdx = 0; varIdx < numVars; ++varIdx)
        result.derivatives[varIdx] = df_dg*b.derivatives[varIdx];

    return result;
}

template <class RhsValueType, class ValueType, int numVars>
Evaluation<ValueType, numVars> operator*(const RhsValueType& a, const Evaluation<ValueType, numVars> &b)
{
    Evaluation<ValueType, numVars> result;

    result.value = a*b.value;
    for (unsigned varIdx = 0; varIdx < numVars; ++varIdx)
        result.derivatives[varIdx] = a*b.derivatives[varIdx];

    return result;
}

template <class ValueType, int numVars>
std::ostream& operator<<(std::ostream& os, const Evaluation<ValueType, numVars>& eval)
{
    os << eval.value;
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

#endif
