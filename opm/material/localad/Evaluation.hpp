// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2015 by Andreas Lauser

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
*/
/*!
 * \file
 *
 * \brief Representation of an evaluation of a function and its derivatives w.r.t. a set
 *        of variables in the localized OPM automatic differentiation (AD) framework.
 */
#ifndef OPM_LOCAL_AD_EVALUATION_HPP
#define OPM_LOCAL_AD_EVALUATION_HPP

#include <iostream>
#include <array>
#include <cassert>
#include <opm/material/common/Valgrind.hpp>

namespace Opm {
namespace LocalAd {
/*!
 * \brief Represents a function evaluation and its derivatives w.r.t. a fixed set of
 *        variables.
 */
template <class ScalarT, class VarSetTag, int numVars>
class Evaluation
{
public:
    typedef ScalarT Scalar;

    enum { size = numVars };

    Evaluation()
    {};

    // copy other function evaluation
    Evaluation(const Evaluation& other)
    {
        // copy evaluated function value and derivatives
        value = other.value;
        std::copy(other.derivatives.begin(), other.derivatives.end(), derivatives.begin());
    };

    // create an evaluation which represents a constant function
    //
    // i.e., f(x) = c. this implies an evaluation with the given value and all
    // derivatives being zero.
    Evaluation(Scalar c)
    {
        value = c;
        std::fill(derivatives.begin(), derivatives.end(), 0.0);
    };

    // create a function evaluation for a "naked" depending variable (i.e., f(x) = x)
    static Evaluation createVariable(Scalar value, int varPos)
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
    static Evaluation createConstant(Scalar value)
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
        for (int varIdx= 0; varIdx < size; ++varIdx)
            this->derivatives[varIdx] += other.derivatives[varIdx];

        return *this;
    }

    Evaluation& operator+=(Scalar other)
    {
        // value is added, derivatives stay the same
        this->value += other;

        return *this;
    }

    Evaluation& operator-=(const Evaluation& other)
    {
        // value and derivatives are subtracted
        this->value -= other.value;
        for (int varIdx= 0; varIdx < size; ++varIdx)
            this->derivatives[varIdx] -= other.derivatives[varIdx];

        return *this;
    }

    Evaluation& operator-=(Scalar other)
    {
        // for constants, values are subtracted, derivatives stay the same
        this->value -= other;

        return *this;
    }

    Evaluation& operator*=(const Evaluation& other)
    {
        // while the values are multiplied, the derivatives follow the product rule,
        // i.e., (u*v)' = (v'u + u'v).
        Scalar u = this->value;
        Scalar v = other.value;
        this->value *= v;
        for (int varIdx= 0; varIdx < size; ++varIdx) {
            Scalar uPrime = this->derivatives[varIdx];
            Scalar vPrime = other.derivatives[varIdx];

            this->derivatives[varIdx] = (v*uPrime + u*vPrime);
        }

        return *this;
    }

    Evaluation& operator*=(Scalar other)
    {
        // values and derivatives are multiplied
        this->value *= other;
        for (int varIdx= 0; varIdx < size; ++varIdx)
            this->derivatives[varIdx] *= other;

        return *this;
    }

    Evaluation& operator/=(const Evaluation& other)
    {
        // values are divided, derivatives follow the rule for division, i.e., (u/v)' = (v'u -
        // u'v)/v^2.
        Scalar u = this->value;
        Scalar v = other.value;
        this->value /= v;
        for (int varIdx= 0; varIdx < size; ++varIdx) {
            Scalar uPrime = this->derivatives[varIdx];
            Scalar vPrime = other.derivatives[varIdx];

            this->derivatives[varIdx] = (v*uPrime - u*vPrime)/(v*v);
        }

        return *this;
    }

    Evaluation& operator/=(Scalar other)
    {
        // values and derivatives are divided
        other = 1.0/other;
        this->value *= other;
        for (int varIdx= 0; varIdx < size; ++varIdx)
            this->derivatives[varIdx] *= other;

        return *this;
    }

    Evaluation operator+(const Evaluation& other) const
    {
        Evaluation result(*this);
        result += other;
        return result;
    }

    Evaluation operator+(Scalar other) const
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

    Evaluation operator-(Scalar other) const
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
        for (int varIdx= 0; varIdx < size; ++varIdx)
            result.derivatives[varIdx] = - this->derivatives[varIdx];

        return result;
    }

    Evaluation operator*(const Evaluation& other) const
    {
        Evaluation result(*this);
        result *= other;
        return result;
    }

    Evaluation operator*(Scalar other) const
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

    Evaluation operator/(Scalar other) const
    {
        Evaluation result(*this);
        result /= other;
        return result;
    }

    Evaluation& operator=(Scalar other)
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

    bool operator==(Scalar other) const
    { return this->value == other; }

    bool operator==(const Evaluation& other) const
    {
        if (this->value != other.value)
            return false;

        for (int varIdx= 0; varIdx < size; ++varIdx)
            if (this->derivatives[varIdx] != other.derivatives[varIdx])
                return false;

        return true;
    }

    bool isSame(const Evaluation& other, Scalar tolerance) const
    {
        Scalar tmp = other.value - other.value;
        if (std::abs(tmp) > tolerance && std::abs(tmp)/tolerance > 1.0)
            return false;

        for (int varIdx= 0; varIdx < size; ++varIdx) {
            Scalar tmp = other.derivatives[varIdx] - this->derivatives[varIdx];
            if (std::abs(tmp) > tolerance && std::abs(tmp)/tolerance > 1.0)
                return false;
        }

        return true;
    }

    bool operator!=(const Evaluation& other) const
    { return !operator==(other); }

    bool operator>(Scalar other) const
    { return this->value > other; }

    bool operator>(const Evaluation& other) const
    { return this->value > other.value; }

    bool operator<(Scalar other) const
    { return this->value < other; }

    bool operator<(const Evaluation& other) const
    { return this->value < other.value; }

    bool operator>=(Scalar other) const
    { return this->value >= other; }

    bool operator>=(const Evaluation& other) const
    { return this->value >= other.value; }

    bool operator<=(Scalar other) const
    { return this->value <= other; }

    bool operator<=(const Evaluation& other) const
    { return this->value <= other.value; }

    // maybe this should be made 'private'...
    Scalar value;
    std::array<Scalar, size> derivatives;
};

template <class ScalarA, class Scalar, class VarSetTag, int numVars>
bool operator<(const ScalarA& a, const Evaluation<Scalar, VarSetTag, numVars> &b)
{ return b > a; }

template <class ScalarA, class Scalar, class VarSetTag, int numVars>
bool operator>(const ScalarA& a, const Evaluation<Scalar, VarSetTag, numVars> &b)
{ return b < a; }

template <class ScalarA, class Scalar, class VarSetTag, int numVars>
bool operator<=(const ScalarA& a, const Evaluation<Scalar, VarSetTag, numVars> &b)
{ return b >= a; }

template <class ScalarA, class Scalar, class VarSetTag, int numVars>
bool operator>=(const ScalarA& a, const Evaluation<Scalar, VarSetTag, numVars> &b)
{ return b <= a; }

template <class ScalarA, class Scalar, class VarSetTag, int numVars>
bool operator!=(const ScalarA& a, const Evaluation<Scalar, VarSetTag, numVars> &b)
{ return a != b.value; }

template <class ScalarA, class Scalar, class VarSetTag, int numVars>
Evaluation<Scalar, VarSetTag, numVars> operator+(const ScalarA& a, const Evaluation<Scalar, VarSetTag, numVars> &b)
{
    Evaluation<Scalar, VarSetTag, numVars> result(b);

    result += a;

    return result;
}

template <class ScalarA, class Scalar, class VarSetTag, int numVars>
Evaluation<Scalar, VarSetTag, numVars> operator-(const ScalarA& a, const Evaluation<Scalar, VarSetTag, numVars> &b)
{
    Evaluation<Scalar, VarSetTag, numVars> result;

    result.value = a - b.value;
    for (int varIdx= 0; varIdx < numVars; ++varIdx)
        result.derivatives[varIdx] = - b.derivatives[varIdx];

    return result;
}

template <class ScalarA, class Scalar, class VarSetTag, int numVars>
Evaluation<Scalar, VarSetTag, numVars> operator/(const ScalarA& a, const Evaluation<Scalar, VarSetTag, numVars> &b)
{
    Evaluation<Scalar, VarSetTag, numVars> result;

    result.value = a/b.value;

    // outer derivative
    Scalar df_dg = - a/(b.value*b.value);
    for (int varIdx= 0; varIdx < numVars; ++varIdx)
        result.derivatives[varIdx] = df_dg*b.derivatives[varIdx];

    return result;
}

template <class ScalarA, class Scalar, class VarSetTag, int numVars>
Evaluation<Scalar, VarSetTag, numVars> operator*(const ScalarA& a, const Evaluation<Scalar, VarSetTag, numVars> &b)
{
    Evaluation<Scalar, VarSetTag, numVars> result;

    result.value = a*b.value;
    for (int varIdx= 0; varIdx < numVars; ++varIdx)
        result.derivatives[varIdx] = a*b.derivatives[varIdx];

    return result;
}

template <class Scalar, class VarSetTag, int numVars>
std::ostream& operator<<(std::ostream& os, const Evaluation<Scalar, VarSetTag, numVars>& eval)
{
    os << eval.value;
    return os;
}

} // namespace LocalAd
} // namespace Opm

// make the Dune matrix/vector classes happy. Obviously, this is not very elegant...

#ifdef DUNE_DENSEMATRIX_HH
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
//
#error "Due to some C++ peculiarity regarding function template specialization,"
#error "the 'evaluation.hh' header must be included before Dune's 'densematrix.hh'!"
#endif

#include <dune/common/ftraits.hh>

namespace Dune {
template <class Scalar, class VarSetTag, int numVars>
struct FieldTraits<Opm::LocalAd::Evaluation<Scalar, VarSetTag, numVars> >
{
public:
    typedef Opm::LocalAd::Evaluation<Scalar, VarSetTag, numVars> field_type;
    typedef Scalar real_type;
};

namespace fvmeta {
template <class Scalar, class VarSetTag, int numVars>
inline Scalar absreal(const Opm::LocalAd::Evaluation<Scalar, VarSetTag, numVars>& k)
{
    return std::abs(k.value);
}

}} // namespace fvmeta, Dune

#endif
