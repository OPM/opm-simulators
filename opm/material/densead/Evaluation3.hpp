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

#ifndef OPM_LOCAL_AD_EVALUATION_3_HPP
#define OPM_LOCAL_AD_EVALUATION_3_HPP

#include "Math.hpp"

#include <opm/common/Valgrind.hpp>

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
template <class ValueT>
class Evaluation< ValueT, 3 >
{
    static constexpr int numVars = 3;
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
    Evaluation(const Evaluation& other) //: data_( other.data_ )
    {
        data_[ 0 ] = other.data_[ 0 ];
        data_[ 1 ] = other.data_[ 1 ];
        data_[ 2 ] = other.data_[ 2 ];
        data_[ 3 ] = other.data_[ 3 ];
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
        data_[ 1 ] = 0;
        data_[ 2 ] = 0;
        data_[ 3 ] = 0;
    }


    // create a function evaluation for a "naked" depending variable (i.e., f(x) = x)
    template <class RhsValueType>
    static Evaluation devide(const RhsValueType& a, const Evaluation& b )
    {
        Evaluation<ValueType, numVars> result;
        result.setValue( a/b.value() );
        const ValueType df_dg = - result.value()/b.value();
        result.data_[ 1 ] = df_dg*b.data_[ 1 ];
        result.data_[ 2 ] = df_dg*b.data_[ 2 ];
        result.data_[ 3 ] = df_dg*b.data_[ 3 ];
        return result;
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
       data_[ 1 ] = other.data_[ 1 ];
       data_[ 2 ] = other.data_[ 2 ];
       data_[ 3 ] = other.data_[ 3 ];
    }


    // add value and derivatives from other to this values and derivatives
    Evaluation& operator+=(const Evaluation& other)
    {
        data_[ 0 ] += other.data_[ 0 ];
        data_[ 1 ] += other.data_[ 1 ];
        data_[ 2 ] += other.data_[ 2 ];
        data_[ 3 ] += other.data_[ 3 ];

        return *this;
    }

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
    {
        // value and derivatives are subtracted
        data_[ 0 ] -= other.data_[ 0 ];
        data_[ 1 ] -= other.data_[ 1 ];
        data_[ 2 ] -= other.data_[ 2 ];
        data_[ 3 ] -= other.data_[ 3 ];

        return *this;
    }

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
    {
        // while the values are multiplied, the derivatives follow the product rule,
        // i.e., (u*v)' = (v'u + u'v).
        const ValueType u = value();
        const ValueType v = other.value();

        data_[ 0 ] *= v ;
        data_[ 1 ]  = data_[ 1 ] * v + other.data_[ 1 ] * u;
        data_[ 2 ]  = data_[ 2 ] * v + other.data_[ 2 ] * u;
        data_[ 3 ]  = data_[ 3 ] * v + other.data_[ 3 ] * u;

        return *this;
    }

    // m(u*v)' = (v'u + u'v)
    template <class RhsValueType>
    Evaluation& operator*=(RhsValueType other)
    {
        data_[ 0 ] *= other;
        data_[ 1 ] *= other;
        data_[ 2 ] *= other;
        data_[ 3 ] *= other;

        return *this;
    }

    // m(u*v)' = (v'u + u'v)
    Evaluation& operator/=(const Evaluation& other)
    {
        // values are divided, derivatives follow the rule for division, i.e., (u/v)' = (v'u -
        // u'v)/v^2.
        const ValueType v_vv = 1.0 / other.value();
        const ValueType u_vv = value() * v_vv * v_vv;

        data_[ 0 ]  *= v_vv;
        data_[ 1 ]   = data_[ 1 ] * v_vv - other.data_[ 1 ] * u_vv ;
        data_[ 2 ]   = data_[ 2 ] * v_vv - other.data_[ 2 ] * u_vv ;
        data_[ 3 ]   = data_[ 3 ] * v_vv - other.data_[ 3 ] * u_vv ;

        return *this;
    }

    // multiply value and derivatives by value of other
    template <class RhsValueType>
    Evaluation& operator/=(const RhsValueType& other)
    {
        // values and derivatives are divided
        const ValueType factor = (1.0/other);
        return this->operator*=( factor );
    }

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
        result -= other;
        return result;
    }

    // subtract constant from evaluation object
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
        // set value and derivatives to negative
        for (int idx = 0; idx < length_; ++idx)
            result.data_[idx] = - data_[idx];

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
        data_ = other.data_;
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
    bool operator>(const RhsValueType& other) const
    { return value() > other; }

    bool operator>(const Evaluation& other) const
    { return value() > other.value(); }

    template <class RhsValueType>
    bool operator<(const RhsValueType& other) const
    { return value() < other; }

    bool operator<(const Evaluation& other) const
    { return value() < other.value(); }

    template <class RhsValueType>
    bool operator>=(const RhsValueType& other) const
    { return value() >= other; }

    bool operator>=(const Evaluation& other) const
    { return value() >= other.value(); }

    template <class RhsValueType>
    bool operator<=(const RhsValueType& other) const
    { return value() <= other; }

    bool operator<=(const Evaluation& other) const
    { return value() <= other.value(); }

    // return value of variable
    const ValueType& value() const
    { return data_[valuepos_]; }

    // set value of variable
    void setValue(const ValueType& val)
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
    std::array<ValueType, length_> data_;
};

} // namespace DenseAD
} // namespace Dune

// #include <opm/material/densead/EvaluationSIMD.hpp>

#endif
