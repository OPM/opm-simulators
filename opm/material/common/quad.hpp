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
 * \brief This file provides the infrastructure to use quad-precision
 *        floating point values in the numerical models.
 */
#if !defined OPM_COMMON_QUAD_HPP && HAVE_QUAD
#define OPM_COMMON_QUAD_HPP

#include <opm/material/common/Unused.hpp>

#include <cmath>
#include <complex>
#include <string>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <type_traits>

extern "C" {
#include <quadmath.h>
}

typedef __float128 quad;

namespace std {

// provide the numeric limits for the quad precision type
template <>
class numeric_limits<quad>
{
public:
    static constexpr bool is_specialized = true;

    static constexpr quad min() throw()
    { return FLT128_MIN; }
    static constexpr quad max() throw()
    { return FLT128_MAX; }

    // number of bits in mantissa
    static constexpr int digits = FLT128_MANT_DIG;
    // number of decimal digits
    static constexpr int digits10 = FLT128_DIG;
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = false;
    static constexpr int radix = 0;
    static constexpr quad epsilon() throw()
    { return FLT128_EPSILON; }
    static constexpr quad round_error() throw()
    { return 0.5; }

    static constexpr int min_exponent = FLT128_MIN_EXP;
    static constexpr int min_exponent10 = FLT128_MIN_10_EXP;
    static constexpr int max_exponent = FLT128_MAX_EXP;
    static constexpr int max_exponent10 = FLT128_MAX_10_EXP;

    static constexpr bool has_infinity = true;
    static constexpr bool has_quiet_NaN = true;
    static constexpr bool has_signaling_NaN = true;
    static constexpr float_denorm_style has_denorm = denorm_present;
    static constexpr bool has_denorm_loss = false;
    static constexpr quad infinity() throw()
    { return __builtin_huge_valq(); }
    static constexpr quad quiet_NaN() throw()
    { return __builtin_nan(""); }
    static constexpr quad signaling_NaN() throw()
    { return __builtin_nans(""); }
    static constexpr quad denorm_min() throw()
    { return FLT128_DENORM_MIN; }

    static constexpr bool is_iec559 = true;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo = false;

    static constexpr bool traps = std::numeric_limits<double>::traps;
    static constexpr bool tinyness_before = std::numeric_limits<double>::tinyness_before;
    static constexpr float_round_style round_style = round_to_nearest;
};

// provide some type traits for the quadruple precision type
template <>
struct is_floating_point<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_arithmetic<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_fundamental<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_scalar<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_pod<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_signed<quad>
    : public integral_constant<bool, true>
{};


template <>
struct is_standard_layout<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_trivial<quad>
    : public integral_constant<bool, true>
{};

/*
template <>
struct is_trivially_copyable<quad>
    : public integral_constant<bool, true>
{};
*/

template <class OtherType>
struct is_assignable<quad, OtherType>
    : public integral_constant<bool, is_arithmetic<OtherType>::value>
{};

template <class OtherType>
struct is_nothrow_assignable<quad, OtherType>
    : public is_assignable<quad, OtherType>
{};

/*
template <class OtherType>
struct is_trivially_assignable<quad, OtherType>
    : public integral_constant<bool, is_arithmetic<OtherType>::value>
{};
*/

template <>
struct is_copy_assignable<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_nothrow_copy_assignable<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_move_assignable<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_nothrow_move_assignable<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_constructible<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_nothrow_constructible<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_default_constructible<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_nothrow_default_constructible<quad>
    : public integral_constant<bool, true>
{};

/*
template <>
struct is_trivially_default_constructible<quad>
    : public integral_constant<bool, true>
{};
*/

template <>
struct is_copy_constructible<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_move_constructible<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_nothrow_move_constructible<quad>
    : public integral_constant<bool, true>
{};


template <>
struct is_destructible<quad>
    : public integral_constant<bool, true>
{};

template <>
struct is_nothrow_destructible<quad>
    : public integral_constant<bool, true>
{};

template <class OtherType>
struct is_convertible<quad, OtherType>
    : public is_arithmetic<OtherType>
{ };

inline std::ostream& operator<<(std::ostream& os, const quad& val)
{
    if (os.precision() > std::numeric_limits<double>::digits10)
        throw std::runtime_error("The precision requested for output cannot "
                                 "be represented by a double precision floating "
                                 "point object");

    return os << static_cast<double>(val);
}

inline std::istream& operator>>(std::istream& is, quad& val)
{
    double tmp;
    std::istream& ret = (is >> tmp);
    val = tmp;
    return ret;
}

inline quad real(quad val)
{ return val; }

inline quad real(const std::complex<quad>& val)
{ return val.real(); }

inline quad imag(quad val OPM_UNUSED)
{ return 0.0; }

inline quad imag(const std::complex<quad>& val)
{ return val.imag(); }

inline quad abs(quad val)
{ return (val < 0) ? -val : val; }

inline quad floor(quad val)
{ return floorq(val); }

inline quad ceil(quad val)
{ return ceilq(val); }

inline quad max(quad a, quad b)
{ return (a > b) ? a : b; }

inline quad min(quad a, quad b)
{ return (a < b) ? a : b; }

inline quad sqrt(quad val)
{ return sqrtq(val); }

template <class ExpType>
inline quad pow(quad base, ExpType exp)
{ return powq(base, static_cast<quad>(exp)); }

template <class BaseType>
inline quad pow(BaseType base, quad exp)
{ return powq(static_cast<quad>(base), exp); }

inline quad pow(quad base, quad exp)
{ return powq(base, exp); }

inline quad exp(quad val)
{ return expq(val); }

inline quad log(quad val)
{ return logq(val); }

inline quad sin(quad val)
{ return sinq(val); }

inline quad cos(quad val)
{ return cosq(val); }

inline quad tan(quad val)
{ return tanq(val); }

inline quad atan(quad val)
{ return atanq(val); }

inline quad atan2(quad a, quad b)
{ return atan2q(a, b); }

inline quad round(quad val)
{ return roundq(val); }

inline bool isfinite(quad val)
{ return finiteq(val); }

inline bool isnan(quad val)
{ return isnanq(val); }

inline bool isinf(quad val)
{ return isinfq(val); }

} // namespace std

// specialize Dune::className for __float128 since it former does not work properly with
// __float128 (this is mainly the fault of GCC/libstdc++)
#include <dune/common/classname.hh>

namespace Dune {
template <>
inline std::string className<__float128>()
{ return "quad"; }
} // namespace Dune

#endif // OPM_COMMON_QUAD_HPP
