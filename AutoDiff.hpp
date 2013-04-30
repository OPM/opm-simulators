/*===========================================================================
//
// File: AutoDiff.hpp
//
// Created: 2013-04-29 10:51:23+0200
//
// Authors: Knut-Andreas Lie      <Knut-Andreas.Lie@sintef.no>
//          Halvor M. Nilsen      <HalvorMoll.Nilsen@sintef.no>
//          Atgeirr F. Rasmussen  <atgeirr@sintef.no>
//          Xavier Raynaud        <Xavier.Raynaud@sintef.no>
//          Bård Skaflestad       <Bard.Skaflestad@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2013 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include <cmath>

#ifndef OPM_AUTODIFF_HPP_HEADER
#define OPM_AUTODIFF_HPP_HEADER

namespace AutoDiff {
    template <typename Scalar>
    class Forward {
    public:
        explicit Forward(const Scalar& x)
            : val_(x), der_(Scalar(1))
        {}

        Forward(const Scalar x, const Scalar dx)
            : val_(x), der_(dx)
        {}

        void
        operator +=(const Scalar& rhs)
        {
            val_ += rhs;
        }

        void
        operator +=(const Forward& rhs)
        {
            val_ += rhs.val_;
            der_ += rhs.der_;
        }

        void
        operator -=(const Scalar& rhs)
        {
            val_ -= rhs;
        }

        void
        operator -=(const Forward& rhs)
        {
            val_ -= rhs.val_;
            der_ -= rhs.der_;
        }

        void
        operator *=(const Scalar& rhs)
        {
            val_ *= rhs;
            der_ *= rhs;
        }

        void
        operator *=(const Forward& rhs)
        {
            der_   = der_*rhs.val_ + val_*rhs.der_;
            val_  *= rhs.val_;
        }

        void
        operator /=(const Scalar& rhs)
        {
            val_ /= rhs;
            der_ /= rhs;
        }

        void
        operator /=(const Forward& rhs)
        {
            der_   = (der_*rhs.val_ - val_*rhs.der_) / (rhs.val_ * rhs.val_);
            val_  /= rhs.val_;
        }

        template <class Ostream>
        Ostream&
        print(Ostream& os) const
        {
            os << "(x,dx) = (" << val_ << ',' << der_ << ")";

            return os;
        }

        const Scalar val() const { return val_; }
        const Scalar der() const { return der_; }

    private:
        Scalar val_ ;
        Scalar der_;

    };

    template <class Ostream, typename Scalar>
    Ostream&
    operator<<(Ostream& os, const Forward<Scalar>& fw)
    {
        return fw.print(os);
    }

    template <typename Scalar>
    Forward<Scalar>
    operator +(const Forward<Scalar>& lhs,
               const Forward<Scalar>& rhs)
    {
        Forward<Scalar> ret = lhs;
        ret += rhs;

        return ret;
    }

    template <typename Scalar, typename T>
    Forward<Scalar>
    operator +(const T                lhs,
               const Forward<Scalar>& rhs)
    {
        Forward<Scalar> ret = rhs;
        ret += Scalar(lhs);

        return ret;
    }

    template <typename Scalar, typename T>
    Forward<Scalar>
    operator +(const Forward<Scalar>& lhs,
               const T                rhs)
    {
        Forward<Scalar> ret = lhs;
        ret += Scalar(rhs);

        return ret;
    }

    template <typename Scalar>
    Forward<Scalar>
    operator -(const Forward<Scalar>& lhs,
               const Forward<Scalar>& rhs)
    {
        Forward<Scalar> ret = lhs;
        ret -= rhs;

        return ret;
    }

    template <typename Scalar, typename T>
    Forward<Scalar>
    operator -(const T                lhs,
               const Forward<Scalar>& rhs)
    {
        Forward<Scalar> ret(Scalar(lhs) - rhs.val(), -rhs.der());

        return ret;
    }

    template <typename Scalar, typename T>
    Forward<Scalar>
    operator -(const Forward<Scalar>& lhs,
               const T                rhs)
    {
        Forward<Scalar> ret = lhs;
        ret -= Scalar(rhs);

        return ret;
    }

    template <typename Scalar>
    Forward<Scalar>
    operator *(const Forward<Scalar>& lhs,
               const Forward<Scalar>& rhs)
    {
        Forward<Scalar> ret = lhs;
        ret *= rhs;

        return ret;
    }

    template <typename Scalar, typename T>
    Forward<Scalar>
    operator *(const T                lhs,
               const Forward<Scalar>& rhs)
    {
        Forward<Scalar> ret = rhs;
        ret *= Scalar(lhs);

        return ret;
    }

    template <typename Scalar, typename T>
    Forward<Scalar>
    operator *(const Forward<Scalar>& lhs,
               const T                rhs)
    {
        Forward<Scalar> ret = lhs;
        ret *= Scalar(rhs);

        return ret;
    }

    template <typename Scalar>
    Forward<Scalar>
    operator /(const Forward<Scalar>& lhs,
               const Forward<Scalar>& rhs)
    {
        Forward<Scalar> ret = lhs;
        ret /= rhs;

        return ret;
    }

    template <typename Scalar, typename T>
    Forward<Scalar>
    operator /(const T                lhs,
               const Forward<Scalar>& rhs)
    {
        Scalar a =  Scalar(lhs) / rhs.val();
        Scalar b = -Scalar(lhs) / (rhs.val() * rhs.val());

        Forward<Scalar> ret(a, b);

        return ret;
    }

    template <typename Scalar, typename T>
    Forward<Scalar>
    operator /(const Forward<Scalar>& lhs,
               const T                rhs)
    {
        Scalar a = lhs.val() / Scalar(rhs);
        Scalar b = lhs.der() / Scalar(rhs);

        Forward<Scalar> ret(a, b);

        return ret;
    }

    template <typename Scalar>
    Forward<Scalar>
    cos(const Forward<Scalar>& x)
    {
        Forward<Scalar> ret( std::cos(x.val()),
                            -std::sin(x.val()) * x.der());

        return ret;
    }

    template <typename Scalar>
    Forward<Scalar>
    sqrt(const Forward<Scalar>& x)
    {
        Scalar a = std::sqrt(x.val());
        Scalar b = Scalar(1.0) / (Scalar(2.0) * a);
        Forward<Scalar> ret(a, b * x.der());

        return ret;
    }
} // namespace AutoDiff

namespace std {
    using AutoDiff::cos;
    using AutoDiff::sqrt;
}

#endif  /* OPM_AUTODIFF_HPP_HEADER */
