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
            : x_(x), dx_(Scalar(1))
        {}

        Forward(const Scalar x, const Scalar dx)
            : x_(x), dx_(dx)
        {}

        Forward&
        operator +=(const Scalar& rhs)
        {
            x_ += rhs;

            return *this;
        }

        Forward&
        operator +=(const Forward& rhs)
        {
            x_  += rhs.x_;
            dx_ += rhs.dx_;

            return *this;
        }

        Forward&
        operator -=(const Scalar& rhs)
        {
            x_ -= rhs;

            return *this;
        }

        Forward&
        operator -=(const Forward& rhs)
        {
            x_  -= rhs.x_;
            dx_ -= rhs.dx_;

            return *this;
        }

        Forward&
        operator *=(const Scalar& rhs)
        {
            x_  *= rhs;
            dx_ *= rhs;

            return *this;
        }

        Forward&
        operator *=(const Forward& rhs)
        {
            x_  *= rhs.x_;
            dx_ *= dx_*rhs.x_ + x_*rhs.dx_;

            return *this;
        }

        Forward&
        operator /=(const Scalar& rhs)
        {
            x_  /= rhs;
            dx_ /= rhs;

            return *this;
        }

        Forward&
        operator /=(const Forward& rhs)
        {
            x_  /= rhs.x_;
            dx_  = (dx_*rhs.x_ - x_*rhs.dx_) / (rhs.x_ * rhs.x_);

            return *this;
        }

        template <class Ostream>
        Ostream&
        print(Ostream& os) const
        {
            os << "(x,dx) = (" << x_ << ',' << dx_ << ")";

            return os;
        }

        const Scalar val() const { return x_ ; }
        const Scalar der() const { return dx_; }

    private:
        Scalar x_ ;
        Scalar dx_;

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
        Scalar a = rhs.val() / Scalar(lhs);
        Scalar b = rhs.der() / Scalar(lhs);

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
