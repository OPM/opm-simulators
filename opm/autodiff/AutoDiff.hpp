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

#ifndef OPM_AUTODIFF_HPP_HEADER
#define OPM_AUTODIFF_HPP_HEADER

#include <cmath>

namespace Opm
{

    /// A simple class for forward-mode automatic differentiation.
    ///
    /// The class represents a single value and a single derivative.
    /// Only basic arithmetic operators and a few functions are
    /// implemented for it, it is mostly intended for simple
    /// experimentation.
    template <typename Scalar>
    class AutoDiff
    {
    public:
        /// Create an AutoDiff object representing a constant, that
        /// is, its derivative is zero.
        static AutoDiff
        constant(const Scalar x)
        {
            return function(x, Scalar(0));
        }

        /// Create an AutoDiff object representing a primary variable,
        /// that is, its derivative is one.
        static AutoDiff
        variable(const Scalar x)
        {
            return function(x, Scalar(1));
        }

        /// Create an AutoDiff object representing a function value
        /// and its derivative.
        static AutoDiff
        function(const Scalar x, const Scalar dx)
        {
            return AutoDiff(x, dx);
        }

        void
        operator +=(const Scalar& rhs)
        {
            val_ += rhs;
        }

        void
        operator +=(const AutoDiff& rhs)
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
        operator -=(const AutoDiff& rhs)
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
        operator *=(const AutoDiff& rhs)
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
        operator /=(const AutoDiff& rhs)
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
        AutoDiff(const Scalar x, const Scalar dx)
            : val_(x), der_(dx)
        {}

        Scalar val_;
        Scalar der_;
    };


    template <class Ostream, typename Scalar>
    Ostream&
    operator<<(Ostream& os, const AutoDiff<Scalar>& fw)
    {
        return fw.print(os);
    }

    template <typename Scalar>
    AutoDiff<Scalar>
    operator +(const AutoDiff<Scalar>& lhs,
               const AutoDiff<Scalar>& rhs)
    {
        AutoDiff<Scalar> ret = lhs;
        ret += rhs;

        return ret;
    }

    template <typename Scalar, typename T>
    AutoDiff<Scalar>
    operator +(const T                lhs,
               const AutoDiff<Scalar>& rhs)
    {
        AutoDiff<Scalar> ret = rhs;
        ret += Scalar(lhs);

        return ret;
    }

    template <typename Scalar, typename T>
    AutoDiff<Scalar>
    operator +(const AutoDiff<Scalar>& lhs,
               const T                rhs)
    {
        AutoDiff<Scalar> ret = lhs;
        ret += Scalar(rhs);

        return ret;
    }

    template <typename Scalar>
    AutoDiff<Scalar>
    operator -(const AutoDiff<Scalar>& lhs,
               const AutoDiff<Scalar>& rhs)
    {
        AutoDiff<Scalar> ret = lhs;
        ret -= rhs;

        return ret;
    }

    template <typename Scalar, typename T>
    AutoDiff<Scalar>
    operator -(const T                lhs,
               const AutoDiff<Scalar>& rhs)
    {
        return AutoDiff<Scalar>::function(Scalar(lhs) - rhs.val(), -rhs.der());
    }

    template <typename Scalar, typename T>
    AutoDiff<Scalar>
    operator -(const AutoDiff<Scalar>& lhs,
               const T                rhs)
    {
        AutoDiff<Scalar> ret = lhs;
        ret -= Scalar(rhs);

        return ret;
    }

    template <typename Scalar>
    AutoDiff<Scalar>
    operator *(const AutoDiff<Scalar>& lhs,
               const AutoDiff<Scalar>& rhs)
    {
        AutoDiff<Scalar> ret = lhs;
        ret *= rhs;

        return ret;
    }

    template <typename Scalar, typename T>
    AutoDiff<Scalar>
    operator *(const T                lhs,
               const AutoDiff<Scalar>& rhs)
    {
        AutoDiff<Scalar> ret = rhs;
        ret *= Scalar(lhs);

        return ret;
    }

    template <typename Scalar, typename T>
    AutoDiff<Scalar>
    operator *(const AutoDiff<Scalar>& lhs,
               const T                rhs)
    {
        AutoDiff<Scalar> ret = lhs;
        ret *= Scalar(rhs);

        return ret;
    }

    template <typename Scalar>
    AutoDiff<Scalar>
    operator /(const AutoDiff<Scalar>& lhs,
               const AutoDiff<Scalar>& rhs)
    {
        AutoDiff<Scalar> ret = lhs;
        ret /= rhs;

        return ret;
    }

    template <typename Scalar, typename T>
    AutoDiff<Scalar>
    operator /(const T                lhs,
               const AutoDiff<Scalar>& rhs)
    {
        Scalar a =  Scalar(lhs) / rhs.val();
        Scalar b = (-Scalar(lhs) / (rhs.val() * rhs.val())) * rhs.der();

        return AutoDiff<Scalar>::function(a, b);
    }

    template <typename Scalar, typename T>
    AutoDiff<Scalar>
    operator /(const AutoDiff<Scalar>& lhs,
               const T                rhs)
    {
        Scalar a = lhs.val() / Scalar(rhs);
        Scalar b = lhs.der() / Scalar(rhs);

        return AutoDiff<Scalar>::function(a, b);
    }

    template <typename Scalar>
    AutoDiff<Scalar>
    cos(const AutoDiff<Scalar>& x)
    {
        Scalar a = std::cos(x.val());
        Scalar b = -std::sin(x.val()) * x.der();

        return AutoDiff<Scalar>::function(a, b);
    }

    template <typename Scalar>
    AutoDiff<Scalar>
    sqrt(const AutoDiff<Scalar>& x)
    {
        Scalar a = std::sqrt(x.val());
        Scalar b = (Scalar(1.0) / (Scalar(2.0) * a)) * x.der();

        return AutoDiff<Scalar>::function(a, b);
    }

} // namespace Opm

namespace std {
    using Opm::cos;
    using Opm::sqrt;
}

#endif  /* OPM_AUTODIFF_HPP_HEADER */
