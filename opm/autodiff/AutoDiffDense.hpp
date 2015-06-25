/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Statoil ASA.

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

#ifndef OPM_AUTODIFFDENSE_HEADER_INCLUDED
#define OPM_AUTODIFFDENSE_HEADER_INCLUDED


#include <opm/core/utility/platform_dependent/disable_warnings.h>

#include <Eigen/Eigen>

#include <opm/core/utility/platform_dependent/reenable_warnings.h>


namespace Opm
{

    /// A simple class for forward-mode automatic differentiation.
    ///
    /// The class represents a single value and a small, dense vector
    /// of derivatives.  Only basic arithmetic operators are
    /// implemented for it.
    template <typename Scalar, int NumDerivs>
    class AutoDiffDense
    {
    public:
        typedef Scalar Value;
        typedef Eigen::Array<Scalar, NumDerivs, 1> Derivative;

        /// Create an AutoDiffDense object representing a constant, that
        /// is, its derivative is zero.
        static AutoDiffDense
        constant(const Scalar x)
        {
            return function(x, Derivative::Zero());
        }

        /// Create an AutoDiffDense object representing a primary variable,
        /// that is, its derivative is one.
        static AutoDiffDense
        variable(const int index, const Scalar x)
        {
            Derivative dx = Derivative::Zero();
            dx(index) = 1.0;
            return function(x, dx);
        }

        /// Create an AutoDiffDense object representing a function value
        /// and its derivative.
        static AutoDiffDense
        function(const Scalar x, const Derivative& dx)
        {
            return AutoDiffDense(x, dx);
        }

        void
        operator +=(const Scalar& rhs)
        {
            val_ += rhs;
        }

        void
        operator +=(const AutoDiffDense& rhs)
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
        operator -=(const AutoDiffDense& rhs)
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
        operator *=(const AutoDiffDense& rhs)
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
        operator /=(const AutoDiffDense& rhs)
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
        const Derivative& der() const { return der_; }

    private:
        AutoDiffDense(const Scalar x, const Derivative& dx)
            : val_(x), der_(dx)
        {}

        Scalar val_;
        Derivative der_;
    };






    template <class Ostream, typename Scalar, int NumDerivs>
    Ostream&
    operator<<(Ostream& os, const AutoDiffDense<Scalar, NumDerivs>& fw)
    {
        return fw.print(os);
    }

    template <typename Scalar, int NumDerivs>
    AutoDiffDense<Scalar, NumDerivs>
    operator +(const AutoDiffDense<Scalar, NumDerivs>& lhs,
               const AutoDiffDense<Scalar, NumDerivs>& rhs)
    {
        AutoDiffDense<Scalar, NumDerivs> ret = lhs;
        ret += rhs;

        return ret;
    }

    template <typename Scalar, int NumDerivs, typename T>
    AutoDiffDense<Scalar, NumDerivs>
    operator +(const T                lhs,
               const AutoDiffDense<Scalar, NumDerivs>& rhs)
    {
        AutoDiffDense<Scalar, NumDerivs> ret = rhs;
        ret += Scalar(lhs);

        return ret;
    }

    template <typename Scalar, int NumDerivs, typename T>
    AutoDiffDense<Scalar, NumDerivs>
    operator +(const AutoDiffDense<Scalar, NumDerivs>& lhs,
               const T                rhs)
    {
        AutoDiffDense<Scalar, NumDerivs> ret = lhs;
        ret += Scalar(rhs);

        return ret;
    }

    template <typename Scalar, int NumDerivs>
    AutoDiffDense<Scalar, NumDerivs>
    operator -(const AutoDiffDense<Scalar, NumDerivs>& lhs,
               const AutoDiffDense<Scalar, NumDerivs>& rhs)
    {
        AutoDiffDense<Scalar, NumDerivs> ret = lhs;
        ret -= rhs;

        return ret;
    }

    template <typename Scalar, int NumDerivs, typename T>
    AutoDiffDense<Scalar, NumDerivs>
    operator -(const T                lhs,
               const AutoDiffDense<Scalar, NumDerivs>& rhs)
    {
        return AutoDiffDense<Scalar, NumDerivs>::function(Scalar(lhs) - rhs.val(), -rhs.der());
    }

    template <typename Scalar, int NumDerivs, typename T>
    AutoDiffDense<Scalar, NumDerivs>
    operator -(const AutoDiffDense<Scalar, NumDerivs>& lhs,
               const T                rhs)
    {
        AutoDiffDense<Scalar, NumDerivs> ret = lhs;
        ret -= Scalar(rhs);

        return ret;
    }

    template <typename Scalar, int NumDerivs>
    AutoDiffDense<Scalar, NumDerivs>
    operator *(const AutoDiffDense<Scalar, NumDerivs>& lhs,
               const AutoDiffDense<Scalar, NumDerivs>& rhs)
    {
        AutoDiffDense<Scalar, NumDerivs> ret = lhs;
        ret *= rhs;

        return ret;
    }

    template <typename Scalar, int NumDerivs, typename T>
    AutoDiffDense<Scalar, NumDerivs>
    operator *(const T                lhs,
               const AutoDiffDense<Scalar, NumDerivs>& rhs)
    {
        AutoDiffDense<Scalar, NumDerivs> ret = rhs;
        ret *= Scalar(lhs);

        return ret;
    }

    template <typename Scalar, int NumDerivs, typename T>
    AutoDiffDense<Scalar, NumDerivs>
    operator *(const AutoDiffDense<Scalar, NumDerivs>& lhs,
               const T                rhs)
    {
        AutoDiffDense<Scalar, NumDerivs> ret = lhs;
        ret *= Scalar(rhs);

        return ret;
    }

    template <typename Scalar, int NumDerivs>
    AutoDiffDense<Scalar, NumDerivs>
    operator /(const AutoDiffDense<Scalar, NumDerivs>& lhs,
               const AutoDiffDense<Scalar, NumDerivs>& rhs)
    {
        AutoDiffDense<Scalar, NumDerivs> ret = lhs;
        ret /= rhs;

        return ret;
    }

    template <typename Scalar, int NumDerivs, typename T>
    AutoDiffDense<Scalar, NumDerivs>
    operator /(const T                lhs,
               const AutoDiffDense<Scalar, NumDerivs>& rhs)
    {
        Scalar a =  Scalar(lhs) / rhs.val();
        typedef typename AutoDiffDense<Scalar, NumDerivs>::Derivative Derivative;
        Derivative b = (-Scalar(lhs) / (rhs.val() * rhs.val())) * rhs.der();

        return AutoDiffDense<Scalar, NumDerivs>::function(a, b);
    }

    template <typename Scalar, int NumDerivs, typename T>
    AutoDiffDense<Scalar, NumDerivs>
    operator /(const AutoDiffDense<Scalar, NumDerivs>& lhs,
               const T                rhs)
    {
        Scalar a = lhs.val() / Scalar(rhs);
        auto b = lhs.der() / Scalar(rhs);

        return AutoDiffDense<Scalar, NumDerivs>::function(a, b);
    }

} // namespace Opm


#endif // OPM_AUTODIFFDENSE_HEADER_INCLUDED
