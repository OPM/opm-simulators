/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_AUTODIFFVEC_HEADER_INCLUDED
#define OPM_AUTODIFFVEC_HEADER_INCLUDED

#include <opm/autodiff/AutoDiff.hpp>
#include <Eigen/Eigen>
#include <Eigen/Sparse>

namespace AutoDiff
{

    template <typename Scalar>
    class ForwardVec
    {
    public:
        /// Underlying types for scalar vectors and jacobians.
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> V;
        typedef Eigen::SparseMatrix<Scalar> M;

        /// Named constructor pattern used here.
        static ForwardVec constant(const V& val)
        {
            return ForwardVec(val);
        }

        static ForwardVec variable(const V& val)
        {
            ForwardVec ret(val);

            // typedef Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D;
            // D ones = V::Ones(val.size()).matrix().asDiagonal();
            // ret.jac_ = ones;
            ret.jac_.reserve(Eigen::VectorXi::Constant(val.size(), 1));
            for (typename M::Index row = 0; row < val.size(); ++row) {
                ret.jac_.insert(row, row) = Scalar(1.0);
            }
            ret.jac_.makeCompressed();
            return ret;
        }

        static ForwardVec function(const V& val, const M& jac)
        {
            return ForwardVec(val, jac);
        }

        /// Operators.
        ForwardVec operator+(const ForwardVec& rhs)
        {
            return function(val_ + rhs.val_, jac_ + rhs.jac_);
        }

        /// Operators.
        ForwardVec operator-(const ForwardVec& rhs)
        {
            return function(val_ - rhs.val_, jac_ - rhs.jac_);
        }

        /// Operators.
        ForwardVec operator*(const ForwardVec& rhs)
        {
            typedef Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D;
            D D1 = val_.matrix().asDiagonal();
            D D2 = rhs.val_.matrix().asDiagonal();
            return function(val_ * rhs.val_, D2*jac_ + D1*rhs.jac_);
        }

        /// Operators.
        ForwardVec operator/(const ForwardVec& rhs)
        {
            typedef Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D;
            D D1 = val_.matrix().asDiagonal();
            D D2 = rhs.val_.matrix().asDiagonal();
            D D3 = std::pow(rhs.val_, -2).matrix().asDiagonal();
            return function(val_ / rhs.val_, D3 * (D2*jac_ - D1*rhs.jac_));
        }

        /// I/O.
        template <class Ostream>
        Ostream&
        print(Ostream& os) const
        {
            os << "val =\n" << val_ << "\n\njac =\n" << jac_ << "\n";

            return os;
        }

    private:
        explicit ForwardVec(const V& val)
            : val_(val), jac_(val.size(), val.size())
        {
        }

        ForwardVec(const V& val, const M& jac)
        : val_(val), jac_(jac)
        {
        }

        V val_;
        M jac_;
    };


    template <class Ostream, typename Scalar>
    Ostream&
    operator<<(Ostream& os, const ForwardVec<Scalar>& fw)
    {
        return fw.print(os);
    }



} // namespace Autodiff



#endif // OPM_AUTODIFFVEC_HEADER_INCLUDED
