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

#include "AutoDiff.hpp"
#include <Eigen/Eigen>
#include <Eigen/Sparse>

namespace AutoDiff
{

    template <typename Scalar>
    class ForwardVec
    {
    public:
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> V;
        typedef Eigen::SparseMatrix<Scalar> M;

        static ForwardVec constant(const V& val)
        {
            return ForwardVec(val);
        }

        static ForwardVec variable(const V& val)
        {
            ForwardVec ret(val);

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

        explicit ForwardVec(const V& val)
            : val_(val), jac_(val.size(), val.size())
        {
        }

        ForwardVec(const V& val, const M& jac)
        : val_(val), jac_(jac)
        {
        }

        template <class Ostream>
        Ostream&
        print(Ostream& os) const
        {
            os << "val =\n" << val_ << "\n\njac =\n" << jac_ << "\n";

            return os;
        }

        
    private:
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
