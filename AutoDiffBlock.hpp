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

#ifndef OPM_AUTODIFFBLOCK_HEADER_INCLUDED
#define OPM_AUTODIFFBLOCK_HEADER_INCLUDED

#include "AutoDiff.hpp"
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <vector>
#include <cassert>

namespace AutoDiff
{

    template <typename Scalar>
    class ForwardBlock
    {
    public:
        /// Underlying types for scalar vectors and jacobians.
        typedef Eigen::Array<Scalar, Eigen::Dynamic, 1> V;
        typedef Eigen::SparseMatrix<Scalar> M;

        /// Named constructor pattern used here.
        static ForwardBlock constant(const int index, const V& val, const std::vector<int>& blocksizes)
        {
            std::vector<M> jac;
            const int num_elem = val.size();
            const int num_blocks = blocksizes.size();
            // For constants, all jacobian blocks are zero.
            jac.resize(num_blocks);
            for (int i = 0; i < num_blocks; ++i) {
                jac[i] = M(num_elem, blocksizes[i]);
            }
            return ForwardBlock(index, val, jac);
        }

        static ForwardBlock variable(const int index, const V& val, const std::vector<int>& blocksizes)
        {

            std::vector<M> jac;
            const int num_elem = val.size();
            const int num_blocks = blocksizes.size();
            // First, set all jacobian blocks to zero...
            jac.resize(num_blocks);
            for (int i = 0; i < num_blocks; ++i) {
                jac[i] = M(num_elem, blocksizes[i]);
            }
            // ... then set the one corrresponding to this variable to identity.
            jac[index].reserve(Eigen::VectorXi::Constant(val.size(), 1));
            for (typename M::Index row = 0; row < val.size(); ++row) {
                jac[index].insert(row, row) = Scalar(1.0);
            }
            return ForwardBlock(index, val, jac);
        }

        static ForwardBlock function(const int index, const V& val, const std::vector<M>& jac)
        {
            return ForwardBlock(index, val, jac);
        }

#if 0
        /// Operators.
        ForwardBlock operator+(const ForwardBlock& rhs)
        {
            return function(val_ + rhs.val_, jac_ + rhs.jac_);
        }

        /// Operators.
        ForwardBlock operator-(const ForwardBlock& rhs)
        {
            return function(val_ - rhs.val_, jac_ - rhs.jac_);
        }

        /// Operators.
        ForwardBlock operator*(const ForwardBlock& rhs)
        {
            typedef Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D;
            D D1 = val_.matrix().asDiagonal();
            D D2 = rhs.val_.matrix().asDiagonal();
            return function(val_ * rhs.val_, D2*jac_ + D1*rhs.jac_);
        }

        /// Operators.
        ForwardBlock operator/(const ForwardBlock& rhs)
        {
            typedef Eigen::DiagonalMatrix<Scalar, Eigen::Dynamic> D;
            D D1 = val_.matrix().asDiagonal();
            D D2 = rhs.val_.matrix().asDiagonal();
            D D3 = std::pow(rhs.val_, -2).matrix().asDiagonal();
            return function(val_ / rhs.val_, D3 * (D2*jac_ - D1*rhs.jac_));
        }
#endif
        /// I/O.
        template <class Ostream>
        Ostream&
        print(Ostream& os) const
        {
            int num_blocks = jac_.size();
            os << "Variable index " << index_ << ":\nValue =\n" << val_ << "\n\nJacobian =\n";
            for (int i = 0; i < num_blocks; ++i) {
                os << "Sub Jacobian #" << i << '\n' << jac_[i] << "\n";
            }
            return os;
        }

    private:
        ForwardBlock(const int index,
                     const V& val,
                     const std::vector<M>& jac)
            : index_(index), val_(val), jac_(jac)
        {
            const int num_elem = val_.size();
            const int num_blocks = jac_.size();
            assert(index_ < num_blocks);
            assert(num_elem == jac_[index].cols());
            for (int block = 0; block < num_blocks; ++block) {
                assert(num_elem == jac_[block].rows());
            }
        }

        int index_;
        V val_;
        std::vector<M> jac_;
    };


    template <class Ostream, typename Scalar>
    Ostream&
    operator<<(Ostream& os, const ForwardBlock<Scalar>& fw)
    {
        return fw.print(os);
    }



} // namespace Autodiff



#endif // OPM_AUTODIFFBLOCK_HEADER_INCLUDED
