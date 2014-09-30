/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_AUTODIFFMATRIX_HEADER_INCLUDED
#define OPM_AUTODIFFMATRIX_HEADER_INCLUDED

#include <opm/core/utility/platform_dependent/disable_warnings.h>

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include <opm/core/utility/platform_dependent/reenable_warnings.h>


namespace Opm
{

    /// Implementation details for class AutoDiffMatrix.
    namespace AutoDiffMatrixDetail
    {


        class Zero;
        class Identity;
        class Diagonal;
        class Sparse;
        typedef std::shared_ptr<Zero> ZeroMat;
        typedef std::shared_ptr<Identity> IdentityMat;
        typedef std::shared_ptr<Diagonal> DiagonalMat;
        typedef std::shared_ptr<Sparse> SparseMat;






        class Interface
        {
        public:
            typedef std::shared_ptr<Interface> Mat;

            virtual ~Interface()
            {
            }

            virtual Mat operator+(const Mat& rhs) = 0;
            virtual Mat addIdentity(const IdentityMat& rhs) = 0;
            virtual Mat addDiagonal(const DiagonalMat& rhs) = 0;
            virtual Mat addSparse(const SparsMate& rhs) = 0;

            virtual Mat operator*(const Mat& rhs) = 0;
            virtual Mat leftMulDiagonal(const DiagonalMat& rhs) = 0;
            virtual Mat leftMulSparse(const SparseMat& rhs) = 0;
        };






        class Zero : public Interface
        {
        public:
            virtual Mat operator+(const Mat& rhs)
            {
                return rhs;
            }

            virtual Mat addIdentity(const IdentityMat& rhs)
            {
                return rhs;
            }

            virtual Mat addDiagonal(const DiagonalMat& rhs)
            {
                return rhs;
            }

            virtual Mat addSparse(const SparseMat& rhs)
            {
                return rhs;
            }


            virtual Mat operator*(const Mat& rhs)
            {
                return std::make_shared<Zero>();
            }

            virtual Mat leftMulDiagonal(const DiagonalMat& rhs)
            {
                return std::make_shared<Zero>();
            }

            virtual Mat leftMulSparse(const SparseMat& rhs)
            {
                return std::make_shared<Zero>();
            }
        };





        class Identity : public Interface
        {
        public:
            virtual Mat operator+(const Mat& rhs)
            {
                return rhs->addIdentity(*this);
            }

            virtual Mat addIdentity(const IdentityMat& rhs)
            {
                return rhs;
            }

            virtual Mat addDiagonal(const DiagonalMat& rhs)
            {
                return rhs;
            }

            virtual Mat addSparse(const SparseMat& rhs)
            {
                return rhs;
            }


            virtual Mat operator*(const Mat& rhs)
            {
                return std::make_shared<Zero>();
            }

            virtual Mat leftMulDiagonal(const DiagonalMat& rhs)
            {
                return std::make_shared<Zero>();
            }

            virtual Mat leftMulSparse(const SparseMat& rhs)
            {
                return std::make_shared<Zero>();
            }
        };





        class Diagonal : public Interface
        {
        public:
            virtual Mat operator+(const Mat& rhs)
            {
                return (*rhs) + (*this);
            }
            operator+(const IdentityMat& rhs)
            {
                // TODO return Diagnonal(...);
            }
            operator+(const DiagonalMat& rhs)
            {
                // TODO return Diagonal(...);
            }
            operator+(const SparseMat& rhs)
            {
                // TODO return Sparse(...);
            }

            virtual Mat operator*(const Mat& rhs)
            {
                return (*rhs) * (;
            }
            Mat operator*(const IdentityMat& rhs)
            {
                return *this;
            }
            Mat operator*(const DiagonalMat& rhs)
            {
                // TODO return Diagonal(...);
            }
            Mat operator*(const SparseMat& rhs)
            {
                // TODO return Sparse(...);
            }
        };





        class Sparse : public Interface
        {
            virtual Mat operator+(const Mat& rhs)
            {
                return (*rhs) + (*this);
            }
            operator+(const IdentityMat& rhs)
            {
                // TODO return Sparse(...);
            }
            operator+(const DiagonalMat& rhs)
            {
                // TODO return Sparse(...);
            }
            operator+(const SparseMat& rhs)
            {
                // TODO return Sparse(...);
            }

            virtual Mat operator*(const Mat& rhs)
            {
                return rhs;
            }
            Mat operator*(const IdentityMat& rhs)
            {
                return *this;
            }
            Mat operator*(const DiagonalMat& rhs)
            {
                // TODO return Sparse(...);
            }
            Mat operator*(const SparseMat& rhs)
            {
                // TODO return Sparse(...);
            }
        };

    } // namespace AutoDiffMatrixDetail
} // namespace Opm


#endif // OPM_AUTODIFFMATRIX_HEADER_INCLUDED
