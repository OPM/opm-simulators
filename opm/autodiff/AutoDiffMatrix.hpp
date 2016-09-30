/*
  Copyright 2014, 2015 SINTEF ICT, Applied Mathematics.

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

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/common/ErrorMacros.hpp>
#include <opm/autodiff/fastSparseOperations.hpp>
#include <vector>


namespace Opm
{

    /**
     * AutoDiffMatrix is a wrapper class that optimizes matrix operations.
     * Internally, an AutoDiffMatrix can be either Zero, Identity, Diagonal,
     * or Sparse, and we utilize this to perform faster matrix operations.
     */
    class AutoDiffMatrix
    {
    public:
        typedef std::vector<double> DiagRep;
        typedef Eigen::SparseMatrix<double> SparseRep;


        /**
         * Creates an empty zero matrix
         */
        AutoDiffMatrix()
            : type_(Zero),
              rows_(0),
              cols_(0),
              diag_(),
              sparse_()
        {
        }


        /**
         * Creates a zero matrix with num_rows x num_cols entries
         */
        AutoDiffMatrix(const int num_rows, const int num_cols)
            : type_(Zero),
              rows_(num_rows),
              cols_(num_cols),
              diag_(),
              sparse_()
        {
        }




        /**
         * Creates an identity matrix with num_rows_cols x num_rows_cols entries
         */
        static AutoDiffMatrix createIdentity(const int num_rows_cols)
        {
            return AutoDiffMatrix(Identity, num_rows_cols, num_rows_cols);
        }



        /**
         * Creates a diagonal matrix from an Eigen diagonal matrix
         */
        explicit AutoDiffMatrix(const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& d)
            : type_(Diagonal),
              rows_(d.rows()),
              cols_(d.cols()),
              diag_(d.diagonal().array().data(), d.diagonal().array().data() + d.rows()),
              sparse_()
        {
            assert(rows_ == cols_);
        }



        /**
         * Creates a sparse matrix from an Eigen sparse matrix
         */
        explicit AutoDiffMatrix(const Eigen::SparseMatrix<double>& s)
            : type_(Sparse),
              rows_(s.rows()),
              cols_(s.cols()),
              diag_(),
              sparse_(s)
        {
        }



        AutoDiffMatrix(const AutoDiffMatrix& other) = default;
        AutoDiffMatrix& operator=(const AutoDiffMatrix& other) = default;



        AutoDiffMatrix(AutoDiffMatrix&& other)
            : type_(Zero),
              rows_(0),
              cols_(0),
              diag_(),
              sparse_()
        {
            swap(other);
        }



        AutoDiffMatrix& operator=(AutoDiffMatrix&& other)
        {
            swap(other);
            return *this;
        }



        void swap(AutoDiffMatrix& other)
        {
            std::swap(type_, other.type_);
            std::swap(rows_, other.rows_);
            std::swap(cols_, other.cols_);
            diag_.swap(other.diag_);
            sparse_.swap(other.sparse_);
        }



        /**
         * Adds two AutoDiffMatrices. Internally, this function optimizes
         * the addition operation based on the structure of the matrix, e.g.,
         * adding two zero matrices will obviously yield a zero matrix, and
         * so on.
         */
        AutoDiffMatrix operator+(const AutoDiffMatrix& rhs) const
        {
            assert(rows_ == rhs.rows_);
            assert(cols_ == rhs.cols_);
            switch (type_) {
            case Zero:
                return rhs;
            case Identity:
                switch (rhs.type_) {
                case Zero:
                    return *this;
                case Identity:
                    return addII(*this, rhs);
                case Diagonal:
                    return addDI(rhs, *this);
                case Sparse:
                    return addSI(rhs, *this);
                default:
                    OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << rhs.type_);
                }
            case Diagonal:
                switch (rhs.type_) {
                case Zero:
                    return *this;
                case Identity:
                    return addDI(*this, rhs);
                case Diagonal:
                    return addDD(*this, rhs);
                case Sparse:
                    return addSD(rhs, *this);
                default:
                    OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << rhs.type_);
                }
            case Sparse:
                switch (rhs.type_) {
                case Zero:
                    return *this;
                case Identity:
                    return addSI(*this, rhs);
                case Diagonal:
                    return addSD(*this, rhs);
                case Sparse:
                    return addSS(*this, rhs);
                default:
                    OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << rhs.type_);
                }
            default:
                OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << rhs.type_);
            }
        }






        /**
         * Multiplies two AutoDiffMatrices. Internally, this function optimizes
         * the multiplication operation based on the structure of the matrix, e.g.,
         * multiplying M with a zero matrix will obviously yield a zero matrix.
         */
        AutoDiffMatrix operator*(const AutoDiffMatrix& rhs) const
        {
            assert(cols_ == rhs.rows_);
            switch (type_) {
            case Zero:
                return AutoDiffMatrix(rows_, rhs.cols_);
            case Identity:
                return rhs;
            case Diagonal:
                switch (rhs.type_) {
                case Zero:
                    return AutoDiffMatrix(rows_, rhs.cols_);
                case Identity:
                    return *this;
                case Diagonal:
                    return mulDD(*this, rhs);
                case Sparse:
                    return mulDS(*this, rhs);
                default:
                    OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << rhs.type_);
                }
            case Sparse:
                switch (rhs.type_) {
                case Zero:
                    return AutoDiffMatrix(rows_, rhs.cols_);
                case Identity:
                    return *this;
                case Diagonal:
                    return mulSD(*this, rhs);
                case Sparse:
                    return mulSS(*this, rhs);
                default:
                    OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << rhs.type_);
                }
            default:
                OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << rhs.type_);
            }
        }







        AutoDiffMatrix& operator+=(const AutoDiffMatrix& rhs)
        {
            if( type_ == Sparse && rhs.type_ == Sparse )
            {
                fastSparseAdd( sparse_, rhs.sparse_ );
            }
            else {
                *this = *this + rhs;
            }
            return *this;
        }






        AutoDiffMatrix& operator-=(const AutoDiffMatrix& rhs)
        {
            if( type_ == Sparse && rhs.type_ == Sparse )
            {
                fastSparseSubstract( sparse_, rhs.sparse_ );
            }
            else {
                *this = *this + (rhs * -1.0);
            }
            return *this;
        }






        /**
         * Multiplies an AutoDiffMatrix with a scalar. Optimizes internally
         * by exploiting that e.g., an identity matrix multiplied by a scalar x
         * yields a diagonal matrix with x the diagonal.
         */
        AutoDiffMatrix operator*(const double rhs) const
        {
            switch (type_) {
            case Zero:
                return *this;
            case Identity:
                {
                    AutoDiffMatrix retval(*this);
                    retval.type_ = Diagonal;
                    retval.diag_.assign(rows_, rhs);
                    return retval;
                }
            case Diagonal:
                {
                    AutoDiffMatrix retval(*this);
                    for (double& elem : retval.diag_) {
                        elem *= rhs;
                    }
                    return retval;
                }
            case Sparse:
                {
                    AutoDiffMatrix retval(*this);
                    retval.sparse_ *= rhs;
                    return retval;
                }
            default:
                OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << type_);
            }
        }






        /**
         * Divides an AutoDiffMatrix by a scalar. Optimizes internally
         * by exploiting that e.g., an identity matrix divided by a scalar x
         * yields a diagonal matrix with 1/x on the diagonal.
         */
        AutoDiffMatrix operator/(const double rhs) const
        {
            switch (type_) {
            case Zero:
                return *this;
            case Identity:
                {
                    AutoDiffMatrix retval(*this);
                    retval.type_ = Diagonal;
                    retval.diag_.assign(rows_, 1.0/rhs);
                    return retval;
                }
            case Diagonal:
                {
                    AutoDiffMatrix retval(*this);
                    for (double& elem : retval.diag_) {
                        elem /= rhs;
                    }
                    return retval;
                }
            case Sparse:
                {
                    AutoDiffMatrix retval(*this);
                    retval.sparse_ /= rhs;
                    return retval;
                }
            default:
                OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << type_);
            }
        }






        /**
         * Multiplies an AutoDiffMatrix with a vector. Optimizes internally
         * by exploiting that e.g., an identity matrix multiplied by a vector
         * yields the vector itself.
         */
        Eigen::VectorXd operator*(const Eigen::VectorXd& rhs) const
        {
            assert(cols_ == rhs.size());
            switch (type_) {
            case Zero:
                return Eigen::VectorXd::Zero(rows_);
            case Identity:
                return rhs;
            case Diagonal:
                {
                    const Eigen::VectorXd d = Eigen::Map<const Eigen::VectorXd>(diag_.data(), rows_);
                    return d.cwiseProduct(rhs);
                }
            case Sparse:
                return sparse_ * rhs;
            default:
                OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << type_);
            }
        }





        // Add identity to identity
        static AutoDiffMatrix addII(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == Identity);
            assert(rhs.type_ == Identity);
            AutoDiffMatrix retval;
            retval.type_ = Diagonal;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.diag_.assign(lhs.rows_, 2.0);
            return retval;
        }

        // Add diagonal to identity
        static AutoDiffMatrix addDI(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            static_cast<void>(rhs); // Silence release-mode warning.
            assert(lhs.type_ == Diagonal);
            assert(rhs.type_ == Identity);
            AutoDiffMatrix retval = lhs;
            for (int r = 0; r < lhs.rows_; ++r) {
                retval.diag_[r] += 1.0;
            }
            return retval;
        }

        // Add diagonal to diagonal
        static AutoDiffMatrix addDD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == Diagonal);
            assert(rhs.type_ == Diagonal);
            AutoDiffMatrix retval = lhs;
            for (int r = 0; r < lhs.rows_; ++r) {
                retval.diag_[r] += rhs.diag_[r];
            }
            return retval;
        }

        // Add sparse to identity
        static AutoDiffMatrix addSI(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            static_cast<void>(rhs); // Silence release-mode warning.
            assert(lhs.type_ == Sparse);
            assert(rhs.type_ == Identity);
            AutoDiffMatrix retval = lhs;
            retval.sparse_ += spdiag(Eigen::VectorXd::Ones(lhs.rows_));;
            return retval;
        }

        // Add sparse to diagonal
        static AutoDiffMatrix addSD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == Sparse);
            assert(rhs.type_ == Diagonal);
            AutoDiffMatrix retval = lhs;
            retval.sparse_ += spdiag(rhs.diag_);
            return retval;
        }

        // Add sparse to sparse
        static AutoDiffMatrix addSS(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == Sparse);
            assert(rhs.type_ == Sparse);
            AutoDiffMatrix retval = lhs;
            retval.sparse_ += rhs.sparse_;
            return retval;
        }




        // Multiply diagonal with diagonal
        static AutoDiffMatrix mulDD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == Diagonal);
            assert(rhs.type_ == Diagonal);
            AutoDiffMatrix retval = lhs;
            for (int r = 0; r < lhs.rows_; ++r) {
                retval.diag_[r] *= rhs.diag_[r];
            }
            return retval;
        }

        // Multiply diagonal with sparse
        static AutoDiffMatrix mulDS(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == Diagonal);
            assert(rhs.type_ == Sparse);
            AutoDiffMatrix retval;
            retval.type_ = Sparse;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            fastDiagSparseProduct(lhs.diag_, rhs.sparse_, retval.sparse_);
            return retval;
        }

        // Multiply sparse with diagonal
        static AutoDiffMatrix mulSD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == Sparse);
            assert(rhs.type_ == Diagonal);
            AutoDiffMatrix retval;
            retval.type_ = Sparse;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            fastSparseDiagProduct(lhs.sparse_, rhs.diag_, retval.sparse_);
            return retval;
        }

        // Multiply sparse with sparse
        static AutoDiffMatrix mulSS(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == Sparse);
            assert(rhs.type_ == Sparse);
            AutoDiffMatrix retval;
            retval.type_ = Sparse;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            fastSparseProduct(lhs.sparse_, rhs.sparse_, retval.sparse_);
            return retval;
        }




        /**
         * Converts the AutoDiffMatrix to an Eigen SparseMatrix.This might be
         * an expensive operation to perform for e.g., an identity matrix or a
         * diagonal matrix.
         */
        template<class Scalar, int Options, class Index>
        void toSparse(Eigen::SparseMatrix<Scalar, Options, Index>& s) const
        {
            switch (type_) {
            case Zero:
                s = Eigen::SparseMatrix<Scalar, Options, Index>(rows_, cols_);
                return;
            case Identity:
                s = spdiag(Eigen::VectorXd::Ones(rows_));
                return;
            case Diagonal:
                s = spdiag(diag_);
                return;
            case Sparse:
                s = sparse_;
                return;
            }
        }


        /**
         * Converts the AutoDiffMatrix to an Eigen SparseMatrix.This might be
         * an expensive operation to perform for e.g., an identity matrix or a
         * diagonal matrix.
         */
        template<class Scalar, int Options, class Index>
        void assign(const Eigen::SparseMatrix<Scalar, Options, Index>& s)
        {
            (*this) = AutoDiffMatrix( s.rows(), s.cols() );
            type_ = Sparse;
            sparse_ = s;
        }


        template<class Scalar, int Options, class Index>
        void assign(const Eigen::SparseMatrix<Scalar, Options, Index>&& s)
        {
            (*this) = AutoDiffMatrix( s.rows(), s.cols() );
            type_ = Sparse;
            sparse_ = std::move(s);
        }


        /**
         * Returns number of rows in the matrix
         */
        int rows() const
        {
            return rows_;
        }



        /**
         * Returns number of columns in the matrix
         */
        int cols() const
        {
            return cols_;
        }



        /**
         * Returns number of non-zero elements in the matrix. Optimizes internally
         * by exploiting that e.g., an n*n identity matrix has n non-zeros.
         * Note that an n*n diagonal matrix is defined to have n non-zeros, even though
         * several diagonal elements might be 0.0.
         */
        int nonZeros() const
        {
            switch (type_) {
            case Zero:
                return 0;
            case Identity:
                return rows_;
            case Diagonal:
                return rows_;
            case Sparse:
                return sparse_.nonZeros();
            default:
                OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << type_);
            }
        }




        /**
         * Returns element (row, col) in the matrix
         */
        double coeff(const int row, const int col) const
        {
            switch (type_) {
            case Zero:
                return 0.0;
            case Identity:
                return (row == col) ? 1.0 : 0.0;
            case Diagonal:
                return (row == col) ? diag_[row] : 0.0;
            case Sparse:
                return sparse_.coeff(row, col);
            default:
                OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << type_);
            }
        }





        /**
         * Returns the sparse representation of this matrix. Note that this might
         * be an expensive operation to perform if the internal structure is not
         * sparse.
         */
        const SparseRep& getSparse() const {
            if (type_ != Sparse) {
                /**
                 * If we are not a sparse matrix, our internal variable sparse_
                 * is undefined, and hence changing it so that it happens to be
                 * a sparse representation of our true data does not change our
                 * true data, and hence justifies that we do not really violate
                 * the const qualifier.
                 */
                SparseRep& mutable_sparse = const_cast<SparseRep&>(sparse_);
                toSparse(mutable_sparse);
            }
            return sparse_;
        }


    private:
        enum AudoDiffMatrixType { Zero, Identity, Diagonal, Sparse };

        AudoDiffMatrixType type_;  //<  Type of matrix
        int rows_;                 //<  Number of rows
        int cols_;                 //<  Number of columns
        DiagRep diag_;             //<  Diagonal representation (only if type==Diagonal)
        SparseRep sparse_;         //<  Sparse representation (only if type==Sparse)



        /**
         * Constructor which sets all members
         */
        AutoDiffMatrix(AudoDiffMatrixType type, int rows_arg, int cols_arg,
                DiagRep diag=DiagRep(), SparseRep sparse=SparseRep())
            : type_(type),
              rows_(rows_arg),
              cols_(cols_arg),
              diag_(diag),
              sparse_(sparse)
        {
        }





        /**
         * Creates a sparse diagonal matrix from d.
         * Typical use is to convert a standard vector to an
         * Eigen sparse matrix.
         */
        template <class V>
        static inline
        SparseRep
        spdiag(const V& d)
        {
            const int n = d.size();
            SparseRep mat(n, n);
            mat.reserve(Eigen::ArrayXi::Ones(n, 1));
            for (SparseRep::Index i = 0; i < n; ++i) {
                if (d[i] != 0.0) {
                    mat.insert(i, i) = d[i];
                }
            }

            return mat;
        }

    };



    /**
     * Utility function to lessen code changes required elsewhere.
     */
    inline void fastSparseProduct(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs, AutoDiffMatrix& res)
    {
        res = lhs * rhs;
    }


    /**
     * Utility function to lessen code changes required elsewhere.
     */
    inline void fastSparseProduct(const Eigen::SparseMatrix<double>& lhs, const AutoDiffMatrix& rhs, AutoDiffMatrix& res)
    {
        res = AutoDiffMatrix(lhs) * rhs;
    }


    /**
     * Multiplies an Eigen sparse matrix with an AutoDiffMatrix.
     */
    inline AutoDiffMatrix operator*(const Eigen::SparseMatrix<double>& lhs, const AutoDiffMatrix& rhs)
    {
        AutoDiffMatrix retval;
        fastSparseProduct(lhs, rhs, retval);
        return retval;
    }

} // namespace Opm


#endif // OPM_AUTODIFFMATRIX_HEADER_INCLUDED
