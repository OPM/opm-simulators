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

#include <opm/core/utility/platform_dependent/disable_warnings.h>

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include <opm/core/utility/platform_dependent/reenable_warnings.h>

#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/autodiff/fastSparseProduct.hpp>
#include <vector>


namespace Opm
{

    class AutoDiffMatrix
    {
    public:
        typedef std::vector<double> DiagRep;
        typedef Eigen::SparseMatrix<double> SparseRep;


        AutoDiffMatrix()
            : type_(Zero),
              rows_(0),
              cols_(0),
              diag_(),
              sparse_()
        {
        }



        AutoDiffMatrix(const int num_rows, const int num_cols)
            : type_(Zero),
              rows_(num_rows),
              cols_(num_cols),
              diag_(),
              sparse_()
        {
        }



        enum CreationType { ZeroMatrix, IdentityMatrix };


        AutoDiffMatrix(const CreationType t, const int num_rows)
            : type_(t == ZeroMatrix ? Zero : Identity),
              rows_(num_rows),
              cols_(num_rows),
              diag_(),
              sparse_()
        {
        }



        explicit AutoDiffMatrix(const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& d)
            : type_(Diagonal),
              rows_(d.rows()),
              cols_(d.cols()),
              diag_(d.diagonal().array().data(), d.diagonal().array().data() + d.rows()),
              sparse_()
        {
            assert(rows_ == cols_);
        }



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
            *this = *this + rhs;
            return *this;
        }






        AutoDiffMatrix& operator-=(const AutoDiffMatrix& rhs)
        {
            *this = *this + (rhs * -1.0);
            return *this;
        }






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

        static AutoDiffMatrix addSI(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            static_cast<void>(rhs); // Silence release-mode warning.
            assert(lhs.type_ == Sparse);
            assert(rhs.type_ == Identity);
            AutoDiffMatrix retval = lhs;
            retval.sparse_ += spdiag(Eigen::VectorXd::Ones(lhs.rows_));;
            return retval;
        }

        static AutoDiffMatrix addSD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == Sparse);
            assert(rhs.type_ == Diagonal);
            AutoDiffMatrix retval = lhs;
            retval.sparse_ += spdiag(rhs.diag_);
            return retval;
        }

        static AutoDiffMatrix addSS(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == Sparse);
            assert(rhs.type_ == Sparse);
            AutoDiffMatrix retval = lhs;
            retval.sparse_ += rhs.sparse_;
            return retval;
        }





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


        int rows() const
        {
            return rows_;
        }

        int cols() const
        {
            return cols_;
        }

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
        AudoDiffMatrixType type_;
        int rows_;
        int cols_;
        DiagRep diag_;
        SparseRep sparse_;

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




    inline void fastSparseProduct(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs, AutoDiffMatrix& res)
    {
        res = lhs * rhs;
    }


    inline void fastSparseProduct(const Eigen::SparseMatrix<double>& lhs, const AutoDiffMatrix& rhs, AutoDiffMatrix& res)
    {
        res = AutoDiffMatrix(lhs) * rhs;
    }


    inline AutoDiffMatrix operator*(const Eigen::SparseMatrix<double>& lhs, const AutoDiffMatrix& rhs)
    {
        AutoDiffMatrix retval;
        fastSparseProduct(lhs, rhs, retval);
        return retval;
    }

} // namespace Opm


#endif // OPM_AUTODIFFMATRIX_HEADER_INCLUDED
