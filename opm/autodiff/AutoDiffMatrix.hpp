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
        typedef std::vector<double> Diag;
        typedef Eigen::SparseMatrix<double> Sparse;

        AutoDiffMatrix()
            : type_(Z),
              rows_(0),
              cols_(0)
        {
        }



        AutoDiffMatrix(const int num_rows, const int num_cols)
            : type_(Z),
              rows_(num_rows),
              cols_(num_cols)
        {
        }



        enum CreationType { ZeroMatrix, IdentityMatrix };


        AutoDiffMatrix(const CreationType t, const int num_rows)
            : type_(t == ZeroMatrix ? Z : I),
              rows_(num_rows),
              cols_(num_rows)
        {
        }



        explicit AutoDiffMatrix(const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& d)
            : type_(D),
              rows_(d.rows()),
              cols_(d.cols()),
              diag_(d.diagonal().array().data(), d.diagonal().array().data() + d.rows())
        {
        }



        explicit AutoDiffMatrix(const Eigen::SparseMatrix<double>& s)
            : type_(S),
              rows_(s.rows()),
              cols_(s.cols())
        {
            sparse_[0] = s;
        }



        AutoDiffMatrix(const AutoDiffMatrix& other)
        {
            *this = other;
        }

        AutoDiffMatrix& operator=(const AutoDiffMatrix& other)
        {
            type_ = other.type_;
            rows_ = other.rows_;
            cols_ = other.cols_;
            switch(type_) {
                case D:
                    diag_ = other.diag_;
                    break;
                case S:
                    sparse_[0] = other.sparse_[0];
                    break;
                default:
                    break;
            }
        }



        AutoDiffMatrix(AutoDiffMatrix&& other)
            : AutoDiffMatrix()
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
            std::swap(sparse_, other.sparse_);
        }



        AutoDiffMatrix operator+(const AutoDiffMatrix& rhs) const
        {
            assert(rows_ == rhs.rows_);
            assert(cols_ == rhs.cols_);
            switch (type_) {
            case Z:
                return rhs;
            case I:
                switch (rhs.type_) {
                case Z:
                    return *this;
                case I:
                    return addII(*this, rhs);
                case D:
                    return rhs + (*this);
                case S:
                    return rhs + (*this);
                default:
                    OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << rhs.type_);
                }
            case D:
                switch (rhs.type_) {
                case Z:
                    return *this;
                case I:
                    return addDI(*this, rhs);
                case D:
                    return addDD(*this, rhs);
                case S:
                    return rhs + (*this);
                default:
                    OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << rhs.type_);
                }
            case S:
                switch (rhs.type_) {
                case Z:
                    return *this;
                case I:
                    return addSI(*this, rhs);
                case D:
                    return addSD(*this, rhs);
                case S:
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
            case Z:
                return AutoDiffMatrix(rows_, rhs.cols_);
            case I:
                switch (rhs.type_) {
                case Z:
                    return rhs;
                case I:
                    return rhs;
                case D:
                    return rhs;
                case S:
                    return rhs;
                default:
                    OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << rhs.type_);
                }
            case D:
                switch (rhs.type_) {
                case Z:
                    return AutoDiffMatrix(rows_, rhs.cols_);
                case I:
                    return *this;
                case D:
                    return mulDD(*this, rhs);
                case S:
                    return mulDS(*this, rhs);
                default:
                    OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << rhs.type_);
                }
            case S:
                switch (rhs.type_) {
                case Z:
                    return AutoDiffMatrix(rows_, rhs.cols_);
                case I:
                    return *this;
                case D:
                    return mulSD(*this, rhs);
                case S:
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
            *this = *this + rhs * -1.0;
            return *this;
        }






        AutoDiffMatrix operator*(const double rhs) const
        {
            switch (type_) {
            case Z:
                return *this;
            case I:
                {
                    AutoDiffMatrix retval(*this);
                    retval.type_ = D;
                    retval.diag_.assign(rows_, rhs);
                    return retval;
                }
            case D:
                {
                    AutoDiffMatrix retval(*this);
                    for (double& elem : retval.diag_) {
                        elem *= rhs;
                    }
                    return retval;
                }
            case S:
                {
                    AutoDiffMatrix retval(*this);
                    retval.sparse_[0] *= rhs;
                    return retval;
                }
            default:
                OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << type_);
            }
        }






        AutoDiffMatrix operator/(const double rhs) const
        {
            switch (type_) {
            case Z:
                return *this;
            case I:
                {
                    AutoDiffMatrix retval(*this);
                    retval.type_ = D;
                    retval.diag_.assign(rows_, 1.0/rhs);
                    return retval;
                }
            case D:
                {
                    AutoDiffMatrix retval(*this);
                    for (double& elem : retval.diag_) {
                        elem /= rhs;
                    }
                    return retval;
                }
            case S:
                {
                    AutoDiffMatrix retval(*this);
                    retval.sparse_[0] /= rhs;
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
            case Z:
                return Eigen::VectorXd::Zero(rows_);
            case I:
                return rhs;
            case D:
                return Eigen::Map<const Eigen::VectorXd>(diag_.data(), rows_) * rhs;
            case S:
                return sparse_[0] * rhs;
            default:
                OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << type_);
            }
        }






        static AutoDiffMatrix addII(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == I);
            assert(rhs.type_ == I);
            AutoDiffMatrix retval;
            retval.type_ = D;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.diag_.assign(lhs.rows_, 2.0);
            return retval;
        }

        static AutoDiffMatrix addDI(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            static_cast<void>(rhs); // Silence release-mode warning.
            assert(lhs.type_ == D);
            assert(rhs.type_ == I);
            AutoDiffMatrix retval = lhs;
            for (int r = 0; r < lhs.rows_; ++r) {
                retval.diag_[r] += 1.0;
            }
            return retval;
        }

        static AutoDiffMatrix addDD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == D);
            assert(rhs.type_ == D);
            AutoDiffMatrix retval = lhs;
            for (int r = 0; r < lhs.rows_; ++r) {
                retval.diag_[r] += rhs.diag_[r];
            }
            return retval;
        }

        static AutoDiffMatrix addSI(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == S);
            assert(rhs.type_ == I);
            AutoDiffMatrix retval;
            retval.type_ = S;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.sparse_[0] = lhs.sparse_[0];
            retval.sparse_[0] += spdiag(Eigen::VectorXd::Ones(lhs.rows_));
            return retval;
        }

        static AutoDiffMatrix addSD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == S);
            assert(rhs.type_ == D);
            AutoDiffMatrix retval;
            retval.type_ = S;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.sparse_[0] = lhs.sparse_[0];
            retval.sparse_[0] += spdiag(rhs.diag_);
            return retval;
        }

        static AutoDiffMatrix addSS(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == S);
            assert(rhs.type_ == S);
            AutoDiffMatrix retval = lhs;
            retval.sparse_[0] += rhs.sparse_[0];
            return retval;
        }





        static AutoDiffMatrix mulDD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == D);
            assert(rhs.type_ == D);
            AutoDiffMatrix retval = lhs;
            for (int r = 0; r < lhs.rows_; ++r) {
                retval.diag_[r] *= rhs.diag_[r];
            }
            return retval;
        }

        static AutoDiffMatrix mulDS(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == D);
            assert(rhs.type_ == S);
            AutoDiffMatrix retval;
            retval.type_ = S;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.sparse_[0] = std::move(fastDiagSparseProduct(lhs.diag_, rhs.sparse_[0]));
            return retval;
        }

        static AutoDiffMatrix mulSD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == S);
            assert(rhs.type_ == D);
            AutoDiffMatrix retval;
            retval.type_ = S;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.sparse_[0] = std::move(fastSparseDiagProduct(lhs.sparse_[0], rhs.diag_));
            return retval;
        }

        static AutoDiffMatrix mulSS(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == S);
            assert(rhs.type_ == S);
            AutoDiffMatrix retval;
            retval.type_ = S;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.sparse_[0] = std::move(fastSparseProduct<Sparse>(lhs.sparse_[0], rhs.sparse_[0]));
            return retval;
        }


        template<class Scalar, int Options, class Index>
        void toSparse(Eigen::SparseMatrix<Scalar, Options, Index>& s) const
        {
            switch (type_) {
            case Z:
                s = Eigen::SparseMatrix<Scalar, Options, Index>(rows_, cols_);
                return;
            case I:
                s = spdiag(Eigen::VectorXd::Ones(rows_));
                return;
            case D:
                s = spdiag(diag_);
                return;
            case S:
                s = sparse_[0];
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
            case Z:
                return 0;
            case I:
                return rows_;
            case D:
                return rows_;
            case S:
                return sparse_[0].nonZeros();
            default:
                OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << type_);
            }
        }


        double coeff(const int row, const int col) const
        {
            switch (type_) {
            case Z:
                return 0.0;
            case I:
                return (row == col) ? 1.0 : 0.0;
            case D:
                return (row == col) ? diag_[row] : 0.0;
            case S:
                return sparse_[0].coeff(row, col);
            default:
                OPM_THROW(std::logic_error, "Invalid AutoDiffMatrix type encountered: " << type_);
            }
        }

    private:
        enum MatrixType { Z, I, D, S };
        MatrixType type_;
        int rows_;
        int cols_;
        Diag diag_;

        /**
         * Eigen uses memory allocation within the default constructor, so that
         * Sparse a; actually calls malloc. To prevent this, we here use
         * Sparse a[1], and only construct the object when needed.
         */
        Sparse sparse_[1];

        template <class V>
        static inline
        Sparse
        spdiag(const V& d)
        {
            const int n = d.size();
            Sparse mat(n, n);
            mat.reserve(Eigen::ArrayXi::Ones(n, 1));
            for (Sparse::Index i = 0; i < n; ++i) {
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
