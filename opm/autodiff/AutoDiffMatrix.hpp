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

#include <opm/autodiff/fastSparseProduct.hpp>
#include <vector>


namespace Opm
{

    class AutoDiffMatrix
    {
    public:
        AutoDiffMatrix()
            : type_(Z),
              rows_(0),
              cols_(0)
        {
        }



        AutoDiffMatrix(const int rows, const int cols)
            : type_(Z),
              rows_(rows),
              cols_(cols)
        {
        }



        enum CreationType { ZeroMatrix, IdentityMatrix };


        AutoDiffMatrix(const CreationType t, const int rows)
            : type_(t == ZeroMatrix ? Z : I),
              rows_(rows),
              cols_(rows)
        {
        }



        explicit AutoDiffMatrix(const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& d)
            : type_(D),
              rows_(d.rows()),
              cols_(d.cols()),
              d_(d.diagonal().array().data(), d.diagonal().array().data() + d.rows())
        {
        }



        explicit AutoDiffMatrix(const Eigen::SparseMatrix<double>& s)
            : type_(S),
              rows_(s.rows()),
              cols_(s.cols()),
              s_(s)
        {
        }



        AutoDiffMatrix(const AutoDiffMatrix& other) = default;
        AutoDiffMatrix& operator=(const AutoDiffMatrix& other) = default;



        AutoDiffMatrix(AutoDiffMatrix&& other)
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
            d_.swap(other.d_);
            s_.swap(other.s_);
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
                    return sumII(*this, rhs);
                case D:
                    return rhs + (*this);
                case S:
                    return rhs + (*this);
                }
            case D:
                switch (rhs.type_) {
                case Z:
                    return *this;
                case I:
                    return sumDI(*this, rhs);
                case D:
                    return sumDD(*this, rhs);
                case S:
                    return rhs + (*this);
                }
            case S:
                switch (rhs.type_) {
                case Z:
                    return *this;
                case I:
                    return sumSI(*this, rhs);
                case D:
                    return sumSD(*this, rhs);
                case S:
                    return sumSS(*this, rhs);
                }
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
                }
            case D:
                switch (rhs.type_) {
                case Z:
                    return AutoDiffMatrix(rows_, rhs.cols_);
                case I:
                    return *this;
                case D:
                    return prodDD(*this, rhs);
                case S:
                    return prodDS(*this, rhs);
                }
            case S:
                switch (rhs.type_) {
                case Z:
                    return AutoDiffMatrix(rows_, rhs.cols_);
                case I:
                    return *this;
                case D:
                    return prodSD(*this, rhs);
                case S:
                    return prodSS(*this, rhs);
                }
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
                    retval.d_.assign(rows_, rhs);
                    return retval;
                }
            case D:
                {
                    AutoDiffMatrix retval(*this);
                    for (double& elem : retval.d_) {
                        elem *= rhs;
                    }
                    return retval;
                }
            case S:
                {
                    AutoDiffMatrix retval(*this);
                    retval.s_ *= rhs;
                    return retval;
                }
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
                return Eigen::Map<const Eigen::VectorXd>(d_.data(), rows_) * rhs;
            case S:
                return s_ * rhs;
            }
        }






        static AutoDiffMatrix sumII(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == I);
            assert(rhs.type_ == I);
            AutoDiffMatrix retval;
            retval.type_ = D;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.d_.assign(lhs.rows_, 2.0);
            return retval;
        }

        static AutoDiffMatrix sumDI(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            static_cast<void>(rhs); // Silence release-mode warning.
            assert(lhs.type_ == D);
            assert(rhs.type_ == I);
            AutoDiffMatrix retval = lhs;
            for (int r = 0; r < lhs.rows_; ++r) {
                retval.d_[r] += 1.0;
            }
            return retval;
        }

        static AutoDiffMatrix sumDD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == D);
            assert(rhs.type_ == D);
            AutoDiffMatrix retval = lhs;
            for (int r = 0; r < lhs.rows_; ++r) {
                retval.d_[r] += rhs.d_[r];
            }
            return retval;
        }

        static AutoDiffMatrix sumSI(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == S);
            assert(rhs.type_ == I);
            AutoDiffMatrix retval;
            Eigen::SparseMatrix<double> ident = spdiag(Eigen::VectorXd::Ones(lhs.rows_));
            retval.type_ = S;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.s_ = lhs.s_ + ident;
            return retval;
        }

        static AutoDiffMatrix sumSD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == S);
            assert(rhs.type_ == D);
            AutoDiffMatrix retval;
            Eigen::SparseMatrix<double> diag = spdiag(rhs.d_);
            retval.type_ = S;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.s_ = lhs.s_ + diag;
            return retval;
        }

        static AutoDiffMatrix sumSS(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == S);
            assert(rhs.type_ == S);
            AutoDiffMatrix retval;
            retval.type_ = S;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.s_ = lhs.s_ + rhs.s_;
            return retval;
        }





        static AutoDiffMatrix prodDD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == D);
            assert(rhs.type_ == D);
            AutoDiffMatrix retval;
            retval.type_ = D;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            retval.d_.resize(lhs.rows_);
            for (int r = 0; r < lhs.rows_; ++r) {
                retval.d_[r] = lhs.d_[r] * rhs.d_[r];
            }
            return retval;
        }

        static AutoDiffMatrix prodDS(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == D);
            assert(rhs.type_ == S);
            AutoDiffMatrix retval;
            retval.type_ = S;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            fastDiagSparseProduct(lhs.d_, rhs.s_, retval.s_);
            return retval;
        }

        static AutoDiffMatrix prodSD(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == S);
            assert(rhs.type_ == D);
            AutoDiffMatrix retval;
            retval.type_ = S;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            fastSparseDiagProduct(lhs.s_, rhs.d_, retval.s_);
            return retval;
        }

        static AutoDiffMatrix prodSS(const AutoDiffMatrix& lhs, const AutoDiffMatrix& rhs)
        {
            assert(lhs.type_ == S);
            assert(rhs.type_ == S);
            AutoDiffMatrix retval;
            retval.type_ = S;
            retval.rows_ = lhs.rows_;
            retval.cols_ = rhs.cols_;
            fastSparseProduct(lhs.s_, rhs.s_, retval.s_);
            return retval;
        }



        void toSparse(Eigen::SparseMatrix<double>& s) const
        {
            switch (type_) {
            case Z:
                s = Eigen::SparseMatrix<double>(rows_, cols_);
                return;
            case I:
                s = spdiag(Eigen::VectorXd::Ones(rows_));
                return;
            case D:
                s = spdiag(d_);
                return;
            case S:
                s = s_;
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
                return s_.nonZeros();
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
                return (row == col) ? d_[row] : 0.0;
            case S:
                return s_.coeff(row, col);
            }
        }

    private:
        enum MatrixType { Z, I, D, S };
        MatrixType type_;
        int rows_;
        int cols_;
        std::vector<double> d_;
        Eigen::SparseMatrix<double> s_;

        template <class V>
        static inline
        Eigen::SparseMatrix<double>
        spdiag(const V& d)
        {
            typedef Eigen::SparseMatrix<double> M;
            const int n = d.size();
            M mat(n, n);
            mat.reserve(Eigen::ArrayXi::Ones(n, 1));
            for (M::Index i = 0; i < n; ++i) {
                mat.insert(i, i) = d[i];
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
