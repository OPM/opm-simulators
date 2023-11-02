/*
  Copyright 2022-2023 SINTEF AS
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
#define BOOST_TEST_MODULE TestSeqDILU

#include <config.h>
#include <opm/simulators/linalg/DILU.hpp>

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>

using NumericTypes = boost::mpl::list<double, float>;


BOOST_AUTO_TEST_CASE_TEMPLATE(SeqDILUDiagIsCorrect2x2NoZeros, T, NumericTypes)
{
    /*
        Tests that the dilu decomposition mathces the expected result
        for a 2x2 matrix with no zero blocks, with block size 2x2.
 
                 A         
        | | 3  1| | 1  0| |
        | | 2  1| | 0  1| |
        |                 |
        | | 2  0| |-1  0| |
        | | 0  2| | 0 -1| |
    */

    const int N = 2;
    constexpr int bz = 2;
    const int nonZeroes = 4;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;

    Matrix A(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = A.createbegin(); row != A.createend(); ++row) {
        row.insert(0);
        row.insert(1);
    }

    A[0][0][0][0]=3.0;
    A[0][0][0][1]=1.0;
    A[0][0][1][0]=2.0;
    A[0][0][1][1]=1.0;

    A[0][1][0][0]=1.0;
    A[0][1][1][1]=1.0;

    A[1][0][0][0]=2.0;
    A[1][0][1][1]=2.0;


    A[1][1][0][0]=-1.0;
    A[1][1][1][1]=-1.0;


    auto D_00 = A[0][0];
    auto D_00_inv = D_00;
    D_00_inv.invert();
    // D_11 = A_11 - L_10 D_00_inv U_01
    auto D_11 = A[1][1] - A[1][0]*D_00_inv*A[0][1];

    Dune::SeqDilu<Matrix, Vector, Vector> seqdilu(A);

    auto Dinv = seqdilu.getDiagonal();

    // diagonal stores inverse
    auto D_00_dilu = Dinv[0];
    D_00_dilu.invert();
    auto D_11_dilu = Dinv[1];
    D_11_dilu.invert();


    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            BOOST_CHECK_CLOSE(D_00_dilu[i][j], D_00[i][j], 1e-7);
            BOOST_CHECK_CLOSE(D_11_dilu[i][j], D_11[i][j], 1e-7);
        }
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(SeqDILUDiagIsCorrect2x2, T, NumericTypes)
{
    /*
        Tests that the dilu decomposition mathces the expected result
        for a 2x2 matrix, with block size 2x2.
        
                 A         
        | | 3  1| | 1  0| |
        | | 2  1| | 0  1| |
        |                 |
        | | 0  0| |-1  0| |
        | | 0  0| | 0 -1| |
    */

    const int N = 2;
    constexpr int bz = 2;
    const int nonZeroes = 3;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;

    Matrix A(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = A.createbegin(); row != A.createend(); ++row) {
        row.insert(row.index());
        if (row.index() == 0) {
            row.insert(row.index() + 1);
        }
    }

    A[0][0][0][0]=3.0;
    A[0][0][0][1]=1.0;
    A[0][0][1][0]=2.0;
    A[0][0][1][1]=1.0;

    A[0][1][0][0]=1.0;
    A[0][1][1][1]=1.0;

    A[1][1][0][0]=-1.0;
    A[1][1][1][1]=-1.0;


    auto D_00 = A[0][0];
    auto D_00_inv = D_00;
    D_00_inv.invert();
    // D_11 = A_11 - L_10 D_00_inv U_01 = A_11
    auto D_11 = A[1][1];

    Dune::SeqDilu<Matrix, Vector, Vector> seqdilu(A);


    auto Dinv = seqdilu.getDiagonal();

    // diagonal stores inverse
    auto D_00_dilu = Dinv[0];
    D_00_dilu.invert();
    auto D_11_dilu = Dinv[1];
    D_11_dilu.invert();

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            BOOST_CHECK_CLOSE(D_00_dilu[i][j], D_00[i][j], 1e-7);
            BOOST_CHECK_CLOSE(D_11_dilu[i][j], D_11[i][j], 1e-7);
        }
    }
}




BOOST_AUTO_TEST_CASE_TEMPLATE(SeqDILUApplyIsCorrectNoZeros, T, NumericTypes)
{
    /*
        Tests that applying the dilu preconditioner mathces the expected result
        for a 2x2 matrix with no zero blocks, with block size 2x2.

                 A               x     =     b
        | | 3  1| | 1  0| |   | |1| |     | |2| |
        | | 2  1| | 0  1| |   | |2| |     | |1| |
        |                 | x |     |  =  |     |
        | | 2  0| |-1  0| |   | |1| |     | |3| |
        | | 0  2| | 0 -1| |   | |1| |     | |4| |
    */

    const int N = 2;
    constexpr int bz = 2;
    const int nonZeroes = 4;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;

    Matrix A(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = A.createbegin(); row != A.createend(); ++row) {
        row.insert(0);
        row.insert(1);
    }

    A[0][0][0][0]=3.0;
    A[0][0][0][1]=1.0;
    A[0][0][1][0]=2.0;
    A[0][0][1][1]=1.0;

    A[0][1][0][0]=1.0;
    A[0][1][1][1]=1.0;

    A[1][0][0][0]=2.0;
    A[1][0][1][1]=2.0;


    A[1][1][0][0]=-1.0;
    A[1][1][1][1]=-1.0;


    Vector x(2);
    x[0][0] = 1.0;
    x[0][1] = 2.0;
    x[1][0] = 1.0;
    x[1][1] = 1.0;

    Vector b(2);
    b[0][0] = 2.0;
    b[0][1] = 1.0;
    b[1][0] = 3.0;
    b[1][1] = 4.0;


    auto D_00 = A[0][0];
    auto D_00_inv = D_00;
    D_00_inv.invert();
    // D_11= A_11 - L_10 D_00_inv U_01
    auto D_11 = A[1][1] - A[1][0]*D_00_inv*A[0][1];
    auto D_11_inv = D_11;
    D_11_inv.invert(); 

    // define: z = M^-1(b - Ax)
    // where: M = (D + L_A) A^-1 (D + U_A)
    // lower triangular solve: (E + L) y = b - Ax
    // y_0 = D_00_inv*[b_0 - (A_00*x_0 + A_01*x_1)]
    Vector y(2);

    auto rhs = b[0];
    A[0][0].mmv(x[0], rhs);
    A[0][1].mmv(x[1], rhs);
    D_00_inv.mv(rhs, y[0]);

    // y_1 = D_11_inv*(b_1 - (A_10*x_0 + A_11*x_1) - A_10*y_0)
    rhs = b[1];
    A[1][0].mmv(x[0], rhs);
    A[1][1].mmv(x[1], rhs);
    A[1][0].mmv(y[0], rhs);
    D_11_inv.mv(rhs, y[1]);


    // upper triangular solve: (E + U) z = Ey
    // z_1 = y_1
    Vector z(2);
    z[1] = y[1];

    // z_0 = y_0 - D_00_inv*A_01*z_1
    z[0] = y[0];
    auto temp = D_00_inv*A[0][1];
    temp.mmv(z[1], z[0]); 

    // x_k+1 = x_k + z
    Vector new_x = x;
    new_x += z;
    
    Dune::SeqDilu<Matrix, Vector, Vector> seqdilu(A);
    seqdilu.apply(x, b);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            BOOST_CHECK_CLOSE(x[i][j], new_x[i][j], 1e-7);
        }
    }
}



BOOST_AUTO_TEST_CASE_TEMPLATE(SeqDILUApplyIsCorrect1, T, NumericTypes)
{
    /*
        Tests that applying the dilu preconditioner mathces the expected result
        for a 2x2 matrix, with block size 2x2.

                 A               x     =     b
        | | 3  1| | 1  0| |   | |1| |     | |2| |
        | | 2  1| | 0  1| |   | |2| |     | |1| |
        |                 | x |     |  =  |     |
        | | 0  0| |-1  0| |   | |1| |     | |3| |
        | | 0  0| | 0 -1| |   | |1| |     | |4| |
    */


    const int N = 2;
    constexpr int bz = 2;
    const int nonZeroes = 3;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;

    Matrix A(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = A.createbegin(); row != A.createend(); ++row) {
        row.insert(row.index());
        if (row.index() == 0) {
            row.insert(row.index() + 1);
        }
    }

    A[0][0][0][0]=3.0;
    A[0][0][0][1]=1.0;
    A[0][0][1][0]=2.0;
    A[0][0][1][1]=1.0;

    A[0][1][0][0]=1.0;
    A[0][1][1][1]=1.0;

    A[1][1][0][0]=-1.0;
    A[1][1][1][1]=-1.0;


    Vector x(2);
    x[0][0] = 1.0;
    x[0][1] = 2.0;
    x[1][0] = 1.0;
    x[1][1] = 1.0;

    Vector b(2);
    b[0][0] = 2.0;
    b[0][1] = 1.0;
    b[1][0] = 3.0;
    b[1][1] = 4.0;


    auto D_00 = A[0][0];
    auto D_00_inv = D_00;
    D_00_inv.invert();
    // D_11 = A_11 - L_10 D_0_inv U_01 = A_11
    auto D_11 = A[1][1];
    auto D_11_inv = D_11;
    D_11_inv.invert(); 

    // define: z = M^-1(b - Ax)
    // where: M = (D + L_A) A^-1 (D + U_A)
    // lower triangular solve: (E + L) y = b - Ax
    // y_0 = D_00_inv*[b_0 - (A_00*x_0 + A_01*x_1)]
    Vector y(2);

    auto rhs = b[0];
    A[0][0].mmv(x[0], rhs);
    A[0][1].mmv(x[1], rhs);
    D_00_inv.mv(rhs, y[0]);

    // y_1 = D_11_inv*(b_1 - (A_10*x_0 + A_11*x_1) - A_10*y_0) = D_11_inv*(b_1 - A_11*x_1)
    rhs = b[1];
    A[1][1].mmv(x[1], rhs);
    D_11_inv.mv(rhs, y[1]);


    // upper triangular solve: (E + U) z = Ey
    // z_1 = y_1
    Vector z(2);
    z[1] = y[1];

    // z_0 = y_0 - D_00_inv*A_01*z_1
    z[0] = y[0];
    auto temp = D_00_inv*A[0][1];
    temp.mmv(z[1], z[0]); 

    // x_k+1 = x_k + z
    Vector new_x = x;
    new_x += z;
    
    Dune::SeqDilu<Matrix, Vector, Vector> seqdilu(A);
    seqdilu.apply(x, b);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            BOOST_CHECK_CLOSE(x[i][j], new_x[i][j], 1e-7);
        }
    }
}



BOOST_AUTO_TEST_CASE_TEMPLATE(SeqDILUApplyIsCorrect2, T, NumericTypes)
{
    /*
        Tests that applying the dilu preconditioner mathces the expected result
        for a 2x2 matrix, with block size 2x2.

                 A               x     =     b
        | | 3  1| | 0  0| |   | |1| |     | |2| |
        | | 2  1| | 0  0| |   | |2| |     | |1| |
        |                 | x |     |  =  |     |
        | | 2  0| |-1  0| |   | |1| |     | |3| |
        | | 0  2| | 0 -1| |   | |1| |     | |4| |
    */
    const int N = 2;
    constexpr int bz = 2;
    const int nonZeroes = 3;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;

    Matrix A(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = A.createbegin(); row != A.createend(); ++row) {
        row.insert(row.index());
        if (row.index() == 1) {
            row.insert(row.index() - 1);
        }
    }

    A[0][0][0][0]=3.0;
    A[0][0][0][1]=1.0;
    A[0][0][1][0]=2.0;
    A[0][0][1][1]=1.0;

    A[1][1][0][0]=2.0;
    A[1][1][1][1]=2.0;

    A[1][1][0][0]=-1.0;
    A[1][1][1][1]=-1.0;

    Vector x(2);
    x[0][0] = 1.0;
    x[0][1] = 2.0;
    x[1][0] = 1.0;
    x[1][1] = 1.0;

    Vector b(2);
    b[0][0] = 2.0;
    b[0][1] = 1.0;
    b[1][0] = 3.0;
    b[1][1] = 4.0;


    auto D_00 = A[0][0];
    auto D_00_inv = D_00;
    D_00_inv.invert();
    // D_11 = A_11 - L_10 D_0_inv U_01 = A_11
    auto D_11 = A[1][1];
    auto D_11_inv = D_11;
    D_11_inv.invert(); 

    // define: z = M^-1(b - Ax)
    // where: M = (D + L_A) A^-1 (D + U_A)
    // lower triangular solve: (E + L) y = b - Ax
    // y_0 = D_00_inv*[b_0 - (A_00*x_0 + A_01*x_1)] = D_00_inv*[b_0 - A_00*x_0]
    Vector y(2);
    auto rhs = b[0];

    A[0][0].mmv(x[0], rhs);
    D_00_inv.mv(rhs, y[0]);

    // y_1 = D_11_inv*(b_1 - (A_10*x_0 + A_11*x_1) - A_10*y_0)
    rhs = b[1];
    A[1][1].mmv(x[1], rhs);
    D_11_inv.mv(rhs, y[1]);


    // upper triangular solve: (E + U) z = Ey
    // z_1 = y_1
    Vector z(2);
    z[1] = y[1];

    // z_0 = y_0 - D_00_inv*A_01*z_1 = y_0
    z[0] = y[0];

    // x_k+1 = x_k + z
    Vector new_x = x;
    new_x += z;
    
    Dune::SeqDilu<Matrix, Vector, Vector> seqdilu(A);
    seqdilu.apply(x, b);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            BOOST_CHECK_CLOSE(x[i][j], new_x[i][j], 1e-7);
        }
    }
}


BOOST_AUTO_TEST_CASE_TEMPLATE(SeqDILUDiagIsCorrect3x3, T, NumericTypes)
{
    /*
        Tests that the dilu decomposition mathces the expected result
        for a 3x3 matrix, with block size 3x3.

                          A                 
        | | 3  1  2| | 0  0  0| | 0  0  0| |
        | | 2  3  1| | 0  0  0| | 0  0  0| |
        | | 2  1  0| | 0  0  0| | 0  0  0| |
        |                                  |
        | | 0  0  0| | 1  0  1| | 1  0  2| |
        | | 0  0  0| | 4  1  0| | 0  1  1| |  
        | | 0  0  0| | 3  1  3| | 0  1  3| |
        |                                  |
        | | 0  0  0| | 1  0  2| | 1  3  2| |
        | | 0  0  0| | 0  1  4| | 2  1  3| |
        | | 0  0  0| | 5  1  1| | 3  1  2| |
    */


    const int N = 3;
    constexpr int bz = 3;
    const int nonZeroes = 5;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;

    Matrix A(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = A.createbegin(); row != A.createend(); ++row) {
        if (row.index() == 0) {
            row.insert(row.index());
        }
        else if (row.index() == 1) {
            row.insert(row.index());
            row.insert(row.index() + 1);
        }
        else if (row.index() == 2) {
            row.insert(row.index() - 1);
            row.insert(row.index());
        }
    }

    A[0][0][0][0]=3.0;    A[1][1][0][0]=1.0;    A[1][2][0][0]=1.0;
    A[0][0][0][1]=1.0;    A[1][1][0][1]=0.0;    A[1][2][0][1]=0.0;
    A[0][0][0][2]=2.0;    A[1][1][0][2]=1.0;    A[1][2][0][2]=2.0;
    A[0][0][1][0]=2.0;    A[1][1][1][0]=4.0;    A[1][2][1][0]=0.0;
    A[0][0][1][1]=3.0;    A[1][1][1][1]=1.0;    A[1][2][1][1]=1.0;
    A[0][0][1][2]=1.0;    A[1][1][1][2]=0.0;    A[1][2][1][2]=1.0;
    A[0][0][2][0]=2.0;    A[1][1][2][0]=3.0;    A[1][2][2][0]=0.0;
    A[0][0][2][1]=1.0;    A[1][1][2][1]=1.0;    A[1][2][2][1]=1.0;
    A[0][0][2][2]=0.0;    A[1][1][2][2]=3.0;    A[1][2][2][2]=3.0;

    A[2][1][0][0]=1.0;    A[2][2][0][0]=1.0;
    A[2][1][0][1]=0.0;    A[2][2][0][1]=3.0;
    A[2][1][0][2]=2.0;    A[2][2][0][2]=2.0;
    A[2][1][1][0]=0.0;    A[2][2][1][0]=2.0;
    A[2][1][1][1]=1.0;    A[2][2][1][1]=1.0;
    A[2][1][1][2]=4.0;    A[2][2][1][2]=3.0;
    A[2][1][2][0]=5.0;    A[2][2][2][0]=3.0;
    A[2][1][2][1]=1.0;    A[2][2][2][1]=1.0;
    A[2][1][2][2]=1.0;    A[2][2][2][2]=2.0;


    auto D_00 = A[0][0];
    auto D_00_inv = D_00;
    D_00_inv.invert();
    // D_11 = A_11 - L_10 D_00_inv U_01 = A_11
    auto D_11 = A[1][1];
    auto D_11_inv = D_11;
    D_11_inv.invert(); 
    // D_22 = A_22 - A_20 D_00_inv A_02 - A_21 D_11_inv A_12 =  A_22 - A_21 D_11_inv A_12
    auto D_22 = A[2][2] - A[2][1]*D_11_inv*A[1][2];
    auto D_22_inv = D_22;
    D_22_inv.invert(); 

    Dune::SeqDilu<Matrix, Vector, Vector> seqdilu(A);
    auto Dinv = seqdilu.getDiagonal();

    // diagonal stores inverse
    auto D_00_dilu = Dinv[0];
    D_00_dilu.invert();
    auto D_11_dilu = Dinv[1];
    D_11_dilu.invert();
    auto D_22_dilu = Dinv[2];
    D_22_dilu.invert();


    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            BOOST_CHECK_CLOSE(D_00_dilu[i][j], D_00[i][j], 1e-7);
            BOOST_CHECK_CLOSE(D_11_dilu[i][j], D_11[i][j], 1e-7);
            BOOST_CHECK_CLOSE(D_22_dilu[i][j], D_22[i][j], 1e-7);
        }
    }
}



BOOST_AUTO_TEST_CASE_TEMPLATE(SeqDILUApplyIsCorrect3, T, NumericTypes)
{
    /*
        Tests that applying the dilu preconditioner mathces the expected result
        for a 3x3 matrix, with block size 3x3.

                          A                       x    =     b
        | | 3  1  2| | 0  0  0| | 0  0  0| |   | |1| |    | |2| |
        | | 2  3  1| | 0  0  0| | 0  0  0| |   | |2| |    | |1| |
        | | 2  1  0| | 0  0  0| | 0  0  0| |   | |3| |    | |2| |
        |                                  |   |     |    |     |
        | | 0  0  0| | 1  0  1| | 1  0  2| |   | |1| |    | |2| |
        | | 0  0  0| | 4  1  0| | 0  1  1| | x | |3| |  = | |3| |  
        | | 0  0  0| | 3  1  3| | 0  1  3| |   | |2| |    | |2| |
        |                                  |   |     |    |     |
        | | 0  0  0| | 1  0  2| | 1  3  2| |   | |1| |    | |0| |
        | | 0  0  0| | 0  1  4| | 2  1  3| |   | |0| |    | |2| |
        | | 0  0  0| | 5  1  1| | 3  1  2| |   | |2| |    | |1| |
    
    */

    const int N = 3;
    constexpr int bz = 3;
    const int nonZeroes = 5;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;

    Matrix A(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = A.createbegin(); row != A.createend(); ++row) {
        if (row.index() == 0) {
            row.insert(row.index());
        }
        else if (row.index() == 1) {
            row.insert(row.index());
            row.insert(row.index() + 1);
        }
        else if (row.index() == 2) {
            row.insert(row.index() - 1);
            row.insert(row.index());
        }
    }

    A[0][0][0][0]=3.0;    A[1][1][0][0]=1.0;    A[1][2][0][0]=1.0;
    A[0][0][0][1]=1.0;    A[1][1][0][1]=0.0;    A[1][2][0][1]=0.0;
    A[0][0][0][2]=2.0;    A[1][1][0][2]=1.0;    A[1][2][0][2]=2.0;
    A[0][0][1][0]=2.0;    A[1][1][1][0]=4.0;    A[1][2][1][0]=0.0;
    A[0][0][1][1]=3.0;    A[1][1][1][1]=1.0;    A[1][2][1][1]=1.0;
    A[0][0][1][2]=1.0;    A[1][1][1][2]=0.0;    A[1][2][1][2]=1.0;
    A[0][0][2][0]=2.0;    A[1][1][2][0]=3.0;    A[1][2][2][0]=0.0;
    A[0][0][2][1]=1.0;    A[1][1][2][1]=1.0;    A[1][2][2][1]=1.0;
    A[0][0][2][2]=0.0;    A[1][1][2][2]=3.0;    A[1][2][2][2]=3.0;

    A[2][1][0][0]=1.0;    A[2][2][0][0]=1.0;
    A[2][1][0][1]=0.0;    A[2][2][0][1]=3.0;
    A[2][1][0][2]=2.0;    A[2][2][0][2]=2.0;
    A[2][1][1][0]=0.0;    A[2][2][1][0]=2.0;
    A[2][1][1][1]=1.0;    A[2][2][1][1]=1.0;
    A[2][1][1][2]=4.0;    A[2][2][1][2]=3.0;
    A[2][1][2][0]=5.0;    A[2][2][2][0]=3.0;
    A[2][1][2][1]=1.0;    A[2][2][2][1]=1.0;
    A[2][1][2][2]=1.0;    A[2][2][2][2]=2.0;

    Vector x(3);
    x[0][0] = 1.0;    x[1][0] = 1.0;    x[2][0] = 1.0;
    x[0][1] = 2.0;    x[1][1] = 3.0;    x[2][1] = 0.0;
    x[0][2] = 3.0;    x[1][2] = 2.0;    x[2][2] = 2.0;

    Vector b(3);
    b[0][0] = 2.0;    b[1][0] = 2.0;    b[2][0] = 0.0;
    b[0][1] = 1.0;    b[1][1] = 3.0;    b[2][1] = 2.0;
    b[0][2] = 2.0;    b[1][2] = 2.0;    b[2][2] = 1.0;


    // D_00 = A_00
    auto D_00 = A[0][0];
    auto D_00_inv = D_00;
    D_00_inv.invert();
    // D_11 = A_11 - L_10 D_00_inv U_01
    //      = A_11
    auto D_11 = A[1][1];
    auto D_11_inv = D_11;
    D_11_inv.invert(); 
    // D_22 = A_22 - A_20 D_00_inv A_02 - A_21 D_11_inv A_12
    //      = A_22 - A_21 D_11_inv A_12
    auto D_22 = A[2][2] - A[2][1]*D_11_inv*A[1][2];
    auto D_22_inv = D_22;
    D_22_inv.invert(); 

    // define: z = M^-1(b - Ax)
    // where: M = (D + L_A) A^-1 (D + U_A)
    // lower triangular solve: (E + L) y = b - Ax
    
    Vector y(3);
    // y_0 = D_00_inv*[b_0 - (A_00*x_0 + A_01*x_1)]
    //     = D_00_inv*[b_0 - A_00*x_0]
    auto rhs = b[0];
    A[0][0].mmv(x[0], rhs);
    D_00_inv.mv(rhs, y[0]);

    // y_1 = D_11_inv*(b_1 - (A_10*x_0 + A_11*x_1 + A_12*x_2) - A_10*y_0)
    //     = D_11_inv*(b_1 - A_11*x_1)
    rhs = b[1];
    A[1][1].mmv(x[1], rhs);
    A[1][2].mmv(x[2], rhs);
    D_11_inv.mv(rhs, y[1]);

    // y_2 = D_22_inv*(b_2 - (A_20*x_0 + A_21*x_1 + A_22*x_2) - (A_20*y_0 + A_21*y_1))
    //     = D_22_inv*(b_2 - (A_21*x_1 + A_22*x_2) - (A_21*y_1))
    rhs = b[2];
    A[2][1].mmv(x[1], rhs);
    A[2][2].mmv(x[2], rhs);
    A[2][1].mmv(y[1], rhs);
    D_22_inv.mv(rhs, y[2]);


    // upper triangular solve: (E + U) z = Ey
    Vector z(3);
    // z_2 = y_2
    z[2] = y[2];

    // z_1 = y_1 - D_11_inv*A_12*z_2
    z[1] = y[1];
    auto temp = D_11_inv*A[1][2];
    temp.mmv(z[2], z[1]); 

    // z_0 = y_0 - D_00_inv(A_01*z_1 + A_02*z_2)
    // z_0 = y_0
    z[0] = y[0];

    // x_k+1 = x_k + z
    Vector new_x = x;
    new_x += z;
    
    Dune::SeqDilu<Matrix, Vector, Vector> seqdilu(A);
    seqdilu.apply(x, b);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            BOOST_CHECK_CLOSE(x[i][j], new_x[i][j], 1e-7);
        }
    }
}



BOOST_AUTO_TEST_CASE_TEMPLATE(SeqDILUApplyIsEqualToDuneSeqILUApply, T, NumericTypes)
{
    /*
        Tests that applying the DILU preconditioner is equivalent to applying a ILU preconditioner
        for a block diagonal matrix.

                 A               x     =     b
        | | 3  1| | 0  0| |   | |1| |     | |2| |
        | | 2  1| | 0  0| |   | |2| |     | |1| |
        |                 | x |     |  =  |     |
        | | 0  0| |-1  0| |   | |1| |     | |3| |
        | | 0  0| | 0 -1| |   | |1| |     | |4| |
    */

    const int N = 2;
    constexpr int bz = 2;
    const int nonZeroes = 2;
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, bz, bz>>;
    using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;

    Matrix A(N, N, nonZeroes, Matrix::row_wise);
    for (auto row = A.createbegin(); row != A.createend(); ++row) {
        row.insert(row.index());
    }

    A[0][0][0][0]=3.0;
    A[0][0][0][1]=1.0;
    A[0][0][1][0]=2.0;
    A[0][0][1][1]=1.0;

    A[1][1][0][0]=-1.0;
    A[1][1][1][1]=-1.0;


    Dune::SeqDilu<Matrix, Vector, Vector> seqdilu(A);
    Dune::SeqILU<Matrix, Vector, Vector> seqilu(A, 1.0);

    Vector dilu_x(2);
    dilu_x[0][0] = 1.0;
    dilu_x[0][1] = 2.0;
    dilu_x[1][0] = 1.0;
    dilu_x[1][1] = 1.0;

    Vector dilu_b(2);
    dilu_b[0][0] = 2.0;
    dilu_b[0][1] = 1.0;
    dilu_b[1][0] = 3.0;
    dilu_b[1][1] = 4.0;

    Vector ilu_x(2);
    ilu_x[0][0] = 1.0;
    ilu_x[0][1] = 2.0;
    ilu_x[1][0] = 1.0;
    ilu_x[1][1] = 1.0;

    Vector ilu_b(2);
    ilu_b[0][0] = 2.0;
    ilu_b[0][1] = 1.0;
    ilu_b[1][0] = 3.0;
    ilu_b[1][1] = 4.0;

    seqdilu.apply(dilu_x, dilu_b);
    seqilu.apply(ilu_x, ilu_b);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            BOOST_CHECK_CLOSE(dilu_x[i][j], ilu_x[i][j], 1e-7);
        }
    }
}
