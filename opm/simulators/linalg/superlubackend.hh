// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::Linear::SuperLUBackend
 */
#ifndef EWOMS_SUPER_LU_BACKEND_HH
#define EWOMS_SUPER_LU_BACKEND_HH

#if HAVE_SUPERLU

#include <opm/models/linear/istlsparsematrixbackend.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/simulators/linalg/linalgproperties.hh>

#include <opm/material/common/Unused.hpp>

#include <dune/istl/superlu.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/version.hh>

namespace Opm::Properties::TTag {
struct SuperLULinearSolver {};
} // namespace Opm::Properties::TTag

namespace Opm {
namespace Linear {
template <class Scalar, class TypeTag, class Matrix, class Vector>
class SuperLUSolve_;

/*!
 * \ingroup Linear
 * \brief A linear solver backend for the SuperLU sparse matrix library.
 */
template <class TypeTag>
class SuperLUBackend
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using Matrix = typename SparseMatrixAdapter::block_type;

    static_assert(std::is_same<SparseMatrixAdapter, IstlSparseMatrixAdapter<MatrixBlock>::value,
                  "The SuperLU linear solver backend requires the IstlSparseMatrixAdapter");

public:
    SuperLUBackend(Simulator& simulator OPM_UNUSED)
    {}

    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverVerbosity,
                             "The verbosity level of the linear solver");
    }

    /*!
     * \brief Causes the solve() method to discared the structure of the linear system of
     *        equations the next time it is called.
     *
     * Since the SuperLU backend does not create any internal matrices, this is a no-op.
     */
    void eraseMatrix()
    { }

    void prepare(const SparseMatrixAdapter& M, const Vector& b)
    { }

    void setResidual(const Vector& b)
    { b_ = &b; }

    void getResidual(Vector& b) const
    { b = *b_; }

    void setMatrix(const SparseMatrixAdapter& M)
    { M_ = &M; }

    bool solve(Vector& x)
    { return SuperLUSolve_<Scalar, TypeTag, Matrix, Vector>::solve_(*M_, x, *b_); }

private:
    const Matrix* M_;
    Vector* b_;
};

template <class Scalar, class TypeTag, class Matrix, class Vector>
class SuperLUSolve_
{
public:
    static bool solve_(const Matrix& A, Vector& x, const Vector& b)
    {
        Vector bTmp(b);

        int verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);
        Dune::InverseOperatorResult result;
        Dune::SuperLU<Matrix> solver(A, verbosity > 0);
        solver.apply(x, bTmp, result);

        if (result.converged) {
            // make sure that the result only contains finite values.
            Scalar tmp = 0;
            for (unsigned i = 0; i < x.size(); ++i) {
                const auto& xi = x[i];
                for (unsigned j = 0; j < Vector::block_type::dimension; ++j)
                    tmp += xi[j];
            }
            result.converged = std::isfinite(tmp);
        }

        return result.converged;
    }
};

// the following is required to make the SuperLU adapter of dune-istl happy with
// quadruple precision math on Dune 2.4. this is because the most which SuperLU can
// handle is double precision (i.e., the linear systems of equations are always solved
// with at most double precision if chosing SuperLU as the linear solver...)
#if HAVE_QUAD
template <class TypeTag, class Matrix, class Vector>
class SuperLUSolve_<__float128, TypeTag, Matrix, Vector>
{
public:
    static bool solve_(const Matrix& A,
                       Vector& x,
                       const Vector& b)
    {
        static const int numEq = getPropValue<TypeTag, Properties::NumEq>();
        using DoubleEqVector = Dune::FieldVector<double, numEq>;
        using DoubleEqMatrix = Dune::FieldMatrix<double, numEq, numEq>;
        using DoubleVector = Dune::BlockVector<DoubleEqVector>;
        using DoubleMatrix = Dune::BCRSMatrix<DoubleEqMatrix>;

        // copy the inputs into the double precision data structures
        DoubleVector bDouble(b);
        DoubleVector xDouble(x);
        DoubleMatrix ADouble(A);

        bool res =
            SuperLUSolve_<double, TypeTag, Matrix, Vector>::solve_(ADouble,
                                                                   xDouble,
                                                                   bDouble);

        // copy the result back into the quadruple precision vector.
        x = xDouble;

        return res;
    }
};
#endif

} // namespace Linear
} // namespace Opm

namespace Opm::Properties {

template<class TypeTag>
struct LinearSolverVerbosity<TypeTag, TTag::SuperLULinearSolver> { static constexpr int value = 0; };
template<class TypeTag>
struct LinearSolverBackend<TypeTag, TTag::SuperLULinearSolver> { using type = Opm::Linear::SuperLUBackend<TypeTag>; };

} // namespace Opm::Properties

#endif // HAVE_SUPERLU

#endif
