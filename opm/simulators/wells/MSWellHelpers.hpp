/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.

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


#ifndef OPM_MSWELLHELPERS_HEADER_INCLUDED
#define OPM_MSWELLHELPERS_HEADER_INCLUDED

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/istl/solvers.hh>
#if HAVE_UMFPACK
#include <dune/istl/umfpack.hh>
#endif // HAVE_UMFPACK
#include <cmath>

namespace Opm {

namespace mswellhelpers
{
    // obtain y = D^-1 * x with a direct solver
    template <typename MatrixType, typename VectorType>
    VectorType
    invDXDirect(const MatrixType& D, VectorType x)
    {
#if HAVE_UMFPACK
        VectorType y(x.size());
        y = 0.;

        Dune::UMFPack<MatrixType> linsolver(D, 0);

        // Object storing some statistics about the solving process
        Dune::InverseOperatorResult res;

        // Solve
        linsolver.apply(y, x, res);

        // Checking if there is any inf or nan in y
        // it will be the solution before we find a way to catch the singularity of the matrix
        for (size_t i_block = 0; i_block < y.size(); ++i_block) {
            for (size_t i_elem = 0; i_elem < y[i_block].size(); ++i_elem) {
                if (std::isinf(y[i_block][i_elem]) || std::isnan(y[i_block][i_elem]) ) {
                    OPM_THROW(Opm::NumericalIssue, "nan or inf value found in invDXDirect due to singular matrix");
                }
            }
        }

        return y;
#else
        // this is not thread safe
        OPM_THROW(std::runtime_error, "Cannot use invDXDirect() without UMFPACK. "
                  "Reconfigure opm-simulator with SuiteSparse/UMFPACK support and recompile.");
#endif // HAVE_UMFPACK
    }





    // obtain y = D^-1 * x with a BICSSTAB iterative solver
    template <typename MatrixType, typename VectorType>
    VectorType
    invDX(const MatrixType& D, VectorType x, Opm::DeferredLogger& deferred_logger)
    {
        // the function will change the value of x, so we should not use reference of x here.

        // TODO: store some of the following information to avoid to call it again and again for
        // efficiency improvement.
        // Bassically, only the solve / apply step is different.

        VectorType y(x.size());
        y = 0.;

        Dune::MatrixAdapter<MatrixType, VectorType, VectorType> linearOperator(D);

        // Sequential incomplete LU decomposition as the preconditioner
        Dune::SeqILU0<MatrixType, VectorType, VectorType> preconditioner(D, 1.0);
        // Dune::SeqILUn<MatrixType, VectorType, VectorType> preconditioner(D, 1, 0.92);
        // Dune::SeqGS<MatrixType, VectorType, VectorType> preconditioner(D, 1, 1);
        // Dune::SeqJac<MatrixType, VectorType, VectorType> preconditioner(D, 1, 1);

        // Preconditioned BICGSTAB solver
        Dune::BiCGSTABSolver<VectorType> linsolver(linearOperator,
                                                   preconditioner,
                                                   1.e-8, // desired residual reduction factor
                                                   250, // maximum number of iterations
                                                   0); // verbosity of the solver */

        // Object storing some statistics about the solving process
        Dune::InverseOperatorResult res;

        // Solve
        linsolver.apply(y, x, res);

        if ( !res.converged ) {
            OPM_DEFLOG_THROW(Opm::NumericalIssue, "the invDX does not get converged! ", deferred_logger);
        }

        return y;
    }




    template <typename ValueType>
    inline ValueType haalandFormular(const ValueType& re, const double diameter, const double roughness)
    {
        const ValueType value = -3.6 * Opm::log10(6.9 / re + std::pow(roughness / (3.7 * diameter), 10. / 9.) );

        // sqrt(1/f) should be non-positive
        assert(value >= 0.0);

        return 1. / (value * value);
    }




    template <typename ValueType>
    inline ValueType calculateFrictionFactor(const double area, const double diameter,
                                          const ValueType& w, const double roughness, const ValueType& mu)
    {

        ValueType f = 0.;
        // Reynolds number
        const ValueType re = Opm::abs( diameter * w / (area * mu));

        if ( re == 0.0 ) {
            // make sure it is because the mass rate is zero
            assert(w == 0.);
            return 0.0;
        }

        const ValueType re_value1 = 2000.;
        const ValueType re_value2 = 4000.;

        if (re < re_value1) {
            f = 16. / re;
        } else if (re > re_value2){
            f = haalandFormular(re, diameter, roughness);
        } else { // in between
            const ValueType f1 = 16. / re_value1;
            const ValueType f2 = haalandFormular(re_value2, diameter, roughness);

            f = (f2 - f1) / (re_value2 - re_value1) * (re - re_value1) + f1;
        }
        return f;
    }






    // calculating the friction pressure loss
    // l is the segment length
    // area is the segment cross area
    // diameter is the segment inner diameter
    // w is mass flow rate through the segment
    // density is density
    // roughness is the absolute roughness
    // mu is the average phase viscosity
    template <typename ValueType>
    ValueType frictionPressureLoss(const double l, const double diameter, const double area, const double roughness,
                                   const ValueType& density, const ValueType& w, const ValueType& mu)
    {
        const ValueType f = calculateFrictionFactor(area, diameter, w, roughness, mu);
        // \Note: a factor of 2 needs to be here based on the dimensional analysis
        return 2. * f * l * w * w / (area * area * diameter * density);
    }





    template <typename ValueType>
    ValueType velocityHead(const double area, const ValueType& mass_rate, const ValueType& density)
    {
        return (0.5 * mass_rate * mass_rate / (area * area * density));
    }


} // namespace mswellhelpers

}

#endif
