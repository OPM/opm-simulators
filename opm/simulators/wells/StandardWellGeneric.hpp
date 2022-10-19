/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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


#ifndef OPM_STANDARDWELL_GENERIC_HEADER_INCLUDED
#define OPM_STANDARDWELL_GENERIC_HEADER_INCLUDED

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <opm/simulators/wells/WellHelpers.hpp>

#include <optional>
#include <vector>

namespace Opm
{

class ConvergenceReport;
class DeferredLogger;
class ParallelWellInfo;
class Schedule;
class SummaryState;
class WellInterfaceGeneric;
class WellState;

template<class Scalar>
class StandardWellGeneric
{
protected:
    // sparsity pattern for the matrices
    //[A C^T    [x       =  [ res
    // B  D ]   x_well]      res_well]

    // the vector type for the res_well and x_well
    using VectorBlockWellType = Dune::DynamicVector<Scalar>;
    using BVectorWell = Dune::BlockVector<VectorBlockWellType>;

    // the matrix type for the diagonal matrix D
    using DiagMatrixBlockWellType = Dune::DynamicMatrix<Scalar>;
    using DiagMatWell = Dune::BCRSMatrix<DiagMatrixBlockWellType>;

    // the matrix type for the non-diagonal matrix B and C^T
    using OffDiagMatrixBlockWellType = Dune::DynamicMatrix<Scalar>;
    using OffDiagMatWell = Dune::BCRSMatrix<OffDiagMatrixBlockWellType>;

public:
    /// get the number of blocks of the C and B matrices, used to allocate memory in a WellContributions object
    unsigned int getNumBlocks() const;

protected:
    StandardWellGeneric(const WellInterfaceGeneric& baseif);

    // calculate a relaxation factor to avoid overshoot of total rates
    static double relaxationFactorRate(const std::vector<double>& primary_variables,
                                       const BVectorWell& dwells);

    // relaxation factor considering only one fraction value
    static double relaxationFactorFraction(const double old_value,
                                           const double dx);

    void computeConnectionPressureDelta();

    // Base interface reference
    const WellInterfaceGeneric& baseif_;

    // residuals of the well equations
    BVectorWell resWell_;

    // densities of the fluid in each perforation
    std::vector<double> perf_densities_;
    // pressure drop between different perforations
    std::vector<double> perf_pressure_diffs_;

    // two off-diagonal matrices
    OffDiagMatWell duneB_;
    OffDiagMatWell duneC_;
    // diagonal matrix for the well
    DiagMatWell invDuneD_;
    DiagMatWell duneD_;

    // Wrapper for the parallel application of B for distributed wells
    wellhelpers::ParallelStandardWellB<Scalar> parallelB_;

    // several vector used in the matrix calculation
    mutable BVectorWell Bx_;
    mutable BVectorWell invDrw_;

    double getRho() const { return perf_densities_[0]; }
};

}

#endif // OPM_STANDARDWELL_GENERIC_HEADER_INCLUDED
