/*
  Copyright 2024 Equinor ASA

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

#ifndef OPM_CPRCREATION_HPP
#define OPM_CPRCREATION_HPP


#include <dune/istl/paamg/matrixhierarchy.hh>
#include <dune/istl/umfpack.hh>

#include <opm/simulators/linalg/gpubridge/Matrix.hpp>
#include <opm/simulators/linalg/gpubridge/Preconditioner.hpp>

#include <type_traits>

namespace Opm::Accelerator {

template<class Scalar> class BlockedMatrix;

/// This class implements a Constrained Pressure Residual (CPR) preconditioner
template <class Scalar, unsigned int block_size>
class CprCreation
{
    int cprN;
    int cprNb;
    int cprnnz;
    int cprnnzb;

public:     
    CprCreation();

protected:
    
    int num_levels;
    std::vector<Scalar> weights, coarse_vals, coarse_x, coarse_y;
    std::vector<Matrix<Scalar>> Amatrices, Rmatrices; // scalar matrices that represent the AMG hierarchy
    std::vector<std::vector<int> > PcolIndices; // prolongation does not need a full matrix, only store colIndices
    std::vector<std::vector<Scalar> > invDiags; // inverse of diagonal of Amatrices
    
    BlockedMatrix<Scalar> *mat = nullptr;    // input matrix, blocked

    using DuneMat = Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1> >;
    using DuneVec = Dune::BlockVector<Dune::FieldVector<Scalar, 1> >;
    using MatrixOperator = Dune::MatrixAdapter<DuneMat, DuneVec, DuneVec>;
    using DuneAmg = Dune::Amg::MatrixHierarchy<MatrixOperator, Dune::Amg::SequentialInformation>;
    std::unique_ptr<DuneAmg> dune_amg;
    std::unique_ptr<DuneMat> dune_coarse;       // extracted pressure matrix, finest level in AMG hierarchy
    std::shared_ptr<MatrixOperator> dune_op;    // operator, input to Dune AMG
    std::vector<int> level_sizes;               // size of each level in the AMG hierarchy
    std::vector<std::vector<int> > diagIndices; // index of diagonal value for each level
    std::conditional_t<std::is_same_v<Scalar,double>,
                       Dune::UMFPack<DuneMat>, int> umfpack; // dune/istl/umfpack object used to solve the coarsest level of AMG
    bool always_recalculate_aggregates = false; // OPM always reuses the aggregates by default
    bool recalculate_aggregates = true;         // only rerecalculate if true
    const int pressure_idx = 1;                 // hardcoded to mimic OPM
    unsigned num_pre_smooth_steps;              // number of Jacobi smooth steps before restriction
    unsigned num_post_smooth_steps;             // number of Jacobi smooth steps after prolongation

    // Analyze the AMG hierarchy build by Dune
    void analyzeHierarchy();

    // Analyze the aggregateMaps from the AMG hierarchy
    // These can be reused, so only use when recalculate_aggregates is true
    void analyzeAggregateMaps();

    void create_preconditioner_amg(BlockedMatrix<Scalar> *mat);
};

} // namespace Opm

#endif

