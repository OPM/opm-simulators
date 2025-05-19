/*
  Copyright 2025 Equinor ASA.

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

#ifndef OPM_GPU_PRESSURE_TRANSFER_POLICY_HEADER_INCLUDED
#define OPM_GPU_PRESSURE_TRANSFER_POLICY_HEADER_INCLUDED


#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/WellOperators.hpp>
#include <opm/simulators/linalg/gpuistl/GpuSparseMatrix.hpp>
#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl/detail/cpr_amg_operations.hpp>

#include <cstddef>

namespace Opm { namespace Details {
    template<class Scalar> using GpuPressureMatrixType = gpuistl::GpuSparseMatrix<Scalar>;
    template<class Scalar> using GpuPressureVectorType = gpuistl::GpuVector<Scalar>;
    template<class Scalar> using GpuSeqCoarseOperatorType = Dune::MatrixAdapter<GpuPressureMatrixType<Scalar>,
                                                                                GpuPressureVectorType<Scalar>,
                                                                                GpuPressureVectorType<Scalar>>;

    template<class Scalar, class Comm>
    using GpuCoarseOperatorType = GpuSeqCoarseOperatorType<Scalar>;
} // namespace Details
} // namespace Opm

namespace Opm::gpuistl
{
// Note: This is not a subclass of LevelTransferPolicyCpr anymore, but a wrapper that implements the same interface.
// This is because the constructor of LevelTransferPolicyCpr requires default constructors which we do not want with gpu types.
template <class FineOperator, class Communication, class Scalar, bool transpose = false>
class GpuPressureTransferPolicy
{
public:
    using CoarseOperator = typename Details::GpuCoarseOperatorType<Scalar,Communication>;
    using ParallelInformation = Communication;
    using GpuMatrixType = GpuSparseMatrix<Scalar>;

    using FineRangeType = typename FineOperator::range_type;
    using FineDomainType = typename FineOperator::domain_type;
    using CoarseRangeType = typename CoarseOperator::range_type;
    using CoarseDomainType = typename CoarseOperator::domain_type;

public:
    GpuPressureTransferPolicy(const Communication& comm,
                           const FineRangeType& weights,
                           [[maybe_unused]] const PropertyTree& prm,
                           int pressure_var_index)
        : communication_(&const_cast<Communication&>(comm)),
          weights_(weights),
          pressure_var_index_(pressure_var_index)
    {
    }

    void createCoarseLevelSystem(const FineOperator& fineOperator)
    {
        const auto& fineLevelMatrix = fineOperator.getmat();
        // Create coarse level matrix on GPU with block size 1 but same sparsity pattern
        coarseLevelMatrix_ = std::make_shared<GpuMatrixType>(
            fineLevelMatrix.getRowIndices(),
            fineLevelMatrix.getColumnIndices(),
            1 // block size 1
        );

        // Calculate entries for coarse matrix
        calculateCoarseEntries(fineOperator);

        if (!lhs_)
            // We also set lhs to .N() (and not .M()) because we assume a square matrix
            lhs_ = std::make_shared<FineDomainType>(coarseLevelMatrix_->N());
        if (!rhs_)
            rhs_ = std::make_shared<FineRangeType>(coarseLevelMatrix_->N());

        // Create a MatrixAdapter that wraps the matrix without copying it
        operator_ = std::make_shared<CoarseOperator>(*coarseLevelMatrix_);
    }

    void calculateCoarseEntries(const FineOperator& fineOperator)
    {
        const auto& fineLevelMatrix = fineOperator.getmat();

        // Set coarse matrix to zero
        coarseLevelMatrix_->getNonZeroValues() = 0.0;

        // Calculate coarse entries on GPU
        detail::calculateCoarseEntries(
            fineLevelMatrix,
            *coarseLevelMatrix_,
            weights_,
            pressure_var_index_,
            transpose
        );
    }

    void moveToCoarseLevel(const FineRangeType& fine)
    {
        // Set coarse vector to zero
        *rhs_ = 0;

        // Restrict vector on GPU
        detail::restrictVector(
            fine,
            *rhs_,
            weights_,
            pressure_var_index_,
            transpose
        );

        *lhs_ = 0;
    }

    void moveToFineLevel(FineDomainType& fine)
    {
        // Prolongate vector on GPU
        detail::prolongateVector(
            *lhs_,
            fine,
            weights_,
            pressure_var_index_,
            transpose
        );
    }

    GpuPressureTransferPolicy* clone() const
    {
        return new GpuPressureTransferPolicy(*this);
    }

    std::size_t getPressureIndex() const
    {
        return pressure_var_index_;
    }


    std::shared_ptr<CoarseOperator>& getCoarseLevelOperator()
    {
        return operator_;
    }

    CoarseRangeType& getCoarseLevelRhs()
    {
        return *rhs_;
    }

    CoarseDomainType& getCoarseLevelLhs()
    {
        return *lhs_;
    }

    Communication& getCoarseLevelCommunication()
    {
        return *communication_;
    }

    const Communication& getCoarseLevelCommunication() const
    {
        return *communication_;
    }

private:
    Communication* communication_;
    const FineRangeType& weights_;
    const std::size_t pressure_var_index_;
    std::shared_ptr<GpuMatrixType> coarseLevelMatrix_;
    std::shared_ptr<CoarseOperator> operator_;
    std::shared_ptr<FineDomainType> lhs_;
    std::shared_ptr<FineRangeType> rhs_;
};

} // namespace Opm::gpuistl

#endif // OPM_GPU_PRESSURE_TRANSFER_POLICY_HEADER_INCLUDED
