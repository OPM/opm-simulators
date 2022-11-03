/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2019 Dr. Blatt - HPC-Simulation-Software & Services.

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

#ifndef OPM_PRESSURE_TRANSFER_POLICY_HEADER_INCLUDED
#define OPM_PRESSURE_TRANSFER_POLICY_HEADER_INCLUDED


#include <opm/simulators/linalg/twolevelmethodcpr.hh>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/WellOperators.hpp>

#include <cstddef>

namespace Opm { namespace Details {
    template<class Scalar> using PressureMatrixType = Dune::BCRSMatrix<MatrixBlock<Scalar, 1, 1>>;
    template<class Scalar> using PressureVectorType = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;
    template<class Scalar> using SeqCoarseOperatorType = Dune::MatrixAdapter<PressureMatrixType<Scalar>,
                                                                             PressureVectorType<Scalar>,
                                                                             PressureVectorType<Scalar>>;
    template<class Scalar, class Comm>
    using ParCoarseOperatorType
        = Opm::GhostLastMatrixAdapter<PressureMatrixType<Scalar>,
                                       PressureVectorType<Scalar>,
                                       PressureVectorType<Scalar>,
                                       Comm>;
    template<class Scalar, class Comm>
    using CoarseOperatorType = std::conditional_t<std::is_same<Comm, Dune::Amg::SequentialInformation>::value,
                                                  SeqCoarseOperatorType<Scalar>,
                                                  ParCoarseOperatorType<Scalar,Comm>>;
} // namespace Details



template <class FineOperator, class Communication, class Scalar, bool transpose = false>
class PressureTransferPolicy
    : public Dune::Amg::LevelTransferPolicyCpr<FineOperator, Details::CoarseOperatorType<Scalar,Communication>>
{
public:
    using CoarseOperator = typename Details::CoarseOperatorType<Scalar,Communication>;
    using ParentType = Dune::Amg::LevelTransferPolicyCpr<FineOperator, CoarseOperator>;
    using ParallelInformation = Communication;
    using FineVectorType = typename FineOperator::domain_type;

public:
    PressureTransferPolicy(const Communication& comm,
                           const FineVectorType& weights,
                           const PropertyTree& /*prm*/,
                           int pressure_var_index)
        : communication_(&const_cast<Communication&>(comm))
        , weights_(weights)
        , pressure_var_index_(pressure_var_index)
    {
    }

    void createCoarseLevelSystem(const FineOperator& fineOperator) override
    {
        using CoarseMatrix = typename CoarseOperator::matrix_type;
        const auto& fineLevelMatrix = fineOperator.getmat();
        coarseLevelMatrix_.reset(new CoarseMatrix(fineLevelMatrix.N(), fineLevelMatrix.M(), CoarseMatrix::row_wise));
        auto createIter = coarseLevelMatrix_->createbegin();

        for (const auto& row : fineLevelMatrix) {
            for (auto col = row.begin(), cend = row.end(); col != cend; ++col) {
                createIter.insert(col.index());
            }
            ++createIter;
        }

        calculateCoarseEntries(fineOperator);
        coarseLevelCommunication_.reset(communication_, [](Communication*) {});

        this->lhs_.resize(this->coarseLevelMatrix_->M());
        this->rhs_.resize(this->coarseLevelMatrix_->N());
        using OperatorArgs = typename Dune::Amg::ConstructionTraits<CoarseOperator>::Arguments;
        OperatorArgs oargs(coarseLevelMatrix_, *coarseLevelCommunication_);
        this->operator_ = Dune::Amg::ConstructionTraits<CoarseOperator>::construct(oargs);
    }

    void calculateCoarseEntries(const FineOperator& fineOperator) override
    {
        const auto& fineMatrix = fineOperator.getmat();
        *coarseLevelMatrix_ = 0;
        auto rowCoarse = coarseLevelMatrix_->begin();
        for (auto row = fineMatrix.begin(), rowEnd = fineMatrix.end(); row != rowEnd; ++row, ++rowCoarse) {
            assert(row.index() == rowCoarse.index());
            auto entryCoarse = rowCoarse->begin();
            for (auto entry = row->begin(), entryEnd = row->end(); entry != entryEnd; ++entry, ++entryCoarse) {
                assert(entry.index() == entryCoarse.index());
                Scalar matrix_el = 0;
                if (transpose) {
                    const auto& bw = weights_[entry.index()];
                    for (std::size_t i = 0; i < bw.size(); ++i) {
                        matrix_el += (*entry)[pressure_var_index_][i] * bw[i];
                    }
                } else {
                    const auto& bw = weights_[row.index()];
                    for (std::size_t i = 0; i < bw.size(); ++i) {
                        matrix_el += (*entry)[i][pressure_var_index_] * bw[i];
                    }
                }
                (*entryCoarse) = matrix_el;
            }
        }
        assert(rowCoarse == coarseLevelMatrix_->end());
    }

    void moveToCoarseLevel(const typename ParentType::FineRangeType& fine) override
    {
        // Set coarse vector to zero
        this->rhs_ = 0;

        auto end = fine.end(), begin = fine.begin();

        for (auto block = begin; block != end; ++block) {
            const auto& bw = weights_[block.index()];
            Scalar rhs_el = 0.0;
            if (transpose) {
                rhs_el = (*block)[pressure_var_index_];
            } else {
                for (std::size_t i = 0; i < block->size(); ++i) {
                    rhs_el += (*block)[i] * bw[i];
                }
            }
            this->rhs_[block - begin] = rhs_el;
        }

        this->lhs_ = 0;
    }

    void moveToFineLevel(typename ParentType::FineDomainType& fine) override
    {
        auto end = fine.end(), begin = fine.begin();

        for (auto block = begin; block != end; ++block) {
            if (transpose) {
                const auto& bw = weights_[block.index()];
                for (std::size_t i = 0; i < block->size(); ++i) {
                    (*block)[i] = this->lhs_[block - begin] * bw[i];
                }
            } else {
                (*block)[pressure_var_index_] = this->lhs_[block - begin];
            }
        }
    }

    PressureTransferPolicy* clone() const override
    {
        return new PressureTransferPolicy(*this);
    }

    const Communication& getCoarseLevelCommunication() const
    {
        return *coarseLevelCommunication_;
    }

    std::size_t getPressureIndex() const
    {
        return pressure_var_index_;
    }

private:
    Communication* communication_;
    const FineVectorType& weights_;
    const std::size_t pressure_var_index_;
    std::shared_ptr<Communication> coarseLevelCommunication_;
    std::shared_ptr<typename CoarseOperator::matrix_type> coarseLevelMatrix_;
};

} // namespace Opm

#endif // OPM_PRESSURE_TRANSFER_POLICY_HEADER_INCLUDED
