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

#pragma once

#include <opm/common/TimingMacros.hpp>

#include <opm/simulators/linalg/matrixblock.hh>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/twolevelmethodcpr.hh>

#include <dune/istl/paamg/pinfo.hh>

#include <cstddef>

namespace Opm
{
    template <class Communication>
    void extendCommunicatorWithWells(const Communication& comm,
                                     std::shared_ptr<Communication>& commRW,
                                     const int nw)
    {
        OPM_TIMEBLOCK(extendCommunicatorWithWells);
        // used for extending the coarse communicator pattern
        using IndexSet = typename Communication::ParallelIndexSet;
        using LocalIndex = typename IndexSet::LocalIndex;
        const IndexSet& indset = comm.indexSet();
        IndexSet& indset_rw = commRW->indexSet();
        const int max_nw = comm.communicator().max(nw);
        const int rank = comm.communicator().rank();
        int glo_max = 0;
        std::size_t loc_max = 0;
        indset_rw.beginResize();
        for (auto ind = indset.begin(), indend = indset.end(); ind != indend; ++ind) {
            indset_rw.add(ind->global(), LocalIndex(ind->local(), ind->local().attribute(), true));
            const int glo = ind->global();
            const std::size_t loc = ind->local().local();
            glo_max = std::max(glo_max, glo);
            loc_max = std::max(loc_max, loc);
        }
        const int global_max = comm.communicator().max(glo_max);
        // used append the welldofs at the end
        assert(loc_max + 1 == indset.size());
        std::size_t local_ind = loc_max + 1;
        for (int i = 0; i < nw; ++i) {
            // need to set unique global number
            const std::size_t v = global_max + max_nw * rank + i + 1;
            // set to true to not have problems with higher levels if growing of domains is used
            indset_rw.add(v, LocalIndex(local_ind, Dune::OwnerOverlapCopyAttributeSet::owner, true));
            ++local_ind;
        }
        indset_rw.endResize();
        assert(indset_rw.size() == indset.size() + nw);
        // assume same communication pattern
        commRW->remoteIndices().setNeighbours(comm.remoteIndices().getNeighbours());
        commRW->remoteIndices().template rebuild<true>();
        //commRW->clearInterfaces(); may need this for correct rebuild
    }

    namespace Details
    {
        using PressureMatrixType = Dune::BCRSMatrix<Opm::MatrixBlock<double, 1, 1>>;
        using PressureVectorType = Dune::BlockVector<Dune::FieldVector<double, 1>>;
        using SeqCoarseOperatorType = Dune::MatrixAdapter<PressureMatrixType, PressureVectorType, PressureVectorType>;
        template <class Comm>
        using ParCoarseOperatorType
            = Dune::OverlappingSchwarzOperator<PressureMatrixType, PressureVectorType, PressureVectorType, Comm>;
        template <class Comm>
        using CoarseOperatorType = std::conditional_t<std::is_same<Comm, Dune::Amg::SequentialInformation>::value,
                                                      SeqCoarseOperatorType,
                                                      ParCoarseOperatorType<Comm>>;
    } // namespace Details

    template <class FineOperator, class Communication, bool transpose = false>
    class PressureBhpTransferPolicy : public Dune::Amg::LevelTransferPolicyCpr<FineOperator, Details::CoarseOperatorType<Communication>>
    {
    public:
        typedef typename Details::CoarseOperatorType<Communication> CoarseOperator;
        typedef Dune::Amg::LevelTransferPolicyCpr<FineOperator, CoarseOperator> ParentType;
        typedef Communication ParallelInformation;
        typedef typename FineOperator::domain_type FineVectorType;

    public:
        PressureBhpTransferPolicy(const Communication& comm,
                                  const FineVectorType& weights,
                                  const Opm::PropertyTree& prm,
                                  const std::size_t pressureIndex)
            : communication_(&const_cast<Communication&>(comm))
            , weights_(weights)
            , prm_(prm)
            , pressure_var_index_(pressureIndex)
        {
        }

        virtual void createCoarseLevelSystem(const FineOperator& fineOperator) override
        {
            OPM_TIMEBLOCK(createCoarseLevelSystem);
            using CoarseMatrix = typename CoarseOperator::matrix_type;
            const auto& fineLevelMatrix = fineOperator.getmat();
            const auto& nw = fineOperator.getNumberOfExtraEquations();
            if (prm_.get<bool>("add_wells")) {
                const std::size_t average_elements_per_row
                    = static_cast<std::size_t>(std::ceil(fineLevelMatrix.nonzeroes() / fineLevelMatrix.N()));
                const double overflow_fraction = 1.2;
                coarseLevelMatrix_.reset(new CoarseMatrix(fineLevelMatrix.N() + nw,
                                                          fineLevelMatrix.M() + nw,
                                                          average_elements_per_row,
                                                          overflow_fraction,
                                                          CoarseMatrix::implicit));
                int rownum = 0;
                for (const auto& row : fineLevelMatrix) {
                    for (auto col = row.begin(), cend = row.end(); col != cend; ++col) {
                        coarseLevelMatrix_->entry(rownum, col.index()) = 0.0;
                    }
                    ++rownum;
                }
            } else {
                coarseLevelMatrix_.reset(
                    new CoarseMatrix(fineLevelMatrix.N(), fineLevelMatrix.M(), CoarseMatrix::row_wise));
                auto createIter = coarseLevelMatrix_->createbegin();
                for (const auto& row : fineLevelMatrix) {
                    for (auto col = row.begin(), cend = row.end(); col != cend; ++col) {
                        createIter.insert(col.index());
                    }
                    ++createIter;
                }
            }
        if constexpr (std::is_same_v<Communication, Dune::Amg::SequentialInformation>) {
            coarseLevelCommunication_ = std::make_shared<Communication>();
        } else {
            coarseLevelCommunication_ = std::make_shared<Communication>(
                communication_->communicator(), communication_->category(), false);
        }
        if (prm_.get<bool>("add_wells")) {
            fineOperator.addWellPressureEquationsStruct(*coarseLevelMatrix_);
            coarseLevelMatrix_->compress(); // all elemenst should be set
            if constexpr (!std::is_same_v<Communication, Dune::Amg::SequentialInformation>) {
                extendCommunicatorWithWells(*communication_, coarseLevelCommunication_, nw);
            }
        }
        calculateCoarseEntries(fineOperator);

        this->lhs_.resize(this->coarseLevelMatrix_->M());
        this->rhs_.resize(this->coarseLevelMatrix_->N());
        using OperatorArgs = typename Dune::Amg::ConstructionTraits<CoarseOperator>::Arguments;
        OperatorArgs oargs(coarseLevelMatrix_, *coarseLevelCommunication_);
        this->operator_ = Dune::Amg::ConstructionTraits<CoarseOperator>::construct(oargs);
    }

    virtual void calculateCoarseEntries(const FineOperator& fineOperator) override
    {
        OPM_TIMEBLOCK(calculateCoarseEntries);
        const auto& fineMatrix = fineOperator.getmat();
        *coarseLevelMatrix_ = 0;
        auto rowCoarse = coarseLevelMatrix_->begin();
        for (auto row = fineMatrix.begin(), rowEnd = fineMatrix.end(); row != rowEnd; ++row, ++rowCoarse) {
            assert(row.index() == rowCoarse.index());
            auto entryCoarse = rowCoarse->begin();
            for (auto entry = row->begin(), entryEnd = row->end(); entry != entryEnd; ++entry, ++entryCoarse) {
                assert(entry.index() == entryCoarse.index());
                double matrix_el = 0;
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
        if (prm_.get<bool>("add_wells")) {
            OPM_TIMEBLOCK(cprwAddWellEquation);
            assert(transpose == false); // not implemented
            bool use_well_weights = prm_.get<bool>("use_well_weights");
            fineOperator.addWellPressureEquations(*coarseLevelMatrix_, weights_, use_well_weights);
#ifndef NDEBUG
            std::advance(rowCoarse, fineOperator.getNumberOfExtraEquations());
            assert(rowCoarse == coarseLevelMatrix_->end());
#endif

        }
    }

    virtual void moveToCoarseLevel(const typename ParentType::FineRangeType& fine) override
    {
        OPM_TIMEBLOCK(moveToCoarseLevel);
        //NB we iterate over fine assumming welldofs is at the end
        // Set coarse vector to zero
        this->rhs_ = 0;

        auto end = fine.end(), begin = fine.begin();

        for (auto block = begin; block != end; ++block) {
            const auto& bw = weights_[block.index()];
            double rhs_el = 0.0;
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

    virtual void moveToFineLevel(typename ParentType::FineDomainType& fine) override
    {
        OPM_TIMEBLOCK(moveToFineLevel);
        //NB we iterate over fine assumming welldofs is at the end
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

    virtual PressureBhpTransferPolicy* clone() const override
    {
        return new PressureBhpTransferPolicy(*this);
    }

    const Communication& getCoarseLevelCommunication() const
    {
        return *coarseLevelCommunication_;
    }

private:
    Communication* communication_;
    const FineVectorType& weights_;
    PropertyTree prm_;
    const int pressure_var_index_;
    std::shared_ptr<Communication> coarseLevelCommunication_;
    std::shared_ptr<typename CoarseOperator::matrix_type> coarseLevelMatrix_;
};

} // namespace Opm

