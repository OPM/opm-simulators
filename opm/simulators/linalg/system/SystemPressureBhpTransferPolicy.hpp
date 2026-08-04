/*
  Copyright Equinor ASA 2026

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
#ifndef OPM_SYSTEMPRESSUREBHPTRANSFERPOLICY_HEADER_INCLUDED
#define OPM_SYSTEMPRESSUREBHPTRANSFERPOLICY_HEADER_INCLUDED

#include <opm/simulators/linalg/system/SystemCprwPressureStage.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>

#include <opm/simulators/linalg/PressureBhpTransferPolicy.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/twolevelmethodcpr.hh>

#include <dune/istl/paamg/pinfo.hh>

#include <cstddef>
#include <memory>

namespace Opm
{

// --------------------------------------------------------------------------
// The CPRW pressure system of the coupled (reservoir, well) system, presented
// as a Dune level transfer policy so that it can drive TwoLevelMethodCpr.
//
// The assembly and the transfers themselves live in SystemCprwPressureStage;
// this only adapts them to the interface the two-level machinery expects. The
// fine level is the whole system (a SystemVector), the coarse level is the
// scalar pressure system of dimension nCells + nWells.
//
// This is the counterpart of PressureBhpTransferPolicy for the system solver,
// with the essential difference that it needs no well model: the coarse system
// comes from the B/C/D blocks and the weights alone.
// --------------------------------------------------------------------------
template <class FineOperator, class Comm, class Scalar>
class SystemPressureBhpTransferPolicy
    : public Dune::Amg::LevelTransferPolicyCpr<FineOperator, Details::CoarseOperatorType<Scalar, Comm>>
{
public:
    using CoarseOperator = Details::CoarseOperatorType<Scalar, Comm>;
    using ParentType = Dune::Amg::LevelTransferPolicyCpr<FineOperator, CoarseOperator>;
    using ParallelInformation = Comm;
    using Stage = SystemCprwPressureStage<Scalar, Comm>;

    static constexpr auto _0 = Dune::Indices::_0;
    static constexpr auto _1 = Dune::Indices::_1;

    SystemPressureBhpTransferPolicy(const SystemMatrix<Scalar>& S,
                                    const SystemVector<Scalar>& weights,
                                    const PropertyTree& prm,
                                    const int pressureIndex,
                                    const WellTransfer wellTransfer = WellTransfer::Full,
                                    const Comm* comm = nullptr)
        : weights_(weights)
        , stage_(std::make_shared<Stage>(S, prm, pressureIndex, wellTransfer, comm))
    {
    }

    void createCoarseLevelSystem(const FineOperator&) override
    {
        stage_->buildCoarseSystem(weights_);
        const auto& coarse = stage_->coarseMatrixPtr();
        this->lhs_.resize(coarse->M());
        this->rhs_.resize(coarse->N());
        using OperatorArgs = typename Dune::Amg::ConstructionTraits<CoarseOperator>::Arguments;
        OperatorArgs oargs(coarse, stage_->coarseCommunication());
        this->operator_ = Dune::Amg::ConstructionTraits<CoarseOperator>::construct(oargs);
    }

    void calculateCoarseEntries(const FineOperator&) override
    {
        stage_->assembleCoarseEntries(weights_);
    }

    void moveToCoarseLevel(const typename ParentType::FineRangeType& fine) override
    {
        stage_->moveToCoarseLevel(fine[_0], fine[_1], weights_, this->rhs_);
        this->lhs_ = 0;
    }

    void moveToFineLevel(typename ParentType::FineDomainType& fine) override
    {
        stage_->moveToFineLevel(this->lhs_, fine[_0], fine[_1]);
    }

    SystemPressureBhpTransferPolicy* clone() const override
    {
        return new SystemPressureBhpTransferPolicy(*this);
    }

    const Comm& getCoarseLevelCommunication() const
    {
        return stage_->coarseCommunication();
    }

private:
    const SystemVector<Scalar>& weights_;
    // Shared so that clone() is a shallow copy, as it is for
    // PressureBhpTransferPolicy.
    std::shared_ptr<Stage> stage_;
};

} // namespace Opm

#endif // OPM_SYSTEMPRESSUREBHPTRANSFERPOLICY_HEADER_INCLUDED
