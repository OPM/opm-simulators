// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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
#ifndef OPM_TPSA_PRECONDITIONER_IMPL_HPP
#define OPM_TPSA_PRECONDITIONER_IMPL_HPP

#ifndef OPM_TPSA_PRECONDITIONER_HPP
#include <config.h>
#include <opm/simulators/linalg/tpsa/TpsaPreconditioner.hpp>
#endif

namespace Opm
{

template <class Scalar, class DispOp, class RotOp, class SPresOp, class Comm>
TpsaPreconditioner<Scalar, DispOp, RotOp, SPresOp, Comm>::
TpsaPreconditioner(const Linear::TpsaMatrixView<Scalar>& S, const PropertyTree& prm)
    requires (!isParallel)
    : S_(S)
{
    initSubSolvers_(prm);
}

template <class Scalar, class DispOp, class RotOp, class SPresOp, class Comm>
TpsaPreconditioner<Scalar, DispOp, RotOp, SPresOp, Comm>::
TpsaPreconditioner(const Linear::TpsaMatrixView<Scalar>& S,
                   const PropertyTree& prm,
                   const Comm& comm)
    requires (isParallel)
    : S_(S)
    , comm_(&comm)
{
    initSubSolvers_(prm);
}

template <class Scalar, class DispOp, class RotOp, class SPresOp, class Comm>
void
TpsaPreconditioner<Scalar, DispOp, RotOp, SPresOp, Comm>::
pre(MultiVector&, MultiVector&)
{
}

template <class Scalar, class DispOp, class RotOp, class SPresOp, class Comm>
void
TpsaPreconditioner<Scalar, DispOp, RotOp, SPresOp, Comm>::
post(MultiVector&)
{
}

template <class Scalar, class DispOp, class RotOp, class SPresOp, class Comm>
Dune::SolverCategory::Category
TpsaPreconditioner<Scalar, DispOp, RotOp, SPresOp, Comm>::
category() const
{
    if constexpr (isParallel) {
        return Dune::SolverCategory::overlapping;
    } else {
        return Dune::SolverCategory::sequential;
    }
}

template <class Scalar, class DispOp, class RotOp, class SPresOp, class Comm>
void
TpsaPreconditioner<Scalar, DispOp, RotOp, SPresOp, Comm>::
update()
{
    dispSolver0_->preconditioner().update();
    dispSolver1_->preconditioner().update();
    dispSolver2_->preconditioner().update();
    rotSolver_->preconditioner().update();
    sPresSolver_->preconditioner().update();
}

template <class Scalar, class DispOp, class RotOp, class SPresOp, class Comm>
bool
TpsaPreconditioner<Scalar, DispOp, RotOp, SPresOp, Comm>::
hasPerfectUpdate() const
{
    return true;
}

template <class Scalar, class DispOp, class RotOp, class SPresOp, class Comm>
void
TpsaPreconditioner<Scalar, DispOp, RotOp, SPresOp, Comm>::
apply(MultiVector& v, const MultiVector& d)
{
    Dune::InverseOperatorResult result;

    // The defects of the coupled fields are updated as the sweep proceeds,
    // so they need their own copies.
    auto d0 = d[_0];
    auto d1 = d[_1];
    auto d2 = d[_2];
    auto d3 = d[_3];
    auto d4 = d[_4];

    // Ensure that v is zero-initialized for apply()
    v = 0.0;

    dispSolver0_->apply(v[_0], d0, result);
    dispSolver1_->apply(v[_1], d1, result);
    dispSolver2_->apply(v[_2], d2, result);

    S_[_3][_0].mmv(v[_0], d3);
    S_[_3][_1].mmv(v[_1], d3);
    S_[_3][_2].mmv(v[_2], d3);
    rotSolver_->apply(v[_3], d3, result);

    S_[_4][_0].mmv(v[_0], d4);
    S_[_4][_1].mmv(v[_1], d4);
    S_[_4][_2].mmv(v[_2], d4);
    sPresSolver_->apply(v[_4], d4, result);
}

template <class Scalar, class DispOp, class RotOp, class SPresOp, class Comm>
void
TpsaPreconditioner<Scalar, DispOp, RotOp, SPresOp, Comm>::
initSubSolvers_(const PropertyTree& prm)
{
    const auto dispPrm = prm.get_child("disp_disp_solver");
    const auto rotPrm = prm.get_child("rot_rot_solver");
    const auto sPresPrm = prm.get_child("spres_spres_solver");

    std::function<Linear::DispVector0T<Scalar>()> dispWeightCalc;
    std::function<Linear::RotVectorT<Scalar>()> rotWeightCalc;
    std::function<Linear::SPresVectorT<Scalar>()> sPresWeightCalc;

    if constexpr (isParallel) {
        dispOp0_ = std::make_unique<DispOp>(S_[_0][_0], *comm_);
        dispSolver0_ = std::make_unique<DispSolver>(*dispOp0_,
                                                    *comm_,
                                                    dispPrm,
                                                    dispWeightCalc,
                                                    pressureIdx);

        dispOp1_ = std::make_unique<DispOp>(S_[_1][_1], *comm_);
        dispSolver1_ = std::make_unique<DispSolver>(*dispOp1_,
                                                    *comm_,
                                                    dispPrm,
                                                    dispWeightCalc,
                                                    pressureIdx);

        dispOp2_ = std::make_unique<DispOp>(S_[_2][_2], *comm_);
        dispSolver2_ = std::make_unique<DispSolver>(*dispOp2_,
                                                    *comm_,
                                                    dispPrm,
                                                    dispWeightCalc,
                                                    pressureIdx);

        rotOp_ = std::make_unique<RotOp>(S_[_3][_3], *comm_);
        rotSolver_ = std::make_unique<RotSolver>(*rotOp_,
                                                 *comm_,
                                                 rotPrm,
                                                 rotWeightCalc,
                                                 pressureIdx);

        sPresOp_ = std::make_unique<SPresOp>(S_[_4][_4], *comm_);
        sPresSolver_ = std::make_unique<SPresSolver>(*sPresOp_,
                                                     *comm_,
                                                     sPresPrm,
                                                     sPresWeightCalc,
                                                     pressureIdx);
    } else {
        dispOp0_ = std::make_unique<DispOp>(S_[_0][_0]);
        dispSolver0_ = std::make_unique<DispSolver>(*dispOp0_,
                                                    dispPrm,
                                                    dispWeightCalc,
                                                    pressureIdx);

        dispOp1_ = std::make_unique<DispOp>(S_[_1][_1]);
        dispSolver1_ = std::make_unique<DispSolver>(*dispOp1_,
                                                    dispPrm,
                                                    dispWeightCalc,
                                                    pressureIdx);

        dispOp2_ = std::make_unique<DispOp>(S_[_2][_2]);
        dispSolver2_ = std::make_unique<DispSolver>(*dispOp2_,
                                                    dispPrm,
                                                    dispWeightCalc,
                                                    pressureIdx);

        rotOp_ = std::make_unique<RotOp>(S_[_3][_3]);
        rotSolver_ = std::make_unique<RotSolver>(*rotOp_,
                                                 rotPrm,
                                                 rotWeightCalc,
                                                 pressureIdx);

        sPresOp_ = std::make_unique<SPresOp>(S_[_4][_4]);
        sPresSolver_ = std::make_unique<SPresSolver>(*sPresOp_,
                                                     sPresPrm,
                                                     sPresWeightCalc,
                                                     pressureIdx);
    }
}

} // namespace Opm

#endif // OPM_TPSA_PRECONDITIONER_IMPL_HPP
