/*
  Copyright 2024, SINTEF Digital

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

#ifndef OPM_COMP_WELL_PRIMARY_VARIABLES_HPP
#define OPM_COMP_WELL_PRIMARY_VARIABLES_HPP

#include <opm/material/densead/Evaluation.hpp>

#include <flowexperimental/comp/wells/CompWellEquations.hpp>

#include <array>

namespace Opm {

template <typename FluidSystem, typename Scalar>
class SingleCompWellState;

template <typename FluidSystem, typename Indices>
class CompWellPrimaryVariables
{
public:
    static constexpr int numWellConservationEq = FluidSystem::numComponents;
    static constexpr int numWellControlEq = 1;
    static constexpr int numWellEq = numWellConservationEq + numWellControlEq;
    // the indices for the primary variables
    // the first primary variable will be the total surface rate
    // the last primary variable will be the BHP
    // the one in the middle with will the mole fractions for the numWellEq - 1 components
    // this can be changed based on the implementation itself
    static constexpr int QTotal = 0; // TODO: for now, it is the total surface rate, but later, we might make it total mass rate
    static constexpr int Bhp = numWellEq - numWellControlEq;

    using Scalar = typename FluidSystem::Scalar;

    // this is the number of the equations for the reservoir conservation equations
    static constexpr int numResEq = Indices::numEq;
    // we use hard-coded Evaluation type for now
    // TODO: we can try to use DyanmicEvaluation here
    using EvalWell = DenseAd::Evaluation<Scalar, numWellEq + numResEq>;
    using Eval = DenseAd::Evaluation<Scalar, numResEq>;

    using BVectorWell = typename CompWellEquations<Scalar, numWellEq, numResEq>::BVectorWell;
    using SingleWellState = SingleCompWellState<FluidSystem, Scalar>;

    template <typename T>
    using FluidState = CompositionalFluidState<T, FluidSystem>;

    template <typename T>
    FluidState<T> toFluidState() const;

    void update(const SingleWellState& well_state);

    void updateEvaluation();

    EvalWell getBhp() const;

    EvalWell getTotalRate() const;

    static EvalWell extendEval(const Eval& in);

    static Eval restrictEval(const EvalWell& in);

    void updateNewton(const BVectorWell& dwells);

private:
    std::array<Scalar, numWellEq> value_;
    std::array<EvalWell, numWellEq> evaluation_;
    // temperature for now is constant, so it is not part of the primary variables
    // we need it for the flash calculation
    Scalar temperature_{0.};

    template <typename T>
    T getValue_(int index) const;

};

} // end of namespace Opm

#include "CompWellPrimaryVariables_impl.hpp"

#endif // OPM_COMP_WELL_PRIMARY_VARIABLES_HPP
