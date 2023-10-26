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


#ifndef OPM_STANDARDWELL_PRIMARY_VARIABLES_HEADER_INCLUDED
#define OPM_STANDARDWELL_PRIMARY_VARIABLES_HEADER_INCLUDED

#include <opm/material/densead/Evaluation.hpp>

#include <opm/simulators/wells/StandardWellEquations.hpp>

#include <vector>

namespace Opm
{

class DeferredLogger;
template<class FluidSystem, class Indices, class Scalar> class WellInterfaceIndices;
class WellState;

//! \brief Class holding primary variables for StandardWell.
template<class FluidSystem, class Indices, class Scalar>
class StandardWellPrimaryVariables {
protected:
    // the positions of the primary variables for StandardWell
    // the first one is the weighted total rate (WQ_t), the second and the third ones are F_w and F_g,
    // which represent the fraction of Water and Gas based on the weighted total rate, the last one is BHP.
    // correspondingly, we have four well equations for blackoil model, the first three are mass
    // converstation equations, and the last one is the well control equation.
    // primary variables related to other components, will be before the Bhp and after F_g.
    // well control equation is always the last well equation.
    // TODO: in the current implementation, we use the well rate as the first primary variables for injectors,
    // instead of G_t.

    // Table showing the primary variable indices, depending on what phases are present:
    //
    //         WOG     OG     WG     WO    W/O/G (single phase)
    // WQTotal   0      0      0      0                       0
    // WFrac     1  -1000     -1000   1                   -1000
    // GFrac     2      1      1  -1000                   -1000
    // Spres     3      2      2      2                       1

    //! \brief Number of the well control equations.
    static constexpr int numWellControlEq = 1;

public:
    //! \brief Number of the conservation equations.
    static constexpr int numWellConservationEq = Indices::numPhases + Indices::numSolvents;

    //! \brief Number of the well equations that will always be used.
    //! \details Based on the solution strategy, there might be other well equations be introduced.
    static constexpr int numStaticWellEq = numWellConservationEq + numWellControlEq;

    static constexpr int WQTotal = 0; //!< The index for the weighted total rate

    //! \brief The index for Bhp in primary variables and the index of well control equation.
    //! \details They both will be the last one in their respective system.
    //! \todo: We should have indices for the well equations and well primary variables separately.
    static constexpr int Bhp = numStaticWellEq - numWellControlEq;

    static constexpr bool has_wfrac_variable = Indices::waterEnabled && Indices::oilEnabled;
    static constexpr bool has_gfrac_variable = Indices::gasEnabled && Indices::numPhases > 1;
    static constexpr int WFrac = has_wfrac_variable ? 1 : -1000;
    static constexpr int GFrac = has_gfrac_variable ? has_wfrac_variable + 1 : -1000;
    static constexpr int SFrac = !Indices::enableSolvent ? -1000 : has_wfrac_variable+has_gfrac_variable+1;

    //! \brief Evaluation for the well equations.
    using EvalWell = DenseAd::DynamicEvaluation<Scalar, numStaticWellEq + Indices::numEq + 1>;
    using BVectorWell = typename StandardWellEquations<Scalar,Indices::numEq>::BVectorWell;

    //! \brief Constructor initializes reference to well interface.
    StandardWellPrimaryVariables(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well)
        : well_(well)
    {}

    //! \brief Initialize evaluations from values.
    void init();

    //! \brief Resize values and evaluations.
    void resize(const int numWellEq);

    //! \brief Returns number of well equations.
    int numWellEq() const { return numWellEq_; }

    //! \brief Copy values from well state.
    void update(const WellState& well_state,
                const bool stop_or_zero_rate_target,
                DeferredLogger& deferred_logger);

    //! \brief Copy polymer molecular weigt values from well state.
    void updatePolyMW(const WellState& well_state);

    //! \brief Update values from newton update vector.
    void updateNewton(const BVectorWell& dwells,
                      const bool stop_or_zero_rate_target,
                      const double dFLimit,
                      const double dBHPLimit,
                      DeferredLogger& deferred_logger);

    //! \brief Update polymer molecular weight values from newton update vector.
    void updateNewtonPolyMW(const BVectorWell& dwells);

    //! \brief Check that all values are finite.
    void checkFinite(DeferredLogger& deferred_logger) const;

    //! \brief Copy values to well state.
    void copyToWellState(WellState& well_state, DeferredLogger& deferred_logger) const;

    //! \brief Copy polymer molecular weight values to well state.
    void copyToWellStatePolyMW(WellState& well_state) const;

    //! \brief Returns scaled volume fraction for a component.
    EvalWell volumeFractionScaled(const int compIdx) const;

    //! \brief Returns surface volume fraction for a component.
    EvalWell surfaceVolumeFraction(const int compIdx) const;

    //! \brief Returns scaled rate for a component.
    EvalWell getQs(const int compIdx) const;

    //! \brief Returns a value.
    Scalar value(const int idx) const
    { return value_[idx]; }

    //! \brief Returns a const ref to an evaluation.
    const EvalWell& eval(const int idx) const
    { return evaluation_[idx]; }

    //! \brief Set a value. Note that this does not also set the corresponding evaluation.
    void setValue(const int idx, const Scalar val)
    { value_[idx] = val; }

private:
    //! \brief Calculate a relaxation factor for producers.
    //! \details To avoid overshoot of the fractions which might result in negative rates.
    double relaxationFactorFractionsProducer(const BVectorWell& dwells, DeferredLogger& deferred_logger) const;

    //! \brief Returns volume fraction for a component.
    EvalWell volumeFraction(const unsigned compIdx) const;

    //! \brief Handle non-reasonable fractions due to numerical overshoot.
    void processFractions();

    //! \brief The values for the primary variables.
    //! \details Based on different solution strategies, the wells can have different primary variables.
    std::vector<Scalar> value_;

    //! \brief The Evaluation for the well primary variables.
    //! \details Contain derivatives and are used in AD calculation
    std::vector<EvalWell> evaluation_;

    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& well_; //!< Reference to well interface

    //! \brief Total number of the well equations and primary variables.
    //! \details There might be extra equations be used, numWellEq will be updated during the initialization
    int numWellEq_ = numStaticWellEq;
};

}

#endif // OPM_STANDARDWELL_PRIMARY_VARIABLES_HEADER_INCLUDED
