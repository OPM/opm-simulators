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


#ifndef OPM_STANDARDWELL_EVAL_HEADER_INCLUDED
#define OPM_STANDARDWELL_EVAL_HEADER_INCLUDED

#include <opm/simulators/wells/StandardWellEquations.hpp>
#include <opm/simulators/wells/StandardWellGeneric.hpp>

#include <opm/material/densead/DynamicEvaluation.hpp>

#include <optional>
#include <vector>

namespace Opm
{

class ConvergenceReport;
class DeferredLogger;
class GroupState;
class Schedule;
class SummaryState;
class WellContributions;
template<class FluidSystem, class Indices, class Scalar> class WellInterfaceIndices;
class WellState;

template<class FluidSystem, class Indices, class Scalar>
class StandardWellEval : public StandardWellGeneric<Scalar>
{
protected:
    using BVectorWell = typename StandardWellEquations<Indices,Scalar>::BVectorWell;
    static constexpr int numStaticWellEq = StandardWellEquations<Indices,Scalar>::numStaticWellEq;
    static constexpr int Bhp = StandardWellEquations<Indices,Scalar>::Bhp;

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

    static const int WQTotal = 0; 
    // In the implementation, one should use has_wfrac_variable
    // rather than has_water to check if you should do something
    // with the variable at the WFrac location, similar for GFrac.
    // (following implementation MultisegmentWellEval.hpp)
    static const bool waterEnabled = Indices::waterEnabled;
    static const bool gasEnabled = Indices::gasEnabled;
    static const bool oilEnabled = Indices::oilEnabled;

    static constexpr bool has_wfrac_variable = Indices::waterEnabled && Indices::oilEnabled;
    static constexpr bool has_gfrac_variable = Indices::gasEnabled && Indices::numPhases > 1;
    static constexpr int WFrac = has_wfrac_variable ? 1 : -1000;
    static constexpr int GFrac = has_gfrac_variable ? has_wfrac_variable + 1 : -1000;
    static constexpr int SFrac = !Indices::enableSolvent ? -1000 : 3;

public:
    using EvalWell = DenseAd::DynamicEvaluation<Scalar, numStaticWellEq + Indices::numEq + 1>;
    using Eval = DenseAd::Evaluation<Scalar, Indices::numEq>;

    /// get the number of blocks of the C and B matrices, used to allocate memory in a WellContributions object
    unsigned int getNumBlocks() const
    {
        return this->linSys_.getNumBlocks();
    }

    /// add the contribution (C, D^-1, B matrices) of this Well to the WellContributions object
    void addWellContribution(WellContributions& wellContribs) const;

protected:
    StandardWellEval(const WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif);

    const WellInterfaceIndices<FluidSystem,Indices,Scalar>& baseif_;

    void initPrimaryVariablesEvaluation() const;

    const EvalWell& getBhp() const
    {
        return primary_variables_evaluation_[Bhp];
    }

    const EvalWell& getWQTotal() const
    {
        return primary_variables_evaluation_[WQTotal];
    }

    EvalWell extendEval(const Eval& in) const;
    EvalWell getQs(const int compIdx) const;
    EvalWell wellSurfaceVolumeFraction(const int compIdx) const;
    EvalWell wellVolumeFraction(const unsigned compIdx) const;
    EvalWell wellVolumeFractionScaled(const int phase) const;

    // calculate a relaxation factor to avoid overshoot of the fractions for producers
    // which might result in negative rates
    static double relaxationFactorFractionsProducer(const std::vector<double>& primary_variables,
                                                    const BVectorWell& dwells);

    // computing the accumulation term for later use in well mass equations
    void computeAccumWell();

    // TODO: not total sure whether it is a good idea to put this function here
    // the major reason to put here is to avoid the usage of Wells struct
    void computeConnectionDensities(const std::vector<double>& perfComponentRates,
                                    const std::vector<double>& b_perf,
                                    const std::vector<double>& rsmax_perf,
                                    const std::vector<double>& rvmax_perf,
                                    const std::vector<double>& rvwmax_perf,
                                    const std::vector<double>& surf_dens_perf,
                                    DeferredLogger& deferred_logger);

    ConvergenceReport getWellConvergence(const WellState& well_state,
                                         const std::vector<double>& B_avg,
                                         const double maxResidualAllowed,
                                         const double tol_wells,
                                         const double relaxed_tolerance_flow,
                                         const bool relax_tolerance,
                                         std::vector<double>& res,
                                         DeferredLogger& deferred_logger) const;

    void init(std::vector<double>& perf_depth,
              const std::vector<double>& depth_arg,
              const int num_cells,
              const bool has_polymermw);

    // handle the non reasonable fractions due to numerical overshoot
    void processFractions() const;

    void updatePrimaryVariables(const WellState& well_state,
                                DeferredLogger& deferred_logger) const;

    void updatePrimaryVariablesPolyMW(const BVectorWell& dwells) const;

    void updateWellStateFromPrimaryVariables(WellState& well_state,
                                             DeferredLogger& deferred_logger) const;

    void updatePrimaryVariablesNewton(const BVectorWell& dwells,
                                      const double dFLimit,
                                      const double dBHPLimit) const;

    void updateWellStateFromPrimaryVariablesPolyMW(WellState& well_state) const;

    // total number of the well equations and primary variables
    // there might be extra equations be used, numWellEq will be updated during the initialization
    int numWellEq_ = numStaticWellEq;

    // the values for the primary varibles
    // based on different solutioin strategies, the wells can have different primary variables
    mutable std::vector<double> primary_variables_;

    // the Evaluation for the well primary variables, which contain derivativles and are used in AD calculation
    mutable std::vector<EvalWell> primary_variables_evaluation_;

    // the saturations in the well bore under surface conditions at the beginning of the time step
    std::vector<double> F0_;

    StandardWellEquations<Indices,Scalar> linSys_; //!< Equation system
};

}

#endif // OPM_STANDARDWELL_EVAL_HEADER_INCLUDED
