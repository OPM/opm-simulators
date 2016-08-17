/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015 IRIS AS

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

#ifndef OPM_BLACKOILMODELEBOS_HEADER_INCLUDED
#define OPM_BLACKOILMODELEBOS_HEADER_INCLUDED

#include <applications/ebos/eclproblem.hh>
#include <ewoms/common/start.hh>

#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/StandardWellsDense.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/WellHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/autodiff/VFPProperties.hpp>
#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/autodiff/VFPInjProperties.hpp>
#include <opm/autodiff/DefaultBlackoilSolutionState.hpp>
#include <opm/autodiff/BlackoilDetails.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>

#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/well_controls.h>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/common/data/SimulationDataContainer.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <algorithm>
//#include <fstream>



namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(EclFlowProblem, INHERITS_FROM(BlackOilModel, EclBaseProblem));
SET_BOOL_PROP(EclFlowProblem, DisableWells, true);
SET_BOOL_PROP(EclFlowProblem, EnableDebuggingChecks, false);
}}

namespace Opm {


    namespace parameter { class ParameterGroup; }
    class DerivedGeology;
    class RockCompressibility;
    class NewtonIterationBlackoilInterface;
    class VFPProperties;
    class SimulationDataContainer;




    /// A model implementation for three-phase black oil.
    ///
    /// The simulator is capable of handling three-phase problems
    /// where gas can be dissolved in oil and vice versa. It
    /// uses an industry-standard TPFA discretization with per-phase
    /// upwind weighting of mobilities.
    class BlackoilModelEbos
    {
    public:
        // ---------  Types and enums  ---------
        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;

        typedef BlackoilState ReservoirState;
        typedef WellStateFullyImplicitBlackoil WellState;
        typedef BlackoilModelParameters ModelParameters;
        typedef DefaultBlackoilSolutionState SolutionState;

        typedef typename TTAG(EclFlowProblem) TypeTag;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator)         Simulator ;
        typedef typename GET_PROP_TYPE(TypeTag, Grid)              Grid;
        typedef typename GET_PROP_TYPE(TypeTag, SolutionVector)    SolutionVector ;
        typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables)  PrimaryVariables ;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem)       FluidSystem;
        typedef typename GET_PROP_TYPE(TypeTag, Indices)           BlackoilIndices;
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw)       MaterialLaw;
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
        //typedef typename SolutionVector :: value_type            PrimaryVariables ;

        // ---------  Public methods  ---------

        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] fluid            fluid properties
        /// \param[in] geo              rock properties
        /// \param[in] rock_comp_props  if non-null, rock compressibility properties
        /// \param[in] wells            well structure
        /// \param[in] vfp_properties   Vertical flow performance tables
        /// \param[in] linsolver        linear solver
        /// \param[in] eclState         eclipse state
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilModelEbos(Simulator& ebosSimulator,
                          const ModelParameters&          param,
                          const BlackoilPropsAdInterface& fluid,
                          const DerivedGeology&           geo  ,
                          const RockCompressibility*      rock_comp_props,
                          const StandardWellsDense&                well_model,
                          const NewtonIterationBlackoilInterface& linsolver,
                          const bool terminal_output)
        : ebosSimulator_(ebosSimulator)
        , grid_(ebosSimulator_.gridManager().grid())
        , fluid_ (fluid)
        , geo_   (geo)
        , rock_comp_props_(rock_comp_props)
        , vfp_properties_(
            eclState().getTableManager().getVFPInjTables(),
            eclState().getTableManager().getVFPProdTables())
        , linsolver_ (linsolver)
        , active_(detail::activePhases(fluid.phaseUsage()))
        , canph_ (detail::active2Canonical(fluid.phaseUsage()))
        , cells_ (detail::buildAllCells(Opm::AutoDiffGrid::numCells(grid_)))
        , ops_   (grid_, geo.nnc())
        , has_disgas_(FluidSystem::enableDissolvedGas())
        , has_vapoil_(FluidSystem::enableVaporizedOil())
        , param_( param )
        , use_threshold_pressure_(false)
        , rq_    (fluid.numPhases())
        , phaseCondition_(AutoDiffGrid::numCells(grid_))
        , well_model_ (well_model)
        , isRs_(V::Zero(AutoDiffGrid::numCells(grid_)))
        , isRv_(V::Zero(AutoDiffGrid::numCells(grid_)))
        , isSg_(V::Zero(AutoDiffGrid::numCells(grid_)))
        , residual_ ( { std::vector<ADB>(fluid.numPhases(), ADB::null()),
                    ADB::null(),
                    ADB::null(),
                    { 1.1169, 1.0031, 0.0031 }, // the default magic numbers
                    false } )
        , terminal_output_ (terminal_output)
        , current_relaxation_(1.0)
        , isBeginReportStep_(false)
        , wellVariables_(wells().number_of_wells * wells().number_of_phases)
        , F0_(wells().number_of_wells * wells().number_of_phases)
        {
            const double gravity = detail::getGravity(geo_.gravity(), UgGridHelpers::dimensions(grid_));
            const V depth = Opm::AutoDiffGrid::cellCentroidsZToEigen(grid_);

            well_model_.init(&fluid_, &active_, &phaseCondition_, &vfp_properties_, gravity, depth);

            wellModel().setWellsActive( localWellsActive() );
            global_nc_    =  Opm::AutoDiffGrid::numCells(grid_);
        }

        const EclipseState& eclState() const
        { return *ebosSimulator_.gridManager().eclState(); }

        /// Called once before each time step.
        /// \param[in] timer                  simulation timer
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void prepareStep(const SimulatorTimerInterface& timer,
                         const ReservoirState& reservoir_state,
                         const WellState& /* well_state */)
        {
            if (active_[Gas]) {
                updatePrimalVariableFromState(reservoir_state);
            }
        }


        /// Called once per nonlinear iteration.
        /// This model will perform a Newton-Raphson update, changing reservoir_state
        /// and well_state. It will also use the nonlinear_solver to do relaxation of
        /// updates if necessary.
        /// \param[in] iteration              should be 0 for the first call of a new timestep
        /// \param[in] timer                  simulation timer
        /// \param[in] nonlinear_solver       nonlinear solver used (for oscillation/relaxation control)
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        template <class NonlinearSolverType>
        IterationReport nonlinearIteration(const int iteration,
                                           const SimulatorTimerInterface& timer,
                                           NonlinearSolverType& nonlinear_solver,
                                           ReservoirState& reservoir_state,
                                           WellState& well_state)
        {
            const double dt = timer.currentStepLength();
            if (iteration == 0) {
                // For each iteration we store in a vector the norms of the residual of
                // the mass balance for each active phase, the well flux and the well equations.
                residual_norms_history_.clear();
                current_relaxation_ = 1.0;
                dx_old_ = V::Zero(sizeNonLinear());
            }
            IterationReport iter_report = assemble(timer, iteration, reservoir_state, well_state);
            residual_norms_history_.push_back(computeResidualNorms());
            const bool converged = getConvergence(timer, iteration);
            const bool must_solve = (iteration < nonlinear_solver.minIter()) || (!converged);
            if (must_solve) {
                // enable single precision for solvers when dt is smaller then 20 days
                residual_.singlePrecision = (unit::convert::to(dt, unit::day) < 20.) ;

                // Compute the nonlinear update.
                V dx = solveJacobianSystem();

                // Stabilize the nonlinear update.
                bool isOscillate = false;
                bool isStagnate = false;
                nonlinear_solver.detectOscillations(residual_norms_history_, iteration, isOscillate, isStagnate);
                if (isOscillate) {
                    current_relaxation_ -= nonlinear_solver.relaxIncrement();
                    current_relaxation_ = std::max(current_relaxation_, nonlinear_solver.relaxMax());
                    if (terminalOutputEnabled()) {
                        std::string msg = " Oscillating behavior detected: Relaxation set to "
                            + std::to_string(current_relaxation_);
                        OpmLog::info(msg);
                    }
                }
                nonlinear_solver.stabilizeNonlinearUpdate(dx, dx_old_, current_relaxation_);

                // Apply the update, applying model-dependent
                // limitations and chopping of the update.
                updateState(dx, reservoir_state, well_state);
            }
            const bool failed = false; // Not needed in this model.
            const int linear_iters = must_solve ? linearIterationsLastSolve() : 0;
            return IterationReport{ failed, converged, linear_iters, iter_report.well_iterations };
        }


        /// Called once after each time step.
        /// In this class, this function does nothing.
        /// \param[in] timer                  simulation timer
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void afterStep(const SimulatorTimerInterface& timer,
                       const ReservoirState& reservoir_state,
                       WellState& well_state)
        {
        }

        /// Assemble the residual and Jacobian of the nonlinear system.
        /// \param[in]      reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        /// \param[in]      initial_assembly  pass true if this is the first call to assemble() in this timestep
        IterationReport assemble(const SimulatorTimerInterface& timer,
                                 const int iterationIdx,
                                 const ReservoirState& reservoir_state,
                                 WellState& well_state)
        {
            using namespace Opm::AutoDiffGrid;

            // Possibly switch well controls and updating well state to
            // get reasonable initial conditions for the wells
            wellModel().updateWellControls(well_state);

            // Create the primary variables.
            setWellVariables(well_state);

            //SolutionState state(/*numPhases=*/3);
            //setupLegacyState(state, reservoir_state, well_state);

            // -------- Mass balance equations --------
            assembleMassBalanceEq(timer, iterationIdx, reservoir_state);

            // -------- Well equations ----------

            if (iterationIdx == 0) {
                // Create the (constant, derivativeless) initial state.
                // Compute initial accumulation contributions
                // and well connection pressures.
                computeWellConnectionPressures(well_state);

                computeAccumWells();
                // Create the (constant, derivativeless) initial state.
                //SolutionState state0 = state;
                //makeConstantState(state0);
                // Compute initial accumulation contributions
                // and well connection pressures.
                //computeWellConnectionPressures(state0, well_state);
                //wellModel().computeAccumWells(state0);
            }

            IterationReport iter_report = {false, false, 0, 0};
            if ( ! wellsActive() ) {
                return iter_report;
            }


            double dt = timer.currentStepLength();
            std::vector<ADB> mob_perfcells;
            std::vector<ADB> b_perfcells;
            //wellModel().extractWellPerfProperties(state, rq_, mob_perfcells, b_perfcells);
            if (param_.solve_welleq_initially_ && iterationIdx == 0) {
                // solve the well equations as a pre-processing step
                iter_report = solveWellEq(mob_perfcells, b_perfcells, dt, well_state);
            }
            std::vector<ADB> cq_s;                
            computeWellFluxDense(cq_s, 4);
            updatePerfPhaseRatesAndPressures(cq_s, well_state);

            const int np = wells().number_of_phases;
            const int nw = wells().number_of_wells;
                for (int p = 0; p < np; ++p) {
                    for (int w = 0; w < nw; ++w) {
                    //std::cout << p << " " << w << " " << getQs(w,p).value<<std::endl;
                    //print( getQs(w, p));
                }
            }
            //std::cout << state.qs.value() << std::endl;

            //wellModel().addWellFluxEq(cq_s, state, dt, residual_);
            addWellFluxEq(cq_s, dt, 4);
            addWellContributionToMassBalanceEq(cq_s);
            //wellModel().addWellControlEq(state, well_state, aliveWells, residual_);

            if (param_.compute_well_potentials_) {
                //SolutionState state0 = state;
                //makeConstantState(state0);
                //wellModel().computeWellPotentials(mob_perfcells, b_perfcells, state0, well_state);
            }
            return iter_report;
        }

        typedef DenseAd::Evaluation<double, /*size=*/3> Eval;
        EvalWell extendEval(Eval in) const {
            EvalWell out = 0.0;
            out.value = in.value;
            for(int i = 0;i<3;++i) {
                out.derivatives[i] = in.derivatives[flowToEbosPvIdx(i)];
            }
            return out;
        }

        void print(EvalWell in) const {
            std::cout << in.value << std::endl;
            for (int i = 0; i < in.derivatives.size(); ++i) {
                std::cout << in.derivatives[i] << std::endl;
            }
        }

        void
        setWellVariables(const WellState& xw) {
            const int np = wells().number_of_phases;
            const int nw = wells().number_of_wells;
            for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                for (int w = 0; w < nw; ++w) {
                    wellVariables_[w + nw*phaseIdx] = 0.0;
                    wellVariables_[w + nw*phaseIdx].value = xw.wellSolutions()[w + nw* phaseIdx];
                    wellVariables_[w + nw*phaseIdx].derivatives[np + phaseIdx] = 1.0;
                }
            }
        }

        void
        computeAccumWells() {
            const int np = wells().number_of_phases;
            const int nw = wells().number_of_wells;
            for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                for (int w = 0; w < nw; ++w) {
                    F0_[w + nw * phaseIdx] = wellVolumeFraction(w,phaseIdx).value;
                }
            }
        }


        template <class WellState>
        void
        updatePerfPhaseRatesAndPressures(const std::vector<ADB>& cq_s,
                                         WellState& xw) const
        {
            if ( !localWellsActive() )
            {
                // If there are no wells in the subdomain of the proces then
                // cq_s has zero size and will cause a segmentation fault below.
                return;
            }

            // Update the perforation phase rates (used to calculate the pressure drop in the wellbore).
            const int np = wells().number_of_phases;
            const int nw = wells().number_of_wells;
            const int nperf = wells().well_connpos[nw];

            V cq = superset(cq_s[0].value(), Span(nperf, np, 0), nperf*np);
            for (int phase = 1; phase < np; ++phase) {
                cq += superset(cq_s[phase].value(), Span(nperf, np, phase), nperf*np);
            }
            xw.perfPhaseRates().assign(cq.data(), cq.data() + nperf*np);

            // Update the perforation pressures.
            const V& cdp = wellModel().wellPerforationPressureDiffs();
            for (int w = 0; w < nw; ++w  ) {
                for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {
                    xw.perfPress()[perf] = cdp[perf] + xw.bhp()[w];
                }
            }

        }

        void
        addWellFluxEq(std::vector<ADB> cq_s,
                      const double dt,
                      const int numBlocks)
        {
            if( !localWellsActive() )
            {
                // If there are no wells in the subdomain of the proces then
                // cq_s has zero size and will cause a segmentation fault below.
                return;
            }

            const int np = wells().number_of_phases;
            const int nw = wells().number_of_wells;

            double volume = 0.002831684659200; // 0.1 cu ft;
            //std::vector<ADB> F = wellVolumeFractions(state);
            //std::cout << F0_[0] << std::endl;
            //std::cout << F[0] << std::endl;
            //std::cout << "før Ebos" <<residual_.well_flux_eq << std::endl;
            ADB qs = ADB::constant(ADB::V::Zero(np*nw));
            for (int p = 0; p < np; ++p) {

                std::vector<EvalWell> res_vec(nw);
                for (int w = 0; w < nw; ++w) {

                    EvalWell res = (wellVolumeFraction(w, p) - F0_[w + nw*p]) * volume / dt;
                    res += getQs(w, p);
                    //for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {
                    //    res -= cq_s[perf*np + p];
                    //}
                    res_vec[w] = res;

                }
                ADB tmp = convertToADBWell(res_vec, numBlocks);
                qs += superset(tmp,Span(nw,1,p*nw), nw*np);
            }
            //std::cout << residual_.well_flux_eq << std::endl;


            //wellModel().convertToADB(res_vec, well_cells, nc, well_id, nw, numBlocks);
            //ADB qs = state.qs;
            for (int phase = 0; phase < np; ++phase) {
                qs -= superset(wellModel().wellOps().p2w * cq_s[phase], Span(nw, 1, phase*nw), nw*np);
                //qs += superset((F[phase]-F0_[phase]) * vol_dt, Span(nw,1,phase*nw), nw*np);
            }

            residual_.well_flux_eq = qs;

            //std::cout << "etter Ebos" << residual_.well_flux_eq << std::endl;

        }
        const AutoDiffBlock<double> convertToADBWell(const std::vector<EvalWell>& local, const int numVars) const
        {
            typedef typename ADB::M  M;
            const int nLocal = local.size();
            typename ADB::V value( nLocal );
            //const int numVars = 5;
            const int np = wells().number_of_phases;
            const int nw = wells().number_of_wells;

            Eigen::SparseMatrix<double> matFlux(nLocal,np*nw);

            for( int i=0; i<nLocal; ++i )
            {
                value[ i ] = local[ i ].value;
                for (int phase = 0; phase < np; ++phase) {
                    matFlux.insert(i, nw*phase + i) = local[i].derivatives[np + phase];
                }
            }

            std::vector< M > jacs( numVars );
            if (numVars == 4) {
                for( int d=0; d<np; ++d ) {
                    //Eigen::DiagonalMatrix<double>(deri[d]);
                    //jacs[ d ] = M(mat[d]);
                }

                jacs[3] = M(matFlux);
                //jacs[4] = M(matBHP);
            }
            else if (numVars == 1) {
                jacs[0] = M(matFlux);
                //jacs[1] = M(matBHP);
            }
            //std::cout << numVars << std::endl;

            return ADB::function( std::move( value ), std::move( jacs ));
        }


        void
        computeWellFluxDense(std::vector<ADB>& cq_s, const int numBlocks) const
        {
            if( ! localWellsActive() ) return ;
            const int np = wells().number_of_phases;
            const int nw = wells().number_of_wells;
            const int nperf = wells().well_connpos[nw];
            const Opm::PhaseUsage& pu = fluid_.phaseUsage();
            V Tw = Eigen::Map<const V>(wellModel().wells().WI, nperf);
            const std::vector<int>& well_cells = wellModel().wellOps().well_cells;
            std::vector<int> well_id(nperf);
            std::vector<std::vector<EvalWell>> cq_s_dense(np, std::vector<EvalWell>(nperf,0.0));


            // pressure diffs computed already (once per step, not changing per iteration)
            const V& cdp = wellModel().wellPerforationPressureDiffs();

            //std::vector<std::vector<EvalWell>> cq_s_dense(np, std::vector<EvalWell>(nperf,0.0));
            //std::vector<ADB> cmix_s_ADB = wellModel().wellVolumeFractions(state);

            for (int w = 0; w < nw; ++w) {

                EvalWell bhp = getBhp(w); //wellModel().extractDenseADWell(state.bhp,w);
//                std::cout << "well " << w << std::endl;
//                std::cout << "bhpF " << std::endl;
//                print(bhp);
//                std::cout << "state.bhp " << std::endl;
//                print(wellModel().extractDenseADWell(state.bhp,w));


                // TODO: fix for 2-phase case
                std::vector<EvalWell> cmix_s(np,0.0);
                for (int phase = 0; phase < np; ++phase) {
                    //int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                    cmix_s[phase] = wellVolumeFraction(w,phase);
                }

                //std::cout <<"cmix gas "<< w<< " "<<cmix_s[Gas] << std::endl;


                for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {
                    const int cell_idx = well_cells[perf];
                    const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                    const auto& fs = intQuants.fluidState();
                    well_id[perf] = w;
                    EvalWell pressure = extendEval(fs.pressure(FluidSystem::oilPhaseIdx)); //wellModel().extractDenseAD(state.pressure, cell_idx, cell_idx);
                    EvalWell rs = extendEval(fs.Rs()); //wellModel().extractDenseAD(state.rs, cell_idx, cell_idx);
                    EvalWell rv = extendEval(fs.Rv()); //wellModel().extractDenseAD(state.rv, cell_idx, cell_idx);
                    std::vector<EvalWell> b_perfcells_dense(np, 0.0);
                    std::vector<EvalWell> mob_perfcells_dense(np, 0.0);
                    for (int phase = 0; phase < np; ++phase) {
                        int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                        b_perfcells_dense[phase] = extendEval(fs.invB(ebosPhaseIdx));
                        mob_perfcells_dense[phase] = extendEval(intQuants.mobility(ebosPhaseIdx));
                    }

                    // Pressure drawdown (also used to determine direction of flow)
                    EvalWell well_pressure = bhp + cdp[perf];
                    EvalWell drawdown = pressure - well_pressure;

                    // injection perforations
                    if ( drawdown.value > 0 )  {

                        //Do nothing if crossflow is not allowed
                        if (!wells().allow_cf[w] && wells().type[w] == INJECTOR)
                            continue;
                        // compute phase volumetric rates at standard conditions
                        std::vector<EvalWell> cq_ps(np, 0.0);
                        for (int phase = 0; phase < np; ++phase) {
                            const EvalWell cq_p = - Tw[perf] * (mob_perfcells_dense[phase] * drawdown);
                            cq_ps[phase] = b_perfcells_dense[phase] * cq_p;
                        }


                        if ((active_)[Oil] && (active_)[Gas]) {
                            const int oilpos = pu.phase_pos[Oil];
                            const int gaspos = pu.phase_pos[Gas];
                            const EvalWell cq_psOil = cq_ps[oilpos];
                            const EvalWell cq_psGas = cq_ps[gaspos];
                            cq_ps[gaspos] += rs * cq_psOil;
                            cq_ps[oilpos] += rv * cq_psGas;
                        }


                        // map to ADB
                        for (int phase = 0; phase < np; ++phase) {
                            cq_s_dense[phase][perf] = cq_ps[phase];

                        }

                    } else {
                        //Do nothing if crossflow is not allowed
                        if (!wells().allow_cf[w] && wells().type[w] == PRODUCER)
                            continue;

                        // Using total mobilities
                        EvalWell total_mob_dense = mob_perfcells_dense[0];
                        for (int phase = 1; phase < np; ++phase) {
                            total_mob_dense += mob_perfcells_dense[phase];
                        }
                        // injection perforations total volume rates
                        const EvalWell cqt_i = - Tw[perf] * (total_mob_dense * drawdown);

                        // compute volume ratio between connection at standard conditions
                        EvalWell volumeRatio = 0.0;
                        if ((active_)[Water]) {
                            const int watpos = pu.phase_pos[Water];
                            volumeRatio += cmix_s[watpos] / b_perfcells_dense[watpos];
                        }

                        if ((active_)[Oil] && (active_)[Gas]) {
                            EvalWell well_temperature = extendEval(fs.temperature(FluidSystem::oilPhaseIdx));
                            EvalWell rsSatEval = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), well_temperature, well_pressure);
                            EvalWell rvSatEval = FluidSystem::gasPvt().saturatedOilVaporizationFactor(fs.pvtRegionIndex(), well_temperature, well_pressure);

                            const int oilpos = pu.phase_pos[Oil];
                            const int gaspos = pu.phase_pos[Gas];
                            EvalWell rvPerf = 0.0;
                            if (cmix_s[gaspos] > 0)
                                rvPerf = cmix_s[oilpos] / cmix_s[gaspos];

                            if (rvPerf.value > rvSatEval.value) {
                                rvPerf = rvSatEval;
                                //rvPerf.value = rvSatEval.value;
                            }

                            EvalWell rsPerf = 0.0;
                            if (cmix_s[oilpos] > 0)
                                rsPerf = cmix_s[gaspos] / cmix_s[oilpos];

                            if (rsPerf.value > rsSatEval.value) {
                                //rsPerf = 0.0;
                                rsPerf= rsSatEval;
                            }

                            // Incorporate RS/RV factors if both oil and gas active
                            const EvalWell d = 1.0 - rvPerf * rsPerf;

                            const EvalWell tmp_oil = (cmix_s[oilpos] - rvPerf * cmix_s[gaspos]) / d;
                            //std::cout << "tmp_oil " <<tmp_oil << std::endl;
                            volumeRatio += tmp_oil / b_perfcells_dense[oilpos];

                            const EvalWell tmp_gas = (cmix_s[gaspos] - rsPerf * cmix_s[oilpos]) / d;
                            //std::cout << "tmp_gas " <<tmp_gas << std::endl;
                            volumeRatio += tmp_gas / b_perfcells_dense[gaspos];
                        }
                        else {
                            if ((active_)[Oil]) {
                                const int oilpos = pu.phase_pos[Oil];
                                volumeRatio += cmix_s[oilpos] / b_perfcells_dense[oilpos];
                            }
                            if ((active_)[Gas]) {
                                const int gaspos = pu.phase_pos[Gas];
                                volumeRatio += cmix_s[gaspos] / b_perfcells_dense[gaspos];
                            }
                        }
                        // injecting connections total volumerates at standard conditions
                        EvalWell cqt_is = cqt_i/volumeRatio;
                        //std::cout << "volrat " << volumeRatio << " " << volrat_perf_[perf] << std::endl;
                        for (int phase = 0; phase < np; ++phase) {
                            cq_s_dense[phase][perf] = cmix_s[phase] * cqt_is; // * b_perfcells_dense[phase];
                        }
                    }
                }
            }
            const int nc = Opm::AutoDiffGrid::numCells(grid_);
            cq_s.resize(np, ADB::null());
            for (int phase = 0; phase < np; ++phase) {
                cq_s[phase] = wellModel().convertToADB(cq_s_dense[phase], well_cells, nc, well_id, nw, numBlocks);
                //std::cout << "cq_s " <<cq_s[phase] << std::endl;
            }
            //std::cout << aliveWells << std::endl;
            //std::vector<ADB> cq_s2;
            //Vector aliveWells;
            //computeWellFlux(state,mob_perfcells,b_perfcells, aliveWells,cq_s2);

            //for (int phase = 0; phase < np; ++phase) {
                //if( !(((cq_s[phase].value() - cq_s2[phase].value()).abs()<1e-10).all()) ) {
                    //std::cout << "phase " << phase << std::endl;
                    //std::cout << cq_s2[phase].value() << std::endl;
                    //std::cout << cq_s[phase].value() << std::endl;
                //}
            //}
        }
        void
        computeWellConnectionPressures(const WellState& xw)
        {
            if( ! localWellsActive() ) return ;
            // 1. Compute properties required by computeConnectionPressureDelta().
            //    Note that some of the complexity of this part is due to the function
            //    taking std::vector<double> arguments, and not Eigen objects.
            std::vector<double> b_perf;
            std::vector<double> rsmax_perf;
            std::vector<double> rvmax_perf;
            std::vector<double> surf_dens_perf;
            computePropertiesForWellConnectionPressures(xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

            const double gravity = detail::getGravity(geo_.gravity(), UgGridHelpers::dimensions(grid_));
            const V depth = Opm::AutoDiffGrid::cellCentroidsZToEigen(grid_);
            const V& pdepth = subset(depth, wellModel().wellOps().well_cells);
            const int nperf = wells().well_connpos[wells().number_of_wells];
            const std::vector<double> depth_perf(pdepth.data(), pdepth.data() + nperf);

            wellModel().computeWellConnectionDensitesPressures(xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf, depth_perf, gravity);

        }

        template<class WellState>
        void
        computePropertiesForWellConnectionPressures(const WellState& xw,
                                                    std::vector<double>& b_perf,
                                                    std::vector<double>& rsmax_perf,
                                                    std::vector<double>& rvmax_perf,
                                                    std::vector<double>& surf_dens_perf)
        {
            const int nperf = wells().well_connpos[wells().number_of_wells];
            const int nw = wells().number_of_wells;
            const std::vector<int>& well_cells = wellModel().wellOps().well_cells;
            const PhaseUsage& pu = fluid_.phaseUsage();
            const int np = fluid_.numPhases();
            b_perf.resize(nperf*np);
            rsmax_perf.resize(nperf);
            rvmax_perf.resize(nperf);

            std::vector<PhasePresence> perf_cond(nperf);
            for (int perf = 0; perf < nperf; ++perf) {
                perf_cond[perf] = phaseCondition_[well_cells[perf]];
            }

            // Compute the average pressure in each well block
            const V perf_press = Eigen::Map<const V>(xw.perfPress().data(), nperf);
            //V avg_press = perf_press*0;
            for (int w = 0; w < nw; ++w) {
                for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
                    const int cell_idx = well_cells[perf];
                    const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                    const auto& fs = intQuants.fluidState();

                    const double p_above = perf == wells().well_connpos[w] ? xw.bhp()[w] : perf_press[perf - 1];
                    const double p_avg = (perf_press[perf] + p_above)/2;
                    double temperature = fs.temperature(FluidSystem::oilPhaseIdx).value;

                    if (pu.phase_used[BlackoilPhases::Aqua]) {
                        b_perf[ pu.phase_pos[BlackoilPhases::Aqua] + perf * pu.num_phases] =
                                FluidSystem::waterPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                    }


                    if (pu.phase_used[BlackoilPhases::Vapour]) {
                        int gaspos = pu.phase_pos[BlackoilPhases::Vapour] + perf * pu.num_phases;
                        if (perf_cond[perf].hasFreeOil()) {
                            b_perf[gaspos] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                        }
                        else {
                            double rv = fs.Rv().value;
                            b_perf[gaspos] = FluidSystem::gasPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg, rv);
                        }
                    }

                    if (pu.phase_used[BlackoilPhases::Liquid]) {
                        int oilpos = pu.phase_pos[BlackoilPhases::Liquid] + perf * pu.num_phases;
                        if (perf_cond[perf].hasFreeGas()) {
                            b_perf[oilpos] = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                        }
                        else {
                            double rs = fs.Rs().value;
                            b_perf[oilpos] = FluidSystem::oilPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg, rs);
                        }
                    }

                    if (pu.phase_used[BlackoilPhases::Liquid] && pu.phase_used[BlackoilPhases::Vapour]) {
                        rsmax_perf[perf] = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), temperature, p_avg);
                        rvmax_perf[perf] = FluidSystem::gasPvt().saturatedOilVaporizationFactor(fs.pvtRegionIndex(), temperature, p_avg);
                    }

                }
            }


//            // Use cell values for the temperature as the wells don't knows its temperature yet.
//            const ADB perf_temp = subset(state.temperature, well_cells);

//            // Compute b, rsmax, rvmax values for perforations.
//            // Evaluate the properties using average well block pressures
//            // and cell values for rs, rv, phase condition and temperature.
//            const ADB avg_press_ad = ADB::constant(avg_press);
//            // const std::vector<PhasePresence>& pc = phaseCondition();

//            DataBlock b(nperf, pu.num_phases);
//            if (pu.phase_used[BlackoilPhases::Aqua]) {
//                const V bw = fluid_.bWat(avg_press_ad, perf_temp, well_cells).value();
//                b.col(pu.phase_pos[BlackoilPhases::Aqua]) = bw;
//            }
//            assert((*active_)[Oil]);
//            const V perf_so =  subset(state.saturation[pu.phase_pos[Oil]].value(), well_cells);
//            if (pu.phase_used[BlackoilPhases::Liquid]) {
//                const ADB perf_rs = (state.rs.size() > 0) ? subset(state.rs, well_cells) : ADB::null();
//                const V bo = fluid_.bOil(avg_press_ad, perf_temp, perf_rs, perf_cond, well_cells).value();
//                b.col(pu.phase_pos[BlackoilPhases::Liquid]) = bo;
//            }
//            if (pu.phase_used[BlackoilPhases::Vapour]) {
//                const ADB perf_rv = (state.rv.size() > 0) ? subset(state.rv, well_cells) : ADB::null();
//                const V bg = fluid_.bGas(avg_press_ad, perf_temp, perf_rv, perf_cond, well_cells).value();
//                b.col(pu.phase_pos[BlackoilPhases::Vapour]) = bg;
//            }
//            if (pu.phase_used[BlackoilPhases::Liquid] && pu.phase_used[BlackoilPhases::Vapour]) {
//                const V rssat = fluid_.rsSat(ADB::constant(avg_press), ADB::constant(perf_so), well_cells).value();
//                //rsmax_perf.assign(rssat.data(), rssat.data() + nperf);

//                const V rvsat = fluid_.rvSat(ADB::constant(avg_press), ADB::constant(perf_so), well_cells).value();
//                //rvmax_perf.assign(rvsat.data(), rvsat.data() + nperf);
//            }

//            // b is row major, so can just copy data.
//            //b_perf.assign(b.data(), b.data() + nperf * pu.num_phases);

            // Surface density.
            // The compute density segment wants the surface densities as
            // an np * number of wells cells array
            V rho = superset(fluid_.surfaceDensity(0 , well_cells), Span(nperf, pu.num_phases, 0), nperf*pu.num_phases);
            for (int phase = 1; phase < pu.num_phases; ++phase) {
                rho += superset(fluid_.surfaceDensity(phase , well_cells), Span(nperf, pu.num_phases, phase), nperf*pu.num_phases);
            }
            surf_dens_perf.assign(rho.data(), rho.data() + nperf * pu.num_phases);

        }




        /// \brief Compute the residual norms of the mass balance for each phase,
        /// the well flux, and the well equation.
        /// \return a vector that contains for each phase the norm of the mass balance
        /// and afterwards the norm of the residual of the well flux and the well equation.
        std::vector<double> computeResidualNorms() const
        {
            std::vector<double> residualNorms;

            std::vector<ADB>::const_iterator massBalanceIt = residual_.material_balance_eq.begin();
            const std::vector<ADB>::const_iterator endMassBalanceIt = residual_.material_balance_eq.end();

            for (; massBalanceIt != endMassBalanceIt; ++massBalanceIt) {
                const double massBalanceResid = detail::infinityNorm( (*massBalanceIt),
                                                                      linsolver_.parallelInformation() );
                if (!std::isfinite(massBalanceResid)) {
                    OPM_THROW(Opm::NumericalProblem,
                              "Encountered a non-finite residual");
                }
                residualNorms.push_back(massBalanceResid);
            }

            // the following residuals are not used in the oscillation detection now
            const double wellFluxResid = detail::infinityNormWell( residual_.well_flux_eq,
                                                                   linsolver_.parallelInformation() );
            if (!std::isfinite(wellFluxResid)) {
                OPM_THROW(Opm::NumericalProblem,
                          "Encountered a non-finite residual");
            }
            residualNorms.push_back(wellFluxResid);

            const double wellResid = detail::infinityNormWell( residual_.well_eq,
                                                               linsolver_.parallelInformation() );
            if (!std::isfinite(wellResid)) {
                OPM_THROW(Opm::NumericalProblem,
                          "Encountered a non-finite residual");
            }
            residualNorms.push_back(wellResid);

            return residualNorms;
        }


        /// \brief compute the relative change between to simulation states
        //  \return || u^n+1 - u^n || / || u^n+1 ||
        double relativeChange( const SimulationDataContainer& previous, const SimulationDataContainer& current ) const
        {
            std::vector< double > p0  ( previous.pressure() );
            std::vector< double > sat0( previous.saturation() );

            const std::size_t pSize = p0.size();
            const std::size_t satSize = sat0.size();

            // compute u^n - u^n+1
            for( std::size_t i=0; i<pSize; ++i ) {
                p0[ i ] -= current.pressure()[ i ];
            }

            for( std::size_t i=0; i<satSize; ++i ) {
                sat0[ i ] -= current.saturation()[ i ];
            }

            // compute || u^n - u^n+1 ||
            const double stateOld  = detail::euclidianNormSquared( p0.begin(),   p0.end(), 1, linsolver_.parallelInformation() ) +
                detail::euclidianNormSquared( sat0.begin(), sat0.end(),
                                              current.numPhases(),
                                              linsolver_.parallelInformation() );

            // compute || u^n+1 ||
            const double stateNew  = detail::euclidianNormSquared( current.pressure().begin(),   current.pressure().end(), 1, linsolver_.parallelInformation() ) +
                detail::euclidianNormSquared( current.saturation().begin(), current.saturation().end(),
                                              current.numPhases(),
                                              linsolver_.parallelInformation() );

            if( stateNew > 0.0 ) {
                return stateOld / stateNew ;
            }
            else {
                return 0.0;
            }
        }


        /// The size (number of unknowns) of the nonlinear system of equations.
        int sizeNonLinear() const
        {
            return residual_.sizeNonLinear();
        }

        /// Number of linear iterations used in last call to solveJacobianSystem().
        int linearIterationsLastSolve() const
        {
            return linsolver_.iterations();
        }


        /// Solve the Jacobian system Jx = r where J is the Jacobian and
        /// r is the residual.
        V solveJacobianSystem() const
        {
            return linsolver_.computeNewtonIncrement(residual_);
        }

        /// Apply an update to the primary variables, chopped if appropriate.
        /// \param[in]      dx                updates to apply to primary variables
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void updateState(const V& dx,
                         ReservoirState& reservoir_state,
                         WellState& well_state)
        {
            using namespace Opm::AutoDiffGrid;
            const int np = fluid_.numPhases();
            const int nc = numCells(grid_);
            const V null;
            assert(null.size() == 0);
            const V zero = V::Zero(nc);

            // Extract parts of dx corresponding to each part.
            const V dp = subset(dx, Span(nc));
            int varstart = nc;
            const V dsw = active_[Water] ? subset(dx, Span(nc, 1, varstart)) : null;
            varstart += dsw.size();

            const V dxvar = active_[Gas] ? subset(dx, Span(nc, 1, varstart)): null;
            varstart += dxvar.size();

            // Extract well parts np phase rates + bhp
            const V dwells = subset(dx, Span(wellModel().numWellVars(), 1, varstart));
            varstart += dwells.size();

            assert(varstart == dx.size());

            // Pressure update.
            const double dpmaxrel = dpMaxRel();
            const V p_old = Eigen::Map<const V>(&reservoir_state.pressure()[0], nc, 1);
            const V absdpmax = dpmaxrel*p_old.abs();
            const V dp_limited = sign(dp) * dp.abs().min(absdpmax);
            const V p = (p_old - dp_limited).max(zero);
            std::copy(&p[0], &p[0] + nc, reservoir_state.pressure().begin());


            // Saturation updates.
            const Opm::PhaseUsage& pu = fluid_.phaseUsage();
            const DataBlock s_old = Eigen::Map<const DataBlock>(& reservoir_state.saturation()[0], nc, np);
            const double dsmax = dsMax();

            V so;
            V sw;
            V sg;

            {
                V maxVal = zero;
                V dso = zero;
                if (active_[Water]){
                    maxVal = dsw.abs().max(maxVal);
                    dso = dso - dsw;
                }

                V dsg;
                if (active_[Gas]){
                    dsg = isSg_ * dxvar - isRv_ * dsw;
                    maxVal = dsg.abs().max(maxVal);
                    dso = dso - dsg;
                }

                maxVal = dso.abs().max(maxVal);

                V step = dsmax/maxVal;
                step = step.min(1.);

                if (active_[Water]) {
                    const int pos = pu.phase_pos[ Water ];
                    const V sw_old = s_old.col(pos);
                    sw = sw_old - step * dsw;
                }

                if (active_[Gas]) {
                    const int pos = pu.phase_pos[ Gas ];
                    const V sg_old = s_old.col(pos);
                    sg = sg_old - step * dsg;
                }

                assert(active_[Oil]);
                const int pos = pu.phase_pos[ Oil ];
                const V so_old = s_old.col(pos);
                so = so_old - step * dso;
            }

            // Appleyard chop process.
            if (active_[Gas]) {
                auto ixg = sg < 0;
                for (int c = 0; c < nc; ++c) {
                    if (ixg[c]) {
                        if (active_[Water]) {
                            sw[c] = sw[c] / (1-sg[c]);
                        }
                        so[c] = so[c] / (1-sg[c]);
                        sg[c] = 0;
                    }
                }
            }

            if (active_[Oil]) {
                auto ixo = so < 0;
                for (int c = 0; c < nc; ++c) {
                    if (ixo[c]) {
                        if (active_[Water]) {
                            sw[c] = sw[c] / (1-so[c]);
                        }
                        if (active_[Gas]) {
                            sg[c] = sg[c] / (1-so[c]);
                        }
                        so[c] = 0;
                    }
                }
            }

            if (active_[Water]) {
                auto ixw = sw < 0;
                for (int c = 0; c < nc; ++c) {
                    if (ixw[c]) {
                        so[c] = so[c] / (1-sw[c]);
                        if (active_[Gas]) {
                            sg[c] = sg[c] / (1-sw[c]);
                        }
                        sw[c] = 0;
                    }
                }
            }

            //const V sumSat = sw + so + sg;
            //sw = sw / sumSat;
            //so = so / sumSat;
            //sg = sg / sumSat;

            // Update the reservoir_state
            if (active_[Water]) {
                for (int c = 0; c < nc; ++c) {
                    reservoir_state.saturation()[c*np + pu.phase_pos[ Water ]] = sw[c];
                }
            }

            if (active_[Gas]) {
                for (int c = 0; c < nc; ++c) {
                    reservoir_state.saturation()[c*np + pu.phase_pos[ Gas ]] = sg[c];
                }
            }

            if (active_[ Oil ]) {
                for (int c = 0; c < nc; ++c) {
                    reservoir_state.saturation()[c*np + pu.phase_pos[ Oil ]] = so[c];
                }
            }

            // Update rs and rv
            const double drmaxrel = drMaxRel();
            V rs;
            if (has_disgas_) {
                const V rs_old = Eigen::Map<const V>(&reservoir_state.gasoilratio()[0], nc);
                const V drs = isRs_ * dxvar;
                const V drs_limited = sign(drs) * drs.abs().min(rs_old.abs()*drmaxrel);
                rs = rs_old - drs_limited;
            }
            V rv;
            if (has_vapoil_) {
                const V rv_old = Eigen::Map<const V>(&reservoir_state.rv()[0], nc);
                const V drv = isRv_ * dxvar;
                const V drv_limited = sign(drv) * drv.abs().min(rv_old.abs()*drmaxrel);
                rv = rv_old - drv_limited;
            }

            // Sg is used as primal variable for water only cells.
            const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
            auto watOnly = sw >  (1 - epsilon);

            // phase translation sg <-> rs
            std::vector<HydroCarbonState>& hydroCarbonState = reservoir_state.hydroCarbonState();
            std::fill(hydroCarbonState.begin(), hydroCarbonState.end(), HydroCarbonState::GasAndOil);

            if (has_disgas_) {
                const V rsSat0 = fluidRsSat(p_old, s_old.col(pu.phase_pos[Oil]), cells_);
                const V rsSat = fluidRsSat(p, so, cells_);
                // The obvious case
                auto hasGas = (sg > 0 && isRs_ == 0);

                // Set oil saturated if previous rs is sufficiently large
                const V rs_old = Eigen::Map<const V>(&reservoir_state.gasoilratio()[0], nc);
                auto gasVaporized =  ( (rs > rsSat * (1+epsilon) && isRs_ == 1 ) && (rs_old > rsSat0 * (1-epsilon)) );
                auto useSg = watOnly || hasGas || gasVaporized;
                for (int c = 0; c < nc; ++c) {
                    if (useSg[c]) {
                        rs[c] = rsSat[c];
                    } else {
                        hydroCarbonState[c] = HydroCarbonState::OilOnly;
                    }
                }

            }

            // phase transitions so <-> rv
            if (has_vapoil_) {

                // The gas pressure is needed for the rvSat calculations
                const V gaspress_old = computeGasPressure(p_old, s_old.col(Water), s_old.col(Oil), s_old.col(Gas));
                const V gaspress = computeGasPressure(p, sw, so, sg);
                const V rvSat0 = fluidRvSat(gaspress_old, s_old.col(pu.phase_pos[Oil]), cells_);
                const V rvSat = fluidRvSat(gaspress, so, cells_);

                // The obvious case
                auto hasOil = (so > 0 && isRv_ == 0);

                // Set oil saturated if previous rv is sufficiently large
                const V rv_old = Eigen::Map<const V>(&reservoir_state.rv()[0], nc);
                auto oilCondensed = ( (rv > rvSat * (1+epsilon) && isRv_ == 1) && (rv_old > rvSat0 * (1-epsilon)) );
                auto useSg = watOnly || hasOil || oilCondensed;
                for (int c = 0; c < nc; ++c) {
                    if (useSg[c]) {
                        rv[c] = rvSat[c];
                    } else {
                        hydroCarbonState[c] = HydroCarbonState::GasOnly;
                    }
                }

            }

            // Update the reservoir_state
            if (has_disgas_) {
                std::copy(&rs[0], &rs[0] + nc, reservoir_state.gasoilratio().begin());
            }

            if (has_vapoil_) {
                std::copy(&rv[0], &rv[0] + nc, reservoir_state.rv().begin());
            }


            wellModel().updateWellState(dwells, dpMaxRel(), well_state);

            // Update phase conditions used for property calculations.
            updatePhaseCondFromPrimalVariable(reservoir_state);
        }

        /// Return true if output to cout is wanted.
        bool terminalOutputEnabled() const
        {
            return terminal_output_;
        }


        /// Compute convergence based on total mass balance (tol_mb) and maximum
        /// residual mass balance (tol_cnv).
        /// \param[in]   timer       simulation timer
        /// \param[in]   dt          timestep length
        /// \param[in]   iteration   current iteration number
        bool getConvergence(const SimulatorTimerInterface& timer, const int iteration)
        {
            const double dt = timer.currentStepLength();
            const double tol_mb    = param_.tolerance_mb_;
            const double tol_cnv   = param_.tolerance_cnv_;
            const double tol_wells = param_.tolerance_wells_;

            const int nc = Opm::AutoDiffGrid::numCells(grid_);
            const int np = numPhases();
            assert(int(rq_.size()) == np);

            const V& pv = geo_.poreVolume();

            std::vector<double> R_sum(np);
            std::vector<double> B_avg(np);
            std::vector<double> maxCoeff(np);
            std::vector<double> maxNormWell(np);
            Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> B(nc, np);
            Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> R(nc, np);
            Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> tempV(nc, np);

            for ( int idx = 0; idx < np; ++idx )
            {
                V b(nc);
                for (int cell_idx = 0; cell_idx < nc; ++cell_idx) {
                    const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                    const auto& fs = intQuants.fluidState();

                    int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(idx);

                    b[cell_idx] = 1 / fs.invB(ebosPhaseIdx).value;
                }
                B.col(idx) = b;
            }


            for ( int idx = 0; idx < np; ++idx )
            {
                //const ADB& tempB = rq_[idx].b;
                //B.col(idx)       = 1./tempB.value();
                R.col(idx)       = residual_.material_balance_eq[idx].value();
                tempV.col(idx)   = R.col(idx).abs()/pv;
            }

            const double pvSum = convergenceReduction(B, tempV, R,
                                                      R_sum, maxCoeff, B_avg, maxNormWell,
                                                      nc);

            std::vector<double> CNV(np);
            std::vector<double> mass_balance_residual(np);
            std::vector<double> well_flux_residual(np);

            bool converged_MB = true;
            bool converged_CNV = true;
            bool converged_Well = true;
            // Finish computation
            for ( int idx = 0; idx < np; ++idx )
            {
                CNV[idx]                    = B_avg[idx] * dt * maxCoeff[idx];
                mass_balance_residual[idx]  = std::abs(B_avg[idx]*R_sum[idx]) * dt / pvSum;
                converged_MB                = converged_MB && (mass_balance_residual[idx] < tol_mb);
                converged_CNV               = converged_CNV && (CNV[idx] < tol_cnv);
                // Well flux convergence is only for fluid phases, not other materials
                // in our current implementation.
                assert(np >= np);
                if (idx < np) {
                    well_flux_residual[idx] = B_avg[idx] * maxNormWell[idx];
                    converged_Well = converged_Well && (well_flux_residual[idx] < tol_wells);
                }
            }

            const double residualWell     = detail::infinityNormWell(residual_.well_eq,
                                                                     linsolver_.parallelInformation());
            converged_Well = converged_Well && (residualWell < Opm::unit::barsa);
            const bool converged = converged_MB && converged_CNV && converged_Well;

            // Residual in Pascal can have high values and still be ok.
            const double maxWellResidualAllowed = 1000.0 * maxResidualAllowed();

            if ( terminal_output_ )
            {
                // Only rank 0 does print to std::cout
                if (iteration == 0) {
                    std::string msg = "Iter";
                    for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                        const std::string& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));
                        msg += "   MB(" + phaseName + ") ";
                    }
                    for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                        const std::string& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));
                        msg += "    CNV(" + phaseName + ") ";
                    }
                    for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                        const std::string& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));
                        msg += "  W-FLUX(" + phaseName + ")";
                    }
                    // std::cout << "  WELL-CONT ";
                    OpmLog::note(msg);
                }
                std::ostringstream ss;
                const std::streamsize oprec = ss.precision(3);
                const std::ios::fmtflags oflags = ss.setf(std::ios::scientific);
                ss << std::setw(4) << iteration;
                for (int idx = 0; idx < np; ++idx) {
                    ss << std::setw(11) << mass_balance_residual[idx];
                }
                for (int idx = 0; idx < np; ++idx) {
                    ss << std::setw(11) << CNV[idx];
                }
                for (int idx = 0; idx < np; ++idx) {
                    ss << std::setw(11) << well_flux_residual[idx];
                }
                // std::cout << std::setw(11) << residualWell;
                ss.precision(oprec);
                ss.flags(oflags);
                OpmLog::note(ss.str());
            }

            for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                const auto& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));

                if (std::isnan(mass_balance_residual[phaseIdx])
                    || std::isnan(CNV[phaseIdx])
                    || (phaseIdx < np && std::isnan(well_flux_residual[phaseIdx]))) {
                    OPM_THROW(Opm::NumericalProblem, "NaN residual for phase " << phaseName);
                }
                if (mass_balance_residual[phaseIdx] > maxResidualAllowed()
                    || CNV[phaseIdx] > maxResidualAllowed()
                    || (phaseIdx < np && well_flux_residual[phaseIdx] > maxResidualAllowed())) {
                    OPM_THROW(Opm::NumericalProblem, "Too large residual for phase " << phaseName);
                }
            }
            if (std::isnan(residualWell) || residualWell > maxWellResidualAllowed) {
                OPM_THROW(Opm::NumericalProblem, "NaN or too large residual for well control equation");
            }

            return converged;
        }


        /// The number of active fluid phases in the model.
        int numPhases() const
        {
            return fluid_.numPhases();
        }

        /// Update the scaling factors for mass balance equations
        void updateEquationsScaling()
        {
            ADB::V B;
            const Opm::PhaseUsage& pu = fluid_.phaseUsage();
            for ( int idx=0; idx<MaxNumPhases; ++idx )
            {
                if (active_[idx]) {
                    const int pos    = pu.phase_pos[idx];
                    const ADB& temp_b = rq_[pos].b;
                    B = 1. / temp_b.value();
#if HAVE_MPI
                    if ( linsolver_.parallelInformation().type() == typeid(ParallelISTLInformation) )
                    {
                        const ParallelISTLInformation& real_info =
                            boost::any_cast<const ParallelISTLInformation&>(linsolver_.parallelInformation());
                        double B_global_sum = 0;
                        real_info.computeReduction(B, Reduction::makeGlobalSumFunctor<double>(), B_global_sum);
                        residual_.matbalscale[idx] = B_global_sum / global_nc_;
                    }
                    else
#endif
                    {
                        residual_.matbalscale[idx] = B.mean();
                    }
                }
            }
        }


    protected:

        // ---------  Types and enums  ---------

        typedef Eigen::Array<double,
                             Eigen::Dynamic,
                             Eigen::Dynamic,
                             Eigen::RowMajor> DataBlock;

        struct ReservoirResidualQuant {
            ReservoirResidualQuant()
                : b    (   ADB::null())
                , dh   (   ADB::null())
                , mob  (   ADB::null())
            {
            }

            ADB              b;     // Reciprocal FVF
            ADB              dh;    // Pressure drop across int. interfaces
            ADB              mob;   // Phase mobility (per cell)
        };

        // ---------  Data members  ---------

        Simulator& ebosSimulator_;
        const Grid&         grid_;
        const BlackoilPropsAdInterface& fluid_;
        const DerivedGeology&           geo_;
        const RockCompressibility*      rock_comp_props_;
        VFPProperties                   vfp_properties_;
        const NewtonIterationBlackoilInterface&    linsolver_;
        // For each canonical phase -> true if active
        const std::vector<bool>         active_;
        // Size = # active phases. Maps active -> canonical phase indices.
        const std::vector<int>          canph_;
        const std::vector<int>          cells_;  // All grid cells
        HelperOps                       ops_;
        const bool has_disgas_;
        const bool has_vapoil_;

        ModelParameters                 param_;
        bool use_threshold_pressure_;
        V threshold_pressures_by_connection_;

        std::vector<ReservoirResidualQuant> rq_;
        std::vector<PhasePresence> phaseCondition_;

        // Well Model
        StandardWellsDense                       well_model_;

        V isRs_;
        V isRv_;
        V isSg_;

        LinearisedBlackoilResidual residual_;

        /// \brief Whether we print something to std::cout
        bool terminal_output_;
        /// \brief The number of cells of the global grid.
        int global_nc_;

        std::vector<std::vector<double>> residual_norms_history_;
        double current_relaxation_;
        V dx_old_;

        std::vector<EvalWell> wellVariables_;
        std::vector<double> F0_;

        // ---------  Protected methods  ---------

    public:

        EvalWell getBhp(const int wellIdx) const {
            const WellControls* wc = wells().ctrls[wellIdx];
            if (well_controls_get_current_type(wc) == BHP) {
                EvalWell bhp = 0.0;
                const double target_rate = well_controls_get_current_target(wc);
                bhp.value = target_rate;
                return bhp;
            }
            return wellVariables_[wellIdx];
        }

        EvalWell getQs(const int wellIdx, const int phaseIdx) const {
            EvalWell qs = 0.0;
            const WellControls* wc = wells().ctrls[wellIdx];
            const int np = fluid_.numPhases();
            const double target_rate = well_controls_get_current_target(wc);
//            std::cout << "well info " << std::endl;
//            std::cout << wellIdx << " " << wells().type[wellIdx] << " " <<target_rate << std::endl;
//            std::cout << phaseIdx << " " << wells().comp_frac[np*wellIdx + phaseIdx] << std::endl;
//            std::cout << well_controls_get_current_type(wc) << std::endl;
//            std::cout << "well info end" << std::endl;


            if (wells().type[wellIdx] == INJECTOR) {
                const double comp_frac = wells().comp_frac[np*wellIdx + phaseIdx];
                if (comp_frac == 0.0)
                    return qs;

                if (well_controls_get_current_type(wc) == BHP) {
                    return wellVariables_[wellIdx];
                }
                qs.value = target_rate;
                return qs;
            }

            // Producers
            if (well_controls_get_current_type(wc) == BHP) {
                return wellVariables_[wellIdx] * wellVolumeFractionScaled(wellIdx,phaseIdx);
            }
            if (well_controls_get_current_type(wc) == SURFACE_RATE) {
                const double comp_frac = wells().comp_frac[np*wellIdx + phaseIdx];

                if (comp_frac == 1.0) {
                    qs.value = target_rate;
                    return qs;
                }
                int currentControlIdx = 0;
                for (int i = 0; i < np; ++i) {
                    currentControlIdx += wells().comp_frac[np*wellIdx + i] * i;
                }

                if (wellVolumeFractionScaled(wellIdx,currentControlIdx) == 0) {
                    return qs;
                }
//                std::cout << "phase Idx " <<phaseIdx <<std::endl;
//                std::cout << "currentcontrollidx " <<currentControlIdx <<std::endl;
//                std::cout << "fraction " <<wellVolumeFractionScaled(wellIdx,phaseIdx) <<std::endl;
//                std::cout << "controll fraction " <<wellVolumeFractionScaled(wellIdx,currentControlIdx) <<std::endl;
//                std::cout << "wellvariable " <<  wellVolumeFraction(wellIdx,phaseIdx) <<std::endl;




                return (target_rate * wellVolumeFractionScaled(wellIdx,phaseIdx) / wellVolumeFractionScaled(wellIdx,currentControlIdx));
            }
            // ReservoirRate
            return target_rate * wellVolumeFractionScaled(wellIdx,phaseIdx);
        }

        EvalWell wellVolumeFraction(const int wellIdx, const int phaseIdx) const {
            assert(fluid_.numPhases() == 3);
            const int nw = wells().number_of_wells;
            if (phaseIdx == Water) {
               return wellVariables_[nw + wellIdx];
            }

            if (phaseIdx == Gas) {
               return wellVariables_[2*nw + wellIdx];
            }

            // Oil
            return 1.0 - wellVariables_[nw + wellIdx] - wellVariables_[2 * nw + wellIdx];
        }

        EvalWell wellVolumeFractionScaled(const int wellIdx, const int phaseIdx) const {
            const WellControls* wc = wells().ctrls[wellIdx];
            if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {
                const double* distr = well_controls_get_current_distr(wc);
                return wellVolumeFraction(wellIdx, phaseIdx) / distr[phaseIdx];
            }
            std::vector<double> g = {1,1,0.01};
            return (wellVolumeFraction(wellIdx, phaseIdx) / g[phaseIdx]);
        }




        /// return the StandardWells object
        StandardWellsDense& wellModel() { return well_model_; }
        const StandardWellsDense& wellModel() const { return well_model_; }

        /// return the Well struct in the StandardWells
        const Wells& wells() const { return well_model_.wells(); }

        /// return true if wells are available in the reservoir
        bool wellsActive() const { return well_model_.wellsActive(); }

        /// return true if wells are available on this process
        bool localWellsActive() const { return well_model_.localWellsActive(); }

        void
        makeConstantState(SolutionState& state) const
        {
            // HACK: throw away the derivatives. this may not be the most
            // performant way to do things, but it will make the state
            // automatically consistent with variableState() (and doing
            // things automatically is all the rage in this module ;)
            state.pressure = ADB::constant(state.pressure.value());
            state.temperature = ADB::constant(state.temperature.value());
            state.rs = ADB::constant(state.rs.value());
            state.rv = ADB::constant(state.rv.value());
            const int num_phases = state.saturation.size();
            for (int phaseIdx = 0; phaseIdx < num_phases; ++ phaseIdx) {
                state.saturation[phaseIdx] = ADB::constant(state.saturation[phaseIdx].value());
            }
            state.qs = ADB::constant(state.qs.value());
            state.bhp = ADB::constant(state.bhp.value());
            assert(state.canonical_phase_pressures.size() == static_cast<std::size_t>(Opm::BlackoilPhases::MaxNumPhases));
            for (int canphase = 0; canphase < Opm::BlackoilPhases::MaxNumPhases; ++canphase) {
                ADB& pp = state.canonical_phase_pressures[canphase];
                pp = ADB::constant(pp.value());
            }
            state.wellVariables = ADB::constant(state.wellVariables.value());
        }

        void setupLegacyState(SolutionState& state,
                              const ReservoirState& x,
                              const WellState& xw) const
        {
            const int nc = Opm::AutoDiffGrid::numCells(grid_);
            const int np = x.numPhases();

            std::vector<V> vars0(np, V::Ones(nc, 1));
            wellModel().variableWellStateInitials(xw, vars0);
            std::vector<ADB> vars = ADB::variables(vars0);

            std::vector<int> indices = {{Pressure, Sw, Xvar}};
            int foo = indices.size();
            indices.resize(4);
            wellModel().variableStateWellIndices(indices, foo);

            const ADB& ones = ADB::constant(V::Ones(nc, 1));


            // temperature cannot be a variable at this time (only constant).
            state.temperature = ones;

            // saturations
            state.saturation[Water] = std::move(vars[indices[Sw]]);

            const ADB& sw = state.saturation[Water];
            const ADB& xvar = vars[indices[Xvar]];
            state.saturation[Gas] = xvar;
            state.saturation[Oil] = sw + xvar;

            // pressures
            state.pressure = std::move(vars[indices[Pressure]]);
            const ADB& po = state.pressure;

            const ADB& tmp = po + sw + xvar;
            state.canonical_phase_pressures[Gas] = tmp;
            state.canonical_phase_pressures[Water] = tmp;
            state.canonical_phase_pressures[Oil] = po;

            if (has_disgas_) {
                state.rs = po + xvar;
            } else {
                state.rs = po;
            }
            if (has_vapoil_) {
                state.rv = po + xvar;
            } else {
                state.rv = po;
            }

            // Note that so is never a primary variable.
            //state.saturation[Oil] = std::move(so);

            // wells
            wellModel().variableStateExtractWellsVars(indices, vars, state, xw);
        }

        V
        fluidRsSat(const V&                p,
                   const V&                satOil,
                   const std::vector<int>& cells) const
        {
            return fluid_.rsSat(ADB::constant(p), ADB::constant(satOil), cells).value();
        }

        V
        fluidRvSat(const V&                p,
                   const V&                satOil,
                   const std::vector<int>& cells) const
        {
            return fluid_.rvSat(ADB::constant(p), ADB::constant(satOil), cells).value();
        }

        void convertInput( const int iterationIdx,
                           const ReservoirState& reservoirState,
                           Simulator& simulator ) const
        {
            SolutionVector& solution = simulator.model().solution( 0 /* timeIdx */ );
            const Opm::PhaseUsage pu = fluid_.phaseUsage();

            const int numCells = reservoirState.numCells();
            const int numPhases = fluid_.numPhases();
            const auto& oilPressure = reservoirState.pressure();
            const auto& saturations = reservoirState.saturation();
            const auto& rs          = reservoirState.gasoilratio();
            const auto& rv          = reservoirState.rv();
            for( int cellIdx = 0; cellIdx<numCells; ++cellIdx )
            {
                // set non-switching primary variables
                PrimaryVariables& cellPv = solution[ cellIdx ];
                // set water saturation
                cellPv[BlackoilIndices::waterSaturationIdx] = saturations[cellIdx*numPhases + pu.phase_pos[Water]];

                // set switching variable and interpretation
                if( isRs_[ cellIdx ] && has_disgas_ )
                {
                    cellPv[BlackoilIndices::compositionSwitchIdx] = rs[cellIdx];
                    cellPv[BlackoilIndices::pressureSwitchIdx] = oilPressure[cellIdx];
                    cellPv.setPrimaryVarsMeaning( PrimaryVariables::Sw_po_Rs );
                }
                else if( isRv_[ cellIdx ] && has_vapoil_ )
                {
                    // this case (-> gas only with vaporized oil in the gas) is
                    // relatively expensive as it requires to compute the capillary
                    // pressure in order to get the gas phase pressure. (the reason why
                    // ebos uses the gas pressure here is that it makes the common case
                    // of the primary variable switching code fast because to determine
                    // whether the oil phase appears one needs to compute the Rv value
                    // for the saturated gas phase and if this is not available as a
                    // primary variable, it needs to be computed.) luckily for here, the
                    // gas-only case is not too common, so the performance impact of this
                    // is limited.
                    typedef Opm::SimpleModularFluidState<double,
                                                         /*numPhases=*/3,
                                                         /*numComponents=*/3,
                                                         FluidSystem,
                                                         /*storePressure=*/false,
                                                         /*storeTemperature=*/false,
                                                         /*storeComposition=*/false,
                                                         /*storeFugacity=*/false,
                                                         /*storeSaturation=*/true,
                                                         /*storeDensity=*/false,
                                                         /*storeViscosity=*/false,
                                                         /*storeEnthalpy=*/false> SatOnlyFluidState;
                    SatOnlyFluidState fluidState;
                    fluidState.setSaturation(FluidSystem::waterPhaseIdx, saturations[cellIdx*numPhases + pu.phase_pos[Water]]);
                    fluidState.setSaturation(FluidSystem::oilPhaseIdx, saturations[cellIdx*numPhases + pu.phase_pos[Oil]]);
                    fluidState.setSaturation(FluidSystem::gasPhaseIdx, saturations[cellIdx*numPhases + pu.phase_pos[Gas]]);

                    double pC[/*numPhases=*/3] = { 0.0, 0.0, 0.0 };
                    const MaterialLawParams& matParams = simulator.problem().materialLawParams(cellIdx);
                    MaterialLaw::capillaryPressures(pC, matParams, fluidState);
                    double pg = oilPressure[cellIdx] + (pC[FluidSystem::gasPhaseIdx] - pC[FluidSystem::oilPhaseIdx]);

                    cellPv[BlackoilIndices::compositionSwitchIdx] = rv[cellIdx];
                    cellPv[BlackoilIndices::pressureSwitchIdx] = pg;
                    cellPv.setPrimaryVarsMeaning( PrimaryVariables::Sw_pg_Rv );
                }
                else
                {
                    assert(isSg_[cellIdx]);
                    cellPv[BlackoilIndices::compositionSwitchIdx] = saturations[cellIdx*numPhases + pu.phase_pos[Gas]];
                    cellPv[BlackoilIndices::pressureSwitchIdx] = oilPressure[ cellIdx ];
                    cellPv.setPrimaryVarsMeaning( PrimaryVariables::Sw_po_Sg );
                }
            }

            if( iterationIdx == 0 )
            {
                simulator.model().solution( 1 /* timeIdx */ ) = solution;
            }
        }

    public:
        /*
        int ebosCompToFlowPhaseIdx( const int compIdx ) const
        {
            const int compToPhase[ 3 ] = { Oil, Water, Gas };
            return compToPhase[ compIdx ];
        }
        */

        int flowToEbosPvIdx( const int flowPv ) const
        {
            const int flowToEbos[ 3 ] = {
                                          BlackoilIndices::pressureSwitchIdx,
                                          BlackoilIndices::waterSaturationIdx,
                                          BlackoilIndices::compositionSwitchIdx
                                        };
            return flowToEbos[ flowPv ];
        }

        int flowPhaseToEbosCompIdx( const int phaseIdx ) const
        {
            const int phaseToComp[ 3 ] = { FluidSystem::waterCompIdx, FluidSystem::oilCompIdx, FluidSystem::gasCompIdx };
            return phaseToComp[ phaseIdx ];
        }

        int flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
        {
            const int flowToEbos[ 3 ] = { FluidSystem::waterPhaseIdx, FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx };
            return flowToEbos[ phaseIdx ];
        }

    private:
        void convertResults(const Simulator& simulator)
        {
            const auto& ebosJac = simulator.model().linearizer().matrix();
            const auto& ebosResid = simulator.model().linearizer().residual();

            const int numPhases = wells().number_of_phases;
            const int numCells = ebosJac.N();
            const int cols = ebosJac.M();
            assert( numCells == cols );

            // create the matrices and the right hand sides in a format which is more
            // appropriate for the conversion from what's produced ebos to the flow stuff
            typedef Eigen::SparseMatrix<double, Eigen::RowMajor> M;
            typedef ADB::V  V;
            std::vector< std::vector< M > > jacs( numPhases );
            std::vector< V > resid (numPhases);
            for( int eqIdx = 0; eqIdx < numPhases; ++eqIdx )
            {
                jacs[ eqIdx ].resize( numPhases );
                resid[ eqIdx ].resize( numCells );
                for( int pvIdx = 0; pvIdx < numPhases; ++pvIdx )
                {
                    jacs[ eqIdx ][ pvIdx ] = M( numCells, cols );
                    jacs[ eqIdx ][ pvIdx ].reserve( ebosJac.nonzeroes() );
                }
            }

            // write the right-hand-side values from the ebosJac into the objects
            // allocated above.
            const auto endrow = ebosJac.end();
            for( int cellIdx = 0; cellIdx < numCells; ++cellIdx )
            {
                const double cellVolume = simulator.model().dofTotalVolume(cellIdx);
                const auto& cellRes = ebosResid[ cellIdx ];

                for( int flowPhaseIdx = 0; flowPhaseIdx < numPhases; ++flowPhaseIdx )
                {
                    const double refDens = FluidSystem::referenceDensity( flowPhaseToEbosPhaseIdx( flowPhaseIdx ), 0 );
                    double ebosVal = cellRes[ flowPhaseToEbosCompIdx( flowPhaseIdx ) ] / refDens * cellVolume;

                    resid[ flowPhaseIdx ][ cellIdx ] = ebosVal;
                }
            }

            for( auto row = ebosJac.begin(); row != endrow; ++row )
            {
                const int rowIdx = row.index();
                const double cellVolume = simulator.model().dofTotalVolume(rowIdx);

                for( int flowPhaseIdx = 0; flowPhaseIdx < numPhases; ++flowPhaseIdx )
                {
                    for( int pvIdx = 0; pvIdx < numPhases; ++pvIdx )
                    {
                        jacs[flowPhaseIdx][pvIdx].startVec(rowIdx);
                    }
                }

                // translate the Jacobian of the residual from the format used by ebos to
                // the one expected by flow
                const auto endcol = row->end();
                for( auto col = row->begin(); col != endcol; ++col )
                {
                    const int colIdx = col.index();
                    for( int flowPhaseIdx = 0; flowPhaseIdx < numPhases; ++flowPhaseIdx )
                    {
                        const double refDens = FluidSystem::referenceDensity( flowPhaseToEbosPhaseIdx( flowPhaseIdx ), 0 );
                        for( int pvIdx=0; pvIdx<numPhases; ++pvIdx )
                        {
                            double ebosVal = (*col)[flowPhaseToEbosCompIdx(flowPhaseIdx)][flowToEbosPvIdx(pvIdx)]/refDens*cellVolume;
                            if (ebosVal != 0.0)
                                jacs[flowPhaseIdx][pvIdx].insertBackByOuterInnerUnordered(rowIdx, colIdx) = ebosVal;
                        }
                    }
                }
            }

            // convert the resulting matrices from/ row major ordering to colum major.
            typedef typename ADB::M ADBJacobianMatrix;
            std::vector< std::vector< ADBJacobianMatrix > > adbJacs( numPhases );

            for( int flowPhaseIdx = 0; flowPhaseIdx < numPhases; ++flowPhaseIdx )
            {
                adbJacs[ flowPhaseIdx ].resize( numPhases + 1 );
                for( int pvIdx = 0; pvIdx < numPhases; ++pvIdx )
                {
                    jacs[ flowPhaseIdx ][ pvIdx ].finalize();
                    adbJacs[ flowPhaseIdx ][ pvIdx ].assign( std::move(jacs[ flowPhaseIdx ][ pvIdx ]) );
                }
                // add two "dummy" matrices for the well primary variables
                const int numWells = wells().number_of_wells;
                for( int pvIdx = numPhases; pvIdx < numPhases + 1; ++pvIdx ) {
                    adbJacs[ flowPhaseIdx ][ pvIdx ] = ADBJacobianMatrix( numPhases*numWells, numPhases*numWells );
                        //sparsityPattern.derivative()[pvIdx];
                }
            }

            for( int eqIdx = 0; eqIdx < numPhases; ++eqIdx )
            {
                residual_.material_balance_eq[ eqIdx ] =
                    ADB::function(std::move(resid[eqIdx]),
                                  std::move(adbJacs[eqIdx]));
            }
        }

        void updateLegacyState(const Simulator& simulator, SolutionState& legacyState)
        {
            const int numPhases = 3;
            const int numCells = simulator.model().numGridDof();

            typedef Eigen::SparseMatrix<double, Eigen::ColMajor> EigenMatrix;

            ///////
            // create the value vectors for the legacy state
            ///////
            V poVal;
            V TVal;
            std::vector<V> SVal(numPhases);
            std::vector<V> mobVal(numPhases);
            std::vector<V> bVal(numPhases);
            std::vector<V> pVal(numPhases);
            V RsVal;
            V RvVal;

            poVal.resize(numCells);
            TVal.resize(numCells);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                SVal[phaseIdx].resize(numCells);
                mobVal[phaseIdx].resize(numCells);
                bVal[phaseIdx].resize(numCells);
                pVal[phaseIdx].resize(numCells);
            }
            RsVal.resize(numCells);
            RvVal.resize(numCells);

            ///////
            // create the Jacobian matrices for the legacy state. here we assume that the
            // sparsity pattern of the inputs is already correct
            ///////
            std::vector<EigenMatrix> poJac(numPhases + 1);
            //std::vector<EigenMatrix> TJac(numPhases + 2);
            std::vector<std::vector<EigenMatrix>> SJac(numPhases);
            std::vector<std::vector<EigenMatrix>> mobJac(numPhases);
            std::vector<std::vector<EigenMatrix>> bJac(numPhases);
            std::vector<std::vector<EigenMatrix>> pJac(numPhases);
            std::vector<EigenMatrix> RsJac(numPhases + 1);
            std::vector<EigenMatrix> RvJac(numPhases + 1);

            // reservoir stuff
            for (int pvIdx = 0; pvIdx < numPhases; ++ pvIdx) {
                poJac[pvIdx].resize(numCells, numCells);
                //TJac[pvIdx].resize(numCells, numCells);
                RsJac[pvIdx].resize(numCells, numCells);
                RvJac[pvIdx].resize(numCells, numCells);

                poJac[pvIdx].reserve(numCells);
                //TJac[pvIdx].reserve(numCells);
                RsJac[pvIdx].reserve(numCells);
                RvJac[pvIdx].reserve(numCells);
            }

            // auxiliary equations
            for (int pvIdx = numPhases; pvIdx < numPhases + 1; ++ pvIdx) {
                legacyState.pressure.derivative()[pvIdx].toSparse(poJac[pvIdx]);
                //legacyState.temperature.derivative()[pvIdx].toSparse(TJac[pvIdx]);
                legacyState.rs.derivative()[pvIdx].toSparse(RsJac[pvIdx]);
                legacyState.rv.derivative()[pvIdx].toSparse(RvJac[pvIdx]);
            }

            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                SJac[phaseIdx].resize(numPhases + 1);
                mobJac[phaseIdx].resize(numPhases + 1);
                bJac[phaseIdx].resize(numPhases + 1);
                pJac[phaseIdx].resize(numPhases + 1);
                for (int pvIdx = 0; pvIdx < numPhases; ++ pvIdx) {
                    SJac[phaseIdx][pvIdx].resize(numCells, numCells);
                    SJac[phaseIdx][pvIdx].reserve(numCells);

                    mobJac[phaseIdx][pvIdx].resize(numCells, numCells);
                    mobJac[phaseIdx][pvIdx].reserve(numCells);

                    bJac[phaseIdx][pvIdx].resize(numCells, numCells);
                    bJac[phaseIdx][pvIdx].reserve(numCells);

                    pJac[phaseIdx][pvIdx].resize(numCells, numCells);
                    pJac[phaseIdx][pvIdx].reserve(numCells);
                }

                // auxiliary equations for the saturations and pressures
                for (int pvIdx = numPhases; pvIdx < numPhases + 1; ++ pvIdx) {
                    legacyState.saturation[phaseIdx].derivative()[pvIdx].toSparse(SJac[phaseIdx][pvIdx]);
                    legacyState.saturation[phaseIdx].derivative()[pvIdx].toSparse(mobJac[phaseIdx][pvIdx]);
                    legacyState.saturation[phaseIdx].derivative()[pvIdx].toSparse(bJac[phaseIdx][pvIdx]);
                    legacyState.canonical_phase_pressures[phaseIdx].derivative()[pvIdx].toSparse(pJac[phaseIdx][pvIdx]);
                }
            }

            ///////
            // write the values and the derivatives into the data structures for the
            // legacy state.
            ///////
            for( int cellIdx = 0; cellIdx < numCells; ++cellIdx )
            {
                const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cellIdx, /*timeIdx=*/0));
                const auto& fs = intQuants.fluidState();

                poVal[cellIdx] = fs.pressure(FluidSystem::oilPhaseIdx).value;
                TVal[cellIdx] = fs.temperature(0).value;
                RsVal[cellIdx] = fs.Rs().value;
                RvVal[cellIdx] = fs.Rv().value;

                for (int pvIdx = 0; pvIdx < numPhases; ++pvIdx) {
                    poJac[pvIdx].startVec(cellIdx);
                    //TJac[pvIdx].startVec(cellIdx);
                    RsJac[pvIdx].startVec(cellIdx);
                    RvJac[pvIdx].startVec(cellIdx);

                    poJac[pvIdx].insertBackByOuterInnerUnordered(cellIdx, cellIdx) = fs.pressure(FluidSystem::oilPhaseIdx).derivatives[flowToEbosPvIdx(pvIdx)];
                    //TJac[pvIdx].insertBackByOuterInnerUnordered(cellIdx, cellIdx) = fs.temperature(FluidSystem::oilPhaseIdx).derivatives[flowToEbosPvIdx(pvIdx)];
                    RsJac[pvIdx].insertBackByOuterInnerUnordered(cellIdx, cellIdx) = fs.Rs().derivatives[flowToEbosPvIdx(pvIdx)];
                    RvJac[pvIdx].insertBackByOuterInnerUnordered(cellIdx, cellIdx) = fs.Rv().derivatives[flowToEbosPvIdx(pvIdx)];
                }

                for( int flowPhaseIdx = 0; flowPhaseIdx < numPhases; ++flowPhaseIdx )
                {
                    int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(flowPhaseIdx);
                    SVal[flowPhaseIdx][cellIdx] = fs.saturation(ebosPhaseIdx).value;
                    mobVal[flowPhaseIdx][cellIdx] = intQuants.mobility(ebosPhaseIdx).value;
                    bVal[flowPhaseIdx][cellIdx] = fs.invB(ebosPhaseIdx).value;
                    pVal[flowPhaseIdx][cellIdx] = fs.pressure(ebosPhaseIdx).value;

                    for (int pvIdx = 0; pvIdx < numPhases; ++pvIdx) {
                        SJac[flowPhaseIdx][pvIdx].startVec(cellIdx);
                        mobJac[flowPhaseIdx][pvIdx].startVec(cellIdx);
                        bJac[flowPhaseIdx][pvIdx].startVec(cellIdx);
                        pJac[flowPhaseIdx][pvIdx].startVec(cellIdx);
                        SJac[flowPhaseIdx][pvIdx].insertBackByOuterInnerUnordered(cellIdx, cellIdx) = fs.saturation(ebosPhaseIdx).derivatives[flowToEbosPvIdx(pvIdx)];
                        mobJac[flowPhaseIdx][pvIdx].insertBackByOuterInnerUnordered(cellIdx, cellIdx) = intQuants.mobility(ebosPhaseIdx).derivatives[flowToEbosPvIdx(pvIdx)];
                        bJac[flowPhaseIdx][pvIdx].insertBackByOuterInnerUnordered(cellIdx, cellIdx) = fs.invB(ebosPhaseIdx).derivatives[flowToEbosPvIdx(pvIdx)];
                        pJac[flowPhaseIdx][pvIdx].insertBackByOuterInnerUnordered(cellIdx, cellIdx) = fs.pressure(ebosPhaseIdx).derivatives[flowToEbosPvIdx(pvIdx)];
                    }
                }
            }

            // finalize all Jacobian matrices
            for (int pvIdx = 0; pvIdx < numPhases; ++pvIdx) {
                poJac[pvIdx].finalize();
                //TJac[pvIdx].finalize();
                RsJac[pvIdx].finalize();
                RvJac[pvIdx].finalize();

                for (int phaseIdx = 0; phaseIdx < 3; ++ phaseIdx) {
                    SJac[phaseIdx][pvIdx].finalize();
                    mobJac[phaseIdx][pvIdx].finalize();
                    bJac[phaseIdx][pvIdx].finalize();
                    pJac[phaseIdx][pvIdx].finalize();
                }
            }

            ///////
            // create Opm::AutoDiffMatrix objects from Eigen::SparseMatrix
            // objects. (Opm::AutoDiffMatrix is not directly assignable, wtf?)
            ///////
            typedef typename ADB::M ADBJacobianMatrix;

            std::vector<ADBJacobianMatrix> poAdbJacs;
            std::vector<ADBJacobianMatrix> RsAdbJacs;
            std::vector<ADBJacobianMatrix> RvAdbJacs;

            poAdbJacs.resize(numPhases + 1);
            RsAdbJacs.resize(numPhases + 1);
            RvAdbJacs.resize(numPhases + 1);
            for(int pvIdx = 0; pvIdx < numPhases + 1; ++pvIdx)
            {
                poAdbJacs[pvIdx].assign(poJac[pvIdx]);
                RsAdbJacs[pvIdx].assign(RsJac[pvIdx]);
                RvAdbJacs[pvIdx].assign(RvJac[pvIdx]);
            }

            std::vector<std::vector<ADBJacobianMatrix>> SAdbJacs(numPhases);
            std::vector<std::vector<ADBJacobianMatrix>> mobAdbJacs(numPhases);
            std::vector<std::vector<ADBJacobianMatrix>> bAdbJacs(numPhases);
            std::vector<std::vector<ADBJacobianMatrix>> pAdbJacs(numPhases);
            for(int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                SAdbJacs[phaseIdx].resize(numPhases + 1);
                mobAdbJacs[phaseIdx].resize(numPhases + 1);
                bAdbJacs[phaseIdx].resize(numPhases + 1);
                pAdbJacs[phaseIdx].resize(numPhases + 1);
                for(int pvIdx = 0; pvIdx < numPhases + 1; ++pvIdx)
                {
                    SAdbJacs[phaseIdx][pvIdx].assign(SJac[phaseIdx][pvIdx]);
                    mobAdbJacs[phaseIdx][pvIdx].assign(mobJac[phaseIdx][pvIdx]);
                    bAdbJacs[phaseIdx][pvIdx].assign(bJac[phaseIdx][pvIdx]);
                    pAdbJacs[phaseIdx][pvIdx].assign(pJac[phaseIdx][pvIdx]);
                }
            }

            ///////
            // create the ADB objects in the legacy state
            ///////
            legacyState.pressure =
                ADB::function(std::move(poVal),
                              std::move(poAdbJacs));
            legacyState.temperature =
                ADB::constant(std::move(TVal));
            legacyState.rs =
                ADB::function(std::move(RsVal),
                              std::move(RsAdbJacs));
            legacyState.rv =
                ADB::function(std::move(RvVal),
                              std::move(RvAdbJacs));

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            {
                legacyState.saturation[phaseIdx] =
                    ADB::function(std::move(SVal[phaseIdx]),
                                  std::move(SAdbJacs[phaseIdx]));
                legacyState.canonical_phase_pressures[phaseIdx] =
                    ADB::function(std::move(pVal[phaseIdx]),
                                  std::move(pAdbJacs[phaseIdx]));
                rq_[phaseIdx].mob =
                    ADB::function(std::move(mobVal[phaseIdx]),
                                 std::move(mobAdbJacs[phaseIdx]));
                rq_[phaseIdx].b =
                    ADB::function(std::move(bVal[phaseIdx]),
                                  std::move(bAdbJacs[phaseIdx]));
            }
        }

    public:
        void beginReportStep()
        {
            isBeginReportStep_ = true;
        }

        void endReportStep()
        {
            ebosSimulator_.problem().endEpisode();
        }

    private:
        void assembleMassBalanceEq(const SimulatorTimerInterface& timer,
                                   const int iterationIdx,
                                   const ReservoirState& reservoirState)
        {
            convertInput( iterationIdx, reservoirState, ebosSimulator_ );

            ebosSimulator_.startNextEpisode( timer.currentStepLength() );
            ebosSimulator_.setEpisodeIndex( timer.reportStepNum() );
            ebosSimulator_.setTimeStepIndex( timer.reportStepNum() );
            ebosSimulator_.model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);
            ebosSimulator_.model().newtonMethod().setIterationIndex(iterationIdx);

            static int prevEpisodeIdx = 10000;

            // notify ebos about the end of the previous episode and time step if applicable
            if (isBeginReportStep_) {
                isBeginReportStep_ = false;
                ebosSimulator_.problem().beginEpisode();
            }

            // doing the notifactions here is conceptually wrong and also causes the
            // endTimeStep() and endEpisode() methods to be not called for the
            // simulation's last time step and episode.
            if (ebosSimulator_.model().newtonMethod().numIterations() == 0
                && prevEpisodeIdx < timer.reportStepNum())
            {
                ebosSimulator_.problem().endTimeStep();
            }

            ebosSimulator_.setTimeStepSize( timer.currentStepLength() );
            if (ebosSimulator_.model().newtonMethod().numIterations() == 0)
            {
                ebosSimulator_.problem().beginTimeStep();
            }

            ebosSimulator_.problem().beginIteration();
            ebosSimulator_.model().linearizer().linearize();
            ebosSimulator_.problem().endIteration();

            prevEpisodeIdx = ebosSimulator_.episodeIndex();

            convertResults(ebosSimulator_);
            //updateLegacyState(ebosSimulator_, state);

            if (param_.update_equations_scaling_) {
                updateEquationsScaling();
            }

        }

        IterationReport solveWellEq(const std::vector<ADB>& mob_perfcells,
                                    const std::vector<ADB>& b_perfcells,
                                    const double dt,
                                    WellState& well_state)
        {
            const int np = wells().number_of_phases;
            //std::vector<ADB> cq_s(np, ADB::null());
            const int nw = wells().number_of_wells;
            const int nperf = wells().well_connpos[nw];
            std::vector<ADB> cq_s(np, ADB::null());
            WellState well_state0 = well_state;

            int it  = 0;
            bool converged;
            do {
                // bhp and Q for the wells
                computeWellFluxDense(cq_s, 1);
                updatePerfPhaseRatesAndPressures(cq_s, well_state);
                addWellFluxEq(cq_s, dt, 1);
                converged = getWellConvergence(it);

                if (converged) {
                    break;
                }

                ++it;
                if( localWellsActive() )
                {
                    std::vector<ADB> eqs;
                    eqs.reserve(1);
                    eqs.push_back(residual_.well_flux_eq);
                    //eqs.push_back(residual_.well_eq);
                    ADB total_residual = vertcatCollapseJacs(eqs);
                    const std::vector<M>& Jn = total_residual.derivative();
                    typedef Eigen::SparseMatrix<double> Sp;
                    Sp Jn0;
                    Jn[0].toSparse(Jn0);
                    const Eigen::SparseLU< Sp > solver(Jn0);
                    ADB::V total_residual_v = total_residual.value();
                    const Eigen::VectorXd& dx = solver.solve(total_residual_v.matrix());
                    assert(dx.size() == total_residual_v.size());
                    wellModel().updateWellState(dx.array(), dpMaxRel(), well_state);
                    wellModel().updateWellControls(well_state);
                    setWellVariables(well_state);
                }
            } while (it < 15);

            if (!converged) {
                well_state = well_state0;
            }

            const bool failed = false; // Not needed in this method.
            const int linear_iters = 0; // Not needed in this method
            return IterationReport{failed, converged, linear_iters, it};
        }


        void
        addWellContributionToMassBalanceEq(const std::vector<ADB>& cq_s)
        {
            if ( !localWellsActive() )
            {
                // If there are no wells in the subdomain of the proces then
                // cq_s has zero size and will cause a segmentation fault below.
                return;
            }

            // Add well contributions to mass balance equations
            const int nc = Opm::AutoDiffGrid::numCells(grid_);
            const int np = numPhases();
            for (int phase = 0; phase < np; ++phase) {
                residual_.material_balance_eq[phase] -= superset(cq_s[phase], wellModel().wellOps().well_cells, nc);
                }
        }


        bool getWellConvergence(const int iteration)
        {
            const double tol_wells = param_.tolerance_wells_;

            const int nc = Opm::AutoDiffGrid::numCells(grid_);
            const int np = numPhases();

            const V& pv = geo_.poreVolume();
            std::vector<double> R_sum(np);
            std::vector<double> B_avg(np);
            std::vector<double> maxCoeff(np);
            std::vector<double> maxNormWell(np);
            Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> B(nc, np);
            Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> R(nc, np);
            Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> tempV(nc, np);

            for ( int idx = 0; idx < np; ++idx )
            {
                V b(nc);
                for (int cell_idx = 0; cell_idx < nc; ++cell_idx) {
                    const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                    const auto& fs = intQuants.fluidState();

                    int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(idx);

                    b[cell_idx] = 1 / fs.invB(ebosPhaseIdx).value;
                }
                B.col(idx) = b;
            }

            for ( int idx = 0; idx < np; ++idx )
            {
                //const ADB& tempB = rq_[idx].b;
                //B.col(idx)       = 1./tempB.value();
                R.col(idx)       = residual_.material_balance_eq[idx].value();
                tempV.col(idx)   = R.col(idx).abs()/pv;
            }

            convergenceReduction(B, tempV, R, R_sum, maxCoeff, B_avg, maxNormWell, nc);

            std::vector<double> well_flux_residual(np);
            bool converged_Well = true;
            // Finish computation
            for ( int idx = 0; idx < np; ++idx )
            {
                well_flux_residual[idx] = B_avg[idx] * maxNormWell[idx];
                converged_Well = converged_Well && (well_flux_residual[idx] < tol_wells);
            }

            const double residualWell     = detail::infinityNormWell(residual_.well_eq,
                                                                     linsolver_.parallelInformation());
            converged_Well = converged_Well && (residualWell < Opm::unit::barsa);
            const bool converged = converged_Well;

            // if one of the residuals is NaN, throw exception, so that the solver can be restarted
            for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                const auto& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));

                if (std::isnan(well_flux_residual[phaseIdx])) {
                    OPM_THROW(Opm::NumericalProblem, "NaN residual for phase " << phaseName);
                }
                if (well_flux_residual[phaseIdx] > maxResidualAllowed()) {
                    OPM_THROW(Opm::NumericalProblem, "Too large residual for phase " << phaseName);
                }
            }

            if ( terminal_output_ )
            {
                // Only rank 0 does print to std::cout
                if (iteration == 0) {
                    std::string msg;
                    msg = "Iter";
                    for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                        const std::string& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));
                        msg += "  W-FLUX(" + phaseName + ")";
                    }
                    OpmLog::note(msg);
                }
                std::ostringstream ss;
                const std::streamsize oprec = ss.precision(3);
                const std::ios::fmtflags oflags = ss.setf(std::ios::scientific);
                ss << std::setw(4) << iteration;
                for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                    ss << std::setw(11) << well_flux_residual[phaseIdx];
                }
                ss.precision(oprec);
                ss.flags(oflags);
                OpmLog::note(ss.str());
            }
            return converged;
        }

        std::vector<ADB>
        computePressures(const ADB& po,
                         const ADB& sw,
                         const ADB& so,
                         const ADB& sg) const
        {
            // convert the pressure offsets to the capillary pressures
            std::vector<ADB> pressure = fluid_.capPress(sw, so, sg, cells_);
            for (int phaseIdx = 0; phaseIdx < BlackoilPhases::MaxNumPhases; ++phaseIdx) {
                // The reference pressure is always the liquid phase (oil) pressure.
                if (phaseIdx == BlackoilPhases::Liquid)
                    continue;
                if (active_[phaseIdx]) {
                    pressure[phaseIdx] = pressure[phaseIdx] - pressure[BlackoilPhases::Liquid];
                }
            }

            // Since pcow = po - pw, but pcog = pg - po,
            // we have
            //   pw = po - pcow
            //   pg = po + pcgo
            // This is an unfortunate inconsistency, but a convention we must handle.
            for (int phaseIdx = 0; phaseIdx < BlackoilPhases::MaxNumPhases; ++phaseIdx) {
                if (active_[phaseIdx]) {
                    if (phaseIdx == BlackoilPhases::Aqua) {
                        pressure[phaseIdx] = po - pressure[phaseIdx];
                    } else {
                        pressure[phaseIdx] += po;
                    }
                }
            }

            return pressure;
        }

        V computeGasPressure(const V& po,
                           const V& sw,
                           const V& so,
                           const V& sg) const
        {
            assert (active_[Gas]);
            std::vector<ADB> cp = fluid_.capPress(ADB::constant(sw),
                                                  ADB::constant(so),
                                                  ADB::constant(sg),
                                                  cells_);
            return cp[Gas].value() + po;
        }

        const std::vector<PhasePresence>
        phaseCondition() const {return phaseCondition_;}

        /// update the primal variable for Sg, Rv or Rs. The Gas phase must
        /// be active to call this method.
        void
        updatePrimalVariableFromState(const ReservoirState& state)
        {
            updatePhaseCondFromPrimalVariable(state);
        }


        /// Update the phaseCondition_ member based on the primalVariable_ member.
        /// Also updates isRs_, isRv_ and isSg_;
        void
        updatePhaseCondFromPrimalVariable(const ReservoirState& state)
        {
            const int nc = Opm::AutoDiffGrid::numCells(grid_);
            isRs_ = V::Zero(nc);
            isRv_ = V::Zero(nc);
            isSg_ = V::Zero(nc);

            if (! (active_[Gas] && active_[Oil])) {
                // updatePhaseCondFromPrimarVariable() logic requires active gas and oil phase.
                phaseCondition_.assign(nc, PhasePresence());
                return;
            }
            for (int c = 0; c < nc; ++c) {
                phaseCondition_[c] = PhasePresence(); // No free phases.
                phaseCondition_[c].setFreeWater(); // Not necessary for property calculation usage.
                switch (state.hydroCarbonState()[c]) {
                case HydroCarbonState::GasAndOil:
                    phaseCondition_[c].setFreeOil();
                    phaseCondition_[c].setFreeGas();
                    isSg_[c] = 1;
                    break;
                case HydroCarbonState::OilOnly:
                    phaseCondition_[c].setFreeOil();
                    isRs_[c] = 1;
                    break;
                case HydroCarbonState::GasOnly:
                    phaseCondition_[c].setFreeGas();
                    isRv_[c] = 1;
                    break;
                default:
                    OPM_THROW(std::logic_error, "Unknown primary variable enum value in cell " << c << ": " << state.hydroCarbonState()[c]);
                }
            }
        }


        // TODO: added since the interfaces of the function are different
        // TODO: for StandardWells and MultisegmentWells
//        void
//        computeWellConnectionPressures(const SolutionState& state,
//                                       const WellState& well_state)
//        {
//            wellModel().computeWellConnectionPressures(state, well_state);
//        }


        /// \brief Compute the reduction within the convergence check.
        /// \param[in] B     A matrix with MaxNumPhases columns and the same number rows
        ///                  as the number of cells of the grid. B.col(i) contains the values
        ///                  for phase i.
        /// \param[in] tempV A matrix with MaxNumPhases columns and the same number rows
        ///                  as the number of cells of the grid. tempV.col(i) contains the
        ///                   values
        ///                  for phase i.
        /// \param[in] R     A matrix with MaxNumPhases columns and the same number rows
        ///                  as the number of cells of the grid. B.col(i) contains the values
        ///                  for phase i.
        /// \param[out] R_sum An array of size MaxNumPhases where entry i contains the sum
        ///                   of R for the phase i.
        /// \param[out] maxCoeff An array of size MaxNumPhases where entry i contains the
        ///                   maximum of tempV for the phase i.
        /// \param[out] B_avg An array of size MaxNumPhases where entry i contains the average
        ///                   of B for the phase i.
        /// \param[out] maxNormWell The maximum of the well flux equations for each phase.
        /// \param[in]  nc    The number of cells of the local grid.
        /// \return The total pore volume over all cells.
        double
        convergenceReduction(const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& B,
                             const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& tempV,
                             const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& R,
                             std::vector<double>& R_sum,
                             std::vector<double>& maxCoeff,
                             std::vector<double>& B_avg,
                             std::vector<double>& maxNormWell,
                             int nc) const
        {
            const int np = numPhases();
            const int nw = residual_.well_flux_eq.size() / np;
            assert(nw * np == int(residual_.well_flux_eq.size()));

            // Do the global reductions
#if HAVE_MPI
            if ( linsolver_.parallelInformation().type() == typeid(ParallelISTLInformation) )
            {
                const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>(linsolver_.parallelInformation());

                // Compute the global number of cells and porevolume
                std::vector<int> v(nc, 1);
                auto nc_and_pv = std::tuple<int, double>(0, 0.0);
                auto nc_and_pv_operators = std::make_tuple(Opm::Reduction::makeGlobalSumFunctor<int>(),
                                                           Opm::Reduction::makeGlobalSumFunctor<double>());
                auto nc_and_pv_containers  = std::make_tuple(v, geo_.poreVolume());
                info.computeReduction(nc_and_pv_containers, nc_and_pv_operators, nc_and_pv);

                for ( int idx = 0; idx < np; ++idx )
                {
                    auto values     = std::tuple<double,double,double>(0.0 ,0.0 ,0.0);
                    auto containers = std::make_tuple(B.col(idx),
                                                      tempV.col(idx),
                                                      R.col(idx));
                    auto operators  = std::make_tuple(Opm::Reduction::makeGlobalSumFunctor<double>(),
                                                      Opm::Reduction::makeGlobalMaxFunctor<double>(),
                                                      Opm::Reduction::makeGlobalSumFunctor<double>());
                    info.computeReduction(containers, operators, values);
                    B_avg[idx]       = std::get<0>(values)/std::get<0>(nc_and_pv);
                    maxCoeff[idx]    = std::get<1>(values);
                    R_sum[idx]       = std::get<2>(values);
                    assert(np >= np);
                    if (idx < np) {
                        maxNormWell[idx] = 0.0;
                        for ( int w = 0; w < nw; ++w ) {
                            maxNormWell[idx]  = std::max(maxNormWell[idx], std::abs(residual_.well_flux_eq.value()[nw*idx + w]));
                        }
                    }
                }
                info.communicator().max(maxNormWell.data(), np);
                // Compute pore volume
                return std::get<1>(nc_and_pv);
            }
            else
#endif
            {
                B_avg.resize(np);
                maxCoeff.resize(np);
                R_sum.resize(np);
                maxNormWell.resize(np);
                for ( int idx = 0; idx < np; ++idx )
                {
                    B_avg[idx] = B.col(idx).sum()/nc;
                    maxCoeff[idx] = tempV.col(idx).maxCoeff();
                    R_sum[idx] = R.col(idx).sum();

                    assert(np >= np);
                    if (idx < np) {
                        maxNormWell[idx] = 0.0;
                        for ( int w = 0; w < nw; ++w ) {
                            maxNormWell[idx] = std::max(maxNormWell[idx], std::abs(residual_.well_flux_eq.value()[nw*idx + w]));
                        }
                    }
                }
                // Compute total pore volume
                return geo_.poreVolume().sum();
            }
        }


        double dpMaxRel() const { return param_.dp_max_rel_; }
        double dsMax() const { return param_.ds_max_; }
        double drMaxRel() const { return param_.dr_max_rel_; }
        double maxResidualAllowed() const { return param_.max_residual_allowed_; }

    public:
        bool isBeginReportStep_;

    };
} // namespace Opm

#endif // OPM_BLACKOILMODELBASE_IMPL_HEADER_INCLUDED
