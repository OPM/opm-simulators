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

#include <dune/istl/solvers.hh>

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

        typedef double Scalar;
        typedef Dune::FieldVector<Scalar, 3    >       VectorBlockType;
        typedef Dune::FieldMatrix<Scalar, 3, 3 >      MatrixBlockType;
        typedef Dune::BCRSMatrix <MatrixBlockType>      Mat;
        typedef Dune::BlockVector<VectorBlockType>      BVector;
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
                          const StandardWellsDense<FluidSystem, BlackoilIndices>&                well_model,
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
                //residual_.singlePrecision = (unit::convert::to(dt, unit::day) < 20.) ;

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

            // -------- Mass balance equations --------
            assembleMassBalanceEq(timer, iterationIdx, reservoir_state);
            // -------- Well equations ----------
            double dt = timer.currentStepLength();

            IterationReport iter_report;
            try
            {
                iter_report = wellModel().assemble(ebosSimulator_, iterationIdx, dt, well_state, residual_);
            }
            catch ( const Dune::FMatrixError& e  )
            {
                OPM_THROW(Opm::NumericalProblem,"no convergence");
            }

            typedef double Scalar;
            typedef Dune::FieldVector<Scalar, 3    >       VectorBlockType;
            typedef Dune::FieldMatrix<Scalar, 3, 3 >      MatrixBlockType;
            typedef Dune::BCRSMatrix <MatrixBlockType>      Mat;
            typedef Dune::BlockVector<VectorBlockType>      BVector;

            typedef Dune::MatrixAdapter<Mat,BVector,BVector> Operator;
            auto& ebosJac = ebosSimulator_.model().linearizer().matrix();
            auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            convertResults2(ebosResid, ebosJac);
            wellModel().addRhs(ebosResid, ebosJac);

            return iter_report;
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

            typedef double Scalar;
            typedef Dune::FieldVector<Scalar, 3    >       VectorBlockType;
            typedef Dune::FieldMatrix<Scalar, 3, 3 >      MatrixBlockType;
            typedef Dune::BCRSMatrix <MatrixBlockType>      Mat;
            typedef Dune::BlockVector<VectorBlockType>      BVector;

            typedef Dune::MatrixAdapter<Mat,BVector,BVector> Operator;
            auto& ebosJac = ebosSimulator_.model().linearizer().matrix();
            auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            Operator opA(ebosJac);
            const double relax = 0.9;
            typedef Dune::SeqILU0<Mat, BVector, BVector> SeqPreconditioner;
            SeqPreconditioner precond(opA.getmat(), relax);
            Dune::SeqScalarProduct<BVector> sp;


            wellModel().apply(ebosJac, ebosResid);
            Dune::BiCGSTABSolver<BVector> linsolve(opA, sp, precond,
                                                   0.001,
                                                   100,
                                                   false);


            const int np = numPhases();
            const int nc = AutoDiffGrid::numCells(grid_);
//            std::cout << "ebosResid" << std::endl;
//            for( int p=0; p<np; ++p) {
//                for (int cell = 0 ; cell < 1; ++cell) {
//                    std::cout << ebosResid[cell][p] << std::endl;
//                }
//            }



            // Solve system.
            Dune::InverseOperatorResult result;
            BVector x(ebosJac.M());
            x = 0.0;
            linsolve.apply(x, ebosResid, result);

            const int nw = wellModel().wells().number_of_wells;
            BVector xw(nw);
            wellModel().recoverVariable(x, xw);

            // convert to V;


            V dx( (nc + nw) * np);
            for( int p=0; p<np; ++p) {
                for (int i = 0; i < nc; ++i) {
                    int idx = i + nc*ebosCompToFlowPhaseIdx(p);
                    dx(idx) = x[i][p];
                }
                for (int w = 0; w < nw; ++w) {
                    int idx = w + nw*p + nc*np;
                    dx(idx) = xw[w][flowPhaseToEbosCompIdx(p)];
                }
            }

            //V dx2 = linsolver_.computeNewtonIncrement(residual_);
            //std::cout << "------dx------- " << std::endl;
            //std::cout << dx << std::endl;
            //std::cout << "------dx2------- " << std::endl;
            //std::cout << dx2 << std::endl;

            //return dx;
            return dx;
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


            wellModel().updateWellState(dwells, well_state);

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
            Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> R2(nc, np);
            Eigen::Array<V::Scalar, Eigen::Dynamic, Eigen::Dynamic> tempV(nc, np);

            auto ebosResid = ebosSimulator_.model().linearizer().residual();
            for ( int idx = 0; idx < np; ++idx )
            {
                V b(nc);
                V r(nc);
                for (int cell_idx = 0; cell_idx < nc; ++cell_idx) {
                    const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                    const auto& fs = intQuants.fluidState();

                    int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(idx);
                    int ebosCompIdx = flowPhaseToEbosCompIdx(idx);

                    b[cell_idx] = 1 / fs.invB(ebosPhaseIdx).value;
                    r[cell_idx] = ebosResid[cell_idx][ebosCompIdx];

                }
                R2.col(idx) = r;
                B.col(idx) = b;
            }

            for ( int idx = 0; idx < np; ++idx )
            {
                //const ADB& tempB = rq_[idx].b;
                //B.col(idx)       = 1./tempB.value();
                //R.col(idx)       = residual_.material_balance_eq[idx].value();
                tempV.col(idx)   = R2.col(idx).abs()/pv;
//                std::cout << "------R------- " << idx << std::endl;
//                std::cout << R.col(idx)[0] << std::endl;
//                std::cout << "------R2------- " << idx << std::endl;
//                std::cout << R2.col(idx)[0] << std::endl;
            }

            std::vector<double> pv_vector (geo_.poreVolume().data(), geo_.poreVolume().data() + geo_.poreVolume().size());
            const double pvSum = detail::convergenceReduction(B, tempV, R2,
                                                      R_sum, maxCoeff, B_avg, maxNormWell,
                                                      nc, np, pv_vector, wellModel().residual());

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
        StandardWellsDense<FluidSystem, BlackoilIndices> well_model_;

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



        // ---------  Protected methods  ---------

    public:


        /// return the StandardWells object
        StandardWellsDense<FluidSystem, BlackoilIndices>& wellModel() { return well_model_; }
        const StandardWellsDense<FluidSystem, BlackoilIndices>& wellModel() const { return well_model_; }

        /// return the Well struct in the StandardWells
        const Wells& wells() const { return well_model_.wells(); }

        /// return true if wells are available in the reservoir
        bool wellsActive() const { return well_model_.wellsActive(); }

        /// return true if wells are available on this process
        bool localWellsActive() const { return well_model_.localWellsActive(); }


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
        int ebosCompToFlowPhaseIdx( const int compIdx ) const
        {
            const int compToPhase[ 3 ] = { Oil, Water, Gas };
            return compToPhase[ compIdx ];
        }

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
        void convertResults2(BVector& ebosResid, Mat& ebosJac) const
        {
            const int numPhases = wells().number_of_phases;
            const int numCells = ebosJac.N();
            const int cols = ebosJac.M();
            assert( numCells == cols );

            // write the right-hand-side values from the ebosJac into the objects
            // allocated above.
            const auto endrow = ebosJac.end();
            for( int cellIdx = 0; cellIdx < numCells; ++cellIdx )
            {
                const double cellVolume = ebosSimulator_.model().dofTotalVolume(cellIdx);
                auto& cellRes = ebosResid[ cellIdx ];

                for( int flowPhaseIdx = 0; flowPhaseIdx < numPhases; ++flowPhaseIdx )
                {
                    const double refDens = FluidSystem::referenceDensity( flowPhaseToEbosPhaseIdx( flowPhaseIdx ), 0 );
                    cellRes[ flowPhaseToEbosCompIdx( flowPhaseIdx ) ] /= refDens;
                    cellRes[ flowPhaseToEbosCompIdx( flowPhaseIdx ) ] *= cellVolume;
                }
            }

            for( auto row = ebosJac.begin(); row != endrow; ++row )
            {
                const int rowIdx = row.index();
                const double cellVolume = ebosSimulator_.model().dofTotalVolume(rowIdx);


                // translate the Jacobian of the residual from the format used by ebos to
                // the one expected by flow
                const auto endcol = row->end();
                for( auto col = row->begin(); col != endcol; ++col )
                {
                    for( int flowPhaseIdx = 0; flowPhaseIdx < numPhases; ++flowPhaseIdx )
                    {
                        const double refDens = FluidSystem::referenceDensity( flowPhaseToEbosPhaseIdx( flowPhaseIdx ), 0 );
                        for( int pvIdx=0; pvIdx<numPhases; ++pvIdx )
                        {
                            (*col)[flowPhaseToEbosCompIdx(flowPhaseIdx)][flowToEbosPvIdx(pvIdx)] /= refDens;
                            (*col)[flowPhaseToEbosCompIdx(flowPhaseIdx)][flowToEbosPvIdx(pvIdx)] *= cellVolume;
                        }
                    }
                }
            }
        }

        int flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
        {
            const int flowToEbos[ 3 ] = { FluidSystem::waterPhaseIdx, FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx };
            return flowToEbos[ phaseIdx ];
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

            //convertResults(ebosSimulator_);

            if (param_.update_equations_scaling_) {
                std::cout << "scaling" << std::endl;
                //updateEquationsScaling();
            }

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




        double dpMaxRel() const { return param_.dp_max_rel_; }
        double dsMax() const { return param_.ds_max_; }
        double drMaxRel() const { return param_.dr_max_rel_; }
        double maxResidualAllowed() const { return param_.max_residual_allowed_; }

    public:
        bool isBeginReportStep_;

    };
} // namespace Opm

#endif // OPM_BLACKOILMODELBASE_IMPL_HEADER_INCLUDED
