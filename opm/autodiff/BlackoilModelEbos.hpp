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
#include <opm/autodiff/NewtonIterationBlackoilInterface.hpp>

#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/core/well_controls.h>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/simulator/SimulatorTimer.hpp>
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
        typedef BlackoilModelEbos ThisType;
    public:
        // ---------  Types and enums  ---------
        typedef BlackoilState ReservoirState;
        typedef WellStateFullyImplicitBlackoilDense WellState;
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
        typedef Dune::FieldVector<Scalar, 3    >        VectorBlockType;
        typedef Dune::FieldMatrix<Scalar, 3, 3 >        MatrixBlockType;
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
        , vfp_properties_(
            eclState().getTableManager().getVFPInjTables(),
            eclState().getTableManager().getVFPProdTables())
        , linsolver_ (linsolver)
        , active_(detail::activePhases(fluid.phaseUsage()))
        , has_disgas_(FluidSystem::enableDissolvedGas())
        , has_vapoil_(FluidSystem::enableVaporizedOil())
        , param_( param )
        , well_model_ (well_model)
        , terminal_output_ (terminal_output)
        , current_relaxation_(1.0)
        , dx_old_(AutoDiffGrid::numCells(grid_))
        , isBeginReportStep_(false)
        , isRestart_(false)
        {
            const double gravity = detail::getGravity(geo_.gravity(), UgGridHelpers::dimensions(grid_));
            const std::vector<double> pv(geo_.poreVolume().data(), geo_.poreVolume().data() + geo_.poreVolume().size());
            const std::vector<double> depth(geo_.z().data(), geo_.z().data() + geo_.z().size());
            well_model_.init(&fluid_, &active_, &vfp_properties_, gravity, depth, pv);
            wellModel().setWellsActive( localWellsActive() );
            global_nc_    =  Opm::AutoDiffGrid::numCells(grid_);
        }

        const EclipseState& eclState() const
        { return *ebosSimulator_.gridManager().eclState(); }

        /// Called once before each time step.
        /// \param[in] timer                  simulation timer
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void prepareStep(const SimulatorTimerInterface& /*timer*/,
                         const ReservoirState& /*reservoir_state*/,
                         const WellState& /* well_state */)
        {
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
            if (iteration == 0) {
                // For each iteration we store in a vector the norms of the residual of
                // the mass balance for each active phase, the well flux and the well equations.
                residual_norms_history_.clear();
                current_relaxation_ = 1.0;
                dx_old_ = 0.0;
            }
            IterationReport iter_report = assemble(timer, iteration, reservoir_state, well_state);
            std::vector<double> residual_norms;
            const bool converged = getConvergence(timer, iteration,residual_norms);
            residual_norms_history_.push_back(residual_norms);
            bool must_solve = (iteration < nonlinear_solver.minIter()) || (!converged);
            // is first set to true if a linear solve is needed, but then it is set to false if the solver succeed.
            isRestart_ = must_solve && (iteration == nonlinear_solver.maxIter());
            // don't solve if we have reached the maximum number of iteration.
            must_solve = must_solve && (iteration < nonlinear_solver.maxIter());
            Dune::InverseOperatorResult result;
            if (must_solve) {
                // enable single precision for solvers when dt is smaller then 20 days
                //residual_.singlePrecision = (unit::convert::to(dt, unit::day) < 20.) ;

                // Compute the nonlinear update.
                const int nc = AutoDiffGrid::numCells(grid_);
                const int nw = wellModel().wells().number_of_wells;
                BVector x(nc);
                BVector xw(nw);
                solveJacobianSystem(result, x, xw);

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
                nonlinear_solver.stabilizeNonlinearUpdate(x, dx_old_, current_relaxation_);

                // Apply the update, applying model-dependent
                // limitations and chopping of the update.
                updateState(x,reservoir_state);
                wellModel().updateWellState(xw, well_state);

                // since the solution was changed, the cache for the intensive quantities
                // are invalid
                ebosSimulator_.model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);

                // solver has succeed i.e. no need for restart.
                isRestart_ = false;
            }
            const bool failed = false; // Not needed in this model.
            const int linear_iters = must_solve ? result.iterations : 0;
            return IterationReport{ failed, converged, linear_iters, iter_report.well_iterations };
        }
        void printIf(int c, double x, double y, double eps, std::string type) {
            if (std::abs(x-y) > eps) {
                std::cout << type << " " <<c << ": "<<x << " " << y << std::endl;
            }
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
                iter_report = wellModel().assemble(ebosSimulator_, iterationIdx, dt, well_state);
            }
            catch ( const Dune::FMatrixError& e  )
            {
                isRestart_ = true;
                OPM_THROW(Opm::NumericalProblem,"Well equation did not converge");
            }

            return iter_report;
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
            const int nc = Opm::AutoDiffGrid::numCells(grid_);
            const int nw = wellModel().wells().number_of_wells;
            return numPhases() * (nc + nw);
        }

        /// Number of linear iterations used in last call to solveJacobianSystem().
        int linearIterationsLastSolve() const
        {
            return linsolver_.iterations();
        }

        template <class X, class Y>
        void applyWellModel(const X& x, Y& y )
        {
           wellModel().apply(x, y);
        }


        /// Solve the Jacobian system Jx = r where J is the Jacobian and
        /// r is the residual.
        void solveJacobianSystem(Dune::InverseOperatorResult& result, BVector& x, BVector& xw) const
        {

            typedef double Scalar;
            typedef Dune::FieldVector<Scalar, 3    >       VectorBlockType;
            typedef Dune::FieldMatrix<Scalar, 3, 3 >      MatrixBlockType;
            typedef Dune::BCRSMatrix <MatrixBlockType>      Mat;
            typedef Dune::BlockVector<VectorBlockType>      BVector;

            const auto& ebosJac = ebosSimulator_.model().linearizer().matrix();
            auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            typedef WellModelMatrixAdapter<Mat,BVector,BVector, ThisType> Operator;
            Operator opA(ebosJac, const_cast< ThisType&  > (*this));
            const double relax = 0.9;
            typedef Dune::SeqILU0<Mat, BVector, BVector> SeqPreconditioner;
            SeqPreconditioner precond(opA.getmat(), relax);
            Dune::SeqScalarProduct<BVector> sp;

            // apply well residual to the residual.
            wellModel().apply(ebosResid);

            Dune::BiCGSTABSolver<BVector> linsolve(opA, sp, precond,
                                                   0.01,
                                                   100,
                                                   false);
            // Solve system.
            x = 0.0;
            linsolve.apply(x, ebosResid, result);

            // recover wells.
            xw = 0.0;
            wellModel().recoverVariable(x, xw);
        }

        //=====================================================================
        // Implementation for ISTL-matrix based operator
        //=====================================================================

        /*!
           \brief Adapter to turn a matrix into a linear operator.

           Adapts a matrix to the assembled linear operator interface
         */
        template<class M, class X, class Y, class WellModel>
        class WellModelMatrixAdapter : public Dune::MatrixAdapter<M,X,Y>
        {
          typedef Dune::MatrixAdapter<M,X,Y> BaseType;
        public:
          //! export types
          typedef M matrix_type;
          typedef X domain_type;
          typedef Y range_type;
          typedef typename X::field_type field_type;

          //! constructor: just store a reference to a matrix
          explicit WellModelMatrixAdapter (const M& A, WellModel& wellMod ) : BaseType( A ), wellMod_( wellMod ) {}

          //! apply operator to x:  \f$ y = A(x) \f$
          virtual void apply (const X& x, Y& y) const
          {
            BaseType::apply( x, y );

            wellMod_.applyWellModel(x, y );
          }

        private:
          WellModel& wellMod_;
        };

        /** @} end documentation */


        /// Apply an update to the primary variables, chopped if appropriate.
        /// \param[in]      dx                updates to apply to primary variables
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void updateState(const BVector& dx,
                         ReservoirState& reservoir_state)
        {
            using namespace Opm::AutoDiffGrid;
            const int np = fluid_.numPhases();
            const int nc = numCells(grid_);

            for (int cell_idx = 0; cell_idx < nc; ++cell_idx) {
                double dp = dx[cell_idx][flowPhaseToEbosCompIdx(0)];
                //reservoir_state.pressure()[cell_idx] -= dp;
                double& p = reservoir_state.pressure()[cell_idx];
                p -= dp;
                p = std::max(p, 1e5);

                // Saturation updates.
                const double dsw = active_[Water] ? dx[cell_idx][flowPhaseToEbosCompIdx(1)] : 0.0;
                const int xvar_ind = active_[Water] ? 2 : 1;
                const double dxvar = active_[Gas] ? dx[cell_idx][flowPhaseToEbosCompIdx(xvar_ind)] : 0.0;

                double dso = 0.0;
                double dsg = 0.0;
                double drs = 0.0;
                double drv = 0.0;

                double maxVal = 0.0;
                // water phase
                maxVal = std::max(std::abs(dsw),maxVal);
                dso -= dsw;
                // gas phase
                switch (reservoir_state.hydroCarbonState()[cell_idx]) {
                case HydroCarbonState::GasAndOil:
                    dsg = dxvar;
                    break;
                case HydroCarbonState::OilOnly:
                    drs = dxvar;
                    break;
                case HydroCarbonState::GasOnly:
                    dsg -= dsw;
                    drv = dxvar;
                    break;
                default:
                    OPM_THROW(std::logic_error, "Unknown primary variable enum value in cell " << cell_idx << ": " << reservoir_state.hydroCarbonState()[cell_idx]);
                }
                dso -= dsg;

                // Appleyard chop process.
                maxVal = std::max(std::abs(dsg),maxVal);
                double step = dsMax()/maxVal;
                step = std::min(step, 1.0);


                const Opm::PhaseUsage& pu = fluid_.phaseUsage();
                if (active_[Water]) {
                    double& sw = reservoir_state.saturation()[cell_idx*np + pu.phase_pos[ Water ]];
                    sw -= step * dsw;
                }
                if (active_[Gas]) {
                    double& sg = reservoir_state.saturation()[cell_idx*np + pu.phase_pos[ Gas ]];
                    sg -= step * dsg;
                }
                double& so = reservoir_state.saturation()[cell_idx*np + pu.phase_pos[ Oil ]];
                so -= step * dso;

                // const double drmaxrel = drMaxRel();
                // Update rs and rv
                if (has_disgas_) {
                    double& rs = reservoir_state.gasoilratio()[cell_idx];
                    rs -= drs;
                    rs = std::max(rs, 0.0);

                }
                if (has_vapoil_) {
                    double& rv = reservoir_state.rv()[cell_idx];
                    rv -= drv;
                    rv = std::max(rv, 0.0);
                }

                // Sg is used as primal variable for water only cells.
                const double epsilon = 1e-4; //std::sqrt(std::numeric_limits<double>::epsilon());
                double& sw = reservoir_state.saturation()[cell_idx*np + pu.phase_pos[ Water ]];
                double& sg = reservoir_state.saturation()[cell_idx*np + pu.phase_pos[ Gas ]];
                double& rs = reservoir_state.gasoilratio()[cell_idx];
                double& rv = reservoir_state.rv()[cell_idx];


                // phase translation sg <-> rs
                const HydroCarbonState hydroCarbonState = reservoir_state.hydroCarbonState()[cell_idx];
                const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                const auto& fs = intQuants.fluidState();
                switch (hydroCarbonState) {
                case HydroCarbonState::GasAndOil: {

                    if (sw > (1.0 - epsilon)) // water only i.e. do nothing
                        break;

                    if (sg <= 0.0 && has_disgas_) {
                        reservoir_state.hydroCarbonState()[cell_idx] = HydroCarbonState::OilOnly; // sg --> rs
                        sg = 0;
                        so = 1.0 - sw - sg;
                        double rsSat = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), reservoir_state.temperature()[cell_idx], reservoir_state.pressure()[cell_idx]);
                        double& rs = reservoir_state.gasoilratio()[cell_idx];
                        rs = rsSat*(1-epsilon);
                    } else if (so <= 0.0 && has_vapoil_) {
                        reservoir_state.hydroCarbonState()[cell_idx] = HydroCarbonState::GasOnly; // sg --> rv
                        so = 0;
                        sg = 1.0 - sw - so;
                        double& rv = reservoir_state.rv()[cell_idx];
                        // use gas pressure?
                        double rvSat = FluidSystem::gasPvt().saturatedOilVaporizationFactor(fs.pvtRegionIndex(), reservoir_state.temperature()[cell_idx], reservoir_state.pressure()[cell_idx]);
                        rv = rvSat*(1-epsilon);
                    }
                    break;
                }
                case HydroCarbonState::OilOnly: {
                    if (sw > (1.0 - epsilon)) {
                        // water only change to Sg
                        rs = 0;
                        rv = 0;
                        reservoir_state.hydroCarbonState()[cell_idx] = HydroCarbonState::GasAndOil;
                        //std::cout << "watonly rv -> sg" << cell_idx << std::endl;
                        break;
                    }




                    double rsSat = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), reservoir_state.temperature()[cell_idx], reservoir_state.pressure()[cell_idx]);
                    if (rs > ( rsSat * (1+epsilon) ) ) {
                        reservoir_state.hydroCarbonState()[cell_idx] = HydroCarbonState::GasAndOil;
                        sg = epsilon;
                        so -= epsilon;
                        rs = rsSat;
                    }
                    break;
                }
                case HydroCarbonState::GasOnly: {
                    if (sw > (1.0 - epsilon)) {
                        // water only change to Sg
                        rs = 0;
                        rv = 0;
                        reservoir_state.hydroCarbonState()[cell_idx] = HydroCarbonState::GasAndOil;
                        //std::cout << "watonly rv -> sg" << cell_idx << std::endl;
                        break;
                    }
                    double rvSat = FluidSystem::gasPvt().saturatedOilVaporizationFactor(fs.pvtRegionIndex(), reservoir_state.temperature()[cell_idx], reservoir_state.pressure()[cell_idx]);
                    if (rv > rvSat * (1+epsilon) ) {
                        reservoir_state.hydroCarbonState()[cell_idx] = HydroCarbonState::GasAndOil;
                        so = epsilon;
                        rv = rvSat;
                        sg -= epsilon;
                    }
                    break;
                }

                default:
                    OPM_THROW(std::logic_error, "Unknown primary variable enum value in cell " << cell_idx << ": " << hydroCarbonState);
                }
            }

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
        bool getConvergence(const SimulatorTimerInterface& timer, const int iteration, std::vector<double>& residual_norms)
        {
            typedef std::vector< double > Vector;

            const double dt = timer.currentStepLength();
            const double tol_mb    = param_.tolerance_mb_;
            const double tol_cnv   = param_.tolerance_cnv_;
            const double tol_wells = param_.tolerance_wells_;

            const int nc = Opm::AutoDiffGrid::numCells(grid_);
            const int np = numPhases();

            const auto& pv = geo_.poreVolume();

            Vector R_sum(np);
            Vector B_avg(np);
            Vector maxCoeff(np);
            Vector maxNormWell(np);

            std::vector< Vector > B( np, Vector( nc ) );
            std::vector< Vector > R( np, Vector( nc ) );
            std::vector< Vector > R2( np, Vector( nc ) );
            std::vector< Vector > tempV( np, Vector( nc ) );

            const auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            for ( int idx = 0; idx < np; ++idx )
            {
                Vector& R2_idx = R2[ idx ];
                Vector& B_idx  = B[ idx ];
                const int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(idx);
                const int ebosCompIdx = flowPhaseToEbosCompIdx(idx);

                for (int cell_idx = 0; cell_idx < nc; ++cell_idx) {
                    const auto& intQuants = *(ebosSimulator_.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
                    const auto& fs = intQuants.fluidState();

                    B_idx [cell_idx] = 1 / fs.invB(ebosPhaseIdx).value();
                    R2_idx[cell_idx] = ebosResid[cell_idx][ebosCompIdx];
                }
            }

            for ( int idx = 0; idx < np; ++idx )
            {
                //tempV.col(idx)   = R2.col(idx).abs()/pv;
                Vector& tempV_idx = tempV[ idx ];
                Vector& R2_idx    = R2[ idx ];
                for( int cell_idx = 0; cell_idx < nc; ++cell_idx )
                {
                    tempV_idx[ cell_idx ] = std::abs( R2_idx[ cell_idx ] ) / pv[ cell_idx ];
                }
            }

            Vector pv_vector (geo_.poreVolume().data(), geo_.poreVolume().data() + geo_.poreVolume().size());
            Vector wellResidual =  wellModel().residual();
            const double pvSum = detail::convergenceReduction(B, tempV, R2,
                                                      R_sum, maxCoeff, B_avg, maxNormWell,
                                                      nc, np, pv_vector, wellResidual );

            Vector CNV(np);
            Vector mass_balance_residual(np);
            Vector well_flux_residual(np);

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
                residual_norms.push_back(CNV[idx]);
            }

            const bool converged = converged_MB && converged_CNV && converged_Well;

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
                ss.precision(oprec);
                ss.flags(oflags);
                OpmLog::note(ss.str());
            }

            for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                const auto& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));

                if (std::isnan(mass_balance_residual[phaseIdx])
                    || std::isnan(CNV[phaseIdx])
                    || (phaseIdx < np && std::isnan(well_flux_residual[phaseIdx]))) {
                    isRestart_ = true;
                    OPM_THROW(Opm::NumericalProblem, "NaN residual for phase " << phaseName);
                }
                if (mass_balance_residual[phaseIdx] > maxResidualAllowed()
                    || CNV[phaseIdx] > maxResidualAllowed()
                    || (phaseIdx < np && well_flux_residual[phaseIdx] > maxResidualAllowed())) {
                    isRestart_ = true;
                    OPM_THROW(Opm::NumericalProblem, "Too large residual for phase " << phaseName);
                }
            }

            return converged;
        }


        /// The number of active fluid phases in the model.
        int numPhases() const
        {
            return fluid_.numPhases();
        }

        std::vector<std::vector<double> >
        computeFluidInPlace(const ReservoirState& x,
                            const std::vector<int>& fipnum) const
        {
            OPM_THROW(std::logic_error,
                      "computeFluidInPlace() not implemented by BlackoilModelEbos!");
        }

        const Simulator& ebosSimulator() const
        { return ebosSimulator_; }

    protected:

        // ---------  Data members  ---------

        Simulator& ebosSimulator_;
        const Grid&         grid_;
        const BlackoilPropsAdInterface& fluid_;
        const DerivedGeology&           geo_;
        VFPProperties                   vfp_properties_;
        const NewtonIterationBlackoilInterface&    linsolver_;
        // For each canonical phase -> true if active
        const std::vector<bool>         active_;
        // Size = # active phases. Maps active -> canonical phase indices.
        const std::vector<int>          cells_;  // All grid cells
        const bool has_disgas_;
        const bool has_vapoil_;

        ModelParameters                 param_;

        // Well Model
        StandardWellsDense<FluidSystem, BlackoilIndices> well_model_;

        /// \brief Whether we print something to std::cout
        bool terminal_output_;
        /// \brief The number of cells of the global grid.
        int global_nc_;

        std::vector<std::vector<double>> residual_norms_history_;
        double current_relaxation_;
        BVector dx_old_;



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
                if( reservoirState.hydroCarbonState()[cellIdx] == HydroCarbonState::OilOnly && has_disgas_ )
                {
                    cellPv[BlackoilIndices::compositionSwitchIdx] = rs[cellIdx];
                    cellPv[BlackoilIndices::pressureSwitchIdx] = oilPressure[cellIdx];
                    cellPv.setPrimaryVarsMeaning( PrimaryVariables::Sw_po_Rs );
                }
                else if( reservoirState.hydroCarbonState()[cellIdx] == HydroCarbonState::GasOnly && has_vapoil_ )
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
                    assert( reservoirState.hydroCarbonState()[cellIdx] == HydroCarbonState::GasAndOil);
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

        void convertResults(BVector& ebosResid, Mat& ebosJac) const
        {
            const int numPhases = wells().number_of_phases;
            const int numCells = ebosJac.N();
            assert( numCells == static_cast<int>(ebosJac.M()) );

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
            // if the last step failed we want to recalculate the IntesiveQuantities.
            if (isRestart_) {
                ebosSimulator_.model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);
            }

            ebosSimulator_.problem().beginIteration();
            ebosSimulator_.model().linearizer().linearize();
            ebosSimulator_.problem().endIteration();

            prevEpisodeIdx = ebosSimulator_.episodeIndex();

            auto& ebosJac = ebosSimulator_.model().linearizer().matrix();
            auto& ebosResid = ebosSimulator_.model().linearizer().residual();
            convertResults(ebosResid, ebosJac);

            if (param_.update_equations_scaling_) {
                std::cout << "equation scaling not suported yet" << std::endl;
                //updateEquationsScaling();
            }

        }


        double dpMaxRel() const { return param_.dp_max_rel_; }
        double dsMax() const { return param_.ds_max_; }
        double drMaxRel() const { return param_.dr_max_rel_; }
        double maxResidualAllowed() const { return param_.max_residual_allowed_; }

    public:
        bool isBeginReportStep_;
        bool isRestart_;

    };
} // namespace Opm

#endif // OPM_BLACKOILMODELBASE_IMPL_HEADER_INCLUDED
