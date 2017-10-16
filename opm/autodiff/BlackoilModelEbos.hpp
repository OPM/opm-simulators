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

#include <ebos/eclproblem.hh>
#include <ewoms/common/start.hh>

#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/BlackoilWellModel.hpp>
#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/WellHelpers.hpp>
#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/autodiff/VFPProperties.hpp>
#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/autodiff/VFPInjProperties.hpp>
#include <opm/autodiff/BlackoilDetails.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>
#include <opm/autodiff/NewtonIterationBlackoilInterface.hpp>
#include <opm/autodiff/RateConverter.hpp>

#include <opm/core/grid.h>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/core/well_controls.h>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/autodiff/ISTLSolver.hpp>
#include <opm/common/data/SimulationDataContainer.hpp>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/timer.hh>
#include <dune/common/unused.hh>

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
SET_BOOL_PROP(EclFlowProblem, ExportGlobalTransmissibility, true);
// default in flow is to formulate the equations in surface volumes
SET_BOOL_PROP(EclFlowProblem, BlackoilConserveSurfaceVolume, true);
SET_BOOL_PROP(EclFlowProblem, UseVolumetricResidual, false);


// SWATINIT is done by the flow part of flow_ebos. this can be removed once the legacy
// code for fluid and satfunc handling gets fully retired.
SET_BOOL_PROP(EclFlowProblem, EnableSwatinit, false);
}}

namespace Opm {
    /// A model implementation for three-phase black oil.
    ///
    /// The simulator is capable of handling three-phase problems
    /// where gas can be dissolved in oil and vice versa. It
    /// uses an industry-standard TPFA discretization with per-phase
    /// upwind weighting of mobilities.
    template <class TypeTag>
    class BlackoilModelEbos
    {
    public:
        // ---------  Types and enums  ---------
        typedef BlackoilState ReservoirState;
        typedef WellStateFullyImplicitBlackoil WellState;
        typedef BlackoilModelParameters ModelParameters;

        typedef typename GET_PROP_TYPE(TypeTag, Simulator)         Simulator;
        typedef typename GET_PROP_TYPE(TypeTag, Grid)              Grid;
        typedef typename GET_PROP_TYPE(TypeTag, ElementContext)    ElementContext;
        typedef typename GET_PROP_TYPE(TypeTag, SolutionVector)    SolutionVector ;
        typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables)  PrimaryVariables ;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem)       FluidSystem;
        typedef typename GET_PROP_TYPE(TypeTag, Indices)           Indices;
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw)       MaterialLaw;
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

        typedef double Scalar;
        static const int numEq = Indices::numEq;
        static const int contiSolventEqIdx = Indices::contiSolventEqIdx;
        static const int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
        static const int solventSaturationIdx = Indices::solventSaturationIdx;
        static const int polymerConcentrationIdx = Indices::polymerConcentrationIdx;

        typedef Dune::FieldVector<Scalar, numEq >        VectorBlockType;
        typedef Dune::FieldMatrix<Scalar, numEq, numEq >        MatrixBlockType;
        typedef Dune::BCRSMatrix <MatrixBlockType>      Mat;
        typedef Dune::BlockVector<VectorBlockType>      BVector;

        typedef ISTLSolver< MatrixBlockType, VectorBlockType, Indices::pressureSwitchIdx >  ISTLSolverType;
        //typedef typename SolutionVector :: value_type            PrimaryVariables ;

        // For the conversion between the surface volume rate and resrevoir voidage rate
        using RateConverterType = RateConverter::
            SurfaceToReservoirVoidage<BlackoilPropsAdFromDeck::FluidSystem, std::vector<int> >;

        typedef Opm::FIPData FIPDataType;

        // ---------  Public methods  ---------

        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] wells            well structure
        /// \param[in] vfp_properties   Vertical flow performance tables
        /// \param[in] linsolver        linear solver
        /// \param[in] eclState         eclipse state
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilModelEbos(Simulator& ebosSimulator,
                          const ModelParameters& param,
                          BlackoilWellModel<TypeTag>& well_model,
                          RateConverterType& rate_converter,
                          const NewtonIterationBlackoilInterface& linsolver,
                          const bool terminal_output
                          )
        : ebosSimulator_(ebosSimulator)
        , grid_(ebosSimulator_.gridManager().grid())
        , istlSolver_( dynamic_cast< const ISTLSolverType* > (&linsolver) )
        , phaseUsage_(phaseUsageFromDeck(eclState()))
        , vfp_properties_(
            eclState().getTableManager().getVFPInjTables(),
            eclState().getTableManager().getVFPProdTables())
        , active_(detail::activePhases(phaseUsage_))
        , has_disgas_(FluidSystem::enableDissolvedGas())
        , has_vapoil_(FluidSystem::enableVaporizedOil())
        , has_solvent_(GET_PROP_VALUE(TypeTag, EnableSolvent))
        , has_polymer_(GET_PROP_VALUE(TypeTag, EnablePolymer))
        , param_( param )
        , well_model_ (well_model)
        , terminal_output_ (terminal_output)
        , rate_converter_(rate_converter)
        , current_relaxation_(1.0)
        , dx_old_(AutoDiffGrid::numCells(grid_))
        , isBeginReportStep_(false)
        {
            // Wells are active if they are active wells on at least
            // one process.
            int wellsActive = localWellsActive() ? 1 : 0;
            wellsActive = grid_.comm().max(wellsActive);
            wellModel().setWellsActive( wellsActive );
            // compute global sum of number of cells
            global_nc_ = detail::countGlobalCells(grid_);
            wellModel().setVFPProperties(&vfp_properties_);
            if (!istlSolver_)
            {
                OPM_THROW(std::logic_error,"solver down cast to ISTLSolver failed");
            }
        }

        bool isParallel() const
        { return  grid_.comm().size() > 1; }

        const EclipseState& eclState() const
        { return ebosSimulator_.gridManager().eclState(); }

        /// Called once before each time step.
        /// \param[in] timer                  simulation timer
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void prepareStep(const SimulatorTimerInterface& /*timer*/,
                         const ReservoirState& /*reservoir_state*/,
                         const WellState& /* well_state */)
        {
            if ( wellModel().wellCollection()->havingVREPGroups() ) {
                updateRateConverter();
            }

            unsigned numDof = ebosSimulator_.model().numGridDof();
            wasSwitched_.resize(numDof);
            std::fill(wasSwitched_.begin(), wasSwitched_.end(), false);
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
        SimulatorReport nonlinearIteration(const int iteration,
                                           const SimulatorTimerInterface& timer,
                                           NonlinearSolverType& nonlinear_solver,
                                           ReservoirState& /*reservoir_state*/,
                                           WellState& well_state)
        {
            SimulatorReport report;
            failureReport_ = SimulatorReport();
            Dune::Timer perfTimer;

            perfTimer.start();
            if (iteration == 0) {
                // For each iteration we store in a vector the norms of the residual of
                // the mass balance for each active phase, the well flux and the well equations.
                residual_norms_history_.clear();
                current_relaxation_ = 1.0;
                dx_old_ = 0.0;
            }

            report.total_linearizations = 1;

            try {
                report += assemble(timer, iteration, well_state);
                report.assemble_time += perfTimer.stop();
            }
            catch (...) {
                report.assemble_time += perfTimer.stop();
                failureReport_ += report;
                // todo (?): make the report an attribute of the class
                throw; // continue throwing the stick
            }

            std::vector<double> residual_norms;
            perfTimer.reset();
            perfTimer.start();
            // the step is not considered converged until at least minIter iterations is done
            report.converged = getConvergence(timer, iteration,residual_norms) && iteration > nonlinear_solver.minIter();

             // checking whether the group targets are converged
             if (wellModel().wellCollection()->groupControlActive()) {
                  report.converged = report.converged && wellModel().wellCollection()->groupTargetConverged(well_state.wellRates());
             }

            report.update_time += perfTimer.stop();
            residual_norms_history_.push_back(residual_norms);
            if (!report.converged) {
                perfTimer.reset();
                perfTimer.start();
                report.total_newton_iterations = 1;

                // enable single precision for solvers when dt is smaller then 20 days
                //residual_.singlePrecision = (unit::convert::to(dt, unit::day) < 20.) ;

                // Compute the nonlinear update.
                const int nc = AutoDiffGrid::numCells(grid_);
                const int nw = numWells();
                BVector x(nc);

                try {
                    solveJacobianSystem(x);
                    report.linear_solve_time += perfTimer.stop();
                    report.total_linear_iterations += linearIterationsLastSolve();
                }
                catch (...) {
                    report.linear_solve_time += perfTimer.stop();
                    report.total_linear_iterations += linearIterationsLastSolve();

                    failureReport_ += report;
                    throw; // re-throw up
                }

                perfTimer.reset();
                perfTimer.start();

                // handling well state update before oscillation treatment is a decision based
                // on observation to avoid some big performance degeneration under some circumstances.
                // there is no theorectical explanation which way is better for sure.

                if( nw > 0 )
                {
                    wellModel().recoverWellSolutionAndUpdateWellState(x, well_state);
                }

                if (param_.use_update_stabilization_) {
                    // Stabilize the nonlinear update.
                    bool isOscillate = false;
                    bool isStagnate = false;
                    nonlinear_solver.detectOscillations(residual_norms_history_, iteration, isOscillate, isStagnate);
                    if (isOscillate) {
                        current_relaxation_ -= nonlinear_solver.relaxIncrement();
                        current_relaxation_ = std::max(current_relaxation_, nonlinear_solver.relaxMax());
                        if (terminalOutputEnabled()) {
                            std::string msg = "    Oscillating behavior detected: Relaxation set to "
                                    + std::to_string(current_relaxation_);
                            OpmLog::info(msg);
                        }
                    }
                    nonlinear_solver.stabilizeNonlinearUpdate(x, dx_old_, current_relaxation_);
                }

                // Apply the update, with considering model-dependent limitations and
                // chopping of the update.
                updateState(x,iteration);

                report.update_time += perfTimer.stop();
            }

            return report;
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
            DUNE_UNUSED_PARAMETER(timer);
            DUNE_UNUSED_PARAMETER(reservoir_state);
            DUNE_UNUSED_PARAMETER(well_state);
        }

        /// Assemble the residual and Jacobian of the nonlinear system.
        /// \param[in]      reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        /// \param[in]      initial_assembly  pass true if this is the first call to assemble() in this timestep
        SimulatorReport assemble(const SimulatorTimerInterface& timer,
                                 const int iterationIdx,
                                 WellState& well_state)
        {
            using namespace Opm::AutoDiffGrid;

            SimulatorReport report;

            // when having VREP group control, update the rate converter based on reservoir state
            if ( wellModel().wellCollection()->havingVREPGroups() ) {
                updateRateConverter();
            }

            // -------- Mass balance equations --------
            assembleMassBalanceEq(timer, iterationIdx);

            // -------- Well equations ----------
            double dt = timer.currentStepLength();

            try
            {
                report = wellModel().assemble(ebosSimulator_, iterationIdx, dt, well_state);

                // apply well residual to the residual.
                auto& ebosResid = ebosSimulator_.model().linearizer().residual();
                wellModel().apply(ebosResid);
            }
            catch ( const Dune::FMatrixError& e  )
            {
                OPM_THROW(Opm::NumericalProblem,"Error encounted when solving well equations");
            }

            return report;
        }

        // compute the "relative" change of the solution between time steps
        template <class Dummy>
        double relativeChange(const Dummy&, const Dummy&) const
        {
            Scalar resultDelta = 0.0;
            Scalar resultDenom = 0.0;

            const auto& elemMapper = ebosSimulator_.model().elementMapper();
            const auto& gridView = ebosSimulator_.gridView();
            auto elemIt = gridView.template begin</*codim=*/0>();
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const auto& elem = *elemIt;
                if (elem.partitionType() != Dune::InteriorEntity)
                    continue;

		unsigned globalElemIdx = elemMapper.index(elem);
                const auto& priVarsNew = ebosSimulator_.model().solution(/*timeIdx=*/0)[globalElemIdx];

                Scalar pressureNew;
		pressureNew = priVarsNew[Indices::pressureSwitchIdx];

		Scalar saturationsNew[FluidSystem::numPhases] = { 0.0 };
                Scalar oilSaturationNew = 1.0;
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    saturationsNew[FluidSystem::waterPhaseIdx] = priVarsNew[Indices::waterSaturationIdx];
                    oilSaturationNew -= saturationsNew[FluidSystem::waterPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && priVarsNew.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
                    saturationsNew[FluidSystem::gasPhaseIdx] = priVarsNew[Indices::compositionSwitchIdx];
                    oilSaturationNew -= saturationsNew[FluidSystem::gasPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    saturationsNew[FluidSystem::oilPhaseIdx] = oilSaturationNew;
                }

                const auto& priVarsOld = ebosSimulator_.model().solution(/*timeIdx=*/1)[globalElemIdx];

                Scalar pressureOld;
                pressureOld = priVarsOld[Indices::pressureSwitchIdx];

                Scalar saturationsOld[FluidSystem::numPhases] = { 0.0 };
                Scalar oilSaturationOld = 1.0;
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    saturationsOld[FluidSystem::waterPhaseIdx] = priVarsOld[Indices::waterSaturationIdx];
                    oilSaturationOld -= saturationsOld[FluidSystem::waterPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && priVarsOld.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
                    saturationsOld[FluidSystem::gasPhaseIdx] = priVarsOld[Indices::compositionSwitchIdx];
                    oilSaturationOld -= saturationsOld[FluidSystem::gasPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    saturationsOld[FluidSystem::oilPhaseIdx] = oilSaturationOld;
                }

                Scalar tmp = pressureNew - pressureOld;
                resultDelta += tmp*tmp;
                resultDenom += pressureNew*pressureNew;

                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++ phaseIdx) {
                    Scalar tmp = saturationsNew[phaseIdx] - saturationsOld[phaseIdx];
		    resultDelta += tmp*tmp;
                    resultDenom += saturationsNew[phaseIdx]*saturationsNew[phaseIdx];
                }
            }

            resultDelta = gridView.comm().sum(resultDelta);
            resultDenom = gridView.comm().sum(resultDenom);

	    if (resultDenom > 0.0)
	      return resultDelta/resultDenom;
	    return 0.0;
        }


        /// The size (number of unknowns) of the nonlinear system of equations.
        int sizeNonLinear() const
        {
            const int nc = Opm::AutoDiffGrid::numCells(grid_);
            const int nw = numWells();
            return numComponents() * (nc + nw);
        }

        /// Number of linear iterations used in last call to solveJacobianSystem().
        int linearIterationsLastSolve() const
        {
            return istlSolver().iterations();
        }

        /// Solve the Jacobian system Jx = r where J is the Jacobian and
        /// r is the residual.
        void solveJacobianSystem(BVector& x) const
        {
            const auto& ebosJac = ebosSimulator_.model().linearizer().matrix();
            auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            // set initial guess
            x = 0.0;

            // Solve system.
            if( isParallel() )
            {
                typedef WellModelMatrixAdapter< Mat, BVector, BVector, BlackoilWellModel<TypeTag>, true > Operator;
                Operator opA(ebosJac, well_model_, istlSolver().parallelInformation() );
                assert( opA.comm() );
                istlSolver().solve( opA, x, ebosResid, *(opA.comm()) );
            }
            else
            {
                typedef WellModelMatrixAdapter< Mat, BVector, BVector, BlackoilWellModel<TypeTag>, false > Operator;
                Operator opA(ebosJac, well_model_);
                istlSolver().solve( opA, x, ebosResid );
            }
        }

        //=====================================================================
        // Implementation for ISTL-matrix based operator
        //=====================================================================

        /*!
           \brief Adapter to turn a matrix into a linear operator.

           Adapts a matrix to the assembled linear operator interface
         */
        template<class M, class X, class Y, class WellModel, bool overlapping >
        class WellModelMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
        {
          typedef Dune::AssembledLinearOperator<M,X,Y> BaseType;

        public:
          typedef M matrix_type;
          typedef X domain_type;
          typedef Y range_type;
          typedef typename X::field_type field_type;

#if HAVE_MPI
          typedef Dune::OwnerOverlapCopyCommunication<int,int> communication_type;
#else
          typedef Dune::CollectiveCommunication< Grid > communication_type;
#endif

          enum {
            //! \brief The solver category.
            category = overlapping ?
                Dune::SolverCategory::overlapping :
                Dune::SolverCategory::sequential
          };

          //! constructor: just store a reference to a matrix
          WellModelMatrixAdapter (const M& A, const WellModel& wellMod, const boost::any& parallelInformation = boost::any() )
              : A_( A ), wellMod_( wellMod ), comm_()
          {
#if HAVE_MPI
            if( parallelInformation.type() == typeid(ParallelISTLInformation) )
            {
              const ParallelISTLInformation& info =
                  boost::any_cast<const ParallelISTLInformation&>( parallelInformation);
              comm_.reset( new communication_type( info.communicator() ) );
            }
#endif
          }

          virtual void apply( const X& x, Y& y ) const
          {
            A_.mv( x, y );
            // add well model modification to y
            wellMod_.apply(x, y );

#if HAVE_MPI
            if( comm_ )
              comm_->project( y );
#endif
          }

          // y += \alpha * A * x
          virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
          {
            A_.usmv(alpha,x,y);
            // add scaled well model modification to y
            wellMod_.applyScaleAdd( alpha, x, y );

#if HAVE_MPI
            if( comm_ )
              comm_->project( y );
#endif
          }

          virtual const matrix_type& getmat() const { return A_; }

          communication_type* comm()
          {
              return comm_.operator->();
          }

        protected:
          const matrix_type& A_ ;
          const WellModel& wellMod_;
          std::unique_ptr< communication_type > comm_;
        };

        /// Apply an update to the primary variables, chopped if appropriate.
        /// \param[in]      dx                updates to apply to primary variables
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        void updateState(const BVector& dx,
                         const int iterationIdx)
        {
            using namespace Opm::AutoDiffGrid;

            const auto& ebosProblem = ebosSimulator_.problem();

            unsigned numSwitched = 0;
            ElementContext elemCtx( ebosSimulator_ );
            const auto& gridView = ebosSimulator_.gridView();
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            SolutionVector& solution = ebosSimulator_.model().solution( 0 /* timeIdx */ );

            // Store the initial solution.
            if( iterationIdx == 0 )
            {
                ebosSimulator_.model().solution( 1 /* timeIdx */ ) = solution;
            }

            for (auto elemIt = gridView.template begin</*codim=*/0>();
                 elemIt != elemEndIt;
                 ++elemIt)
            {
                const auto& elem = *elemIt;

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);

                PrimaryVariables& priVars = solution[ cell_idx ];

                const double& dp = dx[cell_idx][Indices::pressureSwitchIdx];
                double& p = priVars[Indices::pressureSwitchIdx];
                const double& dp_rel_max = dpMaxRel();
                const int sign_dp = dp > 0 ? 1: -1;
                p -= sign_dp * std::min(std::abs(dp), std::abs(p)*dp_rel_max);
                p = std::max(p, 0.0);

                // Saturation updates.
                const double dsw = active_[Water] ? dx[cell_idx][Indices::waterSaturationIdx] : 0.0;
                const double dxvar = active_[Gas] ? dx[cell_idx][Indices::compositionSwitchIdx] : 0.0;

                double dso = 0.0;
                double dsg = 0.0;
                double drs = 0.0;
                double drv = 0.0;

                // determine the saturation delta values
                if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
                    dsg = dxvar;
                }
                else if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Rs) {
                    drs = dxvar;
                }
                else {
                    assert(priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv);
                    drv = dxvar;
                    dsg = 0.0;
                }

                // solvent
                const double dss = has_solvent_ ? dx[cell_idx][Indices::solventSaturationIdx] : 0.0;

                // polymer
                const double dc = has_polymer_ ? dx[cell_idx][Indices::polymerConcentrationIdx] : 0.0;

                // oil
                dso = - (dsw + dsg + dss);

                // compute a scaling factor for the saturation update so that the maximum
                // allowed change of saturations between iterations is not exceeded
                double maxVal = 0.0;
                maxVal = std::max(std::abs(dsw),maxVal);
                maxVal = std::max(std::abs(dsg),maxVal);
                maxVal = std::max(std::abs(dso),maxVal);
                maxVal = std::max(std::abs(dss),maxVal);

                double satScaleFactor = 1.0;
                if (maxVal > dsMax()) {
                    satScaleFactor = dsMax()/maxVal;
                }

                if (active_[Water]) {
                    double& sw = priVars[Indices::waterSaturationIdx];
                    sw -= satScaleFactor * dsw;
                }

                if (active_[Gas]) {
                     if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
                           double& sg = priVars[Indices::compositionSwitchIdx];
                           sg -= satScaleFactor * dsg;
                     }
                }

                if (has_solvent_) {
                    double& ss = priVars[Indices::solventSaturationIdx];
                    ss -= satScaleFactor * dss;
                    ss = std::min(std::max(ss, 0.0),1.0);
                }
                if (has_polymer_) {
                    double& c = priVars[Indices::polymerConcentrationIdx];
                    c -= satScaleFactor * dc;
                    c = std::max(c, 0.0);
                }

                // Update rs and rv
                if (active_[Gas] && active_[Oil] ) {
                    unsigned pvtRegionIdx = ebosSimulator_.problem().pvtRegionIndex(cell_idx);
                    const double drmaxrel = drMaxRel();
                    if (has_disgas_) {
                        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Rs) {
                            Scalar RsSat =
                                FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx, 300.0, p);

                            double& rs = priVars[Indices::compositionSwitchIdx];
                            rs -= ((drs<0)?-1:1)*std::min(std::abs(drs), RsSat*drmaxrel);
                            rs = std::max(rs, 0.0);
                        }

                    }
                    if (has_vapoil_) {
                        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv) {
                            Scalar RvSat =
                                FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx, 300.0, p);

                            double& rv = priVars[Indices::compositionSwitchIdx];
                            rv -= ((drv<0)?-1:1)*std::min(std::abs(drv), RvSat*drmaxrel);
                            rv = std::max(rv, 0.0);
                        }
                    }
                }

               // Add an epsilon to make it harder to switch back immediately after the primary variable was changed.
                if (wasSwitched_[cell_idx])
                    wasSwitched_[cell_idx] = priVars.adaptPrimaryVariables(ebosProblem, cell_idx, 1e-5);
                else
                    wasSwitched_[cell_idx] = priVars.adaptPrimaryVariables(ebosProblem, cell_idx);

                if (wasSwitched_[cell_idx])
                    ++numSwitched;
            }

            // if the solution is updated the intensive Quantities need to be recalculated
            ebosSimulator_.model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);

        }

        /// Return true if output to cout is wanted.
        bool terminalOutputEnabled() const
        {
            return terminal_output_;
        }

        template <class CollectiveCommunication>
        double convergenceReduction(const CollectiveCommunication& comm,
                                    const double pvSumLocal,
                                    std::vector< Scalar >& R_sum,
                                    std::vector< Scalar >& maxCoeff,
                                    std::vector< Scalar >& B_avg)
        {
            // Compute total pore volume (use only owned entries)
            double pvSum = pvSumLocal;

            if( comm.size() > 1 )
            {
                // global reduction
                std::vector< Scalar > sumBuffer;
                std::vector< Scalar > maxBuffer;
                const int numComp = B_avg.size();
                sumBuffer.reserve( 2*numComp + 1 ); // +1 for pvSum
                maxBuffer.reserve( numComp );
                for( int compIdx = 0; compIdx < numComp; ++compIdx )
                {
                    sumBuffer.push_back( B_avg[ compIdx ] );
                    sumBuffer.push_back( R_sum[ compIdx ] );
                    maxBuffer.push_back( maxCoeff[ compIdx ] );
                }

                // Compute total pore volume
                sumBuffer.push_back( pvSum );

                // compute global sum
                comm.sum( sumBuffer.data(), sumBuffer.size() );

                // compute global max
                comm.max( maxBuffer.data(), maxBuffer.size() );

                // restore values to local variables
                for( int compIdx = 0, buffIdx = 0; compIdx < numComp; ++compIdx, ++buffIdx )
                {
                    B_avg[ compIdx ]    = sumBuffer[ buffIdx ];
                    ++buffIdx;

                    R_sum[ compIdx ]       = sumBuffer[ buffIdx ];
                }

                for( int compIdx = 0; compIdx < numComp; ++compIdx )
                {
                    maxCoeff[ compIdx ] = maxBuffer[ compIdx ];
                }

                // restore global pore volume
                pvSum = sumBuffer.back();
            }

            // return global pore volume
            return pvSum;
        }

        /// Compute convergence based on total mass balance (tol_mb) and maximum
        /// residual mass balance (tol_cnv).
        /// \param[in]   timer       simulation timer
        /// \param[in]   dt          timestep length
        /// \param[in]   iteration   current iteration number
        bool getConvergence(const SimulatorTimerInterface& timer, const int iteration, std::vector<double>& residual_norms)
        {
            typedef std::vector< Scalar > Vector;

            const double dt = timer.currentStepLength();
            const double tol_mb    = param_.tolerance_mb_;
            const double tol_cnv   = param_.tolerance_cnv_;

            const int np = numPhases();
            const int numComp = numComponents();

            Vector R_sum(numComp, 0.0 );
            Vector B_avg(numComp, 0.0 );
            Vector maxCoeff(numComp, std::numeric_limits< Scalar >::lowest() );

            const auto& ebosModel = ebosSimulator_.model();
            const auto& ebosProblem = ebosSimulator_.problem();

            const auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            ElementContext elemCtx(ebosSimulator_);
            const auto& gridView = ebosSimulator().gridView();
            const auto& elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();

            double pvSumLocal = 0.0;
            for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
                 elemIt != elemEndIt;
                 ++elemIt)
            {
                const auto& elem = *elemIt;
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();

                const double pvValue = ebosProblem.porosity(cell_idx) * ebosModel.dofTotalVolume( cell_idx );
                pvSumLocal += pvValue;

                for ( int phaseIdx = 0; phaseIdx < np; ++phaseIdx )
                {
                    const int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phaseIdx);
                    const int ebosCompIdx = flowPhaseToEbosCompIdx(phaseIdx);

                    B_avg[ phaseIdx ] += 1.0 / fs.invB(ebosPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][ebosCompIdx];

                    R_sum[ phaseIdx ] += R2;
                    maxCoeff[ phaseIdx ] = std::max( maxCoeff[ phaseIdx ], std::abs( R2 ) / pvValue );
                }

                if ( has_solvent_ ) {
                    B_avg[ contiSolventEqIdx ] += 1.0 / intQuants.solventInverseFormationVolumeFactor().value();
                    const auto R2 = ebosResid[cell_idx][contiSolventEqIdx];
                    R_sum[ contiSolventEqIdx ] += R2;
                    maxCoeff[ contiSolventEqIdx ] = std::max( maxCoeff[ contiSolventEqIdx ], std::abs( R2 ) / pvValue );
                }
                if (has_polymer_ ) {
                    B_avg[ contiPolymerEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][contiPolymerEqIdx];
                    R_sum[ contiPolymerEqIdx ] += R2;
                    maxCoeff[ contiPolymerEqIdx ] = std::max( maxCoeff[ contiPolymerEqIdx ], std::abs( R2 ) / pvValue );
                }

            }

            // compute local average in terms of global number of elements
            const int bSize = B_avg.size();
            for ( int i = 0; i<bSize; ++i )
            {
                B_avg[ i ] /= Scalar( global_nc_ );
            }

            // TODO: we remove the maxNormWell for now because the convergence of wells are on a individual well basis.
            // Anyway, we need to provide some infromation to help debug the well iteration process.


            // compute global sum and max of quantities
            const double pvSum = convergenceReduction(grid_.comm(), pvSumLocal,
                                                      R_sum, maxCoeff, B_avg);

            Vector CNV(numComp);
            Vector mass_balance_residual(numComp);

            bool converged_MB = true;
            bool converged_CNV = true;
            // Finish computation
            for ( int compIdx = 0; compIdx < numComp; ++compIdx )
            {
                CNV[compIdx]                    = B_avg[compIdx] * dt * maxCoeff[compIdx];
                mass_balance_residual[compIdx]  = std::abs(B_avg[compIdx]*R_sum[compIdx]) * dt / pvSum;
                converged_MB                = converged_MB && (mass_balance_residual[compIdx] < tol_mb);
                converged_CNV               = converged_CNV && (CNV[compIdx] < tol_cnv);

                residual_norms.push_back(CNV[compIdx]);
            }

            const bool converged_Well = wellModel().getWellConvergence(ebosSimulator_, B_avg);

            bool converged = converged_MB && converged_Well;

            // do not care about the cell based residual in the last two Newton
            // iterations
            if (iteration < param_.max_strict_iter_)
                converged = converged && converged_CNV;

            if ( terminal_output_ )
            {
                // Only rank 0 does print to std::cout
                if (iteration == 0) {
                    std::string msg = "Iter";

                    std::vector< std::string > key( numComp );
                    for (int phaseIdx = 0; phaseIdx < numPhases(); ++phaseIdx) {
                        const std::string& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));
                        key[ phaseIdx ] = std::toupper( phaseName.front() );
                    }
                    if (has_solvent_) {
                        key[ solventSaturationIdx ] = "S";
                    }

                    if (has_polymer_) {
                        key[ polymerConcentrationIdx ] = "P";
                    }

                    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                        msg += "    MB(" + key[ compIdx ] + ")  ";
                    }
                    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                        msg += "    CNV(" + key[ compIdx ] + ") ";
                    }
                    OpmLog::note(msg);
                }
                std::ostringstream ss;
                const std::streamsize oprec = ss.precision(3);
                const std::ios::fmtflags oflags = ss.setf(std::ios::scientific);
                ss << std::setw(4) << iteration;
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    ss << std::setw(11) << mass_balance_residual[compIdx];
                }
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    ss << std::setw(11) << CNV[compIdx];
                }
                ss.precision(oprec);
                ss.flags(oflags);
                OpmLog::note(ss.str());
            }

            for (int phaseIdx = 0; phaseIdx < numPhases(); ++phaseIdx) {
                const auto& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));

                if (std::isnan(mass_balance_residual[phaseIdx])
                    || std::isnan(CNV[phaseIdx])) {
                    OPM_THROW(Opm::NumericalProblem, "NaN residual for phase " << phaseName);
                }
                if (mass_balance_residual[phaseIdx] > maxResidualAllowed()
                    || CNV[phaseIdx] > maxResidualAllowed()) {
                    OPM_THROW(Opm::NumericalProblem, "Too large residual for phase " << phaseName);
                }
            }

            return converged;
        }


        /// The number of active fluid phases in the model.
        int numPhases() const
        {
            return phaseUsage_.num_phases;
        }

        int numComponents() const
        {
            if (numPhases() == 2) {
                return 2;
            }
            int numComp = FluidSystem::numComponents;
            if (has_solvent_)
                numComp ++;

            if (has_polymer_)
                numComp ++;

            return numComp;
        }

        /// Wrapper required due to not following generic API
        template<class T>
        std::vector<std::vector<double> >
        computeFluidInPlace(const T&, const std::vector<int>& fipnum) const
        {
            return computeFluidInPlace(fipnum);
        }

        std::vector<std::vector<double> >
        computeFluidInPlace(const std::vector<int>& fipnum) const
        {
            const auto& comm = grid_.comm();
            const auto& gridView = ebosSimulator().gridView();
            const int nc = gridView.size(/*codim=*/0);
            const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
            int ntFip = *std::max_element(fipnum.begin(), fipnum.end());
            ntFip = comm.max(ntFip);

            std::vector<double> tpv(ntFip, 0.0);
            std::vector<double> hcpv(ntFip, 0.0);

            std::vector<std::vector<double> > regionValues(ntFip, std::vector<double>(FIPDataType::fipValues,0.0));

            for (int i = 0; i<FIPDataType::fipValues; i++) {
                fip_.fip[i].resize(nc,0.0);
            }

            ElementContext elemCtx(ebosSimulator_);
            const auto& elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();
            for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
                 elemIt != elemEndIt;
                 ++elemIt)
            {
                elemCtx.updatePrimaryStencil(*elemIt);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                const unsigned cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();

                const int regionIdx = fipnum[cellIdx] - 1;
                if (regionIdx < 0) {
                    // the given cell is not attributed to any region
                    continue;
                }

                // calculate the pore volume of the current cell. Note that the porosity
                // returned by the intensive quantities is defined as the ratio of pore
                // space to total cell volume and includes all pressure dependent (->
                // rock compressibility) and static modifiers (MULTPV, MULTREGP, NTG,
                // PORV, MINPV and friends). Also note that because of this, the porosity
                // returned by the intensive quantities can be outside of the physical
                // range [0, 1] in pathetic cases.
                const double pv =
                    ebosSimulator_.model().dofTotalVolume(cellIdx)
                    * intQuants.porosity().value();

                for (int phase = 0; phase < maxnp; ++phase) {
                    const double b = fs.invB(flowPhaseToEbosPhaseIdx(phase)).value();
                    const double s = fs.saturation(flowPhaseToEbosPhaseIdx(phase)).value();

                    fip_.fip[phase][cellIdx] = b * s * pv;

                    if (active_[ phase ]) {
                        regionValues[regionIdx][phase] += fip_.fip[phase][cellIdx];
                    }
                }

                if (active_[ Oil ] && active_[ Gas ]) {
                    // Account for gas dissolved in oil and vaporized oil
                    fip_.fip[FIPDataType::FIP_DISSOLVED_GAS][cellIdx] = fs.Rs().value() * fip_.fip[FIPDataType::FIP_LIQUID][cellIdx];
                    fip_.fip[FIPDataType::FIP_VAPORIZED_OIL][cellIdx] = fs.Rv().value() * fip_.fip[FIPDataType::FIP_VAPOUR][cellIdx];

                    regionValues[regionIdx][FIPData::FIP_DISSOLVED_GAS] += fip_.fip[FIPData::FIP_DISSOLVED_GAS][cellIdx];
                    regionValues[regionIdx][FIPData::FIP_VAPORIZED_OIL] += fip_.fip[FIPData::FIP_VAPORIZED_OIL][cellIdx];
                }

                const double hydrocarbon = fs.saturation(FluidSystem::oilPhaseIdx).value() + fs.saturation(FluidSystem::gasPhaseIdx).value();

                tpv[regionIdx] += pv;
                hcpv[regionIdx] += pv * hydrocarbon;
            }

            // sum tpv (-> total pore volume of the regions) and hcpv (-> pore volume of the
            // the regions that is occupied by hydrocarbons)
            comm.sum(tpv.data(), tpv.size());
            comm.sum(hcpv.data(), hcpv.size());

            for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
                 elemIt != elemEndIt;
                 ++elemIt)
            {
                const auto& elem = *elemIt;

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                unsigned cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const int regionIdx = fipnum[cellIdx] - 1;
                if (regionIdx < 0) {
                    // the cell is not attributed to any region. ignore it!
                    continue;
                }

                const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();

                // calculate the pore volume of the current cell. Note that the
                // porosity returned by the intensive quantities is defined as the
                // ratio of pore space to total cell volume and includes all pressure
                // dependent (-> rock compressibility) and static modifiers (MULTPV,
                // MULTREGP, NTG, PORV, MINPV and friends). Also note that because of
                // this, the porosity returned by the intensive quantities can be
                // outside of the physical range [0, 1] in pathetic cases.
                const double pv =
                    ebosSimulator_.model().dofTotalVolume(cellIdx)
                    * intQuants.porosity().value();

                fip_.fip[FIPDataType::FIP_PV][cellIdx] = pv;
                const double hydrocarbon = fs.saturation(FluidSystem::oilPhaseIdx).value() + fs.saturation(FluidSystem::gasPhaseIdx).value();

                //Compute hydrocarbon pore volume weighted average pressure.
                //If we have no hydrocarbon in region, use pore volume weighted average pressure instead
                if (hcpv[regionIdx] > 1e-10) {
                    fip_.fip[FIPDataType::FIP_WEIGHTED_PRESSURE][cellIdx] = pv * fs.pressure(FluidSystem::oilPhaseIdx).value() * hydrocarbon / hcpv[regionIdx];
                } else {
                    fip_.fip[FIPDataType::FIP_WEIGHTED_PRESSURE][cellIdx] = pv * fs.pressure(FluidSystem::oilPhaseIdx).value() / tpv[regionIdx];
                }

                regionValues[regionIdx][FIPDataType::FIP_PV] += fip_.fip[FIPDataType::FIP_PV][cellIdx];
                regionValues[regionIdx][FIPDataType::FIP_WEIGHTED_PRESSURE] += fip_.fip[FIPDataType::FIP_WEIGHTED_PRESSURE][cellIdx];
            }

            // sum the results over all processes
            for(int regionIdx=0; regionIdx < ntFip; ++regionIdx) {
                comm.sum(regionValues[regionIdx].data(), regionValues[regionIdx].size());
            }

            return regionValues;
        }

        SimulationDataContainer getSimulatorData ( const SimulationDataContainer& localState) const
        {
            typedef std::vector<double> VectorType;

            const auto& ebosModel = ebosSimulator().model();
            const auto& phaseUsage = phaseUsage_;

            // extract everything which can possibly be written to disk
            const int numCells   = ebosModel.numGridDof();
            const int num_phases = numPhases();

            SimulationDataContainer simData( numCells, 0, num_phases );

            //Get shorthands for water, oil, gas
            const int aqua_active = phaseUsage.phase_used[Opm::PhaseUsage::Aqua];
            const int liquid_active = phaseUsage.phase_used[Opm::PhaseUsage::Liquid];
            const int vapour_active = phaseUsage.phase_used[Opm::PhaseUsage::Vapour];

            const int aqua_pos   = phaseUsage.phase_pos[ Opm::PhaseUsage::Aqua ];
            const int liquid_pos = phaseUsage.phase_pos[ Opm::PhaseUsage::Liquid ];
            const int vapour_pos = phaseUsage.phase_pos[ Opm::PhaseUsage::Vapour ];

            VectorType zero;

            VectorType& pressureOil = simData.pressure();
            VectorType& temperature = simData.temperature();
            VectorType& saturation = simData.saturation();

            // WATER
            if( aqua_active ) {
                simData.registerCellData( "1OVERBW", 1 );
                simData.registerCellData( "WAT_DEN", 1 );
                simData.registerCellData( "WAT_VISC", 1 );
                simData.registerCellData( "WATKR", 1 );
            }

            VectorType& bWater   = aqua_active ? simData.getCellData( "1OVERBW" ) : zero;
            VectorType& rhoWater = aqua_active ? simData.getCellData( "WAT_DEN" ) : zero;
            VectorType& muWater  = aqua_active ? simData.getCellData( "WAT_VISC" ) : zero;
            VectorType& krWater  = aqua_active ? simData.getCellData( "WATKR" ) : zero;

            // OIL
            if( liquid_active ) {
                simData.registerCellData( "1OVERBO", 1 );
                simData.registerCellData( "OIL_DEN", 1 );
                simData.registerCellData( "OIL_VISC", 1 );
                simData.registerCellData( "OILKR", 1 );
            }

            VectorType& bOil   = liquid_active ? simData.getCellData( "1OVERBO" ) : zero;
            VectorType& rhoOil = liquid_active ? simData.getCellData( "OIL_DEN" ) : zero;
            VectorType& muOil  = liquid_active ? simData.getCellData( "OIL_VISC" ) : zero;
            VectorType& krOil  = liquid_active ? simData.getCellData( "OILKR" ) : zero;

            // GAS
            if( vapour_active ) {
                simData.registerCellData( "1OVERBG", 1 );
                simData.registerCellData( "GAS_DEN", 1 );
                simData.registerCellData( "GAS_VISC", 1 );
                simData.registerCellData( "GASKR", 1 );
            }

            VectorType& bGas   = vapour_active ? simData.getCellData( "1OVERBG" ) : zero;
            VectorType& rhoGas = vapour_active ? simData.getCellData( "GAS_DEN" ) : zero;
            VectorType& muGas  = vapour_active ? simData.getCellData( "GAS_VISC" ) : zero;
            VectorType& krGas  = vapour_active ? simData.getCellData( "GASKR" ) : zero;

            simData.registerCellData( BlackoilState::GASOILRATIO, 1 );
            simData.registerCellData( BlackoilState::RV, 1 );
            simData.registerCellData( "RSSAT", 1 );
            simData.registerCellData( "RVSAT", 1 );

            VectorType& Rs    = simData.getCellData( BlackoilState::GASOILRATIO );
            VectorType& Rv    = simData.getCellData( BlackoilState::RV );
            VectorType& RsSat = simData.getCellData( "RSSAT" );
            VectorType& RvSat = simData.getCellData( "RVSAT" );

            simData.registerCellData( "PBUB", 1 );
            simData.registerCellData( "PDEW", 1 );

            VectorType& Pb = simData.getCellData( "PBUB" );
            VectorType& Pd = simData.getCellData( "PDEW" );

            simData.registerCellData( "SOMAX", 1 );
            VectorType& somax = simData.getCellData( "SOMAX" );

            // Two components for hysteresis parameters
            // pcSwMdc/krnSwMdc, one for oil-water and one for gas-oil
            simData.registerCellData( "PCSWMDC_GO", 1 );
            simData.registerCellData( "KRNSWMDC_GO", 1 );

            simData.registerCellData( "PCSWMDC_OW", 1 );
            simData.registerCellData( "KRNSWMDC_OW", 1 );

            VectorType& pcSwMdc_go = simData.getCellData( "PCSWMDC_GO" );
            VectorType& krnSwMdc_go = simData.getCellData( "KRNSWMDC_GO" );

            VectorType& pcSwMdc_ow = simData.getCellData( "PCSWMDC_OW" );
            VectorType& krnSwMdc_ow = simData.getCellData( "KRNSWMDC_OW" );

            if (has_solvent_) {
                simData.registerCellData( "SSOL", 1 );
            }
            VectorType& ssol  = has_solvent_ ? simData.getCellData( "SSOL" ) : zero;

            if (has_polymer_) {
                simData.registerCellData( "POLYMER", 1 );
            }
            VectorType& cpolymer  = has_polymer_ ? simData.getCellData( "POLYMER" ) : zero;

            std::vector<int> failed_cells_pb;
            std::vector<int> failed_cells_pd;
            const auto& gridView = ebosSimulator().gridView();
            auto elemIt = gridView.template begin</*codim=*/ 0, Dune::Interior_Partition>();
            const auto& elemEndIt = gridView.template end</*codim=*/ 0, Dune::Interior_Partition>();
            ElementContext elemCtx(ebosSimulator());

            for (; elemIt != elemEndIt; ++elemIt) {
                const auto& elem = *elemIt;

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                const unsigned cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);

                const auto& fs = intQuants.fluidState();

                const int satIdx = cellIdx * num_phases;

                pressureOil[cellIdx] = fs.pressure(FluidSystem::oilPhaseIdx).value();

                temperature[cellIdx] = fs.temperature(FluidSystem::oilPhaseIdx).value();

                somax[cellIdx] = ebosSimulator().model().maxOilSaturation(cellIdx);

                const auto& matLawManager = ebosSimulator().problem().materialLawManager();
                if (matLawManager->enableHysteresis()) {
                    matLawManager->oilWaterHysteresisParams(
                            pcSwMdc_ow[cellIdx],
                            krnSwMdc_ow[cellIdx],
                            cellIdx);
                    matLawManager->gasOilHysteresisParams(
                            pcSwMdc_go[cellIdx],
                            krnSwMdc_go[cellIdx],
                            cellIdx);
                }

                if (aqua_active) {
                    saturation[ satIdx + aqua_pos ] = fs.saturation(FluidSystem::waterPhaseIdx).value();
                    bWater[cellIdx] = fs.invB(FluidSystem::waterPhaseIdx).value();
                    rhoWater[cellIdx] = fs.density(FluidSystem::waterPhaseIdx).value();
                    muWater[cellIdx] = fs.viscosity(FluidSystem::waterPhaseIdx).value();
                    krWater[cellIdx] = intQuants.relativePermeability(FluidSystem::waterPhaseIdx).value();
                }
                if (vapour_active) {
                    saturation[ satIdx + vapour_pos ]  = fs.saturation(FluidSystem::gasPhaseIdx).value();
                    bGas[cellIdx] = fs.invB(FluidSystem::gasPhaseIdx).value();
                    rhoGas[cellIdx] = fs.density(FluidSystem::gasPhaseIdx).value();
                    muGas[cellIdx] = fs.viscosity(FluidSystem::gasPhaseIdx).value();
                    krGas[cellIdx] = intQuants.relativePermeability(FluidSystem::gasPhaseIdx).value();
                    Rs[cellIdx] = fs.Rs().value();
                    Rv[cellIdx] = fs.Rv().value();
                    RsSat[cellIdx] = FluidSystem::saturatedDissolutionFactor(fs,
                                                                             FluidSystem::oilPhaseIdx,
                                                                             intQuants.pvtRegionIndex(),
                                                                             /*maxOilSaturation=*/1.0).value();
                    RvSat[cellIdx] = FluidSystem::saturatedDissolutionFactor(fs,
                                                                             FluidSystem::gasPhaseIdx,
                                                                             intQuants.pvtRegionIndex(),
                                                                             /*maxOilSaturation=*/1.0).value();
                    try {
                        Pb[cellIdx] = FluidSystem::bubblePointPressure(fs, intQuants.pvtRegionIndex()).value();
                    }
                    catch (const NumericalProblem& e) {
                        const auto globalIdx = ebosSimulator_.gridManager().grid().globalCell()[cellIdx];
                        failed_cells_pb.push_back(globalIdx);
                    }
                    try {
                        Pd[cellIdx] = FluidSystem::dewPointPressure(fs, intQuants.pvtRegionIndex()).value();
                    }
                    catch (const NumericalProblem& e) {
                        const auto globalIdx = ebosSimulator_.gridManager().grid().globalCell()[cellIdx];
                        failed_cells_pd.push_back(globalIdx);
                    }
                }
                if( liquid_active )
                {
                    saturation[ satIdx + liquid_pos ] = fs.saturation(FluidSystem::oilPhaseIdx).value();
                    bOil[cellIdx] = fs.invB(FluidSystem::oilPhaseIdx).value();
                    rhoOil[cellIdx] = fs.density(FluidSystem::oilPhaseIdx).value();
                    muOil[cellIdx] = fs.viscosity(FluidSystem::oilPhaseIdx).value();
                    krOil[cellIdx] = intQuants.relativePermeability(FluidSystem::oilPhaseIdx).value();
                }

                if (has_solvent_)
                {
                    ssol[cellIdx] = intQuants.solventSaturation().value();
                }

                if (has_polymer_)
                {
                    cpolymer[cellIdx] = intQuants.polymerConcentration().value();
                }

                // hack to make the intial output of rs and rv Ecl compatible.
                // For cells with swat == 1 Ecl outputs; rs = rsSat and rv=rvSat, in all but the initial step
                // where it outputs rs and rv values calculated by the initialization. To be compatible we overwrite
                // rs and rv with the values passed by the localState.
                // Volume factors, densities and viscosities need to be recalculated with the updated rs and rv values.
                if (ebosSimulator_.episodeIndex() < 0 && vapour_active && liquid_active ) {

                    Rs[cellIdx] = localState.getCellData( BlackoilState::GASOILRATIO )[cellIdx];
                    Rv[cellIdx] = localState.getCellData( BlackoilState::RV)[cellIdx];

                    // copy the fluidstate and set the new rs and rv values
                    auto fs_updated = fs;
                    auto rs_eval = fs_updated.Rs();
                    rs_eval.setValue( Rs[cellIdx] );
                    fs_updated.setRs(rs_eval);
                    auto rv_eval = fs_updated.Rv();
                    rv_eval.setValue( Rv[cellIdx] );
                    fs_updated.setRv(rv_eval);

                    //re-compute the volume factors, viscosities and densities.
                    rhoOil[cellIdx] = FluidSystem::density(fs_updated,
                                                           FluidSystem::oilPhaseIdx,
                                                           intQuants.pvtRegionIndex()).value();
                    rhoGas[cellIdx] = FluidSystem::density(fs_updated,
                                                           FluidSystem::gasPhaseIdx,
                                                           intQuants.pvtRegionIndex()).value();

                    bOil[cellIdx] = FluidSystem::inverseFormationVolumeFactor(fs_updated,
                                                           FluidSystem::oilPhaseIdx,
                                                           intQuants.pvtRegionIndex()).value();
                    bGas[cellIdx] = FluidSystem::inverseFormationVolumeFactor(fs_updated,
                                                           FluidSystem::gasPhaseIdx,
                                                           intQuants.pvtRegionIndex()).value();

                    muOil[cellIdx] = FluidSystem::viscosity(fs_updated,
                                                           FluidSystem::oilPhaseIdx,
                                                           intQuants.pvtRegionIndex()).value();
                    muGas[cellIdx] = FluidSystem::viscosity(fs_updated,
                                                           FluidSystem::gasPhaseIdx,
                                                           intQuants.pvtRegionIndex()).value();

                }
            }

            const size_t max_num_cells_faillog = 20;

            int pb_size = failed_cells_pb.size(), pd_size = failed_cells_pd.size();
            std::vector<int> displ_pb, displ_pd, recv_len_pb, recv_len_pd;
            const auto& comm = grid_.comm();

            if ( comm.rank() == 0 )
            {
                displ_pb.resize(comm.size()+1, 0);
                displ_pd.resize(comm.size()+1, 0);
                recv_len_pb.resize(comm.size());
                recv_len_pd.resize(comm.size());
            }

            comm.gather(&pb_size, recv_len_pb.data(), 1, 0);
            comm.gather(&pd_size, recv_len_pd.data(), 1, 0);
            std::partial_sum(recv_len_pb.begin(), recv_len_pb.end(), displ_pb.begin()+1);
            std::partial_sum(recv_len_pd.begin(), recv_len_pd.end(), displ_pd.begin()+1);
            std::vector<int> global_failed_cells_pb, global_failed_cells_pd;

            if ( comm.rank() == 0 )
            {
                global_failed_cells_pb.resize(displ_pb.back());
                global_failed_cells_pd.resize(displ_pd.back());
            }

            comm.gatherv(failed_cells_pb.data(), static_cast<int>(failed_cells_pb.size()),
                         global_failed_cells_pb.data(), recv_len_pb.data(),
                         displ_pb.data(), 0);
            comm.gatherv(failed_cells_pd.data(), static_cast<int>(failed_cells_pd.size()),
                         global_failed_cells_pd.data(),  recv_len_pd.data(),
                         displ_pd.data(), 0);
            std::sort(global_failed_cells_pb.begin(), global_failed_cells_pb.end());
            std::sort(global_failed_cells_pd.begin(), global_failed_cells_pd.end());

            if (global_failed_cells_pb.size() > 0) {
                std::stringstream errlog;
                errlog << "Finding the bubble point pressure failed for " << global_failed_cells_pb.size() << " cells [";
                errlog << global_failed_cells_pb[0];
                const size_t max_elems = std::min(max_num_cells_faillog, failed_cells_pb.size());
                for (size_t i = 1; i < max_elems; ++i) {
                    errlog << ", " << global_failed_cells_pb[i];
                }
                if (global_failed_cells_pb.size() > max_num_cells_faillog) {
                    errlog << ", ...";
                }
                errlog << "]";
                OpmLog::warning("Bubble point numerical problem", errlog.str());
            }
            if (global_failed_cells_pd.size() > 0) {
                std::stringstream errlog;
                errlog << "Finding the dew point pressure failed for " << global_failed_cells_pd.size() << " cells [";
                errlog << global_failed_cells_pd[0];
                const size_t max_elems = std::min(max_num_cells_faillog, global_failed_cells_pd.size());
                for (size_t i = 1; i < max_elems; ++i) {
                    errlog << ", " << global_failed_cells_pd[i];
                }
                if (global_failed_cells_pd.size() > max_num_cells_faillog) {
                    errlog << ", ...";
                }
                errlog << "]";
                OpmLog::warning("Dew point numerical problem", errlog.str());
            }

            return simData;
        }

        const FIPDataType& getFIPData() const {
            return fip_;
        }

        const Simulator& ebosSimulator() const
        { return ebosSimulator_; }

        /// return the statistics if the nonlinearIteration() method failed
        const SimulatorReport& failureReport() const
        { return failureReport_; }

    protected:
        const ISTLSolverType& istlSolver() const
        {
            assert( istlSolver_ );
            return *istlSolver_;
        }

        // ---------  Data members  ---------

        Simulator& ebosSimulator_;
        const Grid&            grid_;
        const ISTLSolverType*  istlSolver_;
        const PhaseUsage phaseUsage_;
        VFPProperties                   vfp_properties_;
        // For each canonical phase -> true if active
        const std::vector<bool>         active_;
        // Size = # active phases. Maps active -> canonical phase indices.
        const std::vector<int>          cells_;  // All grid cells
        const bool has_disgas_;
        const bool has_vapoil_;
        const bool has_solvent_;
        const bool has_polymer_;

        ModelParameters                 param_;
        SimulatorReport failureReport_;

        // Well Model
        BlackoilWellModel<TypeTag>& well_model_;

        /// \brief Whether we print something to std::cout
        bool terminal_output_;
        /// \brief The number of cells of the global grid.
        long int global_nc_;

        // rate converter between the surface volume rates and reservoir voidage rates
        RateConverterType& rate_converter_;

        std::vector<std::vector<double>> residual_norms_history_;
        double current_relaxation_;
        BVector dx_old_;
        mutable FIPDataType fip_;

    public:
        /// return the StandardWells object
        BlackoilWellModel<TypeTag>&
        wellModel() { return well_model_; }

        const BlackoilWellModel<TypeTag>&
        wellModel() const { return well_model_; }

        int numWells() const { return well_model_.numWells(); }

        /// return true if wells are available on this process
        bool localWellsActive() const { return well_model_.localWellsActive(); }




    public:
        int flowPhaseToEbosCompIdx( const int phaseIdx ) const
        {
            const auto& pu = phaseUsage_;
            if (active_[Water] && pu.phase_pos[Water] == phaseIdx)
                return Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            if (active_[Oil] && pu.phase_pos[Oil] == phaseIdx)
                return Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
            if (active_[Gas] && pu.phase_pos[Gas] == phaseIdx)
                return Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);

            // for other phases return the index
            return phaseIdx;
        }

        int flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
        {
            const auto& pu = phaseUsage_;
            if (active_[Water] && pu.phase_pos[Water] == phaseIdx)
                return FluidSystem::waterPhaseIdx;
            if (active_[Oil] && pu.phase_pos[Oil] == phaseIdx)
                return FluidSystem::oilPhaseIdx;
            if (active_[Gas] && pu.phase_pos[Gas] == phaseIdx)
                return FluidSystem::gasPhaseIdx;

            assert(phaseIdx < 3);
            // for other phases return the index
            return phaseIdx;
        }

    private:

        void updateRateConverter()
        {
            rate_converter_.defineState<ElementContext>(ebosSimulator_);
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
                                   const int iterationIdx)
        {
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
            // if the last time step failed we need to update the solution varables in ebos
            // and recalculate the Intesive Quantities.
            if ( timer.lastStepFailed() && iterationIdx == 0  ) {
                ebosSimulator_.model().solution( 0 /* timeIdx */ ) = ebosSimulator_.model().solution( 1 /* timeIdx */ );
                ebosSimulator_.model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);
            }

            ebosSimulator_.problem().beginIteration();
            ebosSimulator_.model().linearizer().linearize();
            ebosSimulator_.problem().endIteration();

            prevEpisodeIdx = ebosSimulator_.episodeIndex();

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
        std::vector<bool> wasSwitched_;
    };
} // namespace Opm

#endif // OPM_BLACKOILMODELBASE_IMPL_HEADER_INCLUDED
