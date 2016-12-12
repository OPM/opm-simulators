/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2015 Statoil ASA.

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

#ifndef OPM_NONLINEARSOLVER_HEADER_INCLUDED
#define OPM_NONLINEARSOLVER_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/simulator/SimulatorTimerInterface.hpp>
#include <opm/autodiff/DuneMatrix.hpp>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <memory>

namespace Opm {


    /// A nonlinear solver class suitable for general fully-implicit models,
    /// as well as pressure, transport and sequential models.
    template <class PhysicalModel>
    class NonlinearSolver
    {
    public:
        // ---------  Types and enums  ---------
        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;

        // Available relaxation scheme types.
        enum RelaxType { DAMPEN, SOR };

        // Solver parameters controlling nonlinear process.
        struct SolverParameters
        {
            enum RelaxType relax_type_;
            double         relax_max_;
            double         relax_increment_;
            double         relax_rel_tol_;
            int            max_iter_; // max nonlinear iterations
            int            min_iter_; // min nonlinear iterations

            explicit SolverParameters( const parameter::ParameterGroup& param );
            SolverParameters();

            void reset();
        };

        // Forwarding types from PhysicalModel.
        typedef typename PhysicalModel::ReservoirState ReservoirState;
        typedef typename PhysicalModel::WellState WellState;

        // ---------  Public methods  ---------

        /// Construct solver for a given model.
        ///
        /// The model is a std::unique_ptr because the object to which model points to is
        /// not allowed to be deleted as long as the NonlinearSolver object exists.
        ///
        /// \param[in]      param   parameters controlling nonlinear process
        /// \param[in, out] model   physical simulation model.
        explicit NonlinearSolver(const SolverParameters& param,
                              std::unique_ptr<PhysicalModel> model);

        /// Take a single forward step, after which the states will be modified
        /// according to the physical model.
        /// \param[in] timer                  simulation timer
        /// \param[in, out] reservoir_state     reservoir state variables
        /// \param[in, out] well_state          well state variables
        SimulatorReport
        step(const SimulatorTimerInterface& timer,
             ReservoirState& reservoir_state,
             WellState& well_state);

        /// Take a single forward step, after which the states will be modified
        /// according to the physical model. This version allows for the
        /// states passed as in/out arguments to be different from the initial
        /// states.
        /// \param[in] timer                  simulation timer
        /// \param[in] initial_reservoir_state  reservoir state variables at start of timestep
        /// \param[in] initial_well_state       well state variables at start of timestep
        /// \param[in, out] reservoir_state     reservoir state variables
        /// \param[in, out] well_state          well state variables
        /// \return                             number of linear iterations used
        SimulatorReport
        step(const SimulatorTimerInterface& timer,
             const ReservoirState& initial_reservoir_state,
             const WellState& initial_well_state,
             ReservoirState& reservoir_state,
             WellState& well_state);


        /// Number of linearizations used in all calls to step().
        int linearizations() const;

        /// Number of full nonlinear solver iterations used in all calls to step().
        int nonlinearIterations() const;

        /// Number of linear solver iterations used in all calls to step().
        int linearIterations() const;

        /// Number of well iterations used in all calls to step().
        int wellIterations() const;

        /// Number of nonlinear solver iterations used in the last call to step().
        int nonlinearIterationsLastStep() const;

        /// Number of linear solver iterations used in the last call to step().
        int linearIterationsLastStep() const;

        /// Number of well iterations used in all calls to step().
        int wellIterationsLastStep() const;

        /// Compute fluid in place.
        /// \param[in]    ReservoirState
        /// \param[in]    FIPNUM for active cells not global cells.
        /// \return fluid in place, number of fip regions, each region contains 5 values which are liquid, vapour, water, free gas and dissolved gas.
        std::vector<V>
        computeFluidInPlace(const ReservoirState& x,
                            const std::vector<int>& fipnum) const
        {
            return model_->computeFluidInPlace(x, fipnum);
        }

        std::vector<std::vector<double>>
        computeFluidInPlace(const std::vector<int>& fipnum) const
        {
            return model_->computeFluidInPlace(fipnum);
        }


        /// Reference to physical model.
        const PhysicalModel& model() const;

        /// Mutable reference to physical model.
        PhysicalModel& model();

        /// Detect oscillation or stagnation in a given residual history.
        void detectOscillations(const std::vector<std::vector<double>>& residual_history,
                                const int it, bool& oscillate, bool& stagnate) const;

        /// Apply a stabilization to dx, depending on dxOld and relaxation parameters.
        void stabilizeNonlinearUpdate(V& dx, V& dxOld, const double omega) const;

        /// Apply a stabilization to dx, depending on dxOld and relaxation parameters.
        /// Implemention for Dune block vectors.
        template <class BVector>
        void stabilizeNonlinearUpdate(BVector& dx, BVector& dxOld, const double omega) const;

        /// The greatest relaxation factor (i.e. smallest factor) allowed.
        double relaxMax() const          { return param_.relax_max_; }

        /// The step-change size for the relaxation factor.
        double relaxIncrement() const    { return param_.relax_increment_; }

        /// The relaxation type (DAMPEN or SOR).
        enum RelaxType relaxType() const { return param_.relax_type_; }

        /// The relaxation relative tolerance.
        double relaxRelTol() const       { return param_.relax_rel_tol_; }

        /// The maximum number of nonlinear iterations allowed.
        int maxIter() const           { return param_.max_iter_; }

        /// The minimum number of nonlinear iterations allowed.
        double minIter() const           { return param_.min_iter_; }

    private:
        // ---------  Data members  ---------
        SolverParameters param_;
        std::unique_ptr<PhysicalModel> model_;
        int linearizations_;
        int nonlinearIterations_;
        int linearIterations_;
        int wellIterations_;
        int nonlinearIterationsLast_;
        int linearIterationsLast_;
        int wellIterationsLast_;
    };
} // namespace Opm

#include "NonlinearSolver_impl.hpp"

#endif // OPM_NONLINEARSOLVER_HEADER_INCLUDED
