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
#include <opm/core/utility/parameters/ParameterGroup.hpp>
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
        /// \param[in] dt                time step size
        /// \param[in] reservoir_state   reservoir state variables
        /// \param[in] well_state        well state variables
        /// \return                      number of linear iterations used
        int
        step(const double dt,
             ReservoirState& reservoir_state,
             WellState& well_state);

        /// Number of nonlinear solver iterations used in all calls to step().
        unsigned int nonlinearIterations() const;

        /// Number of linear solver iterations used in all calls to step().
        unsigned int linearIterations() const;

        /// Number of nonlinear solver iterations used in the last call to step().
        unsigned int nonlinearIterationsLastStep() const;

        /// Number of linear solver iterations used in the last call to step().
        unsigned int linearIterationsLastStep() const;

        /// return reference to physical model
        const PhysicalModel& model() const;

    private:
        // ---------  Data members  ---------
        SolverParameters param_;
        std::unique_ptr<PhysicalModel> model_;
        unsigned int nonlinearIterations_;
        unsigned int linearIterations_;
        unsigned int nonlinearIterationsLast_;
        unsigned int linearIterationsLast_;

        // ---------  Private methods  ---------
        enum RelaxType relaxType() const { return param_.relax_type_; }
        double relaxMax() const          { return param_.relax_max_; }
        double relaxIncrement() const    { return param_.relax_increment_; }
        double relaxRelTol() const       { return param_.relax_rel_tol_; }
        double maxIter() const           { return param_.max_iter_; }
        double minIter() const           { return param_.min_iter_; }
        void detectOscillations(const std::vector<std::vector<double>>& residual_history,
                                const int it, const double relaxRelTol,
                                bool& oscillate, bool& stagnate) const;
        void stabilizeNonlinearUpdate(V& dx, V& dxOld, const double omega, const RelaxType relax_type) const;
    };
} // namespace Opm

#include "NonlinearSolver_impl.hpp"

#endif // OPM_NONLINEARSOLVER_HEADER_INCLUDED
