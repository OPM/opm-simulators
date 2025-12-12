// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef TPSA_NEWTON_METHOD_HPP
#define TPSA_NEWTON_METHOD_HPP

#include <opm/common/Exceptions.hpp>

#include <opm/models/tpsa/tpsabaseproperties.hpp>
#include <opm/models/tpsa/tpsanewtonmethodparams.hpp>
#include <opm/models/utils/timer.hpp>
#include <opm/models/utils/timerguard.hh>

#include <opm/simulators/linalg/PropertyTree.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <unistd.h>


namespace Opm {

/*!
* \brief Newton method solving for generic TPSA model.
*
* Generates the Jacobian matrix, J(u^n) and residual vector, R(u^n), with a solution vector, u^n, at iteration n.
* Subsequently the linear system J(u^n)\delta u^n = -R(u^n) is solved to get u^{n+1} = u^n + \Delta u^n.
*/
template <class TypeTag>
class TpsaNewtonMethod
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Model = GetPropType<TypeTag, Properties::ModelTPSA>;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVectorTPSA>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVectorTPSA>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariablesTPSA>;
    using Constraints = GetPropType<TypeTag, Properties::Constraints>;
    using EqVector = GetPropType<TypeTag, Properties::EqVectorTPSA>;
    using Linearizer = GetPropType<TypeTag, Properties::LinearizerTPSA>;
    using LinearSolverBackend = GetPropType<TypeTag, Properties::LinearSolverBackendTPSA>;
    using ConvergenceWriter = GetPropType<TypeTag, Properties::NewtonConvergenceWriterTPSA>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapterTPSA>;

    using IstlMatrix = typename SparseMatrixAdapter::IstlMatrix;

public:
    /*!
    * \brief Constructor
    *
    * \param simulator Simulator object
    */
    explicit TpsaNewtonMethod(Simulator& simulator)
        : simulator_(simulator)
        , linearSolver_(simulator)
        , error_(1e100)
        , lastError_(1e100)
        , numIterations_(0)
        , numLinearizations_(0)
        , convergenceWriter_()
    {
        // Read runtime/default Newton parameters
        params_.read();
    }

    /*!
    * \brief Register all run-time parameters for the Newton method.
    */
    static void registerParameters()
    {
        LinearSolverBackend::registerParameters();
        TpsaNewtonMethodParams<Scalar>::registerParameters();
    }

    /*!
    * \brief Finialize the construction of the object.
    *
    * At this point, it can be assumed that all objects featured by the simulator have been allocated. (But not that
    * they have been fully initialized yet.)
    */
    void finishInit()
    { }

    /*!
    * \brief Run the Newton method.
    *
    * \returns Bool indicating if Newton converged
    */
    bool apply()
    {
        // Make sure all timers are prestine
        prePostProcessTimer_.halt();
        linearizeTimer_.halt();
        solveTimer_.halt();
        updateTimer_.halt();

        // Get vectors and linearizer
        SolutionVector& nextSolution = model().solution(/*historyIdx=*/0);
        SolutionVector currentSolution(nextSolution);
        GlobalEqVector solutionUpdate(nextSolution.size());

        Linearizer& linearizer = model().linearizer();

        TimerGuard prePostProcessTimerGuard(prePostProcessTimer_);

        // Tell the implementation that we begin solving
        prePostProcessTimer_.start();
        begin_(nextSolution);
        prePostProcessTimer_.stop();

        try {
            TimerGuard innerPrePostProcessTimerGuard(prePostProcessTimer_);
            TimerGuard linearizeTimerGuard(linearizeTimer_);
            TimerGuard updateTimerGuard(updateTimer_);
            TimerGuard solveTimerGuard(solveTimer_);

            // Execute the method as long as the implementation thinks that we should do another iteration
            while (proceed_()) {
                // Notify the implementation that we're about to start a new iteration
                prePostProcessTimer_.start();
                beginIteration_();
                prePostProcessTimer_.stop();

                // Make the current solution to the old one
                currentSolution = nextSolution;

                // Do the actual linearization
                linearizeTimer_.start();
                linearizeDomain_();
                linearizeAuxiliaryEquations_();
                linearizeTimer_.stop();

                // Get residual and Jacobian for convergence check and preparation of linear solver
                solveTimer_.start();
                auto& residual = linearizer.residual();
                const auto& jacobian = linearizer.jacobian();
                linearSolver_.prepare(jacobian, residual);
                linearSolver_.getResidual(residual);
                solveTimer_.stop();

                // The preSolve_() method usually computes the errors, but it can do something else in addition.
                // TODO: should its costs be counted to the linearization or to the update?
                updateTimer_.start();
                preSolve_(currentSolution, residual);
                updateTimer_.stop();

                // Check convergence criteria
                if (converged()) {
                    // Tell the implementation that we're done with this iteration
                    prePostProcessTimer_.start();
                    endIteration_(nextSolution, currentSolution);
                    prePostProcessTimer_.stop();

                    break;
                }

                // Solve A x = b, where b is the residual, A is its Jacobian and x is the update of the solution
                solveTimer_.start();
                solutionUpdate = 0.0;
                const bool conv = linearSolver_.solve(solutionUpdate);
                solveTimer_.stop();

                if (!conv) {
                    solveTimer_.stop();
                    if (verbose_()) {
                        std::cout << "TPSA: Linear solver did not converge!" << std::endl;
                    }

                    prePostProcessTimer_.start();
                    failed_();
                    prePostProcessTimer_.stop();

                    return false;
                }

                // Update the current solution with the delta
                updateTimer_.start();
                postSolve_(currentSolution, residual, solutionUpdate);
                update_(nextSolution, currentSolution, solutionUpdate, residual);
                updateTimer_.stop();

                // End of iteration calculations
                prePostProcessTimer_.start();
                endIteration_(nextSolution, currentSolution);
                prePostProcessTimer_.stop();
            }
        }
        catch (const Dune::Exception& e)
        {
            if (verbose_()) {
                std::cout << "TPSA: Newton method caught exception: \""
                          << e.what() << "\"\n" << std::flush;
            }

            prePostProcessTimer_.start();
            failed_();
            prePostProcessTimer_.stop();

            return false;
        }
        catch (const NumericalProblem& e)
        {
            if (verbose_()) {
                std::cout << "TPSA: Newton method caught exception: \""
                          << e.what() << "\"\n" << std::flush;
            }

            prePostProcessTimer_.start();
            failed_();
            prePostProcessTimer_.stop();

            return false;
        }

        // Tell the implementation that we're done
        prePostProcessTimer_.start();
        end_();
        prePostProcessTimer_.stop();

        // print the timing summary of the time step
        if (verbose_()) {
            Scalar elapsedTot =
                linearizeTimer_.realTimeElapsed() +
                solveTimer_.realTimeElapsed() +
                updateTimer_.realTimeElapsed();
            const auto default_precision{std::cout.precision()};
            std::cout << std::setprecision(2)
                      << "TPSA: "
                      << "Newton iter = " << numIterations() << " (error="
                      << error_ << ") | "
                      << "linearization = "
                      << linearizeTimer_.realTimeElapsed() << "s ("
                      << 100 * linearizeTimer_.realTimeElapsed() / elapsedTot << "%) | "
                      << "solve = "
                      << solveTimer_.realTimeElapsed() << "s ("
                      << 100 * solveTimer_.realTimeElapsed() / elapsedTot << "%) | "
                      << "update = "
                      << updateTimer_.realTimeElapsed() << "s ("
                      << 100 * updateTimer_.realTimeElapsed() / elapsedTot << "%)"
                      << "\n" << std::flush;
            std::cout << std::setprecision(default_precision); // restore default output width
        }

        // if we're not converged, tell the implementation that we've failed
        if (!converged()) {
            prePostProcessTimer_.start();
            failed_();
            std::cout << "TPSA: Newton iterations did not converge!" << std::endl;
            prePostProcessTimer_.stop();
            return false;
        }

        // if we converged, tell the implementation that we've succeeded
        prePostProcessTimer_.start();
        succeeded_();
        prePostProcessTimer_.stop();

        return true;
    }

    /*!
    * \brief Set the index of current iteration.
    *
    * \param value New iteration index
    *
    * Normally this does not need to be called, but if the non-linear solver is implemented externally, it needs to be
    * set in order for the model to do the Right Thing (TM) while linearizing.
    */
    void setIterationIndex(int value)
    {
        numIterations_ = value;
    }

    /*!
    * \brief Set the current tolerance at which the Newton method considers itself to be converged.
    *
    * \param value New tolerance value
    */
    void setTolerance(Scalar value)
    {
        params_.tolerance_ = value;
    }

    /*!
    * \brief Returns true if the error of the solution is below the tolerance.
    *
    * \returns Bool indicating if convergence has been achived
    */
    bool converged() const
    { return error_ <= tolerance(); }

    /*!
    * \brief Returns a reference to the object describing the current physical problem.
    *
    * \returns Reference to problem object
    */
    Problem& problem()
    { return simulator_.problem(); }

    /*!
    * \brief Returns a reference to the object describing the current physical problem.
    *
    * \returns Reference to problem object
    */
    const Problem& problem() const
    { return simulator_.problem(); }

    /*!
    * \brief Returns a reference to the geomechanics model.
    *
    * \returns Reference to geomechanics model object
    */
    Model& model()
    { return simulator_.problem().geoMechModel(); }

    /*!
    * \brief Returns a reference to the geomechanics model.
    *
    * \returns Reference to geomechanics model object
    */
    const Model& model() const
    { return simulator_.problem().geoMechModel(); }

    /*!
    * \brief Returns the linear solver backend object for external use.
    *
    * \returns Reference to linear solver object
    */
    LinearSolverBackend& linearSolver()
    { return linearSolver_; }

    /*!
    * \brief Returns the linear solver backend object for external use.
    *
    * \returns Reference to linear solver object
    */
    const LinearSolverBackend& linearSolver() const
    { return linearSolver_; }

    /*!
    * \brief Returns the number of iterations done since the Newton method was invoked.
    *
    * \returns Number of Newton iteratinos
    */
    int numIterations() const
    { return numIterations_; }

    /*!
    * \brief Returns the number of linearizations that has done since the Newton method was invoked.
    *
    * \returns Number of linearizations
    */
    int numLinearizations() const
    { return numLinearizations_; }

    /*!
    * \brief Return the current tolerance at which the Newton method considers itself to be converged.
    *
    * \returns Tolerance for Newton error
    */
    Scalar tolerance() const
    { return params_.tolerance_; }

    /*!
    * \brief Returns minimum number of Newton iterations used
    *
    * \returns Minimum number of Newton iterations
    */
    Scalar minIterations() const
    { return params_.minIterations_; }

    /*!
    * \brief Return post-process timer
    *
    * \returns Reference to post-process timer object
    */
    const Timer& prePostProcessTimer() const
    { return prePostProcessTimer_; }

    /*!
    * \brief Return linearization timer
    *
    * \returns Reference to linearization timer object
    */
    const Timer& linearizeTimer() const
    { return linearizeTimer_; }

    /*!
    * \brief Return linear solver timer
    *
    * \returns Reference to linear solver timer object
    */
    const Timer& solveTimer() const
    { return solveTimer_; }

    /*!
    * \brief Return solution update timer
    *
    * \returns Reference to solution update timer object
    */
    const Timer& updateTimer() const
    { return updateTimer_; }

protected:
    /*!
    * \brief Returns true if the Newton method ought to be chatty.
    *
    * \returns Bool indicating if the Newton procedure should be verbose
    */
    bool verbose_() const
    { return params_.verbose_ && (simulator_.gridView().comm().rank() == 0); }

    /*!
    * \brief Called before the Newton method is applied to an non-linear system of equations.
    *
    * \param nextSolution The initial solution vector
    */
    void begin_(const SolutionVector& /*nextSolution*/)
    {
        numIterations_ = 0;
        numLinearizations_ = 0;
        error_ = 1e100;

        if (params_.writeConvergence_) {
            convergenceWriter_.beginTimeStep();
        }
    }

    /*!
    * \brief Indicates the beginning of a Newton iteration.
    */
    void beginIteration_()
    {
        // start with a clean message stream
        const auto& comm = simulator_.gridView().comm();
        bool succeeded = true;
        succeeded = comm.min(succeeded);

        if (!succeeded) {
            throw NumericalProblem("TPSA: Pre-processing of the problem failed!");
        }

        lastError_ = error_;
    }

    /*!
    * \brief Linearize the global non-linear system of equations associated with the spatial domain.
    */
    void linearizeDomain_()
    {
        model().linearizer().linearizeDomain();
        ++numLinearizations_;
    }

    /*!
    * \brief Linearize auxillary equations
    */
    void linearizeAuxiliaryEquations_()
    {
        model().linearizer().linearizeAuxiliaryEquations();
        model().linearizer().finalize();
    }

    /*!
    * \brief Compute error before a Newton iteration
    *
    * \param currentSolution Current solution vector
    * \param currentResidual Current residual vector
    */
    void preSolve_(const SolutionVector& /*currentSolution*/,
                   const GlobalEqVector& currentResidual)
    {
        const auto& constraintsMap = model().linearizer().constraintsMap();
        lastError_ = error_;
        Scalar newtonMaxError = params_.maxError_;

        // Calculate the error as the maximum weighted tolerance of the solution's residual
        error_ = 0;
        const auto& elemMapper = simulator_.model().elementMapper();
        for (const auto& elem : elements(simulator_.gridView(), Dune::Partitions::interior)) {
            unsigned dofIdx = elemMapper.index(elem);

            // Do not consider auxiliary DOFs for the error
            if (dofIdx >= model().numGridDof() || model().dofTotalVolume(dofIdx) <= 0.0) {
                continue;
            }

            // Also do not consider DOFs which are constraint
            if (enableConstraints_()) {
                if (constraintsMap.count(dofIdx) > 0) {
                    continue;
                }
            }

            const auto& r = currentResidual[dofIdx];
            for (unsigned eqIdx = 0; eqIdx < r.size(); ++eqIdx) {
                error_ = max(std::abs(r[eqIdx] * model().eqWeight(dofIdx, eqIdx)), error_);
            }
        }

        // Take the other processes into account
        error_ = simulator_.gridView().comm().max(error_);

        // Make sure that the error never grows beyond the maximum allowed one
        if (error_ > newtonMaxError) {
            throw NumericalProblem("TPSA: Newton error " + std::to_string(double(error_)) +
                                   " is larger than maximum allowed error of " +
                                   std::to_string(double(newtonMaxError)));
        }
    }

    /*!
    * \brief Update the error of the solution given the previous iteration.
    *
    * \param currentSolution The solution at the beginning the current iteration
    * \param currentResidual The residual (i.e., right-hand-side) of the current iteration's solution.
    * \param solutionUpdate The difference between the current and the next solution
    *
    * For our purposes, the error of a solution is defined as the maximum of the weighted residual of a given solution.
    */
    void postSolve_(const SolutionVector& /*currentSolution*/,
                    const GlobalEqVector& /*currentResidual*/,
                    GlobalEqVector& /*solutionUpdate*/)
    { }

    /*!
    * \brief Update the current solution with a delta vector.
    *
    * \param nextSolution The solution vector after the current iteration
    * \param currentSolution The solution vector after the last iteration
    * \param solutionUpdate The delta vector as calculated by solving the linear system of equations
    * \param currentResidual The residual vector of the current Newton-Raphson iteraton
    *
    * Different update strategies, such as chopped updates can be implemented by overriding this method. The default
    * behavior is use the standard Newton-Raphson update strategy, i.e. u^{n+1} = u^n - \Delta u^n
    */
    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual)
    {
        const auto& constraintsMap = model().linearizer().constraintsMap();

        // first, write out the current solution to make convergence
        // analysis possible
        writeConvergence_(currentSolution, solutionUpdate);

        // make sure not to swallow non-finite values at this point
        if (!std::isfinite(solutionUpdate.one_norm())) {
            throw NumericalProblem("TPSA: Non-finite update in Newton!");
        }

        std::size_t numGridDof = model().numGridDof();
        for (unsigned dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
            if (enableConstraints_()) {
                if (constraintsMap.count(dofIdx) > 0) {
                    const auto& constraints = constraintsMap.at(dofIdx);
                    updateConstraintDof_(dofIdx,
                                         nextSolution[dofIdx],
                                         constraints);
                }
                else {
                    updatePrimaryVariables_(dofIdx,
                                            nextSolution[dofIdx],
                                            currentSolution[dofIdx],
                                            solutionUpdate[dofIdx],
                                            currentResidual[dofIdx]);
                }
            }
            else {
                updatePrimaryVariables_(dofIdx,
                                        nextSolution[dofIdx],
                                        currentSolution[dofIdx],
                                        solutionUpdate[dofIdx],
                                        currentResidual[dofIdx]);
            }
        }

        // update the DOFs of the auxiliary equations
        // std::size_t numDof = model().numTotalDof();
        // for (std::size_t dofIdx = numGridDof; dofIdx < numDof; ++dofIdx) {
        //     nextSolution[dofIdx] = currentSolution[dofIdx];
        //     nextSolution[dofIdx] -= solutionUpdate[dofIdx];
        // }

        // Update material state
        model().updateMaterialState(/*timeIdx=*/0);
    }

    /*!
    * \brief Update the primary variables for a degree of freedom which is constraint.
    *
    * \param dofIdx Degree-of-freedom index
    * \param nextValue The solution vector after the current iteration
    * \param constraints Constraints for solution
    *
    */
    void updateConstraintDof_(unsigned /*dofIdx*/,
                              PrimaryVariables& /*nextValue*/,
                              const Constraints& /*constraints*/)
    { }

    /*!
    * \brief Update a single primary variables object.
    *
    * \param nextValue The solution vector after the current iteration
    * \param currentValue The solution vector after the last iteration
    * \param update The delta vector as calculated by solving the linear system of equations
    * \param currentResidual The residual vector of the current Newton-Raphson iteraton
    */
    void updatePrimaryVariables_(unsigned /*dofIdx*/,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& /*currentResidual*/)
    {
        nextValue = currentValue;
        nextValue -= update;
    }

    /*!
    * \brief Write the convergence behaviour of the newton method to disk.
    *
    * \param currentSolution The solution vector after the last iteration
    * \param solutionUpdate The delta vector as calculated by solving the linear system of equations
    *
    * This method is called as part of the update proceedure.
    */
    void writeConvergence_(const SolutionVector& currentSolution,
                           const GlobalEqVector& solutionUpdate)
    {
        if (params_.writeConvergence_) {
            convergenceWriter_.beginIteration();
            convergenceWriter_.writeFields(currentSolution, solutionUpdate);
            convergenceWriter_.endIteration();
        }
    }

    /*!
    * \brief Indicates that one Newton iteration was finished.
    *
    * \param nextSolution The solution after the current Newton iteration
    * \param currentSolution The solution at the beginning of the current Newton iteration
    */
    void endIteration_(const SolutionVector& /*nextSolution*/,
                       const SolutionVector& /*currentSolution*/)
    {
        ++numIterations_;

        const auto& comm = simulator_.gridView().comm();
        bool succeeded = true;
        succeeded = comm.min(succeeded);

        if (!succeeded) {
            throw NumericalProblem("TPSA: Post-processing of the problem failed!");
        }

        // if (verbose_()) {
        //     std::cout << "TPSA: End Newton iteration " << numIterations_ << ""
        //               << " with error = " << error_
        //               << std::endl;
        // }
    }

    /*!
    * \brief Returns true iff another Newton iteration should be done.
    *
    * \returns Bool indicating if Newton iterations should continue
    */
    bool proceed_() const
    {
        if (numIterations() < params_.minIterations_) {
            return true;
        }
        else if (converged()) {
            // we are below the specified tolerance, so we don't have to
            // do more iterations
            return false;
        }
        else if (numIterations() >= params_.maxIterations_) {
            // we have exceeded the allowed number of steps.  If the
            // error was reduced by a factor of at least 4,
            // in the last iterations we proceed even if we are above
            // the maximum number of steps
            return error_ * 4.0 < lastError_;
        }

        return true;
    }

    /*!
    * \brief Indicates that we're done solving the non-linear system of equations.
    */
    void end_()
    {
        if (params_.writeConvergence_) {
            convergenceWriter_.endTimeStep();
        }
    }

    /*!
    * \brief Called if the Newton method broke down.
    *
    * \note This method is called _after_ end_()
    */
    void failed_()
    { numIterations_ = params_.targetIterations_ * 2; }

    /*!
    * \brief Called if the Newton method was successful.
    *
    * \note This method is called _after_ end_()
    */
    void succeeded_()
    { }

    /*!
    * \brief Check if constraints are enabled
    *
    * \returns Bool indicating if constraints are enabled
    */
    static bool enableConstraints_()
    { return getPropValue<TypeTag, Properties::EnableConstraintsTPSA>(); }

    Simulator& simulator_;
    LinearSolverBackend linearSolver_;

    Timer prePostProcessTimer_;
    Timer linearizeTimer_;
    Timer solveTimer_;
    Timer updateTimer_;

    Scalar error_;
    Scalar lastError_;
    TpsaNewtonMethodParams<Scalar> params_;

    int numIterations_;
    int numLinearizations_;

    ConvergenceWriter convergenceWriter_;
};  // class TpsaNewtonMethod

}  // namespace Opm

#endif
