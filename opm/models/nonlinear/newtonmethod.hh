// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
/*!
 * \file
 * \copydoc Opm::NewtonMethod
 */
#ifndef EWOMS_NEWTON_METHOD_HH
#define EWOMS_NEWTON_METHOD_HH

#include <dune/istl/istlexception.hh>
#include <dune/common/classname.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <opm/common/Exceptions.hpp>

#include <opm/material/densead/Math.hpp>

#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <opm/models/nonlinear/newtonmethodparams.hpp>
#include <opm/models/nonlinear/newtonmethodproperties.hh>
#include <opm/models/nonlinear/nullconvergencewriter.hh>

#include <opm/models/utils/timer.hpp>
#include <opm/models/utils/timerguard.hh>

#include <opm/simulators/linalg/linalgproperties.hh>

#include <iostream>
#include <sstream>

#include <unistd.h>

namespace Opm {
// forward declaration of classes
template <class TypeTag>
class NewtonMethod;
}

namespace Opm {
// forward declaration of property tags
} // namespace Opm

namespace Opm::Properties {

namespace TTag {

//! The type tag on which the default properties for the Newton method
//! are attached
struct NewtonMethod {};

} // namespace TTag

// set default values for the properties
template<class TypeTag>
struct NewtonMethod<TypeTag, TTag::NewtonMethod> { using type = ::Opm::NewtonMethod<TypeTag>; };
template<class TypeTag>
struct NewtonConvergenceWriter<TypeTag, TTag::NewtonMethod> { using type = NullConvergenceWriter<TypeTag>; };

} // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup Newton
 * \brief The multi-dimensional Newton method.
 *
 * This class uses static polymorphism to allow implementations to
 * implement different update/convergence strategies.
 */
template <class TypeTag>
class NewtonMethod
{
    using Implementation = GetPropType<TypeTag, Properties::NewtonMethod>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Model = GetPropType<TypeTag, Properties::Model>;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Constraints = GetPropType<TypeTag, Properties::Constraints>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Linearizer = GetPropType<TypeTag, Properties::Linearizer>;
    using LinearSolverBackend = GetPropType<TypeTag, Properties::LinearSolverBackend>;
    using ConvergenceWriter = GetPropType<TypeTag, Properties::NewtonConvergenceWriter>;

    using Communicator = typename Dune::MPIHelper::MPICommunicator;
    using CollectiveCommunication = typename Dune::Communication<typename Dune::MPIHelper::MPICommunicator>;

public:
    NewtonMethod(Simulator& simulator)
        : simulator_(simulator)
        , endIterMsgStream_(std::ostringstream::out)
        , linearSolver_(simulator)
        , comm_(Dune::MPIHelper::getCommunicator())
        , convergenceWriter_(asImp_())
    {
        lastError_ = 1e100;
        error_ = 1e100;

        numIterations_ = 0;
        params_.read();
    }

    /*!
     * \brief Register all run-time parameters for the Newton method.
     */
    static void registerParameters()
    {
        LinearSolverBackend::registerParameters();
        NewtonMethodParams<Scalar>::registerParameters();
    }

    /*!
     * \brief Finialize the construction of the object.
     *
     * At this point, it can be assumed that all objects featured by the simulator have
     * been allocated. (But not that they have been fully initialized yet.)
     */
    void finishInit()
    { }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool converged() const
    { return error_ <= tolerance(); }

    /*!
     * \brief Returns a reference to the object describing the current physical problem.
     */
    Problem& problem()
    { return simulator_.problem(); }

    /*!
     * \brief Returns a reference to the object describing the current physical problem.
     */
    const Problem& problem() const
    { return simulator_.problem(); }

    /*!
     * \brief Returns a reference to the numeric model.
     */
    Model& model()
    { return simulator_.model(); }

    /*!
     * \brief Returns a reference to the numeric model.
     */
    const Model& model() const
    { return simulator_.model(); }

    /*!
     * \brief Returns the number of iterations done since the Newton method
     *        was invoked.
     */
    int numIterations() const
    { return numIterations_; }

    /*!
     * \brief Set the index of current iteration.
     *
     * Normally this does not need to be called, but if the non-linear solver is
     * implemented externally, it needs to be set in order for the model to do the Right
     * Thing (TM) while linearizing.
     */
    void setIterationIndex(int value)
    { numIterations_ = value; }

    /*!
     * \brief Return the current tolerance at which the Newton method considers itself to
     *        be converged.
     */
    Scalar tolerance() const
    { return params_.tolerance_; }

    /*!
     * \brief Set the current tolerance at which the Newton method considers itself to
     *        be converged.
     */
    void setTolerance(Scalar value)
    { params_.tolerance_ = value; }

    /*!
     * \brief Run the Newton method.
     *
     * The actual implementation can influence all the strategic
     * decisions via callbacks using static polymorphism.
     */
    bool apply()
    {
        // Clear the current line using an ansi escape
        // sequence.  For an explanation see
        // http://en.wikipedia.org/wiki/ANSI_escape_code
        const char *clearRemainingLine = "\n";
        if (isatty(fileno(stdout))) {
            static const char blubb[] = { 0x1b, '[', 'K', '\r', 0 };
            clearRemainingLine = blubb;
        }

        // make sure all timers are prestine
        prePostProcessTimer_.halt();
        linearizeTimer_.halt();
        solveTimer_.halt();
        updateTimer_.halt();

        SolutionVector& nextSolution = model().solution(/*historyIdx=*/0);
        SolutionVector currentSolution(nextSolution);
        GlobalEqVector solutionUpdate(nextSolution.size());

        Linearizer& linearizer = model().linearizer();

        TimerGuard prePostProcessTimerGuard(prePostProcessTimer_);

        // tell the implementation that we begin solving
        prePostProcessTimer_.start();
        asImp_().begin_(nextSolution);
        prePostProcessTimer_.stop();

        try {
            TimerGuard innerPrePostProcessTimerGuard(prePostProcessTimer_);
            TimerGuard linearizeTimerGuard(linearizeTimer_);
            TimerGuard updateTimerGuard(updateTimer_);
            TimerGuard solveTimerGuard(solveTimer_);

            // execute the method as long as the implementation thinks
            // that we should do another iteration
            while (asImp_().proceed_()) {
                // linearize the problem at the current solution

                // notify the implementation that we're about to start
                // a new iteration
                prePostProcessTimer_.start();
                asImp_().beginIteration_();
                prePostProcessTimer_.stop();

                // make the current solution to the old one
                currentSolution = nextSolution;

                if (asImp_().verbose_()) {
                    std::cout << "Linearize: r(x^k) = dS/dt + div F - q;   M = grad r"
                              << clearRemainingLine
                              << std::flush;
                }

                // do the actual linearization
                linearizeTimer_.start();
                asImp_().linearizeDomain_();
                asImp_().linearizeAuxiliaryEquations_();
                linearizeTimer_.stop();

                solveTimer_.start();
                auto& residual = linearizer.residual();
                const auto& jacobian = linearizer.jacobian();
                linearSolver_.prepare(jacobian, residual);
                linearSolver_.setResidual(residual);
                linearSolver_.getResidual(residual);
                solveTimer_.stop();

                // The preSolve_() method usually computes the errors, but it can do
                // something else in addition. TODO: should its costs be counted to
                // the linearization or to the update?
                updateTimer_.start();
                asImp_().preSolve_(currentSolution, residual);
                updateTimer_.stop();

                if (!asImp_().proceed_()) {
                    if (asImp_().verbose_() && isatty(fileno(stdout)))
                        std::cout << clearRemainingLine
                                  << std::flush;

                    // tell the implementation that we're done with this iteration
                    prePostProcessTimer_.start();
                    asImp_().endIteration_(nextSolution, currentSolution);
                    prePostProcessTimer_.stop();

                    break;
                }

                // solve the resulting linear equation system
                if (asImp_().verbose_()) {
                    std::cout << "Solve: M deltax^k = r"
                              << clearRemainingLine
                              << std::flush;
                }

                solveTimer_.start();
                // solve A x = b, where b is the residual, A is its Jacobian and x is the
                // update of the solution
                linearSolver_.setMatrix(jacobian);
                solutionUpdate = 0.0;
                bool converged = linearSolver_.solve(solutionUpdate);
                solveTimer_.stop();

                if (!converged) {
                    solveTimer_.stop();
                    if (asImp_().verbose_())
                        std::cout << "Newton: Linear solver did not converge\n" << std::flush;

                    prePostProcessTimer_.start();
                    asImp_().failed_();
                    prePostProcessTimer_.stop();

                    return false;
                }

                // update the solution
                if (asImp_().verbose_()) {
                    std::cout << "Update: x^(k+1) = x^k - deltax^k"
                              << clearRemainingLine
                              << std::flush;
                }

                // update the current solution (i.e. uOld) with the delta
                // (i.e. u). The result is stored in u
                updateTimer_.start();
                asImp_().postSolve_(currentSolution,
                                    residual,
                                    solutionUpdate);
                asImp_().update_(nextSolution, currentSolution, solutionUpdate, residual);
                updateTimer_.stop();

                if (asImp_().verbose_() && isatty(fileno(stdout)))
                    // make sure that the line currently holding the cursor is prestine
                    std::cout << clearRemainingLine
                              << std::flush;

                // tell the implementation that we're done with this iteration
                prePostProcessTimer_.start();
                asImp_().endIteration_(nextSolution, currentSolution);
                prePostProcessTimer_.stop();
            }
        }
        catch (const Dune::Exception& e)
        {
            if (asImp_().verbose_())
                std::cout << "Newton method caught exception: \""
                          << e.what() << "\"\n" << std::flush;

            prePostProcessTimer_.start();
            asImp_().failed_();
            prePostProcessTimer_.stop();

            return false;
        }
        catch (const NumericalProblem& e)
        {
            if (asImp_().verbose_())
                std::cout << "Newton method caught exception: \""
                          << e.what() << "\"\n" << std::flush;

            prePostProcessTimer_.start();
            asImp_().failed_();
            prePostProcessTimer_.stop();

            return false;
        }

        // clear current line on terminal
        if (asImp_().verbose_() && isatty(fileno(stdout)))
            std::cout << clearRemainingLine
                      << std::flush;

        // tell the implementation that we're done
        prePostProcessTimer_.start();
        asImp_().end_();
        prePostProcessTimer_.stop();

        // print the timing summary of the time step
        if (asImp_().verbose_()) {
            Scalar elapsedTot =
                linearizeTimer_.realTimeElapsed()
                + solveTimer_.realTimeElapsed()
                + updateTimer_.realTimeElapsed();
            std::cout << "Linearization/solve/update time: "
                      << linearizeTimer_.realTimeElapsed() << "("
                      << 100 * linearizeTimer_.realTimeElapsed()/elapsedTot << "%)/"
                      << solveTimer_.realTimeElapsed() << "("
                      << 100 * solveTimer_.realTimeElapsed()/elapsedTot << "%)/"
                      << updateTimer_.realTimeElapsed() << "("
                      << 100 * updateTimer_.realTimeElapsed()/elapsedTot << "%)"
                      << "\n" << std::flush;
        }


        // if we're not converged, tell the implementation that we've failed
        if (!asImp_().converged()) {
            prePostProcessTimer_.start();
            asImp_().failed_();
            prePostProcessTimer_.stop();
            return false;
        }

        // if we converged, tell the implementation that we've succeeded
        prePostProcessTimer_.start();
        asImp_().succeeded_();
        prePostProcessTimer_.stop();

        return true;
    }

    /*!
     * \brief Suggest a new time-step size based on the old time-step
     *        size.
     *
     * The default behavior is to suggest the old time-step size
     * scaled by the ratio between the target iterations and the
     * iterations required to actually solve the last time-step.
     */
    Scalar suggestTimeStepSize(Scalar oldDt) const
    {
        // be aggressive reducing the time-step size but
        // conservative when increasing it. the rationale is
        // that we want to avoid failing in the next time
        // integration which would be quite expensive
        if (numIterations_ > params_.targetIterations_) {
            Scalar percent = Scalar(numIterations_ - params_.targetIterations_) / params_.targetIterations_;
            Scalar nextDt = std::max(problem().minTimeStepSize(),
                                     oldDt / (Scalar{1.0} + percent));
            return nextDt;
        }

        Scalar percent = Scalar(params_.targetIterations_ - numIterations_) / params_.targetIterations_;
        Scalar nextDt = std::max(problem().minTimeStepSize(),
                                 oldDt*(Scalar{1.0} + percent / Scalar{1.2}));
        return nextDt;
    }

    /*!
     * \brief Message that should be printed for the user after the
     *        end of an iteration.
     */
    std::ostringstream& endIterMsg()
    { return endIterMsgStream_; }

    /*!
     * \brief Causes the solve() method to discared the structure of the linear system of
     *        equations the next time it is called.
     */
    void eraseMatrix()
    { linearSolver_.eraseMatrix(); }

    /*!
     * \brief Returns the linear solver backend object for external use.
     */
    LinearSolverBackend& linearSolver()
    { return linearSolver_; }

    /*!
     * \copydoc linearSolver()
     */
    const LinearSolverBackend& linearSolver() const
    { return linearSolver_; }

    const Timer& prePostProcessTimer() const
    { return prePostProcessTimer_; }

    const Timer& linearizeTimer() const
    { return linearizeTimer_; }

    const Timer& solveTimer() const
    { return solveTimer_; }

    const Timer& updateTimer() const
    { return updateTimer_; }

protected:
    /*!
     * \brief Returns true if the Newton method ought to be chatty.
     */
    bool verbose_() const
    {
        return params_.verbose_ && (comm_.rank() == 0);
    }

    /*!
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param u The initial solution
     */
    void begin_(const SolutionVector&)
    {
        numIterations_ = 0;

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
        endIterMsgStream_.str("");
        const auto& comm = simulator_.gridView().comm();
        bool succeeded = true;
        try {
            problem().beginIteration();
        }
        catch (const std::exception& e) {
            succeeded = false;

            std::cout << "rank " << simulator_.gridView().comm().rank()
                      << " caught an exception while pre-processing the problem:" << e.what()
                      << "\n"  << std::flush;
        }

        succeeded = comm.min(succeeded);

        if (!succeeded)
            throw NumericalProblem("pre processing of the problem failed");

        lastError_ = error_;
    }

    /*!
     * \brief Linearize the global non-linear system of equations associated with the
     *        spatial domain.
     */
    void linearizeDomain_()
    {
        model().linearizer().linearizeDomain();
    }

    void linearizeAuxiliaryEquations_()
    {
        model().linearizer().linearizeAuxiliaryEquations();
        model().linearizer().finalize();
    }

    void preSolve_(const SolutionVector&,
                   const GlobalEqVector& currentResidual)
    {
        const auto& constraintsMap = model().linearizer().constraintsMap();
        lastError_ = error_;
        Scalar newtonMaxError = params_.maxError_;

        // calculate the error as the maximum weighted tolerance of
        // the solution's residual
        error_ = 0;
        for (unsigned dofIdx = 0; dofIdx < currentResidual.size(); ++dofIdx) {
            // do not consider auxiliary DOFs for the error
            if (dofIdx >= model().numGridDof() || model().dofTotalVolume(dofIdx) <= 0.0)
                continue;

            // also do not consider DOFs which are constraint
            if (enableConstraints_()) {
                if (constraintsMap.count(dofIdx) > 0)
                    continue;
            }

            const auto& r = currentResidual[dofIdx];
            for (unsigned eqIdx = 0; eqIdx < r.size(); ++eqIdx)
                error_ = max(std::abs(r[eqIdx] * model().eqWeight(dofIdx, eqIdx)), error_);
        }

        // take the other processes into account
        error_ = comm_.max(error_);

        // make sure that the error never grows beyond the maximum
        // allowed one
        if (error_ > newtonMaxError)
            throw NumericalProblem("Newton: Error "+std::to_string(double(error_))
                                   + " is larger than maximum allowed error of "
                                   + std::to_string(double(newtonMaxError)));
    }

    /*!
     * \brief Update the error of the solution given the previous
     *        iteration.
     *
     * For our purposes, the error of a solution is defined as the
     * maximum of the weighted residual of a given solution.
     *
     * \param currentSolution The solution at the beginning the current iteration
     * \param currentResidual The residual (i.e., right-hand-side) of the current
     *                        iteration's solution.
     * \param solutionUpdate The difference between the current and the next solution
     */
    void postSolve_(const SolutionVector&,
                    const GlobalEqVector&,
                    GlobalEqVector& solutionUpdate)
    {
        // loop over the auxiliary modules and ask them to post process the solution
        // vector.
        auto& model = simulator_.model();
        const auto& comm = simulator_.gridView().comm();
        for (unsigned i = 0; i < model.numAuxiliaryModules(); ++i) {
            auto& auxMod = *model.auxiliaryModule(i);

            bool succeeded = true;
            try {
                auxMod.postSolve(solutionUpdate);
            }
            catch (const std::exception& e) {
                succeeded = false;

                std::cout << "rank " << simulator_.gridView().comm().rank()
                          << " caught an exception while post processing an auxiliary module:" << e.what()
                          << "\n"  << std::flush;
            }

            succeeded = comm.min(succeeded);

            if (!succeeded)
                throw NumericalProblem("post processing of an auxilary equation failed");
        }
    }

    /*!
     * \brief Update the current solution with a delta vector.
     *
     * Different update strategies, such as chopped updates can be
     * implemented by overriding this method. The default behavior is
     * use the standard Newton-Raphson update strategy, i.e.
     * \f[ u^{k+1} = u^k - \Delta u^k \f]
     *
     * \param nextSolution The solution vector after the current iteration
     * \param currentSolution The solution vector after the last iteration
     * \param solutionUpdate The delta vector as calculated by solving the linear system
     *                       of equations
     * \param currentResidual The residual vector of the current Newton-Raphson iteraton
     */
    void update_(SolutionVector& nextSolution,
                 const SolutionVector& currentSolution,
                 const GlobalEqVector& solutionUpdate,
                 const GlobalEqVector& currentResidual)
    {
        const auto& constraintsMap = model().linearizer().constraintsMap();

        // first, write out the current solution to make convergence
        // analysis possible
        asImp_().writeConvergence_(currentSolution, solutionUpdate);

        // make sure not to swallow non-finite values at this point
        if (!std::isfinite(solutionUpdate.one_norm()))
            throw NumericalProblem("Non-finite update!");

        size_t numGridDof = model().numGridDof();
        for (unsigned dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
            if (enableConstraints_()) {
                if (constraintsMap.count(dofIdx) > 0) {
                    const auto& constraints = constraintsMap.at(dofIdx);
                    asImp_().updateConstraintDof_(dofIdx,
                                                  nextSolution[dofIdx],
                                                  constraints);
                }
                else
                    asImp_().updatePrimaryVariables_(dofIdx,
                                                     nextSolution[dofIdx],
                                                     currentSolution[dofIdx],
                                                     solutionUpdate[dofIdx],
                                                     currentResidual[dofIdx]);
            }
            else
                asImp_().updatePrimaryVariables_(dofIdx,
                                                 nextSolution[dofIdx],
                                                 currentSolution[dofIdx],
                                                 solutionUpdate[dofIdx],
                                                 currentResidual[dofIdx]);
        }

        // update the DOFs of the auxiliary equations
        size_t numDof = model().numTotalDof();
        for (size_t dofIdx = numGridDof; dofIdx < numDof; ++dofIdx) {
            nextSolution[dofIdx] = currentSolution[dofIdx];
            nextSolution[dofIdx] -= solutionUpdate[dofIdx];
        }
    }

    /*!
     * \brief Update the primary variables for a degree of freedom which is constraint.
     */
    void updateConstraintDof_(unsigned,
                              PrimaryVariables& nextValue,
                              const Constraints& constraints)
    { nextValue = constraints; }

    /*!
     * \brief Update a single primary variables object.
     */
    void updatePrimaryVariables_(unsigned,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector&)
    {
        nextValue = currentValue;
        nextValue -= update;
    }

    /*!
     * \brief Write the convergence behaviour of the newton method to
     *        disk.
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
    void endIteration_(const SolutionVector&,
                       const SolutionVector&)
    {
        ++numIterations_;

        const auto& comm = simulator_.gridView().comm();
        bool succeeded = true;
        try {
            problem().endIteration();
        }
        catch (const std::exception& e) {
            succeeded = false;

            std::cout << "rank " << simulator_.gridView().comm().rank()
                      << " caught an exception while letting the problem post-process:" << e.what()
                      << "\n"  << std::flush;
        }

        succeeded = comm.min(succeeded);

        if (!succeeded)
            throw NumericalProblem("post processing of the problem failed");

        if (asImp_().verbose_()) {
            std::cout << "Newton iteration " << numIterations_ << ""
                      << " error: " << error_
                      << endIterMsg().str() << "\n" << std::flush;
        }
    }

    /*!
     * \brief Returns true iff another Newton iteration should be done.
     */
    bool proceed_() const
    {
        if (asImp_().numIterations() < 1)
            return true; // we always do at least one full iteration
        else if (asImp_().converged()) {
            // we are below the specified tolerance, so we don't have to
            // do more iterations
            return false;
        }
        else if (asImp_().numIterations() >= params_.maxIterations_) {
            // we have exceeded the allowed number of steps.  If the
            // error was reduced by a factor of at least 4,
            // in the last iterations we proceed even if we are above
            // the maximum number of steps
            return error_ * 4.0 < lastError_;
        }

        return true;
    }

    /*!
     * \brief Indicates that we're done solving the non-linear system
     *        of equations.
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
     * This method is called _after_ end_()
     */
    void failed_()
    { numIterations_ = params_.targetIterations_ * 2; }

    /*!
     * \brief Called if the Newton method was successful.
     *
     * This method is called _after_ end_()
     */
    void succeeded_()
    {}

    static bool enableConstraints_()
    { return getPropValue<TypeTag, Properties::EnableConstraints>(); }

    Simulator& simulator_;

    Timer prePostProcessTimer_;
    Timer linearizeTimer_;
    Timer solveTimer_;
    Timer updateTimer_;

    std::ostringstream endIterMsgStream_;

    Scalar error_;
    Scalar lastError_;
    NewtonMethodParams<Scalar> params_;

    // actual number of iterations done so far
    int numIterations_;

    // the linear solver
    LinearSolverBackend linearSolver_;

    // the collective communication used by the simulation (i.e. fake
    // or MPI)
    CollectiveCommunication comm_;

    // the object which writes the convergence behaviour of the Newton
    // method to disk
    ConvergenceWriter convergenceWriter_;

private:
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }
    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // namespace Opm

#endif
