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

#include "nullconvergencewriter.hh"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/timer.hh>
#include <opm/models/utils/timerguard.hh>

#include <opm/material/densead/Math.hpp>
#include <opm/material/common/Unused.hpp>

#include <opm/material/common/Exceptions.hpp>

#include <dune/istl/istlexception.hh>
#include <dune/common/classname.hh>
#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

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

BEGIN_PROPERTIES

//! The type tag on which the default properties for the Newton method
//! are attached
NEW_TYPE_TAG(NewtonMethod);

//! The simulation management class of the simulation
NEW_PROP_TAG(Simulator);

//! The physical model which we would like to solve
NEW_PROP_TAG(Problem);

//! The model describing the PDEs for the conservation quantities
NEW_PROP_TAG(Model);

//! The type of scalar values
NEW_PROP_TAG(Scalar);

//! Specifies the type of the actual Newton method
NEW_PROP_TAG(NewtonMethod);

//! Specifies the type of a solution
NEW_PROP_TAG(SolutionVector);

//! Specifies the type of a solution for a single degee of freedom
NEW_PROP_TAG(PrimaryVariables);

//! Specifies whether the problem to be simulated exhibits contraint degrees of freedom
NEW_PROP_TAG(EnableConstraints);

//! Specifies the type of objects which specify constraints for a single degee of freedom
NEW_PROP_TAG(Constraints);

//! Vector containing a quantity of for equation on the whole grid
NEW_PROP_TAG(GlobalEqVector);

//! Vector containing a quantity of for equation for a single degee of freedom
NEW_PROP_TAG(EqVector);

//! The class which linearizes the non-linear system of equations
NEW_PROP_TAG(Linearizer);

//! Specifies the type of a global Jacobian matrix
NEW_PROP_TAG(SparseMatrixAdapter);

//! Specifies the type of the linear solver to be used
NEW_PROP_TAG(LinearSolverBackend);

//! Specifies whether the Newton method should print messages or not
NEW_PROP_TAG(NewtonVerbose);

//! Specifies the type of the class which writes out the Newton convergence
NEW_PROP_TAG(NewtonConvergenceWriter);

//! Specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration
NEW_PROP_TAG(NewtonWriteConvergence);

//! Specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration
NEW_PROP_TAG(ConvergenceWriter);

/*!
 * \brief The value for the error below which convergence is declared
 *
 * This value can (and for the porous media models will) be changed to account for grid
 * scaling and other effects.
 */
NEW_PROP_TAG(NewtonTolerance);

//! The maximum error which may occur in a simulation before the
//! Newton method for the time step is aborted
NEW_PROP_TAG(NewtonMaxError);

/*!
 * \brief The number of iterations at which the Newton method
 *        should aim at.
 *
 * This is used to control the time-step size. The heuristic used
 * is to scale the last time-step size by the deviation of the
 * number of iterations used from the target steps.
 */
NEW_PROP_TAG(NewtonTargetIterations);

//! Number of maximum iterations for the Newton method.
NEW_PROP_TAG(NewtonMaxIterations);

// set default values for the properties
SET_TYPE_PROP(NewtonMethod, NewtonMethod, Opm::NewtonMethod<TypeTag>);
SET_TYPE_PROP(NewtonMethod, NewtonConvergenceWriter, Opm::NullConvergenceWriter<TypeTag>);
SET_BOOL_PROP(NewtonMethod, NewtonWriteConvergence, false);
SET_BOOL_PROP(NewtonMethod, NewtonVerbose, true);
SET_SCALAR_PROP(NewtonMethod, NewtonTolerance, 1e-8);
// set the abortion tolerace to some very large value. if not
// overwritten at run-time this basically disables abortions
SET_SCALAR_PROP(NewtonMethod, NewtonMaxError, 1e100);
SET_INT_PROP(NewtonMethod, NewtonTargetIterations, 10);
SET_INT_PROP(NewtonMethod, NewtonMaxIterations, 18);

END_PROPERTIES

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
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, Linearizer) Linearizer;
    typedef typename GET_PROP_TYPE(TypeTag, LinearSolverBackend) LinearSolverBackend;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonConvergenceWriter) ConvergenceWriter;

    typedef typename Dune::MPIHelper::MPICommunicator Communicator;
    typedef Dune::CollectiveCommunication<Communicator> CollectiveCommunication;

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
        tolerance_ = EWOMS_GET_PARAM(TypeTag, Scalar, NewtonTolerance);

        numIterations_ = 0;
    }

    /*!
     * \brief Register all run-time parameters for the Newton method.
     */
    static void registerParameters()
    {
        LinearSolverBackend::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, NewtonVerbose,
                             "Specify whether the Newton method should inform "
                             "the user about its progress or not");
        EWOMS_REGISTER_PARAM(TypeTag, bool, NewtonWriteConvergence,
                             "Write the convergence behaviour of the Newton "
                             "method to a VTK file");
        EWOMS_REGISTER_PARAM(TypeTag, int, NewtonTargetIterations,
                             "The 'optimum' number of Newton iterations per "
                             "time step");
        EWOMS_REGISTER_PARAM(TypeTag, int, NewtonMaxIterations,
                             "The maximum number of Newton iterations per time "
                             "step");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, NewtonTolerance,
                             "The maximum raw error tolerated by the Newton"
                             "method for considering a solution to be "
                             "converged");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, NewtonMaxError,
                             "The maximum error tolerated by the Newton "
                             "method to which does not cause an abort");
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
    { return tolerance_; }

    /*!
     * \brief Set the current tolerance at which the Newton method considers itself to
     *        be converged.
     */
    void setTolerance(Scalar value)
    { tolerance_ = value; }

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

        Opm::TimerGuard prePostProcessTimerGuard(prePostProcessTimer_);

        // tell the implementation that we begin solving
        prePostProcessTimer_.start();
        asImp_().begin_(nextSolution);
        prePostProcessTimer_.stop();

        try {
            Opm::TimerGuard innerPrePostProcessTimerGuard(prePostProcessTimer_);
            Opm::TimerGuard linearizeTimerGuard(linearizeTimer_);
            Opm::TimerGuard updateTimerGuard(updateTimer_);
            Opm::TimerGuard solveTimerGuard(solveTimer_);

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
        catch (const Opm::NumericalIssue& e)
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
        if (numIterations_ > targetIterations_()) {
            Scalar percent = Scalar(numIterations_ - targetIterations_())/targetIterations_();
            Scalar nextDt = std::max(problem().minTimeStepSize(),
                                     oldDt/(1.0 + percent));
            return nextDt;
        }

        Scalar percent = Scalar(targetIterations_() - numIterations_)/targetIterations_();
        Scalar nextDt = std::max(problem().minTimeStepSize(),
                                 oldDt*(1.0 + percent/1.2));
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

    const Opm::Timer& prePostProcessTimer() const
    { return prePostProcessTimer_; }

    const Opm::Timer& linearizeTimer() const
    { return linearizeTimer_; }

    const Opm::Timer& solveTimer() const
    { return solveTimer_; }

    const Opm::Timer& updateTimer() const
    { return updateTimer_; }

protected:
    /*!
     * \brief Returns true if the Newton method ought to be chatty.
     */
    bool verbose_() const
    {
        return EWOMS_GET_PARAM(TypeTag, bool, NewtonVerbose) && (comm_.rank() == 0);
    }

    /*!
     * \brief Called before the Newton method is applied to an
     *        non-linear system of equations.
     *
     * \param u The initial solution
     */
    void begin_(const SolutionVector& u  OPM_UNUSED)
    {
        numIterations_ = 0;

        if (EWOMS_GET_PARAM(TypeTag, bool, NewtonWriteConvergence))
            convergenceWriter_.beginTimeStep();
    }

    /*!
     * \brief Indicates the beginning of a Newton iteration.
     */
    void beginIteration_()
    {
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
#if ! DUNE_VERSION_NEWER(DUNE_COMMON, 2,5)
        catch (const Dune::Exception& e)
        {
            succeeded = false;

            std::cout << "rank " << simulator_.gridView().comm().rank()
                      << " caught an exception while pre-processing the problem:" << e.what()
                      << "\n"  << std::flush;
        }
#endif

        succeeded = comm.min(succeeded);

        if (!succeeded)
            throw Opm::NumericalIssue("pre processing of the problem failed");

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

    void preSolve_(const SolutionVector& currentSolution  OPM_UNUSED,
                   const GlobalEqVector& currentResidual)
    {
        const auto& constraintsMap = model().linearizer().constraintsMap();
        lastError_ = error_;
        Scalar newtonMaxError = EWOMS_GET_PARAM(TypeTag, Scalar, NewtonMaxError);

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
                error_ = Opm::max(std::abs(r[eqIdx] * model().eqWeight(dofIdx, eqIdx)), error_);
        }

        // take the other processes into account
        error_ = comm_.max(error_);

        // make sure that the error never grows beyond the maximum
        // allowed one
        if (error_ > newtonMaxError)
            throw Opm::NumericalIssue("Newton: Error "+std::to_string(double(error_))
                                        +" is larger than maximum allowed error of "
                                        +std::to_string(double(newtonMaxError)));
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
    void postSolve_(const SolutionVector& currentSolution OPM_UNUSED,
                    const GlobalEqVector& currentResidual OPM_UNUSED,
                    GlobalEqVector& solutionUpdate OPM_UNUSED)
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
#if ! DUNE_VERSION_NEWER(DUNE_COMMON, 2,5)
            catch (const Dune::Exception& e)
            {
                succeeded = false;

                std::cout << "rank " << simulator_.gridView().comm().rank()
                          << " caught an exception while post processing an auxiliary module:" << e.what()
                          << "\n"  << std::flush;
            }
#endif

            succeeded = comm.min(succeeded);

            if (!succeeded)
                throw Opm::NumericalIssue("post processing of an auxilary equation failed");
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
            throw Opm::NumericalIssue("Non-finite update!");

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
    void updateConstraintDof_(unsigned globalDofIdx  OPM_UNUSED,
                              PrimaryVariables& nextValue,
                              const Constraints& constraints)
    { nextValue = constraints; }

    /*!
     * \brief Update a single primary variables object.
     */
    void updatePrimaryVariables_(unsigned globalDofIdx  OPM_UNUSED,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& currentResidual  OPM_UNUSED)
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
        if (EWOMS_GET_PARAM(TypeTag, bool, NewtonWriteConvergence)) {
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
    void endIteration_(const SolutionVector& nextSolution  OPM_UNUSED,
                       const SolutionVector& currentSolution  OPM_UNUSED)
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
#if ! DUNE_VERSION_NEWER(DUNE_COMMON, 2,5)
        catch (const Dune::Exception& e)
        {
            succeeded = false;

            std::cout << "rank " << simulator_.gridView().comm().rank()
                      << " caught an exception while letting the problem post process:" << e.what()
                      << "\n"  << std::flush;
        }
#endif

        succeeded = comm.min(succeeded);

        if (!succeeded)
            throw Opm::NumericalIssue("post processing of the problem failed");

        if (asImp_().verbose_()) {
            std::cout << "Newton iteration " << numIterations_ << ""
                      << " error: " << error_
                      << endIterMsg().str() << "\n" << std::flush;
        }

        endIterMsgStream_.str("");
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
        else if (asImp_().numIterations() >= asImp_().maxIterations_()) {
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
        if (EWOMS_GET_PARAM(TypeTag, bool, NewtonWriteConvergence))
            convergenceWriter_.endTimeStep();
    }

    /*!
     * \brief Called if the Newton method broke down.
     *
     * This method is called _after_ end_()
     */
    void failed_()
    { numIterations_ = targetIterations_() * 2; }

    /*!
     * \brief Called if the Newton method was successful.
     *
     * This method is called _after_ end_()
     */
    void succeeded_()
    {}

    // optimal number of iterations we want to achieve
    int targetIterations_() const
    { return EWOMS_GET_PARAM(TypeTag, int, NewtonTargetIterations); }
    // maximum number of iterations we do before giving up
    int maxIterations_() const
    { return EWOMS_GET_PARAM(TypeTag, int, NewtonMaxIterations); }

    static bool enableConstraints_()
    { return GET_PROP_VALUE(TypeTag, EnableConstraints); }

    Simulator& simulator_;

    Opm::Timer prePostProcessTimer_;
    Opm::Timer linearizeTimer_;
    Opm::Timer solveTimer_;
    Opm::Timer updateTimer_;

    std::ostringstream endIterMsgStream_;

    Scalar error_;
    Scalar lastError_;
    Scalar tolerance_;

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
