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
 *
 * \copydoc Opm::Simulator
 */
#ifndef EWOMS_SIMULATOR_HH
#define EWOMS_SIMULATOR_HH

#include <opm/models/io/restart.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/timer.hh>
#include <opm/models/utils/timerguard.hh>
#include <opm/models/parallel/mpiutil.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <memory>

#define EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(code)                      \
    {                                                                   \
        const auto& comm = Dune::MPIHelper::getCollectiveCommunication(); \
        bool exceptionThrown = false;                                   \
        try { code; }                                                   \
        catch (const Dune::Exception& e) {                              \
            exceptionThrown = true;                                     \
            std::cerr << "Process " << comm.rank() << " threw a fatal exception: " \
                      << e.what() << ". Abort!" << std::endl;           \
        }                                                               \
        catch (const std::exception& e) {                               \
            exceptionThrown = true;                                     \
            std::cerr << "Process " << comm.rank() << " threw a fatal exception: " \
                      << e.what() << ". Abort!" << std::endl;           \
        }                                                               \
        catch (...) {                                                   \
            exceptionThrown = true;                                     \
            std::cerr << "Process " << comm.rank() << " threw a fatal exception. " \
                      <<" Abort!" << std::endl;                         \
        }                                                               \
                                                                        \
        if (comm.max(exceptionThrown))                                  \
            std::abort();                                               \
    }

namespace Opm {

/*!
 * \ingroup Common
 *
 * \brief Manages the initializing and running of time dependent
 *        problems.
 *
 * This class instantiates the grid, the model and the problem to be
 * simlated and runs the simulation loop. The time axis is treated as
 * a sequence of "episodes" which are defined as time intervals for
 * which the problem exhibits boundary conditions and source terms
 * that do not depend on time.
 */
template <class TypeTag>
class Simulator
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Vanguard = GetPropType<TypeTag, Properties::Vanguard>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

public:
    // do not allow to copy simulators around
    Simulator(const Simulator& ) = delete;

    Simulator(bool verbose = true)
    {
        TimerGuard setupTimerGuard(setupTimer_);

        setupTimer_.start();

        const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
        verbose_ = verbose && comm.rank() == 0;

        timeStepIdx_ = 0;
        startTime_ = 0.0;
        time_ = 0.0;
        endTime_ = EWOMS_GET_PARAM(TypeTag, Scalar, EndTime);
        timeStepSize_ = EWOMS_GET_PARAM(TypeTag, Scalar, InitialTimeStepSize);
        assert(timeStepSize_ > 0);
        const std::string& predetTimeStepFile =
            EWOMS_GET_PARAM(TypeTag, std::string, PredeterminedTimeStepsFile);
        if (!predetTimeStepFile.empty()) {
            std::ifstream is(predetTimeStepFile);
            while (!is.eof()) {
                Scalar dt;
                is >> dt;
                forcedTimeSteps_.push_back(dt);
            }
        }

        episodeIdx_ = 0;
        episodeStartTime_ = 0;
        episodeLength_ = std::numeric_limits<Scalar>::max();

        finished_ = false;

        if (verbose_)
            std::cout << "Allocating the simulation vanguard\n" << std::flush;

        int exceptionThrown = 0;
        std::string what;
        try
        { vanguard_.reset(new Vanguard(*this)); }
        catch (const std::exception& e) {
            exceptionThrown = 1;
            what = e.what();
            if (comm.size() > 1) {
                what += " (on rank " + std::to_string(comm.rank()) + ")";
            }
            if (verbose_)
                std::cerr << "Rank " << comm.rank() << " threw an exception: " << e.what() << std::endl;
        }

        if (comm.max(exceptionThrown)) {
            auto all_what = gatherStrings(what);
            assert(!all_what.empty());
            throw std::runtime_error("Allocating the simulation vanguard failed: " + all_what.front());
        }

        if (verbose_)
            std::cout << "Distributing the vanguard's data\n" << std::flush;

        try
        { vanguard_->loadBalance(); }
        catch (const std::exception& e) {
            exceptionThrown = 1;
            what = e.what();
            if (comm.size() > 1) {
                what += " (on rank " + std::to_string(comm.rank()) + ")";
            }
            if (verbose_)
                std::cerr << "Rank " << comm.rank() << " threw an exception: " << e.what() << std::endl;
        }

        if (comm.max(exceptionThrown)) {
            auto all_what = gatherStrings(what);
            assert(!all_what.empty());
            throw std::runtime_error("Could not distribute the vanguard data: " + all_what.front());
        }

        if (verbose_)
            std::cout << "Allocating the model\n" << std::flush;
        model_.reset(new Model(*this));

        if (verbose_)
            std::cout << "Allocating the problem\n" << std::flush;
        problem_.reset(new Problem(*this));

        if (verbose_)
            std::cout << "Initializing the model\n" << std::flush;

        try
        { model_->finishInit(); }
        catch (const std::exception& e) {
            exceptionThrown = 1;
            what = e.what();
            if (comm.size() > 1) {
                what += " (on rank " + std::to_string(comm.rank()) + ")";
            }
            if (verbose_)
                std::cerr << "Rank " << comm.rank() << " threw an exception: " << e.what() << std::endl;
        }

        if (comm.max(exceptionThrown)) {
            auto all_what = gatherStrings(what);
            assert(!all_what.empty());
            throw std::runtime_error("Could not initialize the model: " + all_what.front());
        }

        if (verbose_)
            std::cout << "Initializing the problem\n" << std::flush;

        try
        { problem_->finishInit(); }
        catch (const std::exception& e) {
            exceptionThrown = 1;
            what = e.what();
            if (comm.size() > 1) {
                what += " (on rank " + std::to_string(comm.rank()) + ")";
            }
            if (verbose_)
                std::cerr << "Rank " << comm.rank() << " threw an exception: " << e.what() << std::endl;
        }

        if (comm.max(exceptionThrown)) {
            auto all_what = gatherStrings(what);
            assert(!all_what.empty());
            throw std::runtime_error("Could not initialize the problem: " + all_what.front());
        }

        setupTimer_.stop();

        if (verbose_)
            std::cout << "Simulator successfully set up\n" << std::flush;
    }

    /*!
     * \brief Registers all runtime parameters used by the simulation.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EndTime,
                             "The simulation time at which the simulation is finished [s]");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, InitialTimeStepSize,
                             "The size of the initial time step [s]");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, RestartTime,
                             "The simulation time at which a restart should be attempted [s]");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, PredeterminedTimeStepsFile,
                             "A file with a list of predetermined time step sizes (one "
                             "time step per line)");

        Vanguard::registerParameters();
        Model::registerParameters();
        Problem::registerParameters();
    }

    /*!
     * \brief Return a reference to the grid manager of simulation
     */
    Vanguard& vanguard()
    { return *vanguard_; }

    /*!
     * \brief Return a reference to the grid manager of simulation
     */
    const Vanguard& vanguard() const
    { return *vanguard_; }

    /*!
     * \brief Return the grid view for which the simulation is done
     */
    const GridView& gridView() const
    { return vanguard_->gridView(); }

    /*!
     * \brief Return the physical model used in the simulation
     */
    Model& model()
    { return *model_; }

    /*!
     * \brief Return the physical model used in the simulation
     */
    const Model& model() const
    { return *model_; }

    /*!
     * \brief Return the object which specifies the pysical setup of
     *        the simulation
     */
    Problem& problem()
    { return *problem_; }

    /*!
     * \brief Return the object which specifies the pysical setup of
     *        the simulation
     */
    const Problem& problem() const
    { return *problem_; }

    /*!
     * \brief Set the time of the start of the simulation.
     *
     * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
     */
    void setStartTime(Scalar t)
    { startTime_ = t; }

    /*!
     * \brief Return the time of the start of the simulation.
     */
    Scalar startTime() const
    { return startTime_; }

    /*!
     * \brief Set the current simulated time, don't change the current
     *        time step index.
     *
     * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
     */
    void setTime(Scalar t)
    { time_ = t; }

    /*!
     * \brief Set the current simulated time and the time step index.
     *
     * \param t The time \f$\mathrm{[s]}\f$ which should be jumped to
     * \param stepIdx The new time step index
     */
    void setTime(Scalar t, unsigned stepIdx)
    {
        time_ = t;
        timeStepIdx_ = stepIdx;
    }

    /*!
     * \brief Return the number of seconds of simulated time which have elapsed since the
     *        start time.
     *
     * To get the time after the time integration, you have to add
     * timeStepSize() to time().
     */
    Scalar time() const
    { return time_; }

    /*!
     * \brief Set the time of simulated seconds at which the simulation runs.
     *
     * \param t The time \f$\mathrm{[s]}\f$ at which the simulation is finished
     */
    void setEndTime(Scalar t)
    { endTime_ = t; }

    /*!
     * \brief Returns the number of (simulated) seconds which the simulation
     *        runs.
     */
    Scalar endTime() const
    { return endTime_; }

    /*!
     * \brief Returns a reference to the timer object which measures the time needed to
     *        set up and initialize the simulation
     */
    const Timer& setupTimer() const
    { return setupTimer_; }

    /*!
     * \brief Returns a reference to the timer object which measures the time needed to
     *        run the simulation
     */
    const Timer& executionTimer() const
    { return executionTimer_; }
    Timer& executionTimer()
    { return executionTimer_; }

    /*!
     * \brief Returns a reference to the timer object which measures the time needed for
     *        pre- and postprocessing of the solutions.
     */
    const Timer& prePostProcessTimer() const
    { return prePostProcessTimer_; }

    /*!
     * \brief Returns a reference to the timer object which measures the time needed for
     *        linarizing the solutions.
     */
    const Timer& linearizeTimer() const
    { return linearizeTimer_; }

    /*!
     * \brief Returns a reference to the timer object which measures the time needed by
     *        the solver.
     */
    const Timer& solveTimer() const
    { return solveTimer_; }

    /*!
     * \brief Returns a reference to the timer object which measures the time needed to
     *        the solutions of the non-linear system of equations.
     */
    const Timer& updateTimer() const
    { return updateTimer_; }

    /*!
     * \brief Returns a reference to the timer object which measures the time needed to
     *        write the visualization output
     */
    const Timer& writeTimer() const
    { return writeTimer_; }

    /*!
     * \brief Set the current time step size to a given value.
     *
     * If the step size would exceed the length of the current
     * episode, the timeStep() method will take care that the step
     * size won't exceed the episode or the end of the simulation,
     * though.
     *
     * \param timeStepSize The new value for the time step size \f$\mathrm{[s]}\f$
     */
    void setTimeStepSize(Scalar value)
    {
         timeStepSize_ = value;
    }

    /*!
     * \brief Set the current time step index to a given value.
     *
     * \param timeStepIndex The new value for the time step index
     */
    void setTimeStepIndex(unsigned value)
    { timeStepIdx_ = value; }

    /*!
     * \brief Returns the time step length \f$\mathrm{[s]}\f$ so that we
     *        don't miss the beginning of the next episode or cross
     *        the end of the simlation.
     */
    Scalar timeStepSize() const
    { return timeStepSize_; }

    /*!
     * \brief Returns number of time steps which have been
     *        executed since the beginning of the simulation.
     */
    int timeStepIndex() const
    { return timeStepIdx_; }

    /*!
     * \brief Specify whether the simulation is finished
     *
     * \param yesno If true the simulation is considered finished
     *              before the end time is reached, else it is only
     *              considered finished if the end time is reached.
     */
    void setFinished(bool yesno = true)
    { finished_ = yesno; }

    /*!
     * \brief Returns true if the simulation is finished.
     *
     * This is the case if either setFinished(true) has been called or
     * if the end time is reached.
     */
    bool finished() const
    {
        assert(timeStepSize_ >= 0.0);
        Scalar eps =
            std::max(Scalar(std::abs(this->time())), timeStepSize())
            *std::numeric_limits<Scalar>::epsilon()*1e3;
        return finished_ || (this->time()*(1.0 + eps) >= endTime());
    }

    /*!
     * \brief Returns true if the simulation is finished after the
     *        time level is incremented by the current time step size.
     */
    bool willBeFinished() const
    {
        static const Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e3;

        return finished_ || (this->time() + timeStepSize_)*(1.0 + eps) >= endTime();
    }

    /*!
     * \brief Aligns the time step size to the episode boundary and to
     *        the end time of the simulation.
     */
    Scalar maxTimeStepSize() const
    {
        if (finished())
            return 0.0;

        return std::min(episodeMaxTimeStepSize(),
                        std::max<Scalar>(0.0, endTime() - this->time()));
    }

    /*!
     * \brief Change the current episode of the simulation.
     *
     * \param episodeStartTime Time when the episode began \f$\mathrm{[s]}\f$
     * \param episodeLength Length of the episode \f$\mathrm{[s]}\f$
     */
    void startNextEpisode(Scalar episodeStartTime, Scalar episodeLength)
    {
        ++episodeIdx_;
        episodeStartTime_ = episodeStartTime;
        episodeLength_ = episodeLength;
    }

    /*!
     * \brief Start the next episode, but don't change the episode
     *        identifier.
     *
     * \param len Length of the episode \f$\mathrm{[s]}\f$, infinite if not
     *            specified.
     */
    void startNextEpisode(Scalar len = std::numeric_limits<Scalar>::max())
    {
        ++episodeIdx_;
        episodeStartTime_ = startTime_ + time_;
        episodeLength_ = len;
    }

    /*!
     * \brief Sets the index of the current episode.
     *
     * Use this method with care!
     */
    void setEpisodeIndex(int episodeIdx)
    { episodeIdx_ = episodeIdx; }

    /*!
     * \brief Returns the index of the current episode.
     *
     * The first episode has the index 0.
     */
    int episodeIndex() const
    { return episodeIdx_; }

    /*!
     * \brief Returns the absolute time when the current episode
     *        started \f$\mathrm{[s]}\f$.
     */
    Scalar episodeStartTime() const
    { return episodeStartTime_; }

    /*!
     * \brief Sets the length in seconds of the current episode.
     *
     * Use this method with care!
     */
    void setEpisodeLength(Scalar dt)
    { episodeLength_ = dt; }

    /*!
     * \brief Returns the length of the current episode in
     *        simulated time \f$\mathrm{[s]}\f$.
     */
    Scalar episodeLength() const
    { return episodeLength_; }

    /*!
     * \brief Returns true if the current episode has just been started at the
     *        current time.
     */
    bool episodeStarts() const
    {
        static const Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e3;

        return this->time() <= (episodeStartTime_ - startTime())*(1 + eps);
    }

    /*!
     * \brief Returns true if the current episode is finished at the
     *        current time.
     */
    bool episodeIsOver() const
    {
        static const Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e3;

        return this->time() >= (episodeStartTime_ - startTime() + episodeLength())*(1 - eps);
    }

    /*!
     * \brief Returns true if the current episode will be finished
     *        after the current time step.
     */
    bool episodeWillBeOver() const
    {
        static const Scalar eps = std::numeric_limits<Scalar>::epsilon()*1e3;

        return this->time() + timeStepSize()
            >=  (episodeStartTime_ - startTime() + episodeLength())*(1 - eps);
    }

    /*!
     * \brief Aligns the time step size to the episode boundary if the
     *        current time step crosses the boundary of the current episode.
     */
    Scalar episodeMaxTimeStepSize() const
    {
        // if the current episode is over and the simulation
        // wants to give it some extra time, we will return
        // the time step size it suggested instead of trying
        // to align it to the end of the episode.
        if (episodeIsOver())
            return 0.0;

        // make sure that we don't exceed the end of the
        // current episode.
        return std::max<Scalar>(0.0,
                                (episodeStartTime() + episodeLength())
                                - (this->time() + this->startTime()));
    }

    /*
     * \}
     */

    /*!
     * \brief Runs the simulation using a given problem class.
     *
     * This method makes sure that time steps sizes are aligned to
     * episode boundaries, amongst other stuff.
     */
    void run()
    {
        // create TimerGuard objects to hedge for exceptions
        TimerGuard setupTimerGuard(setupTimer_);
        TimerGuard executionTimerGuard(executionTimer_);
        TimerGuard prePostProcessTimerGuard(prePostProcessTimer_);
        TimerGuard writeTimerGuard(writeTimer_);

        setupTimer_.start();
        Scalar restartTime = EWOMS_GET_PARAM(TypeTag, Scalar, RestartTime);
        if (restartTime > -1e30) {
            // try to restart a previous simulation
            time_ = restartTime;

            Restart res;
            EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(res.deserializeBegin(*this, time_));
            if (verbose_)
                std::cout << "Deserialize from file '" << res.fileName() << "'\n" << std::flush;
            EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(this->deserialize(res));
            EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->deserialize(res));
            EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(model_->deserialize(res));
            EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(res.deserializeEnd());
            if (verbose_)
                std::cout << "Deserialization done."
                          << " Simulator time: " << time() << humanReadableTime(time())
                          << " Time step index: " << timeStepIndex()
                          << " Episode index: " << episodeIndex()
                          << "\n" << std::flush;
        }
        else {
            // if no restart is done, apply the initial solution
            if (verbose_)
                std::cout << "Applying the initial solution of the \"" << problem_->name()
                          << "\" problem\n" << std::flush;

            Scalar oldTimeStepSize = timeStepSize_;
            int oldTimeStepIdx = timeStepIdx_;
            timeStepSize_ = 0.0;
            timeStepIdx_ = -1;

            EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(model_->applyInitialSolution());

            // write initial condition
            if (problem_->shouldWriteOutput())
                EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->writeOutput());

            timeStepSize_ = oldTimeStepSize;
            timeStepIdx_ = oldTimeStepIdx;
        }
        setupTimer_.stop();

        executionTimer_.start();
        bool episodeBegins = episodeIsOver() || (timeStepIdx_ == 0);
        // do the time steps
        while (!finished()) {
            prePostProcessTimer_.start();
            if (episodeBegins) {
                // notify the problem that a new episode has just been
                // started.
                EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->beginEpisode());

                if (finished()) {
                    // the problem can chose to terminate the simulation in
                    // beginEpisode(), so we have handle this case.
                    EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->endEpisode());
                    prePostProcessTimer_.stop();

                    break;
                }
            }
            episodeBegins = false;

            if (verbose_) {
                std::cout << "Begin time step " << timeStepIndex() + 1 << ". "
                          << "Start time: " << this->time() << " seconds" << humanReadableTime(this->time())
                          << ", step size: " << timeStepSize() << " seconds" << humanReadableTime(timeStepSize())
                          << "\n";
            }

            // pre-process the current solution
            EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->beginTimeStep());

            if (finished()) {
                // the problem can chose to terminate the simulation in
                // beginTimeStep(), so we have handle this case.
                EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->endTimeStep());
                EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->endEpisode());
                prePostProcessTimer_.stop();

                break;
            }
            prePostProcessTimer_.stop();

            try {
                // execute the time integration scheme
                problem_->timeIntegration();
            }
            catch (...) {
                // exceptions in the time integration might be recoverable. clean up in
                // case they are
                const auto& model = problem_->model();
                prePostProcessTimer_ += model.prePostProcessTimer();
                linearizeTimer_ += model.linearizeTimer();
                solveTimer_ += model.solveTimer();
                updateTimer_ += model.updateTimer();

                throw;
            }

            const auto& model = problem_->model();
            prePostProcessTimer_ += model.prePostProcessTimer();
            linearizeTimer_ += model.linearizeTimer();
            solveTimer_ += model.solveTimer();
            updateTimer_ += model.updateTimer();

            // post-process the current solution
            prePostProcessTimer_.start();
            EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->endTimeStep());
            prePostProcessTimer_.stop();

            // write the result to disk
            writeTimer_.start();
            if (problem_->shouldWriteOutput())
                EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->writeOutput());
            writeTimer_.stop();

            // do the next time integration
            Scalar oldDt = timeStepSize();
            EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->advanceTimeLevel());

            if (verbose_) {
                std::cout << "Time step " << timeStepIndex() + 1 << " done. "
                          << "CPU time: " << executionTimer_.realTimeElapsed() << " seconds" << humanReadableTime(executionTimer_.realTimeElapsed())
                          << ", end time: " << this->time() + oldDt << " seconds" << humanReadableTime(this->time() + oldDt)
                          << ", step size: " << oldDt << " seconds" << humanReadableTime(oldDt)
                          << "\n" << std::flush;
            }

            // advance the simulated time by the current time step size
            time_ += oldDt;
            ++timeStepIdx_;

            prePostProcessTimer_.start();
            // notify the problem if an episode is finished
            if (episodeIsOver()) {
                // Notify the problem about the end of the current episode...
                EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->endEpisode());
                episodeBegins = true;
            }
            else {
                Scalar dt;
                if (timeStepIdx_ < static_cast<int>(forcedTimeSteps_.size()))
                    // use the next time step size from the input file
                    dt = forcedTimeSteps_[timeStepIdx_];
                else
                    // ask the problem to provide the next time step size
                    dt = std::min(maxTimeStepSize(), problem_->nextTimeStepSize());
                assert(finished() || dt > 0);
                setTimeStepSize(dt);
            }
            prePostProcessTimer_.stop();

            // write restart file if mandated by the problem
            writeTimer_.start();
            if (problem_->shouldWriteRestartFile())
                EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(serialize());
            writeTimer_.stop();
        }
        executionTimer_.stop();

        EWOMS_CATCH_PARALLEL_EXCEPTIONS_FATAL(problem_->finalize());
    }

    /*!
     * \brief Given a time step size in seconds, return it in a format which is more
     *        easily parsable by humans.
     *
     * e.g. 874000.0 will become "10.12 days"
     */
    static std::string humanReadableTime(Scalar timeInSeconds, bool isAmendment=true)
    {
        std::ostringstream oss;
        oss << std::setprecision(4);
        if (isAmendment)
            oss << " (";
        if (timeInSeconds >= 365.25*24*60*60) {
            int years = static_cast<int>(timeInSeconds/(365.25*24*60*60));
            int days = static_cast<int>((timeInSeconds - years*(365.25*24*60*60))/(24*60*60));

            double accuracy = 1e-2;
            double hours =
                std::round(1.0/accuracy*
                           (timeInSeconds
                            - years*(365.25*24*60*60)
                            - days*(24*60*60))/(60*60))
                *accuracy;

            oss << years << " years, " << days << " days, "  << hours << " hours";
        }
        else if (timeInSeconds >= 24.0*60*60) {
            int days = static_cast<int>(timeInSeconds/(24*60*60));
            int hours = static_cast<int>((timeInSeconds - days*(24*60*60))/(60*60));

            double accuracy = 1e-2;
            double minutes =
                std::round(1.0/accuracy*
                           (timeInSeconds
                            - days*(24*60*60)
                            - hours*(60*60))/60)
                *accuracy;

            oss << days << " days, " << hours << " hours, " << minutes << " minutes";
        }
        else if (timeInSeconds >= 60.0*60) {
            int hours = static_cast<int>(timeInSeconds/(60*60));
            int minutes = static_cast<int>((timeInSeconds - hours*(60*60))/60);

            double accuracy = 1e-2;
            double seconds =
                std::round(1.0/accuracy*
                           (timeInSeconds
                            - hours*(60*60)
                            - minutes*60))
                *accuracy;

            oss << hours << " hours, " << minutes << " minutes, "  << seconds << " seconds";
        }
        else if (timeInSeconds >= 60.0) {
            int minutes = static_cast<int>(timeInSeconds/60);

            double accuracy = 1e-3;
            double seconds =
                std::round(1.0/accuracy*
                           (timeInSeconds
                            - minutes*60))
                *accuracy;

            oss << minutes << " minutes, "  << seconds << " seconds";
        }
        else if (!isAmendment)
            oss << timeInSeconds << " seconds";
        else
            return "";
        if (isAmendment)
            oss << ")";

        return oss.str();
    }

    /*!
     * \name Saving/restoring the simulation state
     * \{
     */

    /*!
     * \brief This method writes the complete state of the simulation
     *        to the harddisk.
     *
     * The file will start with the prefix returned by the name()
     * method, has the current time of the simulation clock in it's
     * name and uses the extension <tt>.ers</tt>. (Ewoms ReStart
     * file.)  See Opm::Restart for details.
     */
    void serialize()
    {
        using Restarter = Restart;
        Restarter res;
        res.serializeBegin(*this);
        if (gridView().comm().rank() == 0)
            std::cout << "Serialize to file '" << res.fileName() << "'"
                      << ", next time step size: " << timeStepSize()
                      << "\n" << std::flush;

        this->serialize(res);
        problem_->serialize(res);
        model_->serialize(res);
        res.serializeEnd();
    }

    /*!
     * \brief Write the time manager's state to a restart file.
     *
     * \tparam Restarter The type of the object which takes care to serialize
     *                   data
     * \param restarter The serializer object
     */
    template <class Restarter>
    void serialize(Restarter& restarter)
    {
        restarter.serializeSectionBegin("Simulator");
        restarter.serializeStream()
            << episodeIdx_ << " "
            << episodeStartTime_ << " "
            << episodeLength_ << " "
            << startTime_ << " "
            << time_ << " "
            << timeStepIdx_ << " ";
        restarter.serializeSectionEnd();
    }

    /*!
     * \brief Read the time manager's state from a restart file.
     *
     * \tparam Restarter The type of the object which takes care to deserialize
     *                   data
     * \param restarter The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter& restarter)
    {
        restarter.deserializeSectionBegin("Simulator");
        restarter.deserializeStream()
            >> episodeIdx_
            >> episodeStartTime_
            >> episodeLength_
            >> startTime_
            >> time_
            >> timeStepIdx_;
        restarter.deserializeSectionEnd();
    }

private:
    std::unique_ptr<Vanguard> vanguard_;
    std::unique_ptr<Model> model_;
    std::unique_ptr<Problem> problem_;

    int episodeIdx_;
    Scalar episodeStartTime_;
    Scalar episodeLength_;

    Timer setupTimer_;
    Timer executionTimer_;
    Timer prePostProcessTimer_;
    Timer linearizeTimer_;
    Timer solveTimer_;
    Timer updateTimer_;
    Timer writeTimer_;

    std::vector<Scalar> forcedTimeSteps_;
    Scalar startTime_;
    Scalar time_;
    Scalar endTime_;

    Scalar timeStepSize_;
    int timeStepIdx_;

    bool finished_;
    bool verbose_;
};

namespace Properties {
template<class TypeTag>
struct Simulator<TypeTag, TTag::NumericModel> { using type = ::Opm::Simulator<TypeTag>; };
}

} // namespace Opm

#endif
