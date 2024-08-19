/*
  Copyright 2024 Equinor ASA.

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

#ifndef OPM_ADAPTIVE_TIME_STEPPING_IMPL_HPP
#define OPM_ADAPTIVE_TIME_STEPPING_IMPL_HPP
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#endif


namespace Opm {
/*********************************************
 * Public methods of AdaptiveTimeStepping
 * ******************************************/


//! \brief contructor taking parameter object
template<class TypeTag>
AdaptiveTimeStepping<TypeTag>::AdaptiveTimeStepping(
                                 const UnitSystem& unit_system,
                                 const double max_next_tstep,
                                 const bool terminal_output
)
    : time_step_control_{}
    , restart_factor_{Parameters::Get<Parameters::SolverRestartFactor<Scalar>>()} // 0.33
    , growth_factor_{Parameters::Get<Parameters::SolverGrowthFactor<Scalar>>()} // 2.0
    , max_growth_{Parameters::Get<Parameters::SolverMaxGrowth<Scalar>>()} // 3.0
    , max_time_step_{
        Parameters::Get<Parameters::SolverMaxTimeStepInDays<Scalar>>() * 24 * 60 * 60} // 365.25
    , min_time_step_{
        unit_system.to_si(UnitSystem::measure::time,
                          Parameters::Get<Parameters::SolverMinTimeStep<Scalar>>())} // 1e-12;
    , ignore_convergence_failure_{
        Parameters::Get<Parameters::SolverContinueOnConvergenceFailure>()} // false;
    , solver_restart_max_{Parameters::Get<Parameters::SolverMaxRestarts>()} // 10
    , solver_verbose_{Parameters::Get<Parameters::SolverVerbosity>() > 0 && terminal_output} // 2
    , timestep_verbose_{Parameters::Get<Parameters::TimeStepVerbosity>() > 0 && terminal_output} // 2
    , suggested_next_timestep_{
        (max_next_tstep <= 0 ? Parameters::Get<Parameters::InitialTimeStepInDays>()
                             : max_next_tstep) * 24 * 60 * 60} // 1.0
    , full_timestep_initially_{Parameters::Get<Parameters::FullTimeStepInitially>()} // false
    , timestep_after_event_{
        Parameters::Get<Parameters::TimeStepAfterEventInDays<Scalar>>() * 24 * 60 * 60} // 1e30
    , use_newton_iteration_{false}
    , min_time_step_before_shutting_problematic_wells_{
        Parameters::Get<Parameters::MinTimeStepBeforeShuttingProblematicWellsInDays>() * unit::day}

{
    init_(unit_system);
}

//! \brief contructor
//! \param tuning Pointer to ecl TUNING keyword
template<class TypeTag>
AdaptiveTimeStepping<TypeTag>::AdaptiveTimeStepping(double max_next_tstep,
                                                 const Tuning& tuning,
                                                 const UnitSystem& unit_system,
                                                 const bool terminal_output
)
    : time_step_control_{}
    , restart_factor_{tuning.TSFCNV}
    , growth_factor_{tuning.TFDIFF}
    , max_growth_{tuning.TSFMAX}
    , max_time_step_{tuning.TSMAXZ} // 365.25
    , min_time_step_{tuning.TSFMIN} // 0.1;
    , ignore_convergence_failure_{true}
    , solver_restart_max_{Parameters::Get<Parameters::SolverMaxRestarts>()} // 10
    , solver_verbose_{Parameters::Get<Parameters::SolverVerbosity>() > 0 && terminal_output} // 2
    , timestep_verbose_{Parameters::Get<Parameters::TimeStepVerbosity>() > 0 && terminal_output} // 2
    , suggested_next_timestep_{
        max_next_tstep <= 0 ? Parameters::Get<Parameters::InitialTimeStepInDays>() * 24 * 60 * 60
                            : max_next_tstep} // 1.0
    , full_timestep_initially_{Parameters::Get<Parameters::FullTimeStepInitially>()} // false
    , timestep_after_event_{tuning.TMAXWC} // 1e30
    , use_newton_iteration_{false}
    , min_time_step_before_shutting_problematic_wells_{
        Parameters::Get<Parameters::MinTimeStepBeforeShuttingProblematicWellsInDays>() * unit::day}
{
    init_(unit_system);
}

template<class TypeTag>
bool
AdaptiveTimeStepping<TypeTag>::
operator==(const AdaptiveTimeStepping<TypeTag>& rhs)
{
    if (this->time_step_control_type_ != rhs.time_step_control_type_ ||
        (this->time_step_control_ && !rhs.time_step_control_) ||
        (!this->time_step_control_ && rhs.time_step_control_)) {
        return false;
    }

    bool result = false;
    switch (this->time_step_control_type_) {
    case TimeStepControlType::HardCodedTimeStep:
        result = castAndComp<HardcodedTimeStepControl>(rhs);
        break;
    case TimeStepControlType::PIDAndIterationCount:
        result = castAndComp<PIDAndIterationCountTimeStepControl>(rhs);
        break;
    case TimeStepControlType::SimpleIterationCount:
        result = castAndComp<SimpleIterationCountTimeStepControl>(rhs);
        break;
    case TimeStepControlType::PID:
        result = castAndComp<PIDTimeStepControl>(rhs);
        break;
    }

    return result &&
           this->restart_factor_ == rhs.restart_factor_ &&
           this->growth_factor_ == rhs.growth_factor_ &&
           this->max_growth_ == rhs.max_growth_ &&
           this->max_time_step_ == rhs.max_time_step_ &&
           this->min_time_step_ == rhs.min_time_step_ &&
           this->ignore_convergence_failure_ == rhs.ignore_convergence_failure_ &&
           this->solver_restart_max_== rhs.solver_restart_max_ &&
           this->solver_verbose_ == rhs.solver_verbose_ &&
           this->full_timestep_initially_ == rhs.full_timestep_initially_ &&
           this->timestep_after_event_ == rhs.timestep_after_event_ &&
           this->use_newton_iteration_ == rhs.use_newton_iteration_ &&
           this->min_time_step_before_shutting_problematic_wells_ ==
               rhs.min_time_step_before_shutting_problematic_wells_;
}

template<class TypeTag>
void
AdaptiveTimeStepping<TypeTag>::
registerParameters()
{
    registerEclTimeSteppingParameters<Scalar>();
    detail::registerAdaptiveParameters();
}

template<class TypeTag>
void
AdaptiveTimeStepping<TypeTag>::
setReservoirCouplingMaster(ReservoirCouplingMaster *reservoir_coupling_master)
{
    this->reservoir_coupling_master_ = reservoir_coupling_master;
}

template<class TypeTag>
void
AdaptiveTimeStepping<TypeTag>::
setReservoirCouplingSlave(ReservoirCouplingSlave *reservoir_coupling_slave)
{
    this->reservoir_coupling_slave_ = reservoir_coupling_slave;
}

/** \brief  step method that acts like the solver::step method
            in a sub cycle of time steps
    \param tuningUpdater Function used to update TUNING parameters before each
                         time step. ACTIONX might change tuning.
*/
template<class TypeTag>
template <class Solver>
SimulatorReport
AdaptiveTimeStepping<TypeTag>::
step(const SimulatorTimer& simulator_timer,
                                    Solver& solver,
                                    const bool is_event,
                                    const std::vector<int>* fipnum,
                                    const std::function<bool()> tuning_updater)
{
    SubStepper sub_stepper{
        *this, simulator_timer, solver, is_event, fipnum, tuning_updater,
    };
    return sub_stepper.run();
}

template<class TypeTag>
template<class Serializer>
void
AdaptiveTimeStepping<TypeTag>::
serializeOp(Serializer& serializer)
{
    serializer(this->time_step_control_type_);
    switch (this->time_step_control_type_) {
    case TimeStepControlType::HardCodedTimeStep:
        allocAndSerialize<HardcodedTimeStepControl>(serializer);
        break;
    case TimeStepControlType::PIDAndIterationCount:
        allocAndSerialize<PIDAndIterationCountTimeStepControl>(serializer);
        break;
    case TimeStepControlType::SimpleIterationCount:
        allocAndSerialize<SimpleIterationCountTimeStepControl>(serializer);
        break;
    case TimeStepControlType::PID:
        allocAndSerialize<PIDTimeStepControl>(serializer);
        break;
    }
    serializer(this->restart_factor_);
    serializer(this->growth_factor_);
    serializer(this->max_growth_);
    serializer(this->max_time_step_);
    serializer(this->min_time_step_);
    serializer(this->ignore_convergence_failure_);
    serializer(this->solver_restart_max_);
    serializer(this->solver_verbose_);
    serializer(this->timestep_verbose_);
    serializer(this->suggested_next_timestep_);
    serializer(this->full_timestep_initially_);
    serializer(this->timestep_after_event_);
    serializer(this->use_newton_iteration_);
    serializer(this->min_time_step_before_shutting_problematic_wells_);
}

template<class TypeTag>
AdaptiveTimeStepping<TypeTag>
AdaptiveTimeStepping<TypeTag>::
serializationTestObjectHardcoded()
{
    return serializationTestObject_<HardcodedTimeStepControl>();
}

template<class TypeTag>
AdaptiveTimeStepping<TypeTag>
AdaptiveTimeStepping<TypeTag>::
serializationTestObjectPID()
{
    return serializationTestObject_<PIDTimeStepControl>();
}

template<class TypeTag>
AdaptiveTimeStepping<TypeTag>
AdaptiveTimeStepping<TypeTag>::
serializationTestObjectPIDIt()
{
    return serializationTestObject_<PIDAndIterationCountTimeStepControl>();
}

template<class TypeTag>
AdaptiveTimeStepping<TypeTag>
AdaptiveTimeStepping<TypeTag>::
serializationTestObjectSimple()
{
    return serializationTestObject_<SimpleIterationCountTimeStepControl>();
}


template<class TypeTag>
void
AdaptiveTimeStepping<TypeTag>::
setSuggestedNextStep(const double x)
{
    this->suggested_next_timestep_ = x;
}

template<class TypeTag>
double
AdaptiveTimeStepping<TypeTag>::
suggestedNextStep() const
{
    return this->suggested_next_timestep_;
}


template<class TypeTag>
void
AdaptiveTimeStepping<TypeTag>::
updateNEXTSTEP(double max_next_tstep)
{
     // \Note Only update next suggested step if TSINIT was explicitly
     //  set in TUNING or NEXTSTEP is active.
    if (max_next_tstep > 0) {
        this->suggested_next_timestep_ = max_next_tstep;
    }
}

template<class TypeTag>
void
AdaptiveTimeStepping<TypeTag>::
updateTUNING(double max_next_tstep, const Tuning& tuning)
{
    this->restart_factor_ = tuning.TSFCNV;
    this->growth_factor_ = tuning.TFDIFF;
    this->max_growth_ = tuning.TSFMAX;
    this->max_time_step_ = tuning.TSMAXZ;
    updateNEXTSTEP(max_next_tstep);
    this->timestep_after_event_ = tuning.TMAXWC;
}

/*********************************************
 * Private methods of AdaptiveTimeStepping
 * ******************************************/

template<class TypeTag>
template<class T, class Serializer>
void
AdaptiveTimeStepping<TypeTag>::
allocAndSerialize(Serializer& serializer)
{
    if (!serializer.isSerializing()) {
        this->time_step_control_ = std::make_unique<T>();
    }
    serializer(*static_cast<T*>(this->time_step_control_.get()));
}

template<class TypeTag>
template<class T>
bool
AdaptiveTimeStepping<TypeTag>::
castAndComp(const AdaptiveTimeStepping<TypeTag>& Rhs) const
{
    const T* lhs = static_cast<const T*>(this->time_step_control_.get());
    const T* rhs = static_cast<const T*>(Rhs.time_step_control_.get());
    return *lhs == *rhs;
}

template<class TypeTag>
void
AdaptiveTimeStepping<TypeTag>::
maybeModifySuggestedTimeStepAtBeginningOfReportStep_(const double original_time_step,
                                                     bool is_event)
{
    // init last time step as a fraction of the given time step
    if (this->suggested_next_timestep_ < 0) {
        this->suggested_next_timestep_ = this->restart_factor_ * original_time_step;
    }

    if (this->full_timestep_initially_) {
        this->suggested_next_timestep_ = original_time_step;
    }

    // use seperate time step after event
    if (is_event && this->timestep_after_event_ > 0) {
        this->suggested_next_timestep_ = this->timestep_after_event_;
    }
}

template<class TypeTag>
template<class Controller>
AdaptiveTimeStepping<TypeTag>
AdaptiveTimeStepping<TypeTag>::
serializationTestObject_()
{
    AdaptiveTimeStepping<TypeTag> result;

    result.restart_factor_ = 1.0;
    result.growth_factor_ = 2.0;
    result.max_growth_ = 3.0;
    result.max_time_step_ = 4.0;
    result.min_time_step_ = 5.0;
    result.ignore_convergence_failure_ = true;
    result.solver_restart_max_ = 6;
    result.solver_verbose_ = true;
    result.timestep_verbose_ = true;
    result.suggested_next_timestep_ = 7.0;
    result.full_timestep_initially_ = true;
    result.use_newton_iteration_ = true;
    result.min_time_step_before_shutting_problematic_wells_ = 9.0;
    result.time_step_control_type_ = Controller::Type;
    result.time_step_control_ =
                       std::make_unique<Controller>(Controller::serializationTestObject());

    return result;
}

/*********************************************
 * Protected methods of AdaptiveTimeStepping
 * ******************************************/

template<class TypeTag>
void AdaptiveTimeStepping<TypeTag>::
init_(const UnitSystem& unitSystem)
{
    // valid are "pid" and "pid+iteration"
    std::string control = Parameters::Get<Parameters::TimeStepControl>(); // "pid"

    const double tol =  Parameters::Get<Parameters::TimeStepControlTolerance>(); // 1e-1
    if (control == "pid") {
        this->time_step_control_ = std::make_unique<PIDTimeStepControl>(tol);
        this->time_step_control_type_ = TimeStepControlType::PID;
    }
    else if (control == "pid+iteration") {
        const int iterations =  Parameters::Get<Parameters::TimeStepControlTargetIterations>(); // 30
        const double decay_damping_factor =
                                Parameters::Get<Parameters::TimeStepControlDecayDampingFactor>(); // 1.0
        const double growth_damping_factor =
                                Parameters::Get<Parameters::TimeStepControlGrowthDampingFactor>(); // 3.2
        this->time_step_control_ = std::make_unique<PIDAndIterationCountTimeStepControl>(
                                iterations, decay_damping_factor, growth_damping_factor, tol);
        this->time_step_control_type_ = TimeStepControlType::PIDAndIterationCount;
    }
    else if (control == "pid+newtoniteration") {
        const int iterations =
            Parameters::Get<Parameters::TimeStepControlTargetNewtonIterations>(); // 8
        const double decay_damping_factor =
            Parameters::Get<Parameters::TimeStepControlDecayDampingFactor>(); // 1.0
        const double growth_damping_factor =
            Parameters::Get<Parameters::TimeStepControlGrowthDampingFactor>(); // 3.2
        const double non_dimensional_min_time_step_iterations =
            Parameters::Get<Parameters::MinTimeStepBasedOnNewtonIterations>(); // 0.0 by default
        // the min time step can be reduced by the newton iteration numbers
        double min_time_step_reduced_by_iterations =
            unitSystem.to_si(UnitSystem::measure::time, non_dimensional_min_time_step_iterations);
        this->time_step_control_ =
            std::make_unique<PIDAndIterationCountTimeStepControl>(
                iterations, decay_damping_factor, growth_damping_factor,
                    tol, min_time_step_reduced_by_iterations);
        this->time_step_control_type_ = TimeStepControlType::PIDAndIterationCount;
        this->use_newton_iteration_ = true;
    }
    else if (control == "iterationcount") {
        const int iterations =  Parameters::Get<Parameters::TimeStepControlTargetIterations>(); // 30
        const double decay_rate = Parameters::Get<Parameters::TimeStepControlDecayRate>(); // 0.75
        const double growth_rate = Parameters::Get<Parameters::TimeStepControlGrowthRate>(); // 1.25
        this->time_step_control_ =
            std::make_unique<SimpleIterationCountTimeStepControl>(iterations, decay_rate, growth_rate);
        this->time_step_control_type_ = TimeStepControlType::SimpleIterationCount;
    }
    else if (control == "newtoniterationcount") {
        const int iterations =  Parameters::Get<Parameters::TimeStepControlTargetNewtonIterations>(); // 8
        const double decay_rate = Parameters::Get<Parameters::TimeStepControlDecayRate>(); // 0.75
        const double growth_rate = Parameters::Get<Parameters::TimeStepControlGrowthRate>(); // 1.25
        this->time_step_control_ =
            std::make_unique<SimpleIterationCountTimeStepControl>(iterations, decay_rate, growth_rate);
        this->use_newton_iteration_ = true;
        this->time_step_control_type_ = TimeStepControlType::SimpleIterationCount;
    }
    else if (control == "hardcoded") {
        const std::string filename = Parameters::Get<Parameters::TimeStepControlFileName>(); // "timesteps"
        this->time_step_control_ = std::make_unique<HardcodedTimeStepControl>(filename);
        this->time_step_control_type_ = TimeStepControlType::HardCodedTimeStep;
    }
    else
        OPM_THROW(std::runtime_error,
                  "Unsupported time step control selected " + control);

    // make sure growth factor is something reasonable
    assert(this->growth_factor_ >= 1.0);
}



/************************************************
 * Private class SubStepper public methods
 ************************************************/

template<class TypeTag>
template<class Solver>
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
SubStepper(
    AdaptiveTimeStepping<TypeTag>& adaptive_time_stepping,
    const SimulatorTimer& simulator_timer,
    Solver& solver,
    const bool is_event,
    const std::vector<int>* fipnum,
    const std::function<bool()>& tuning_updater
)
    : adaptive_time_stepping_{adaptive_time_stepping}
    , simulator_timer_{simulator_timer}
    , solver_{solver}
    , is_event_{is_event}
    , fipnum_{fipnum}
    , tuning_updater_{tuning_updater}
{
}

template<class TypeTag>
template<class Solver>
AdaptiveTimeStepping<TypeTag>&
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
getAdaptiveTimerStepper()
{
    return adaptive_time_stepping_;
}

template<class TypeTag>
template<class Solver>
SimulatorReport
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
run()
{
#if HAVE_MPI
    if (isReservoirCouplingSlave_()) {
        return runStepReservoirCouplingSlave_();
    }
    else if (isReservoirCouplingMaster_()) {
        return runStepReservoirCouplingMaster_();
    }
    else {
        return runStepOriginal_();
    }
#else
    return runStepOriginal_();
#endif
}

/************************************************
 * Private class SubStepper private methods
 ************************************************/


template<class TypeTag>
template<class Solver>
bool
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
isReservoirCouplingMaster_() const
{
    return adaptive_time_stepping_.reservoir_coupling_master_ != nullptr;
}

template<class TypeTag>
template<class Solver>
bool
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
isReservoirCouplingSlave_() const
{
    return adaptive_time_stepping_.reservoir_coupling_slave_ != nullptr;
}

template<class TypeTag>
template<class Solver>
void
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
maybeModifySuggestedTimeStepAtBeginningOfReportStep_(const double original_time_step)
{
    this->adaptive_time_stepping_.maybeModifySuggestedTimeStepAtBeginningOfReportStep_(
        original_time_step, this->is_event_
    );
}

template<class TypeTag>
template<class Solver>
bool
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
maybeUpdateTuning_() const
{
    return this->tuning_updater_();
}

template<class TypeTag>
template<class Solver>
double
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
maxTimeStep_() const
{
    return this->adaptive_time_stepping_.max_time_step_;
}

template <class TypeTag>
template <class Solver>
SimulatorReport
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
runStepOriginal_()
{
    maybeUpdateTuning_();
    const double original_time_step = this->simulator_timer_.currentStepLength();
    maybeModifySuggestedTimeStepAtBeginningOfReportStep_(original_time_step);

    AdaptiveSimulatorTimer substep_timer{
        this->simulator_timer_.startDateTime(),
        original_time_step,
        this->simulator_timer_.simulationTimeElapsed(),
        suggestedNextTimestep_(),
        this->simulator_timer_.reportStepNum(),
        maxTimeStep_()
    };
    SubStepIteration substepIteration{*this, substep_timer, original_time_step, /*final_step=*/true};
    return substepIteration.run();
}

#if HAVE_MPI
template <class TypeTag>
template <class Solver>
ReservoirCouplingMaster&
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
reservoirCouplingMaster_()
{
    return *adaptive_time_stepping_.reservoir_coupling_master_;
}
#endif

#if HAVE_MPI
template <class TypeTag>
template <class Solver>
ReservoirCouplingSlave&
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
reservoirCouplingSlave_()
{
    return *adaptive_time_stepping_.reservoir_coupling_slave_;
}
#endif

#if HAVE_MPI
template <class TypeTag>
template <class Solver>
SimulatorReport
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
runStepReservoirCouplingMaster_()
{
    bool substep_done = false;
    int iteration = 0;
    const double original_time_step = this->simulator_timer_.currentStepLength();
    TimePoint current_time{this->simulator_timer_.simulationTimeElapsed()};
    TimePoint step_end_time = current_time + original_time_step;
    auto current_step_length = original_time_step;
    SimulatorReport report;
    while(!substep_done) {
        reservoirCouplingMaster_().receiveNextReportDateFromSlaves();
        if (iteration == 0) {
            maybeUpdateTuning_();
        }
        current_step_length = reservoirCouplingMaster_().maybeChopSubStep(
                                          current_step_length, current_time.getTime());
        reservoirCouplingMaster_().sendNextTimeStepToSlaves(current_step_length);
        if (iteration == 0) {
            maybeModifySuggestedTimeStepAtBeginningOfReportStep_(current_step_length);
        }
        AdaptiveSimulatorTimer substep_timer{
            this->simulator_timer_.startDateTime(),
            /*stepLength=*/current_step_length,
            /*elapsedTime=*/current_time.getTime(),
            /*timeStepEstimate=*/suggestedNextTimestep_(),
            this->simulator_timer_.reportStepNum(),
            maxTimeStep_()
        };
        bool final_step = (current_time + current_step_length) >= step_end_time;
        SubStepIteration substepIteration{*this, substep_timer, current_step_length, final_step};
        auto sub_steps_report = substepIteration.run();
        report += sub_steps_report;
        current_time += current_step_length;
        if (final_step) {
            break;
        }
        iteration++;
    }
    return report;
}
#endif

#if HAVE_MPI
template <class TypeTag>
template <class Solver>
SimulatorReport
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
runStepReservoirCouplingSlave_()
{
    bool substep_done = false;
    int iteration = 0;
    const double original_time_step = this->simulator_timer_.currentStepLength();
    TimePoint current_time{this->simulator_timer_.simulationTimeElapsed()};
    TimePoint step_end_time = current_time + original_time_step;
    SimulatorReport report;
    while(!substep_done) {
        reservoirCouplingSlave_().sendNextReportDateToMasterProcess();
        auto timestep = reservoirCouplingSlave_().receiveNextTimeStepFromMaster();
        if (iteration == 0) {
            maybeUpdateTuning_();
            maybeModifySuggestedTimeStepAtBeginningOfReportStep_(timestep);
        }
        AdaptiveSimulatorTimer substep_timer{
            this->simulator_timer_.startDateTime(),
            /*step_length=*/timestep,
            /*elapsed_time=*/current_time.getTime(),
            /*time_step_estimate=*/suggestedNextTimestep_(),
            this->simulator_timer_.reportStepNum(),
            maxTimeStep_()
        };
        bool final_step = (current_time + timestep) >= step_end_time;
        SubStepIteration substepIteration{*this, substep_timer, timestep, final_step};
        auto sub_steps_report = substepIteration.run();
        report += sub_steps_report;
        current_time += timestep;
        if (final_step) {
            break;
        }
        iteration++;
    }
    return report;
}

#endif

template <class TypeTag>
template <class Solver>
double
AdaptiveTimeStepping<TypeTag>::SubStepper<Solver>::
suggestedNextTimestep_() const
{
    return this->adaptive_time_stepping_.suggestedNextStep();
}



/************************************************
 * Private class SubStepIteration public methods
 ************************************************/

template<class TypeTag>
template<class Solver>
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
SubStepIteration(
    SubStepper<Solver>& substepper,
    AdaptiveSimulatorTimer& substep_timer,
    const double original_time_step,
    bool final_step
)
    : substepper_{substepper}
    , substep_timer_{substep_timer}
    , original_time_step_{original_time_step}
    , final_step_{final_step}
    , adaptive_time_stepping_{substepper.getAdaptiveTimerStepper()}
{
}

template <class TypeTag>
template <class Solver>
SimulatorReport
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
run()
{
    auto& simulator = solver_().model().simulator();
    auto& problem = simulator.problem();
    // counter for solver restarts
    int restarts = 0;
    SimulatorReport report;

    // sub step time loop
    while (!this->substep_timer_.done()) {
        auto old_value = suggestedNextTimestep_();
        if (maybeUpdateTuning_()) {
            setTimeStep_(suggestedNextTimestep_());
            // maybeUpdateTuning() might change the suggested time step
            setSuggestedNextStep_(old_value);
        }
        const double dt = this->substep_timer_.currentStepLength();
        if (timeStepVerbose_()) {
            detail::logTimer(this->substep_timer_);
        }

        auto substep_report = runSubStep_();

        //Pass substep to eclwriter for summary output
        problem.setSubStepReport(substep_report);

        report += substep_report;

        if (substep_report.converged || checkContinueOnUnconvergedSolution_(dt)) {
            ++this->substep_timer_;   // advance by current dt

            const int iterations = getNumIterations_(substep_report);
            auto dt_estimate = timeStepControlComputeEstimate_(
                                     dt, iterations, this->substep_timer_.simulationTimeElapsed());

            assert(dt_estimate > 0);
            dt_estimate = maybeRestrictTimeStepGrowth_(dt, dt_estimate, restarts);
            restarts = 0;         // solver converged, reset restarts counter

            maybeReportSubStep_(substep_report);
            if (this->final_step_ && this->substep_timer_.done()) {
                // if the time step is done we do not need to write it as this will be done
                // by the simulator anyway.
            }
            else {
                report.success.output_write_time += writeOutput_();
            }

            // set new time step length
            setTimeStep_(dt_estimate);

            report.success.converged = this->substep_timer_.done();
            this->substep_timer_.setLastStepFailed(false);
        }
        else { // in case of no convergence
            this->substep_timer_.setLastStepFailed(true);
            checkTimeStepMaxRestartLimit_(restarts);

            const double new_time_step = restartFactor_() * dt;
            checkTimeStepMinLimit_(new_time_step);
            bool wells_shut = false;
            if (new_time_step > minTimeStepBeforeClosingWells_()) {
                chopTimeStep_(new_time_step);
            } else {
                wells_shut = chopTimeStepOrCloseFailingWells_(new_time_step);
            }
            if (wells_shut) {
                setTimeStep_(dt);  // retry the old timestep
            }
            else {
                restarts++;   // only increase if no wells were shut
            }
        }
        problem.setNextTimeStepSize(this->substep_timer_.currentStepLength());
    }
    updateSuggestedNextStep_();
    return report;
}


/************************************************
 * Private class SubStepIteration private methods
 ************************************************/


template<class TypeTag>
template<class Solver>
bool
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
checkContinueOnUnconvergedSolution_(double dt) const
{
    bool continue_on_uncoverged_solution = ignoreConvergenceFailure_() && dt <= minTimeStep_();
    if (continue_on_uncoverged_solution && solverVerbose_()) {
        // NOTE: This method is only called if the solver failed to converge.
        const auto msg = fmt::format(
            "Solver failed to converge but timestep {} is smaller or equal to {}\n"
            "which is the minimum threshold given by option --solver-min-time-step\n",
            dt, minTimeStep_()
        );
        OpmLog::problem(msg);
    }
    return continue_on_uncoverged_solution;
}

template<class TypeTag>
template<class Solver>
void
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
checkTimeStepMaxRestartLimit_(const int restarts) const
{
   // If we have restarted (i.e. cut the timestep) too
    // many times, we have failed and throw an exception.
    if (restarts >= solverRestartMax_()) {
        const auto msg = fmt::format(
            "Solver failed to converge after cutting timestep {} times.", restarts
        );
        if (solverVerbose_()) {
            OpmLog::error(msg);
        }
        // Use throw directly to prevent file and line
        throw TimeSteppingBreakdown{msg};
    }
}

template<class TypeTag>
template<class Solver>
void
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
checkTimeStepMinLimit_(const int new_time_step) const
{
    // If we have restarted (i.e. cut the timestep) too
    // much, we have failed and throw an exception.
    if (new_time_step < minTimeStep_()) {
        const auto msg = fmt::format(
            "Solver failed to converge after cutting timestep to {}\n"
            "which is the minimum threshold given by option --solver-min-time-step\n",
            minTimeStep_()
        );
        if (solverVerbose_()) {
            OpmLog::error(msg);
        }
        // Use throw directly to prevent file and line
        throw TimeSteppingBreakdown{msg};
    }
}

template<class TypeTag>
template<class Solver>
void
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
chopTimeStep_(const double new_time_step)
{
    setTimeStep_(new_time_step);
    if (solverVerbose_()) {
        const auto msg = fmt::format("{}\nTimestep chopped to {} days\n",
                    this->cause_of_failure_,
                    std::to_string(unit::convert::to(
                        this->substep_timer_.currentStepLength(), unit::day)));
        OpmLog::problem(msg);
    }
}

template<class TypeTag>
template<class Solver>
bool
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
chopTimeStepOrCloseFailingWells_(const int new_time_step)
{
    bool wells_shut = false;
    // We are below the threshold, and will check if there are any wells we should close
    // rather than chopping again.
    std::set<std::string> failing_wells = detail::consistentlyFailingWells(
                                         solver_().model().stepReports());
    if (failing_wells.empty()) {
        // Found no wells to close, chop the timestep
        chopTimeStep_(new_time_step);
    } else {
        // Close all consistently failing wells.
        int num_shut_wells = 0;
        for (const auto& well : failing_wells) {
            bool was_shut = solver_().model().wellModel().forceShutWellByName(
                                              well, this->substep_timer_.simulationTimeElapsed());
            if (was_shut) {
                ++num_shut_wells;
            }
        }
        if (num_shut_wells == 0) {
            // None of the problematic wells were shut.
            // We must fall back to chopping again.
            chopTimeStep_(new_time_step);
        } else {
            wells_shut = true;
            if (solverVerbose_()) {
                std::string msg;
                msg = "\nProblematic well(s) were shut: ";
                for (const auto& well : failing_wells) {
                    msg += well;
                    msg += " ";
                }
                msg += "(retrying timestep)\n";
                OpmLog::problem(msg);
            }
        }
    }
    return wells_shut;
}

template<class TypeTag>
template<class Solver>
boost::posix_time::ptime
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
currentDateTime_() const
{
    return simulatorTimer_().currentDateTime();
}

template<class TypeTag>
template<class Solver>
const std::vector<int> *
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
fipnum_() const
{
    return this->substepper_.fipnum_;
}

template<class TypeTag>
template<class Solver>
int
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
getNumIterations_(const SimulatorReportSingle &substep_report) const
{
    if (useNewtonIteration_()) {
        return substep_report.total_newton_iterations;
    }
    else {
        return substep_report.total_linear_iterations;
    }
}

template<class TypeTag>
template<class Solver>
double
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
growthFactor_() const
{
    return this->adaptive_time_stepping_.growth_factor_;
}

template<class TypeTag>
template<class Solver>
bool
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
ignoreConvergenceFailure_() const
{
    return adaptive_time_stepping_.ignore_convergence_failure_;
}

template<class TypeTag>
template<class Solver>
double
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
maxGrowth_() const
{
    return this->adaptive_time_stepping_.max_growth_;
}

template<class TypeTag>
template<class Solver>
void
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
maybeReportSubStep_(SimulatorReportSingle substep_report) const
{
    if (timeStepVerbose_()) {
        std::ostringstream ss;
        substep_report.reportStep(ss);
        OpmLog::info(ss.str());
    }
}

template<class TypeTag>
template<class Solver>
double
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
maybeRestrictTimeStepGrowth_(const double dt, double dt_estimate, const int restarts) const
{
    // limit the growth of the timestep size by the growth factor
    dt_estimate = std::min(dt_estimate, double(maxGrowth_() * dt));
    assert(dt_estimate > 0);
    // further restrict time step size growth after convergence problems
    if (restarts > 0) {
        dt_estimate = std::min(growthFactor_() * dt, dt_estimate);
    }

    // Further restrict time step size if we are in
    // prediction mode with THP constraints.
    if (solver_().model().wellModel().hasTHPConstraints()) {
        const double max_prediction_thp_timestep = 16.0 * unit::day;
        dt_estimate = std::min(dt_estimate, max_prediction_thp_timestep);
    }
    assert(dt_estimate > 0);
    return dt_estimate;
}

template<class TypeTag>
template<class Solver>
bool
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
maybeUpdateTuning_() const
{
    return this->substepper_.maybeUpdateTuning_();
}

template<class TypeTag>
template<class Solver>
double
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
minTimeStepBeforeClosingWells_() const
{
    return this->adaptive_time_stepping_.min_time_step_before_shutting_problematic_wells_;
}

template<class TypeTag>
template<class Solver>
double
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
minTimeStep_() const
{
    return this->adaptive_time_stepping_.min_time_step_;
}

template<class TypeTag>
template<class Solver>
double
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
restartFactor_() const
{
    return this->adaptive_time_stepping_.restart_factor_;
}

template<class TypeTag>
template<class Solver>
SimulatorReportSingle
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
runSubStep_()
{
    SimulatorReportSingle substep_report;

    auto handleFailure = [this, &substep_report]
            (const std::string& failure_reason, const std::exception& e, bool log_exception = true)
    {
        substep_report = solver_().failureReport();
        this->cause_of_failure_ = failure_reason;
        if (log_exception && solverVerbose_()) {
            std::string message;
            message = "Caught Exception: ";
            message += e.what();
            OpmLog::debug(message);
        }
    };

    try {
        substep_report = solver_().step(this->substep_timer_);
        if (solverVerbose_()) {
            // report number of linear iterations
            OpmLog::debug("Overall linear iterations used: "
                          + std::to_string(substep_report.total_linear_iterations));
        }
    }
    catch (const TooManyIterations& e) {
        handleFailure("Solver convergence failure - Iteration limit reached", e);
    }
    catch (const ConvergenceMonitorFailure& e) {
        handleFailure("Convergence monitor failure", e, /*log_exception=*/false);
    }
    catch (const LinearSolverProblem& e) {
        handleFailure("Linear solver convergence failure", e);
    }
    catch (const NumericalProblem& e) {
        handleFailure("Solver convergence failure - Numerical problem encountered", e);
    }
    catch (const std::runtime_error& e) {
        handleFailure("Runtime error encountered", e);
    }
    catch (const Dune::ISTLError& e) {
        handleFailure("ISTL error - Time step too large", e);
    }
    catch (const Dune::MatrixBlockError& e) {
        handleFailure("Matrix block error", e);
    }

    return substep_report;
}

template<class TypeTag>
template<class Solver>
void
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
setTimeStep_(double dt_estimate)
{
    this->substep_timer_.provideTimeStepEstimate(dt_estimate);
}

template<class TypeTag>
template<class Solver>
Solver&
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
solver_() const
{
    return this->substepper_.solver_;
}


template<class TypeTag>
template<class Solver>
int
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
solverRestartMax_() const
{
    return this->adaptive_time_stepping_.solver_restart_max_;
}

template<class TypeTag>
template<class Solver>
void
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
setSuggestedNextStep_(double step)
{
    this->adaptive_time_stepping_.setSuggestedNextStep(step);
}

template <class TypeTag>
template <class Solver>
const SimulatorTimer&
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
simulatorTimer_() const
{
    return this->substepper_.simulator_timer_;
}

template <class TypeTag>
template <class Solver>
bool
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
solverVerbose_() const
{
    return this->adaptive_time_stepping_.solver_verbose_;
}

template<class TypeTag>
template<class Solver>
boost::posix_time::ptime
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
startDateTime_() const
{
    return simulatorTimer_().startDateTime();
}

template <class TypeTag>
template <class Solver>
double
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
suggestedNextTimestep_() const
{
    return this->adaptive_time_stepping_.suggestedNextStep();
}

template <class TypeTag>
template <class Solver>
double
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
timeStepControlComputeEstimate_(const double dt, const int iterations, double elapsed) const
{
    // create object to compute the time error, simply forwards the call to the model
    SolutionTimeErrorSolverWrapper<Solver> relative_change{solver_()};
    return this->adaptive_time_stepping_.time_step_control_->computeTimeStepSize(
        dt, iterations, relative_change, elapsed);
}

template <class TypeTag>
template <class Solver>
bool
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
timeStepVerbose_() const
{
    return this->adaptive_time_stepping_.timestep_verbose_;
}

template <class TypeTag>
template <class Solver>
void
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
updateSuggestedNextStep_()
{
    auto suggested_next_step = this->substep_timer_.currentStepLength();
    if (! std::isfinite(suggested_next_step)) { // check for NaN
        suggested_next_step = this->original_time_step_;
    }
    if (timeStepVerbose_()) {
        std::ostringstream ss;
        this->substep_timer_.report(ss);
        ss << "Suggested next step size = "
           << unit::convert::to(suggested_next_step, unit::day) << " (days)" << std::endl;
        OpmLog::debug(ss.str());
    }
    setSuggestedNextStep_(suggested_next_step);
}

template <class TypeTag>
template <class Solver>
bool
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
useNewtonIteration_() const
{
    return this->adaptive_time_stepping_.use_newton_iteration_;
}

template <class TypeTag>
template <class Solver>
double
AdaptiveTimeStepping<TypeTag>::SubStepIteration<Solver>::
writeOutput_() const
{
    if (fipnum_()) {
        solver_().computeFluidInPlace(*(fipnum_()));
    }
    time::StopWatch perf_timer;
    perf_timer.start();
    auto& problem = solver_().model().simulator().problem();
    problem.writeOutput();
    return perf_timer.secsSinceStart();
}

/************************************************
 * Private class SolutionTimeErrorSolverWrapper
 * **********************************************/

template<class TypeTag>
template<class Solver>
AdaptiveTimeStepping<TypeTag>::SolutionTimeErrorSolverWrapper<Solver>::SolutionTimeErrorSolverWrapper(
    const Solver& solver
)
    : solver_{solver}
{}

template<class TypeTag>
template<class Solver>
double AdaptiveTimeStepping<TypeTag>::SolutionTimeErrorSolverWrapper<Solver>::relativeChange() const
{
    // returns:   || u^n+1 - u^n || / || u^n+1 ||
    return solver_.model().relativeChange();
}

} // namespace Opm
