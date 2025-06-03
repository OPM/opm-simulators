/*
  Copyright 2014 IRIS AS
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 Statoil AS

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
#include <config.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>
#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>

#include <fmt/format.h>

namespace Opm
{
    ////////////////////////////////////////////////////////
    //
    //  InterationCountTimeStepControl Implementation
    //
    ////////////////////////////////////////////////////////

    SimpleIterationCountTimeStepControl::
    SimpleIterationCountTimeStepControl( const int target_iterations,
                                         const double decayrate,
                                         const double growthrate,
                                         const bool verbose)
        : target_iterations_( target_iterations )
        , decayrate_( decayrate )
        , growthrate_( growthrate )
        , verbose_( verbose )
    {
        if( decayrate_  > 1.0 ) {
            OPM_THROW(std::runtime_error,
                      "SimpleIterationCountTimeStepControl: "
                      "decay should be <= 1 " + std::to_string(decayrate_));
        }
        if( growthrate_ < 1.0 ) {
            OPM_THROW(std::runtime_error,
                      "SimpleIterationCountTimeStepControl: "
                      "growth should be >= 1 " + std::to_string(growthrate_));
        }
    }

    SimpleIterationCountTimeStepControl
    SimpleIterationCountTimeStepControl::serializationTestObject()
    {
        return {1, 1.0, 2.0, true};
    }

    double SimpleIterationCountTimeStepControl::
    computeTimeStepSize( const double dt, const int iterations, const RelativeChangeInterface& /* relativeChange */, const AdaptiveSimulatorTimer& /* substepTimer */) const
    {
        double dtEstimate = dt ;

        // reduce the time step size if we exceed the number of target iterations
        if( iterations > target_iterations_ )
        {
            // scale dtEstimate down with a given rate
            dtEstimate *= decayrate_;
        }
        // increase the time step size if we are below the number of target iterations
        else if ( iterations < target_iterations_ )
        {
            dtEstimate *= growthrate_;
        }

        return dtEstimate;
    }

    bool SimpleIterationCountTimeStepControl::
    operator==(const SimpleIterationCountTimeStepControl& ctrl) const
    {
         return this->target_iterations_ == ctrl.target_iterations_ &&
                this->decayrate_ == ctrl.decayrate_ &&
                this->growthrate_ == ctrl.growthrate_ &&
                this->verbose_ == ctrl.verbose_;
    }

    ////////////////////////////////////////////////////////
    //
    //  HardcodedTimeStepControl Implementation
    //
    ////////////////////////////////////////////////////////

    HardcodedTimeStepControl::
    HardcodedTimeStepControl(const std::string& filename)
    {
        std::ifstream infile (filename);
        if (!infile.is_open()) {
            OPM_THROW(std::runtime_error,"Incorrect or no filename is provided to the hardcodedTimeStep. Use timestep.control.filename=your_file_name");
        }
        std::string::size_type sz;
        std::string line;
        while ( std::getline(infile, line)) {
            if( line[0] != '-') { // ignore lines starting with '-'
                const double time = std::stod(line,&sz); // read the first number i.e. the actual substep time
                subStepTime_.push_back( time * unit::day );
            }

        }
    }

    HardcodedTimeStepControl HardcodedTimeStepControl::serializationTestObject()
    {
        HardcodedTimeStepControl result;
        result.subStepTime_ = {1.0, 2.0};

        return result;
    }

    double HardcodedTimeStepControl::
    computeTimeStepSize( const double /*dt */, const int /*iterations */, const RelativeChangeInterface& /* relativeChange */ , const AdaptiveSimulatorTimer& substepTimer) const
    {
        auto nextTime = std::upper_bound(subStepTime_.begin(), subStepTime_.end(), substepTimer.simulationTimeElapsed());
        return (*nextTime - substepTimer.simulationTimeElapsed());
    }

    bool HardcodedTimeStepControl::operator==(const HardcodedTimeStepControl& ctrl) const
    {
        return this->subStepTime_ == ctrl.subStepTime_;
    }



    ////////////////////////////////////////////////////////
    //
    //  PIDTimeStepControl Implementation
    //
    ////////////////////////////////////////////////////////

    PIDTimeStepControl::PIDTimeStepControl( const double tol,
                                            const bool verbose )
        : tol_( tol )
        , errors_( 3, tol_ )
        , verbose_( verbose )
    {}

    PIDTimeStepControl
    PIDTimeStepControl::serializationTestObject()
    {
        PIDTimeStepControl result(1.0, true);
        result.errors_ = {2.0, 3.0};

        return result;;
    }

    double PIDTimeStepControl::
    computeTimeStepSize( const double dt, const int /* iterations */, const RelativeChangeInterface& relChange, const AdaptiveSimulatorTimer& /* substepTimer */) const
    {
        // shift errors
        for( int i=0; i<2; ++i ) {
            errors_[ i ] = errors_[i+1];
        }

        // store new error
        const double error = relChange.relativeChange();
        errors_[ 2 ] = error;
        for( int i=0; i<2; ++i ) {
            assert(std::isfinite(errors_[i]));
        }

        if( errors_[2] > tol_ )
        {
            // adjust dt by given tolerance
            const double newDt = dt * tol_ / error;
            if ( verbose_ )
                    OpmLog::info(fmt::format("Computed step size (tol): {} days", unit::convert::to( newDt, unit::day )));
            return newDt;
        }
        else if (errors_[0] == 0 || errors_[1] == 0 || errors_[2] == 0.)
        {
            if ( verbose_ )
                OpmLog::info("The solution between time steps does not change, there is no time step constraint from the PID time step control ");
            return std::numeric_limits<double>::max();
        }
        else
        {
            // values taking from turek time stepping paper
            const double kP = 0.075 ;
            const double kI = 0.175 ;
            const double kD = 0.01 ;
            const double newDt = (dt * std::pow( errors_[ 1 ] / errors_[ 2 ], kP ) *
                                 std::pow( tol_         / errors_[ 2 ], kI ) *
                                 std::pow( errors_[1]*errors_[1]/errors_[ 0 ]/errors_[ 2 ], kD ));
            if( verbose_ )
                OpmLog::info(fmt::format("Computed step size (pow): {} days", unit::convert::to( newDt, unit::day )));
            return newDt;
        }
    }

    bool PIDTimeStepControl::operator==(const PIDTimeStepControl& ctrl) const
    {
        return this->tol_ == ctrl.tol_ &&
               this->errors_ == ctrl.errors_ &&
               this->verbose_ == ctrl.verbose_;
    }



    ////////////////////////////////////////////////////////////
    //
    //  PIDAndIterationCountTimeStepControl  Implementation
    //
    ////////////////////////////////////////////////////////////

    PIDAndIterationCountTimeStepControl::
    PIDAndIterationCountTimeStepControl( const int target_iterations,
                                         const double decayDampingFactor,
                                         const double growthDampingFactor,
                                         const double tol,
                                         const double minTimeStepBasedOnIterations,
                                         const bool verbose)
        : PIDTimeStepControl( tol, verbose )
        , target_iterations_( target_iterations )
        , decayDampingFactor_( decayDampingFactor )
        , growthDampingFactor_( growthDampingFactor )
        , minTimeStepBasedOnIterations_(minTimeStepBasedOnIterations)
    {}

    PIDAndIterationCountTimeStepControl
    PIDAndIterationCountTimeStepControl::serializationTestObject()
    {
        return PIDAndIterationCountTimeStepControl{1, 2.0, 3.0, 4.0, 5.0, true};
    }

    double PIDAndIterationCountTimeStepControl::
    computeTimeStepSize( const double dt, const int iterations, const RelativeChangeInterface& relChange,  const AdaptiveSimulatorTimer& substepTimer) const
    {
        double dtEstimatePID = PIDTimeStepControl :: computeTimeStepSize( dt, iterations, relChange, substepTimer);

        // adjust timesteps based on target iteration
        double dtEstimateIter;
        if (iterations > target_iterations_) {
            double off_target_fraction = double(iterations - target_iterations_) / target_iterations_;
            dtEstimateIter = dt / (1.0 + off_target_fraction * decayDampingFactor_);
            if (dtEstimateIter < minTimeStepBasedOnIterations_) {
                dtEstimateIter = minTimeStepBasedOnIterations_;
            }
        } else {
            double off_target_fraction = double(target_iterations_ - iterations) / target_iterations_;
            // Be a bit more careful when increasing.
            dtEstimateIter = dt * (1.0 + off_target_fraction * growthDampingFactor_);
        }

        return std::min(dtEstimatePID, dtEstimateIter);
    }

    bool PIDAndIterationCountTimeStepControl::operator==(const PIDAndIterationCountTimeStepControl& ctrl) const
    {
        return static_cast<const PIDTimeStepControl&>(*this) == ctrl &&
               this->target_iterations_ == ctrl.target_iterations_ &&
               this->decayDampingFactor_ == ctrl.decayDampingFactor_ &&
               this->growthDampingFactor_ == ctrl.growthDampingFactor_ &&
               this->minTimeStepBasedOnIterations_ == ctrl.minTimeStepBasedOnIterations_;
    }



    ////////////////////////////////////////////////////////////
    //
    //  General3rdOrderController  Implementation
    //
    ////////////////////////////////////////////////////////////

    General3rdOrderController::General3rdOrderController( const double tolerance,
                                                          const double safetyFactor,
                                                          const bool rejectCompletedStep,
                                                          const std::string& toleranceTestVersion,
                                                          const double maxReductionTimeStep,
                                                          const std::string& parameters,
                                                          const bool verbose)
        : tolerance_( tolerance )
        , safetyFactor_( safetyFactor )
        , rejectCompletedStep_( rejectCompletedStep )
        , maxReductionTimeStep_( maxReductionTimeStep )
        , verbose_( verbose )
    {
        std::istringstream parameters_stream(parameters);
        std::string value;
        int counter = 0;
        try {
            while (std::getline(parameters_stream, value, ';'))
            {
                if (counter < 3) {
                    beta_[counter] = std::stod(value);
                }
                else {
                    alpha_[counter - 3] = std::stod(value);
                }
                counter++;
            }
        }
        catch (const std::exception& e) {
            OPM_THROW(std::runtime_error, "The parameter string for the time step controller is invalid: " + parameters);
        }

        if (toleranceTestVersion == "standard") {
            toleranceTestVersion_ = ToleranceTestVersions::Standard;
        }
        else if (toleranceTestVersion == "control-error-filtering") {
            toleranceTestVersion_ = ToleranceTestVersions::ControlErrorFiltering;
        }
        else {
            OPM_THROW(std::runtime_error, "Unsupported tolerance test version given: " + toleranceTestVersion);
        }
    }

    General3rdOrderController
    General3rdOrderController::serializationTestObject()
    {
        General3rdOrderController result(1.0, 2.0, false, "standard", 3.0, "0.125;0.25;0.125;0.75;0.25", false);
        result.errors_ = {2.0, 3.0};

        return result;
    }

    double General3rdOrderController::
    computeTimeStepSize(const double dt, const int /* iterations */, const RelativeChangeInterface& /* relChange */, const AdaptiveSimulatorTimer& substepTimer) const
    {
        if (errors_[0] == 0 || errors_[1] == 0 || errors_[2] == 0)
        {
            if ( verbose_ )
            {
                OpmLog::info("The solution between time steps does not change, there is no time step constraint from the controller.");
            }
            return std::numeric_limits<double>::max();
        }

        for( int i = 0; i < 2; ++i )
        {
            assert(std::isfinite(errors_[i]));
        }

        // Use an I controller after report time steps
        if (substepTimer.currentStepNum() < 3)
        {
            controllerVersion_ = InternalControlVersions::IController;
            const double newDt = dt * timeStepFactor(errors_, timeSteps_);
            if (verbose_)
            {
                OpmLog::info(fmt::format("Computed step size (pow): {} days", unit::convert::to( newDt, unit::day )));
            }
            return newDt;
        }
        // Use the general third order controller for all other time steps
        else
        {
            controllerVersion_ = InternalControlVersions::General3rdOrder;
            const double newDt = dt * timeStepFactor(errors_, timeSteps_);
            if (verbose_)
            {
                OpmLog::info(fmt::format("Computed step size (pow): {} days", unit::convert::to( newDt, unit::day )));
            }
            return newDt;
        }
    }

    double General3rdOrderController::
    timeStepFactor(const std::array<double, 3>& errors, const std::array<double, 3>& timeSteps) const
    {
        if (controllerVersion_ == InternalControlVersions::General3rdOrder)
        {
            return std::pow(safetyFactor_ * tolerance_ / errors[2], beta_[0]) *
                   std::pow(safetyFactor_ * tolerance_ / errors[1], beta_[1]) *
                   std::pow(safetyFactor_ * tolerance_ / errors[0], beta_[2]) *
                   std::pow(timeSteps[2] / timeSteps[1], -alpha_[0]) *
                   std::pow(timeSteps[1] / timeSteps[0], -alpha_[1]);
        }
        // controllerVersion_ == InternalControlVersions::IController
        return std::pow(safetyFactor_ * tolerance_ / errors[2], 0.35);
    }

    bool General3rdOrderController::
    timeStepAccepted(const double error, const double timeStepJustCompleted) const
    {
        bool acceptTimeStep = true;

        // Reject time step if chosen tolerance test version fails
        if (toleranceTestVersion_ == ToleranceTestVersions::Standard)
        {
            if (rejectCompletedStep_ && error > tolerance_)
            {
                acceptTimeStep = false;
            }
        }
        else if (toleranceTestVersion_ == ToleranceTestVersions::ControlErrorFiltering)
        {
            const std::array<double, 3> tempErrors{errors_[0], errors_[1], error};
            const std::array<double, 3> tempTimeSteps{timeSteps_[0], timeSteps_[1], timeStepJustCompleted};
            const double stepFactor = timeStepFactor(tempErrors, tempTimeSteps);
            if (rejectCompletedStep_ && stepFactor < maxReductionTimeStep_)
            {
                acceptTimeStep = false;
            }
        }

        if (acceptTimeStep)
        {
            // Shift errors and time steps
            for( int i = 0; i < 2; ++i )
            {
                errors_[i] = errors_[i+1];
                timeSteps_[i] = timeSteps_[i+1];
            }

            // Store new error and time step
            errors_[2] = error;
            timeSteps_[2] = timeStepJustCompleted;
        }
        
        return acceptTimeStep;
    }

    bool General3rdOrderController::operator==(const General3rdOrderController& ctrl) const
    {
        return this->tolerance_ == ctrl.tolerance_ &&
               this->safetyFactor_ == ctrl.safetyFactor_ &&
               this->rejectCompletedStep_ == ctrl.rejectCompletedStep_ &&
               this->errors_ == ctrl.errors_ &&
               this->timeSteps_ == ctrl.timeSteps_ &&
               this->beta_ == ctrl.beta_ &&
               this->alpha_ == ctrl.alpha_ &&
               this->controllerVersion_ == ctrl.controllerVersion_ &&
               this->toleranceTestVersion_ == ctrl.toleranceTestVersion_ &&
               this->maxReductionTimeStep_ == ctrl.maxReductionTimeStep_ &&
               this->verbose_ == ctrl.verbose_;
    }

} // end namespace Opm
