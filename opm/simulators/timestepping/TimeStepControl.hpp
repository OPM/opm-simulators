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
#ifndef OPM_TIMESTEPCONTROL_HEADER_INCLUDED
#define OPM_TIMESTEPCONTROL_HEADER_INCLUDED

#include <opm/simulators/timestepping/TimeStepControlInterface.hpp>
#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>

#include <string>
#include <vector>

namespace Opm
{
    enum class TimeStepControlType {
      SimpleIterationCount,
      PID,
      PIDAndIterationCount,
      HardCodedTimeStep,
      General3rdOrder,
    };

    enum class ToleranceTestVersions {
        Standard,
        ControlErrorFiltering,
    };

    enum class InternalControlVersions {
        IController,
        General3rdOrder,
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    ///  A simple iteration count based adaptive time step control.
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class SimpleIterationCountTimeStepControl : public TimeStepControlInterface
    {
    public:
        static constexpr TimeStepControlType Type = TimeStepControlType::SimpleIterationCount;

        SimpleIterationCountTimeStepControl() = default;

        /// \brief constructor
        /// \param target_iterations  number of desired iterations (e.g. Newton iterations) per time step in one time step
        /// \param decayrate          decayrate of time step when target iterations are not met (should be <= 1)
        /// \param growthrate         growthrate of time step when target iterations are not met (should be >= 1)
        /// \param verbose            if true, get some output
        SimpleIterationCountTimeStepControl(const int target_iterations,
                                            const double decayrate,
                                            const double growthrate,
                                            const bool verbose);

        static SimpleIterationCountTimeStepControl serializationTestObject();

        /// \brief \copydoc TimeStepControlInterface::computeTimeStepSize
        double computeTimeStepSize(const double dt,
                                   const int iterations,
                                   const RelativeChangeInterface& /* relativeChange */,
                                   const AdaptiveSimulatorTimer& /* substepTimer */ ) const override;

        bool timeStepAccepted(const double /* error */,
                              const double /* timeStepJustCompleted */) const override { return true; }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(target_iterations_);
            serializer(decayrate_);
            serializer(growthrate_);
            serializer(verbose_);
        }

        bool operator==(const SimpleIterationCountTimeStepControl&) const;

    protected:
        const int     target_iterations_ = 0;
        const double  decayrate_ = 0.0;
        const double  growthrate_ = 0.0;
        const bool    verbose_ = false;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    ///  PID controller based adaptive time step control as suggested in:
    ///     Turek and Kuzmin. Algebraic Flux Correction III. Incompressible Flow Problems. Uni Dortmund.
    ///
    ///  See also:
    ///     D. Kuzmin and S.Turek. Numerical simulation of turbulent bubbly flows. Techreport Uni Dortmund. 2004
    ///
    ///  and the original article:
    ///     Valli, Coutinho, and Carey. Adaptive Control for Time Step Selection in Finite Element
    ///     Simulation of Coupled Viscous Flow and Heat Transfer. Proc of the 10th
    ///     International Conference on Numerical Methods in Fluids. 1998.
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class PIDTimeStepControl : public TimeStepControlInterface
    {
    public:
        static constexpr TimeStepControlType Type = TimeStepControlType::PID;

        PIDTimeStepControl() = default;

        /// \brief constructor
        /// \param tol      tolerance for the relative changes of the numerical solution to be accepted
        ///                 in one time step
        /// \param verbose  if true, get some output
        PIDTimeStepControl(const double tol,
                           const bool verbose);

        static PIDTimeStepControl serializationTestObject();

        /// \brief \copydoc TimeStepControlInterface::computeTimeStepSize
        /// \param dt Time step length
        /// \param relativeChange Relative change handler
        double computeTimeStepSize(const double dt,
                                   const int /* iterations */,
                                   const RelativeChangeInterface& relativeChange,
                                   const AdaptiveSimulatorTimer& /* substepTimer */ ) const override;

        bool timeStepAccepted(const double /* error */,
                              const double /* timeStepJustCompleted */) const override { return true; }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(tol_);
            serializer(errors_);
            serializer(verbose_);
        }

        bool operator==(const PIDTimeStepControl&) const;

    protected:
        const double tol_ = 0.1;
        mutable std::vector< double > errors_{};
        const bool verbose_ = false;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    ///  PID controller based adaptive time step control as above that also takes
    ///  target iterations into account.
    //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class PIDAndIterationCountTimeStepControl : public PIDTimeStepControl
    {
        typedef PIDTimeStepControl BaseType;
    public:
        static constexpr TimeStepControlType Type = TimeStepControlType::PIDAndIterationCount;

        PIDAndIterationCountTimeStepControl() = default;

        /// \brief constructor
        /// \param target_iterations            number of desired iterations per time step
        /// \param decayDampingFactor           limiting the decrease in time step if iteration count was high
        /// \param growthDampingFactor          limiting the increase in time step if iteration count was low
        /// \param tol                          tolerance for the relative changes of the numerical solution to be
        ///                                     accepted in one time step
        /// \param minTimeStepBasedOnIterations time step suggestion from target iterations should not be below this
        /// \param verbose                      if true, get some output
        PIDAndIterationCountTimeStepControl(const int target_iterations,
                                            const double decayDampingFactor,
                                            const double growthDampingFactor,
                                            const double tol,
                                            const double minTimeStepBasedOnIterations,
                                            const bool verbose);

        static PIDAndIterationCountTimeStepControl serializationTestObject();

        /// \brief \copydoc TimeStepControlInterface::computeTimeStepSize
        /// \param dt Time step length
        /// \param iterations Number of iterations used
        /// \param relativeChange Relative change handler
        double computeTimeStepSize(const double dt,
                                   const int iterations,
                                   const RelativeChangeInterface& relativeChange,
                                   const AdaptiveSimulatorTimer& /* substepTimer */ ) const override;

        bool timeStepAccepted(const double /* error */,
                              const double /* timeStepJustCompleted */) const override { return true; }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(static_cast<PIDTimeStepControl&>(*this));
            serializer(target_iterations_);
            serializer(decayDampingFactor_);
            serializer(growthDampingFactor_);
            serializer(minTimeStepBasedOnIterations_);
            serializer(verbose_);
        }

        bool operator==(const PIDAndIterationCountTimeStepControl&) const;

    protected:
        const int target_iterations_ = 8;
        const double decayDampingFactor_ = 1.0;
        const double growthDampingFactor_ = 3.2;
        const double minTimeStepBasedOnIterations_ = 0.0;
        const bool verbose_ = false;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    ///  General 3rd order controller
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class General3rdOrderController : public TimeStepControlInterface
    {
    public:
        static constexpr TimeStepControlType Type = TimeStepControlType::General3rdOrder;

        General3rdOrderController() = default;

        /// \brief constructor
        /// \param tolerance                    tolerance for the relative changes of the numerical solution to be
        ///                                     accepted in one time step
        /// \param safetyFactor                 multiplied with tolerance to ensure that target relative change is
        ///                                     lower than tolerance
        /// \param rejectCompletedStep          if true, discard the recently completed time step and try again
        /// \param toleranceTestVersion         the test used to decide if a time step should be rejected
        /// \param maxReductionTimeStep         limits the reduction in time step size for control error filtering
        /// \param parameters                   parameter values for the controller formula
        /// \param verbose                      if true, get some output
        General3rdOrderController(const double tolerance,
                                  const double safetyFactor,
                                  const bool rejectCompletedStep,
                                  const std::string& toleranceTestVersion,
                                  const double maxReductionTimeStep,
                                  const std::string& parameters,
                                  const bool verbose);

        static General3rdOrderController serializationTestObject();

        double computeTimeStepSize(const double dt,
                                   const int /* iterations */,
                                   const RelativeChangeInterface& /* relativeChange */,
                                   const AdaptiveSimulatorTimer& substepTimer) const override;

        double timeStepFactor(const std::array<double, 3>& errors, const std::array<double, 3>& timeSteps) const;

        bool timeStepAccepted(const double error,
                              const double timeStepJustCompleted) const override;

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(tolerance_);
            serializer(safetyFactor_);
            serializer(rejectCompletedStep_);
            serializer(errors_);
            serializer(timeSteps_);
            serializer(beta_);
            serializer(alpha_);
            serializer(controllerVersion_);
            serializer(toleranceTestVersion_);
            serializer(maxReductionTimeStep_);
            serializer(verbose_);
        }

        bool operator==(const General3rdOrderController&) const;


    protected:
        const double tolerance_ = 0.1;
        const double safetyFactor_ = 0.8;
        const bool rejectCompletedStep_ = false;
        mutable std::array<double, 3> errors_{};
        mutable std::array<double, 3> timeSteps_{};
        mutable std::array<double, 3> beta_{0.125, 0.25, 0.125};
        mutable std::array<double, 2> alpha_{0.75, 0.25};
        mutable InternalControlVersions controllerVersion_{InternalControlVersions::IController};
        ToleranceTestVersions toleranceTestVersion_{ToleranceTestVersions::Standard};
        const double maxReductionTimeStep_ = 0.1;
        const bool verbose_ = false;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    ///  HardcodedTimeStepControl
    ///  Input generated from summary file using the ert application:
    ///
    ///  ecl_summary DECK TIME > filename
    ///
    ///  Assumes time is given in days
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class HardcodedTimeStepControl : public TimeStepControlInterface
    {
    public:
        static constexpr TimeStepControlType Type = TimeStepControlType::HardCodedTimeStep;

        HardcodedTimeStepControl() = default;

        /// \brief constructor
        /// \param filename   filename contaning the timesteps
        explicit HardcodedTimeStepControl(const std::string& filename);

        static HardcodedTimeStepControl serializationTestObject();

        /// \brief \copydoc TimeStepControlInterface::computeTimeStepSize
        /// \param dt Time step length
        /// \param substepTimer Sub step timer
        double computeTimeStepSize(const double dt,
                                   const int /* iterations */,
                                   const RelativeChangeInterface& /*relativeChange */,
                                   const AdaptiveSimulatorTimer& substepTimer) const override;

        bool timeStepAccepted(const double /* error */,
                              const double /* timeStepJustCompleted */) const override { return true; }

        template<class Serializer>
        void serializeOp(Serializer& serializer)
        {
            serializer(subStepTime_);
        }

        bool operator==(const HardcodedTimeStepControl&) const;

    protected:
        // store the time (in days) of the substeps the simulator should use
        std::vector<double> subStepTime_;
    };


} // end namespace Opm
#endif
