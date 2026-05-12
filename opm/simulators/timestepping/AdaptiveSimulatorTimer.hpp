/*
  Copyright 2014 IRIS AS

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
#ifndef OPM_ADAPTIVESIMULATORTIMER_HEADER_INCLUDED
#define OPM_ADAPTIVESIMULATORTIMER_HEADER_INCLUDED

#include <cassert>
#include <iosfwd>
#include <vector>
#include <limits>
#include <algorithm>
#include <memory>
#include <numeric>
#include <optional>

#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>

namespace Opm
{

    /////////////////////////////////////////////////////////
    ///
    /// \brief Simulation timer for adaptive time stepping
    ///
    /////////////////////////////////////////////////////////
    class AdaptiveSimulatorTimer : public SimulatorTimerInterface
    {
    public:
        /// \brief constructor taking a simulator timer to determine start and end time
        ///  \param simulation_start_time Start time for simulation
        ///  \param step_length Time step length
        ///  \param elapsed_time Current elapsed time
        ///  \param last_step_taken Last suggested time step
        ///  \param report_step Current report step index
        ///  \param max_time_step    maximum time step allowed
        AdaptiveSimulatorTimer(const boost::posix_time::ptime simulation_start_time,
                               const double step_length,
                               const double elapsed_time,
                               const double last_step_taken,
                               const int report_step,
                               const double max_time_step = std::numeric_limits<double>::max());

        /// \brief advance time by currentStepLength
        AdaptiveSimulatorTimer& operator++ ();

        /// \brief advance time by currentStepLength
        void advance() override { this->operator++ (); }

        /// \brief provide and estimate for new time step size
        void provideTimeStepEstimate( const double dt_estimate );

        /// \brief Whether this is the first step
        bool initialStep () const override;

        /// \brief \copydoc SimulationTimer::currentStepNum
        int currentStepNum () const override;

        /// \brief return current report step
        int reportStepNum() const override;

        /// \brief \copydoc SimulationTimer::currentStepLength
        double currentStepLength () const override;

        // \brief Set next step length
        void setCurrentStepLength(double dt);

        /// \brief \copydoc SimulationTimer::totalTime
        double totalTime() const;

        /// \brief \copydoc SimulationTimer::simulationTimeElapsed
        double simulationTimeElapsed() const override;

        /// \brief \copydoc SimulationTimer::done
        bool done () const override;

        /// \brief return average step length used so far
        double averageStepLength() const;

        /// \brief return max step length used so far
        double maxStepLength () const;

        /// \brief return min step length used so far
        double minStepLength () const;

        /// \brief Previous step length. This is the length of the step that
        ///        was taken to arrive at this time.
        double stepLengthTaken () const override;

        /// \brief report start and end time as well as used steps so far
        void report(std::ostream& os) const;

        /// \brief start date time of simulation
        boost::posix_time::ptime startDateTime() const override;

        /// \brief Return true if last time step failed
        bool lastStepFailed() const override { return last_step_failed_; }

        /// \brief tell the timestepper whether timestep failed or not
        void setLastStepFailed(bool last_step_failed) { last_step_failed_ = last_step_failed; }

        /// \brief Reservoir coupling constructs a fresh timer per sync chunk,
        /// so `start_time_`, `total_time_`, and `current_step_` describe the
        /// chunk rather than the enclosing report step. The accessors below
        /// return a "report step view" that the rescoup outer loop populates
        /// (via the setters) so that log lines can show the report-step start
        /// and end and a substep counter that increments across chunks. When
        /// the view is not set, they fall through to the per-timer fields, so
        /// the non-rescoup path's behavior and output are unchanged.
        double reportStepStartTime() const;
        double reportStepTotalTime() const;
        int reportStepSubstepNum() const;

        void setReportStepStartTime(double t) { report_step_start_time_ = t; }
        void setReportStepTotalTime(double t) { report_step_total_time_ = t; }
        void setReportStepSubstepOffset(int n) { report_step_substep_offset_ = n; }

        int reportStepSubstepOffset() const { return report_step_substep_offset_; }

        /// return copy of object
        std::unique_ptr<SimulatorTimerInterface> clone() const override;

    protected:
        std::shared_ptr<boost::posix_time::ptime> start_date_time_;
        const double start_time_;
        const double total_time_;
        const int report_step_;
        const double max_time_step_;

        double current_time_;
        double dt_;
        int current_step_;

        std::vector< double > steps_;
        bool last_step_failed_;

        /// \brief Optional report-step start time for the "report step view"
        ///   accessors. Set by the rescoup outer loop on each per-sync-chunk
        ///   timer; unset on a non-rescoup timer.
        std::optional<double> report_step_start_time_;
        /// \brief Optional report-step end time for the "report step view"
        ///   accessors. Same population pattern as `report_step_start_time_`.
        std::optional<double> report_step_total_time_;
        /// \brief Number of substeps already taken in this report step before
        ///   this timer was constructed (i.e. in earlier sync chunks). Added to
        ///   `current_step_` by `reportStepSubstepNum()` to give a counter that
        ///   increments across sync chunks. Zero on a non-rescoup timer.
        int report_step_substep_offset_ = 0;

    };

} // namespace Opm

#endif // OPM_SIMULATORTIMER_HEADER_INCLUDED
