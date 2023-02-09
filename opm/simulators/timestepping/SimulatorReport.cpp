/*
  Copyright 2012, 2020 SINTEF Digital, Mathematics and Cybernetics.

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

#include "config.h"

#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <iomanip>
#include <ostream>
#include <sstream>
#include <fmt/format.h>

namespace Opm
{
    void SimulatorReportSingle::operator+=(const SimulatorReportSingle& sr)
    {
        pressure_time += sr.pressure_time;
        transport_time += sr.transport_time;
        linear_solve_setup_time += sr.linear_solve_setup_time;
        linear_solve_time += sr.linear_solve_time;
        solver_time += sr.solver_time;
        assemble_time += sr.assemble_time;
        pre_post_time += sr.pre_post_time;
        assemble_time_well += sr.assemble_time_well;
        update_time += sr.update_time;
        output_write_time += sr.output_write_time;
        total_time += sr.total_time;
        total_well_iterations += sr.total_well_iterations;
        total_linearizations += sr.total_linearizations;
        total_newton_iterations += sr.total_newton_iterations;
        total_linear_iterations += sr.total_linear_iterations;
        if (sr.total_linear_iterations > 0) {
            min_linear_iterations = std::min(min_linear_iterations, sr.total_linear_iterations);
        }
        max_linear_iterations = std::max(max_linear_iterations, sr.total_linear_iterations);

        // It makes no sense adding time points. Therefore, do not 
        // overwrite the value of global_time which gets set in 
        // NonlinearSolverEbos.hpp by the line:
        //     report.global_time = timer.simulationTimeElapsed();
    }


    void SimulatorReportSingle::reportStep(std::ostringstream& ss) const
    {
        if (total_well_iterations != 0) {
            ss << fmt::format("Well its={:2}", total_well_iterations);
        }
        ss << fmt::format(" Newton its={:2}, linearizations={:2} ({:2.1f}sec), linear its={:3} ({:2.1f}sec)",
                          total_newton_iterations,
                          total_linearizations,
                          assemble_time,
                          total_linear_iterations,
                          linear_solve_time);
    }

    void SimulatorReportSingle::reportFullyImplicit(std::ostream& os, const SimulatorReportSingle* failureReport) const
    {
        os << fmt::format("Total time (seconds):       {:9.2f} \n", total_time);

         os << fmt::format("Solver time (seconds):      {:9.2f} \n",
                          solver_time + (failureReport ? failureReport->solver_time : 0.0));

        if (assemble_time > 0.0 || linear_solve_time > 0.0) {

            double t = assemble_time + (failureReport ? failureReport->assemble_time : 0.0);
            os << fmt::format(" Assembly time (seconds):   {:9.2f}", t);

            if (failureReport) {
              os << fmt::format(" (Failed: {:2.1f}; {:2.1f}%)",
                                failureReport->assemble_time,
                                100*failureReport->assemble_time/t);
             }
            os << std::endl;

            t = assemble_time_well + (failureReport ? failureReport->assemble_time_well : 0.0);
            os << fmt::format("   Well assembly (seconds):   {:7.2f}", t);
            if (failureReport) {
              os << fmt::format(" (Failed: {:2.1f}; {:2.1f}%)",
                                failureReport->assemble_time_well,
                                100*failureReport->assemble_time_well/t);
            }
            os << std::endl;

            t = linear_solve_time + (failureReport ? failureReport->linear_solve_time : 0.0);
            os << fmt::format(" Linear solve time (seconds):{:8.2f}", t);
            if (failureReport) {
              os << fmt::format(" (Failed: {:2.1f}; {:2.1f}%)",
                                failureReport->linear_solve_time,
                                100*failureReport->linear_solve_time/t);
            }
            os << std::endl;

            t = linear_solve_setup_time + (failureReport ? failureReport->linear_solve_setup_time : 0.0);
            os << fmt::format("   Linear setup (seconds):    {:7.2f}", t);
            if (failureReport) {
              os << fmt::format(" (Failed: {:2.1f}; {:2.1f}%)",
                                failureReport->linear_solve_setup_time,
                                100*failureReport->linear_solve_setup_time/t);
            }
            os << std::endl;

            t = update_time + (failureReport ? failureReport->update_time : 0.0);
            os << fmt::format(" Update time (seconds):       {:7.2f}", t);
            if (failureReport) {
              os << fmt::format(" (Failed: {:2.1f}; {:2.1f}%)",
                                failureReport->update_time,
                                100*failureReport->update_time/t);
            }
            os << std::endl;
            t = pre_post_time + (failureReport ? failureReport->pre_post_time : 0.0);
            os << fmt::format(" Pre/post step (seconds):     {:7.2f}", t);
            if (failureReport) {
              os << fmt::format(" (Failed: {:2.1f}; {:2.1f}%)",
                                failureReport->pre_post_time,
                                100*failureReport->pre_post_time/t);
            }
            os << std::endl;

            os << fmt::format(" Output write time (seconds): {:7.2f}", 
                              output_write_time + (failureReport ? failureReport->output_write_time : 0.0));
            os << std::endl;

        }

        int n = total_linearizations + (failureReport ? failureReport->total_linearizations : 0);
        os << fmt::format("Overall Linearizations:    {:7}", n);
        if (failureReport) {
          os << fmt::format("    (Failed: {:3}; {:2.1f}%)",
                            failureReport->total_linearizations,
                            100.0*failureReport->total_linearizations/n);
        }
        os << std::endl;

        n = total_newton_iterations + (failureReport ? failureReport->total_newton_iterations : 0);
        os << fmt::format("Overall Newton Iterations: {:7}", n);
        if (failureReport) {
          os << fmt::format("    (Failed: {:3}; {:2.1f}%)",
                            failureReport->total_newton_iterations,
                            100.0*failureReport->total_newton_iterations/n);
        }
        os << std::endl;

        n = total_linear_iterations + (failureReport ? failureReport->total_linear_iterations : 0);
        os << fmt::format("Overall Linear Iterations: {:7}", n);
        if (failureReport) {
          os << fmt::format("    (Failed: {:3}; {:2.1f}%)",
                            failureReport->total_linear_iterations,
                            100.0*failureReport->total_linear_iterations/n);
        }
        os << std::endl;
    }

    void SimulatorReport::operator+=(const SimulatorReportSingle& sr)
    {
        if (sr.converged) {
            success += sr;
        } else {
            failure += sr;
        }
        stepreports.push_back(sr);
    }

    void SimulatorReport::operator+=(const SimulatorReport& sr)
    {
        success += sr.success;
        failure += sr.failure;
        stepreports.insert(stepreports.end(), sr.stepreports.begin(), sr.stepreports.end());
    }

    void SimulatorReport::reportFullyImplicit(std::ostream& os) const
    {
        success.reportFullyImplicit(os, &failure);
    }

    void SimulatorReport::fullReports(std::ostream& os) const
    {
        os << "  Time(day)  TStep(day)  Assembly    LSetup    LSolve    Update    Output WellIt Lins NewtIt LinIt Conv\n";
        for (size_t i = 0; i < this->stepreports.size(); ++i) {
            const SimulatorReportSingle& sr = this->stepreports[i];
            os.precision(10);
            os << std::defaultfloat;
            os << std::setw(11) << unit::convert::to(sr.global_time, unit::day) << " ";
            os << std::setw(11) << unit::convert::to(sr.timestep_length, unit::day) << " ";
            os.precision(4);
            os << std::fixed;
            os << std::setw(9) << sr.assemble_time << " ";
            os << std::setw(9) << sr.linear_solve_setup_time << " ";
            os << std::setw(9) << sr.linear_solve_time << " ";
            os << std::setw(9) << sr.update_time << " ";
            os << std::setw(9) << sr.output_write_time << " ";
            os.precision(6);
            os << std::defaultfloat;
            os << std::setw(6) << sr.total_well_iterations << " ";
            os << std::setw(4) << sr.total_linearizations << " ";
            os << std::setw(6) << sr.total_newton_iterations << " ";
            os << std::setw(5) << sr.total_linear_iterations << " ";
            os << std::setw(4) << sr.converged << "\n";
        }
    }

} // namespace Opm
