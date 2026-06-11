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

#include <config.h>
#include <opm/simulators/timestepping/SimulatorReport.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <cstddef>
#include <iomanip>
#include <ostream>
#include <fmt/format.h>

namespace Opm
{
    SimulatorReportSingle SimulatorReportSingle::serializationTestObject()
    {
        return SimulatorReportSingle{1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                                     7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
                                     13, 14, 15, 16, 17, 18,
                                     true, false, false, 19, 20.0, 21.0,
                                     22, 23, 24, 25, 26, 27, 28, 29,
                                     30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0,
                                     38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46, 47};
    }

    bool SimulatorReportSingle::operator==(const SimulatorReportSingle& rhs) const
    {
        return this->pressure_time == rhs.pressure_time &&
               this->transport_time == rhs.transport_time &&
               this->total_time == rhs.total_time &&
               this->solver_time == rhs.solver_time &&
               this->assemble_time == rhs.assemble_time &&
               this->pre_post_time == rhs.pre_post_time &&
               this->assemble_time_well == rhs.assemble_time_well &&
               this->linear_solve_setup_time == rhs.linear_solve_setup_time &&
               this->linear_solve_time == rhs.linear_solve_time &&
               this->local_solve_time == rhs.local_solve_time &&
               this->update_time == rhs.update_time &&
               this->output_write_time == rhs.output_write_time &&
               this->total_well_iterations == rhs.total_well_iterations &&
               this->total_linearizations == rhs.total_linearizations &&
               this->total_newton_iterations == rhs.total_newton_iterations &&
               this->total_linear_iterations == rhs.total_linear_iterations &&
               this->min_linear_iterations == rhs.min_linear_iterations &&
               this->max_linear_iterations == rhs.max_linear_iterations &&
               this->converged == rhs.converged &&
               this->time_step_rejected == rhs.time_step_rejected &&
               this->well_group_control_changed == rhs.well_group_control_changed &&
               this->exit_status == rhs.exit_status &&
               this->global_time == rhs.global_time &&
               this->timestep_length == rhs.timestep_length &&
               this->num_domains == rhs.num_domains &&
               this->num_wells == rhs.num_wells &&
               this->num_overlap_cells == rhs.num_overlap_cells &&
               this->num_owned_cells == rhs.num_owned_cells &&
               this->converged_domains == rhs.converged_domains &&
               this->unconverged_domains == rhs.unconverged_domains &&
               this->accepted_unconverged_domains == rhs.accepted_unconverged_domains &&
               this->skipped_domains == rhs.skipped_domains &&
               this->props_time == rhs.props_time &&
               this->convergence_check_time == rhs.convergence_check_time &&
               this->output_eval_time == rhs.output_eval_time &&
               this->tracer_solve_time == rhs.tracer_solve_time &&
               this->temperature_solve_time == rhs.temperature_solve_time &&
               this->output_disk_write_time == rhs.output_disk_write_time &&
               this->linear_solve_apply_time == rhs.linear_solve_apply_time &&
               this->precond_setup_time == rhs.precond_setup_time &&
               this->precond_apply_time == rhs.precond_apply_time &&
               this->well_solve_time == rhs.well_solve_time &&
               this->well_potential_solve_time == rhs.well_potential_solve_time &&
               this->well_solve_assemble_time == rhs.well_solve_assemble_time &&
               this->well_solve_linear_solve_time == rhs.well_solve_linear_solve_time &&
               this->well_control_network_time == rhs.well_control_network_time &&
               this->gaslift_time == rhs.gaslift_time &&
               this->well_facility_time == rhs.well_facility_time &&
               this->total_well_potential_iterations == rhs.total_well_potential_iterations &&
               this->total_network_iterations == rhs.total_network_iterations;
    }

    void SimulatorReportSingle::operator+=(const SimulatorReportSingle& sr)
    {
        pressure_time += sr.pressure_time;
        transport_time += sr.transport_time;
        linear_solve_setup_time += sr.linear_solve_setup_time;
        linear_solve_time += sr.linear_solve_time;
        local_solve_time += sr.local_solve_time;
        solver_time += sr.solver_time;
        assemble_time += sr.assemble_time;
        pre_post_time += sr.pre_post_time;
        assemble_time_well += sr.assemble_time_well;
        update_time += sr.update_time;
        props_time += sr.props_time;
        convergence_check_time += sr.convergence_check_time;
        output_eval_time += sr.output_eval_time;
        tracer_solve_time += sr.tracer_solve_time;
        temperature_solve_time += sr.temperature_solve_time;
        output_disk_write_time += sr.output_disk_write_time;
        linear_solve_apply_time += sr.linear_solve_apply_time;
        precond_setup_time += sr.precond_setup_time;
        precond_apply_time += sr.precond_apply_time;
        well_solve_time += sr.well_solve_time;
        well_potential_solve_time += sr.well_potential_solve_time;
        well_solve_assemble_time += sr.well_solve_assemble_time;
        well_solve_linear_solve_time += sr.well_solve_linear_solve_time;
        well_control_network_time += sr.well_control_network_time;
        gaslift_time += sr.gaslift_time;
        well_facility_time += sr.well_facility_time;
        output_write_time += sr.output_write_time;
        total_time += sr.total_time;
        total_well_iterations += sr.total_well_iterations;
        total_well_potential_iterations += sr.total_well_potential_iterations;
        total_network_iterations += sr.total_network_iterations;
        total_linearizations += sr.total_linearizations;
        total_newton_iterations += sr.total_newton_iterations;
        total_linear_iterations += sr.total_linear_iterations;
        if (sr.total_linear_iterations > 0) {
            min_linear_iterations = std::min(min_linear_iterations, sr.total_linear_iterations);
        }
        max_linear_iterations = std::max(max_linear_iterations, sr.total_linear_iterations);

        converged_domains += sr.converged_domains;
        unconverged_domains += sr.unconverged_domains;
        accepted_unconverged_domains += sr.accepted_unconverged_domains;
        skipped_domains += sr.skipped_domains;
        // It makes no sense adding time points. Therefore, do not
        // overwrite the value of global_time which gets set in
        // NonlinearSolver.hpp by the line:
        //     report.global_time = timer.simulationTimeElapsed();
    }


    void SimulatorReportSingle::reportStep(std::ostream& ss) const
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

    // Helper lambda to avoid division by zero
    auto noZero = [](auto val)
    {
        if (val == decltype(val){0})
            return decltype(val){1};
        return val;
    };

    void SimulatorReportSingle::reportFullyImplicit(std::ostream& os,
                                                    const SimulatorReportSingle* failureReport,
                                                    bool performance_details) const
    {
        // Disabling this output as it is redundant with the solver_time, now
        // output as "Simulation time". Usually very small difference between them.
        // os << fmt::format("Total time:                 {:9.2f} s\n", total_time);

        os << fmt::format("Simulation time:            {:9.2f} s\n",
                          solver_time + (failureReport ? failureReport->solver_time : 0.0));

        if (assemble_time > 0.0 || linear_solve_time > 0.0) {

            double t = assemble_time + (failureReport ? failureReport->assemble_time : 0.0);
            os << fmt::format("  Assembly time:            {:9.2f} s", t);

            if (failureReport) {
              os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                failureReport->assemble_time,
                                100*failureReport->assemble_time/noZero(t));
             }
            os << std::endl;

            t = assemble_time_well + (failureReport ? failureReport->assemble_time_well : 0.0);
            os << fmt::format("    Well assembly:            {:7.2f} s", t);
            if (failureReport) {
              os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                failureReport->assemble_time_well,
                                100*failureReport->assemble_time_well/noZero(t));
            }
            os << std::endl;

            // The facility calculations cover the time step preparation of
            // the wells and the control/network updates. The indented items
            // below are parts of their parent, but do not necessarily add
            // up to it (e.g. well solves also happen during the
            // control/network updates).
            t = well_facility_time + (failureReport ? failureReport->well_facility_time : 0.0);
            if (performance_details && t > 0.0) {
                os << fmt::format("      Facility calculations:  {:7.2f} s", t);
                os << std::endl;
            }

            t = well_control_network_time + (failureReport ? failureReport->well_control_network_time : 0.0);
            if (performance_details && t > 0.0) {
                os << fmt::format("        Control/network:      {:7.2f} s", t);
                os << std::endl;
            }

            t = gaslift_time + (failureReport ? failureReport->gaslift_time : 0.0);
            if (performance_details && t > 0.0) {
                os << fmt::format("          Gas lift optimize:  {:7.2f} s", t);
                os << std::endl;
            }

            t = well_solve_time + (failureReport ? failureReport->well_solve_time : 0.0);
            if (performance_details && t > 0.0) {
                os << fmt::format("        Well solves:          {:7.2f} s", t);
                os << std::endl;
            }

            t = linear_solve_time + (failureReport ? failureReport->linear_solve_time : 0.0);
            os << fmt::format("  Linear solve time:        {:9.2f} s", t);
            if (failureReport) {
              os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                failureReport->linear_solve_time,
                                100*failureReport->linear_solve_time/noZero(t));
            }
            os << std::endl;

            t = linear_solve_setup_time + (failureReport ? failureReport->linear_solve_setup_time : 0.0);
            os << fmt::format("    Linear setup:             {:7.2f} s", t);
            if (failureReport) {
              os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                failureReport->linear_solve_setup_time,
                                100*failureReport->linear_solve_setup_time/noZero(t));
            }
            os << std::endl;

            t = precond_setup_time + (failureReport ? failureReport->precond_setup_time : 0.0);
            if (performance_details && t > 0.0) {
                os << fmt::format("      Precond setup:          {:7.2f} s", t);
                os << std::endl;
            }

            if (performance_details) {
                t = linear_solve_apply_time + (failureReport ? failureReport->linear_solve_apply_time : 0.0);
                os << fmt::format("    Linear apply:             {:7.2f} s", t);
                if (failureReport) {
                  os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                    failureReport->linear_solve_apply_time,
                                    100*failureReport->linear_solve_apply_time/noZero(t));
                }
                os << std::endl;
            }

            t = precond_apply_time + (failureReport ? failureReport->precond_apply_time : 0.0);
            if (performance_details && t > 0.0) {
                os << fmt::format("      Precond apply:          {:7.2f} s", t);
                os << std::endl;
            }

            if (local_solve_time > 0.0) {
                t = local_solve_time + (failureReport ? failureReport->local_solve_time : 0.0);
                os << fmt::format("  Local solve time:           {:7.2f} s", t);
                if (failureReport) {
                  os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                    failureReport->local_solve_time,
                                    100*failureReport->local_solve_time/noZero(t));
                }
                os << std::endl;
            }

            t = update_time + (failureReport ? failureReport->update_time : 0.0);
            os << fmt::format("  Props/update time:          {:7.2f} s", t);
            if (failureReport) {
              os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                failureReport->update_time,
                                100*failureReport->update_time/noZero(t));
            }
            os << std::endl;

            if (performance_details) {
                t = props_time + (failureReport ? failureReport->props_time : 0.0);
                os << fmt::format("    Properties evaluation:    {:7.2f} s", t);
                if (failureReport) {
                  os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                    failureReport->props_time,
                                    100*failureReport->props_time/noZero(t));
                }
                os << std::endl;
            }

            if (performance_details) {
                t = convergence_check_time + (failureReport ? failureReport->convergence_check_time : 0.0);
                os << fmt::format("    Convergence checks:       {:7.2f} s", t);
                if (failureReport) {
                  os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                    failureReport->convergence_check_time,
                                    100*failureReport->convergence_check_time/noZero(t));
                }
                os << std::endl;
            }

            t = pre_post_time + (failureReport ? failureReport->pre_post_time : 0.0);
            os << fmt::format("  Pre/post step:              {:7.2f} s", t);
            if (failureReport) {
              os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                failureReport->pre_post_time,
                                100*failureReport->pre_post_time/noZero(t));
            }
            os << std::endl;

            if (performance_details) {
                t = output_eval_time + (failureReport ? failureReport->output_eval_time : 0.0);
                os << fmt::format("    Output evaluation:        {:7.2f} s", t);
                if (failureReport) {
                  os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                    failureReport->output_eval_time,
                                    100*failureReport->output_eval_time/noZero(t));
                }
                os << std::endl;
            }

            t = tracer_solve_time + (failureReport ? failureReport->tracer_solve_time : 0.0);
            if (performance_details && t > 0.0) {
                os << fmt::format("    Tracer solve:             {:7.2f} s", t);
                if (failureReport) {
                  os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                    failureReport->tracer_solve_time,
                                    100*failureReport->tracer_solve_time/noZero(t));
                }
                os << std::endl;
            }

            t = temperature_solve_time + (failureReport ? failureReport->temperature_solve_time : 0.0);
            if (performance_details && t > 0.0) {
                os << fmt::format("    Temperature solve:        {:7.2f} s", t);
                if (failureReport) {
                  os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                    failureReport->temperature_solve_time,
                                    100*failureReport->temperature_solve_time/noZero(t));
                }
                os << std::endl;
            }

            os << fmt::format("  Output write time:          {:7.2f} s",
                              output_write_time + (failureReport ? failureReport->output_write_time : 0.0));
            os << std::endl;

            t = output_disk_write_time + (failureReport ? failureReport->output_disk_write_time : 0.0);
            if (performance_details && t > 0.0) {
                os << fmt::format("    Actual disk write:        {:7.2f} s", t);
                os << std::endl;
            }

            const double total_well_solve_time =
                well_solve_time + well_potential_solve_time +
                (failureReport ? failureReport->well_solve_time +
                                 failureReport->well_potential_solve_time : 0.0);
            if (performance_details && total_well_solve_time > 0.0) {
                os << fmt::format("  Well solve time:            {:7.2f} s", total_well_solve_time);
                if (failureReport) {
                  os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                                    failureReport->well_solve_time + failureReport->well_potential_solve_time,
                                    100*(failureReport->well_solve_time + failureReport->well_potential_solve_time)
                                       /noZero(total_well_solve_time));
                }
                os << std::endl;

                t = well_solve_assemble_time + (failureReport ? failureReport->well_solve_assemble_time : 0.0);
                os << fmt::format("    Assembly:                 {:7.2f} s", t);
                os << std::endl;

                t = well_solve_linear_solve_time + (failureReport ? failureReport->well_solve_linear_solve_time : 0.0);
                os << fmt::format("    Linear solve:             {:7.2f} s", t);
                os << std::endl;

                t = well_potential_solve_time + (failureReport ? failureReport->well_potential_solve_time : 0.0);
                os << fmt::format("    Potential solves:         {:7.2f} s", t);
                os << std::endl;
            }
        }

        int n = total_linearizations + (failureReport ? failureReport->total_linearizations : 0);
        os << fmt::format("Overall Linearizations:    {:7}", n);
        if (failureReport) {
          os << fmt::format("      (Wasted: {:5}; {:2.1f}%)",
                            failureReport->total_linearizations,
                            100.0*failureReport->total_linearizations/noZero(n));
        }
        os << std::endl;

        n = total_newton_iterations + (failureReport ? failureReport->total_newton_iterations : 0);
        os << fmt::format("Overall Newton Iterations: {:7}", n);
        if (failureReport) {
          os << fmt::format("      (Wasted: {:5}; {:2.1f}%)",
                            failureReport->total_newton_iterations,
                            100.0*failureReport->total_newton_iterations/noZero(n));
        }
        os << std::endl;

        n = total_linear_iterations + (failureReport ? failureReport->total_linear_iterations : 0);
        os << fmt::format("Overall Linear Iterations: {:7}", n);
        if (failureReport) {
          os << fmt::format("      (Wasted: {:5}; {:2.1f}%)",
                            failureReport->total_linear_iterations,
                            100.0*failureReport->total_linear_iterations/noZero(n));
        }
        os << std::endl;

        n = total_well_iterations + (failureReport ? failureReport->total_well_iterations : 0);
        unsigned int np = total_well_potential_iterations +
            (failureReport ? failureReport->total_well_potential_iterations : 0);
        if (performance_details && (n > 0 || np > 0)) {
            os << fmt::format("Overall Well Iterations:   {:7}", n);
            if (failureReport) {
              os << fmt::format("      (Wasted: {:5}; {:2.1f}%)",
                                failureReport->total_well_iterations,
                                100.0*failureReport->total_well_iterations/noZero(n));
            }
            os << std::endl;

            os << fmt::format("Well Potential Iterations: {:7}", np);
            if (failureReport) {
              os << fmt::format("      (Wasted: {:5}; {:2.1f}%)",
                                failureReport->total_well_potential_iterations,
                                100.0*failureReport->total_well_potential_iterations/noZero(np));
            }
            os << std::endl;
        }

        unsigned int nn = total_network_iterations +
            (failureReport ? failureReport->total_network_iterations : 0);
        if (performance_details && nn > 0) {
            os << fmt::format("Network Balance Iterations:{:7}", nn);
            if (failureReport) {
              os << fmt::format("      (Wasted: {:5}; {:2.1f}%)",
                                failureReport->total_network_iterations,
                                100.0*failureReport->total_network_iterations/noZero(nn));
            }
            os << std::endl;
        }
    }


    void SimulatorReportSingle::reportNLDD(std::ostream& os, const SimulatorReportSingle* failureReport) const
    {
        os << fmt::format("Owned + overlap cells:       {:7}\n", num_owned_cells + num_overlap_cells);
        os << fmt::format("Number of wells:             {:7}\n", num_wells);
        os << fmt::format("Number of domains:           {:7}\n", num_domains);
        os << fmt::format("-------------------------------------------------------\n");
        double t = total_time + (failureReport ? failureReport->total_time : 0.0);
        os << fmt::format("Total time:                   {:9.2f} s\n", t);

        t = pre_post_time + (failureReport ? failureReport->pre_post_time : 0.0);
        os << fmt::format("  Pre/post/wait time:           {:7.2f} s\n", t);

        t = solver_time + (failureReport ? failureReport->solver_time : 0.0);
        os << fmt::format("  Solver time:                {:9.2f} s", t);
        if (failureReport) {
            os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                            failureReport->solver_time,
                            100*failureReport->solver_time/noZero(t));
        }
        os << std::endl;

        t = assemble_time + (failureReport ? failureReport->assemble_time : 0.0);
        os << fmt::format("    Assembly time:            {:9.2f} s", t);
        if (failureReport) {
            os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                            failureReport->assemble_time,
                            100*failureReport->assemble_time/noZero(t));
            }
        os << std::endl;

        t = assemble_time_well + (failureReport ? failureReport->assemble_time_well : 0.0);
        os << fmt::format("      Well assembly:            {:7.2f} s", t);
        if (failureReport) {
            os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                            failureReport->assemble_time_well,
                            100*failureReport->assemble_time_well/noZero(t));
        }
        os << std::endl;

        t = linear_solve_time + (failureReport ? failureReport->linear_solve_time : 0.0);
        os << fmt::format("    Linear solve time:        {:9.2f} s", t);
        if (failureReport) {
            os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                            failureReport->linear_solve_time,
                            100*failureReport->linear_solve_time/noZero(t));
        }
        os << std::endl;

        t = linear_solve_setup_time + (failureReport ? failureReport->linear_solve_setup_time : 0.0);
        os << fmt::format("      Linear setup:             {:7.2f} s", t);
        if (failureReport) {
            os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                            failureReport->linear_solve_setup_time,
                            100*failureReport->linear_solve_setup_time/noZero(t));
        }
        os << std::endl;

        t = update_time + (failureReport ? failureReport->update_time : 0.0);
        os << fmt::format("    Props/update time:          {:7.2f} s", t);
        if (failureReport) {
            os << fmt::format(" (Wasted: {:2.1f} s; {:2.1f}%)",
                            failureReport->update_time,
                            100*failureReport->update_time/noZero(t));
        }
        os << std::endl;

        int n = total_linearizations + (failureReport ? failureReport->total_linearizations : 0);
        os << fmt::format("Overall Linearizations:       {:7}", n);
        if (failureReport) {
          os << fmt::format("   (Wasted: {:5}; {:2.1f}%)",
                            failureReport->total_linearizations,
                            100.0*failureReport->total_linearizations/noZero(n));
        }
        os << std::endl;

        n = total_newton_iterations + (failureReport ? failureReport->total_newton_iterations : 0);
        os << fmt::format("Overall Nonlinear Iterations: {:7}", n);
        if (failureReport) {
          os << fmt::format("   (Wasted: {:5}; {:2.1f}%)",
                            failureReport->total_newton_iterations,
                            100.0*failureReport->total_newton_iterations/noZero(n));
        }
        os << std::endl;

        n = total_linear_iterations + (failureReport ? failureReport->total_linear_iterations : 0);
        os << fmt::format("Overall Linear Iterations:    {:7}", n);
        if (failureReport) {
          os << fmt::format("   (Wasted: {:5}; {:2.1f}%)",
                            failureReport->total_linear_iterations,
                            100.0*failureReport->total_linear_iterations/noZero(n));
        }
        os << std::endl;
        os << fmt::format("-------------------------------------------------------\n");
        n = skipped_domains + (failureReport ? failureReport->skipped_domains : 0);
        os << fmt::format("Skipped domain solves:       {:7}", n);
        os << std::endl;
        n = converged_domains + (failureReport ? failureReport->converged_domains : 0);
        os << fmt::format("Converged domain solves:     {:7}", n);
        os << std::endl;
        n = accepted_unconverged_domains + (failureReport ? failureReport->accepted_unconverged_domains : 0);
        os << fmt::format("  Accepted with relaxed tol: {:7}", n);
        os << std::endl;
        n = unconverged_domains + (failureReport ? failureReport->unconverged_domains : 0);
        os << fmt::format("Unconverged domain solves:   {:7}", n);
        os << std::endl << std::endl;

    }

    SimulatorReport SimulatorReport::serializationTestObject()
    {
        return SimulatorReport{SimulatorReportSingle::serializationTestObject(),
                               SimulatorReportSingle::serializationTestObject(),
                               {SimulatorReportSingle::serializationTestObject()}};
    }

    bool SimulatorReport::operator==(const SimulatorReport& rhs) const
    {
        return this->success == rhs.success &&
               this->failure == rhs.failure &&
               this->stepreports == rhs.stepreports;
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

    void SimulatorReport::reportFullyImplicit(std::ostream& os, bool performance_details) const
    {
        os << fmt::format("Number of timesteps:     {:9}\n", stepreports.size());
        success.reportFullyImplicit(os, &failure, performance_details);
    }

    void SimulatorReport::reportNLDD(std::ostream& os) const
    {
        success.reportNLDD(os, &failure);
    }

    void SimulatorReport::fullReports(std::ostream& os) const
    {
        os << "  Time(day)  TStep(day)  Assembly    LSetup    LSolve    LocSol    Update    Output WellIt Lins NewtIt LinIt Conv\n";
        for (std::size_t i = 0; i < this->stepreports.size(); ++i) {
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
            os << std::setw(9) << sr.local_solve_time << " ";
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
