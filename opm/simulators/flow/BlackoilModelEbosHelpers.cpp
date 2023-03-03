/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

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
#include <opm/simulators/flow/BlackoilModelEbosHelpers.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>

#include <cmath>
#include <iomanip>

namespace Opm {
namespace detail {

template<class Scalar>
ConvergenceReport createConvergenceReport(const int iteration,
                                          const double reportTime,
                                          const double tol_mb,
                                          const double tol_cnv,
                                          const double maxResidualAllowed,
                                          const bool terminal_output,
                                          const std::vector<std::string>& names,
                                          const std::vector<Scalar>& CNV,
                                          const std::vector<Scalar>& mass_balance_residual,
                                          const std::string& terminal_header)
{
    // Create convergence report.
    ConvergenceReport report{reportTime};
    using CR = ConvergenceReport;
    const int numComp = CNV.size();
    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
        double res[2] = { mass_balance_residual[compIdx], CNV[compIdx] };
        CR::ReservoirFailure::Type types[2] = { CR::ReservoirFailure::Type::MassBalance,
                                                CR::ReservoirFailure::Type::Cnv };
        double tol[2] = { tol_mb, tol_cnv };
        for (int ii : {0, 1}) {
            if (std::isnan(res[ii])) {
                report.setReservoirFailed({types[ii], CR::Severity::NotANumber, compIdx});
                if (terminal_output) {
                    OpmLog::debug("NaN residual for " + names[compIdx] + " equation.");
                }
            } else if (res[ii] > maxResidualAllowed) {
                report.setReservoirFailed({types[ii], CR::Severity::TooLarge, compIdx});
                if (terminal_output) {
                    OpmLog::debug("Too large residual for " + names[compIdx] + " equation.");
                }
            } else if (res[ii] < 0.0) {
                report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
                if (terminal_output) {
                    OpmLog::debug("Negative residual for " + names[compIdx] + " equation.");
                }
            } else if (res[ii] > tol[ii]) {
                report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
            }
            report.setReservoirConvergenceMetric(types[ii], compIdx, res[ii]);
        }
    }

    // Output of residuals.
    if (terminal_output) {
        // Only rank 0 does print to std::cout
        if (iteration == 0) {
            std::string msg = terminal_header + "Iter";
            for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                msg += "    MB(";
                msg += names[compIdx][0];
                msg += ")  ";
            }
            for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                msg += "    CNV(";
                msg += names[compIdx][0];
                msg += ") ";
            }
            OpmLog::debug(msg);
        }
        std::ostringstream ss;
        if (!terminal_header.empty()) {
            ss << "| ";
        }
        const std::streamsize oprec = ss.precision(3);
        const std::ios::fmtflags oflags = ss.setf(std::ios::scientific);
        ss << std::setw(4) << iteration;
        for (int compIdx = 0; compIdx < numComp; ++compIdx) {
            ss << std::setw(11) << mass_balance_residual[compIdx];
        }
        for (int compIdx = 0; compIdx < numComp; ++compIdx) {
            ss << std::setw(11) << CNV[compIdx];
        }
        ss.precision(oprec);
        ss.flags(oflags);
        OpmLog::debug(ss.str());
    }

    return report;
}

template<class Scalar>
std::tuple<double,double> convergenceReduction(Parallel::Communication comm,
                                               const double pvSumLocal,
                                               const double numAquiferPvSumLocal,
                                               std::vector<Scalar>& R_sum,
                                               std::vector<Scalar>& maxCoeff,
                                               std::vector<Scalar>& B_avg)
{
    OPM_TIMEBLOCK(convergenceReduction);
    // Compute total pore volume (use only owned entries)
    double pvSum = pvSumLocal;
    double numAquiferPvSum = numAquiferPvSumLocal;

    if (comm.size() > 1) {
        // global reduction
        std::vector<Scalar> sumBuffer;
        std::vector<Scalar> maxBuffer;
        const int numComp = B_avg.size();
        sumBuffer.reserve( 2*numComp + 2 ); // +2 for (numAquifer)pvSum
        maxBuffer.reserve( numComp );
        for (int compIdx = 0; compIdx < numComp; ++compIdx) {
            sumBuffer.push_back(B_avg[compIdx]);
            sumBuffer.push_back(R_sum[compIdx]);
            maxBuffer.push_back(maxCoeff[compIdx]);
        }

        // Compute total pore volume
        sumBuffer.push_back(pvSum);
        sumBuffer.push_back(numAquiferPvSum);

        // compute global sum
        comm.sum(sumBuffer.data(), sumBuffer.size());

        // compute global max
        comm.max(maxBuffer.data(), maxBuffer.size());

        // restore values to local variables
        for (int compIdx = 0, buffIdx = 0; compIdx < numComp; ++compIdx, ++buffIdx) {
            B_avg[compIdx] = sumBuffer[buffIdx];
            ++buffIdx;

            R_sum[compIdx] = sumBuffer[buffIdx];
        }

        for (int compIdx = 0; compIdx < numComp; ++compIdx) {
            maxCoeff[compIdx] = maxBuffer[compIdx];
        }

        // restore global pore volume
        pvSum = sumBuffer[sumBuffer.size() - 2];
        numAquiferPvSum = sumBuffer.back();
    }

    // return global pore volume
    return {pvSum, numAquiferPvSum};
}

template ConvergenceReport
createConvergenceReport<double>(const int,
                                const double,
                                const double,
                                const double,
                                const double,
                                const bool,
                                const std::vector<std::string>&,
                                const std::vector<double>&,
                                const std::vector<double>&,
                                const std::string&);

template std::tuple<double,double>
convergenceReduction<double>(Parallel::Communication,
                             const double,
                             const double,
                             std::vector<double>&,
                             std::vector<double>&,
                             std::vector<double>&);

} // namespace detail
} // namespace Opm
