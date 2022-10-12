/*
  Copyright 2020 Equinor ASA.

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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <opm/simulators/wells/WellProdIndexCalculator.hpp>

#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
    void checkSizeCompatibility(const Opm::WellProdIndexCalculator& wellPICalc,
                                const std::vector<double>&          connMobility)
    {
        if (connMobility.size() != wellPICalc.numConnections()) {
            throw std::logic_error {
                "Mobility vector size does not match expected number of connections"
            };
        }
    }

    double logRescale(const double r0, const double rw, const double rd, const double S)
    {
        const auto numerator = std::log(r0 / rw) + S;
        const auto denom     = std::log(rd / rw) + S;

        return numerator / denom;
    }

    void standardConnFactorsExplicitDrainRadius(const Opm::Well&     well,
                                                std::vector<double>& stdConnFact)
    {
        const auto& connections = well.getConnections();
        const auto rdrain = well.getDrainageRadius();

        std::transform(connections.begin(), connections.end(), stdConnFact.begin(),
            [rdrain](const Opm::Connection& conn)
        {
            return conn.CF() * logRescale(conn.r0(), conn.rw(), rdrain, conn.skinFactor());
        });
    }

    void standardConnFactorsDrainIsEquivalent(const Opm::Well&     well,
                                              std::vector<double>& stdConnFact)
    {
        const auto& connections = well.getConnections();

        std::transform(connections.begin(), connections.end(), stdConnFact.begin(),
            [](const Opm::Connection& conn)
        {
            return conn.CF();
        });
    }

    std::vector<double> calculateStandardConnFactors(const Opm::Well& well)
    {
        std::vector<double> stdConnFact(well.getConnections().size());

        if (well.getDrainageRadius() > 0.0) {
            // Well has an explicit drainage radius.  Apply logarithmic
            // scaling to the CTFs.
            standardConnFactorsExplicitDrainRadius(well, stdConnFact);
        }
        else {
            // Unspecified drainage radius.  Standard mobility connection
            // scaling factors are just the regular CTFs.
            standardConnFactorsDrainIsEquivalent(well, stdConnFact);
        }

        return stdConnFact;
    }
} // namespace Anonymous

Opm::WellProdIndexCalculator::WellProdIndexCalculator(const Well& well)
    : standardConnFactors_{ calculateStandardConnFactors(well) }
{}

void Opm::WellProdIndexCalculator::reInit(const Well& well)
{
    this->standardConnFactors_ = calculateStandardConnFactors(well);
}

double
Opm::WellProdIndexCalculator::
connectionProdIndStandard(const std::size_t connIdx,
                          const double      connMobility) const
{
    return this->standardConnFactors_[connIdx] * connMobility;
}

// ===========================================================================

std::vector<double>
Opm::connectionProdIndStandard(const WellProdIndexCalculator& wellPICalc,
                               const std::vector<double>&     connMobility)
{
    checkSizeCompatibility(wellPICalc, connMobility);

    const auto nConn = wellPICalc.numConnections();
    auto connPI = connMobility;
    for (auto connIx = 0*nConn; connIx < nConn; ++connIx) {
        connPI[connIx] = wellPICalc
            .connectionProdIndStandard(connIx, connMobility[connIx]);
    }

    return connPI;
}

double Opm::wellProdIndStandard(const WellProdIndexCalculator& wellPICalc,
                                const std::vector<double>&     connMobility)
{
    const auto connPI = connectionProdIndStandard(wellPICalc, connMobility);

    return std::accumulate(connPI.begin(), connPI.end(), 0.0);
}
