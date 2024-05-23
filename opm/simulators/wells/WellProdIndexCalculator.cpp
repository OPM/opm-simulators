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

template<class Scalar>
void checkSizeCompatibility(const Opm::WellProdIndexCalculator<Scalar>& wellPICalc,
                            const std::vector<Scalar>&                  connMobility)
{
    if (connMobility.size() != wellPICalc.numConnections()) {
        throw std::logic_error {
            "Mobility vector size does not match expected number of connections"
        };
    }
}

template<class Scalar>
Scalar logRescale(const Scalar r0, const Scalar rw, const Scalar rd, const Scalar S)
{
    const auto numerator = std::log(r0 / rw) + S;
    const auto denom     = std::log(rd / rw) + S;

    return numerator / denom;
}

template<class Scalar>
void standardConnFactorsExplicitDrainRadius(const Opm::Well&     well,
                                            std::vector<Scalar>& stdConnFact)
{
    const auto& connections = well.getConnections();
    const auto rdrain = well.getDrainageRadius();

    std::transform(connections.begin(), connections.end(), stdConnFact.begin(),
        [rdrain](const Opm::Connection& conn)
    {
        return conn.CF() * logRescale(conn.r0(), conn.rw(), rdrain, conn.skinFactor());
    });
}

template<class Scalar>
void standardConnFactorsDrainIsEquivalent(const Opm::Well&     well,
                                          std::vector<Scalar>& stdConnFact)
{
    const auto& connections = well.getConnections();

    std::transform(connections.begin(), connections.end(), stdConnFact.begin(),
        [](const Opm::Connection& conn)
    {
        return conn.CF();
    });
}

template<class Scalar>
std::vector<Scalar> calculateStandardConnFactors(const Opm::Well& well)
{
    std::vector<Scalar> stdConnFact(well.getConnections().size());

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

template<class Scalar>
Opm::WellProdIndexCalculator<Scalar>::
WellProdIndexCalculator(const Well& well)
    : standardConnFactors_{ calculateStandardConnFactors<Scalar>(well) }
{}

template<class Scalar>
void Opm::WellProdIndexCalculator<Scalar>::
reInit(const Well& well)
{
    this->standardConnFactors_ = calculateStandardConnFactors<Scalar>(well);
}

template<class Scalar>
Scalar Opm::WellProdIndexCalculator<Scalar>::
connectionProdIndStandard(const std::size_t connIdx,
                          const Scalar      connMobility) const
{
    return this->standardConnFactors_[connIdx] * connMobility;
}

// ===========================================================================

template<class Scalar>
std::vector<Scalar>
Opm::connectionProdIndStandard(const WellProdIndexCalculator<Scalar>& wellPICalc,
                               const std::vector<Scalar>&             connMobility)
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

template<class Scalar>
Scalar Opm::wellProdIndStandard(const WellProdIndexCalculator<Scalar>& wellPICalc,
                                const std::vector<Scalar>&             connMobility)
{
    const auto connPI = connectionProdIndStandard(wellPICalc, connMobility);

    return std::accumulate(connPI.begin(), connPI.end(), 0.0);
}

template class Opm::WellProdIndexCalculator<double>;
template std::vector<double>
Opm::connectionProdIndStandard(const WellProdIndexCalculator<double>&,
                               const std::vector<double>&);
template double
Opm::wellProdIndStandard(const WellProdIndexCalculator<double>&,
                         const std::vector<double>&);
