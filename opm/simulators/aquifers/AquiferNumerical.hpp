/*
  Copyright (C) 2020 Equinor ASA
  Copyright (C) 2020 SINTEF Digital

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

#ifndef OPM_AQUIFERNUMERICAL_HEADER_INCLUDED
#define OPM_AQUIFERNUMERICAL_HEADER_INCLUDED

#include <opm/output/data/Aquifer.hpp>

#include <opm/parser/eclipse/EclipseState/NumericalAquifer.hpp>

namespace Opm
{
template <typename TypeTag>
class AquiferNumerical
{
public:
    // Constructor
    AquiferNumerical(const SingleNumericalAquifer& aquifer)
    : aquifer_(aquifer)
    , init_pressure_(aquifer.initPressure())
    , pressure_(this->init_pressure_)
    , flux_rate_(0.)
    , cumulative_flux_(0.)
    {
    }

    void initFromRestart(const std::vector<data::AquiferData>& aquiferSoln)
    {
        // NOT handling Restart for now
    }

    void beginTimeStep()
    {
    }

    void endTimeStep()
    {
    }

    Opm::data::AquiferData aquiferData() const
    {
        data::AquiferData data;
        data.aquiferID = this->aquifer_.id();
        data.initPressure = this->init_pressure_;
        data.pressure = this->pressure_;
        data.volume = this->cumulative_flux_;
        data.type = Opm::data::AquiferType::Numerical;
        return data;
    }

private:
    const Opm::SingleNumericalAquifer& aquifer_;
    double init_pressure_;
    double pressure_; // aquifer pressure
    double flux_rate_; // aquifer influx rate
    double cumulative_flux_; // cumulative aquifer influx
};
} // namespace Opm
#endif
