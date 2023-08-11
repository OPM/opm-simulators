/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS

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

#ifndef OPM_AQUIFERINTERFACE_HEADER_INCLUDED
#define OPM_AQUIFERINTERFACE_HEADER_INCLUDED

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>

#include <opm/output/data/Aquifer.hpp>

namespace Opm
{

template <typename TypeTag>
class AquiferInterface
{
public:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    // Constructor
    AquiferInterface(int aqID,
                     const Simulator& ebosSimulator)
        : aquiferID_(aqID)
        , ebos_simulator_(ebosSimulator)
    {
    }

    // Destructor
    virtual ~AquiferInterface() = default;

    virtual void initFromRestart(const data::Aquifers& aquiferSoln) = 0;

    virtual void initialSolutionApplied() = 0;

    virtual void beginTimeStep() = 0;
    virtual void endTimeStep() = 0;

    virtual data::AquiferData aquiferData() const = 0;

    virtual void computeFaceAreaFraction(const std::vector<double>& total_face_area) = 0;
    virtual double totalFaceArea() const = 0;

    template <class Context>
    void addToSource(RateVector& rates,
                     const Context& context,
                     const unsigned spaceIdx,
                     const unsigned timeIdx)
    {
        const unsigned cellIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        addToSource(rates, cellIdx, timeIdx);
    }

    virtual void addToSource(RateVector& rates,
                             const unsigned cellIdx,
                             const unsigned timeIdx) = 0;

    int aquiferID() const { return this->aquiferID_; }

protected:
    bool co2store_or_h2store_() const
    {
        const auto& rspec = ebos_simulator_.vanguard().eclState().runspec();
        return rspec.co2Storage() || rspec.h2Storage();
    }

    int phaseIdx_() const
    {
        // If OIL is used to model brine the aquifer should do the same
        if (co2store_or_h2store_() && FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx))
            return FluidSystem::oilPhaseIdx;

        return FluidSystem::waterPhaseIdx;
    }

    const int aquiferID_{};
    const Simulator& ebos_simulator_;
};

} // namespace Opm

#endif
