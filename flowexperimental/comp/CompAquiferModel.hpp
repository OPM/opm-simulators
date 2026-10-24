/*
  Copyright 2026 SINTEF Digital

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

#ifndef OPM_COMP_AQUIFER_MODEL_HPP
#define OPM_COMP_AQUIFER_MODEL_HPP

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/AquiferConfig.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquifers.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <opm/simulators/aquifers/AquiferNumerical.hpp>

#include <memory>
#include <vector>

namespace Opm
{

/// Aquifer model of the compositional simulator, for numerical aquifers only.
/// Their cells are part of the grid, so the aquifers add no source terms and
/// only report the water that crosses their connections.
template <typename TypeTag>
class CompAquiferModel
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

public:
    explicit CompAquiferModel(Simulator& simulator)
    {
        const auto& aquiferConfig = simulator.vanguard().eclState().aquifer();
        if (aquiferConfig.hasAnalyticalAquifer()) {
            OpmLog::warning("Analytic aquifers are not supported in compositional "
                            "simulations and are ignored.");
        }

        for (const auto& [id, aquifer] : aquiferConfig.numericalAquifers().aquifers()) {
            aquifers_.push_back(std::make_unique<AquiferNumerical<TypeTag>>(aquifer, simulator));
        }
    }

    void initialSolutionApplied()
    {
        for (auto& aquifer : aquifers_) {
            aquifer->initialSolutionApplied();
        }
    }

    void beginEpisode()
    {
    }
    void beginTimeStep()
    {
    }
    void beginIteration()
    {
    }

    template <class RateVector, class Context>
    void addToSource(RateVector&, const Context&, unsigned, unsigned) const
    {
    }

    template <class RateVector>
    void addToSource(RateVector&, unsigned, unsigned) const
    {
    }

    void endIteration()
    {
    }

    void endTimeStep()
    {
        for (auto& aquifer : aquifers_) {
            aquifer->endTimeStep();
        }
    }

    void endEpisode()
    {
    }

    data::Aquifers aquiferData() const
    {
        data::Aquifers data;
        for (const auto& aquifer : aquifers_) {
            data.insert_or_assign(aquifer->aquiferID(), aquifer->aquiferData());
        }
        return data;
    }

    template <class Restarter>
    void serialize(Restarter&)
    {
    }

    template <class Restarter>
    void deserialize(Restarter&)
    {
    }

    template <class Serializer>
    void serializeOp(Serializer& serializer)
    {
        for (auto& aquifer : aquifers_) {
            serializer(*aquifer);
        }
    }

private:
    std::vector<std::unique_ptr<AquiferNumerical<TypeTag>>> aquifers_;
};

} // namespace Opm

#endif // OPM_COMP_AQUIFER_MODEL_HPP
