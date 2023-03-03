/*
  Copyright (C) 2023 Equinor

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

#ifndef OPM_AQUIFERCONSTANTFLUX_HPP
#define OPM_AQUIFERCONSTANTFLUX_HPP

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <opm/simulators/aquifers/AquiferInterface.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/Aquancon.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/AquiferFlux.hpp>

namespace Opm {
template<typename TypeTag>
class AquiferConstantFlux : public AquiferInterface<TypeTag> {
public:
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BlackoilIndices = GetPropType<TypeTag, Properties::Indices>;
    static constexpr int numEq = BlackoilIndices::numEq;
    using Eval = DenseAd::Evaluation<double, /*size=*/numEq>;

    AquiferConstantFlux(const SingleAquiferFlux& aquifer,
                        const std::vector<Aquancon::AquancCell>& connections,
                        const Simulator& ebos_simulator)
        : AquiferInterface<TypeTag>(aquifer.id, ebos_simulator)
         , connections_(connections)
         , aquifer_data_(aquifer)
    {
        // init_cumulative_flux is the flux volume from previoius running
        this->initializeConnections();
        connection_flux_.resize(this->connections_.size(), {0});
    }

    static AquiferConstantFlux serializationTestObject(const Simulator& ebos_simulator)
    {
        AquiferConstantFlux<TypeTag> result({}, {}, ebos_simulator);
        result.cumulative_flux_ = 1.0;

        return result;
    }

    virtual ~AquiferConstantFlux() = default;

    void updateAquifer(const SingleAquiferFlux& aquifer) {
        aquifer_data_ = aquifer;
    }

    void initFromRestart(const data::Aquifers& /* aquiferSoln */) override {
    }

    void initialSolutionApplied() override {
    }

    void beginTimeStep() override {
    }

    void endTimeStep() override {
        this->flux_rate_ = 0.;
        for (const auto& q : this->connection_flux_) {
            this->flux_rate_ += Opm::getValue(q);
        }

        this->cumulative_flux_ += this->flux_rate_ * this->ebos_simulator_.timeStepSize();
    }

    data::AquiferData aquiferData() const override
    {
        data::AquiferData data;
        data.aquiferID = this->aquifer_data_.id;
        // pressure for constant flux aquifer is 0
        data.pressure = 0.;
        data.fluxRate = 0.;
        for (const auto& q : this->connection_flux_) {
            data.fluxRate += q.value();
        }
        data.volume = this->cumulative_flux_;
        // not totally sure whether initPressure matters
        data.initPressure = 0.;
        return data;
    }

    void addToSource(RateVector& rates,
                     const unsigned cellIdx,
                     const unsigned timeIdx) override
    {
        const auto& model = this->ebos_simulator_.model();

        const int idx = this->cellToConnectionIdx_[cellIdx];
        if (idx < 0)
            return;

        const auto* intQuantsPtr = model.cachedIntensiveQuantities(cellIdx, timeIdx);
        if (intQuantsPtr == nullptr) {
            throw std::logic_error("Invalid intensive quantities cache detected in AquiferAnalytical::addToSource()");
        }

        const double fw = this->aquifer_data_.flux;
        this->connection_flux_[idx] = fw * this->connections_[idx].effective_facearea;
        rates[BlackoilIndices::conti0EqIdx + compIdx_()]
                += this->connection_flux_[idx] / model.dofTotalVolume(cellIdx);
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(cumulative_flux_);
    }

    bool operator==(const AquiferConstantFlux& rhs) const
    {
        return this->cumulative_flux_ == rhs.cumulative_flux_;
    }

private:
    const std::vector<Aquancon::AquancCell>& connections_;
    SingleAquiferFlux aquifer_data_;
    std::vector<int> cellToConnectionIdx_;
    std::vector<Eval> connection_flux_;
    double flux_rate_ {};
    double cumulative_flux_ = 0.;

    void initializeConnections() {
        this->cellToConnectionIdx_.resize(this->ebos_simulator_.gridView().size(/*codim=*/0), -1);
        for (std::size_t idx = 0; idx < this->connections_.size(); ++idx) {
            const auto global_index = this->connections_[idx].global_index;
            const int cell_index = this->ebos_simulator_.vanguard().compressedIndexForInterior(global_index);

            if (cell_index < 0) {
                continue;
            }

            this->cellToConnectionIdx_[cell_index] = idx;
        }
        // TODO: at the moment, we are using the effective_facearea from the parser. Should we update the facearea here if
        //  the grid changed during the preprocessing?
    }

    // TODO: this is a function from AquiferAnalytical
    int compIdx_() const
    {
        if (this->co2store_())
            return FluidSystem::oilCompIdx;

        return FluidSystem::waterCompIdx;
    }
};
}

#endif //OPM_AQUIFERCONSTANTFLUX_HPP
