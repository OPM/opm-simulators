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


    // TODO: we need to pass in the previous flux volume
    AquiferConstantFlux(const std::shared_ptr<AquiferFlux>& aquifer,
                        const std::vector<Aquancon::AquancCell>& connections,
                        const Simulator& ebos_simulator,
                        const double init_cumulative_flux = 0.)
        : AquiferInterface<TypeTag>(aquifer->id, ebos_simulator)
         , connections_(connections)
         , aquifer_data_(aquifer)
         , cumulative_flux_(init_cumulative_flux)
    {
        // flux_volume is the flux volume from previoius running
        this->initializeConnections();
        connection_flux_.resize(this->connections_.size(), {0});
    }

    virtual ~AquiferConstantFlux() = default;

    /* void updateAquifer(const std::shared_ptr<AquiferFlux>& aquifer) {
        aquifer_data_ = aquifer;
    } */

    void initFromRestart(const data::Aquifers& /* aquiferSoln */) {
    }

    void initialSolutionApplied() {
        // Note: we can not do this here
        // with the current way, we put the AQUFLUX in the first report step of Schedule
        // Maybe it might bring some undesiable consequence to remove it from the solution
    }

    void beginTimeStep() {
    }

    void endTimeStep() {
        this->flux_rate_ = 0.;
        for (const auto& q : this->connection_flux_) {
            this->flux_rate_ += Opm::getValue(q);
        }

        const auto& comm = this->ebos_simulator_.vanguard().grid().comm();
        comm.sum(&this->flux_rate_, 1);

        this->cumulative_flux_ += this->flux_rate_ * this->ebos_simulator_.timeStepSize();
    }

    data::AquiferData aquiferData() const
    {
        data::AquiferData data;
        data.aquiferID = this->aquifer_data_->id;
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
                     const unsigned timeIdx) {
        const auto& model = this->ebos_simulator_.model();

        const int idx = this->cellToConnectionIdx_[cellIdx];
        if (idx < 0)
            return;

        const auto* intQuantsPtr = model.cachedIntensiveQuantities(cellIdx, timeIdx);
        if (intQuantsPtr == nullptr) {
            throw std::logic_error("Invalid intensive quantities cache detected in AquiferAnalytical::addToSource()");
        }

        const double fw = this->aquifer_data_->flux;
        // const double m = this->connections_[idx].influx_coeff;
        this->connection_flux_[idx] = fw * this->connections_[idx].effective_facearea;
        rates[BlackoilIndices::conti0EqIdx + compIdx_()]
                += this->connection_flux_[idx] / model.dofTotalVolume(cellIdx);
    }

    // TODO: repeated function from AquiferAnalytical
    std::size_t size() const
    {
        return this->connections_.size();
    }

private:
    const std::vector<Aquancon::AquancCell> connections_;
    std::shared_ptr<AquiferFlux> aquifer_data_;
    // TODO: for simple case, faceArea_connected_ is not needed here, since it is calculated when parsing
    // But if the grid change, not sure how to handle
    // std::vector<double> faceArea_connected_;
    std::vector<int> cellToConnectionIdx_;
    std::vector<Eval> connection_flux_;
    double flux_rate_;
    double cumulative_flux_ = 0.;


    void initializeConnections() {
        // this->faceArea_connected_.resize(this->size(), {0});

        this->cellToConnectionIdx_.resize(this->ebos_simulator_.gridView().size(/*codim=*/0), -1);
        const auto& gridView = this->ebos_simulator_.vanguard().gridView();
        for (std::size_t idx = 0; idx < this->size(); ++idx) {
            const auto global_index = this->connections_[idx].global_index;
            const int cell_index = this->ebos_simulator_.vanguard().compressedIndex(global_index);
            auto elemIt = gridView.template begin</*codim=*/ 0>();
            if (cell_index > 0)
                std::advance(elemIt, cell_index);

            //the global_index is not part of this grid
            if (cell_index < 0 || elemIt->partitionType() != Dune::InteriorEntity)
                continue;

            this->cellToConnectionIdx_[cell_index] = idx;
        }
        // TODO: at the moment, we are using the effective_facearea from the parser. Should we update the facearea here?
    }
    // TODO: function from AquiferAnalytical
    int compIdx_() const
    {
        if (this->co2store_())
            return FluidSystem::oilCompIdx;

        return FluidSystem::waterCompIdx;
    }

    double cumulativeFlux() const
    {
        return this->cumulative_flux_;
    }
};
}

#endif //OPM_AQUIFERCONSTANTFLUX_HPP
