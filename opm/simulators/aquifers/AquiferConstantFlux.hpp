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

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <cassert>
#include <numeric>
#include <vector>

namespace Opm {

template<typename TypeTag>
class AquiferConstantFlux : public AquiferInterface<TypeTag>
{
public:
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BlackoilIndices = GetPropType<TypeTag, Properties::Indices>;

    static constexpr int numEq = BlackoilIndices::numEq;
    using Eval = DenseAd::Evaluation<double, /*size=*/numEq>;

    AquiferConstantFlux(const std::vector<Aquancon::AquancCell>& connections,
                        const Simulator&                         ebos_simulator,
                        const SingleAquiferFlux&                 aquifer)
        : AquiferInterface<TypeTag>(aquifer.id, ebos_simulator)
        , connections_             (connections)
        , aquifer_data_            (aquifer)
        , connection_flux_         (connections_.size(), Eval{0})
    {
        this->total_face_area_ = this->initializeConnections();
    }

    static AquiferConstantFlux serializationTestObject(const Simulator& ebos_simulator)
    {
        AquiferConstantFlux<TypeTag> result({}, ebos_simulator, {});
        result.cumulative_flux_ = 1.0;

        return result;
    }

    virtual ~AquiferConstantFlux() = default;

    void computeFaceAreaFraction(const std::vector<double>& total_face_area) override
    {
        assert (total_face_area.size() >= static_cast<std::vector<double>::size_type>(this->aquiferID()));

        this->area_fraction_ = this->totalFaceArea()
            / total_face_area[this->aquiferID() - 1];
    }

    double totalFaceArea() const override
    {
        return this->total_face_area_;
    }

    void updateAquifer(const SingleAquiferFlux& aquifer)
    {
        aquifer_data_ = aquifer;
    }

    void initFromRestart(const data::Aquifers& aquiferSoln) override
    {
        auto xaqPos = aquiferSoln.find(this->aquiferID());
        if (xaqPos == aquiferSoln.end()) {
            return;
        }

        this->cumulative_flux_ = this->area_fraction_ * xaqPos->second.volume;
    }

    void initialSolutionApplied() override
    {}

    void beginTimeStep() override
    {}

    void endTimeStep() override
    {
        this->flux_rate_ = this->totalFluxRate();
        this->cumulative_flux_ +=
            this->flux_rate_ * this->ebos_simulator_.timeStepSize();
    }

    data::AquiferData aquiferData() const override
    {
        data::AquiferData data;

        data.aquiferID = this->aquifer_data_.id;

        // Pressure for constant flux aquifer is 0
        data.pressure = 0.0;
        data.fluxRate = this->totalFluxRate();

        data.volume = this->cumulative_flux_;

        // not totally sure whether initPressure matters
        data.initPressure = 0.0;

        return data;
    }

    void addToSource(RateVector& rates,
                     const unsigned cellIdx,
                     [[maybe_unused]] const unsigned timeIdx) override
    {
        const int idx = this->cellToConnectionIdx_[cellIdx];
        if (idx < 0) {
            return;
        }

        const auto& model = this->ebos_simulator_.model();

        const auto fw = this->aquifer_data_.flux;

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
    std::vector<Eval> connection_flux_{};
    std::vector<int> cellToConnectionIdx_{};
    double flux_rate_{};
    double cumulative_flux_{};
    double total_face_area_{0.0};
    double area_fraction_{1.0};

    double initializeConnections()
    {
        auto connected_face_area = 0.0;

        this->cellToConnectionIdx_
            .resize(this->ebos_simulator_.gridView().size(/*codim=*/0), -1);

        for (std::size_t idx = 0; idx < this->connections_.size(); ++idx) {
            const auto global_index = this->connections_[idx].global_index;
            const int cell_index = this->ebos_simulator_.vanguard()
                .compressedIndexForInterior(global_index);

            if (cell_index < 0) {
                continue;
            }

            this->cellToConnectionIdx_[cell_index] = idx;

            connected_face_area += this->connections_[idx].effective_facearea;
        }

        // TODO: At the moment, we are using the effective_facearea from the
        // parser.  Should we update the facearea here if the grid changed
        // during the preprocessing?

        return connected_face_area;
    }

    double computeFaceAreaFraction(const double connected_face_area) const
    {
        const auto tot_face_area = this->ebos_simulator_.vanguard()
            .grid().comm().sum(connected_face_area);

        return (tot_face_area > 0.0)
            ? connected_face_area / tot_face_area
            : 0.0;
    }

    // TODO: this is a function from AquiferAnalytical
    int compIdx_() const
    {
        if (this->co2store_or_h2store_())
            return FluidSystem::oilCompIdx;

        return FluidSystem::waterCompIdx;
    }

    double totalFluxRate() const
    {
        return std::accumulate(this->connection_flux_.begin(),
                               this->connection_flux_.end(), 0.0,
                               [](const double rate, const auto& q)
                               {
                                   return rate + getValue(q);
                               });
    }
};

} // namespace Opm

#endif //OPM_AQUIFERCONSTANTFLUX_HPP
