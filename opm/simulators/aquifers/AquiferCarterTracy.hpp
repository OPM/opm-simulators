/*
  Copyright 2017 TNO - Heat Transfer & Fluid Dynamics, Modelling & Optimization of the Subsurface
  Copyright 2017 Statoil ASA.

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

#ifndef OPM_AQUIFERCT_HEADER_INCLUDED
#define OPM_AQUIFERCT_HEADER_INCLUDED

#include <opm/simulators/aquifers/AquiferInterface.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <exception>
#include <memory>
#include <stdexcept>
#include <utility>

namespace Opm
{

template <typename TypeTag>
class AquiferCarterTracy : public AquiferInterface<TypeTag>
{
public:
    typedef AquiferInterface<TypeTag> Base;

    using typename Base::BlackoilIndices;
    using typename Base::ElementContext;
    using typename Base::Eval;
    using typename Base::FluidState;
    using typename Base::FluidSystem;
    using typename Base::IntensiveQuantities;
    using typename Base::RateVector;
    using typename Base::Scalar;
    using typename Base::Simulator;
    using typename Base::ElementMapper;

    using Base::waterCompIdx;
    using Base::waterPhaseIdx;
    AquiferCarterTracy(const std::vector<Aquancon::AquancCell>& connections,
                       const Simulator& ebosSimulator,
                       const AquiferCT::AQUCT_data& aquct_data)
        : Base(aquct_data.aquiferID, connections, ebosSimulator)
        , aquct_data_(aquct_data)
    {}

    void endTimeStep() override
    {
        for (const auto& q : this->Qai_) {
            this->W_flux_ += q * this->ebos_simulator_.timeStepSize();
        }
        this->fluxValue_ = this->W_flux_.value();
        const auto& comm = this->ebos_simulator_.vanguard().grid().comm();
        comm.sum(&this->fluxValue_, 1);
    }

    data::AquiferData aquiferData() const
    {
        data::AquiferData data;
        data.aquiferID = this->aquiferID();
        // TODO: not sure how to get this pressure value yet
        data.pressure = this->pa0_;
        data.fluxRate = 0.;
        for (const auto& q : this->Qai_) {
            data.fluxRate += q.value();
        }
        data.volume = this->W_flux_.value();
        data.initPressure = this->pa0_;

        auto* aquCT = data.typeData.template create<data::AquiferType::CarterTracy>();
        aquCT->timeConstant = this->aquct_data_.timeConstant();
        aquCT->influxConstant = this->aquct_data_.influxConstant();
        aquCT->waterDensity = this->aquct_data_.waterDensity();
        aquCT->waterViscosity = this->aquct_data_.waterViscosity();
        aquCT->dimensionless_time = this->dimensionless_time_;
        aquCT->dimensionless_pressure = this->dimensionless_pressure_;

        return data;
    }

protected:
    // Variables constants
    AquiferCT::AQUCT_data aquct_data_;

    Scalar beta_; // Influx constant
    // TODO: it is possible it should be a AD variable
    Scalar fluxValue_{0}; // value of flux

    Scalar dimensionless_time_{0};
    Scalar dimensionless_pressure_{0};

    void assignRestartData(const data::AquiferData& xaq) override
    {
        this->fluxValue_ = xaq.volume;
        this->rhow_ = this->aquct_data_.waterDensity();
    }

    std::pair<Scalar, Scalar>
    getInfluenceTableValues(const Scalar td_plus_dt)
    {
        // We use the opm-common numeric linear interpolator
        this->dimensionless_pressure_ =
            linearInterpolation(this->aquct_data_.dimensionless_time,
                                this->aquct_data_.dimensionless_pressure,
                                this->dimensionless_time_);

        const auto PItd =
            linearInterpolation(this->aquct_data_.dimensionless_time,
                                this->aquct_data_.dimensionless_pressure,
                                td_plus_dt);

        const auto PItdprime =
            linearInterpolationDerivative(this->aquct_data_.dimensionless_time,
                                          this->aquct_data_.dimensionless_pressure,
                                          td_plus_dt);

        return std::make_pair(PItd, PItdprime);
    }

    Scalar dpai(const int idx) const
    {
        const auto gdz =
            this->gravity_() * (this->cell_depth_.at(idx) - this->aquiferDepth());

        const auto dp = this->pa0_ + this->rhow_*gdz
            - this->pressure_previous_.at(idx);

        return dp;
    }

    // This function implements Eqs 5.8 and 5.9 of the EclipseTechnicalDescription
    std::pair<Scalar, Scalar>
    calculateEqnConstants(const int idx, const Simulator& simulator)
    {
        const Scalar td_plus_dt = (simulator.timeStepSize() + simulator.time()) / this->Tc_;
        this->dimensionless_time_ = simulator.time() / this->Tc_;

        const auto [PItd, PItdprime] = this->getInfluenceTableValues(td_plus_dt);

        const auto denom = this->Tc_ * (PItd - this->dimensionless_time_*PItdprime);
        const auto a = (this->beta_*dpai(idx) - this->fluxValue_*PItdprime) / denom;
        const auto b = this->beta_ / denom;

        return std::make_pair(a, b);
    }

    // This function implements Eq 5.7 of the EclipseTechnicalDescription
    inline void calculateInflowRate(int idx, const Simulator& simulator) override
    {
        const auto [a, b] = this->calculateEqnConstants(idx, simulator);

        this->Qai_.at(idx) = this->alphai_.at(idx) *
            (a - b*(this->pressure_current_.at(idx) - this->pressure_previous_.at(idx)));
    }

    inline void calculateAquiferConstants() override
    {
        this->Tc_ = this->aquct_data_.timeConstant();
        this->beta_ = this->aquct_data_.influxConstant();
    }

    inline void calculateAquiferCondition() override
    {
        if (this->solution_set_from_restart_) {
            return;
        }

        if (! this->aquct_data_.initial_pressure.has_value()) {
            this->aquct_data_.initial_pressure =
                this->calculateReservoirEquilibrium();

            const auto& tables = this->ebos_simulator_.vanguard()
                .eclState().getTableManager();

            this->aquct_data_.finishInitialisation(tables);
        }

        this->pa0_ = this->aquct_data_.initial_pressure.value();
        this->rhow_ = this->aquct_data_.waterDensity();
    }

    virtual Scalar aquiferDepth() const override
    {
        return this->aquct_data_.datum_depth;
    }
}; // class AquiferCarterTracy
} // namespace Opm

#endif
