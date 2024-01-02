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

#include <opm/common/utility/numeric/linearInterpolation.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/AquiferCT.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <opm/simulators/aquifers/AquiferAnalytical.hpp>

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Opm
{

template <typename TypeTag>
class AquiferCarterTracy : public AquiferAnalytical<TypeTag>
{
public:
    using Base = AquiferAnalytical<TypeTag>;

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

    AquiferCarterTracy(const std::vector<Aquancon::AquancCell>& connections,
                       const Simulator& ebosSimulator,
                       const AquiferCT::AQUCT_data& aquct_data)
        : Base(aquct_data.aquiferID, connections, ebosSimulator)
        , aquct_data_(aquct_data)
    {}

    static AquiferCarterTracy serializationTestObject(const Simulator& ebosSimulator)
    {
        AquiferCarterTracy result({}, ebosSimulator, {});

        result.pressure_previous_ = {1.0, 2.0, 3.0};
        result.pressure_current_ = {4.0, 5.0};
        result.Qai_ = {{6.0}};
        result.rhow_ = 7.0;
        result.W_flux_ = 8.0;
        result.fluxValue_ = 9.0;
        result.dimensionless_time_ = 10.0;
        result.dimensionless_pressure_ = 11.0;

        return result;
    }

    void endTimeStep() override
    {
        for (const auto& q : this->Qai_) {
            this->W_flux_ += q * this->ebos_simulator_.timeStepSize();
        }
        this->fluxValue_ = this->W_flux_.value();
        const auto& comm = this->ebos_simulator_.vanguard().grid().comm();
        comm.sum(&this->fluxValue_, 1);
    }

    data::AquiferData aquiferData() const override
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

        aquCT->dimensionless_time = this->dimensionless_time_;
        aquCT->dimensionless_pressure = this->dimensionless_pressure_;
        aquCT->influxConstant = this->aquct_data_.influxConstant();

        if (!this->co2store_or_h2store_()) {
            aquCT->timeConstant = this->aquct_data_.timeConstant();
            aquCT->waterDensity = this->aquct_data_.waterDensity();
            aquCT->waterViscosity = this->aquct_data_.waterViscosity();
        } else {
            aquCT->waterDensity = this->rhow_;
            aquCT->timeConstant = this->Tc_;
            const auto x = this->aquct_data_.porosity * this->aquct_data_.total_compr * this->aquct_data_.inner_radius * this->aquct_data_.inner_radius;
            aquCT->waterViscosity = this->Tc_ *  this->aquct_data_.permeability / x;
        }

        return data;
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<Base&>(*this));
        serializer(fluxValue_);
        serializer(dimensionless_time_);
        serializer(dimensionless_pressure_);
    }

    bool operator==(const AquiferCarterTracy& rhs) const
    {
        return static_cast<const AquiferAnalytical<TypeTag>&>(*this) == rhs &&
               this->fluxValue_ == rhs.fluxValue_ &&
               this->dimensionless_time_ == rhs.dimensionless_time_ &&
               this->dimensionless_pressure_ == rhs.dimensionless_pressure_;
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

    std::size_t pvtRegionIdx() const
    {
        return this->aquct_data_.pvttableID - 1;
    }

    // This function implements Eq 5.7 of the EclipseTechnicalDescription
    void calculateInflowRate(int idx, const Simulator& simulator) override
    {
        const auto [a, b] = this->calculateEqnConstants(idx, simulator);

        this->Qai_.at(idx) = this->alphai_.at(idx) *
            (a - b*(this->pressure_current_.at(idx) - this->pressure_previous_.at(idx)));
    }

    void calculateAquiferConstants() override
    {
        this->Tc_ = this->co2store_or_h2store_()
            ? this->timeConstantCO2Store()
            : this->aquct_data_.timeConstant();

        this->beta_ = this->aquct_data_.influxConstant();
    }

    void calculateAquiferCondition() override
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
        if (this->aquct_data_.initial_temperature.has_value()) {
            this->Ta0_ = this->aquct_data_.initial_temperature.value();
        }

        this->rhow_ = this->co2store_or_h2store_()
            ? this->waterDensityCO2Store()
            : this->aquct_data_.waterDensity();
    }

    Scalar aquiferDepth() const override
    {
        return this->aquct_data_.datum_depth;
    }

private:
    Scalar timeConstantCO2Store() const
    {
        const auto press = this->aquct_data_.initial_pressure.value();
        const auto temp = this->reservoirTemperatureCO2Store();

        auto waterViscosity = Scalar { 0 };
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const auto rs = Scalar { 0 }; // no dissolved CO2
            waterViscosity = FluidSystem::oilPvt()
                .viscosity(pvtRegionIdx(), temp, press, rs);
        }
        else if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const auto salt = Scalar { 0 };
            const auto rsw = Scalar { 0 };
            waterViscosity = FluidSystem::waterPvt()
                .viscosity(pvtRegionIdx(), temp, press, rsw, salt);
        }
        else {
            OPM_THROW(std::runtime_error, "water or oil phase is needed to run CO2Store.");
        }

        const auto x = this->aquct_data_.porosity * this->aquct_data_.total_compr
            * this->aquct_data_.inner_radius * this->aquct_data_.inner_radius;

        return waterViscosity * x / this->aquct_data_.permeability;
    }

    Scalar waterDensityCO2Store() const
    {
        const auto press = this->aquct_data_.initial_pressure.value();
        const auto temp = this->reservoirTemperatureCO2Store();

        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            const auto& pvt = FluidSystem::oilPvt();
            const auto reg = this->pvtRegionIdx();

            const auto rs = Scalar { 0 }; // no dissolved CO2
            return pvt.inverseFormationVolumeFactor(reg, temp, press, rs)
                * pvt.oilReferenceDensity(reg);
        }
        else if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const auto& pvt = FluidSystem::waterPvt();
            const auto reg = this->pvtRegionIdx();

            const auto salinity = Scalar { 0 };
            const auto rsw = Scalar { 0 };

            return pvt.inverseFormationVolumeFactor(reg, temp, press, rsw, salinity)
                * pvt.waterReferenceDensity(reg);
        }
        else {
            OPM_THROW(std::runtime_error, "water or oil phase is needed to run CO2Store.");
        }
    }

    Scalar reservoirTemperatureCO2Store() const
    {
        return this->aquct_data_.initial_temperature.has_value()
            ? this->aquct_data_.initial_temperature.value()
            : FluidSystem::reservoirTemperature();
    }

}; // class AquiferCarterTracy

} // namespace Opm

#endif
