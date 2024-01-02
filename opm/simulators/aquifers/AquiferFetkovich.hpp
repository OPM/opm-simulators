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

#ifndef OPM_AQUIFETP_HEADER_INCLUDED
#define OPM_AQUIFETP_HEADER_INCLUDED

#include <opm/simulators/aquifers/AquiferAnalytical.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/Aquifetp.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <stdexcept>
#include <vector>

namespace Opm
{

template <typename TypeTag>
class AquiferFetkovich : public AquiferAnalytical<TypeTag>
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

    AquiferFetkovich(const std::vector<Aquancon::AquancCell>& connections,
                     const Simulator& ebosSimulator,
                     const Aquifetp::AQUFETP_data& aqufetp_data)
        : Base(aqufetp_data.aquiferID, connections, ebosSimulator)
        , aqufetp_data_(aqufetp_data)
    {
    }

    static AquiferFetkovich serializationTestObject(const Simulator& ebosSimulator)
    {
        AquiferFetkovich result({}, ebosSimulator, {});

        result.pressure_previous_ = {1.0, 2.0, 3.0};
        result.pressure_current_ = {4.0, 5.0};
        result.Qai_ = {{6.0}};
        result.rhow_ = 7.0;
        result.W_flux_ = 8.0;
        result.aquifer_pressure_ = 9.0;

        return result;
    }

    void endTimeStep() override
    {
        for (const auto& q : this->Qai_) {
            this->W_flux_ += q * this->ebos_simulator_.timeStepSize();
        }
        aquifer_pressure_ = aquiferPressure();
    }

    data::AquiferData aquiferData() const override
    {
        // TODO: how to unify the two functions?
        auto data = data::AquiferData{};

        data.aquiferID = this->aquiferID();
        data.pressure = this->aquifer_pressure_;
        data.fluxRate = std::accumulate(this->Qai_.begin(), this->Qai_.end(), 0.0,
                                        [](const double flux, const auto& q) -> double
                                        {
                                            return flux + q.value();
                                        });
        data.volume = this->W_flux_.value();
        data.initPressure = this->pa0_;

        auto* aquFet = data.typeData.template create<data::AquiferType::Fetkovich>();
        aquFet->initVolume = this->aqufetp_data_.initial_watvolume;
        aquFet->prodIndex = this->aqufetp_data_.prod_index;
        aquFet->timeConstant = this->aqufetp_data_.timeConstant();

        return data;
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(static_cast<Base&>(*this));
        serializer(aquifer_pressure_);
    }

    bool operator==(const AquiferFetkovich& rhs) const
    {
        return static_cast<const Base&>(*this) == rhs &&
               this->aquifer_pressure_ == rhs.aquifer_pressure_;
    }

protected:
    // Aquifer Fetkovich Specific Variables
    Aquifetp::AQUFETP_data aqufetp_data_;
    Scalar aquifer_pressure_; // aquifer

    void assignRestartData(const data::AquiferData& xaq) override
    {
        if (! xaq.typeData.is<data::AquiferType::Fetkovich>()) {
            throw std::invalid_argument {
                "Analytic aquifer data for unexpected aquifer "
                "type passed to Fetkovich aquifer"
            };
        }

        this->aquifer_pressure_ = xaq.pressure;
        this->rhow_ = this->aqufetp_data_.waterDensity();
    }

    inline Eval dpai(int idx)
    {
        const auto gdz =
            this->gravity_() * (this->cell_depth_[idx] - this->aquiferDepth());

        return this->aquifer_pressure_ + this->rhow_*gdz
            - this->pressure_current_.at(idx);
    }

    // This function implements Eq 5.12 of the EclipseTechnicalDescription
    inline Scalar aquiferPressure()
    {
        Scalar Flux = this->W_flux_.value();

        const auto& comm = this->ebos_simulator_.vanguard().grid().comm();
        comm.sum(&Flux, 1);

        const auto denom =
            this->aqufetp_data_.total_compr * this->aqufetp_data_.initial_watvolume;

        return this->pa0_ - (Flux / denom);
    }

    inline void calculateAquiferConstants() override
    {
        this->Tc_ = this->aqufetp_data_.timeConstant();
    }

    // This function implements Eq 5.14 of the EclipseTechnicalDescription
    inline void calculateInflowRate(int idx, const Simulator& simulator) override
    {
        const Scalar td_Tc_ = simulator.timeStepSize() / this->Tc_;
        const Scalar coef = (1 - exp(-td_Tc_)) / td_Tc_;

        this->Qai_.at(idx) = coef * this->alphai_[idx] *
            this->aqufetp_data_.prod_index * dpai(idx);
    }

    inline void calculateAquiferCondition() override
    {
        if (this->solution_set_from_restart_) {
            return;
        }

        if (! this->aqufetp_data_.initial_pressure.has_value()) {
            this->aqufetp_data_.initial_pressure =
                this->calculateReservoirEquilibrium();

            const auto& tables = this->ebos_simulator_.vanguard()
                .eclState().getTableManager();

            this->aqufetp_data_.finishInitialisation(tables);
        }

        this->rhow_ = this->aqufetp_data_.waterDensity();
        this->pa0_ = this->aqufetp_data_.initial_pressure.value();
        if (this->aqufetp_data_.initial_temperature.has_value())
            this->Ta0_ = this->aqufetp_data_.initial_temperature.value();
        this->aquifer_pressure_ = this->pa0_;
    }

    virtual Scalar aquiferDepth() const override
    {
        return this->aqufetp_data_.datum_depth;
    }
}; // Class AquiferFetkovich

} // namespace Opm

#endif
