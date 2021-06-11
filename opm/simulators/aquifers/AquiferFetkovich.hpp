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

#include <opm/simulators/aquifers/AquiferInterface.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <exception>
#include <stdexcept>

namespace Opm
{

template <typename TypeTag>
class AquiferFetkovich : public AquiferInterface<TypeTag>
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

    AquiferFetkovich(const std::vector<Aquancon::AquancCell>& connections,
                     const Simulator& ebosSimulator,
                     const Aquifetp::AQUFETP_data& aqufetp_data)
        : Base(aqufetp_data.aquiferID, connections, ebosSimulator)
        , aqufetp_data_(aqufetp_data)
    {
    }

    void endTimeStep() override
    {
        for (const auto& q : this->Qai_) {
            this->W_flux_ += q * this->ebos_simulator_.timeStepSize();
        }
        aquifer_pressure_ = aquiferPressure();
    }

    data::AquiferData aquiferData() const
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
        data.type = data::AquiferType::Fetkovich;

        data.aquFet = std::make_shared<data::FetkovichData>();

        return data;
    }

protected:
    // Aquifer Fetkovich Specific Variables
    // TODO: using const reference here will cause segmentation fault, which is very strange
    const Aquifetp::AQUFETP_data aqufetp_data_;
    Scalar aquifer_pressure_; // aquifer

    void assignRestartData(const data::AquiferData& xaq) override
    {
        if (xaq.type != data::AquiferType::Fetkovich) {
            throw std::invalid_argument {
                "Analytic aquifer data for unexpected aquifer "
                "type passed to Fetkovich aquifer"
            };
        }

        this->aquifer_pressure_ = xaq.pressure;
    }

    inline Eval dpai(int idx)
    {
        const Eval dp = aquifer_pressure_ - this->pressure_current_.at(idx)
            + this->rhow_[idx] * this->gravity_() * (this->cell_depth_[idx] - this->aquiferDepth());
        return dp;
    }

    // This function implements Eq 5.12 of the EclipseTechnicalDescription
    inline Scalar aquiferPressure()
    {
        Scalar Flux = this->W_flux_.value();
        const auto& comm = this->ebos_simulator_.vanguard().grid().comm();
        comm.sum(&Flux, 1);
        Scalar pa_ = this->pa0_ - Flux / (aqufetp_data_.C_t * aqufetp_data_.V0);
        return pa_;
    }

    inline void calculateAquiferConstants() override
    {
        this->Tc_ = (aqufetp_data_.C_t * aqufetp_data_.V0) / aqufetp_data_.J;
    }
    // This function implements Eq 5.14 of the EclipseTechnicalDescription
    inline void calculateInflowRate(int idx, const Simulator& simulator) override
    {
        const Scalar td_Tc_ = simulator.timeStepSize() / this->Tc_;
        const Scalar coef = (1 - exp(-td_Tc_)) / td_Tc_;
        this->Qai_.at(idx) = this->alphai_[idx] * aqufetp_data_.J * dpai(idx) * coef;
    }

    inline void calculateAquiferCondition() override
    {
        this->rhow_.resize(this->size(), 0.);

        if (this->solution_set_from_restart_) {
            return;
        }

        if (!aqufetp_data_.p0.first) {
            this->pa0_ = this->calculateReservoirEquilibrium();
        } else {
            this->pa0_ = aqufetp_data_.p0.second;
        }
        aquifer_pressure_ = this->pa0_;
    }

    virtual Scalar aquiferDepth() const override
    {
        return aqufetp_data_.d0;
    }
}; // Class AquiferFetkovich
} // namespace Opm
#endif
