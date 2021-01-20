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
#include <stdexcept>

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
    {
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

    Opm::data::AquiferData aquiferData() const
    {
        data::AquiferData data;
        data.aquiferID = this->aquiferID;
        // TODO: not sure how to get this pressure value yet
        data.pressure = this->pa0_;
        data.fluxRate = 0.;
        for (const auto& q : this->Qai_) {
            data.fluxRate += q.value();
        }
        data.volume = this->W_flux_.value();
        data.initPressure = this->pa0_;
        data.type = Opm::data::AquiferType::CarterTracey;
        return data;
    }

protected:
    // Variables constants
    const AquiferCT::AQUCT_data aquct_data_;
    Scalar beta_; // Influx constant
    // TODO: it is possible it should be a AD variable
    Scalar mu_w_{1}; // water viscosity
    Scalar fluxValue_{0}; // value of flux

    void assignRestartData(const data::AquiferData& /* xaq */) override
    {
        throw std::runtime_error {"Restart-based initialization not currently supported "
                                  "for Carter-Tracey analytic aquifers"};
    }

    inline void getInfluenceTableValues(Scalar& pitd, Scalar& pitd_prime, const Scalar& td)
    {
        // We use the opm-common numeric linear interpolator
        pitd = Opm::linearInterpolation(aquct_data_.td, aquct_data_.pi, td);
        pitd_prime = Opm::linearInterpolationDerivative(aquct_data_.td, aquct_data_.pi, td);
    }

    inline Scalar dpai(int idx)
    {
        Scalar dp = this->pa0_
            + this->rhow_.at(idx).value() * this->gravity_() * (this->cell_depth_.at(idx) - this->aquiferDepth())
            - this->pressure_previous_.at(idx);
        return dp;
    }

    // This function implements Eqs 5.8 and 5.9 of the EclipseTechnicalDescription
    inline void calculateEqnConstants(Scalar& a, Scalar& b, const int idx, const Simulator& simulator)
    {
        const Scalar td_plus_dt = (simulator.timeStepSize() + simulator.time()) / this->Tc_;
        const Scalar td = simulator.time() / this->Tc_;
        Scalar PItdprime = 0.;
        Scalar PItd = 0.;
        getInfluenceTableValues(PItd, PItdprime, td_plus_dt);
        a = 1.0 / this->Tc_ * ((beta_ * dpai(idx)) - (this->fluxValue_ * PItdprime)) / (PItd - td * PItdprime);
        b = beta_ / (this->Tc_ * (PItd - td * PItdprime));
    }

    // This function implements Eq 5.7 of the EclipseTechnicalDescription
    inline void calculateInflowRate(int idx, const Simulator& simulator) override
    {
        Scalar a, b;
        calculateEqnConstants(a, b, idx, simulator);
        this->Qai_.at(idx)
            = this->alphai_.at(idx) * (a - b * (this->pressure_current_.at(idx) - this->pressure_previous_.at(idx)));
    }

    inline void calculateAquiferConstants() override
    {
        // We calculate the influx constant
        beta_ = aquct_data_.c2 * aquct_data_.h * aquct_data_.theta * aquct_data_.phi_aq * aquct_data_.C_t
            * aquct_data_.r_o * aquct_data_.r_o;
        // We calculate the time constant
        this->Tc_ = mu_w_ * aquct_data_.phi_aq * aquct_data_.C_t * aquct_data_.r_o * aquct_data_.r_o
            / (aquct_data_.k_a * aquct_data_.c1);
    }

    inline void calculateAquiferCondition() override
    {

        int pvttableIdx = aquct_data_.pvttableID - 1;
        this->rhow_.resize(this->size(), 0.);
        if (!aquct_data_.p0.first) {
            this->pa0_ = this->calculateReservoirEquilibrium();
        } else {
            this->pa0_ = aquct_data_.p0.second;
        }

        // use the thermodynamic state of the first active cell as a
        // reference. there might be better ways to do this...
        ElementContext elemCtx(this->ebos_simulator_);
        auto elemIt = this->ebos_simulator_.gridView().template begin</*codim=*/0>();
        elemCtx.updatePrimaryStencil(*elemIt);
        elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
        const auto& iq0 = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
        // Initialize a FluidState object first
        FluidState fs_aquifer;
        // We use the temperature of the first cell connected to the aquifer
        // Here we copy the fluidstate of the first cell, so we do not accidentally mess up the reservoir fs
        fs_aquifer.assign(iq0.fluidState());
        Eval temperature_aq, pa0_mean, saltConcentration_aq;
        temperature_aq = fs_aquifer.temperature(0);
        saltConcentration_aq = fs_aquifer.saltConcentration();
        pa0_mean = this->pa0_;
        Eval mu_w_aquifer = FluidSystem::waterPvt().viscosity(pvttableIdx, temperature_aq, pa0_mean, saltConcentration_aq);
        mu_w_ = mu_w_aquifer.value();
    }

    virtual Scalar aquiferDepth() const override
    {
        return aquct_data_.d0;
    }
}; // class AquiferCarterTracy
} // namespace Opm

#endif
