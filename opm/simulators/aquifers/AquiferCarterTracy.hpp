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

    using Base::waterCompIdx;
    using Base::waterPhaseIdx;
    AquiferCarterTracy(const std::vector<Aquancon::AquancCell>& connections,
                       const std::unordered_map<int, int>& cartesian_to_compressed,
                       const Simulator& ebosSimulator,
                       const AquiferCT::AQUCT_data& aquct_data)
        : Base(aquct_data.aquiferID, connections, cartesian_to_compressed, ebosSimulator)
        , aquct_data_(aquct_data)
    {
    }

    void endTimeStep() override
    {
        for (const auto& q : this->Qai_) {
            this->W_flux_ += q * this->ebos_simulator_.timeStepSize();
        }
    }

protected:
    // Variables constants
    const AquiferCT::AQUCT_data aquct_data_;
    Scalar beta_; // Influx constant
    // TODO: it is possible it should be a AD variable
    Scalar mu_w_; // water viscosity

    // This function is used to initialize and calculate the alpha_i for each grid connection to the aquifer
    inline void initializeConnections() override
    {
        const auto& eclState = this->ebos_simulator_.vanguard().eclState();
        const auto& ugrid = this->ebos_simulator_.vanguard().grid();
        const auto& grid = eclState.getInputGrid();

        auto globalCellIdx = ugrid.globalCell();

        // We hack the cell depth values for now. We can actually get it from elementcontext pos
        this->cell_depth_.resize(this->size(), this->aquiferDepth());
        this->alphai_.resize(this->size(), 1.0);
        this->faceArea_connected_.resize(this->size(), 0.0);

        auto cell2Faces = Opm::UgGridHelpers::cell2Faces(ugrid);
        auto faceCells = Opm::UgGridHelpers::faceCells(ugrid);

        // Translate the C face tag into the enum used by opm-parser's TransMult class
        Opm::FaceDir::DirEnum faceDirection;

        // denom_face_areas is the sum of the areas connected to an aquifer
        Scalar denom_face_areas = 0.;
        this->cellToConnectionIdx_.resize(this->ebos_simulator_.gridView().size(/*codim=*/0), -1);
        for (size_t idx = 0; idx < this->size(); ++idx) {
            const int cell_index = this->cartesian_to_compressed_.at(this->connections_[idx].global_index);
            this->cellToConnectionIdx_[cell_index] = idx;

            const auto cellFacesRange = cell2Faces[cell_index];
            for (auto cellFaceIter = cellFacesRange.begin(); cellFaceIter != cellFacesRange.end(); ++cellFaceIter) {
                // The index of the face in the compressed grid
                const int faceIdx = *cellFaceIter;

                // the logically-Cartesian direction of the face
                const int faceTag = Opm::UgGridHelpers::faceTag(ugrid, cellFaceIter);

                switch (faceTag) {
                case 0:
                    faceDirection = Opm::FaceDir::XMinus;
                    break;
                case 1:
                    faceDirection = Opm::FaceDir::XPlus;
                    break;
                case 2:
                    faceDirection = Opm::FaceDir::YMinus;
                    break;
                case 3:
                    faceDirection = Opm::FaceDir::YPlus;
                    break;
                case 4:
                    faceDirection = Opm::FaceDir::ZMinus;
                    break;
                case 5:
                    faceDirection = Opm::FaceDir::ZPlus;
                    break;
                default:
                    OPM_THROW(Opm::NumericalIssue,
                              "Initialization of Aquifer Carter Tracy problem. Make sure faceTag is correctly defined");
                }

                if (faceDirection == this->connections_[idx].face_dir) {
                    this->faceArea_connected_.at(idx) = this->getFaceArea(faceCells, ugrid, faceIdx, idx);
                    denom_face_areas += (this->connections_[idx].influx_mult * this->faceArea_connected_.at(idx));
                }
            }
            auto cellCenter = grid.getCellCenter(this->connections_[idx].global_index);
            this->cell_depth_.at(idx) = cellCenter[2];
        }

        const double eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());
        for (size_t idx = 0; idx < this->size(); ++idx) {
            this->alphai_.at(idx) = (denom_face_areas < eps_sqrt)
                ? // Prevent no connection NaNs due to division by zero
                0.
                : (this->connections_[idx].influx_mult * this->faceArea_connected_.at(idx)) / denom_face_areas;
        }
    }

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
        a = 1.0 / this->Tc_ * ((beta_ * dpai(idx)) - (this->W_flux_.value() * PItdprime)) / (PItd - td * PItdprime);
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
