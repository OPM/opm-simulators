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
            aquifer_pressure_ = aquiferPressure();
        }
    }

    Opm::data::AquiferData aquiferData() const
    {
        // TODO: how to unify the two functions?
        data::AquiferData data;
        data.aquiferID = this->aquiferID;
        data.pressure = this->aquifer_pressure_;
        data.fluxRate = 0.;
        for (const auto& q : this->Qai_) {
            data.fluxRate += q.value();
        }
        data.volume = this->W_flux_.value();
        data.initPressure = this->pa0_;
        data.type = Opm::data::AquiferType::Fetkovich;
        // Not handling std::shared_ptr<FetkovichData> aquFet for now,
        // because we do not need it yet
        return data;
    }

protected:
    // Aquifer Fetkovich Specific Variables
    // TODO: using const reference here will cause segmentation fault, which is very strange
    const Aquifetp::AQUFETP_data aqufetp_data_;
    Scalar aquifer_pressure_; // aquifer

    inline void initializeConnections() override
    {
        this->cell_depth_.resize(this->size(), this->aquiferDepth());
        this->alphai_.resize(this->size(), 1.0);
        this->faceArea_connected_.resize(this->size(), 0.0);

        // Translate the C face tag into the enum used by opm-parser's TransMult class
        Opm::FaceDir::DirEnum faceDirection;

        // denom_face_areas is the sum of the areas connected to an aquifer
        Scalar denom_face_areas = 0.;
        this->cellToConnectionIdx_.resize(this->ebos_simulator_.gridView().size(/*codim=*/0), -1);
        for (size_t idx = 0; idx < this->size(); ++idx) {
            const auto global_index = this->connections_[idx].global_index;
            const int cell_index = this->ebos_simulator_.vanguard().compressedIndex(global_index);
            if (cell_index < 0) //the global_index is not part of this grid
                continue;

            this->cellToConnectionIdx_[cell_index] = idx;
            this->cell_depth_.at(idx) = this->ebos_simulator_.vanguard().cellCenterDepth(cell_index);
        }
        // get areas for all connections
        const auto& gridView = this->ebos_simulator_.vanguard().gridView();
        ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());
        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            unsigned cell_index = elemMapper.index(elem);
            int idx = this->cellToConnectionIdx_[cell_index];

            // only deal with connections given by the aquifer
            if( idx < 0)
                continue;

            if (!this->connections_[idx].influx_coeff.first) { // influx_coeff is defaulted
                auto isIt = gridView.ibegin(elem);
                const auto& isEndIt = gridView.iend(elem);
                for (; isIt != isEndIt; ++ isIt) {
                    // store intersection, this might be costly
                    const auto& intersection = *isIt;

                    // only deal with grid boundaries
                    if (!intersection.boundary())
                        continue;

                    int insideFaceIdx  = intersection.indexInInside();
                    switch (insideFaceIdx) {
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
                            "Initialization of Aquifer Fetkovich problem. Make sure faceTag is correctly defined");                }


                    if (faceDirection == this->connections_[idx].face_dir) {
                        this->faceArea_connected_[idx] = this->getFaceArea(intersection, idx);
                        break;
                    }
                }
            } else {
                this->faceArea_connected_.at(idx) = this->connections_[idx].influx_coeff.second;
            }
            denom_face_areas += (this->connections_[idx].influx_mult * this->faceArea_connected_.at(idx));
        }

        const double eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());
        for (size_t idx = 0; idx < this->size(); ++idx) {
            this->alphai_.at(idx) = (denom_face_areas < eps_sqrt)
                ? // Prevent no connection NaNs due to division by zero
                0.
                : (this->connections_[idx].influx_mult * this->faceArea_connected_.at(idx)) / denom_face_areas;
        }
    }

    void assignRestartData(const data::AquiferData& xaq) override
    {
        if (xaq.type != data::AquiferType::Fetkovich) {
            throw std::invalid_argument {"Analytic aquifer data for unexpected aquifer type "
                                         "passed to Fetkovich aquifer"};
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
