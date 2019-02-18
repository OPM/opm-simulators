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

#include <opm/autodiff/AquiferInterface.hpp>

namespace Opm
{

    template<typename TypeTag>
    class AquiferCarterTracy: public AquiferInterface<TypeTag>
    {
        public:
          typedef AquiferInterface<TypeTag> Base;

          using typename Base::Simulator;
          using typename Base::ElementContext;
          using typename Base::FluidSystem;
          using typename Base::BlackoilIndices;
          using typename Base::RateVector;
          using typename Base::IntensiveQuantities;
          using typename Base::Eval;
          using typename Base::Scalar;
          using typename Base::FluidState;

          using Base::waterCompIdx;
          using Base::waterPhaseIdx;
            AquiferCarterTracy( const Aquancon::AquanconOutput& connection,
                                const std::unordered_map<int, int>& cartesian_to_compressed,
                                const Simulator& ebosSimulator,
                                const AquiferCT::AQUCT_data& aquct_data)
            : Base(connection, cartesian_to_compressed, ebosSimulator)
            , aquct_data_(aquct_data)
            {}

            void endTimeStep()
            {
                for (const auto& Qai: Base::Qai_) {
                    Base::W_flux_ += Qai*Base::ebos_simulator_.timeStepSize();
                }
            }

        protected:
            // Variables constants
            const AquiferCT::AQUCT_data aquct_data_;
            Scalar beta_; // Influx constant

            // This function is used to initialize and calculate the alpha_i for each grid connection to the aquifer
            inline void initializeConnections(const Aquancon::AquanconOutput& connection)
            {
                const auto& eclState = Base::ebos_simulator_.vanguard().eclState();
                const auto& ugrid = Base::ebos_simulator_.vanguard().grid();
                const auto& grid = eclState.getInputGrid();

                Base::cell_idx_ = connection.global_index;
                auto globalCellIdx = ugrid.globalCell();

                assert( Base::cell_idx_ == connection.global_index);
                assert( (Base::cell_idx_.size() <= connection.influx_coeff.size()) );
                assert( (connection.influx_coeff.size() == connection.influx_multiplier.size()) );
                assert( (connection.influx_multiplier.size() == connection.reservoir_face_dir.size()) );

                // We hack the cell depth values for now. We can actually get it from elementcontext pos
                Base::cell_depth_.resize(Base::cell_idx_.size(), aquct_data_.d0);
                Base::alphai_.resize(Base::cell_idx_.size(), 1.0);
                Base::faceArea_connected_.resize(Base::cell_idx_.size(),0.0);

                auto cell2Faces = Opm::UgGridHelpers::cell2Faces(ugrid);
                auto faceCells  = Opm::UgGridHelpers::faceCells(ugrid);

                // Translate the C face tag into the enum used by opm-parser's TransMult class
                Opm::FaceDir::DirEnum faceDirection;

                // denom_face_areas is the sum of the areas connected to an aquifer
                Scalar denom_face_areas = 0.;
                Base::cellToConnectionIdx_.resize(Base::ebos_simulator_.gridView().size(/*codim=*/0), -1);
                for (size_t idx = 0; idx < Base::cell_idx_.size(); ++idx)
                {
                    const int cell_index = Base::cartesian_to_compressed_.at(Base::cell_idx_[idx]);
                    Base::cellToConnectionIdx_[cell_index] = idx;

                    const auto cellFacesRange = cell2Faces[cell_index];
                    for(auto cellFaceIter = cellFacesRange.begin(); cellFaceIter != cellFacesRange.end(); ++cellFaceIter)
                    {
                        // The index of the face in the compressed grid
                        const int faceIdx = *cellFaceIter;

                        // the logically-Cartesian direction of the face
                        const int faceTag = Opm::UgGridHelpers::faceTag(ugrid, cellFaceIter);

                        switch(faceTag)
                        {
                            case 0: faceDirection = Opm::FaceDir::XMinus;
                                    break;
                            case 1: faceDirection = Opm::FaceDir::XPlus;
                                    break;
                            case 2: faceDirection = Opm::FaceDir::YMinus;
                                    break;
                            case 3: faceDirection = Opm::FaceDir::YPlus;
                                    break;
                            case 4: faceDirection = Opm::FaceDir::ZMinus;
                                    break;
                            case 5: faceDirection = Opm::FaceDir::ZPlus;
                                    break;
                            default: OPM_THROW(Opm::NumericalIssue,"Initialization of Aquifer Carter Tracy problem. Make sure faceTag is correctly defined");
                        }

                        if (faceDirection == connection.reservoir_face_dir.at(idx))
                        {
                            Base::faceArea_connected_.at(idx) = Base::getFaceArea(faceCells, ugrid, faceIdx, idx, connection);
                            denom_face_areas += ( connection.influx_multiplier.at(idx) * Base::faceArea_connected_.at(idx) );
                        }
                    }
                    auto cellCenter = grid.getCellCenter(Base::cell_idx_.at(idx));
                    Base::cell_depth_.at(idx) = cellCenter[2];
                }

                const double eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());
                for (size_t idx = 0; idx < Base::cell_idx_.size(); ++idx)
                {
                    Base::alphai_.at(idx) = (denom_face_areas < eps_sqrt)? // Prevent no connection NaNs due to division by zero
                    0.
                    : ( connection.influx_multiplier.at(idx) * Base::faceArea_connected_.at(idx) )/denom_face_areas;
                }
            }

            inline void getInfluenceTableValues(Scalar& pitd, Scalar& pitd_prime, const Scalar& td)
            {
                // We use the opm-common numeric linear interpolator
                pitd = Opm::linearInterpolation(aquct_data_.td, aquct_data_.pi, td);
                pitd_prime = Opm::linearInterpolationDerivative(aquct_data_.td, aquct_data_.pi, td);
            }

            inline Scalar dpai(int idx)
            {
                Scalar dp = Base::pa0_ + Base::rhow_.at(idx).value()*Base::gravity_()*(Base::cell_depth_.at(idx) - aquct_data_.d0) - Base::pressure_previous_.at(idx);
                return dp;
            }

            // This function implements Eqs 5.8 and 5.9 of the EclipseTechnicalDescription
            inline void calculateEqnConstants(Scalar& a, Scalar& b, const int idx, const Simulator& simulator)
            {
                const Scalar td_plus_dt = (simulator.timeStepSize() + simulator.time()) / Base::Tc_;
                const Scalar td = simulator.time() / Base::Tc_;
                Scalar PItdprime = 0.;
                Scalar PItd = 0.;
                getInfluenceTableValues(PItd, PItdprime, td_plus_dt);
                a = 1.0/Base::Tc_ * ( (beta_ * dpai(idx)) - (Base::W_flux_.value() * PItdprime) ) / ( PItd - td*PItdprime );
                b = beta_ / (Base::Tc_ * ( PItd - td*PItdprime));
            }

            // This function implements Eq 5.7 of the EclipseTechnicalDescription
            inline void calculateInflowRate(int idx, const Simulator& simulator)
            {
                Scalar a, b;
                calculateEqnConstants(a,b,idx,simulator);
                Base::Qai_.at(idx) = Base::alphai_.at(idx)*( a - b * ( Base::pressure_current_.at(idx) - Base::pressure_previous_.at(idx) ) );
            }

            inline void calculateAquiferConstants()
            {
                // We calculate the influx constant
                beta_ = aquct_data_.c2 * aquct_data_.h
                        * aquct_data_.theta * aquct_data_.phi_aq
                        * aquct_data_.C_t
                        * aquct_data_.r_o * aquct_data_.r_o;
                // We calculate the time constant
                Base::Tc_ = Base::mu_w_ * aquct_data_.phi_aq
                      * aquct_data_.C_t
                      * aquct_data_.r_o * aquct_data_.r_o
                      / ( aquct_data_.k_a * aquct_data_.c1 );
            }

            inline void calculateAquiferCondition()
            {

                int pvttableIdx = aquct_data_.pvttableID - 1;
                Base::rhow_.resize(Base::cell_idx_.size(),0.);
                if (!aquct_data_.p0)
                {
                   Base::pa0_ = calculateReservoirEquilibrium();
                }
                else
                {
                   Base::pa0_ = *(aquct_data_.p0);
                }

                // use the thermodynamic state of the first active cell as a
                // reference. there might be better ways to do this...
                ElementContext elemCtx(Base::ebos_simulator_);
                auto elemIt = Base::ebos_simulator_.gridView().template begin</*codim=*/0>();
                elemCtx.updatePrimaryStencil(*elemIt);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const auto& iq0 = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                // Initialize a FluidState object first
                FluidState fs_aquifer;
                // We use the temperature of the first cell connected to the aquifer
                // Here we copy the fluidstate of the first cell, so we do not accidentally mess up the reservoir fs
                fs_aquifer.assign( iq0.fluidState() );
                Eval temperature_aq, pa0_mean;
                temperature_aq = fs_aquifer.temperature(0);
                pa0_mean = Base::pa0_;
                Eval mu_w_aquifer = FluidSystem::waterPvt().viscosity(pvttableIdx, temperature_aq, pa0_mean);
                Base::mu_w_ = mu_w_aquifer.value();

            }

            // This function is for calculating the aquifer properties from equilibrium state with the reservoir
            inline Scalar calculateReservoirEquilibrium()
            {
                // Since the global_indices are the reservoir index, we just need to extract the fluidstate at those indices
                std::vector<Scalar> pw_aquifer;
                Scalar water_pressure_reservoir;

                ElementContext elemCtx(Base::ebos_simulator_);
                const auto& gridView = Base::ebos_simulator_.gridView();
                auto elemIt = gridView.template begin</*codim=*/0>();
                const auto& elemEndIt = gridView.template end</*codim=*/0>();
                for (; elemIt != elemEndIt; ++elemIt) {
                    const auto& elem = *elemIt;
                    elemCtx.updatePrimaryStencil(elem);

                    size_t cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                    int idx = Base::cellToConnectionIdx_[cellIdx];
                    if (idx < 0)
                    continue;

                    elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                    const auto& iq0 = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                    const auto& fs = iq0.fluidState();

                    water_pressure_reservoir = fs.pressure(waterPhaseIdx).value();
                    Base::rhow_[idx] = fs.density(waterPhaseIdx);
                    pw_aquifer.push_back( (water_pressure_reservoir - Base::rhow_[idx].value()*Base::gravity_()*(Base::cell_depth_[idx] - aquct_data_.d0))*Base::alphai_[idx] );
                }

                // We take the average of the calculated equilibrium pressures.
                Scalar aquifer_pres_avg = std::accumulate(pw_aquifer.begin(), pw_aquifer.end(), 0.)/pw_aquifer.size();
                return aquifer_pres_avg;
            }
    }; // class AquiferCarterTracy
} // namespace Opm

#endif
