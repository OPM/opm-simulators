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

#include <Eigen/QR>
#include <opm/parser/eclipse/EclipseState/AquiferCT.hpp>
#include <opm/autodiff/BlackoilAquiferModel.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/core/props/BlackoilPhases.hpp>


#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <cassert>

namespace Opm
{

    template<typename TypeTag>
    class AquiferCarterTracy
    {

        public:
            typedef BlackoilModelParameters ModelParameters;

            static const int Water = BlackoilPhases::Aqua;
            static const int Oil = BlackoilPhases::Liquid;
            static const int Gas = BlackoilPhases::Vapour;

            typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
            typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
            typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
            typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
            typedef typename GET_PROP_TYPE(TypeTag, Indices) BlackoilIndices;
            typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
            typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
            typedef typename GridView::template Codim<0>::Entity Element;
            typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

            static const int numEq = BlackoilIndices::numEq;
            typedef double Scalar;

            typedef DenseAd::Evaluation<double, /*size=*/numEq> Eval;

            typedef Ewoms::BlackOilPolymerModule<TypeTag> PolymerModule;

            static const bool has_solvent = GET_PROP_VALUE(TypeTag, EnableSolvent);
            static const bool has_polymer = GET_PROP_VALUE(TypeTag, EnablePolymer);
            static const int contiSolventEqIdx = BlackoilIndices::contiSolventEqIdx;
            static const int contiPolymerEqIdx = BlackoilIndices::contiPolymerEqIdx;


            AquiferCarterTracy(const std::vector<int>& cell_id)
            : phi_aq_ (1.0), //
              C_t_ (1.0), //
              r_o_ (1.0), //
              k_a_ (1.0), //
              c1_ (1.0),
              h_ (1.0), //
              theta_ (1.0), //
              c2_ (1.0), //
              d0_ (1.0),
              cell_idx_ (cell_id)
            {
                mu_w_ = 1e-3;
                aqutab_td_.push_back(1.0);
                aqutab_pi_.push_back(1.0);
                aquiferID_ = 1;
                inftableID_ = 1;
                pvttableID_ = 1;
                init_quantities();
            }

            explicit AquiferCarterTracy( const AquiferCT::AQUCT_data& params, const AquiferCT::AQUANCON_data& aquanconParams, 
                                         const int numComponents, const Scalar gravity                                        )
            : phi_aq_ (params.phi_aq), //
              C_t_ (params.C_t), //
              r_o_ (params.r_o), //
              k_a_ (params.k_a), //
              c1_ (params.c1),
              h_ (params.h), //
              theta_ (params.theta), //
              c2_ (params.c2), //
              d0_ (params.d0),
              aqutab_td_ (params.td),
              aqutab_pi_ (params.pi),
              aquiferID_ (params.aquiferID),
              inftableID_ (params.inftableID),
              pvttableID_ (params.pvttableID),
              cell_idx_ (params.cell_id),
              num_components_ (numComponents),
              gravity_ (gravity)
            {
                mu_w_ = 1e-3;
                init_quantities(aquanconParams);
            }

            inline const PhaseUsage&
            phaseUsage() const
            {
                assert(phase_usage_);

                return *phase_usage_;
            }

            inline int
            flowPhaseToEbosCompIdx( const int phaseIdx ) const
            {
                const auto& pu = phaseUsage();
                if (pu.phase_pos[Water] == phaseIdx)
                    return BlackoilIndices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                if (pu.phase_pos[Oil] == phaseIdx)
                    return BlackoilIndices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                if (pu.phase_pos[Gas] == phaseIdx)
                    return BlackoilIndices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);

                // for other phases return the index
                return phaseIdx;
            }

            inline int
            flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
            {
                const auto& pu = phaseUsage();
                if (pu.phase_pos[Water] == phaseIdx) {
                    return FluidSystem::waterPhaseIdx;
                }
                if (pu.phase_pos[Oil] == phaseIdx) {
                    return FluidSystem::oilPhaseIdx;
                }
                if (pu.phase_pos[Gas] == phaseIdx) {
                    return FluidSystem::gasPhaseIdx;
                }

                assert(phaseIdx < 3);
                // for other phases return the index
                return phaseIdx;
            }

            inline void calculateExplicitQuantities(const Simulator& ebosSimulator)
            {
                std::cout << "In CarterTracy<calculateExplicitQuantities>: I am aquifer #" << aquiferID_ << std::endl;
            }

            inline void assembleAquiferEq(Simulator& ebosSimulator, const SimulatorTimerInterface& timer)
            {
                std::cout << "In CarterTracy<assembleAquiferEq>: I am aquifer #" << aquiferID_ << std::endl;
                // resAqui_ = 0.0;
                dt_ = timer.currentStepLength();
                auto& ebosJac = ebosSimulator.model().linearizer().matrix();
                auto& ebosResid = ebosSimulator.model().linearizer().residual();

                // TODO: it probably can be static member for StandardWell
                const double volume = 0.002831684659200; // 0.1 cu ft;

                auto cellID = cell_idx_.begin();
                size_t idx;
                for ( idx = 0; cellID != cell_idx_.end(); ++cellID, ++idx )
                {
                    Eval qinflow = 0.0;
                    // We are dereferencing the value of IntensiveQuantities because cachedIntensiveQuantities return a const pointer to
                    // IntensiveQuantities of that particular cell_id
                    const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(*cellID, /*timeIdx=*/ 0));
                    // This is the pressure at td + dt
                    get_current_Pressure_cell(pressure_current_,idx,intQuants);
                    get_current_density_cell(rhow_,idx,intQuants);
                    calculate_inflow_rate(idx, timer);
                    qinflow = Qai_[idx];
                    ebosResid[*cellID][flowPhaseToEbosCompIdx(FluidSystem::waterPhaseIdx)] -= qinflow.value();

                    for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) 
                    {
                        // also need to consider the efficiency factor when manipulating the jacobians.
                        ebosJac[*cellID][*cellID][flowPhaseToEbosCompIdx(FluidSystem::waterPhaseIdx)][pvIdx] -= qinflow.derivative(pvIdx);
                    }
                    std::cout << "In CarterTracy<assembleAquiferEq>: I am aquifer #" << aquiferID_
                              // << " -> P_wat[t+dt] = " << pressure_current_[idx] << std::endl
                              << " Qai[t+dt] = " << Qai_[idx] << std::endl;
                }
            }

            inline void before_time_step(Simulator& ebosSimulator, const SimulatorTimerInterface& timer)
            {
                auto cellID = cell_idx_.begin();
                size_t idx;
                for ( idx = 0; cellID != cell_idx_.end(); ++cellID, ++idx )
                {
                    const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(*cellID, /*timeIdx=*/ 0));
                    get_current_Pressure_cell(pressure_previous_ ,idx,intQuants);
                }
            }

            inline void after_time_step()
            {
                for (auto Qai = Qai_.begin(); Qai != Qai_.end(); ++Qai)
                {
                    W_flux_ += (*Qai);
                }
                std::cout << "Aquifer # " << aquiferID_ << ": My cumulative flux = " << W_flux_ << std::endl;
            }

            /* Made into public for testing only!!!!!!. Must be protected */
            inline const Scalar time_constant() const
            {
                Scalar Tc = mu_w_*phi_aq_*C_t_*r_o_*r_o_/(k_a_*c1_);
                return Tc;
            }

            /* Made into public for testing only!!!!!!. Must be protected */
            inline const Scalar aquifer_influx_constant() const
            {
                Scalar beta = c2_*h_*theta_*phi_aq_*C_t_*r_o_*r_o_;
                return beta;
            }

            // This is another hack to get the face area only for SPE1.
            // Ideally it should be a map which given a cell_id, it returns the area fraction
            inline const double area_fraction(const int i)
            {
                return 1000.0*20.0*0.092903/(1000.0*1000.0*0.092903*2 + 1000.0*20.0*0.092903*4);
            }

            inline void print_private_members() const
            {
                std::cout << "Aquifer CT #" << aquiferID_ << std::endl;
                auto ita = aqutab_td_.cbegin();
                auto f_lambda = [&ita] (double i) {std::cout << *ita++ << "    " << i << std::endl;};
                std::for_each( aqutab_pi_.cbegin(), aqutab_pi_.cend(), f_lambda );

                for (auto i = coeff_.begin(); i != coeff_.end(); ++i )
                {
                    std::cout << "Coeff = " << *i << std::endl;
                }
            }

            /* Made into public for testing only!!!!!!. Must be protected */
            inline const std::vector<int> cell_id() const
            {
                return cell_idx_;
            }

            inline const int& aquiferID() const
            {
                return aquiferID_;
            }



        protected:
            const PhaseUsage* phase_usage_;


            // Aquifer ID, and other IDs
            int aquiferID_, inftableID_, pvttableID_;
            int num_components_;

            // Grid variables
            std::vector<int> cell_idx_;

            // Quantities at each grid id
            std::vector<Scalar> cell_depth_;
            std::vector<Scalar> pressure_previous_;
            std::vector<Scalar> pressure_current_;
            std::vector<Scalar> Qai_;
            std::vector<Scalar> rhow_;

            // Variables constants
            Scalar mu_w_ , //water viscosity
                   phi_aq_ , //aquifer porosity
                   d0_,     // aquifer datum depth
                   C_t_ , //total compressibility
                   r_o_ , //aquifer inner radius
                   k_a_ , //aquifer permeability
                   c1_, // 0.008527 (METRIC, PVT-M); 0.006328 (FIELD); 3.6 (LAB)
                   h_ , //aquifer thickness
                   theta_ , //angle subtended by the aquifer boundary
                   c2_ ; //6.283 (METRIC, PVT-M); 1.1191 (FIELD); 6.283 (LAB).

            // Variables for influence table
            std::vector<Scalar> aqutab_td_, aqutab_pi_;

            // Cumulative flux
            Scalar W_flux_, dt_, pa0_, gravity_;

            // Also return the polynomial fit
            std::vector<Scalar> coeff_;
            

            // We fit the tabular data using a polynomial fit
            // Modified from Copyright (C) 2014  Clifford Wolf <clifford@clifford.at> 
            // http://svn.clifford.at/handicraft/2014/polyfit/polyfit.cc
            inline void polynomial_fit( const std::vector<Scalar> &X, const std::vector<Scalar> &y, 
                                        std::vector<Scalar> &coeff, int order, bool bias) const
            {
                int colNum = (bias)? order + 1 : order;
                Eigen::MatrixXd A(X.size(), colNum);
                Eigen::VectorXd y_mapped = Eigen::VectorXd::Map(&y.front(), y.size());
                Eigen::VectorXd result;

                assert(X.size() == y.size());
                assert(X.size() >= colNum);

                // create matrix
                for (size_t i = 0; i < X.size(); i++)
                for (size_t j = 0; j < colNum; j++)
                    A(i, j) = (bias)? pow(X.at(i), j) : pow(X.at(i), j+1);

                // solve for linear least squares fit
                result = A.householderQr().solve(y_mapped);

                coeff.resize(colNum);
                for (size_t i = 0; i < colNum; i++)
                    coeff[i] = result[i];
            }

            inline void init_quantities(const AquiferCT::AQUANCON_data& aquanconParams)
            {
                W_flux_ = 0.;
                // pa0_ is the initial aquifer water pressure. Must be calculated from equilibrium if left default,
                // or we get the information from the deck Hacked to make it at 45e6 Pa
                pa0_ = 45e6;

                pressure_previous_.resize(cell_idx_.size(), 0.);
                pressure_current_.resize(cell_idx_.size(), 0.);
                // We hack the cell depth values for now. We can actually get it from elementcontext pos
                cell_depth_.resize(cell_idx_.size(), d0_); 
                rhow_.resize(cell_idx_.size(), 998.0); 
                Qai_.resize(cell_idx_.size(), 0.);

                polynomial_fit(aqutab_td_, aqutab_pi_, coeff_, 2, true);
            }

            inline void get_current_Pressure_cell(std::vector<Scalar>& pressure_water, const int idx, const IntensiveQuantities& intQuants)
            {
                const auto& fs = intQuants.fluidState();
                pressure_water[idx] = fs.pressure(FluidSystem::waterPhaseIdx).value();
            }

            inline void get_current_density_cell(std::vector<Scalar>& rho_water, const int idx, const IntensiveQuantities& intQuants)
            {
                const auto& fs = intQuants.fluidState();
                rho_water[idx] = fs.density(FluidSystem::waterPhaseIdx).value();
            }

            inline Scalar dpai(int idx)
            {
                Scalar dp = pa0_ - rhow_[idx]*gravity_*(cell_depth_[idx] - d0_) - pressure_previous_[idx];
                return dp;
            }

            inline void calculate_a_b_constants(Scalar& a, Scalar& b, const int idx, const SimulatorTimerInterface& timer)
            {
                // This function implements Eqs 5.8 and 5.9 of the EclipseTechnicalDescription
                Scalar beta = aquifer_influx_constant();
                Scalar Tc = time_constant();
                Scalar td_plus_dt = (timer.currentStepLength() + timer.simulationTimeElapsed()) / Tc;
                Scalar td = timer.simulationTimeElapsed() / Tc;
                Scalar PItdprime = coeff_[1] + 2.0*coeff_[2]*(td_plus_dt);
                Scalar PItd = coeff_[0] + coeff_[1]*td_plus_dt + coeff_[2]*td_plus_dt*td_plus_dt;
                a = 1.0/Tc * ( (beta * dpai(idx)) - (W_flux_ * PItdprime) ) / ( PItd - td*PItdprime );
                b = beta / Tc / ( PItd - td*PItdprime);
            }

            inline void calculate_inflow_rate(int idx, const SimulatorTimerInterface& timer)
            {
                Scalar a, b;
                calculate_a_b_constants(a,b,idx,timer);
                // This function implements Eq 5.7 of the EclipseTechnicalDescription
                Qai_[idx] = area_fraction(idx)*( a - b * ( pressure_current_[idx] - pressure_previous_[idx] ) );
            }
            

    }; // class AquiferCarterTracy


} // namespace Opm

#endif