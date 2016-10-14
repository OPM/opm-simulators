/*
  Copyright 2015, 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil AS.

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

#ifndef OPM_BLACKOILPRESSUREMODEL_HEADER_INCLUDED
#define OPM_BLACKOILPRESSUREMODEL_HEADER_INCLUDED


#include <opm/autodiff/BlackoilModelBase.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/core/simulator/SimulatorTimerInterface.hpp>

#include <algorithm>

namespace Opm {

    /// A model implementation for the pressure equation in three-phase black oil.
    ///
    /// The model is based on the normal black oil model.
    /// It uses automatic differentiation via the class AutoDiffBlock
    /// to simplify assembly of the jacobian matrix.
    template<class Grid, class WellModel>
    class BlackoilPressureModel : public BlackoilModelBase<Grid, WellModel, BlackoilPressureModel<Grid, WellModel> >
    {
    public:

        typedef BlackoilModelBase<Grid, WellModel, BlackoilPressureModel<Grid, WellModel> > Base;
        friend Base;

        typedef typename Base::ReservoirState ReservoirState;
        typedef typename Base::WellState WellState;
        typedef typename Base::SolutionState SolutionState;
        typedef typename Base::V V;

        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] fluid            fluid properties
        /// \param[in] geo              rock properties
        /// \param[in] rock_comp_props  if non-null, rock compressibility properties
        /// \param[in] wells_arg        well structure
        /// \param[in] linsolver        linear solver
        /// \param[in] eclState         eclipse state
        /// \param[in] has_disgas       turn on dissolved gas
        /// \param[in] has_vapoil       turn on vaporized oil feature
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilPressureModel(const typename Base::ModelParameters&   param,
                              const Grid&                             grid,
                              const BlackoilPropsAdInterface&         fluid,
                              const DerivedGeology&                   geo,
                              const RockCompressibility*              rock_comp_props,
                              const StandardWells&                    std_wells,
                              const NewtonIterationBlackoilInterface& linsolver,
                              std::shared_ptr< const EclipseState>    eclState,
                              const bool                              has_disgas,
                              const bool                              has_vapoil,
                              const bool                              terminal_output)
            : Base(param, grid, fluid, geo, rock_comp_props, std_wells, linsolver,
                   eclState, has_disgas, has_vapoil, terminal_output),
              state0_(3),
              max_dp_rel_(std::numeric_limits<double>::infinity()),
              scaling_{ ADB::null(), ADB::null(), ADB::null() }
        {
        }

        /// Called once per timestep.
        void prepareStep(const SimulatorTimerInterface& timer,
                         const ReservoirState& reservoir_state,
                         const WellState& well_state)
        {
            asImpl().wellModel().setStoreWellPerforationFluxesFlag(true);
            Base::prepareStep(timer, reservoir_state, well_state);
            max_dp_rel_ = std::numeric_limits<double>::infinity();
            state0_ = asImpl().variableState(reservoir_state, well_state);
            asImpl().makeConstantState(state0_);
        }


        /// Solve the Jacobian system Jx = r where J is the Jacobian and
        /// r is the residual.
        V solveJacobianSystem() const
        {
            // We make a reduced residual object which just
            // contains the pressure residual.
            // TODO: directly contruct that object in residual_ instead.
            const int n1 = residual_.material_balance_eq[0].size();
            const int n2 = residual_.well_flux_eq.size() + residual_.well_eq.size();
            const int n_full = residual_.sizeNonLinear();
            const auto& mb = residual_.material_balance_eq;
            const auto& fe = residual_.well_flux_eq;
            const auto& we = residual_.well_eq;
            LinearisedBlackoilResidual pressure_res = {
                {
                    // TODO: handle general 2-phase etc.
                    ADB::function(mb[0].value(), { mb[0].derivative()[0], mb[0].derivative()[3], mb[0].derivative()[4] })
                },
                ADB::function(fe.value(), { fe.derivative()[0], fe.derivative()[3], fe.derivative()[4] }),
                ADB::function(we.value(), { we.derivative()[0], we.derivative()[3], we.derivative()[4] }),
                residual_.matbalscale,
                residual_.singlePrecision
            };
            assert(pressure_res.sizeNonLinear() == n1 + n2);
            V dx_pressure = linsolver_.computeNewtonIncrement(pressure_res);
            assert(dx_pressure.size() == n1 + n2);
            V dx_full = V::Zero(n_full);
            dx_full.topRows(n1) = dx_pressure.topRows(n1);
            dx_full.bottomRows(n2) = dx_pressure.bottomRows(n2);
            return dx_full;
        }

        using Base::numPhases;
        using Base::numMaterials;
        using Base::wellModel;

    protected:
        using Base::asImpl;
        using Base::linsolver_;
        using Base::residual_;
        using Base::sd_;
        using Base::grid_;
        using Base::ops_;
        using Base::has_vapoil_;
        using Base::has_disgas_;

        SolutionState state0_;
        double max_dp_rel_ = std::numeric_limits<double>::infinity();
        ADB scaling_[3] = { ADB::null(), ADB::null(), ADB::null() };




        IterationReport
        assemble(const ReservoirState& reservoir_state,
                 WellState& well_state,
                 const bool initial_assembly)
        {
            IterationReport iter_report = Base::assemble(reservoir_state, well_state, initial_assembly);

            if (initial_assembly) {
            }

            // Compute pressure residual.
            ADB pressure_residual = ADB::constant(V::Zero(residual_.material_balance_eq[0].size()));
            for (int phase = 0; phase < numPhases(); ++phase) {
                pressure_residual += residual_.material_balance_eq[phase] * scaling_[phase];
            }
            residual_.material_balance_eq[0] = pressure_residual; // HACK

            // Compute total reservoir volume flux.
            const int n = sd_.rq[0].mflux.size();
            V flux = V::Zero(n);
            for (int phase = 0; phase < numPhases(); ++phase) {
                UpwindSelector<double> upwind(grid_, ops_, sd_.rq[phase].dh.value());
                flux += sd_.rq[phase].mflux.value() / upwind.select(sd_.rq[phase].b.value());
            }

            // Storing the fluxes in the assemble() method is a bit of
            // a hack, but alternatives either require a more
            // significant redesign of the base class or are
            // inefficient.
            ReservoirState& s = const_cast<ReservoirState&>(reservoir_state);
            s.faceflux().resize(n);
            std::copy_n(flux.data(), n, s.faceflux().begin());
            if (asImpl().localWellsActive()) {
                const V& wflux = asImpl().wellModel().getStoredWellPerforationFluxes();
                assert(int(well_state.perfRates().size()) == wflux.size());
                std::copy_n(wflux.data(), wflux.size(), well_state.perfRates().begin());
            }
            return iter_report;
        }





        SolutionState
        variableState(const ReservoirState& x,
                      const WellState& xw) const
        {
            // As Base::variableState(), except making Sw and Xvar constants.
            std::vector<V> vars0 = asImpl().variableStateInitials(x, xw);
            std::vector<ADB> vars = ADB::variables(vars0);
            const std::vector<int> indices = asImpl().variableStateIndices();
            vars[indices[Sw]] = ADB::constant(vars[indices[Sw]].value());
            vars[indices[Xvar]] = ADB::constant(vars[indices[Xvar]].value());
            // OpmLog::debug("Injector pressure: " + std::to_string(vars[indices[Bhp]].value()[1]));
            return asImpl().variableStateExtractVars(x, indices, vars);
        }





        void computeAccum(const SolutionState& state,
                          const int            aix  )
        {
            if (aix == 0) {
                Base::computeAccum(state0_, aix);
            } else {
                Base::computeAccum(state, aix);
            }
        }





        void assembleMassBalanceEq(const SolutionState& state)
        {
            Base::assembleMassBalanceEq(state);

            // Compute the scaling factors.
            // Scaling factors are:
            //    1/bw,  1/bo - rs/bg, 1/bg - rv/bo
            assert(numPhases() == 3);
            assert(numMaterials() == 3);
            V one = V::Constant(state.pressure.size(), 1.0);
            scaling_[Water] = one / sd_.rq[Water].b;
            scaling_[Oil] = one / sd_.rq[Oil].b - state.rs / sd_.rq[Gas].b;
            scaling_[Gas] = one / sd_.rq[Gas].b - state.rv / sd_.rq[Oil].b;
            if (has_disgas_ && has_vapoil_) {
                ADB r_factor = one / (one - state.rs * state.rv);
                scaling_[Oil] = r_factor * scaling_[Oil];
                scaling_[Gas] = r_factor * scaling_[Gas];
            }
            // @TODO: multiply the oil and gas scale by 1/(1-rs*rv)?
            // OpmLog::debug("gas scaling: " + std::to_string(scaling_[2].value()[0]));
        }





        void updateState(const V& dx,
                         ReservoirState& reservoir_state,
                         WellState& well_state)
        {
            // Naively, rs and rv can get overwritten, so we
            // avoid that by storing.
            std::vector<double> rs_old = reservoir_state.gasoilratio();
            std::vector<double> rv_old = reservoir_state.rv();
            auto hs_old = reservoir_state.hydroCarbonState();
            auto phasecond_old = Base::phaseCondition_;
            auto isRs_old = Base::isRs_;
            auto isRv_old = Base::isRv_;
            auto isSg_old = Base::isSg_;

            // Compute the pressure range.
            const auto minmax_iters = std::minmax_element(reservoir_state.pressure().begin(),
                                                          reservoir_state.pressure().end());
            const double range = *minmax_iters.second - *minmax_iters.first;

            // Use the base class' updateState().
            Base::updateState(dx, reservoir_state, well_state);

            // Compute relative change.
            max_dp_rel_ = dx.head(reservoir_state.pressure().size()).abs().maxCoeff() / range;

            // Restore rs and rv, also various state flags.
            reservoir_state.gasoilratio() = rs_old;
            reservoir_state.rv() = rv_old;
            reservoir_state.hydroCarbonState() = hs_old;
            Base::phaseCondition_ = phasecond_old;
            Base::isRs_ = isRs_old;
            Base::isRv_ = isRv_old;
            Base::isSg_ = isSg_old;
        }





        bool getConvergence(const SimulatorTimerInterface& /* timer */, const int iteration)
        {
            const double tol_p = 1e-11;
            const double resmax = residual_.material_balance_eq[0].value().abs().maxCoeff();
            if (Base::terminalOutputEnabled()) {
                // Only rank 0 does print to std::cout
                if (iteration == 0) {
                    OpmLog::info("\nIter  Res(p)     Delta(p)\n");
                }
                std::ostringstream os;
                os.precision(3);
                os.setf(std::ios::scientific);
                os << std::setw(4) << iteration;
                os << std::setw(11) << resmax;
                os << std::setw(11) << max_dp_rel_;
                OpmLog::info(os.str());
            }
            return resmax < tol_p;
        }

    };





    /// Providing types by template specialisation of ModelTraits for BlackoilPressureModel.
    template <class Grid, class WellModel>
    struct ModelTraits< BlackoilPressureModel<Grid, WellModel> >
    {
        typedef BlackoilState ReservoirState;
        typedef WellStateFullyImplicitBlackoil WellState;
        typedef BlackoilModelParameters ModelParameters;
        typedef DefaultBlackoilSolutionState SolutionState;
    };

} // namespace Opm



#endif // OPM_BLACKOILPRESSUREMODEL_HEADER_INCLUDED
