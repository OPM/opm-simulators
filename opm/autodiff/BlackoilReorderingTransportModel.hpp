/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_BLACKOILREORDERINGTRANSPORTMODEL_HEADER_INCLUDED
#define OPM_BLACKOILREORDERINGTRANSPORTMODEL_HEADER_INCLUDED

#include <opm/autodiff/BlackoilModelBase.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/DebugTimeReport.hpp>
#include <opm/autodiff/multiPhaseUpwind.hpp>
#include <opm/core/grid.h>
#include <opm/core/transport/reorder/reordersequence.h>
#include <opm/core/simulator/BlackoilState.hpp>

#include <opm/autodiff/BlackoilTransportModel.hpp>

namespace Opm {


    namespace detail
    {

        template <typename Scalar>
        struct CreateVariable
        {
            Scalar operator()(double value, int index)
            {
                return Scalar::createVariable(value, index);
            }
        };



        template <>
        struct CreateVariable<double>
        {
            double operator()(double value, int)
            {
                return value;
            }
        };



        template <typename Scalar>
        struct CreateConstant
        {
            Scalar operator()(double value)
            {
                return Scalar::createConstant(value);
            }
        };



        template <>
        struct CreateConstant<double>
        {
            double operator()(double value)
            {
                return value;
            }
        };




        struct Connection
        {
            Connection(const int ind, const double s) : index(ind), sign(s) {}
            int index;
            double sign;
        };




        class Connections;



        class ConnectivityGraph
        {
        public:
            explicit ConnectivityGraph(const HelperOps& ops)
                : grad_(ops.grad)
                , div_(ops.div)
            {
                grad_ia_ = grad_.outerIndexPtr();
                grad_ja_ = grad_.innerIndexPtr();
                grad_sign_ = grad_.valuePtr();
                div_ia_ = div_.outerIndexPtr();
                div_ja_ = div_.innerIndexPtr();
                div_sign_ = div_.valuePtr();
            }

            Connections cellConnections(const int cell) const;

            std::array<int, 2> connectionCells(const int connection) const
            {
                const int pos = div_ia_[connection];
                assert(div_ia_[connection + 1] == pos + 2);
                const double sign1 = div_sign_[pos];
                assert(div_sign_[pos + 1] == -sign1);
                if (sign1 > 0.0) {
                    return {{ div_ja_[pos], div_ja_[pos + 1] }};
                } else {
                    return {{ div_ja_[pos + 1], div_ja_[pos] }};
                 }
            }

        private:
            friend class Connections;
            typedef HelperOps::M M;
            const M& grad_;
            const M& div_;
            const int* grad_ia_;
            const int* grad_ja_;
            const double* grad_sign_;
            const int* div_ia_;
            const int* div_ja_;
            const double* div_sign_;
        };


        class Connections
        {
        public:
            Connections(const ConnectivityGraph& cg, const int cell) : cg_(cg), cell_(cell) {}
            int size() const
            {
                return cg_.grad_ia_[cell_ + 1] - cg_.grad_ia_[cell_];
            }
            class Iterator
            {
            public:
                Iterator(const Connections& c, const int index) : c_(c), index_(index) {}
                Iterator& operator++()
                {
                    ++index_;
                    return *this;
                }
                bool operator!=(const Iterator& other)
                {
                    assert(&c_ == &other.c_);
                    return index_ != other.index_;
                }
                Connection operator*()
                {
                    assert(index_ >= 0 && index_ < c_.size());
                    const int pos = c_.cg_.grad_ia_[c_.cell_] + index_;
                    return Connection(c_.cg_.grad_ja_[pos], -c_.cg_.grad_sign_[pos]); // Note the minus sign!
                }
            private:
                const Connections& c_;
                int index_;
            };
            Iterator begin() const { return Iterator(*this, 0); }
            Iterator end() const { return Iterator(*this, size()); }
        private:
            friend class Iterator;
            const ConnectivityGraph& cg_;
            const int cell_;
        };


        inline Connections ConnectivityGraph::cellConnections(const int cell) const
        {
            return Connections(*this, cell);
        }



    } // namespace detail







    /// A model implementation for the transport equation in three-phase black oil.
    template<class Grid, class WellModel>
    class BlackoilReorderingTransportModel
        : public BlackoilModelBase<Grid, WellModel, BlackoilReorderingTransportModel<Grid, WellModel> >
    {
    public:
        typedef BlackoilModelBase<Grid, WellModel, BlackoilReorderingTransportModel<Grid, WellModel> > Base;
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
        BlackoilReorderingTransportModel(const typename Base::ModelParameters&   param,
                                         const Grid&                             grid,
                                         const BlackoilPropsAdFromDeck&          fluid,
                                         const DerivedGeology&                   geo,
                                         const RockCompressibility*              rock_comp_props,
                                         const StandardWells&                    std_wells,
                                         const NewtonIterationBlackoilInterface& linsolver,
                                         std::shared_ptr<const EclipseState>     eclState,
                                         const bool                              has_disgas,
                                         const bool                              has_vapoil,
                                         const bool                              terminal_output)
            : Base(param, grid, fluid, geo, rock_comp_props, std_wells, linsolver,
                   eclState, has_disgas, has_vapoil, terminal_output)
            , graph_(Base::ops_)
            , props_(dynamic_cast<const BlackoilPropsAdFromDeck&>(fluid)) // TODO: remove the need for this cast.
            , state0_{ ReservoirState(0, 0, 0), WellState(), V(), V() }
            , state_{ ReservoirState(0, 0, 0), WellState(), V(), V() }
            , tr_model_(param, grid, fluid, geo, rock_comp_props, std_wells, linsolver,
                        eclState, has_disgas, has_vapoil, terminal_output)
        {
            // Set up the common parts of the mass balance equations
            // for each active phase.
            const V transi = subset(geo_.transmissibility(), ops_.internal_faces);
            const V trans_nnc = ops_.nnc_trans;
            trans_all_ = V::Zero(transi.size() + trans_nnc.size());
            trans_all_ << transi, trans_nnc;
            gdz_ = geo_.gravity()[2] * (ops_.grad * geo_.z().matrix());
            rhos_ = DataBlock::Zero(ops_.div.rows(), 3);
            rhos_.col(Water) = props_.surfaceDensity(Water, Base::cells_);
            rhos_.col(Oil) = props_.surfaceDensity(Oil, Base::cells_);
            rhos_.col(Gas) = props_.surfaceDensity(Gas, Base::cells_);
        }





        void prepareStep(const SimulatorTimerInterface& timer,
                         const ReservoirState& reservoir_state,
                         const WellState& well_state)
        {
            tr_model_.prepareStep(timer, reservoir_state, well_state);
            Base::prepareStep(timer, reservoir_state, well_state);
            Base::param_.solve_welleq_initially_ = false;
            state0_.reservoir_state = reservoir_state;
            state0_.well_state = well_state;
            // Since (reference) pressure is constant, porosity and transmissibility multipliers can
            // be computed just once.
            const std::vector<double>& p = reservoir_state.pressure();
            state0_.tr_mult = Base::transMult(ADB::constant(Eigen::Map<const V>(p.data(), p.size()))).value();
            state0_.pv_mult = Base::poroMult(ADB::constant(Eigen::Map<const V>(p.data(), p.size()))).value();
            const int num_cells = p.size();
            cstate0_.resize(num_cells);
            for (int cell = 0; cell < num_cells; ++cell) {
                computeCellState(cell, state0_, cstate0_[cell]);
            }
            cstate_ = cstate0_;
        }





        template <class NonlinearSolverType>
        SimulatorReport nonlinearIteration(const int iteration,
                                           const SimulatorTimerInterface& timer,
                                           NonlinearSolverType& nonlinear_solver,
                                           ReservoirState& reservoir_state,
                                           const WellState& well_state)
        {
            // Extract reservoir and well fluxes and state.
            {
                DebugTimeReport tr("Extracting fluxes");
                extractFluxes(reservoir_state, well_state);
                extractState(reservoir_state, well_state);
            }

            // Compute cell ordering based on total flux.
            {
                DebugTimeReport tr("Topological sort");
                computeOrdering();
            }

            // Solve in every component (cell or block of cells), in order.
            {
                DebugTimeReport tr("Solving all components");
                for (int ii = 0; ii < 5; ++ii) {
                    DebugTimeReport tr2("Solving components single sweep.");
                    solveComponents();
                }
            }

            // Update states for output.
            reservoir_state = state_.reservoir_state;

            // Assemble with other model,
            {
                auto rs = reservoir_state;
                auto ws = well_state;
                tr_model_.nonlinearIteration(/*iteration*/ 0, timer, nonlinear_solver, rs, ws);
            }

            // Create report and exit.
            SimulatorReport report;
            report.converged = true;
            return report;
        }





        void afterStep(const SimulatorTimerInterface& /* timer */,
                       const ReservoirState& /* reservoir_state */,
                       const WellState& /* well_state */)
        {
            // Does nothing in this model.
        }

        using Base::numPhases;


    protected:

        // ============  Types  ============

        using Vec2 = Dune::FieldVector<double, 2>;
        using Mat22 = Dune::FieldMatrix<double, 2, 2>;
        using Eval = DenseAd::Evaluation<double, 2>;




        struct State
        {
            ReservoirState reservoir_state;
            WellState well_state;
            V tr_mult;
            V pv_mult;
        };





        template <typename ScalarT>
        struct CellState
        {
            using Scalar = ScalarT;

            Scalar s[3];
            Scalar rs;
            Scalar rv;
            Scalar p[3];
            Scalar kr[3];
            Scalar pc[3];
            Scalar temperature;
            Scalar mu[3];
            Scalar b[3];
            Scalar lambda[3];
            Scalar rho[3];
            Scalar rssat;
            Scalar rvsat;

            // Implement interface used for opm-material properties.
            const Scalar& saturation(int phaseIdx) const
            {
                return s[phaseIdx];
            }

            template <typename T>
            CellState<T> flatten() const
            {
                return CellState<T>{
                    { s[0].value(), s[1].value(), s[2].value() },
                    rs.value(),
                    rv.value(),
                    { p[0].value(), p[1].value(), p[2].value() },
                    { kr[0].value(), kr[1].value(), kr[2].value() },
                    { pc[0].value(), pc[1].value(), pc[2].value() },
                    temperature.value(),
                    { mu[0].value(), mu[1].value(), mu[2].value() },
                    { b[0].value(), b[1].value(), b[2].value() },
                    { lambda[0].value(), lambda[1].value(), lambda[2].value() },
                    { rho[0].value(), rho[1].value(), rho[2].value() },
                    rssat.value(),
                    rvsat.value()
                };
            }
        };



        // ============  Data members  ============

        using Base::grid_;
        using Base::geo_;
        using Base::ops_;

        const detail::ConnectivityGraph graph_;

        const BlackoilPropsAdFromDeck& props_;

        State state0_;
        State state_;

        std::vector<CellState<double>> cstate0_;
        std::vector<CellState<double>> cstate_;

        V total_flux_;
        V total_wellperf_flux_;
        DataBlock comp_wellperf_flux_;
        V total_wellflux_cell_;
        V oil_wellflux_cell_;
        V gas_wellflux_cell_;
        std::vector<int> sequence_;
        std::vector<int> components_;
        V trans_all_;
        V gdz_;
        DataBlock rhos_;

        std::array<double, 2> max_abs_dx_;
        std::array<int, 2> max_abs_dx_cell_;

        // TODO: remove this, for debug only.
        BlackoilTransportModel<Grid, WellModel> tr_model_;


        // ============  Member functions  ============





        template <typename Scalar>
        void computeCellState(const int cell, const State& state, CellState<Scalar>& cstate) const
        {
            assert(numPhases() == 3); // I apologize for this to my future self, that will have to fix it.

            // Extract from state and props.
            const auto hcstate = state.reservoir_state.hydroCarbonState()[cell];
            const bool is_sg = (hcstate == HydroCarbonState::GasAndOil);
            const bool is_rs = (hcstate == HydroCarbonState::OilOnly);
            const bool is_rv = (hcstate == HydroCarbonState::GasOnly);
            const double swval = state.reservoir_state.saturation()[3*cell + Water];
            const double sgval = state.reservoir_state.saturation()[3*cell + Gas];
            const double rsval = state.reservoir_state.gasoilratio()[cell];
            const double rvval = state.reservoir_state.rv()[cell];
            const double poval = state.reservoir_state.pressure()[cell];
            const int pvt_region = props_.pvtRegions()[cell];

            // Property functions.
            const auto& waterpvt = props_.waterProps();
            const auto& oilpvt = props_.oilProps();
            const auto& gaspvt = props_.gasProps();
            const auto& satfunc = props_.materialLaws();

            // Create saturation and composition variables.
            detail::CreateVariable<Scalar> variable;
            detail::CreateConstant<Scalar> constant;
            cstate.s[Water] = variable(swval, 0);
            cstate.s[Gas] = is_sg ? variable(sgval, 1) : constant(sgval);
            cstate.s[Oil] = 1.0 - cstate.s[Water] - cstate.s[Gas];
            cstate.rs = is_rs ? variable(rsval, 1) : constant(rsval);
            cstate.rv = is_rv ? variable(rvval, 1) : constant(rvval);

            // Compute relative permeabilities amd capillary pressures.
            const auto& params = satfunc.materialLawParams(cell);
            typedef BlackoilPropsAdFromDeck::MaterialLawManager::MaterialLaw MaterialLaw;
            MaterialLaw::relativePermeabilities(cstate.kr, params, cstate);
            MaterialLaw::capillaryPressures(cstate.pc, params, cstate);

            // Compute phase pressures.
            cstate.p[Oil] = constant(poval);
            cstate.p[Water] = cstate.p[Oil] + cstate.pc[Water]; // pcow = pw - po (!) [different from old convention]
            cstate.p[Gas] =   cstate.p[Oil] + cstate.pc[Gas];   // pcog = pg - po

            // Compute PVT properties.
            cstate.temperature = constant(0.0); // Temperature is not used.
            cstate.mu[Water] = waterpvt.viscosity(pvt_region, cstate.temperature, cstate.p[Water]);
            cstate.mu[Oil] = is_sg
                ? oilpvt.saturatedViscosity(pvt_region, cstate.temperature, cstate.p[Oil])
                : oilpvt.viscosity(pvt_region, cstate.temperature, cstate.p[Oil], cstate.rs);
            cstate.mu[Gas] = is_sg
                ? gaspvt.saturatedViscosity(pvt_region, cstate.temperature, cstate.p[Gas])
                : gaspvt.viscosity(pvt_region, cstate.temperature, cstate.p[Gas], cstate.rv);
            cstate.b[Water] = waterpvt.inverseFormationVolumeFactor(pvt_region, cstate.temperature, cstate.p[Water]);
            cstate.b[Oil] = is_sg
                ? oilpvt.saturatedInverseFormationVolumeFactor(pvt_region, cstate.temperature, cstate.p[Oil])
                : oilpvt.inverseFormationVolumeFactor(pvt_region, cstate.temperature, cstate.p[Oil], cstate.rs);
            cstate.b[Gas] = is_sg
                ? gaspvt.saturatedInverseFormationVolumeFactor(pvt_region, cstate.temperature, cstate.p[Gas])
                : gaspvt.inverseFormationVolumeFactor(pvt_region, cstate.temperature, cstate.p[Gas], cstate.rv);

            // Compute mobilities.
            for (int phase = 0; phase < 3; ++phase) {
                cstate.lambda[phase] = cstate.kr[phase] / cstate.mu[phase];
            }

            // Compute densities.
            cstate.rho[Water] = rhos_(cell, Water) * cstate.b[Water];
            cstate.rho[Oil] = (rhos_(cell, Oil) + cstate.rs*rhos_(cell, Gas)) * cstate.b[Oil]; // TODO: check that this is correct
            cstate.rho[Gas] = (rhos_(cell, Gas) + cstate.rv*rhos_(cell, Oil)) * cstate.b[Gas];

            // Compute saturated rs and rv factors.
            cstate.rssat = oilpvt.saturatedGasDissolutionFactor(pvt_region, cstate.temperature, cstate.p[Oil]);
            cstate.rvsat = gaspvt.saturatedOilVaporizationFactor(pvt_region, cstate.temperature, cstate.p[Gas]);
            // TODO: add vaporization controls such as in BlackoilPropsAdFromDeck::applyVap().
        }




        void extractFluxes(const ReservoirState& reservoir_state,
                           const WellState& well_state)
        {
            // Input face fluxes are for interior faces + nncs.
            total_flux_ = Eigen::Map<const V>(reservoir_state.faceflux().data(),
                                              reservoir_state.faceflux().size());
            total_wellperf_flux_ = Eigen::Map<const V>(well_state.perfRates().data(),
                                                       well_state.perfRates().size());
            comp_wellperf_flux_ = Eigen::Map<const DataBlock>(well_state.perfPhaseRates().data(),
                                                              well_state.perfRates().size(),
                                                              numPhases());
            const int num_cells = reservoir_state.pressure().size();
            total_wellflux_cell_ = superset(total_wellperf_flux_, Base::wellModel().wellOps().well_cells, num_cells);
            assert(Base::numPhases() == 3);
            V oilflux = comp_wellperf_flux_.col(1);
            V gasflux = comp_wellperf_flux_.col(2);
            oil_wellflux_cell_ = superset(oilflux, Base::wellModel().wellOps().well_cells, num_cells);
            gas_wellflux_cell_ = superset(gasflux, Base::wellModel().wellOps().well_cells, num_cells);
            assert(numPhases() * well_state.perfRates().size() == well_state.perfPhaseRates().size());
        }





        void extractState(const ReservoirState& reservoir_state,
                          const WellState& well_state)
        {
            state_.reservoir_state = reservoir_state;
            state_.well_state = well_state;
            const std::vector<double>& p = reservoir_state.pressure();
            state_.tr_mult = Base::transMult(ADB::constant(Eigen::Map<const V>(p.data(), p.size()))).value();
            state_.pv_mult = Base::poroMult(ADB::constant(Eigen::Map<const V>(p.data(), p.size()))).value();
        }





        void computeOrdering()
        {
            assert(!geo_.nnc().hasNNC()); // TODO: support compute_sequence() with grid + nnc.
            static_assert(std::is_same<Grid, UnstructuredGrid>::value,
                          "compute_sequence() is written in C and therefore requires an UnstructuredGrid, "
                          "it must be rewritten to use other grid classes such as CpGrid");
            using namespace Opm::AutoDiffGrid;
            const int num_cells = numCells(grid_);
            sequence_.resize(num_cells);
            components_.resize(num_cells + 1); // max possible size
            int num_components = -1;


            using namespace Opm::AutoDiffGrid;
            const int num_faces = numFaces(grid_);
            V flux_on_all_faces = superset(total_flux_, ops_.internal_faces, num_faces);
            compute_sequence(&grid_, flux_on_all_faces.data(), sequence_.data(), components_.data(), &num_components);
            OpmLog::debug(std::string("Number of components: ") + std::to_string(num_components));
            components_.resize(num_components + 1); // resize to fit actually used part
        }




        void solveComponents()
        {
            // Zero the max changed.
            max_abs_dx_[0] = 0.0;
            max_abs_dx_[1] = 0.0;
            max_abs_dx_cell_[0] = -1;
            max_abs_dx_cell_[1] = -1;

            // Solve the equations.
            const int num_components = components_.size() - 1;
            for (int comp = 0; comp < num_components; ++comp) {
                const int comp_size = components_[comp + 1] - components_[comp];
                if (comp_size == 1) {
                    solveSingleCell(sequence_[components_[comp]]);
                } else {
                    solveMultiCell(comp_size, &sequence_[components_[comp]]);
                }
            }

            // Log the max change.
            {
                std::ostringstream os;
                os << "===  Max abs dx[0]: " << max_abs_dx_[0] << " (cell " << max_abs_dx_cell_[0]
                   <<")  dx[1]: " << max_abs_dx_[1] << " (cell " << max_abs_dx_cell_[1] << ")";
                OpmLog::debug(os.str());
            }
        }





        void solveSingleCell(const int cell)
        {

            Vec2 res;
            Mat22 jac;
            assembleSingleCell(cell, res, jac);

            // Newton loop.
            int iter = 0;
            const int max_iter = 100;
            double relaxation = 1.0;
            while (!getConvergence(cell, res) && iter < max_iter) {
                Vec2 dx;
                jac.solve(dx, res);
                dx *= relaxation;
                const auto hcstate_old = state_.reservoir_state.hydroCarbonState()[cell];
                updateState(cell, -dx);
                const auto hcstate = state_.reservoir_state.hydroCarbonState()[cell];
                assembleSingleCell(cell, res, jac);
                ++iter;
                if (iter > 10) {
                    relaxation = 0.85;
                    if (iter > 15) {
                        relaxation = 0.70;
                    }
                    if (iter > 20) {
                        relaxation = 0.55;
                    }
                    if (iter > 25) {
                        relaxation = 0.40;
                    }
                    if (iter > 30) {
                        relaxation = 0.25;
                    }
                    // std::ostringstream os;
                    // os << "Iteration " << iter << " in cell " << cell << ", residual = " << res
                    //    << ", cell values { s = ( " << cstate_[cell].s[Water] << ", " << cstate_[cell].s[Oil] << ", " << cstate_[cell].s[Gas]
                    //    << " ), rs = " << cstate_[cell].rs << ", rv = " << cstate_[cell].rv << "}, dx = " << dx << ", hcstate: " << hcstate_old << " -> " << hcstate;
                    // OpmLog::debug(os.str());
                }
            }
            if (iter == max_iter) {
                std::ostringstream os;
                os << "Failed to converge in cell " << cell << ", residual = " << res
                   << ", cell values { s = ( " << cstate_[cell].s[Water] << ", " << cstate_[cell].s[Oil] << ", " << cstate_[cell].s[Gas]
                   << " ), rs = " << cstate_[cell].rs << ", rv = " << cstate_[cell].rv << " }";
                OpmLog::debug(os.str());
            }
        }





        void solveMultiCell(const int comp_size, const int* cell_array)
        {
            // OpmLog::warning("solveMultiCell", "solveMultiCell() called with component size " + std::to_string(comp_size));
            for (int ii = 0; ii < comp_size; ++ii) {
                solveSingleCell(cell_array[ii]);
            }
        }




        template <typename Scalar>
        Scalar oilAccumulation(const CellState<Scalar>& cs)
        {
            return cs.b[Oil]*cs.s[Oil] + cs.rv*cs.b[Gas]*cs.s[Gas];
        }




        template <typename Scalar>
        Scalar gasAccumulation(const CellState<Scalar>& cs)
        {
            return cs.b[Gas]*cs.s[Gas] + cs.rs*cs.b[Oil]*cs.s[Oil];
        }




        void applyThresholdPressure(const int connection, Eval& dp)
        {
            const double thres_press = Base::threshold_pressures_by_connection_[connection];
            if (std::fabs(dp.value()) < thres_press) {
                dp.setValue(0.0);
            } else {
                dp -= dp.value() > 0.0 ? thres_press : -thres_press;
            }
        }




        void assembleSingleCell(const int cell, Vec2& res, Mat22& jac)
        {
            assert(numPhases() == 3); // I apologize for this to my future self, that will have to fix it.

            CellState<Eval> st;
            computeCellState(cell, state_, st);
            cstate_[cell] = st.template flatten<double>();

            // Accumulation terms.
            const double pvm0 = state0_.pv_mult[cell];
            const double pvm = state_.pv_mult[cell];
            const double ao0 = oilAccumulation(cstate0_[cell]) * pvm0;
            const Eval ao  = oilAccumulation(st) * pvm;
            const double ag0 = gasAccumulation(cstate0_[cell]) * pvm0;
            const Eval ag  = gasAccumulation(st) * pvm;

            // Flux terms.
            Eval div_oilflux = Eval::createConstant(0.0);
            Eval div_gasflux = Eval::createConstant(0.0);
            for (auto conn : graph_.cellConnections(cell)) {
                auto conn_cells = graph_.connectionCells(conn.index);
                const int from = conn_cells[0];
                const int to = conn_cells[1];
                if (from < 0 || to < 0) {
                    continue; // Boundary.
                }
                assert((from == cell) == (conn.sign > 0.0));
                const int other = from == cell ? to : from;
                const double vt = conn.sign * total_flux_[conn.index];
                const double gdz = conn.sign * gdz_[conn.index];

                // From this point, we treat everything about this
                // connection as going from 'cell' to 'other'. Since
                // we don't want derivatives from the 'other' cell to
                // participate in the solution, we use the constant
                // values from cstate_[other].
                Eval dh[3];
                Eval dh_sat[3];
                const Eval grad_oil_press = cstate_[other].p[Oil] - st.p[Oil];
                for (int phase : { Water, Oil, Gas }) {
                    const Eval gradp = cstate_[other].p[phase] - st.p[phase];
                    const Eval rhoavg = 0.5 * (st.rho[phase] + cstate_[other].rho[phase]);
                    dh[phase] = gradp - rhoavg * gdz;
                    if (Base::use_threshold_pressure_) {
                        applyThresholdPressure(conn.index, dh[phase]);
                    }
                    dh_sat[phase] = grad_oil_press - dh[phase];
                }
                const double tran = trans_all_[conn.index]; // TODO: include tr_mult effect.
                const auto& m1 = st.lambda;
                const auto& m2 = cstate_[other].lambda;
                const auto upw = connectionMultiPhaseUpwind({{ dh_sat[Water].value(), dh_sat[Oil].value(), dh_sat[Gas].value() }},
                                                            {{ m1[Water].value(), m1[Oil].value(), m1[Gas].value() }},
                                                            {{ m2[Water], m2[Oil], m2[Gas] }},
                                                            tran, vt);
                // if (upw[0] != upw[1] || upw[1] != upw[2]) {
                //     OpmLog::debug("Detected countercurrent flow between cells " + std::to_string(from) + " and " + std::to_string(to));
                // }
                Eval b[3];
                Eval mob[3];
                Eval tot_mob = Eval::createConstant(0.0);
                for (int phase : { Water, Oil, Gas }) {
                    b[phase] = upw[phase] > 0.0 ? st.b[phase] : cstate_[other].b[phase];
                    mob[phase] = upw[phase] > 0.0 ? m1[phase] : m2[phase];
                    tot_mob += mob[phase];
                }
                Eval rs = upw[Oil] > 0.0 ? st.rs : cstate_[other].rs;
                Eval rv = upw[Gas] > 0.0 ? st.rv : cstate_[other].rv;

                Eval flux[3];
                for (int phase : { Oil, Gas }) {
                    Eval gflux = Eval::createConstant(0.0);
                    for (int other_phase : { Water, Oil, Gas }) {
                        if (phase != other_phase) {
                            gflux += mob[other_phase] * (dh_sat[phase] - dh_sat[other_phase]);
                        }
                    }
                    flux[phase] = b[phase] * (mob[phase] / tot_mob) * (vt + tran*gflux);
                }
                div_oilflux += flux[Oil] + rv*flux[Gas];
                div_gasflux += flux[Gas] + rs*flux[Oil];
            }

            // Well fluxes.
            if (total_wellflux_cell_[cell] > 0.0) {
                // Injecting perforation. Use given phase rates.
                assert(oil_wellflux_cell_[cell] >= 0.0);
                assert(gas_wellflux_cell_[cell] >= 0.0);
                div_oilflux -= oil_wellflux_cell_[cell];
                div_gasflux -= gas_wellflux_cell_[cell];
            } else if (total_wellflux_cell_[cell] < 0.0) {
                // Producing perforation. Use total rate and fractional flow.
                Eval totmob = st.lambda[Water] + st.lambda[Oil] + st.lambda[Gas];
                Eval oilflux = st.b[Oil] * (st.lambda[Oil]/totmob) *  total_wellflux_cell_[cell];
                Eval gasflux = st.b[Gas] * (st.lambda[Gas]/totmob) *  total_wellflux_cell_[cell];
                div_oilflux -= (oilflux + st.rv * gasflux);
                div_gasflux -= (gasflux + st.rs * oilflux);
            }

            const Eval oileq = Base::pvdt_[cell]*(ao - ao0) + div_oilflux;
            const Eval gaseq = Base::pvdt_[cell]*(ag - ag0) + div_gasflux;

            res[0] = oileq.value();
            res[1] = gaseq.value();
            jac[0][0] = oileq.derivative(0);
            jac[0][1] = oileq.derivative(1);
            jac[1][0] = gaseq.derivative(0);
            jac[1][1] = gaseq.derivative(1);
        }





        bool getConvergence(const int cell, const Vec2& res)
        {
            const double tol = 1e-7;
            // Compute scaled residuals (scaled like saturations).
            double sres[] = { res[0] / (cstate_[cell].b[Oil] * Base::pvdt_[cell]),
                              res[1] / (cstate_[cell].b[Gas] * Base::pvdt_[cell]) };
            return std::fabs(sres[0]) < tol && std::fabs(sres[1]) < tol;
        }




        void updateState(const int cell,
                         const Vec2& dx)
        {
            if (std::fabs(dx[0]) > max_abs_dx_[0]) {
                max_abs_dx_cell_[0] = cell;
            }
            if (std::fabs(dx[1]) > max_abs_dx_[1]) {
                max_abs_dx_cell_[1] = cell;
            }
            max_abs_dx_[0] = std::max(max_abs_dx_[0], std::fabs(dx[0]));
            max_abs_dx_[1] = std::max(max_abs_dx_[1], std::fabs(dx[1]));

            // Get saturation updates.
            const double dsw = dx[0];
            double dsg = 0.0;
            auto& hcstate = state_.reservoir_state.hydroCarbonState()[cell];
            if (hcstate == HydroCarbonState::GasAndOil) {
                dsg = dx[1];
            } else if (hcstate == HydroCarbonState::GasOnly) {
                dsg = -dsw;
            }
            const double dso = -(dsw + dsg);

            // Handle too large saturation changes.
            const double maxval = std::max(std::fabs(dsw), std::max(std::fabs(dso), std::fabs(dsg)));
            const double sfactor = std::min(1.0, Base::dsMax() / maxval);
            double* s = state_.reservoir_state.saturation().data() + 3*cell;
            s[Water] += sfactor*dsw;
            s[Gas] += sfactor*dsg;
            s[Oil] = 1.0 - s[Water] - s[Gas];

            // Handle < 0 saturations.
            for (int phase : { Gas, Oil, Water }) { // TODO: check if ordering here is significant
                if (s[phase] < 0.0) {
                    for (int other_phase : { Water, Oil, Gas }) {
                        if (phase != other_phase) {
                            s[other_phase] /= (1.0 - s[phase]);
                        }
                    }
                    s[phase] = 0.0;
                }
            }

            // Update rs.
            double& rs = state_.reservoir_state.gasoilratio()[cell];
            const double rs_old = rs;
            if (hcstate == HydroCarbonState::OilOnly) {
                // const double max_allowed_change = std::fabs(rs_old) * Base::drMaxRel();
                const double drs = dx[1];
                // const double factor = std::min(1.0, max_allowed_change / std::fabs(drs));
                // rs += factor*drs;
                rs += drs;
                rs = std::max(rs, 0.0);
            }

            // Update rv.
            double& rv = state_.reservoir_state.rv()[cell];
            const double rv_old = rv;
            if (hcstate == HydroCarbonState::GasOnly) {
                // const double max_allowed_change = std::fabs(rv_old) * Base::drMaxRel();
                const double drv = dx[1];
                // const double factor = std::min(1.0, max_allowed_change / std::fabs(drv));
                // rv += factor*drv;
                rv += drv;
                rv = std::max(rv, 0.0);
            }

            const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
            const bool water_only = s[Water] > (1 - epsilon);
            const auto old_hcstate = hcstate;
            hcstate = HydroCarbonState::GasAndOil;
            // sg <-> rs transition.
            {
                const double rssat_old = cstate_[cell].rssat;
                const double rssat = rssat_old; // TODO: This is no longer true with vaporization controls
                const bool is_rs = old_hcstate == HydroCarbonState::OilOnly;
                const bool has_gas = (s[Gas] > 0.0 && !is_rs);
                const bool gas_vaporized = ( (rs > rssat * (1+epsilon) && is_rs ) && (rs_old > rssat_old * (1-epsilon)) );
                if (water_only || has_gas || gas_vaporized) {
                    rs = rssat;
                } else {
                    hcstate = HydroCarbonState::OilOnly;
                }
            }

            // sg <-> rv transition.
            {
                const double rvsat_old = cstate_[cell].rvsat;
                const double rvsat = rvsat_old; // TODO: This is no longer true with vaporization controls
                const bool is_rv = old_hcstate == HydroCarbonState::GasOnly;
                const bool has_oil = (s[Oil] > 0.0 && !is_rv);
                const bool oil_condensed = ( (rv > rvsat * (1+epsilon) && is_rv) && (rv_old > rvsat_old * (1-epsilon)) );
                if (water_only || has_oil || oil_condensed) {
                    rv = rvsat;
                } else {
                    hcstate = HydroCarbonState::GasOnly;
                }
            }
        }
    };








    /// Providing types by template specialisation of ModelTraits for BlackoilReorderingTransportModel.
    template <class Grid, class WellModel>
    struct ModelTraits< BlackoilReorderingTransportModel<Grid, WellModel> >
    {
        typedef BlackoilState ReservoirState;
        typedef WellStateFullyImplicitBlackoil WellState;
        typedef BlackoilModelParameters ModelParameters;
        typedef DefaultBlackoilSolutionState SolutionState;
    };

} // namespace Opm




#endif // OPM_BLACKOILREORDERINGTRANSPORTMODEL_HEADER_INCLUDED
