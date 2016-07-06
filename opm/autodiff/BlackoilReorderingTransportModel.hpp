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
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/core/grid.h>
#include <opm/autodiff/DebugTimeReport.hpp>
#include <opm/core/transport/reorder/reordersequence.h>

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
                                         const BlackoilPropsAdInterface&         fluid,
                                         const DerivedGeology&                   geo,
                                         const RockCompressibility*              rock_comp_props,
                                         const StandardWells&                    std_wells,
                                         const NewtonIterationBlackoilInterface& linsolver,
                                         Opm::EclipseStateConstPtr               eclState,
                                         const bool                              has_disgas,
                                         const bool                              has_vapoil,
                                         const bool                              terminal_output)
            : Base(param, grid, fluid, geo, rock_comp_props, std_wells, linsolver,
                   eclState, has_disgas, has_vapoil, terminal_output)
            , props_(dynamic_cast<const BlackoilPropsAdFromDeck&>(fluid)) // TODO: remove the need for this cast.
            , state0_{ ReservoirState(0, 0, 0), WellState(), V() }
            , state_{ ReservoirState(0, 0, 0), WellState(), V() }
        {
            // Set up the common parts of the mass balance equations
            // for each active phase.
            const V transi = subset(geo_.transmissibility(), ops_.internal_faces);
            const V trans_nnc = ops_.nnc_trans;
            trans_all_ = V::Zero(transi.size() + trans_nnc.size());
            trans_all_ << transi, trans_nnc;
            gdz_ = geo_.gravity()[2] * (ops_.grad * geo_.z().matrix());
        }





        void prepareStep(const double dt,
                         const ReservoirState& reservoir_state,
                         const WellState& well_state)
        {
            Base::prepareStep(dt, reservoir_state, well_state);
            Base::param_.solve_welleq_initially_ = false;
            state0_.reservoir_state = reservoir_state;
            state0_.well_state = well_state;
            // Since (reference) pressure is constant, porosity and transmissibility multipliers can
            // be computed just once.
            const std::vector<double>& p = reservoir_state.pressure();
            Base::pvdt_ *= Base::poroMult(ADB::constant(Eigen::Map<const V>(p.data(), p.size()))).value();
            state0_.tr_mult = Base::transMult(ADB::constant(Eigen::Map<const V>(p.data(), p.size()))).value();
        }





        template <class NonlinearSolverType>
        IterationReport nonlinearIteration(const int /* iteration */,
                                           const double /* dt */,
                                           NonlinearSolverType& /* nonlinear_solver */,
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
                solveComponents();
            }

            // Update states for output.

            // Create report and exit.
            const bool failed = false;
            const bool converged = true;
            const int linear_iterations = 0;
            const int well_iterations = std::numeric_limits<int>::min();
            return IterationReport{failed, converged, linear_iterations, well_iterations};
        }





        void afterStep(const double /* dt */,
                       const ReservoirState& /* reservoir_state */,
                       const WellState& /* well_state */)
        {
            // Does nothing in this model.
        }





        using Base::numPhases;


    protected:

        // ============  Data members  ============
        using Base::grid_;
        using Base::geo_;
        using Base::ops_;
        using Vec2 = Dune::FieldVector<double, 2>;
        using Mat22 = Dune::FieldMatrix<double, 2, 2>;

        const BlackoilPropsAdFromDeck& props_;

        struct State
        {
            ReservoirState reservoir_state;
            WellState well_state;
            V tr_mult;
        };

        State state0_;
        State state_;

        V total_flux_;
        V total_wellperf_flux_;
        DataBlock comp_wellperf_flux_;
        std::vector<int> sequence_;
        std::vector<int> components_;
        V trans_all_;
        V gdz_;



        // ============  Member functions  ============


        void extractFluxes(const ReservoirState& reservoir_state,
                           const WellState& well_state)
        {
            // Input face fluxes are for interior faces only, while rest of code deals with all faces.
            const V face_flux = Eigen::Map<const V>(reservoir_state.faceflux().data(),
                                                    reservoir_state.faceflux().size());
            using namespace Opm::AutoDiffGrid;
            const int num_faces = numFaces(grid_);
            assert(face_flux.size() < num_faces); // Expected to be internal only.
            total_flux_ = superset(face_flux, ops_.internal_faces, num_faces);
            total_wellperf_flux_ = Eigen::Map<const V>(well_state.perfRates().data(),
                                                       well_state.perfRates().size());
            comp_wellperf_flux_ = Eigen::Map<const DataBlock>(well_state.perfPhaseRates().data(),
                                                              well_state.perfRates().size(),
                                                              numPhases());
            assert(numPhases() * well_state.perfRates().size() == well_state.perfPhaseRates().size());
        }





        void extractState(const ReservoirState& reservoir_state,
                          const WellState& well_state)
        {
            state_.reservoir_state = reservoir_state;
            state_.well_state = well_state;
            const std::vector<double>& p = reservoir_state.pressure();
            state_.tr_mult = Base::transMult(ADB::constant(Eigen::Map<const V>(p.data(), p.size()))).value();
        }





        void computeOrdering()
        {
            static_assert(std::is_same<Grid, UnstructuredGrid>::value,
                          "compute_sequence() is written in C and therefore requires an UnstructuredGrid, "
                          "it must be rewritten to use other grid classes such as CpGrid");
            using namespace Opm::AutoDiffGrid;
            const int num_cells = numCells(grid_);
            sequence_.resize(num_cells);
            components_.resize(num_cells + 1); // max possible size
            int num_components = -1;
            compute_sequence(&grid_, total_flux_.data(), sequence_.data(), components_.data(), &num_components);
            OpmLog::debug(std::string("Number of components: ") + std::to_string(num_components));
            components_.resize(num_components + 1); // resize to fit actually used part
        }




        void solveComponents()
        {
            const int num_components = components_.size() - 1;
            for (int comp = 0; comp < num_components; ++comp) {
                const int comp_size = components_[comp + 1] - components_[comp];
                if (comp_size == 1) {
                    solveSingleCell(sequence_[components_[comp]]);
                } else {
                    solveMultiCell(comp_size, &sequence_[components_[comp]]);
                }
            }
        }





        void solveSingleCell(const int cell)
        {

            // Vec2 x = getInitialGuess(cell);
            Vec2 x;
            Vec2 res;
            Mat22 jac;
            assembleSingleCell(cell, x, res, jac);

            // Newton loop.
            while (!getConvergence(res)) {
                Vec2 dx;
                jac.solve(dx, res);
                x = x - dx;
                assembleSingleCell(cell, x, res, jac);
            }
        }





        void solveMultiCell(const int comp_size, const int* cell_array)
        {
            OpmLog::warning("solveMultiCell", "solveMultiCell() called with component size " + std::to_string(comp_size));
            for (int ii = 0; ii < comp_size; ++ii) {
                solveSingleCell(cell_array[ii]);
            }
        }




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

            // Implement interface used for opm-material properties.
            const Scalar& saturation(int phaseIdx) const
            {
                return s[phaseIdx];
            }
        };




        template <typename Scalar>
        void computeCellState(const int cell, const State& state, CellState<Scalar>& cstate)
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
            cstate.p[Water] = cstate.p[Oil] - cstate.pc[Water]; // pcow = po - pw
            cstate.p[Gas] =   cstate.p[Oil] + cstate.pc[Gas];   // pcog = pg - po (!)

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
        }




        void assembleSingleCell(const int cell, const Vec2& x, Vec2& res, Mat22& jac)
        {
            assert(numPhases() == 3); // I apologize for this to my future self, that will have to fix it.

            CellState<double> cstate0;
            computeCellState(cell, state0_, cstate0);
            typedef DenseAd::Evaluation<double, 2> Eval;
            CellState<Eval> cstate;
            computeCellState(cell, state_, cstate);


            // const Eval ao0 = 0.0;
            // const Eval oileq = Base::pvdt_[cell];

            res[0] = Base::pvdt_[cell]*(0.0);
            jac[0][0] = 0.0;

        }





        // Vec2 getInitialGuess(const State& state, const int cell)
        // {
        //     const auto& hcs = state.hydroCarbonState();
        //     double xvar = (hcs[cell] == GasAndOil) ? sg_[cell]
        //         : ((hcs[cell] == OilOnly) ? rs_[cell] : rv_[cell]);
        //     return Vec2{sw_[cell], xvar};
        // }





        bool getConvergence(const Vec2& res)
        {
            const double tol = 1e-6;
            return res[0] < tol && res[1] < tol;
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
