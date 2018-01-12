/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS

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


#ifndef OPM_WELLINTERFACE_HEADER_INCLUDED
#define OPM_WELLINTERFACE_HEADER_INCLUDED

#include <opm/common/OpmLog/OpmLog.hpp>


#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/wells/WellsManager.hpp>

#include <opm/autodiff/VFPProperties.hpp>
#include <opm/autodiff/VFPInjProperties.hpp>
#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/autodiff/WellHelpers.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>
#include <opm/autodiff/RateConverter.hpp>

#include <opm/simulators/WellSwitchingLogger.hpp>

#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/matrixmatrix.hh>

#include <opm/material/densead/Math.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <string>
#include <memory>
#include <vector>
#include <cassert>

namespace Opm
{


    template<typename TypeTag>
    class WellInterface
    {
    public:

        using WellState = WellStateFullyImplicitBlackoil;

        typedef BlackoilModelParameters ModelParameters;

        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;

        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
        typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

        static const int numEq = Indices::numEq;
        typedef double Scalar;

        typedef Dune::FieldVector<Scalar, numEq    > VectorBlockType;
        typedef Dune::FieldMatrix<Scalar, numEq, numEq > MatrixBlockType;
        typedef Dune::BCRSMatrix <MatrixBlockType> Mat;
        typedef Dune::BlockVector<VectorBlockType> BVector;
        typedef DenseAd::Evaluation<double, /*size=*/numEq> Eval;

        typedef Ewoms::BlackOilPolymerModule<TypeTag> PolymerModule;

        static const bool has_solvent = GET_PROP_VALUE(TypeTag, EnableSolvent);
        static const bool has_polymer = GET_PROP_VALUE(TypeTag, EnablePolymer);
        static const int contiSolventEqIdx = Indices::contiSolventEqIdx;
        static const int contiPolymerEqIdx = Indices::contiPolymerEqIdx;

        // For the conversion between the surface volume rate and resrevoir voidage rate
        using RateConverterType = RateConverter::
        SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;

        /// Constructor
        WellInterface(const Well* well, const int time_step, const Wells* wells,
                      const ModelParameters& param,
                      const RateConverterType& rate_converter,
                      const int pvtRegionIdx,
                      const int num_components);

        /// Virutal destructor
        virtual ~WellInterface() {}

        /// Well name.
        const std::string& name() const;

        /// Well cells.
        const std::vector<int>& cells() {return well_cells_; }

        /// Well type, INJECTOR or PRODUCER.
        WellType wellType() const;

        /// Well controls
        WellControls* wellControls() const;

        void setVFPProperties(const VFPProperties* vfp_properties_arg);

        virtual void init(const PhaseUsage* phase_usage_arg,
                          const std::vector<double>& depth_arg,
                          const double gravity_arg,
                          const size_t num_cells);

        virtual void initPrimaryVariablesEvaluation() const = 0;

        /// a struct to collect information about the convergence checking
        struct ConvergenceReport {
            struct ProblemWell {
                std::string well_name;
                std::string phase_name;
            };
            bool converged = true;
            bool nan_residual_found = false;
            std::vector<ProblemWell> nan_residual_wells;
            // We consider Inf is large residual here
            bool too_large_residual_found = false;
            std::vector<ProblemWell> too_large_residual_wells;

            ConvergenceReport& operator+=(const ConvergenceReport& rhs) {
                converged = converged && rhs.converged;
                nan_residual_found = nan_residual_found || rhs.nan_residual_found;
                if (rhs.nan_residual_found) {
                    for (const ProblemWell& well : rhs.nan_residual_wells) {
                        nan_residual_wells.push_back(well);
                    }
                }
                too_large_residual_found = too_large_residual_found || rhs.too_large_residual_found;
                if (rhs.too_large_residual_found) {
                    for (const ProblemWell& well : rhs.too_large_residual_wells) {
                        too_large_residual_wells.push_back(well);
                    }
                }
                return *this;
            }
        };

        virtual ConvergenceReport getWellConvergence(const std::vector<double>& B_avg) const = 0;

        virtual void solveEqAndUpdateWellState(WellState& well_state) = 0;

        virtual void assembleWellEq(Simulator& ebosSimulator,
                                    const double dt,
                                    WellState& well_state,
                                    bool only_wells) = 0;

        void updateListEconLimited(const WellState& well_state,
                                   DynamicListEconLimited& list_econ_limited) const;

        void setWellEfficiencyFactor(const double efficiency_factor);

        void computeRepRadiusPerfLength(const Grid& grid, const std::map<int, int>& cartesian_to_compressed);

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        virtual void recoverWellSolutionAndUpdateWellState(const BVector& x,
                                                           WellState& well_state) const = 0;

        /// Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const = 0;

        /// r = r - C D^-1 Rw
        virtual void apply(BVector& r) const = 0;

        // TODO: before we decide to put more information under mutable, this function is not const
        virtual void computeWellPotentials(const Simulator& ebosSimulator,
                                           const WellState& well_state,
                                           std::vector<double>& well_potentials) = 0;

        virtual void updateWellStateWithTarget(WellState& well_state) const = 0;

        void updateWellControl(WellState& well_state,
                               wellhelpers::WellSwitchingLogger& logger) const;

        virtual void updatePrimaryVariables(const WellState& well_state) const = 0;

        virtual void calculateExplicitQuantities(const Simulator& ebosSimulator,
                                                 const WellState& well_state) = 0; // should be const?

    protected:

        // to indicate a invalid connection
        static const int INVALIDCONNECTION = -100000;

        const Well* well_ecl_;

        const int current_step_;

        // the index of well in Wells struct
        int index_of_well_;

        // simulation parameters
        const ModelParameters& param_;

        // well type
        // INJECTOR or PRODUCER
        enum WellType well_type_;

        // number of phases
        int number_of_phases_;

        // component fractions for each well
        // typically, it should apply to injection wells
        std::vector<double> comp_frac_;

        // controls for this well
        struct WellControls* well_controls_;

        // number of the perforations for this well
        int number_of_perforations_;

        // record the index of the first perforation
        // of states of individual well.
        int first_perf_;

        // well index for each perforation
        std::vector<double> well_index_;

        // depth for each perforation
        std::vector<double> perf_depth_;

        // reference depth for the BHP
        double ref_depth_;

        double well_efficiency_factor_;

        // cell index for each well perforation
        std::vector<int> well_cells_;

        // saturation table nubmer for each well perforation
        std::vector<int> saturation_table_number_;

        // representative radius of the perforations, used in shear calculation
        std::vector<double> perf_rep_radius_;

        // length of the perforations, use in shear calculation
        std::vector<double> perf_length_;

        // well bore diameter
        std::vector<double> bore_diameters_;

        const PhaseUsage* phase_usage_;

        bool getAllowCrossFlow() const;

        const VFPProperties* vfp_properties_;

        double gravity_;

        // For the conversion between the surface volume rate and resrevoir voidage rate
        const RateConverterType& rateConverter_;

        // The pvt region of the well. We assume
        // We assume a well to not penetrate more than one pvt region.
        const int pvtRegionIdx_;

        const int num_components_;

        const PhaseUsage& phaseUsage() const;

        int flowPhaseToEbosCompIdx( const int phaseIdx ) const;

        int ebosCompIdxToFlowCompIdx( const unsigned compIdx ) const;

        double wsolvent() const;

        double wpolymer() const;

        bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                                 const WellState& well_state) const;

        bool wellHasTHPConstraints() const;

        // Component fractions for each phase for the well
        const std::vector<double>& compFrac() const;

        double mostStrictBhpFromBhpLimits() const;

        // a tuple type for ratio limit check.
        // first value indicates whether ratio limit is violated, when the ratio limit is not violated, the following three
        // values should not be used.
        // second value indicates whehter there is only one connection left.
        // third value indicates the indx of the worst-offending connection.
        // the last value indicates the extent of the violation for the worst-offending connection, which is defined by
        // the ratio of the actual value to the value of the violated limit.
        using RatioCheckTuple = std::tuple<bool, bool, int, double>;

        RatioCheckTuple checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                                              const WellState& well_state) const;

        RatioCheckTuple checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                                             const WellState& well_state) const;

        double scalingFactor(const int comp_idx) const;


    };

}

#include "WellInterface_impl.hpp"

#endif // OPM_WELLINTERFACE_HEADER_INCLUDED
