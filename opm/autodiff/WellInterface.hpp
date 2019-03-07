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
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleEnums.hpp>


#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/simulator/SimulatorReport.hpp>

#include <opm/autodiff/VFPProperties.hpp>
#include <opm/autodiff/WellHelpers.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilModelParametersEbos.hpp>
#include <opm/autodiff/RateConverter.hpp>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/DeferredLogger.hpp>

#include <ewoms/models/blackoil/blackoilpolymermodules.hh>
#include <ewoms/models/blackoil/blackoilsolventmodules.hh>

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

        typedef BlackoilModelParametersEbos<TypeTag> ModelParameters;

        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;

        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
        typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
        typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
        typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;

        static const int numEq = Indices::numEq;
        typedef double Scalar;

        typedef Dune::FieldVector<Scalar, numEq    > VectorBlockType;
        typedef Dune::FieldMatrix<Scalar, numEq, numEq > MatrixBlockType;
        typedef typename SparseMatrixAdapter::IstlMatrix Mat;
        typedef Dune::BlockVector<VectorBlockType> BVector;
        typedef DenseAd::Evaluation<double, /*size=*/numEq> Eval;

        typedef Ewoms::BlackOilPolymerModule<TypeTag> PolymerModule;

        static const bool has_solvent = GET_PROP_VALUE(TypeTag, EnableSolvent);
        static const bool has_polymer = GET_PROP_VALUE(TypeTag, EnablePolymer);
        static const bool has_energy = GET_PROP_VALUE(TypeTag, EnableEnergy);
        // flag for polymer molecular weight related
        static const bool has_polymermw = GET_PROP_VALUE(TypeTag, EnablePolymerMW);
        static const int contiSolventEqIdx = Indices::contiSolventEqIdx;
        static const int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
        // index for the polymer molecular weight continuity equation
        static const int contiPolymerMWEqIdx = Indices::contiPolymerMWEqIdx;

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

        /// Index of well in the wells struct and wellState
        const int indexOfWell() const;

        /// Well cells.
        const std::vector<int>& cells() const {return well_cells_; }

        /// Well type, INJECTOR or PRODUCER.
        WellType wellType() const;

        /// Well controls
        WellControls* wellControls() const;

        void setVFPProperties(const VFPProperties<VFPInjProperties,VFPProdProperties>* vfp_properties_arg);

        virtual void init(const PhaseUsage* phase_usage_arg,
                          const std::vector<double>& depth_arg,
                          const double gravity_arg,
                          const int num_cells);

        virtual void initPrimaryVariablesEvaluation() const = 0;

        virtual ConvergenceReport getWellConvergence(const std::vector<double>& B_avg) const = 0;

        virtual void solveEqAndUpdateWellState(WellState& well_state) = 0;

        virtual void assembleWellEq(const Simulator& ebosSimulator,
                                    const double dt,
                                    WellState& well_state,
                                    Opm::DeferredLogger& deferred_logger
                                    ) = 0;

        void updateWellTestState(const WellState& well_state,
                                 const double& simulationTime,
                                 const bool& writeMessageToOPMLog,
                                 WellTestState& wellTestState,
                                 Opm::DeferredLogger& deferred_logger) const;

        void setWellEfficiencyFactor(const double efficiency_factor);

        void computeRepRadiusPerfLength(const Grid& grid, const std::vector<int>& cartesian_to_compressed);

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
                                           std::vector<double>& well_potentials,
                                           Opm::DeferredLogger& deferred_logger) = 0;

        virtual void updateWellStateWithTarget(const Simulator& ebos_simulator,
                                               WellState& well_state,
                                               Opm::DeferredLogger& deferred_logger) const = 0;

        void updateWellControl(/* const */ Simulator& ebos_simulator,
                               WellState& well_state,
                               Opm::DeferredLogger& deferred_logger) /* const */;

        virtual void updatePrimaryVariables(const WellState& well_state) const = 0;

        virtual void calculateExplicitQuantities(const Simulator& ebosSimulator,
                                                 const WellState& well_state) = 0; // should be const?

        /// \brief Wether the Jacobian will also have well contributions in it.
        virtual bool jacobianContainsWellContributions() const
        {
            return false;
        }

        // updating the voidage rates in well_state when requested
        void calculateReservoirRates(WellState& well_state) const;

        // Add well contributions to matrix
        virtual void addWellContributions(Mat&) const
        {}

        void addCellRates(RateVector& rates, int cellIdx) const;

        Scalar volumetricSurfaceRateForConnection(int cellIdx, int phaseIdx) const;


        template <class EvalWell>
        Eval restrictEval(const EvalWell& in) const
        {
            Eval out = 0.0;
            out.setValue(in.value());
            for(int eqIdx = 0; eqIdx < numEq;++eqIdx) {
                out.setDerivative(eqIdx, in.derivative(eqIdx));
            }
            return out;
        }

        void closeCompletions(WellTestState& wellTestState);

        const Well* wellEcl() const;

        // TODO: theoretically, it should be a const function
        // Simulator is not const is because that assembleWellEq is non-const Simulator
        void wellTesting(Simulator& simulator, const std::vector<double>& B_avg,
                         const double simulation_time, const int report_step,
                         const WellTestConfig::Reason testing_reason,
                         /* const */ WellState& well_state, WellTestState& welltest_state,
                         Opm::DeferredLogger& deferred_logger);

        void updatePerforatedCell(std::vector<bool>& is_cell_perforated);

        virtual void checkWellOperability(const Simulator& ebos_simulator, const WellState& well_state, Opm::DeferredLogger& deferred_logger) = 0;

        // whether the well is operable
        bool isOperable() const;

        /// Returns true if the well has one or more THP limits/constraints.
        bool wellHasTHPConstraints() const;

        /// Returns true if the well is currently in prediction mode (i.e. not history mode).
        bool underPredictionMode() const;

        // update perforation water throughput based on solved water rate
        virtual void updateWaterThroughput(const double dt, WellState& well_state) const = 0;

    protected:

        // to indicate a invalid completion
        static const int INVALIDCOMPLETION = INT_MAX;

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

        const VFPProperties<VFPInjProperties,VFPProdProperties>* vfp_properties_;

        double gravity_;

        // For the conversion between the surface volume rate and resrevoir voidage rate
        const RateConverterType& rateConverter_;

        // The pvt region of the well. We assume
        // We assume a well to not penetrate more than one pvt region.
        const int pvtRegionIdx_;

        const int num_components_;

        std::vector<RateVector> connectionRates_;

        const PhaseUsage& phaseUsage() const;

        int flowPhaseToEbosCompIdx( const int phaseIdx ) const;

        int ebosCompIdxToFlowCompIdx( const unsigned compIdx ) const;

        double wsolvent() const;

        double wpolymer() const;

        bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                                 const WellState& well_state,
                                 Opm::DeferredLogger& deferred_logger) const;

        double getTHPConstraint() const;

        int getTHPControlIndex() const;

        // Component fractions for each phase for the well
        const std::vector<double>& compFrac() const;

        double mostStrictBhpFromBhpLimits() const;

        // a tuple type for ratio limit check.
        // first value indicates whether ratio limit is violated, when the ratio limit is not violated, the following two
        // values should not be used.
        // second value indicates the index of the worst-offending completion.
        // the last value indicates the extent of the violation for the worst-offending completion, which is defined by
        // the ratio of the actual value to the value of the violated limit.
        using RatioCheckTuple = std::tuple<bool, int, double>;

        RatioCheckTuple checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                                              const WellState& well_state) const;

        RatioCheckTuple checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                                             const WellState& well_state,
                                             Opm::DeferredLogger& deferred_logger) const;

        double scalingFactor(const int comp_idx) const;

        // whether a well is specified with a non-zero and valid VFP table number
        bool isVFPActive() const;

        struct OperabilityStatus;

        OperabilityStatus operability_status_;

        void wellTestingEconomic(Simulator& simulator, const std::vector<double>& B_avg,
                                 const double simulation_time, const int report_step,
                                 const WellState& well_state, WellTestState& welltest_state, Opm::DeferredLogger& deferred_logger);

        virtual void wellTestingPhysical(Simulator& simulator, const std::vector<double>& B_avg,
                                 const double simulation_time, const int report_step,
                                         WellState& well_state, WellTestState& welltest_state, Opm::DeferredLogger& deferred_logger) = 0;

        void updateWellTestStateEconomic(const WellState& well_state,
                                         const double simulation_time,
                                         const bool write_message_to_opmlog,
                                         WellTestState& well_test_state,
                                         Opm::DeferredLogger& deferred_logger) const;

        void updateWellTestStatePhysical(const WellState& well_state,
                                         const double simulation_time,
                                         const bool write_message_to_opmlog,
                                         WellTestState& well_test_state,
                                         Opm::DeferredLogger& deferred_logger) const;

        void  solveWellForTesting(Simulator& ebosSimulator, WellState& well_state,
                                  const std::vector<double>& B_avg,
                                  Opm::DeferredLogger& deferred_logger);

        bool solveWellEqUntilConverged(Simulator& ebosSimulator,
                                       const std::vector<double>& B_avg,
                                       WellState& well_state,
                                       Opm::DeferredLogger& deferred_logger);

        void scaleProductivityIndex(const int perfIdx, double& productivity_index, Opm::DeferredLogger& deferred_logger);

        // count the number of times an output log message is created in the productivity
        // index calculations
        int well_productivity_index_logger_counter_;


    };




    // definition of the struct OperabilityStatus
    template<typename TypeTag>
    struct
    WellInterface<TypeTag>::
    OperabilityStatus {
        bool isOperable() const {
            if (!operable_under_only_bhp_limit) {
                return false;
            } else {
                return ( (isOperableUnderBHPLimit() || isOperableUnderTHPLimit()) );
            }
        }

        bool isOperableUnderBHPLimit() const {
            return operable_under_only_bhp_limit && obey_thp_limit_under_bhp_limit;
        }

        bool isOperableUnderTHPLimit() const {
            return can_obtain_bhp_with_thp_limit && obey_bhp_limit_with_thp_limit;
        }

        void reset() {
            operable_under_only_bhp_limit = true;
            obey_thp_limit_under_bhp_limit = true;
            can_obtain_bhp_with_thp_limit = true;
            obey_bhp_limit_with_thp_limit = true;
        }

        // whether the well can be operated under bhp limit
        // without considering other limits.
        // if it is false, then the well is not operable for sure.
        bool operable_under_only_bhp_limit = true;
        // if the well can be operated under bhp limit, will it obey(not violate)
        // the thp limit when operated under bhp limit
        bool obey_thp_limit_under_bhp_limit = true;
        // whether the well operate under the thp limit only
        bool can_obtain_bhp_with_thp_limit = true;
        // whether the well obey bhp limit when operated under thp limit
        bool obey_bhp_limit_with_thp_limit = true;

    };
    const std::string modestring[4] = { "BHP", "THP", "RESERVOIR_RATE", "SURFACE_RATE" };

}

#include "WellInterface_impl.hpp"

#endif // OPM_WELLINTERFACE_HEADER_INCLUDED
