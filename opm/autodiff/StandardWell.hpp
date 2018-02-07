/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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


#ifndef OPM_STANDARDWELL_HEADER_INCLUDED
#define OPM_STANDARDWELL_HEADER_INCLUDED


#include <opm/autodiff/WellInterface.hpp>
#include <opm/autodiff/ISTLSolver.hpp>
#include <opm/autodiff/RateConverter.hpp>

namespace Opm
{

    template<typename TypeTag>
    class StandardWell: public WellInterface<TypeTag>
    {

    public:
        typedef WellInterface<TypeTag> Base;
        // TODO: some functions working with AD variables handles only with values (double) without
        // dealing with derivatives. It can be beneficial to make functions can work with either AD or scalar value.
        // And also, it can also be beneficial to make these functions hanle different types of AD variables.
        using typename Base::Simulator;
        using typename Base::WellState;
        using typename Base::IntensiveQuantities;
        using typename Base::FluidSystem;
        using typename Base::MaterialLaw;
        using typename Base::ModelParameters;
        using typename Base::Indices;
        using typename Base::PolymerModule;
        using typename Base::RateConverterType;

        using Base::numEq;

        // the positions of the primary variables for StandardWell
        // there are three primary variables, the second and the third ones are F_w and F_g
        // the first one can be total rate (G_t) or bhp, based on the control

        static const bool gasoil = numEq == 2 && (Indices::compositionSwitchIdx >= 0);
        static const int XvarWell = 0;
        static const int WFrac = gasoil? -1000: 1;
        static const int GFrac = gasoil? 1: 2;
        static const int SFrac = 3;

        using typename Base::Scalar;
        using typename Base::ConvergenceReport;


        using Base::has_solvent;
        using Base::has_polymer;
        using Base::name;
        using Base::Water;
        using Base::Oil;
        using Base::Gas;

        // TODO: with flow_ebos，for a 2P deck, // TODO: for the 2p deck, numEq will be 3, a dummy phase is already added from the reservoir side.
        // it will cause problem here without processing the dummy phase.
        static const int numWellEq = GET_PROP_VALUE(TypeTag, EnablePolymer)? numEq-1 : numEq; // number of wellEq is only numEq - 1 for polymer
        using typename Base::Mat;
        using typename Base::BVector;
        using typename Base::Eval;

        // sparsity pattern for the matrices
        //[A C^T    [x       =  [ res
        // B  D ]   x_well]      res_well]

        // the vector type for the res_well and x_well
        typedef Dune::FieldVector<Scalar, numWellEq> VectorBlockWellType;
        typedef Dune::BlockVector<VectorBlockWellType> BVectorWell;

#if  DUNE_VERSION_NEWER_REV(DUNE_ISTL, 2 , 5, 1)
        // 3x3 matrix block inversion was unstable from at least 2.3 until and
        // including 2.5.0
        // the matrix type for the diagonal matrix D
        typedef Dune::FieldMatrix<Scalar, numWellEq, numWellEq > DiagMatrixBlockWellType;
#else
        // the matrix type for the diagonal matrix D
        typedef Dune::MatrixBlock<Scalar, numWellEq, numWellEq > DiagMatrixBlockWellType;
#endif

        typedef Dune::BCRSMatrix <DiagMatrixBlockWellType> DiagMatWell;

        // the matrix type for the non-diagonal matrix B and C^T
        typedef Dune::FieldMatrix<Scalar, numWellEq, numEq>  OffDiagMatrixBlockWellType;
        typedef Dune::BCRSMatrix<OffDiagMatrixBlockWellType> OffDiagMatWell;

        typedef DenseAd::Evaluation<double, /*size=*/numEq + numWellEq> EvalWell;

        using Base::contiSolventEqIdx;
        using Base::contiPolymerEqIdx;


        StandardWell(const Well* well, const int time_step, const Wells* wells,
                     const ModelParameters& param,
                     const RateConverterType& rate_converter,
                     const int pvtRegionIdx,
                     const int num_components);

        virtual void init(const PhaseUsage* phase_usage_arg,
                          const std::vector<double>& depth_arg,
                          const double gravity_arg,
                          const size_t num_cells);


        virtual void initPrimaryVariablesEvaluation() const;

        virtual void assembleWellEq(Simulator& ebosSimulator,
                                    const double dt,
                                    WellState& well_state,
                                    bool only_wells);

        /// updating the well state based the control mode specified with current
        // TODO: later will check wheter we need current
        virtual void updateWellStateWithTarget(WellState& well_state) const;

        /// check whether the well equations get converged for this well
        virtual ConvergenceReport getWellConvergence(const std::vector<double>& B_avg) const;

        /// Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const;
        /// r = r - C D^-1 Rw
        virtual void apply(BVector& r) const;

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        virtual void recoverWellSolutionAndUpdateWellState(const BVector& x,
                                                           WellState& well_state) const;

        /// computing the well potentials for group control
        virtual void computeWellPotentials(const Simulator& ebosSimulator,
                                           const WellState& well_state,
                                           std::vector<double>& well_potentials) /* const */;

        virtual void updatePrimaryVariables(const WellState& well_state) const;

        virtual void solveEqAndUpdateWellState(WellState& well_state);

        virtual void calculateExplicitQuantities(const Simulator& ebosSimulator,
                                                 const WellState& well_state); // should be const?
    protected:

        // protected functions from the Base class
        using Base::getAllowCrossFlow;
        using Base::phaseUsage;
        using Base::flowPhaseToEbosCompIdx;
        using Base::ebosCompIdxToFlowCompIdx;
        using Base::wsolvent;
        using Base::wpolymer;
        using Base::wellHasTHPConstraints;
        using Base::mostStrictBhpFromBhpLimits;
        using Base::scalingFactor;

        // protected member variables from the Base class
        using Base::vfp_properties_;
        using Base::gravity_;
        using Base::param_;
        using Base::well_efficiency_factor_;
        using Base::first_perf_;
        using Base::ref_depth_;
        using Base::perf_depth_;
        using Base::well_cells_;
        using Base::number_of_perforations_;
        using Base::number_of_phases_;
        using Base::saturation_table_number_;
        using Base::comp_frac_;
        using Base::well_index_;
        using Base::index_of_well_;
        using Base::well_controls_;
        using Base::well_type_;
        using Base::num_components_;

        using Base::perf_rep_radius_;
        using Base::perf_length_;
        using Base::bore_diameters_;

        // densities of the fluid in each perforation
        std::vector<double> perf_densities_;
        // pressure drop between different perforations
        std::vector<double> perf_pressure_diffs_;

        // residuals of the well equations
        BVectorWell resWell_;

        // two off-diagonal matrices
        OffDiagMatWell duneB_;
        OffDiagMatWell duneC_;
        // diagonal matrix for the well
        DiagMatWell invDuneD_;

        // several vector used in the matrix calculation
        mutable BVectorWell Bx_;
        mutable BVectorWell invDrw_;

        // the values for the primary varibles
        // based on different solutioin strategies, the wells can have different primary variables
        mutable std::vector<double> primary_variables_;

        // the Evaluation for the well primary variables, which contain derivativles and are used in AD calculation
        mutable std::vector<EvalWell> primary_variables_evaluation_;

        // the saturations in the well bore under surface conditions at the beginning of the time step
        std::vector<double> F0_;

        // TODO: this function should be moved to the base class.
        // while it faces chanllenges for MSWell later, since the calculation of bhp
        // based on THP is never implemented for MSWell yet.
        EvalWell getBhp() const;

        // TODO: it is also possible to be moved to the base class.
        EvalWell getQs(const int comp_idx) const;

        EvalWell wellVolumeFractionScaled(const int phase) const;

        EvalWell wellVolumeFraction(const unsigned compIdx) const;

        EvalWell wellSurfaceVolumeFraction(const int phase) const;

        EvalWell extendEval(const Eval& in) const;

        bool crossFlowAllowed(const Simulator& ebosSimulator) const;

        // xw = inv(D)*(rw - C*x)
        void recoverSolutionWell(const BVector& x, BVectorWell& xw) const;

        // updating the well_state based on well solution dwells
        void updateWellState(const BVectorWell& dwells,
                             WellState& well_state) const;

        // calculate the properties for the well connections
        // to calulate the pressure difference between well connections.
        void computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                         const WellState& well_state,
                                                         std::vector<double>& b_perf,
                                                         std::vector<double>& rsmax_perf,
                                                         std::vector<double>& rvmax_perf,
                                                         std::vector<double>& surf_dens_perf) const;

        // TODO: not total sure whether it is a good idea to put this function here
        // the major reason to put here is to avoid the usage of Wells struct
        void computeConnectionDensities(const std::vector<double>& perfComponentRates,
                                        const std::vector<double>& b_perf,
                                        const std::vector<double>& rsmax_perf,
                                        const std::vector<double>& rvmax_perf,
                                        const std::vector<double>& surf_dens_perf);

        void computeConnectionPressureDelta();

        void computeWellConnectionDensitesPressures(const WellState& well_state,
                                                    const std::vector<double>& b_perf,
                                                    const std::vector<double>& rsmax_perf,
                                                    const std::vector<double>& rvmax_perf,
                                                    const std::vector<double>& surf_dens_perf);

        // computing the accumulation term for later use in well mass equations
        void computeAccumWell();

        void computeWellConnectionPressures(const Simulator& ebosSimulator,
                                                    const WellState& well_state);

        // TODO: to check whether all the paramters are required
        void computePerfRate(const IntensiveQuantities& intQuants,
                             const std::vector<EvalWell>& mob_perfcells_dense,
                             const double Tw, const EvalWell& bhp, const double& cdp,
                             const bool& allow_cf, std::vector<EvalWell>& cq_s) const;

        // TODO: maybe we should provide a light version of computePerfRate, which does not include the
        // calculation of the derivatives
        void computeWellRatesWithBhp(const Simulator& ebosSimulator,
                                     const EvalWell& bhp,
                                     std::vector<double>& well_flux) const;

        std::vector<double> computeWellPotentialWithTHP(const Simulator& ebosSimulator,
                                                        const double initial_bhp, // bhp from BHP constraints
                                                        const std::vector<double>& initial_potential) const;

        template <class ValueType>
        ValueType calculateBhpFromThp(const std::vector<ValueType>& rates, const int control_index) const;

        double calculateThpFromBhp(const std::vector<double>& rates, const int control_index, const double bhp) const;

        // get the mobility for specific perforation
        void getMobility(const Simulator& ebosSimulator,
                         const int perf,
                         std::vector<EvalWell>& mob) const;

        void updateWaterMobilityWithPolymer(const Simulator& ebos_simulator,
                                            const int perf,
                                            std::vector<EvalWell>& mob_water) const;
    };

}

#include "StandardWell_impl.hpp"

#endif // OPM_STANDARDWELL_HEADER_INCLUDED
