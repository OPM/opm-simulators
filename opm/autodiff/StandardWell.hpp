/*
  Copyright 2017 SINTEF ICT, Applied Mathematics.
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


#ifndef OPM_STANDARDWELL_HEADER_INCLUDED
#define OPM_STANDARDWELL_HEADER_INCLUDED


#include <opm/autodiff/WellInterface.hpp>

namespace Opm
{

    template<typename TypeTag>
    class StandardWell: public WellInterface<TypeTag>
    {

    public:
        // TODO: some functions working with AD variables handles only with values (double) without
        // dealing with derivatives. It can be beneficial to make functions can work with either AD or scalar value.
        // And also, it can also be beneficial to make these functions hanle different types of AD variables.
        // TODO: several functions related to polymer and PLYSHLOG are not incorprated yet,
        // like the function wpolymer, setupCompressedToCartesian, computeRepRadiusPerfLength,
        // They are introduced though PR 1220 and will be included later.
        using typename WellInterface<TypeTag>::Simulator;
        using typename WellInterface<TypeTag>::WellState;
        using typename WellInterface<TypeTag>::IntensiveQuantities;
        using typename WellInterface<TypeTag>::FluidSystem;
        using typename WellInterface<TypeTag>::MaterialLaw;
        using typename WellInterface<TypeTag>::ModelParameters;
        using typename WellInterface<TypeTag>::BlackoilIndices;

        // the positions of the primary variables for StandardWell
        // there are three primary variables, the second and the third ones are F_w and F_g
        // the first one can be total rate (G_t) or bhp, based on the control
        enum WellVariablePositions {
            XvarWell = 0,
            WFrac = 1,
            GFrac = 2,
            SFrac = 3
        };

        using typename WellInterface<TypeTag>::VectorBlockType;
        using typename WellInterface<TypeTag>::MatrixBlockType;
        using typename WellInterface<TypeTag>::Mat;
        using typename WellInterface<TypeTag>::BVector;
        using typename WellInterface<TypeTag>::Eval;

        using WellInterface<TypeTag>::numEq;
        static const int numWellEq = GET_PROP_VALUE(TypeTag, EnablePolymer)? numEq-1 : numEq; // //numEq; //number of wellEq is only numEq for polymer
        // TODO: should these go to WellInterface?
        static const int contiSolventEqIdx = BlackoilIndices::contiSolventEqIdx;
        static const int contiPolymerEqIdx = BlackoilIndices::contiPolymerEqIdx;
        static const int solventSaturationIdx = BlackoilIndices::solventSaturationIdx;
        static const int polymerConcentrationIdx = BlackoilIndices::polymerConcentrationIdx;

        typedef DenseAd::Evaluation<double, /*size=*/numEq + numWellEq> EvalWell;

        StandardWell(const Well* well, const int time_step, const Wells* wells);

        /// the densities of the fluid  in each perforation
        virtual const std::vector<double>& perfDensities() const;
        virtual std::vector<double>& perfDensities();

        /// the pressure difference between different perforations
        virtual const std::vector<double>& perfPressureDiffs() const;
        virtual std::vector<double>& perfPressureDiffs();

        virtual void setWellVariables(const WellState& well_state);

        EvalWell wellVolumeFractionScaled(const int phase) const;

        EvalWell wellVolumeFraction(const int phase) const;

        EvalWell wellSurfaceVolumeFraction(const int phase) const;

        EvalWell extendEval(const Eval& in) const;

        // TODO: to check whether all the paramters are required
        void computePerfRate(const IntensiveQuantities& intQuants,
                             const std::vector<EvalWell>& mob_perfcells_dense,
                             const double Tw, const EvalWell& bhp, const double& cdp,
                             const bool& allow_cf, std::vector<EvalWell>& cq_s) const;

        virtual void assembleWellEq(Simulator& ebosSimulator,
                                    const double dt,
                                    WellState& well_state,
                                    bool only_wells);

        virtual bool crossFlowAllowed(const Simulator& ebosSimulator) const;

        void getMobility(const Simulator& ebosSimulator,
                         const int perf,
                         std::vector<EvalWell>& mob) const;

        // TODO: the parameters need to be optimized/adjusted
        virtual void init(const PhaseUsage* phase_usage_arg,
                          const std::vector<bool>* active_arg,
                          const VFPProperties* vfp_properties_arg,
                          const double gravity_arg,
                          const int num_cells);

        // Update the well_state based on solution
        void updateWellState(const BVector& dwells,
                             const BlackoilModelParameters& param,
                             WellState& well_state) const;

        // TODO: later will check wheter we need current
        virtual void updateWellStateWithTarget(const int current,
                                               WellState& xw) const;

        // TODO: this should go to the WellInterface, while updateWellStateWithTarget
        // will need touch different types of well_state, we will see.
        virtual void updateWellControl(WellState& xw) const;

        virtual bool getWellConvergence(Simulator& ebosSimulator,
                                        const std::vector<double>& B_avg,
                                        const ModelParameters& param) const;

        virtual void computeAccumWell();

        virtual void computeWellConnectionPressures(const Simulator& ebosSimulator,
                                                    const WellState& xw);

        // Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const;
        // r = r - C D^-1 Rw
        virtual void apply(BVector& r) const;

        // using the solution x to recover the solution xw for wells and applying
        // xw to update Well State
        virtual void applySolutionWellState(const BVector& x, const ModelParameters& param,
                                            WellState& well_state) const;

        virtual void computeWellPotentials(const Simulator& ebosSimulator,
                                           const WellState& well_state,
                                           std::vector<double>& well_potentials) const;

        using WellInterface<TypeTag>::has_solvent;
        using WellInterface<TypeTag>::has_polymer;

        using WellInterface<TypeTag>::phaseUsage;
        using WellInterface<TypeTag>::active;
        using WellInterface<TypeTag>::numberOfPerforations;
        using WellInterface<TypeTag>::wellCells;
        using WellInterface<TypeTag>::saturationTableNumber;
        using WellInterface<TypeTag>::indexOfWell;
        using WellInterface<TypeTag>::name;
        using WellInterface<TypeTag>::wellType;
        using WellInterface<TypeTag>::wellControls;
        using WellInterface<TypeTag>::compFrac;
        using WellInterface<TypeTag>::numberOfPhases;
        using WellInterface<TypeTag>::perfDepth;
        using WellInterface<TypeTag>::flowToEbosPvIdx;
        using WellInterface<TypeTag>::flowPhaseToEbosPhaseIdx;
        using WellInterface<TypeTag>::flowPhaseToEbosCompIdx;
        using WellInterface<TypeTag>::numComponents;
        using WellInterface<TypeTag>::numPhases;
        using WellInterface<TypeTag>::wellIndex;
        using WellInterface<TypeTag>::wsolvent;
        using WellInterface<TypeTag>::wpolymer;

    protected:

        // TODO: maybe this function can go to some helper file.
        void localInvert(Mat& istlA) const;

        // xw = inv(D)*(rw - C*x)
        void recoverSolutionWell(const BVector& x, BVector& xw) const;

        // TODO: decide wether to use member function to refer to private member later
        using WellInterface<TypeTag>::vfp_properties_;
        using WellInterface<TypeTag>::gravity_;
        using WellInterface<TypeTag>::well_efficiency_factor_;
        using WellInterface<TypeTag>::phase_usage_;
        using WellInterface<TypeTag>::first_perf_;
        using WellInterface<TypeTag>::ref_depth_;
        using WellInterface<TypeTag>::perf_depth_;
        using WellInterface<TypeTag>::allow_cf_;

        // densities of the fluid in each perforation
        std::vector<double> perf_densities_;
        // pressure drop between different perforations
        std::vector<double> perf_pressure_diffs_;

        // TODO: probably, they should be moved to the WellInterface, when
        // we decide the template paramters.
        // two off-diagonal matrices
        Mat duneB_;
        Mat duneC_;
        // diagonal matrix for the well
        Mat invDuneD_;

        // several vector used in the matrix calculation
        mutable BVector Bx_;
        mutable BVector invDrw_;
        mutable BVector scaleAddRes_;

        BVector resWell_;

        std::vector<EvalWell> well_variables_;
        std::vector<double> F0_;

        // TODO: this function should be moved to the base class.
        // while it faces chanllenges for MSWell later, since the calculation of bhp
        // based on THP is never implemented for MSWell yet.
        EvalWell getBhp() const;

        // TODO: it is also possible to be moved to the base class.
        EvalWell getQs(const int comp_idx) const;

        // calculate the properties for the well connections
        // to calulate the pressure difference between well connections.
        void computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                         const WellState& xw,
                                                         std::vector<double>& b_perf,
                                                         std::vector<double>& rsmax_perf,
                                                         std::vector<double>& rvmax_perf,
                                                         std::vector<double>& surf_dens_perf) const;

        // TODO: not total sure whether it is a good idea to put here
        // the major reason to put here is to avoid the usage of Wells struct
        void computeConnectionDensities(const std::vector<double>& perfComponentRates,
                                        const std::vector<double>& b_perf,
                                        const std::vector<double>& rsmax_perf,
                                        const std::vector<double>& rvmax_perf,
                                        const std::vector<double>& surf_dens_perf);

        void computeConnectionPressureDelta();

        void computeWellConnectionDensitesPressures(const WellState& xw,
                                                    const std::vector<double>& b_perf,
                                                    const std::vector<double>& rsmax_perf,
                                                    const std::vector<double>& rvmax_perf,
                                                    const std::vector<double>& surf_dens_perf);

        virtual void wellEqIteration(Simulator& ebosSimulator,
                                     const ModelParameters& param,
                                     WellState& well_state);

        using WellInterface<TypeTag>::wellHasTHPConstraints;
        using WellInterface<TypeTag>::mostStrictBhpFromBhpLimits;

        // TODO: maybe we should provide a light version of computeWellFlux, which does not include the
        // calculation of the derivatives
        void computeWellRatesWithBhp(const Simulator& ebosSimulator,
                                     const EvalWell& bhp,
                                     std::vector<double>& well_flux) const;

        std::vector<double> computeWellPotentialWithTHP(const Simulator& ebosSimulator,
                                                        const double initial_bhp, // bhp from BHP constraints
                                                        const std::vector<double>& initial_potential) const;
    };

}

#include "StandardWell_impl.hpp"

#endif // OPM_STANDARDWELL_HEADER_INCLUDED
