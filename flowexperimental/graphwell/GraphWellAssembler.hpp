/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_GRAPHWELL_ASSEMBLER_HEADER_INCLUDED
#define OPM_GRAPHWELL_ASSEMBLER_HEADER_INCLUDED

#include <flowexperimental/graphwell/GraphWellEquations.hpp>
#include <flowexperimental/graphwell/GraphWellFluidProperties.hpp>
#include <flowexperimental/graphwell/GraphWellHelpers.hpp>
#include <flowexperimental/graphwell/GraphWellPrimaryVariables.hpp>
#include <flowexperimental/graphwell/GraphWellTopology.hpp>

#include <array>
#include <vector>

namespace Opm {

/// \brief Face-based AD assembly for the GraphWell model.
///
/// The assembly is two transparent loops (segments for accumulation/perforations,
/// connections for fluxes and momentum) plus a control equation on the surface
/// connection. Derivatives are scattered through two generic helpers
/// (scatterMass / scatterConn) that map the self/other/flux Eval slots to matrix
/// (row, column) pairs -- the only place layout knowledge lives.
template<class FluidSystem, int NumPhases, int NumResEq>
class GraphWellAssembler
{
public:
    static constexpr int NP = NumPhases;
    using PV = GraphWellPrimaryVariables<FluidSystem, NumPhases, NumResEq>;
    using FP = GraphWellFluidProperties<FluidSystem, NumPhases, NumResEq>;
    using Eqns = GraphWellEquations<typename PV::Scalar, NumPhases, NumResEq>;
    using Topo = GraphWellTopology<typename PV::Scalar>;
    using Scalar = typename PV::Scalar;
    using Eval = typename PV::Eval;
    using Role = typename PV::Role;

    struct PerforationInput
    {
        Eval pressure{Scalar{0}};                 //!< cell pressure (reservoir-slot derivs)
        std::array<Eval, NP> mob{};               //!< phase/component mobilities
        std::array<Eval, NP> b{};                 //!< inverse formation-volume factors
        Eval rs{Scalar{0}};
        Eval rv{Scalar{0}};
        Scalar Tw{0};                             //!< well index * transmissibility mult
        Scalar cell_perf_press_diff{0};           //!< cell->perf depth pressure correction
    };

    struct Controls
    {
        enum class Type { BHP, SurfaceRate };
        Type type{Type::BHP};
        Scalar target{0};
        int rate_comp{-1};   //!< active component for SurfaceRate (-1 => total)
        bool producer{true};
    };

    //! Pressure-drop model options, mirroring the deck's WELSEGS compPressureDrop flag
    //! (H__ / HF- / HFA) and the use-average-density-MS-wells option.
    struct PdropOptions
    {
        bool friction{true};        //!< include friction (HF- and HFA)
        bool acceleration{true};    //!< include acceleration (HFA only)
        bool average_density{false};//!< hydrostatic uses 0.5*(down+up) density
    };

    void assemble(const Topo& topo,
                  const PV& pv,
                  const FP& fp,
                  const std::vector<PerforationInput>& perfs,
                  const Controls& controls,
                  Scalar dt,
                  Scalar gravity,
                  bool allow_crossflow,
                  const std::vector<std::array<Scalar, NP>>& old_surface_volumes,
                  Eqns& eqns,
                  bool make_solver = true,
                  bool assemble_control = true,
                  Scalar regularization_factor = 1.0,
                  const PdropOptions& pdrop = {}) const
    {
        eqns.clear();
        assembleAccumulation(topo, pv, fp, dt, old_surface_volumes, regularization_factor, eqns);
        assemblePerforations(topo, pv, fp, perfs, controls.producer, allow_crossflow, gravity, eqns);
        for (int c = 0; c < topo.numConnections(); ++c) {
            assembleConnectionFlux(topo, pv, fp, c, eqns);
            if (c == topo.surfaceConnection()) {
                // The control equation may instead be supplied externally (e.g. borrowed
                // from the production MultisegmentWell to reuse its control/THP/group logic).
                if (assemble_control)
                    assembleControl(topo, pv, fp, controls, c, eqns);
            } else {
                assemblePressureDrop(topo, pv, fp, gravity, c, pdrop, eqns);
            }
        }
        if (make_solver)
            eqns.createSolver();
    }

private:
    static int waterComp() { return FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
        ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx) : -1; }
    static int oilComp() { return FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
        ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::oilCompIdx) : -1; }
    static int gasComp() { return FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)
        ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::gasCompIdx) : -1; }

    // -- generic scatter helpers ---------------------------------------------

    //! Scatter a mass-equation term (segment row, given component).
    void scatterMass(Eqns& eqns, int row_seg, int comp,
                     int self_col, int other_col, int conn_col,
                     const Eval& term) const
    {
        eqns.resSeg()[row_seg][comp] += term.value();
        auto& dss = eqns.Dss();
        for (int k = 0; k < NP; ++k)
            dss[row_seg][self_col][comp][k] += term.derivative(PV::RES + k);
        if (other_col >= 0) {
            for (int k = 0; k < NP; ++k)
                dss[row_seg][other_col][comp][k] += term.derivative(PV::OTH + k);
        }
        if (conn_col >= 0)
            eqns.Dsc()[row_seg][conn_col][comp][0] += term.derivative(PV::QIDX);
    }

    //! Scatter a connection-equation term (momentum or control).
    void scatterConn(Eqns& eqns, int conn_row,
                     int self_col, int other_col,
                     const Eval& term) const
    {
        eqns.resConn()[conn_row][0] += term.value();
        auto& dcs = eqns.Dcs();
        for (int k = 0; k < NP; ++k)
            dcs[conn_row][self_col][0][k] += term.derivative(PV::RES + k);
        if (other_col >= 0) {
            for (int k = 0; k < NP; ++k)
                dcs[conn_row][other_col][0][k] += term.derivative(PV::OTH + k);
        }
        eqns.Dcc()[conn_row][conn_row][0][0] += term.derivative(PV::QIDX);
    }

    // -- terms ----------------------------------------------------------------

    void assembleAccumulation(const Topo& topo, const PV& pv, const FP& fp, Scalar dt,
                              const std::vector<std::array<Scalar, NP>>& old_sv,
                              Scalar regularization_factor,
                              Eqns& eqns) const
    {
        for (int s = 0; s < topo.numSegments(); ++s) {
            const Scalar volume = topo.segment(s).volume;
            if (volume <= Scalar{0})
                continue;
            const Eval surface_volume = Eval{volume} / fp.volRatio(s);
            for (int comp = 0; comp < NP; ++comp) {
                const Eval frac = pv.surfaceVolumeFraction(s, comp, Role::Self);
                const Eval accum = regularization_factor
                    * (surface_volume * frac - old_sv[s][comp]) / dt;
                scatterMass(eqns, s, comp, /*self*/ s, /*other*/ -1, /*conn*/ -1, accum);
            }
        }
    }

    void assembleConnectionFlux(const Topo& topo, const PV& pv, const FP& /*fp*/, int c,
                                Eqns& eqns) const
    {
        const auto& conn = topo.connection(c);
        const Scalar Qval = pv.connRateValue(c);
        const int u = topo.upwindSegment(c, Qval);   // upwind segment (a real segment)
        const Eval Q = pv.connRate(c);
        for (int comp = 0; comp < NP; ++comp) {
            const Eval F = pv.volumeFractionScaled(u, comp, Role::Self);
            const Eval q = Q * F;   // down->up component rate (>0 for production)
            // residual_s += outflow(s,c)*q ; q leaves the down segment (+) and
            // enters the up segment (-). Derivs: self block = upwind u, conn = c.
            if (conn.down != Topo::surface_node)
                scatterMass(eqns, conn.down, comp, /*self*/ u, -1, /*conn*/ c, q);
            if (conn.up != Topo::surface_node) {
                const Eval qn = -q;
                scatterMass(eqns, conn.up, comp, /*self*/ u, -1, /*conn*/ c, qn);
            }
        }
    }

    //! r_c = p_down - p_up - dp_hydro - dp_friction - dp_acceleration
    void assemblePressureDrop(const Topo& topo, const PV& pv, const FP& fp, Scalar gravity,
                              int c, const PdropOptions& pdrop, Eqns& eqns) const
    {
        const auto& conn = topo.connection(c);
        const int down = conn.down;
        const int up = conn.up;
        const Scalar Qval = pv.connRateValue(c);
        const int u = topo.upwindSegment(c, Qval);
        const Role uprole = (u == down) ? Role::Self : Role::Other;

        // p_down (Self) - p_up (Other)
        Eval r = pv.segPressure(down, Role::Self) - pv.segPressure(up, Role::Other);

        // hydrostatic: density of the down (deeper) segment, or the average of the two
        const Eval rho_down = fp.density(down, Role::Self);
        Eval rho_hydro = rho_down;
        if (pdrop.average_density)
            rho_hydro = Scalar{0.5} * (rho_down + fp.density(up, Role::Other));
        r -= rho_hydro * gravity * conn.depth_diff;

        if (pdrop.friction || pdrop.acceleration) {
            // mass rate = Q * sum_comp F(upwind) * surf_dens
            const Eval Q = pv.connRate(c);
            Eval mass_rate(Scalar{0});
            for (int comp = 0; comp < NP; ++comp) {
                const Eval F = pv.volumeFractionScaled(u, comp, uprole);
                mass_rate += Q * F * fp.surfaceDensity(comp);
            }
            const Eval rho_up = fp.density(u, uprole);
            const Eval mu_up = fp.viscosity(u, uprole);

            if (pdrop.friction) {
                // friction: magnitude * sign(Q)
                const Scalar signQ = (Qval >= Scalar{0}) ? Scalar{1} : Scalar{-1};
                const Eval fric = graphwellhelpers::frictionPressureLoss(
                    conn.length, conn.diameter, conn.area, conn.roughness, rho_up, mass_rate, mu_up);
                r -= signQ * fric;
            }
            if (pdrop.acceleration) {
                // acceleration: velocity-head difference across the connection (local form)
                const Eval rho_a = fp.density(down, Role::Self);
                const Eval rho_b = fp.density(up, Role::Other);
                const Eval accel = graphwellhelpers::velocityHead(conn.area, mass_rate, rho_b)
                                 - graphwellhelpers::velocityHead(conn.area, mass_rate, rho_a);
                r -= accel;
            }
        }

        scatterConn(eqns, c, /*self*/ down, /*other*/ up, r);
    }

    void assembleControl(const Topo& topo, const PV& pv, const FP& /*fp*/,
                         const Controls& ctrl, int c, Eqns& eqns) const
    {
        const int top = topo.connection(c).down;
        Eval r(Scalar{0});
        if (ctrl.type == Controls::Type::BHP) {
            r = pv.segPressure(top, Role::Self) - ctrl.target;
        } else {
            const Eval Q = pv.connRate(c);
            if (ctrl.rate_comp < 0) {
                r = Q - ctrl.target;
            } else {
                const Eval F = pv.surfaceVolumeFraction(top, ctrl.rate_comp, Role::Self);
                r = Q * F - ctrl.target;
            }
        }
        scatterConn(eqns, c, /*self*/ top, /*other*/ -1, r);
    }

    void assemblePerforations(const Topo& topo, const PV& pv, const FP& fp,
                              const std::vector<PerforationInput>& perfs,
                              bool producer, bool allow_cf, Scalar gravity,
                              Eqns& eqns) const
    {
        auto& dss = eqns.Dss();
        auto& B = eqns.B();
        auto& C = eqns.C();
        for (int s = 0; s < topo.numSegments(); ++s) {
            const Eval p_seg = pv.segPressure(s, Role::Self);
            const Eval rho_seg = fp.density(s, Role::Self);
            for (int perf : topo.perforationsOfSegment(s)) {
                std::array<Eval, NP> cmix_s;
                for (int comp = 0; comp < NP; ++comp)
                    cmix_s[comp] = pv.surfaceVolumeFraction(s, comp, Role::Self);
                std::array<Eval, NP> cq_s;
                computePerfRate(perfs[perf], p_seg, rho_seg, cmix_s,
                                topo.perforationDepthDiff(perf), gravity,
                                producer, allow_cf, cq_s);
                for (int comp = 0; comp < NP; ++comp) {
                    eqns.resSeg()[s][comp] += cq_s[comp].value();
                    for (int k = 0; k < NP; ++k) {
                        dss[s][s][comp][k] += cq_s[comp].derivative(PV::RES + k);
                        // C stored transposed: [seg_var][comp]
                        C[s][perf][k][comp] -= cq_s[comp].derivative(PV::RES + k);
                    }
                    for (int k = 0; k < NumResEq; ++k)
                        B[s][perf][comp][k] += cq_s[comp].derivative(k);
                }
            }
        }
    }

public:
    //! Port of MultisegmentWell::computePerfRate (HFA scope, no perf_rates output).
    //! Public so the well can reuse it for well-state bookkeeping.
    void computePerfRate(const PerforationInput& in,
                         const Eval& p_seg, const Eval& rho_seg,
                         const std::array<Eval, NP>& cmix_s,
                         Scalar perf_depth_diff, Scalar gravity,
                         bool producer, bool allow_cf,
                         std::array<Eval, NP>& cq_s) const
    {
        for (int c = 0; c < NP; ++c)
            cq_s[c] = Eval{Scalar{0}};

        const Eval perf_seg_press_diff = gravity * rho_seg * perf_depth_diff;
        const Eval perf_press = p_seg + perf_seg_press_diff;
        const Eval cell_press_at_perf = in.pressure - in.cell_perf_press_diff;
        const Eval drawdown = cell_press_at_perf - perf_press;

        const int wc = waterComp();
        const int oc = oilComp();
        const int gc = gasComp();
        const bool og = (oc >= 0 && gc >= 0);

        if (drawdown.value() > Scalar{0}) {     // producing perforation
            if (!allow_cf && !producer)
                return;
            for (int comp = 0; comp < NP; ++comp) {
                const Eval cq_p = -in.Tw * (in.mob[comp] * drawdown);
                cq_s[comp] = in.b[comp] * cq_p;
            }
            if (og) {
                const Eval cq_s_oil = cq_s[oc];
                const Eval cq_s_gas = cq_s[gc];
                cq_s[gc] += in.rs * cq_s_oil;
                cq_s[oc] += in.rv * cq_s_gas;
            }
        } else {                                 // injecting perforation
            if (!allow_cf && producer)
                return;
            Eval total_mob = in.mob[0];
            for (int comp = 1; comp < NP; ++comp)
                total_mob += in.mob[comp];

            Eval volume_ratio(Scalar{0});
            if (wc >= 0)
                volume_ratio += cmix_s[wc] / in.b[wc];
            if (og) {
                const Eval d = Eval{Scalar{1}} - in.rv * in.rs;
                const Eval tmp_oil = (cmix_s[oc] - in.rv * cmix_s[gc]) / d;
                volume_ratio += tmp_oil / in.b[oc];
                const Eval tmp_gas = (cmix_s[gc] - in.rs * cmix_s[oc]) / d;
                volume_ratio += tmp_gas / in.b[gc];
            } else {
                if (oc >= 0)
                    volume_ratio += cmix_s[oc] / in.b[oc];
                if (gc >= 0)
                    volume_ratio += cmix_s[gc] / in.b[gc];
            }
            for (int comp = 0; comp < NP; ++comp) {
                const Eval cqt_i = -in.Tw * (total_mob * drawdown);
                const Eval cqt_is = cqt_i / volume_ratio;
                cq_s[comp] = cmix_s[comp] * cqt_is;
            }
        }
    }
};

} // namespace Opm

#endif // OPM_GRAPHWELL_ASSEMBLER_HEADER_INCLUDED
