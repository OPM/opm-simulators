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
#include <cmath>
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
        Scalar perf_depth_diff{0};                //!< perfDepth - segmentDepth (well-side)
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

    //! Per-connection flow-control device (WSEGVALV / WSEGSICD / WSEGAICD). When present
    //! the device pressure drop replaces the normal pipe friction on that connection
    //! (hydrostatic and acceleration are still applied). A shut device forces zero flux.
    struct ConnDevice
    {
        enum class Kind { None, Valve, SpiralICD, AutoICD };
        Kind kind{Kind::None};
        bool shut{false};
        // valve (WSEGVALV): pipe-friction geometry + constriction
        Scalar length{0};
        Scalar diameter{0};
        Scalar area{0};
        Scalar roughness{0};
        Scalar area_con{0};   //!< constriction cross-area (UDA-resolved at assembly time)
        Scalar cv{0};         //!< constriction flow coefficient
        // spiral ICD (WSEGSICD) + shared by AutoICD
        Scalar scaling_factor{0};
        Scalar visc_cali{0};
        Scalar dens_cali{0};
        Scalar strength{0};
        Scalar width_transition{0};
        Scalar critical_value{0};
        Scalar max_visco_ratio{0};
        // autonomous ICD (WSEGAICD)
        Scalar flow_rate_exp{0};
        Scalar visc_exp{0};
        Scalar dens_exp{0};
        Scalar wat_visc_exp{0};
        Scalar oil_visc_exp{0};
        Scalar gas_visc_exp{0};
        Scalar wat_dens_exp{0};
        Scalar oil_dens_exp{0};
        Scalar gas_dens_exp{0};
        Scalar unit_volume_rate{1};   //!< to_si(geometric_volume_rate, 1)
    };

    //! Per-connection pressure-drop components (signed, for well-state reporting / SPRD).
    struct PdropReport
    {
        Scalar hydrostatic{0};
        Scalar friction{0};      //!< friction or device drop (MSW reports both in this slot)
        Scalar acceleration{0};
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
                  const PdropOptions& pdrop = {},
                  Scalar efficiency = 1.0,
                  const std::vector<ConnDevice>& devices = {},
                  std::vector<PdropReport>* report = nullptr) const
    {
        eqns.clear();
        if (report)
            report->assign(topo.numConnections(), PdropReport{});
        // The well efficiency factor (WEFAC) scales the flow terms -- perforation inflow and
        // inter-segment / surface flux -- in both the well equation and the reservoir coupling,
        // but NOT the wellbore accumulation term (matching StandardWell and MultisegmentWell).
        assembleAccumulation(topo, pv, fp, dt, old_surface_volumes, regularization_factor, eqns);
        assemblePerforations(topo, pv, fp, perfs, controls.producer, allow_crossflow, gravity,
                             efficiency, eqns);
        for (int c = 0; c < topo.numConnections(); ++c) {
            assembleConnectionFlux(topo, pv, fp, c, efficiency, eqns);
            if (c == topo.surfaceConnection()) {
                // The control row lives on the surface connection. It can be assembled here
                // with the simple Controls struct, or (assemble_control=false) supplied by the
                // caller — GraphMultisegmentWell builds it natively to cover BHP/rate/THP/group.
                if (assemble_control)
                    assembleControl(topo, pv, fp, controls, c, eqns);
            } else {
                const ConnDevice dev = (c < static_cast<int>(devices.size()))
                    ? devices[c] : ConnDevice{};
                PdropReport* rep = report ? &(*report)[c] : nullptr;
                assemblePressureDrop(topo, pv, fp, gravity, c, pdrop, dev, eqns, rep);
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
                                Scalar efficiency, Eqns& eqns) const
    {
        const auto& conn = topo.connection(c);
        const Scalar Qval = pv.connRateValue(c);
        const int u = topo.upwindSegment(c, Qval);   // upwind segment (a real segment)
        const Eval Q = pv.connRate(c);
        for (int comp = 0; comp < NP; ++comp) {
            const Eval F = pv.volumeFractionScaled(u, comp, Role::Self);
            const Eval q = Q * F * efficiency;   // down->up component rate (>0 for production)
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

    //! r_c = p_down - p_up - dp_hydro - dp_friction/device - dp_acceleration
    void assemblePressureDrop(const Topo& topo, const PV& pv, const FP& fp, Scalar gravity,
                              int c, const PdropOptions& pdrop, const ConnDevice& device,
                              Eqns& eqns, PdropReport* report = nullptr) const
    {
        const auto& conn = topo.connection(c);
        const int down = conn.down;
        const int up = conn.up;

        // A shut flow-control device forces zero flux through the connection
        // (mirrors MultisegmentWell's trivial zero-rate equation for a SHUT valve).
        if (device.kind != ConnDevice::Kind::None && device.shut) {
            scatterConn(eqns, c, /*self*/ down, /*other*/ up, pv.connRate(c));
            return;
        }

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
        const Eval hydro_drop = rho_hydro * gravity * conn.depth_diff;
        r -= hydro_drop;
        if (report)
            report->hydrostatic = hydro_drop.value();

        const bool is_valve = (device.kind == ConnDevice::Kind::Valve);
        const bool is_sicd = (device.kind == ConnDevice::Kind::SpiralICD);
        const bool is_aicd = (device.kind == ConnDevice::Kind::AutoICD);
        if (pdrop.friction || pdrop.acceleration || is_valve || is_sicd || is_aicd) {
            // mass rate = Q * sum_comp F(upwind) * surf_dens
            const Eval Q = pv.connRate(c);
            Eval mass_rate(Scalar{0});
            for (int comp = 0; comp < NP; ++comp) {
                const Eval F = pv.volumeFractionScaled(u, comp, uprole);
                mass_rate += Q * F * fp.surfaceDensity(comp);
            }
            const Eval rho_up = fp.density(u, uprole);
            const Eval mu_up = fp.viscosity(u, uprole);
            const Scalar signQ = (Qval >= Scalar{0}) ? Scalar{1} : Scalar{-1};

            if (is_valve) {
                // Device drop replaces pipe friction: valve-pipe friction + constriction
                // (always applied, independent of the deck friction flag), as in
                // MultisegmentWellSegments::pressureDropValve.
                const Eval fric = graphwellhelpers::frictionPressureLoss(
                    device.length, device.diameter, device.area, device.roughness,
                    rho_up, mass_rate, mu_up);
                const Eval con = graphwellhelpers::valveContrictionPressureLoss(
                    mass_rate, rho_up, device.area_con, device.cv);
                const Eval drop = signQ * (fric + con);
                r -= drop;
                if (report) report->friction = drop.value();
            } else if (is_sicd) {
                // Spiral-ICD drop (MultisegmentWellSegments::pressureDropSpiralICD):
                // strength * (rho/rho_cal)^0.75 * (mu_mix/mu_cal)^0.25 * (scale*q_res)^2,
                // using per-phase reservoir volume fractions and an emulsion viscosity.
                const int wc = waterComp(), oc = oilComp(), gc = gasComp();
                auto frac = [&](int comp){ return comp >= 0 ? fp.phaseVolumeFraction(u, comp, uprole) : Eval{Scalar{0}}; };
                auto visc = [&](int comp){ return comp >= 0 ? fp.phaseViscosity(u, comp, uprole) : Eval{Scalar{0}}; };
                const Eval water_f = frac(wc), water_v = visc(wc);
                const Eval oil_f = frac(oc), oil_v = visc(oc);
                const Eval gas_f = frac(gc), gas_v = visc(gc);

                const Eval liquid_fraction = water_f + oil_f;
                const Eval liquid_visc_frac = (liquid_fraction < Scalar{1.e-30})
                    ? (oil_f * oil_v + water_f * water_v)
                    : liquid_fraction * graphwellhelpers::emulsionViscosity(
                          water_f, water_v, oil_f, oil_v,
                          device.width_transition, device.critical_value, device.max_visco_ratio);
                const Eval mixture_viscosity = liquid_visc_frac + gas_f * gas_v;

                const Eval reservoir_rate_icd = (mass_rate / rho_up) * device.scaling_factor;
                const Eval temp1 = (rho_up > Scalar{0})
                    ? pow(rho_up / device.dens_cali, Scalar{0.75}) : Eval{Scalar{0}};
                const Eval temp2 = (mixture_viscosity > Scalar{0})
                    ? pow(mixture_viscosity / device.visc_cali, Scalar{0.25}) : Eval{Scalar{0}};
                const Eval drop = signQ * temp1 * temp2 * device.strength
                     * reservoir_rate_icd * reservoir_rate_icd;
                r -= drop;
                if (report) report->friction = drop.value();
            } else if (is_aicd) {
                // Autonomous-ICD drop (MultisegmentWellSegments::pressureDropAutoICD):
                // per-phase exponent-weighted mixture viscosity/density, then
                // strength * (rho/rho_cal)^x * (visc_cal/visc)^y * rho * |q_icd|^z * unit^(2-z).
                const int wc = waterComp(), oc = oilComp(), gc = gasComp();
                auto frac = [&](int c){ return c >= 0 ? fp.phaseVolumeFraction(u, c, uprole) : Eval{Scalar{0}}; };
                auto visc = [&](int c){ return c >= 0 ? fp.phaseViscosity(u, c, uprole) : Eval{Scalar{0}}; };
                auto dens = [&](int c){ return c >= 0 ? fp.phaseDensity(u, c, uprole) : Eval{Scalar{0}}; };
                auto safe_pow = [](const Eval& a, Scalar b){ return (a > Scalar{0}) ? pow(a, b) : Eval{Scalar{0}}; };

                const Eval mixture_viscosity =
                      safe_pow(frac(wc), device.wat_visc_exp) * visc(wc)
                    + safe_pow(frac(oc), device.oil_visc_exp) * visc(oc)
                    + safe_pow(frac(gc), device.gas_visc_exp) * visc(gc);
                const Eval mixture_density =
                      safe_pow(frac(wc), device.wat_dens_exp) * dens(wc)
                    + safe_pow(frac(oc), device.oil_dens_exp) * dens(oc)
                    + safe_pow(frac(gc), device.gas_dens_exp) * dens(gc);

                const Eval volume_rate_icd = mass_rate * device.scaling_factor / mixture_density;
                const Eval abs_vri = abs(volume_rate_icd);
                const Eval dp = safe_pow(mixture_density / device.dens_cali, device.dens_exp)
                    * safe_pow(device.visc_cali / mixture_viscosity, device.visc_exp)
                    * mixture_density * device.strength
                    * safe_pow(abs_vri, device.flow_rate_exp)
                    * std::pow(device.unit_volume_rate, Scalar{2} - device.flow_rate_exp);
                const Eval drop = signQ * dp;
                r -= drop;
                if (report) report->friction = drop.value();
            } else if (pdrop.friction) {
                // friction: magnitude * sign(Q)
                const Eval fric = graphwellhelpers::frictionPressureLoss(
                    conn.length, conn.diameter, conn.area, conn.roughness, rho_up, mass_rate, mu_up);
                const Eval drop = signQ * fric;
                r -= drop;
                if (report) report->friction = drop.value();
            }
            if (pdrop.acceleration) {
                // acceleration: velocity-head difference across the connection (local form)
                const Eval rho_a = fp.density(down, Role::Self);
                const Eval rho_b = fp.density(up, Role::Other);
                const Eval accel = graphwellhelpers::velocityHead(conn.area, mass_rate, rho_b)
                                 - graphwellhelpers::velocityHead(conn.area, mass_rate, rho_a);
                r -= accel;
                if (report) report->acceleration = accel.value();
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
                              Scalar efficiency, Eqns& eqns) const
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
                                perfs[perf].perf_depth_diff, gravity,
                                producer, allow_cf, cq_s);
                for (int comp = 0; comp < NP; ++comp)
                    cq_s[comp] *= efficiency;   // WEFAC scales the reservoir<->well flow
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
