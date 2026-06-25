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

#ifndef OPM_GRAPHWELL_FLUID_PROPERTIES_HEADER_INCLUDED
#define OPM_GRAPHWELL_FLUID_PROPERTIES_HEADER_INCLUDED

#include <flowexperimental/graphwell/GraphWellPrimaryVariables.hpp>

#include <opm/material/densead/Math.hpp>

#include <algorithm>
#include <array>
#include <vector>

namespace Opm {

/// \brief Per-segment secondary (fluid) properties for the GraphWell model.
///
/// This is a minimal port of MultisegmentWellSegments::calculatePhaseProperties /
/// getSurfaceVolume onto the GraphWell derivative layout. All quantities are stored
/// in the "Self" derivative role of the segment they belong to; use \c shiftToOther
/// to relocate a value to the "Other" block when assembling a face term whose
/// upwind/density segment is the connection's \c up node.
///
/// Black-oil only; emulsion/polymer/solvent terms are dropped (HFA scope).
template<class FluidSystem, int NumPhases, int NumResEq>
class GraphWellFluidProperties
{
public:
    using PV = GraphWellPrimaryVariables<FluidSystem, NumPhases, NumResEq>;
    using Scalar = typename PV::Scalar;
    using Eval = typename PV::Eval;
    using Role = typename PV::Role;

    explicit GraphWellFluidProperties(int pvt_region) : pvt_region_(pvt_region)
    {
        surf_dens_.assign(NumPhases, Scalar{0});
        for (unsigned p = 0; p < FluidSystem::numPhases; ++p) {
            if (!FluidSystem::phaseIsActive(p))
                continue;
            const int comp = FluidSystem::canonicalToActiveCompIdx(
                FluidSystem::solventComponentIndex(p));
            surf_dens_[comp] = FluidSystem::referenceDensity(p, pvt_region_);
        }
    }

    Scalar surfaceDensity(int comp) const { return surf_dens_[comp]; }

    //! Recompute the secondary quantities for every segment from the primary state.
    void update(const PV& pv, Scalar temperature)
    {
        const int nseg = pv.numSegments();
        density_.assign(nseg, Eval{Scalar{0}});
        viscosity_.assign(nseg, Eval{Scalar{0}});
        vol_ratio_.assign(nseg, Eval{Scalar{0}});
        phase_volfrac_.assign(nseg, {});
        phase_visc_.assign(nseg, {});
        phase_dens_.assign(nseg, {});
        for (int s = 0; s < nseg; ++s)
            computeSegment(pv, s, temperature);
    }

    const Eval& mixtureDensity(int s) const { return density_[s]; }
    const Eval& mixtureViscosity(int s) const { return viscosity_[s]; }
    const Eval& volRatio(int s) const { return vol_ratio_[s]; }

    //! Move the derivatives of a Self-role value into the Other-role block.
    static Eval shiftToOther(const Eval& in)
    {
        Eval out(in.value());
        for (int k = 0; k < NumPhases; ++k)
            out.setDerivative(PV::OTH + k, in.derivative(PV::RES + k));
        return out;
    }

    //! Mixture density of segment \c s, in the requested role.
    Eval density(int s, Role role) const
    { return role == Role::Self ? density_[s] : shiftToOther(density_[s]); }

    //! Mixture viscosity of segment \c s, in the requested role.
    Eval viscosity(int s, Role role) const
    { return role == Role::Self ? viscosity_[s] : shiftToOther(viscosity_[s]); }

    //! Reservoir volume fraction of component \c comp (= mix/b/vol_ratio), in the
    //! requested role. Mirrors MultisegmentWellSegments::phase_fractions_.
    Eval phaseVolumeFraction(int s, int comp, Role role) const
    { return role == Role::Self ? phase_volfrac_[s][comp] : shiftToOther(phase_volfrac_[s][comp]); }

    //! Phase viscosity of component \c comp, in the requested role.
    Eval phaseViscosity(int s, int comp, Role role) const
    { return role == Role::Self ? phase_visc_[s][comp] : shiftToOther(phase_visc_[s][comp]); }

    //! Phase mass density of component \c comp, in the requested role.
    Eval phaseDensity(int s, int comp, Role role) const
    { return role == Role::Self ? phase_dens_[s][comp] : shiftToOther(phase_dens_[s][comp]); }

private:
    void computeSegment(const PV& pv, int s, Scalar temperature)
    {
        const Eval T(temperature);
        const Eval salt(Scalar{0});
        const Eval p = pv.segPressure(s, Role::Self);

        const bool wat = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
        const bool oil = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
        const bool gas = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
        const int wc = wat ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx) : -1;
        const int oc = oil ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::oilCompIdx) : -1;
        const int gc = gas ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::gasCompIdx) : -1;

        std::vector<Eval> mix_s(NumPhases, Eval{Scalar{0}});
        for (int c = 0; c < NumPhases; ++c)
            mix_s[c] = pv.surfaceVolumeFraction(s, c, Role::Self);

        std::vector<Eval> b(NumPhases, Eval{Scalar{0}});
        std::vector<Eval> visc(NumPhases, Eval{Scalar{0}});
        std::vector<Eval> dens(NumPhases, Eval{Scalar{0}});

        if (wat) {
            const Eval rsw(Scalar{0});
            b[wc] = FluidSystem::waterPvt().inverseFormationVolumeFactor(pvt_region_, T, p, rsw, salt);
            visc[wc] = FluidSystem::waterPvt().viscosity(pvt_region_, T, p, rsw, salt);
            dens[wc] = b[wc] * surf_dens_[wc];
        }

        Eval rv(Scalar{0});
        if (gas) {
            const Eval rvw(Scalar{0});
            const bool oil_exist = oil && mix_s[oc] > Scalar{0};
            if (oil_exist) {
                const Eval rvmax = max(FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvt_region_, T, p), Scalar{0});
                if (mix_s[gc] > Scalar{0})
                    rv = std::clamp(mix_s[oc] / mix_s[gc], Eval{Scalar{0}}, rvmax);
                b[gc] = FluidSystem::gasPvt().inverseFormationVolumeFactor(pvt_region_, T, p, rv, rvw);
                visc[gc] = FluidSystem::gasPvt().viscosity(pvt_region_, T, p, rv, rvw);
                dens[gc] = b[gc] * surf_dens_[gc] + rv * b[gc] * surf_dens_[oc];
            } else {
                b[gc] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(pvt_region_, T, p);
                visc[gc] = FluidSystem::gasPvt().saturatedViscosity(pvt_region_, T, p);
                dens[gc] = b[gc] * surf_dens_[gc];
            }
        }

        Eval rs(Scalar{0});
        if (oil) {
            const bool gas_exist = gas && mix_s[gc] > Scalar{0};
            if (gas_exist) {
                const Eval rsmax = max(FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvt_region_, T, p), Scalar{0});
                if (mix_s[oc] > Scalar{0})
                    rs = std::clamp(mix_s[gc] / mix_s[oc], Eval{Scalar{0}}, rsmax);
                b[oc] = FluidSystem::oilPvt().inverseFormationVolumeFactor(pvt_region_, T, p, rs);
                visc[oc] = FluidSystem::oilPvt().viscosity(pvt_region_, T, p, rs);
                dens[oc] = b[oc] * surf_dens_[oc] + rs * b[oc] * surf_dens_[gc];
            } else {
                b[oc] = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(pvt_region_, T, p);
                visc[oc] = FluidSystem::oilPvt().saturatedViscosity(pvt_region_, T, p);
                dens[oc] = b[oc] * surf_dens_[oc];
            }
        }

        std::vector<Eval> mix = mix_s;
        if (oil && gas) {
            const Eval d = Eval{Scalar{1}} - rs * rv;
            if (d > Scalar{0}) {
                if (rs > Scalar{0})
                    mix[gc] = (mix_s[gc] - mix_s[oc] * rs) / d;
                if (rv > Scalar{0})
                    mix[oc] = (mix_s[oc] - mix_s[gc] * rv) / d;
            }
        }

        Eval vol_ratio(Scalar{0});
        for (int c = 0; c < NumPhases; ++c)
            vol_ratio += mix[c] / b[c];
        vol_ratio_[s] = vol_ratio;

        // mixture density: sum surf_dens * mix_s, divided by vol_ratio
        Eval rho(Scalar{0});
        for (int c = 0; c < NumPhases; ++c)
            rho += surf_dens_[c] * mix_s[c];
        density_[s] = rho / vol_ratio;

        // mixture viscosity: sum visc_c * fraction_c, fraction_c = mix_c / b_c / vol_ratio
        // also store the per-phase reservoir volume fraction and viscosity (for ICD devices)
        Eval mu(Scalar{0});
        for (int c = 0; c < NumPhases; ++c) {
            const Eval frac = mix[c] / b[c] / vol_ratio;
            phase_volfrac_[s][c] = frac;
            phase_visc_[s][c] = visc[c];
            phase_dens_[s][c] = dens[c];
            mu += visc[c] * frac;
        }
        viscosity_[s] = mu;
    }

    int pvt_region_;
    std::vector<Scalar> surf_dens_;
    std::vector<Eval> density_;
    std::vector<Eval> viscosity_;
    std::vector<Eval> vol_ratio_;
    std::vector<std::array<Eval, NumPhases>> phase_volfrac_;
    std::vector<std::array<Eval, NumPhases>> phase_visc_;
    std::vector<std::array<Eval, NumPhases>> phase_dens_;
};

} // namespace Opm

#endif // OPM_GRAPHWELL_FLUID_PROPERTIES_HEADER_INCLUDED
