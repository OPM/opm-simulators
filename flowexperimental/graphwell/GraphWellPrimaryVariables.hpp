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

#ifndef OPM_GRAPHWELL_PRIMARY_VARIABLES_HEADER_INCLUDED
#define OPM_GRAPHWELL_PRIMARY_VARIABLES_HEADER_INCLUDED

#include <opm/material/densead/Evaluation.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

namespace Opm {

/// \brief Primary variables for the GraphWell reformulation.
///
/// Per segment there are \c NumPhases DOFs: a pressure plus (NumPhases-1) surface
/// volume fractions. Per connection there is one DOF: the total surface-volume flux.
///
/// All quantities are exposed as fixed-size AD \c Eval values whose derivative slots
/// are laid out as
///   [0, NumResEq)              reservoir cell variables (filled by perforation terms)
///   [RES,  RES+NumPhases)      "self" segment   (pressure, fractions)
///   [OTH,  OTH+NumPhases)      "other" segment  (pressure, fractions)
///   [QIDX]                     connection flux
/// so that any face term (which couples two segments and a flux) fits in one Eval,
/// eliminating the reverse-flow "extra derivatives" double-evaluation of the
/// production MultisegmentWell.
template<class FluidSystem, int NumPhases, int NumResEq>
class GraphWellPrimaryVariables
{
public:
    using Scalar = typename FluidSystem::Scalar;

    static constexpr int numPhases = NumPhases;
    static constexpr int numResEq = NumResEq;

    //! number of AD derivative slots
    static constexpr int numDeriv = NumResEq + 2 * NumPhases + 1;

    using Eval = DenseAd::Evaluation<Scalar, numDeriv>;

    //! derivative-slot bases
    static constexpr int RES = NumResEq;                 //!< self segment block base
    static constexpr int OTH = NumResEq + NumPhases;     //!< other segment block base
    static constexpr int QIDX = NumResEq + 2 * NumPhases; //!< connection flux slot

    //! segment-local DOF: index 0 is pressure, 1..NumPhases-1 are fractions
    static constexpr int SegPres = 0;

    //! Which segment the "self"/"other" derivative block refers to in a face term.
    enum class Role { Self, Other };

    GraphWellPrimaryVariables() = default;

    //! Set up the component<->fraction-DOF mapping from active phases.
    //! \c ref_comp is the active component index of the "remainder" phase
    //! (oil if present, else the first active phase); \c frac_comp[k] (k=1..NP-1)
    //! is the active component index whose fraction lives in segment DOF \c k.
    void setPhaseMap(int ref_comp, const std::array<int, NumPhases>& frac_comp)
    {
        ref_comp_ = ref_comp;
        frac_comp_ = frac_comp;
    }

    //! Per-component rate scaling factors (mirrors WellInterface::scalingFactor). The
    //! defaults are 1 so that the scaled and surface fractions coincide.
    void setScale(const std::array<Scalar, NumPhases>& scale) { scale_ = scale; }

    void resize(int nseg, int nconn)
    {
        seg_.assign(nseg, std::array<Scalar, NumPhases>{});
        conn_.assign(nconn, Scalar{0});
    }

    int numSegments() const { return static_cast<int>(seg_.size()); }
    int numConnections() const { return static_cast<int>(conn_.size()); }

    // -- raw value access -----------------------------------------------------
    Scalar& segValue(int s, int dof) { return seg_[s][dof]; }
    Scalar segValue(int s, int dof) const { return seg_[s][dof]; }
    Scalar& connValue(int c) { return conn_[c]; }
    Scalar connValue(int c) const { return conn_[c]; }

    Scalar segPressureValue(int s) const { return seg_[s][SegPres]; }
    Scalar connRateValue(int c) const { return conn_[c]; }

    // -- AD accessors ---------------------------------------------------------
    int base(Role role) const { return role == Role::Self ? RES : OTH; }

    Eval segPressure(int s, Role role) const
    {
        Eval p(seg_[s][SegPres]);
        p.setDerivative(base(role) + SegPres, Scalar{1});
        return p;
    }

    Eval connRate(int c) const
    {
        Eval q(conn_[c]);
        q.setDerivative(QIDX, Scalar{1});
        return q;
    }

    //! Raw volume fraction of active component \c comp (the stored primary variable;
    //! the remainder component is 1 - sum of the others). These sum to 1.
    Eval volumeFraction(int s, int comp, Role role) const
    {
        const int b = base(role);
        if (comp == ref_comp_) {
            Eval f(Scalar{1});
            Scalar sum = 0;
            for (int k = 1; k < NumPhases; ++k) {
                sum += seg_[s][k];
                f.setDerivative(b + k, Scalar{-1});
            }
            f.setValue(Scalar{1} - sum);
            // Match MSW: a reference fraction that turns slightly negative due to
            // round-off / Newton overshoot is clamped to zero in value (derivatives kept).
            if (f.value() < Scalar{0})
                f.setValue(Scalar{0});
            return f;
        }
        for (int k = 1; k < NumPhases; ++k) {
            if (frac_comp_[k] == comp) {
                Eval f(seg_[s][k]);
                f.setDerivative(b + k, Scalar{1});
                return f;
            }
        }
        return Eval{Scalar{0}};
    }

    //! Scaled volume fraction F_p / g_p (used for component rates, matches MSW).
    Eval volumeFractionScaled(int s, int comp, Role role) const
    {
        return volumeFraction(s, comp, role) / scale_[comp];
    }

    //! Surface volume fraction (sums to 1): scaled fraction normalised by the sum.
    Eval surfaceVolumeFraction(int s, int comp, Role role) const
    {
        Eval sum(Scalar{0});
        for (int c = 0; c < NumPhases; ++c)
            sum += volumeFractionScaled(s, c, role);
        return volumeFractionScaled(s, comp, role) / sum;
    }

    // -- Newton update --------------------------------------------------------
    //! Update the state from a multitype solution increment (segment block and
    //! connection block), applying relaxation and basic clamping.
    template<class SegVector, class ConnVector>
    void updateNewton(const SegVector& dseg,
                      const ConnVector& dconn,
                      Scalar relaxation,
                      Scalar max_pressure_change)
    {
        for (int s = 0; s < numSegments(); ++s) {
            // pressure (clamped change)
            {
                Scalar dp = relaxation * dseg[s][SegPres];
                dp = std::clamp(dp, -max_pressure_change, max_pressure_change);
                seg_[s][SegPres] -= dp;
                seg_[s][SegPres] = std::max(seg_[s][SegPres], pressure_floor_);
            }
            // fractions
            for (int k = 1; k < NumPhases; ++k)
                seg_[s][k] -= relaxation * dseg[s][k];
            processFractions(s);
        }
        for (int c = 0; c < numConnections(); ++c)
            conn_[c] -= relaxation * dconn[c][0];
    }

private:
    //! Keep fractions in a physically sensible range, using the exact clamp/renormalize
    //! sequence of MultisegmentWellPrimaryVariables::processFractions (water, then gas,
    //! then oil): if a phase fraction is negative, set it to zero and divide the other
    //! active phase fractions by (1 - negative_value) so the set still sums to one.
    void processFractions(int s)
    {
        // Reconstruct the full active-component fraction vector with the reference
        // component carrying the remainder.
        std::array<Scalar, NumPhases> f{};
        Scalar sum_nonref = 0;
        for (int k = 1; k < NumPhases; ++k) {
            f[frac_comp_[k]] = seg_[s][k];
            sum_nonref += seg_[s][k];
        }
        f[ref_comp_] = Scalar{1} - sum_nonref;

        const bool wat = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
        const bool oil = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
        const bool gas = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
        const int wc = wat ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::waterCompIdx) : -1;
        const int oc = oil ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::oilCompIdx) : -1;
        const int gc = gas ? FluidSystem::canonicalToActiveCompIdx(FluidSystem::gasCompIdx) : -1;

        if (wat && f[wc] < Scalar{0}) {
            if (gas) f[gc] /= (Scalar{1} - f[wc]);
            if (oil) f[oc] /= (Scalar{1} - f[wc]);
            f[wc] = Scalar{0};
        }
        if (gas && f[gc] < Scalar{0}) {
            if (wat) f[wc] /= (Scalar{1} - f[gc]);
            if (oil) f[oc] /= (Scalar{1} - f[gc]);
            f[gc] = Scalar{0};
        }
        if (oil && f[oc] < Scalar{0}) {
            if (wat) f[wc] /= (Scalar{1} - f[oc]);
            if (gas) f[gc] /= (Scalar{1} - f[oc]);
            f[oc] = Scalar{0};
        }

        for (int k = 1; k < NumPhases; ++k)
            seg_[s][k] = f[frac_comp_[k]];
    }

    std::vector<std::array<Scalar, NumPhases>> seg_;
    std::vector<Scalar> conn_;

    int ref_comp_{0};
    std::array<int, NumPhases> frac_comp_{};
    std::array<Scalar, NumPhases> scale_ = []{ std::array<Scalar, NumPhases> s{};
        s.fill(Scalar{1}); return s; }();

    static constexpr Scalar pressure_floor_ = Scalar{1e3}; // ~0.01 bar
};

} // namespace Opm

#endif // OPM_GRAPHWELL_PRIMARY_VARIABLES_HEADER_INCLUDED
