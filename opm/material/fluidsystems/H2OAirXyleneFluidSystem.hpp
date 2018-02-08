// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::FluidSystems::H2OAirXylene
 */
#ifndef OPM_H2O_AIR_XYLENE_FLUID_SYSTEM_HPP
#define OPM_H2O_AIR_XYLENE_FLUID_SYSTEM_HPP

#include <opm/material/IdealGas.hpp>
#include <opm/material/components/Air.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/components/Xylene.hpp>
#include <opm/material/components/TabulatedComponent.hpp>

#include <opm/material/binarycoefficients/H2O_Air.hpp>
#include <opm/material/binarycoefficients/H2O_Xylene.hpp>
#include <opm/material/binarycoefficients/Air_Xylene.hpp>

#include "BaseFluidSystem.hpp"
#include "NullParameterCache.hpp"

namespace Opm {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief A fluid system with water, gas and NAPL as phases and
 *        water, air and NAPL (contaminant) as components.
 */
template <class Scalar>
class H2OAirXylene
    : public BaseFluidSystem<Scalar, H2OAirXylene<Scalar> >
{
    typedef H2OAirXylene<Scalar> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

public:
    template <class Evaluation>
    struct ParameterCache : public Opm::NullParameterCache<Evaluation>
    {};

    //! The type of the water component
    typedef Opm::H2O<Scalar> H2O;
    //! The type of the xylene/napl component
    typedef Opm::Xylene<Scalar> NAPL;
    //! The type of the air component
    typedef Opm::Air<Scalar> Air;

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 3;
    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 3;

    //! The index of the water phase
    static const int waterPhaseIdx = 0;
    //! The index of the NAPL phase
    static const int naplPhaseIdx = 1;
    //! The index of the gas phase
    static const int gasPhaseIdx = 2;

    //! The index of the water component
    static const int H2OIdx = 0;
    //! The index of the NAPL component
    static const int NAPLIdx = 1;
    //! The index of the air pseudo-component
    static const int airIdx = 2;

    //! \copydoc BaseFluidSystem::init
    static void init()
    { }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gasPhaseIdx;
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(unsigned phaseIdx)
    { return phaseIdx == gasPhaseIdx && H2O::gasIsIdeal() && Air::gasIsIdeal() && NAPL::gasIsIdeal(); }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);

        // we assume Henry's and Rault's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(unsigned phaseIdx)
    {
        return
            (phaseIdx == gasPhaseIdx)
            // gases are always compressible
            ? true
            : (phaseIdx == waterPhaseIdx)
            // the water component decides for the water phase...
            ? H2O::liquidIsCompressible()
            // the NAPL component decides for the napl phase...
            : NAPL::liquidIsCompressible();
    }

    //! \copydoc BaseFluidSystem::phaseName
    static const char* phaseName(unsigned phaseIdx)
    {
        switch (phaseIdx) {
        case waterPhaseIdx: return "water";
        case naplPhaseIdx: return "napl";
        case gasPhaseIdx: return "gas";
        };
        throw std::logic_error("Invalid phase index "+std::to_string(phaseIdx));
    }

    //! \copydoc BaseFluidSystem::componentName
    static const char* componentName(unsigned compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::name();
        case airIdx: return Air::name();
        case NAPLIdx: return NAPL::name();
        };
        throw std::logic_error("Invalid component index "+std::to_string(compIdx));
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(unsigned compIdx)
    {
        return
            (compIdx == H2OIdx)
            // gases are always compressible
            ? H2O::molarMass()
            : (compIdx == airIdx)
            // the water component decides for the water comp...
            ? Air::molarMass()
            // the NAPL component decides for the napl comp...
            : (compIdx == NAPLIdx)
            ? NAPL::molarMass()
            : 1e30;
    }

    //! \copydoc BaseFluidSystem::density
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& /*paramCache*/,
                           unsigned phaseIdx)
    {
        if (phaseIdx == waterPhaseIdx) {
            const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
            const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

            // See: Ochs 2008
            // \todo: proper citation
            const LhsEval& rholH2O = H2O::liquidDensity(T, p);
            const LhsEval& clH2O = rholH2O/H2O::molarMass();

            const auto& xwH2O = Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, H2OIdx));
            const auto& xwAir = Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, airIdx));
            const auto& xwNapl = Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, NAPLIdx));
            // this assumes each dissolved molecule displaces exactly one
            // water molecule in the liquid
            return clH2O*(H2O::molarMass()*xwH2O + Air::molarMass()*xwAir + NAPL::molarMass()*xwNapl);
        }
        else if (phaseIdx == naplPhaseIdx) {
            // assume pure NAPL for the NAPL phase
            const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
            return NAPL::liquidDensity(T, LhsEval(1e30));
        }

        assert (phaseIdx == gasPhaseIdx);
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        const LhsEval& pH2O = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, H2OIdx))*p;
        const LhsEval& pAir = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, airIdx))*p;
        const LhsEval& pNAPL = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, NAPLIdx))*p;
        return
            H2O::gasDensity(T, pH2O) +
            Air::gasDensity(T, pAir) +
            NAPL::gasDensity(T, pNAPL);
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval viscosity(const FluidState& fluidState,
                             const ParameterCache<ParamCacheEval>& /*paramCache*/,
                             unsigned phaseIdx)
    {
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if (phaseIdx == waterPhaseIdx) {
            // assume pure water viscosity
            return H2O::liquidViscosity(T, p);
        }
        else if (phaseIdx == naplPhaseIdx) {
            // assume pure NAPL viscosity
            return NAPL::liquidViscosity(T, p);
        }

        assert (phaseIdx == gasPhaseIdx);

        /* Wilke method. See:
         *
         * See: R. Reid, et al.: The Properties of Gases and Liquids,
         * 4th edition, McGraw-Hill, 1987, 407-410
         * 5th edition, McGraw-Hill, 20001, p. 9.21/22
         *
         * in this case, we use a simplified version in order to avoid
         * computationally costly evaluation of sqrt and pow functions and
         * divisions
         * -- compare e.g. with Promo Class p. 32/33
         */
        const LhsEval mu[numComponents] = {
            H2O::gasViscosity(T, H2O::vaporPressure(T)),
            Air::simpleGasViscosity(T, p),
            NAPL::gasViscosity(T, NAPL::vaporPressure(T))
        };
        // molar masses
        const Scalar M[numComponents] = {
            H2O::molarMass(),
            Air::molarMass(),
            NAPL::molarMass()
        };

        const auto& xgAir = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, airIdx));
        const auto& xgH2O = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, H2OIdx));
        const auto& xgNapl = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, NAPLIdx));

        const LhsEval& xgAW = xgAir + xgH2O;
        const LhsEval& muAW = (mu[airIdx]*xgAir + mu[H2OIdx]*xgH2O)/ xgAW;

        const LhsEval& MAW = (xgAir*Air::molarMass() + xgH2O*H2O::molarMass())/xgAW;

        Scalar phiCAW = 0.3; // simplification for this particular system
        /* actually like this
         * Scalar phiCAW = std::pow(1.+std::sqrt(mu[NAPLIdx]/muAW)*std::pow(MAW/M[NAPLIdx],0.25),2)
         *                 / std::sqrt(8.*(1.+M[NAPLIdx]/MAW));
         */
        const LhsEval& phiAWC = phiCAW * muAW*M[NAPLIdx]/(mu[NAPLIdx]*MAW);

        return (xgAW*muAW)/(xgAW+xgNapl*phiAWC) + (xgNapl*mu[NAPLIdx])/(xgNapl + xgAW*phiCAW);
    }

    //! \copydoc BaseFluidSystem::diffusionCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval diffusionCoefficient(const FluidState& fluidState,
                                        const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                        unsigned phaseIdx,
                                        unsigned compIdx)
    {
        if (phaseIdx==gasPhaseIdx) {
            const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
            const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

            const LhsEval& diffAC = Opm::BinaryCoeff::Air_Xylene::gasDiffCoeff(T, p);
            const LhsEval& diffWC = Opm::BinaryCoeff::H2O_Xylene::gasDiffCoeff(T, p);
            const LhsEval& diffAW = Opm::BinaryCoeff::H2O_Air::gasDiffCoeff(T, p);

            const LhsEval& xga = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, airIdx));
            const LhsEval& xgw = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, H2OIdx));
            const LhsEval& xgc = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, NAPLIdx));

            if (compIdx==NAPLIdx) return (1.- xgw)/(xga/diffAW + xgc/diffWC);
            else if (compIdx==H2OIdx) return (1.- xgc)/(xgw/diffWC + xga/diffAC);
            else if (compIdx==airIdx) throw std::logic_error("Diffusivity of air in the gas phase "
                                                             "is constraint by sum of diffusive fluxes = 0 !\n");
        } else if (phaseIdx==waterPhaseIdx){
            Scalar diffACl = 1.e-9; // BinaryCoeff::Air_Xylene::liquidDiffCoeff(temperature, pressure);
            Scalar diffWCl = 1.e-9; // BinaryCoeff::H2O_Xylene::liquidDiffCoeff(temperature, pressure);
            Scalar diffAWl = 1.e-9; // BinaryCoeff::H2O_Air::liquidDiffCoeff(temperature, pressure);

            const LhsEval& xwa = Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, airIdx));
            const LhsEval& xww = Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, H2OIdx));
            const LhsEval& xwc = Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, NAPLIdx));

            switch (compIdx) {
            case NAPLIdx:
                return (1.- xww)/(xwa/diffAWl + xwc/diffWCl);
            case airIdx:
                return (1.- xwc)/(xww/diffWCl + xwa/diffACl);
            case H2OIdx:
                throw std::logic_error("Diffusivity of water in the water phase "
                                       "is constraint by sum of diffusive fluxes = 0 !\n");
            };
        } else if (phaseIdx==naplPhaseIdx) {

            throw std::logic_error("Diffusion coefficients of "
                                   "substances in liquid phase are undefined!\n");
        }
        return 0;
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval fugacityCoefficient(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx,
                                       unsigned compIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        assert(0 <= compIdx && compIdx < numComponents);

        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if (phaseIdx == waterPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            else if (compIdx == airIdx)
                return Opm::BinaryCoeff::H2O_Air::henry(T)/p;
            else if (compIdx == NAPLIdx)
                return Opm::BinaryCoeff::H2O_Xylene::henry(T)/p;
        }

        // for the NAPL phase, we assume currently that nothing is
        // dissolved. this means that the affinity of the NAPL
        // component to the NAPL phase is much higher than for the
        // other components, i.e. the fugacity cofficient is much
        // smaller.
        if (phaseIdx == naplPhaseIdx) {
            const LhsEval& phiNapl = NAPL::vaporPressure(T)/p;
            if (compIdx == NAPLIdx)
                return phiNapl;
            else if (compIdx == airIdx)
                return 1e6*phiNapl;
            else if (compIdx == H2OIdx)
                return 1e6*phiNapl;
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        assert(phaseIdx == gasPhaseIdx);
        return 1.0;
    }

    //! \copydoc BaseFluidSystem::enthalpy
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval enthalpy(const FluidState& fluidState,
                            const ParameterCache<ParamCacheEval>& /*paramCache*/,
                            unsigned phaseIdx)
    {
        const auto& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const auto& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if (phaseIdx == waterPhaseIdx) {
            return H2O::liquidEnthalpy(T, p);
        }
        else if (phaseIdx == naplPhaseIdx) {
            return NAPL::liquidEnthalpy(T, p);
        }
        else if (phaseIdx == gasPhaseIdx) {  // gas phase enthalpy depends strongly on composition
            const LhsEval& hgc = NAPL::gasEnthalpy(T, p);
            const LhsEval& hgw = H2O::gasEnthalpy(T, p);
            const LhsEval& hga = Air::gasEnthalpy(T, p);

            LhsEval result = 0;
            result += hgw * Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, H2OIdx));
            result += hga * Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, airIdx));
            result += hgc * Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, NAPLIdx));

            return result;
        }
        throw std::logic_error("Invalid phase index "+std::to_string(phaseIdx));
    }

private:
    template <class LhsEval>
    static LhsEval waterPhaseDensity_(const LhsEval& T,
                                      const LhsEval& pw,
                                      const LhsEval& xww,
                                      const LhsEval& xwa,
                                      const LhsEval& xwc)
    {
        const LhsEval& rholH2O = H2O::liquidDensity(T, pw);
        const LhsEval& clH2O = rholH2O/H2O::molarMass();

        // this assumes each dissolved molecule displaces exactly one
        // water molecule in the liquid
        return clH2O*(xww*H2O::molarMass() + xwa*Air::molarMass() + xwc*NAPL::molarMass());
    }

    template <class LhsEval>
    static LhsEval gasPhaseDensity_(const LhsEval& T,
                                    const LhsEval& pg,
                                    const LhsEval& xgw,
                                    const LhsEval& xga,
                                    const LhsEval& xgc)
    { return H2O::gasDensity(T, pg*xgw) + Air::gasDensity(T, pg*xga) + NAPL::gasDensity(T, pg*xgc); }

    template <class LhsEval>
    static LhsEval NAPLPhaseDensity_(const LhsEval& T, const LhsEval& pn)
    { return NAPL::liquidDensity(T, pn); }
};
} // namespace FluidSystems
} // namespace Opm

#endif
