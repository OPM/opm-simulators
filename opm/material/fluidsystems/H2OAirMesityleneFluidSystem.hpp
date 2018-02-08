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
 * \copydoc Opm::FluidSystems::H2OAirMesitylene
 */
#ifndef OPM_H2O_AIR_MESITYLENE_FLUID_SYSTEM_HPP
#define OPM_H2O_AIR_MESITYLENE_FLUID_SYSTEM_HPP

#include "BaseFluidSystem.hpp"
#include "NullParameterCache.hpp"

#include <opm/material/IdealGas.hpp>
#include <opm/material/components/N2.hpp>
#include <opm/material/components/Air.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/components/SimpleH2O.hpp>
#include <opm/material/components/Mesitylene.hpp>
#include <opm/material/components/TabulatedComponent.hpp>
#include <opm/material/binarycoefficients/H2O_Air.hpp>
#include <opm/material/binarycoefficients/H2O_N2.hpp>
#include <opm/material/binarycoefficients/H2O_Mesitylene.hpp>
#include <opm/material/binarycoefficients/Air_Mesitylene.hpp>

#include <iostream>

namespace Opm {
namespace FluidSystems {

/*!
 * \ingroup Fluidsystems
 * \brief A fluid system with water, gas and NAPL as phases and
 *        water, air and mesitylene (DNAPL) as components.
 */
template <class Scalar>
class H2OAirMesitylene
    : public BaseFluidSystem<Scalar, H2OAirMesitylene<Scalar> >
{
    typedef H2OAirMesitylene<Scalar> ThisType;
    typedef BaseFluidSystem<Scalar, ThisType> Base;

    typedef Opm::H2O<Scalar> IapwsH2O;
    typedef Opm::TabulatedComponent<Scalar, IapwsH2O, /*alongVaporPressure=*/false> TabulatedH2O;

public:
    template <class Evaluation>
    struct ParameterCache : public Opm::NullParameterCache<Evaluation>
    {};

    //! The type of the mesithylene/napl component
    typedef Opm::Mesitylene<Scalar> NAPL;

    //! The type of the air component
    typedef Opm::Air<Scalar> Air;

    //! The type of the water component
    //typedef SimpleH2O H2O;
    typedef TabulatedH2O H2O;
    //typedef IapwsH2O H2O;

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
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/50,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/50);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water [K]
     * \param tempMax The maximum temperature used for tabulation of water [K]
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water [Pa]
     * \param pressMax The maximum pressure used for tabulation of water [Pa]
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        if (H2O::isTabulated) {
            TabulatedH2O::init(tempMin, tempMax, nTemp,
                               pressMin, pressMax, nPress);
        }
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gasPhaseIdx;
    }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(unsigned phaseIdx)
    { return phaseIdx == gasPhaseIdx && H2O::gasIsIdeal() && Air::gasIsIdeal() && NAPL::gasIsIdeal(); }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(unsigned phaseIdx)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // gases are always compressible
        return (phaseIdx == gasPhaseIdx)
            ? true
            : (phaseIdx == waterPhaseIdx)
            ? H2O::liquidIsCompressible()
            : NAPL::liquidIsCompressible();
    }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(unsigned /*phaseIdx*/)
    {
        //assert(0 <= phaseIdx && phaseIdx < numPhases);
        // we assume Henry's and Rault's laws for the water phase and
        // and no interaction between gas molecules of different
        // components, so all phases are ideal mixtures!
        return true;
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
            ? H2O::molarMass()
            : (compIdx == airIdx)
            ? Air::molarMass()
            : (compIdx == NAPLIdx)
            ? NAPL::molarMass()
            : -1e10;
        //throw std::logic_error("Invalid component index "+std::to_string(compIdx));
    }

    //! \copydoc BaseFluidSystem::density
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval density(const FluidState& fluidState,
                           const ParameterCache<ParamCacheEval>& /*paramCache*/,
                           unsigned phaseIdx)
    {
        const LhsEval& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));

        if (phaseIdx == waterPhaseIdx) {
            // See: Ochs 2008
            const LhsEval& p =
                H2O::liquidIsCompressible()
                ? Opm::decay<LhsEval>(fluidState.pressure(phaseIdx))
                : 1e30;

            const LhsEval& rholH2O = H2O::liquidDensity(T, p);
            const LhsEval& clH2O = rholH2O/H2O::molarMass();

            // this assumes each dissolved molecule displaces exactly one
            // water molecule in the liquid
            return
                clH2O*(H2O::molarMass()*Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, H2OIdx)) +
                       Air::molarMass()*Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, airIdx)) +
                       NAPL::molarMass()*Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, NAPLIdx)));
        }
        else if (phaseIdx == naplPhaseIdx) {
            // assume pure NAPL for the NAPL phase
            const LhsEval& p =
                NAPL::liquidIsCompressible()
                ? Opm::decay<LhsEval>(fluidState.pressure(phaseIdx))
                : 1e30;
            return NAPL::liquidDensity(T, p);
        }

        assert (phaseIdx == gasPhaseIdx);
        const LhsEval& pg = Opm::decay<LhsEval>(fluidState.pressure(gasPhaseIdx));
        const LhsEval& pH2O = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, H2OIdx))*pg;
        const LhsEval& pAir = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, airIdx))*pg;
        const LhsEval& pNAPL = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, NAPLIdx))*pg;
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
        const LhsEval& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if (phaseIdx == waterPhaseIdx) {
            // assume pure water viscosity

            return H2O::liquidViscosity(T,
                                        p);
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
            Air::gasViscosity(T, p),
            NAPL::gasViscosity(T, NAPL::vaporPressure(T))
        };
        // molar masses
        const Scalar M[numComponents] = {
            H2O::molarMass(),
            Air::molarMass(),
            NAPL::molarMass()
        };

        const LhsEval& xgAir = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, airIdx));
        const LhsEval& xgH2O = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, H2OIdx));
        const LhsEval& xgNapl = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, NAPLIdx));
        const LhsEval& xgAW = xgAir + xgH2O;
        const LhsEval& muAW = (mu[airIdx]*xgAir + mu[H2OIdx]*xgH2O)/xgAW;
        const LhsEval& MAW = (xgAir*Air::molarMass() + xgH2O*H2O::molarMass())/xgAW;

        Scalar phiCAW = 0.3; // simplification for this particular system
        /* actually like this
         * Scalar phiCAW = std::pow(1.+std::sqrt(mu[NAPLIdx]/muAW)*std::pow(MAW/M[NAPLIdx],0.25),2)
         *                 / std::sqrt(8.*(1.+M[NAPLIdx]/MAW));
         */
        const LhsEval& phiAWC = phiCAW * muAW*M[NAPLIdx]/(mu[NAPLIdx]*MAW);

        return (xgAW*muAW)/(xgAW + xgNapl*phiAWC) + (xgNapl*mu[NAPLIdx])/(xgNapl + xgAW*phiCAW);
    }

    //! \copydoc BaseFluidSystem::diffusionCoefficient
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval diffusionCoefficient(const FluidState& /*fluidState*/,
                                        const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                        unsigned /*phaseIdx*/,
                                        unsigned /*compIdx*/)
    {
        return 0;
#if 0
        typedef Opm::MathToolbox<typename FluidState::Scalar> FsToolbox;

        const LhsEval& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        LhsEval diffCont;

        if (phaseIdx==gasPhaseIdx) {
            const LhsEval& diffAC = Opm::BinaryCoeff::Air_Mesitylene::gasDiffCoeff(T, p);
            const LhsEval& diffWC = Opm::BinaryCoeff::H2O_Mesitylene::gasDiffCoeff(T, p);
            const LhsEval& diffAW = Opm::BinaryCoeff::H2O_Air::gasDiffCoeff(T, p);

            const LhsEval& xga = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, airIdx));
            const LhsEval& xgw = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, H2OIdx));
            const LhsEval& xgc = Opm::decay<LhsEval>(fluidState.moleFraction(gasPhaseIdx, NAPLIdx));

            if (compIdx==NAPLIdx) return (1 - xgw)/(xga/diffAW + xgc/diffWC);
            else if (compIdx==H2OIdx) return (1 - xgc)/(xgw/diffWC + xga/diffAC);
            else if (compIdx==airIdx) throw std::logic_error("Diffusivity of air in the gas phase "
                                                             "is constraint by sum of diffusive fluxes = 0 !\n");
        }
        else if (phaseIdx==waterPhaseIdx){
            const LhsEval& diffACl = 1.e-9; // BinaryCoeff::Air_Mesitylene::liquidDiffCoeff(temperature, pressure);
            const LhsEval& diffWCl = 1.e-9; // BinaryCoeff::H2O_Mesitylene::liquidDiffCoeff(temperature, pressure);
            const LhsEval& diffAWl = 1.e-9; // BinaryCoeff::H2O_Air::liquidDiffCoeff(temperature, pressure);

            const LhsEval& xwa = Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, airIdx));
            const LhsEval& xww = Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, H2OIdx));
            const LhsEval& xwc = Opm::decay<LhsEval>(fluidState.moleFraction(waterPhaseIdx, NAPLIdx));

            switch (compIdx) {
            case NAPLIdx:
                diffCont = (1.- xww)/(xwa/diffAWl + xwc/diffWCl);
                return diffCont;
            case airIdx:
                diffCont = (1.- xwc)/(xww/diffWCl + xwa/diffACl);
                return diffCont;
            case H2OIdx:
                throw std::logic_error("Diffusivity of water in the water phase "
                                       "is constraint by sum of diffusive fluxes = 0 !\n");
            };
        }
        else if (phaseIdx==naplPhaseIdx) {
            throw std::logic_error("Diffusion coefficients of "
                                   "substances in liquid phase are undefined!\n");
        }
        return 0;
#endif
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

        const LhsEval& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));
        Valgrind::CheckDefined(T);
        Valgrind::CheckDefined(p);

        if (phaseIdx == waterPhaseIdx) {
            if (compIdx == H2OIdx)
                return H2O::vaporPressure(T)/p;
            else if (compIdx == airIdx)
                return Opm::BinaryCoeff::H2O_N2::henry(T)/p;
            else if (compIdx == NAPLIdx)
                return Opm::BinaryCoeff::H2O_Mesitylene::henry(T)/p;
            assert(false);
        }
        // for the NAPL phase, we assume currently that nothing is
        // dissolved. this means that the affinity of the NAPL
        // component to the NAPL phase is much higher than for the
        // other components, i.e. the fugacity cofficient is much
        // smaller.
        else if (phaseIdx == naplPhaseIdx) {
            const LhsEval& phiNapl = NAPL::vaporPressure(T)/p;
            if (compIdx == NAPLIdx)
                return phiNapl;
            else if (compIdx == airIdx)
                return 1e6*phiNapl;
            else if (compIdx == H2OIdx)
                return 1e6*phiNapl;
            assert(false);
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
        const LhsEval& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
        const LhsEval& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

        if (phaseIdx == waterPhaseIdx) {
            return H2O::liquidEnthalpy(T, p);
        }
        else if (phaseIdx == naplPhaseIdx) {
            return NAPL::liquidEnthalpy(T, p);
        }
        else if (phaseIdx == gasPhaseIdx) {
            // gas phase enthalpy depends strongly on composition
            LhsEval result = 0;
            result += H2O::gasEnthalpy(T, p) * Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, H2OIdx));
            result += NAPL::gasEnthalpy(T, p) * Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, airIdx));
            result += Air::gasEnthalpy(T, p) * Opm::decay<LhsEval>(fluidState.massFraction(gasPhaseIdx, NAPLIdx));

            return result;
        }
        throw std::logic_error("Invalid phase index "+std::to_string(phaseIdx));
    }

    //! \copydoc BaseFluidSystem::thermalConductivity
    template <class FluidState, class LhsEval = typename FluidState::Scalar, class ParamCacheEval = LhsEval>
    static LhsEval thermalConductivity(const FluidState& fluidState,
                                       const ParameterCache<ParamCacheEval>& /*paramCache*/,
                                       unsigned phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        if (phaseIdx == waterPhaseIdx){ // water phase
            const LhsEval& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
            const LhsEval& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

            return H2O::liquidThermalConductivity(T, p);
        }
        else if (phaseIdx == gasPhaseIdx) { // gas phase
            const LhsEval& T = Opm::decay<LhsEval>(fluidState.temperature(phaseIdx));
            const LhsEval& p = Opm::decay<LhsEval>(fluidState.pressure(phaseIdx));

            return Air::gasThermalConductivity(T, p);
        }

        assert(phaseIdx == naplPhaseIdx);

        // Taken from:
        //
        // D. K. H. Briggs: "Thermal Conductivity of Liquids",
        // Ind. Eng. Chem., 1957, 49 (3), pp 418–421
        //
        // Convertion to SI units:
        // 344e-6 cal/(s cm K) = 0.0143964 J/(s m K)
        return 0.0143964;
    }
};
} // namespace FluidSystems
} // namespace Opm

#endif
