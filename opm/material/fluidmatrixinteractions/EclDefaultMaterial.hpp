// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011-2012 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Opm::EclDefaultMaterial
 */
#ifndef OPM_ECL_DEFAULT_MATERIAL_HPP
#define OPM_ECL_DEFAULT_MATERIAL_HPP

#include "EclDefaultMaterialParams.hpp"

#include <opm/material/fluidstates/SaturationOnlyFluidState.hpp>
#include <opm/material/Valgrind.hpp>

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <algorithm>

namespace Opm {

/*!
 * \ingroup material
 *
 * \brief Implements the default three phase capillary pressure law
 *        used by the ECLipse simulator.
 *
 * This material law is valid for three fluid phases and only depends
 * on the saturations.
 *
 * The required two-phase relations are supplied by means of template
 * arguments and can be an arbitrary other material laws. (Provided
 * that these only depend on saturation.)
 */
template <class ScalarT,
          int wPhaseIdxV,
          int oPhaseIdxV,
          int gPhaseIdxV,
          class WaterOilMaterialLaw,
          class OilGasMaterialLaw,
          class ParamsT = EclDefaultMaterialParams<OilGasMaterialLaw::Params,
                                                   WaterOilMaterialLaw::Params> >
class EclDefaultMaterial
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;
    enum { numPhases = 3 };

    enum { wPhaseIdx = wPhaseIdxV };
    enum { nPhaseIdx = nPhaseIdxV };
    enum { gPhaseIdx = gPhaseIdxV };

    /*!
     * \brief Implements the default three phase capillary pressure law
     *        used by the ECLipse simulator.
     *
     * This material law is valid for three fluid phases and only
     * depends on the saturations.
     *
     * The required two-phase relations are supplied by means of template
     * arguments and can be an arbitrary other material laws.
     *
     * \param values Container for the return values
     * \param params Parameters
     * \param state The fluid state
     */
    template <class ContainerT, class FluidState>
    static void capillaryPressures(ContainerT &values,
                                   const Params &params,
                                   const FluidState &state)
    {
        values[gPhaseIdx] = pcgn(params, state);
        values[oPhaseIdx] = 0;
        values[wPhaseIdx] = pcnw(params, state);
    }

    /*!
     * \brief Capillary pressure between the gas and the non-wetting
     *        liquid (i.e., oil) phase.
     *
     * This is defined as
     * \f[
     * p_{c,gn} = p_g - p_n
     * \f]
     */
    template <class FluidState>
    static Scalar pcgn(const Params &params,
                       const FluidState &state)
    {
        typedef SaturationOnlyFluidState<Scalar, /*numPhases=*/2> TwoPhaseFluidState;

        TwoPhaseFluidState twoPhaseFs;

        // calculate the relative permeabilities of water phase.
        twoPhaseFs.setSaturation(OilGasMaterial::wPhaseIdx, 1 - fluidState.saturation(gPhaseIdx));
        twoPhaseFs.setSaturation(OilGasMaterial::nPhaseIdx, fluidState.saturation(gPhaseIdx));

        return OilGasMaterialLaw::pcnw(params.gasOilParams(), twoPhaseFs);
    }

    /*!
     * \brief Capillary pressure between the non-wetting liquid (i.e.,
     *        oil) and the wetting liquid (i.e., water) phase.
     *
     * This is defined as
     * \f[
     * p_{c,nw} = p_n - p_w
     * \f]
     */
    template <class FluidState>
    static Scalar pcnw(const Params &params,
                       const FluidState &state)
    {
        typedef SaturationOnlyFluidState<Scalar, /*numPhases=*/2> TwoPhaseFluidState;

        TwoPhaseFluidState twoPhaseFs;

        // calculate the relative permeabilities of water phase.
        twoPhaseFs.setSaturation(WaterOilMaterial::wPhaseIdx, fluidState.saturation(wPhaseIdx));
        twoPhaseFs.setSaturation(WaterOilMaterial::nPhaseIdx, 1 - fluidState.saturation(wPhaseIdx));

        return WaterOilMaterialLaw::pcnw(params.gasOilParams(), twoPhaseFs);
    }

    /*!
     * \brief The inverse of the capillary pressure
     */
    template <class ContainerT, class FluidState>
    static void saturations(ContainerT &values,
                            const Params &params,
                            const FluidState &state)
    {
        OPM_THROW(std::runtime_error, "Not implemented: Stone1Material::saturations()");
    }

    /*!
     * \brief The relative permeability of all phases.
     *
     * The relative permeability of the water phase it uses the same
     * value as the relative permeability for water in the water-oil
     * law with \f$S_o = 1 - S_w\f$. The gas relative permebility is
     * taken from the gas-oil material law, but with \f$S_o = 1 -
     * S_g\f$.  The relative permeability of the oil phase is
     * calculated using the relative permeabilities of the oil phase
     * in the two two-phase systems.
     *
     * A more detailed description can be found in the "Three phase
     * oil relative permeability models" section of the ECLipse
     * technical description.
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params,
                                       const FluidState &fluidState)
    {
        values[wPhaseIdx] = krw(params, fluidState);
        values[nPhaseIdx] = krn(params, fluidState);
        values[gPhaseIdx] = krg(params, fluidState);
    }

    /*!
     * \brief The relative permeability of the gas phase.
     */
    template <class FluidState>
    static Scalar krg(const Params &params,
                      const FluidState &fluidState)
    {
        typedef SaturationOnlyFluidState<Scalar, /*numPhases=*/2> TwoPhaseFluidState;

        TwoPhaseFluidState twoPhaseFs;

        // calculate the relative permeabilities of water phase.
        twoPhaseFs.setSaturation(OilGasMaterial::wPhaseIdx, fluidState.saturation(wPhaseIdx));
        twoPhaseFs.setSaturation(OilGasMaterial::nPhaseIdx, 1 - fluidState.saturation(wPhaseIdx));

        return OilGasMaterial::krw(params.gasOilParams(), twoPhaseFs);
    }

    /*!
     * \brief The relative permeability of the wetting phase.
     */
    template <class FluidState>
    static Scalar krw(const Params &params,
                      const FluidState &fluidState)
    {
        typedef SaturationOnlyFluidState<Scalar, /*numPhases=*/2> TwoPhaseFluidState;

        TwoPhaseFluidState twoPhaseFs;

        // first, calculate the relative permeabilities of gas phase.
        twoPhaseFs.setSaturation(WaterOilMaterial::wPhaseIdx, 1 - fluidState.saturation(gPhaseIdx));
        twoPhaseFs.setSaturation(WaterOilMaterial::nPhaseIdx, fluidState.saturation(gPhaseIdx));

        return WaterOilMaterial::krn(params.gasOilParams(), twoPhaseFs);
    }

    /*!
     * \brief The relative permeability of the non-wetting (i.e., oil) phase.
     */
    template <class FluidState>
    static Scalar krn(const Params &params,
                      const FluidState &fluidState)
    {
        typedef SaturationOnlyFluidState<Scalar, /*numPhases=*/2> TwoPhaseFluidState;

        TwoPhaseFluidState twoPhaseFs;

        Scalar Sw = fluidState.saturation(wPhaseIdx);
        Scalar So = fluidState.saturation(oPhaseIdx);
        Scalar Sg = fluidState.saturation(gPhaseIdx);

        // connate water. According to the Eclipse TD, this is
        // probably only relevant if hysteresis is enabled...
        Scalar Swco = 0; // todo!

        // calculate the relative permeabilities of water phase.
        twoPhaseFs.setSaturation(OilGasMaterial::wPhaseIdx, So  - Swco);
        twoPhaseFs.setSaturation(OilGasMaterial::nPhaseIdx, 1 - (So - Swco) );

        Scalar krog = OilGasMaterial::krw(params.oilGasParams(), twoPhaseFs);

        // calculate the relative permeabilities of water phase.
        twoPhaseFs.setSaturation(WaterOilMaterial::wPhaseIdx, 1 - So);
        twoPhaseFs.setSaturation(WaterOilMaterial::nPhaseIdx, So);

        Scalar krow = WaterOilMaterial::krn(params.oilGasParams(), twoPhaseFs);

        if (Sg + Sw - Swco < 1e-30)
            return 0; // avoid division by zero
        else
            return (So * krog + (Sw - Swco)*krow) / (Sg + Sw - Swco);
    }
};
} // namespace Opm

#endif
