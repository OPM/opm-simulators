// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2008-2012 by Andreas Lauser
  Copyright (C) 2010 by Philipp Nuske
  Copyright (C) 2010 by Bernd Flemisch

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
*/
/*!
 * \file
 * \copydoc Opm::VanGenuchten
 */
#ifndef OPM_VAN_GENUCHTEN_HPP
#define OPM_VAN_GENUCHTEN_HPP

#include "VanGenuchtenParams.hpp"

#include <algorithm>
#include <cmath>
#include <cassert>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Implementation of the van Genuchten capillary pressure -
       saturation relation.
 *
 * This class only implements the "raw" van-Genuchten curves as static
 * members and doesn't concern itself converting absolute to effective
 * saturations and vice versa.
 *
 * The converion from and to effective saturations can be done using,
 * e.g. EffToAbsLaw.
 *
 * \see VanGenuchtenParams
 */
template <class TraitsT, class ParamsT = VanGenuchtenParams<TraitsT> >
class VanGenuchten : public TraitsT
{
public:
    //! The traits class for this material law
    typedef TraitsT Traits;

    //! The type of the parameter objects for this law
    typedef ParamsT Params;

    //! The type of the scalar values for this law
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases
    static const int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The van Genuchten capillary pressure law only "
                  "applies to the case of two fluid phases");

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = true;

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = true;

    //! Specify whether the quantities defined by this material law
    //! are saturation dependent
    static const bool isSaturationDependent = true;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the absolute pressure
    static const bool isPressureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are temperature dependent
    static const bool isTemperatureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the phase composition
    static const bool isCompositionDependent = false;

    /*!
     * \brief The capillary pressure-saturation curves according to van Genuchten.
     *
     * Van Genuchten's empirical capillary pressure <-> saturation
     * function is given by
     * \f[
     * p_{c,wn} = p_n - p_w = ({S_w}^{-1/m} - 1)^{1/n}/\alpha
     * \f]
     *
     * \param values A random access container which stores the
     *               relative pressure of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the capillary pressure
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container &values, const Params &params, const FluidState &fs)
    {
        values[Traits::wPhaseIdx] = 0.0; // reference phase
        values[Traits::nPhaseIdx] = pcnw(params, fs);
    }

    /*!
     * \brief Calculate the saturations of the phases starting from
     *        their pressure differences.
     */
    template <class Container, class FluidState>
    static void saturations(Container &values, const Params &params, const FluidState &fs)
    {
        values[Traits::wPhaseIdx] = Sw(params, fs);
        values[Traits::nPhaseIdx] = 1 - values[Traits::wPhaseIdx];
    }

    /*!
     * \brief The relative permeability-saturation curves according to van Genuchten.
     *
     * \param values A random access container which stores the
     *               relative permeability of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the relative permeabilities
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void relativePermeabilities(Container &values, const Params &params, const FluidState &fs)
    {
        values[Traits::wPhaseIdx] = krw(params, fs);
        values[Traits::nPhaseIdx] = krn(params, fs);
    }


    /*!
     * \brief The derivative of all capillary pressures in regard to
     *        a given phase saturation.
     */
    template <class ContainerT, class FluidState>
    static void dCapillaryPressures_dSaturation(ContainerT &values,
                                                const Params &params,
                                                const FluidState &state,
                                                int satPhaseIdx)
    {
        values[Traits::wPhaseIdx] = 0;
        values[Traits::nPhaseIdx] = 0;
        if (satPhaseIdx == Traits::wPhaseIdx)
            values[Traits::nPhaseIdx] = dPcnw_dSw(params, state);
    }

    /*!
     * \brief The derivative of all capillary pressures in regard to
     *        a given phase pressure.
     */
    template <class ContainerT, class FluidState>
    static void dCapillaryPressures_dPressure(ContainerT &values,
                                              const Params &params,
                                              const FluidState &state,
                                              int pPhaseIdx)
    {
        // -> not pressure dependent
        for (int pcPhaseIdx = 0; pcPhaseIdx < numPhases; ++pcPhaseIdx)
            values[pcPhaseIdx] = 0.0;
    }

    /*!
     * \brief The derivative of all capillary pressures in regard to
     *        temperature.
     */
    template <class ContainerT, class FluidState>
    static void dCapillaryPressures_dTemperature(ContainerT &values,
                                                 const Params &params,
                                                 const FluidState &state)
    {
        // -> not temperature dependent
        for (int pcPhaseIdx = 0; pcPhaseIdx < numPhases; ++pcPhaseIdx)
            values[pcPhaseIdx] = 0.0;
    }

    /*!
     * \brief The derivative of all capillary pressures in regard to
     *        a given mole fraction of a component in a phase.
     */
    template <class ContainerT, class FluidState>
    static void dCapillaryPressures_dMoleFraction(ContainerT &values,
                                                  const Params &params,
                                                  const FluidState &state,
                                                  int phaseIdx,
                                                  int compIdx)
    {
        // -> not composition dependent
        for (int pcPhaseIdx = 0; pcPhaseIdx < numPhases; ++pcPhaseIdx)
            values[pcPhaseIdx] = 0.0;
    }

    /*!
     * \brief The derivative of all relative permeabilities in regard to
     *        a given phase saturation.
     */
    template <class ContainerT, class FluidState>
    static void dRelativePermeabilities_dSaturation(ContainerT &values,
                                                    const Params &params,
                                                    const FluidState &state,
                                                    int satPhaseIdx)
    {
        if (satPhaseIdx == Traits::wPhaseIdx) {
            values[Traits::wPhaseIdx] = twoPhaseSatDKrw_dSw(params, state.saturation(Traits::wPhaseIdx));
            values[Traits::nPhaseIdx] = 0;
        }
        else {
            values[Traits::wPhaseIdx] = 0;
            values[Traits::nPhaseIdx] = - twoPhaseSatDKrn_dSw(params, 1 - state.saturation(Traits::nPhaseIdx));
        }
    }

    /*!
     * \brief The derivative of all relative permeabilities in regard to
     *        a given phase pressure.
     */
    template <class ContainerT, class FluidState>
    static void dRelativePermeabilities_dPressure(ContainerT &values,
                                                  const Params &params,
                                                  const FluidState &state,
                                                  int pPhaseIdx)
    {
        // -> not pressure dependent
        for (int krPhaseIdx = 0; krPhaseIdx < numPhases; ++krPhaseIdx)
            values[krPhaseIdx] = 0.0;
    }

    /*!
     * \brief The derivative of all relative permeabilities in regard to
     *        temperature.
     */
    template <class ContainerT, class FluidState>
    static void dRelativePermeabilities_dTemperature(ContainerT &values,
                                                     const Params &params,
                                                     const FluidState &state)
    {
        // -> not temperature dependent
        for (int krPhaseIdx = 0; krPhaseIdx < numPhases; ++krPhaseIdx)
            values[krPhaseIdx] = 0.0;
    }

    /*!
     * \brief The derivative of all relative permeabilities in regard to
     *        a given mole fraction of a component in a phase.
     */
    template <class ContainerT, class FluidState>
    static void dRelativePermeabilities_dMoleFraction(ContainerT &values,
                                                      const Params &params,
                                                      const FluidState &state,
                                                      int phaseIdx,
                                                      int compIdx)
    {
        // -> not composition dependent
        for (int krPhaseIdx = 0; krPhaseIdx < numPhases; ++krPhaseIdx)
            values[krPhaseIdx] = 0.0;
    }

    /*!
     * \brief The capillary pressure-saturation curve according to van Genuchten.
     *
     * Van Genuchten's empirical capillary pressure <-> saturation
     * function is given by
     * \f[
     * p_{c,wn} = p_n - p_w = ({S_w}^{-1/m} - 1)^{1/n}/\alpha
     * \f]
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the capillary pressure
     *           ought to be calculated
     */
    template <class FluidState>
    static Scalar pcnw(const Params &params, const FluidState &fs)
    {
        Scalar Sw = fs.saturation(Traits::wPhaseIdx);
        assert(0 <= Sw && Sw <= 1);

        return twoPhaseSatPcnw(params, Sw);
    }

    /*!
     * \brief The saturation-capillary pressure curve according to van
     *        Genuchten using a material law specific API.
     *
     * The advantage of this model is that it is simpler to use
     * because the baggage of the fluid state API does not need to be
     * carried along. The disavantage of this is, that it is very
     * specific to the van Genuchten law (i.e., depends only on the
     * wetting phase saturation, assumes two fluid phases, etc)
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param Sw The effective wetting phase saturation
     */
    static Scalar twoPhaseSatPcnw(const Params &params, Scalar Sw)
    { return std::pow(std::pow(Sw, -1.0/params.vgM()) - 1, 1.0/params.vgN())/params.vgAlpha(); }

    /*!
     * \brief The saturation-capillary pressure curve according to van Genuchten.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     * S_w = {p_C}^{-1} = ((\alpha p_C)^n + 1)^{-m}
     * \f]
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state containing valid phase pressures
     */
    template <class FluidState>
    static Scalar Sw(const Params &params, const FluidState &fs)
    {
        Scalar pC = fs.pressure(Traits::nPhaseIdx) - fs.pressure(Traits::wPhaseIdx);
        return twoPhaseSatSw(params, pC);
    }

    static Scalar twoPhaseSatSw(const Params &params, Scalar pC)
    {
        assert(pC >= 0);

        return std::pow(std::pow(params.vgAlpha()*pC, params.vgN()) + 1, -params.vgM());
    }

    /*!
     * \brief Calculate the non-wetting phase saturations depending on
     *        the phase pressures.
     */
    template <class FluidState>
    static Scalar Sn(const Params &params, const FluidState &fs)
    { return 1 - Sw(params, fs); }

    static Scalar twoPhaseSatSn(const Params &params, Scalar pC)
    { return 1 - twoPhaseSatSw(params, pC); }

    /*!
     * \brief The partial derivative of the capillary pressure with
     *        regard to the saturation according to van Genuchten.
     *
     * This is equivalent to
     * \f[
     * \frac{\partial p_C}{\partial S_w} =
     * -\frac{1}{\alpha n} ({S_w}^{-1/m} - 1)^{1/n - 1} {S_w}^{-1/m - 1}
     * \f]
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state containing valid saturations
     */
    template <class FluidState>
    static Scalar dPcnw_dSw(const Params &params, const FluidState &fs)
    { return twoPhaseSatDPcnw_dSw(params, fs.saturation(Traits::wPhaseIdx)); }

    static Scalar twoPhaseSatDPcnw_dSw(const Params &params, Scalar Sw)
    {
        assert(0 < Sw && Sw < 1);

        Scalar powSw = std::pow(Sw, -1/params.vgM());
        return
            - 1/(params.vgAlpha() * params.vgN() * Sw)
            * std::pow(powSw - 1, 1/params.vgN() - 1) * powSw;
    }

    /*!
     * \brief The relative permeability for the wetting phase of the
     *        medium according to van Genuchten's curve with Mualem
     *        parameterization.
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the relative permeability
     *           ought to be calculated
     */
    template <class FluidState>
    static Scalar krw(const Params &params, const FluidState &fs)
    { return twoPhaseSatKrw(params, fs.saturation(Traits::wPhaseIdx)); }

    static Scalar twoPhaseSatKrw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        Scalar r = 1. - std::pow(1 - std::pow(Sw, 1/params.vgM()), params.vgM());
        return std::sqrt(Sw)*r*r;
    }

    /*!
     * \brief The derivative of the relative permeability of the
     *        wetting phase in regard to the wetting saturation of the
     *        medium implied according to the van Genuchten curve with
     *        Mualem parameters.
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the derivative
     *           ought to be calculated
     */
    template <class FluidState>
    static Scalar dKrw_dSw(const Params &params, const FluidState &fs)
    { return twoPhaseSatDkrw_dSw(params, fs.saturation(Traits::wPhaseIdx)); }

    static Scalar twoPhaseSatDKrw_dSw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        const Scalar x = 1 - std::pow(Sw, 1.0/params.vgM());
        const Scalar xToM = std::pow(x, params.vgM());
        return (1 - xToM)/std::sqrt(Sw) * ( (1 - xToM)/2 + 2*xToM*(1-x)/x );
    }


    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium according to van Genuchten.
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the derivative
     *           ought to be calculated
     */
    template <class FluidState>
    static Scalar krn(const Params &params, const FluidState &fs)
    { return twoPhaseSatKrn(params, 1.0 - fs.saturation(Traits::nPhaseIdx)); }

    static Scalar twoPhaseSatKrn(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        return
            std::pow(1 - Sw, 1.0/3) *
            std::pow(1 - std::pow(Sw, 1/params.vgM()), 2*params.vgM());
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the van Genuchten
     *        parameterization.
     *
     * \param params The parameter object expressing the coefficients
     *               required by the van Genuchten law.
     * \param fs The fluid state for which the derivative
     *           ought to be calculated
     */
    template <class FluidState>
    static Scalar dKrn_dSw(const Params &params, const FluidState &fs)
    { return twoPhaseSatDkrn_dSw(params, fs.saturation(Traits::wPhaseIdx)); }

    static Scalar twoPhaseSatDKrn_dSw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        const Scalar x = std::pow(Sw, 1.0/params.vgM());
        return
            -std::pow(1 - x, 2*params.vgM())
            *std::pow(1 - Sw, -2/3)
            *(1.0/3 + 2*x/Sw);
    }
};
} // namespace Opm

#endif
