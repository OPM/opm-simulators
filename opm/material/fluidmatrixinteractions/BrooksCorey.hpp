// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2008-2013 by Andreas Lauser
  Copyright (C) 2010 by Philipp Nuske

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
 * \copydoc Opm::BrooksCorey
 */
#ifndef OPM_BROOKS_COREY_HPP
#define OPM_BROOKS_COREY_HPP

#include "BrooksCoreyParams.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Implementation of the Brooks-Corey capillary pressure <->
       saturation relation.
 *
 * This class provides the "raw" curves as static members and doesn't
 * concern itself converting absolute to effective saturations and
 * vice versa.
 *
 *\see BrooksCoreyParams
 */
template <class TraitsT, class ParamsT = BrooksCoreyParams<TraitsT> >
class BrooksCorey : public TraitsT
{
public:
    typedef TraitsT Traits;
    typedef ParamsT Params;
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases to which this material law applies.
    static const int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The Brooks-Corey capillary pressure law only applies "
                  "to the case of two fluid phases");

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

    static_assert(Traits::numPhases == 2,
                  "The number of fluid phases must be two if you want to use "
                  "this material law!");

    /*!
     * \brief The capillary pressure-saturation curves.
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
     * \brief The relative permeability-saturation curves.
     *
     * \param values A random access container which stores the
     *               relative permeability of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the material law.
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
     * \brief The capillary pressure-saturation curve according to
     *        Brooks and Corey.
     *
     * The empirical Brooks-Corey capillary pressure-saturation
     * function is defined as
     * \f[
     * p_C = p_e\overline{S}_w^{-1/\lambda}
     * \f]
     *
     * \param Sw Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    template <class FluidState>
    static Scalar pcnw(const Params &params, const FluidState &fs)
    {
        Scalar Sw = fs.saturation(Traits::wPhaseIdx);
        return twoPhaseSatPcnw(params, Sw);
    }

    static Scalar twoPhaseSatPcnw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        return params.entryPressure()*std::pow(Sw, -1.0/params.lambda());
    }

    /*!
     * \brief The saturation-capillary pressure curve according to
     *        Brooks & Corey.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     \overline{S}_w = (\frac{p_C}{p_e})^{-\lambda}
     \f]
     *
     * \param pc Capillary pressure \f$[Pa]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    template <class FluidState>
    static Scalar Sw(const Params &params, const FluidState &fs)
    {
        Scalar pc = fs.pressure(Traits::nPhaseIdx) - fs.pressure(Traits::wPhaseIdx);
        return twoPhaseSatSw(params, pc);
    }

    static Scalar twoPhaseSatSw(const Params &params, Scalar pc)
    {
        assert(pc > 0); // if we don't assume that, std::pow will screw up!

        return std::pow(pc/params.entryPressure(), -params.lambda());
    }

    /*!
     * \brief Calculate the non-wetting phase saturations depending on
     *        the phase pressures.
     */
    template <class FluidState>
    static Scalar Sn(const Params &params, const FluidState &fs)
    { return 1 - Sw(params, fs); }

    static Scalar twoPhaseSatSn(const Params &params, Scalar pc)
    { return 1 - twoPhaseSatSw(params, pc); }

    /*!
     * \brief The partial derivative of the capillary
     *        pressure w.r.t. the effective saturation according to Brooks & Corey.
     *
     * This is equivalent to
     * \f[
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{p_e}{\lambda} \overline{S}_w^{-1/\lambda - 1}
     \f]
     *
     * \param Sw Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
    */
    template <class FluidState>
    static Scalar dPcnw_dSw(const Params &params, const FluidState &fs)
    {
        Scalar Sw = fs.saturation(Traits::wPhaseIdx);
        return twoPhaseSatDPcnw_dSw(params, Sw);
    }

    static Scalar twoPhaseSatDPcnw_dSw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);
        return - params.entryPressure()/params.lambda() * std::pow(Sw, -1/params.lambda() - 1);
    }

    /*!
     * \brief The partial derivative of the effective saturation with
     *        regard to the capillary pressure according to Brooks and
     *        Corey.
     *
     * \param pcnw Capillary pressure \f$[Pa]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    template <class FluidState>
    static Scalar dSw_dpcnw(const Params &params, const FluidState &fs)
    {
        Scalar Sw = fs.saturation(Traits::wPhaseIdx);
        return twoPhaseSatDSw_dpcnw(params, Sw);
    }

    static Scalar twoPhaseSatDSw_dpcnw(const Params &params, Scalar Sw)
    {
        assert(pcnw > 0); // required for std::pow
        return -params.lambda()/params.entryPressure() * std::pow(pcnw/params.entryPressure(), - params.lambda() - 1);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Sw Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    template <class FluidState>
    static Scalar krw(const Params &params, const FluidState &fs)
    { return twoPhaseSatKrw(params, fs.saturation(Traits::wPhaseIdx)); }

    static Scalar twoPhaseSatKrw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        return std::pow(Sw, 2.0/params.lambda() + 3);
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        wetting phase with regard to the wetting saturation of the
     *        medium implied by the Brooks-Corey parameterization.
     *
     * \param Sw Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    static Scalar twoPhaseSatDKrw_dSw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        return (2.0/params.lambda() + 3)*std::pow(Sw, 2.0/params.lambda() + 2);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Sw Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    template <class FluidState>
    static Scalar krn(const Params &params, const FluidState &fs)
    { return twoPhaseSatKrn(params, 1.0 - fs.saturation(Traits::nPhaseIdx)); }

    static Scalar twoPhaseSatKrn(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        Scalar exponent = 2.0/params.lambda() + 1;
        Scalar Sn = 1. - Sw;
        return Sn*Sn*(1. - std::pow(Sw, exponent));
    }

    /*!
     * \brief The derivative of the relative permeability for the
     *        non-wetting phase in regard to the wetting saturation of
     *        the medium as implied by the Brooks-Corey
     *        parameterization.
     *
     * \param Sw Effective saturation of the wetting phase \f$[-]\f$
     * \param params The parameters of the capillary pressure curve
     *               (for Brooks-Corey: Entry pressure and shape factor)
     */
    static Scalar twoPhaseSatDKrn_dSw(const Params &params, Scalar Sw)
    {
        assert(0 <= Sw && Sw <= 1);

        Scalar alpha =
            1.0/params.lambda() + 1.0/2 -
            Sw*(1.0/params.lambda() + 1.0/2);
        return 2.0*(Sw - 1)*(1 + std::pow(Sw, 2.0/params.lambda())*alpha);
    }

};
} // namespace Opm

#endif // BROOKS_COREY_HPP
