// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2012 by Andreas Lauser

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
 * \copydoc Opm::LinearMaterial
 */
#ifndef OPM_LINEAR_MATERIAL_HPP
#define OPM_LINEAR_MATERIAL_HPP

#include "LinearMaterialParams.hpp"

#include <opm/material/Valgrind.hpp>

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <algorithm>
#include <type_traits>

namespace Opm {

/*!
 * \ingroup material
 *
 * \brief Implements a linear saturation-capillary pressure relation
 *
 * Implements a linear saturation-capillary pressure relation for
 * M-phase fluid systems.
 *
 * \sa LinearMaterialParams
 */
template <class TraitsT, class ParamsT = LinearMaterialParams<TraitsT> >
class LinearMaterial : public TraitsT
{
public:
    typedef TraitsT Traits;
    typedef ParamsT Params;
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases
    static const int numPhases = Traits::numPhases;

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = (numPhases == 2);

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = (numPhases == 2);

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
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f[
     p_C = (1 - \overline{S}_w) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     \f]
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
        for (int phaseIdx = 0; phaseIdx < Traits::numPhases; ++phaseIdx) {
            Scalar S = state.saturation(phaseIdx);
            Valgrind::CheckDefined(S);

            values[phaseIdx] =
                S*params.pcMaxSat(phaseIdx) +
                (1.0 - S)*params.pcMinSat(phaseIdx);
        }
    }

    /*!
     * \brief The inverse of the capillary pressure
     */
    template <class ContainerT, class FluidState>
    static void saturations(ContainerT &values,
                            const Params &params,
                            const FluidState &state)
    {
        OPM_THROW(std::runtime_error, "Not implemented: LinearMaterial::saturations()");
    }

    /*!
     * \brief The relative permeability of all phases.
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params,
                                       const FluidState &state)
    {
        for (int phaseIdx = 0; phaseIdx < Traits::numPhases; ++phaseIdx) {
            Scalar S = state.saturation(phaseIdx);
            Valgrind::CheckDefined(S);

            values[phaseIdx] = std::max(std::min(S,1.0),0.0);
        }
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
        for (int pcPhaseIdx = 0; pcPhaseIdx < numPhases; ++pcPhaseIdx)
            values[pcPhaseIdx] = 0.0;

        values[satPhaseIdx] = params.pcMaxSat(satPhaseIdx) - params.pcMinSat(satPhaseIdx);
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
        for (int krPhaseIdx = 0; krPhaseIdx < numPhases; ++krPhaseIdx)
            values[krPhaseIdx] = 0.0;

        // -> linear relation between 0 and 1, else constant
        if (state.saturation(satPhaseIdx) >= 0 &&
            state.saturation(satPhaseIdx) <= 1)
        {
            values[satPhaseIdx] = 1.0;
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
     * \brief The difference between the pressures of the non-wetting and wetting phase.
     */
    template <class FluidState>
    static Scalar pcnw(const Params &params, const FluidState &fs)
    {
        Scalar S = fs.saturation(Traits::wPhaseIdx);
        Valgrind::CheckDefined(S);

        Scalar wPhasePressure =
            S*params.pcMaxSat(Traits::wPhaseIdx) +
            (1.0 - S)*params.pcMinSat(Traits::wPhaseIdx);

        S = fs.saturation(Traits::nPhaseIdx);
        Valgrind::CheckDefined(S);

        Scalar nPhasePressure =
            S*params.pcMaxSat(Traits::nPhaseIdx) +
            (1.0 - S)*params.pcMinSat(Traits::nPhaseIdx);

        return nPhasePressure - wPhasePressure;
    }


    template <class ScalarT = Scalar>
    static typename std::enable_if<Traits::numPhases == 2, ScalarT>::type
    twoPhaseSatPcnw(const Params &params, Scalar Sw)
    {
        Scalar wPhasePressure =
            Sw*params.pcMaxSat(Traits::wPhaseIdx) +
            (1.0 - Sw)*params.pcMinSat(Traits::wPhaseIdx);

        Scalar nPhasePressure =
            (1.0 - Sw)*params.pcMaxSat(Traits::nPhaseIdx) +
            Sw*params.pcMinSat(Traits::nPhaseIdx);

        return nPhasePressure - wPhasePressure;
    }

    /*!
     * \brief Calculate wetting phase saturation given that the rest
     *        of the fluid state has been initialized
     */
    template <class FluidState>
    static Scalar Sw(const Params &params, const FluidState &fs)
    { OPM_THROW(std::runtime_error, "Not implemented: Sw()"); }

    template <class ScalarT = Scalar>
    static typename std::enable_if<Traits::numPhases == 2, ScalarT>::type
    twoPhaseSatSw(const Params &params, Scalar Sw)
    { OPM_THROW(std::runtime_error, "Not implemented: twoPhaseSatSw()"); }

    /*!
     * \brief Calculate non-wetting liquid phase saturation given that
     *        the rest of the fluid state has been initialized
     */
    template <class FluidState>
    static Scalar Sn(const Params &params, const FluidState &fs)
    { OPM_THROW(std::runtime_error, "Not implemented: Sn()"); }

    template <class ScalarT = Scalar>
    static typename std::enable_if<Traits::numPhases == 2, ScalarT>::type
    twoPhaseSatSn(const Params &params, Scalar Sw)
    { OPM_THROW(std::runtime_error, "Not implemented: twoPhaseSatSn()"); }

    /*!
     * \brief Calculate gas phase saturation given that the rest of
     *        the fluid state has been initialized
     *
     * This method is only available for at least three fluid phases
     */
    template <class FluidState, class ScalarT = Scalar>
    static typename std::enable_if< (Traits::numPhases > 2), ScalarT>::type
    Sg(const Params &params, const FluidState &fs)
    { OPM_THROW(std::runtime_error, "Not implemented: Sg()"); }

    /*!
     * \brief The relative permability of the wetting phase
     */
    template <class FluidState>
    static Scalar krw(const Params &params, const FluidState &fs)
    { return std::max(0.0, std::min(1.0, fs.saturation(Traits::wPhaseIdx))); }

    template <class ScalarT = Scalar>
    static typename std::enable_if<Traits::numPhases == 2, ScalarT>::type
    twoPhaseSatKrw(const Params &params, Scalar Sw)
    { return std::max(0.0, std::min(1.0, Sw)); }

    /*!
     * \brief The relative permability of the liquid non-wetting phase
     */
    template <class FluidState>
    static Scalar krn(const Params &params, const FluidState &fs)
    { return std::max(0.0, std::min(1.0, fs.saturation(Traits::nPhaseIdx))); }

    template <class ScalarT = Scalar>
    static typename std::enable_if<Traits::numPhases == 2, ScalarT>::type
    twoPhaseSatKrn(const Params &params, Scalar Sw)
    { return std::max(0.0, std::min(1.0, 1 - Sw)); }

    /*!
     * \brief The relative permability of the gas phase
     *
     * This method is only available for at least three fluid phases
     */
    template <class FluidState, class ScalarT=Scalar>
    static typename std::enable_if< (Traits::numPhases > 2), ScalarT>::type
    krg(const Params &params, const FluidState &fs)
    { return std::max(0.0, std::min(1.0, fs.saturation(Traits::gPhaseIdx))); }

    /*!
     * \brief The difference between the pressures of the gas and the non-wetting phase.
     *
     * This method is only available for at least three fluid phases
     */
    template <class FluidState, class ScalarT=Scalar>
    static typename std::enable_if< (Traits::numPhases > 2), ScalarT>::type
    pcgn(const Params &params, const FluidState &fs)
    {
        Scalar S = fs.saturation(Traits::nPhaseIdx);
        Valgrind::CheckDefined(S);

        Scalar nPhasePressure =
            S*params.pcMaxSat(Traits::nPhaseIdx) +
            (1.0 - S)*params.pcMinSat(Traits::nPhaseIdx);

        S = fs.saturation(Traits::gPhaseIdx);
        Valgrind::CheckDefined(S);

        Scalar gPhasePressure =
            S*params.pcMaxSat(Traits::gPhaseIdx) +
            (1.0 - S)*params.pcMinSat(Traits::gPhaseIdx);

        return gPhasePressure - nPhasePressure;
    }
};
} // namespace Opm

#endif
