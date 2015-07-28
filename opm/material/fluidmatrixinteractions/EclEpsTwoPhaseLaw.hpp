// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2015 by Andreas Lauser

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
 * \copydoc Opm::EclEpsTwoPhaseLaw
 */
#ifndef OPM_ECL_EPS_TWO_PHASE_LAW_HPP
#define OPM_ECL_EPS_TWO_PHASE_LAW_HPP

#include "EclEpsTwoPhaseLawParams.hpp"

#include <opm/material/fluidstates/SaturationOverlayFluidState.hpp>
#include <opm/material/common/ErrorMacros.hpp>
#include <opm/material/common/Exceptions.hpp>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief This material law takes a material law defined for effective
 *        saturations and converts it to a material law defined on absolute
 *        saturations.
 *
 * This class provides the endpoint scaling functionality used by the ECL reservoir
 * simulator.  The basic purpose of this class is the same as the one of \a EffToAbsLaw,
 * but it is quite a bit more complex.
 */
template <class EffLawT,
          class ParamsT = EclEpsTwoPhaseLawParams<EffLawT> >
class EclEpsTwoPhaseLaw : public EffLawT::Traits
{
    typedef EffLawT EffLaw;

public:
    typedef typename EffLaw::Traits Traits;
    typedef ParamsT Params;
    typedef typename EffLaw::Scalar Scalar;

    enum { wettingPhaseIdx = Traits::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = Traits::nonWettingPhaseIdx };

    //! The number of fluid phases
    static const int numPhases = EffLaw::numPhases;
    static_assert(numPhases == 2,
                  "The endpoint scaling applies to the nested twophase laws, not to "
                  "the threephase one!");

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = true;

    static_assert(EffLaw::implementsTwoPhaseApi,
                  "The material laws put into EclEpsTwoPhaseLaw must implement the "
                  "two-phase material law API!");

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = true;

    static_assert(EffLaw::implementsTwoPhaseSatApi,
                  "The material laws put into EclEpsTwoPhaseLaw must implement the "
                  "two-phase material law saturation API!");

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
     * \brief The capillary pressure-saturation curves depending on absolute saturations.
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
        OPM_THROW(NotAvailable,
                  "The capillaryPressures(fs) method is not yet implemented");
    }

    /*!
     * \brief The relative permeability-saturation curves depending on absolute saturations.
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
        OPM_THROW(NotAvailable,
                  "The pcnw(fs) method is not yet implemented");
    }

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     *
     * \param params A object that stores the appropriate coefficients
     *                for the respective law.
     *
     * \return Capillary pressure [Pa] calculated by specific
     *         constitutive relation (e.g. Brooks & Corey, van
     *         Genuchten, linear...)
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params &params, const FluidState &fs)
    {
        OPM_THROW(NotAvailable,
                  "The pcnw(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params &params, const Evaluation& SwAbs)
    {
        const Evaluation& SwEff = effectiveSaturationPc(params, SwAbs);
        const Evaluation& pcUnscaled = EffLaw::twoPhaseSatPcnw(params.effectiveLawParams(), SwEff);
        return scaleCapillaryPressure_(params, pcUnscaled);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnwInv(const Params &params, const Evaluation& pcnw)
    {
        Evaluation pcnwUnscaled = scaleCapillaryPressureInv_(params, pcnw);
        Evaluation SwUnscaled = EffLaw::twoPhaseSatPcnwInv(params.effectiveLawParams(), pcnwUnscaled);
        return effectiveSaturationPcInv(params, SwUnscaled);
    }

    /*!
     * \brief The saturation-capillary pressure curves.
     */
    template <class Container, class FluidState>
    static void saturations(Container &values, const Params &params, const FluidState &fs)
    {
        OPM_THROW(NotAvailable,
                  "The saturations(fs) method is not yet implemented");
    }

    /*!
     * \brief Calculate wetting liquid phase saturation given that
     *        the rest of the fluid state has been initialized
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params &params, const FluidState &fs)
    {
        OPM_THROW(NotAvailable,
                  "The Sw(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatSw(const Params &params, const Evaluation& pc)
    {
        OPM_THROW(NotAvailable,
                  "The twoPhaseSatSw(pc) method is not yet implemented");
    }

    /*!
     * \brief Calculate non-wetting liquid phase saturation given that
     *        the rest of the fluid state has been initialized
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params &params, const FluidState &fs)
    {
        OPM_THROW(NotAvailable,
                  "The Sn(pc) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatSn(const Params &params, const Evaluation& pc)
    {
        OPM_THROW(NotAvailable,
                  "The twoPhaseSatSn(pc) method is not yet implemented");
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the wetting phase calculated as implied by EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     *
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params &params, const FluidState &fs)
    {
        OPM_THROW(NotAvailable,
                  "The krw(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params &params, const Evaluation& Sw)
    {
        const Evaluation& Swe = effectiveSaturationKrw(params, Sw);
        const Evaluation& rawKrw = EffLaw::twoPhaseSatKrw(params.effectiveLawParams(), Swe);
        return scaleKrw_(params, rawKrw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrwInv(const Params &params, const Evaluation& krw)
    {
        Evaluation krwUnscaled = scaleKrwInv_(params, krw);
        Evaluation SwUnscaled = EffLaw::twoPhaseSatKrwInv(params.effectiveLawParams(), krwUnscaled);
        return effectiveSaturationKrwInv(params, SwUnscaled);
    }

    /*!
     * \brief The relative permeability of the non-wetting phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params &params, const FluidState &fs)
    {
        OPM_THROW(NotAvailable,
                  "The krn(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrn(const Params &params, const Evaluation& Sw)
    {
        const Evaluation& Swe = effectiveSaturationKrn(params, Sw);
        const Evaluation& rawKrn = EffLaw::twoPhaseSatKrn(params.effectiveLawParams(), Swe);
        return scaleKrn_(params, rawKrn);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrnInv(const Params &params, const Evaluation& krn)
    {
        Evaluation krnUnscaled = scaleKrnInv_(params, krn);
        Evaluation SwUnscaled = EffLaw::twoPhaseSatKrnInv(params.effectiveLawParams(), krnUnscaled);
        return effectiveSaturationKrnInv(params, SwUnscaled);
    }

    /*!
     * \brief Convert an absolute saturation to an effective one for capillary pressure.
     *
     * The effective saturation is then feed into the "raw" capillary pressure law.
     */
    template <class Evaluation>
    static Evaluation effectiveSaturationPc(const Params &params, const Evaluation& Sw)
    {
        if (!params.config().enableSatScaling())
            return Sw;

        // the saturations of capillary pressure are always scaled using two-point
        // scaling
        return scaleSatTwoPoint_(Sw,
                                 params.unscaledPoints().saturationPcPoints(),
                                 params.scaledPoints().saturationPcPoints());
    }

    template <class Evaluation>
    static Evaluation effectiveSaturationPcInv(const Params &params, const Evaluation& Sw)
    {
        if (!params.config().enableSatScaling())
            return Sw;

        // the saturations of capillary pressure are always scaled using two-point
        // scaling
        return scaleSatTwoPointInv_(Sw,
                                    params.unscaledPoints().saturationPcPoints(),
                                    params.scaledPoints().saturationPcPoints());
    }

    /*!
     * \brief Convert an absolute saturation to an effective one for the scaling of the
     *        relperm of the wetting phase.
     */
    template <class Evaluation>
    static Evaluation effectiveSaturationKrw(const Params &params, const Evaluation& Sw)
    {
        if (!params.config().enableSatScaling())
            return Sw;

        if (params.config().enableThreePointKrSatScaling()) {
            return scaleSatThreePoint_(Sw,
                                       params.unscaledPoints().saturationKrwPoints(),
                                       params.scaledPoints().saturationKrwPoints());
        }
        else { // two-point relperm saturation scaling
            return scaleSatTwoPoint_(Sw,
                                     params.unscaledPoints().saturationKrwPoints(),
                                     params.scaledPoints().saturationKrwPoints());
        }
    }

    template <class Evaluation>
    static Evaluation effectiveSaturationKrwInv(const Params &params, const Evaluation& Sw)
    {
        if (!params.config().enableSatScaling())
            return Sw;

        if (params.config().enableThreePointKrSatScaling()) {
            return scaleSatThreePointInv_(Sw,
                                          params.unscaledPoints().saturationKrwPoints(),
                                          params.scaledPoints().saturationKrwPoints());
        }
        else { // two-point relperm saturation scaling
            return scaleSatTwoPointInv_(Sw,
                                        params.unscaledPoints().saturationKrwPoints(),
                                        params.scaledPoints().saturationKrwPoints());
        }
    }

    /*!
     * \brief Convert an absolute saturation to an effective one for the scaling of the
     *        relperm of the non-wetting phase.
     */
    template <class Evaluation>
    static Evaluation effectiveSaturationKrn(const Params &params, const Evaluation& Sw)
    {
        if (!params.config().enableSatScaling())
            return Sw;

        if (params.config().enableThreePointKrSatScaling())
            return scaleSatThreePoint_(Sw,
                                       params.unscaledPoints().saturationKrnPoints(),
                                       params.scaledPoints().saturationKrnPoints());
        else // two-point relperm saturation scaling
            return scaleSatTwoPoint_(Sw,
                                     params.unscaledPoints().saturationKrnPoints(),
                                     params.scaledPoints().saturationKrnPoints());
    }


    template <class Evaluation>
    static Evaluation effectiveSaturationKrnInv(const Params &params, const Evaluation& Sw)
    {
        if (!params.config().enableSatScaling())
            return Sw;

        if (params.config().enableThreePointKrSatScaling()) {
            return scaleSatThreePointInv_(Sw,
                                          params.unscaledPoints().saturationKrnPoints(),
                                          params.scaledPoints().saturationKrnPoints());
        }
        else { // two-point relperm saturation scaling
            return scaleSatTwoPointInv_(Sw,
                                        params.unscaledPoints().saturationKrnPoints(),
                                        params.scaledPoints().saturationKrnPoints());
        }
    }

private:
    template <class Evaluation, class PointsContainer>
    static Evaluation scaleSatTwoPoint_(const Evaluation& S,
                                        const PointsContainer& unscaledSats,
                                        const PointsContainer& scaledSats)
    {
        return
            unscaledSats[0]
            +
            (S - scaledSats[0])*((unscaledSats[1] - unscaledSats[0])/(scaledSats[1] - scaledSats[0]));
    }

    template <class Evaluation, class PointsContainer>
    static Evaluation scaleSatTwoPointInv_(const Evaluation& S,
                                        const PointsContainer& unscaledSats,
                                        const PointsContainer& scaledSats)
    {
        return
            scaledSats[0]
            +
            (S - unscaledSats[0])*((scaledSats[1] - scaledSats[0])/(unscaledSats[1] - unscaledSats[0]));
    }

    template <class Evaluation, class PointsContainer>
    static Evaluation scaleSatThreePoint_(const Evaluation& S,
                                          const PointsContainer& unscaledSats,
                                          const PointsContainer& scaledSats)
    {
        if (unscaledSats[1] >= unscaledSats[2])
            return scaleSatTwoPoint_(S, unscaledSats, scaledSats);

        if (S < scaledSats[1])
            return
                unscaledSats[0]
                +
                (S - scaledSats[0])*((unscaledSats[1] - unscaledSats[0])/(scaledSats[1] - scaledSats[0]));
        else
            return
                unscaledSats[1]
                +
                (S - scaledSats[1])*((unscaledSats[2] - unscaledSats[1])/(scaledSats[2] - scaledSats[1]));
    }

    template <class Evaluation, class PointsContainer>
    static Evaluation scaleSatThreePointInv_(const Evaluation& S,
                                             const PointsContainer& unscaledSats,
                                             const PointsContainer& scaledSats)
    {
        if (unscaledSats[1] >= unscaledSats[2])
            return scaleSatTwoPointInv_(S, unscaledSats, scaledSats);

        if (S < unscaledSats[1])
            return
                scaledSats[0]
                +
                (S - unscaledSats[0])*((scaledSats[1] - scaledSats[0])/(unscaledSats[1] - unscaledSats[0]));
        else
            return
                scaledSats[1]
                +
                (S - unscaledSats[1])*((scaledSats[2] - scaledSats[1])/(unscaledSats[2] - unscaledSats[1]));
    }

    /*!
     * \brief Scale the capillary pressure according to the given parameters
     */
    template <class Evaluation>
    static Evaluation scaleCapillaryPressure_(const Params &params, const Evaluation& pc)
    {
        if (!params.config().enablePcScaling())
            return pc;

        Scalar alpha = params.scaledPoints().maxPcnw()/params.unscaledPoints().maxPcnw();
        return pc*alpha;
    }

    template <class Evaluation>
    static Evaluation scaleCapillaryPressureInv_(const Params &params, const Evaluation& pc)
    {
        if (!params.config().enablePcScaling())
            return pc;

        Scalar alpha = params.unscaledPoints().maxPcnw()/params.scaledPoints().maxPcnw();
        return pc/alpha;
    }

    /*!
     * \brief Scale the wetting phase relative permeability of a phase according to the given parameters
     */
    template <class Evaluation>
    static Evaluation scaleKrw_(const Params &params, const Evaluation& krw)
    {
        if (!params.config().enableKrwScaling())
            return krw;

        // TODO: three point krw y-scaling
        Scalar alpha = params.scaledPoints().maxKrw()/params.unscaledPoints().maxKrw();
        return krw*alpha;
    }

    template <class Evaluation>
    static Evaluation scaleKrwInv_(const Params &params, const Evaluation& krw)
    {
        if (!params.config().enableKrwScaling())
            return krw;

        Scalar alpha = params.unscaledPoints().maxKrw()/params.scaledPoints().maxKrw();
        return krw*alpha;
    }

    /*!
     * \brief Scale the non-wetting phase relative permeability of a phase according to the given parameters
     */
    template <class Evaluation>
    static Evaluation scaleKrn_(const Params &params, const Evaluation& krn)
    {
        if (!params.config().enableKrnScaling())
            return krn;

        //TODO: three point krn y-scaling
        Scalar alpha = params.scaledPoints().maxKrn()/params.unscaledPoints().maxKrn();
        return krn*alpha;
    }

    template <class Evaluation>
    static Evaluation scaleKrnInv_(const Params &params, const Evaluation& krn)
    {
        if (!params.config().enableKrnScaling())
            return krn;

        Scalar alpha = params.unscaledPoints().maxKrn()/params.scaledPoints().maxKrn();
        return krn*alpha;
    }
};
} // namespace Opm

#endif
