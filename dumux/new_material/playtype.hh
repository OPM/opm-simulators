/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file PlayType.hh
 * \brief Template implementing play type capillary pressure hystersis.
 */
#ifndef PLAY_TYPE_HH
#define PLAY_TYPE_HH

#include <math.h>

#include <algorithm>

#include <dumux/auxiliary/apis.hh>
#include <dumux/new_material/playtypestate.hh>

namespace Dune
{
/*
 * \brief Implements the play-type hysteresis model.
 */
template <class StateT, class CapPressureModelT>
class PlayType
{
public:
    typedef StateT State;
    typedef CapPressureModelT CapPressure;
    typedef typename State::Scalar Scalar;

    typedef TwophaseSat<State> TwophaseSat;

    /*!
     * \brief Resets the hysteresis model to the
     *        initial state on the main drainage curve
     */
    static void reset(State &state)
    {
        Api::require<Api::PlayTypeState>(state);

        state.setSweRef(1.0);
        state.setImbib(false);
    };

    /*!
     * \brief Returns true iff the given absolute saturation
     *        is viable (i.e. larger than the residual
     *        saturation of the wetting phase but smaller
     *        than the maximum saturation of the current
     *        PISC)
     */
    static bool isValidSw(const State &state, Scalar Sw)
    {
        Api::require<Api::PlayTypeParams>(state);

        Scalar Swe = TwophaseSat::Swe(state, Sw);
        return 0 <= Swe && Swe <= 1;
    }

    /*!
     * \brief Set the current absolute saturation for the
     *        current timestep
     */
    static void updateState(State &state, Scalar Sw)
    {
        Api::require<Api::PlayTypeState>(state);

        Scalar Swe = TwophaseSat::Swe(state, Sw);
        if (state.isImbib()) {
            if (Swe <= state.SweRef() - state.deltaSwe()) {
                state.setImbib(false);
                state.setSweRef(Swe);
            }
            else if (Swe > state.SweRef())
                state.setSweRef(Swe);
        }
        else if (!state.isImbib()) {
            if (Swe >= state.SweRef() + state.deltaSwe()) {
                state.setImbib(true);
                state.setSweRef(Swe);
            }
            else if (Swe < state.SweRef()) {
                state.setSweRef(Swe);
            }
        };
    };

    /*!
     * \brief Returns the capillary pressure dependend on
     *        the absolute saturation of the wetting phase.
     */
    static Scalar pC(const State &state, Scalar Sw)
    {
        Api::require<Api::PlayTypeParams>(state);

        Scalar Swe = TwophaseSat::Swe(state, Sw);

        Scalar SwMdc_ref;
        Scalar SwMic_ref;
        if (!state.isImbib()) {
            SwMdc_ref = state.SweRef();
            SwMic_ref = state.SweRef() + state.deltaSwe();
        }
        else { // state.isImbib()
            SwMic_ref = state.SweRef();
            SwMdc_ref = state.SweRef() - state.deltaSwe();
        }

        if (Swe <= SwMdc_ref) {
            return CapPressure::pC(state.mdcParams(), Swe);
        }
        else if (Swe >= SwMic_ref) {
            Swe = std::min((Scalar) 1.0, Swe/(1 - state.Snre()));
            return CapPressure::pC(state.micParams(), Swe);
        }

        Scalar pos = (Swe - SwMdc_ref)/(SwMic_ref - SwMdc_ref);
        assert(0.0 <= pos && pos <= 1.0);

        Scalar pcDrain = CapPressure::pC(state.mdcParams(), SwMdc_ref);
        SwMic_ref = std::min((Scalar) 1.0, SwMic_ref/(1-state.Snre()));
        Scalar pcImbib = CapPressure::pC(state.micParams(), SwMic_ref);

        return pcDrain*(1 - pos) + pcImbib*pos;
    }

    /*!
     * \brief Returns the absolute saturation to a given on the
     *        capillary pressure.
     */
    static Scalar Sw(const State &state, Scalar pC)
    {
        Api::require<Api::PlayTypeParams>(state);

        assert(0);
        //DUNE_THROW(Dune::NotImplemented, "PlayType::Sw");
        return 0;
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the absolute saturation on the relevant
     *        hysteresis curve.
     *
     * \param Sw The absolute saturation where the derivative is
     *           supposed to be evaluated.
     */
    static Scalar dpC_dSw(const State &state, Scalar Sw)
    {
        Api::require<Api::PlayTypeParams>(state);

        Scalar Swe = TwophaseSat::Swe(state, Sw);

        Scalar SwMdc_ref;
        Scalar SwMic_ref;
        if (!state.isImbib()) {
            SwMdc_ref = state.SweRef();
            SwMic_ref = state.SweRef() + state.deltaSwe();
        }
        else { // state.isImbib()
            SwMic_ref = state.SweRef();
            SwMdc_ref = state.SweRef() - state.deltaSwe();
        }

        if (Swe <= SwMdc_ref) {
            return CapPressure::dpC_dSw(state.mdcParams(), Swe)*
                TwophaseSat::dSwe_dSw(state);
        }
        else if (Swe >= SwMic_ref) {
            Scalar tmp = std::min((Scalar) 1.0, Swe/(1 - state.Snre()));
            return CapPressure::dpC_dSw(state.micParams(), tmp)
                / (1 - state.Snre())
                * TwophaseSat::dSwe_dSw(state);
        }

        Scalar deltaPC = (CapPressure::pC(state.micParams(), SwMic_ref)
                          - CapPressure::pC(state.mdcParams(), SwMdc_ref));
        return deltaPC/(SwMic_ref - SwMdc_ref)
            * TwophaseSat::dSwe_dSw(state);
    };

    /*!
     * \brief Returns the partial derivative of absolute
     *        saturation to the capillary pressure on the relevant
     *        hysteresis curve.
     *
     * \param Sw The absolute saturation where the derivative is
     *           supposed to be evaluated.
     */
    static Scalar dSw_dpC(const State &state, Scalar pC)
    {
        Api::require<Api::PlayTypeParams>(state);

        assert(0);
        // DUNE_THROW(Dune::NotImplemented, "PlayType::dSw_dpC");
        return 0;
    }


    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium.
     */
    static Scalar krw(const State &state, Scalar Sw)
    {
        Api::require<Api::PlayTypeParams>(state);

        assert(Sw >= 0 && Sw <= 1);
        Scalar Swe = TwophaseSat::Swe(state, Sw);
        assert(0 <= Swe && Swe <= 1);

        if (Swe > state.SweRef()) {
            Swe = std::min((Scalar) 1.0,
                           Swe/(1 - state.Snre()));
            return CapPressure::krw(state.micParams(), Swe);
        }
        else {
            return CapPressure::krw(state.mdcParams(), Swe);
        }
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the params.
     */
    static Scalar krn(const State &state, Scalar Sw)
    {
        Api::require<Api::PlayTypeParams>(state);

        assert(Sw >= 0 && Sw <= 1);
        Scalar Swe = TwophaseSat::Swe(state, Sw);
        assert(0 <= Swe && Swe <= 1);

        if (Swe > state.SweRef()) {
            Swe = std::min((Scalar) 1.0,
                           Swe/(1 - state.Snre()));
            return CapPressure::krn(state.micParams(), Swe);
        }
        else {
            return CapPressure::krn(state.mdcParams(), Swe);
        }
    }
};
}; // namespace Dune

#endif // PLAY_TYPE_HH
