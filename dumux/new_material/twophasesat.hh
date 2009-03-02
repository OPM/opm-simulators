/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: and _at_ poware.org                                              *
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
 * \file TwophaseSat.hh Twophase saturation relations
 */
#ifndef TWOPHASE_SAT_HH
#define TWOPHASE_SAT_HH

#include <dumux/auxiliary/apis.hh>

#include <dumux/new_material/twophasesatstate.hh>

#include <stdlib.h>

namespace Dune
{
namespace Api
{
BEGIN_API_DEF(TwophaseSat)
{
    typedef typename Implementation::Scalar Scalar;
    typedef typename Implementation::State  State;

    const State *cTmpState = NULL;
    State *tmpState = NULL;
    Scalar tmp = 0.5;
    tmp = Implementation::Swe(*cTmpState, tmp);
    tmp = Implementation::Sw(*cTmpState, tmp);
    tmp = Implementation::Swmob(*cTmpState, tmp);
    tmp = Implementation::SwFromSwmob(*cTmpState, tmp);

    tmp = Implementation::dSw_dSwe(*cTmpState);
    tmp = Implementation::dSwe_dSw(*cTmpState);
}
END_API_DEF;
}; // namespace Api

/*!
 * \brief Saturation parameters for a two phase porpous medium.
 */
template <class StateT>
class TwophaseSat
{
public:
    typedef StateT State;
    typedef typename State::Scalar Scalar;

    /*!
     * \brief Returns effective wetting saturation given an absolute one.
     */
    static Scalar Swe(const State &state, Scalar Sw)
    {
        Api::require<Api::TwophaseSatParams>(state);

        return (Sw - state.Swr())/(1 - state.Swr());
    }

    /*!
     * \brief Returns absolute wetting saturation given an effective one.
     */
    static Scalar Sw(const State &state, Scalar Swe)
    {
        Api::require<Api::TwophaseSatParams>(state);
        return Swe*(1 - state.Swr()) + state.Swr();
    }

    /*!
     * \brief Convert an absolute wetting saturation to the mobile saturation.
     */
    static Scalar Swmob(const State &state, Scalar Sw)
    {
        Api::require<Api::TwophaseSatParams>(state);

        return (Sw - state.Swr() - state.Snr())/
            (1 - state.Swr() - state.Snr());
    }

    /*!
     * \brief Convert a mobile wetting saturation to the absolute one
     */
    static Scalar SwFromSwmob(const State &state, Scalar Swmob)
    {
        Api::require<Api::TwophaseSatParams>(state);

        return Swmob*(1 - state.Swr() - state.Snr())
            + state.Swr() + state.Snr();
    }


    /*!
     * \brief Derivative of the absolute wetting saturation regarding
     *        the effective saturation.
     */
    static Scalar dSw_dSwe(const State &state)
    {
        Api::require<Api::TwophaseSatParams>(state);

        return (1 - state.Swr());
    }

    /*!
     * \brief Derivative of the effective wetting saturation regarding
     *        the absolute saturation.
     */
    static Scalar dSwe_dSw(const State &state)
    {
        Api::require<Api::TwophaseSatParams>(state);

        return 1/(1 - state.Swr());
    }

};

/*!
 * \brief Identity twophase saturation relation
 */
template <class StateT>
class TwophaseSatId
{
public:
    typedef StateT State;
    typedef typename StateT::Scalar Scalar;

    /*!
     * \brief Returns effective wetting saturation given an absolute one.
     */
    static Scalar Swe(const State &state, Scalar Sw)
    {
        return Sw;
    }

    /*!
     * \brief Returns absolute wetting saturation given an effective one.
     */
    static Scalar Sw(const State &state, Scalar Swe)
    {
        return Swe;
    }

    /*!
     * \brief Convert an absolute wetting saturation to the mobile saturation.
     */
    static Scalar Swmob(const State &state, Scalar Sw)
    {
        return Sw;
    }

    /*!
     * \brief Convert a mobile wetting saturation to the absolute one
     */
    static Scalar SwFromSwmob(const State &state, Scalar Swmob)
    {
        return Swmob;
    }

    /*!
     * \brief Derivative of the absolute wetting saturation regarding
     *        the effective saturation.
     */
    static Scalar dSw_dSwe(const State &state)
    {
        return 1.0;
    }

    /*!
     * \brief Derivative of the effective wetting saturation regarding
     *        the absolute saturation.
     */
    static Scalar dSwe_dSw(const State &state)
    {
        return 1.0;
    }
};
};

#endif // TWOPHASE_SAT_HH
