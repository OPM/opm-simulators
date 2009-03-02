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
 * \file
 *
 * \brief Implementation of the capillary pressure <-> saturation
 *         relation due to Brooks and Corey.
 */
#ifndef DUMUX_BROOKS_COREY_HH
#define DUMUX_BROOKS_COREY_HH

#include <dumux/new_material/brookscoreystate.hh>

#include <algorithm>

#include <math.h>
#include <assert.h>


namespace Dune
{
/*!\ingroup material
 *
 * \brief Implementation of the Brooks-Corey capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vince versa.
 *
 * \sa BrooksCorey, BrooksCoreyTwophase
 */
template <class StateT>
class BrooksCorey
{
public:
    typedef StateT State;
    typedef typename State::Scalar Scalar;

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     * The Brooks-Corey empirical capillary pressure <-> saturation
     * function is given by
     * \f[
     p_C = p_e\overline{S}_w^{-1/\alpha}
     \f]
     * \param Swe   Effective saturation of of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar pC(const State &state, Scalar Swe)
    {
        Api::require<Api::BrooksCoreyParams>(state);
        assert(0 <= Swe && Swe <= 1);

        return state.pe()*pow(Swe, -1.0/state.alpha());
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     \overline{S}_w = (\frac{p_C}{p_e})^{-\alpha}
     \f]
     *
     * \param pC Capillary pressure \f$\p_C\f$
     * \return The effective saturaion of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar Sw(const State &state, Scalar pC)
    {
        Api::require<Api::BrooksCoreyParams>(state);
        assert(pC >= 0);

        Scalar tmp = pow(pC/state.pe(), -state.alpha());
        return std::min(std::max(tmp, Scalar(0.0)), Scalar(1.0));
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     *
     * This is equivalent to
     * \f[
     \frac{\partial p_C}{\partial \overline{S}_w} =
     -\frac{p_e}{\alpha} \overline{S}_w^{-1/\alpha - 1}
     \f]
    */
    static Scalar dpC_dSw(const State &state, Scalar Swe)
    {
        Api::require<Api::BrooksCoreyParams>(state);
        assert(0 <= Swe && Swe <= 1);

        return - state.pe()/state.alpha() * pow(Swe, -1/state.alpha() - 1);
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const State &state, Scalar pC)
    {
        Api::require<Api::BrooksCoreyParams>(state);
        assert(pC >= 0);

        return -state.alpha()/state.pe() * pow(pC/state.pe(), - state.alpha() - 1);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param Sw_mob The mobile saturation of the wetting phase.
     */
    static Scalar krw(const State &state, Scalar Sw_mob)
    {
        Api::require<Api::BrooksCoreyParams>(state);
        assert(0 <= Sw_mob && Sw_mob <= 1);

        const Scalar krwLambda = 0.5;
        return pow(Sw_mob, (2. + 3*krwLambda) / krwLambda);
    };

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * \param Sw_mob The mobile saturation of the wetting phase.
     */
    static Scalar krn(const State &state, Scalar Sw_mob)
    {
        Api::require<Api::BrooksCoreyParams>(state);
        assert(0 <= Sw_mob && Sw_mob <= 1);

        const Scalar krnLambda = 0.5;
        Scalar exponent = (2. + krnLambda)/krnLambda;
        return pow(1. - Sw_mob, 2)*(1. - pow(Sw_mob, exponent));
    }
};
}

#endif // BROOKS_COREY_HH
