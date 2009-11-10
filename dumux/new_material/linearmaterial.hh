/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
 * \file linearmaterial.hh Implements a linear saturation-capillary
 *                    pressure relation
 */
#ifndef LINEAR_MATERIAL_HH
#define LINEAR_MATERIAL_HH

#include <dumux/new_material/linearmaterialcontext.hh>

#include <algorithm>

#include <math.h>
#include <assert.h>

namespace Dune
{
/*!
 * \ingroup material
 *
 * \brief Implements a linear saturation-capillary pressure relation
 *
 *
 * The entry pressure is reached at \f$S_w = 1\f$, the maximum
 * capillary pressure is observed at \f$S_w = 0\f$.
 *
 * \sa LinearMaterialContext
 */
template <class ContextT>
class LinearMaterial
{
public:
    typedef ContextT Context;
    typedef typename Context::Scalar Scalar;

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * This material law is linear:
     * \f[
     p_C = (1 - \overline{S}_w) (p_{C,max} - p_{C,entry}) + p_{C,entry}
     \f]
     *
     * \param Swe Effective saturation of of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar pC(const Context &context, Scalar Swe)
    {
        return (1 - Swe)*(context.maxPC() - context.entryPC()) + context.entryPC();
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     * This is the inverse of the capillary pressure-saturation curve:
     * \f[
     S_w = 1 - \frac{p_C - p_{C,entry}}{p_{C,max} - p_{C,entry}}
     \f]
     *
     * \param pC Capillary pressure \f$\p_C\f$
     * \return The effective saturaion of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar Sw(const Context &context, Scalar pC)
    {
        return 1 - (pC - context.entryPC())/(context.maxPC() - context.entryPC());
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     *
     * This is equivalent to
     * \f[
     \frac{\partial p_C}{\partial \overline{S}_w} =
     - (p_{C,max} - p_{C,min})
     \f]
    */
    static Scalar dpC_dSw(const Context &context, Scalar Swe)
    {
        return - (context.maxPC() - context.entryPC());
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Context &context, Scalar pC)
    {
        return - 1/(context.maxPC() - context.entryPC());
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param Sw_mob The mobile saturation of the wetting phase.
     */
    static Scalar krw(const Context &context, Scalar Sw_mob)
    {
        return std::min(Scalar(1),
                        std::max(Scalar(0),
                                 Sw_mob));
    };

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param Sw_mob The mobile saturation of the wetting phase.
     */
    static Scalar krn(const Context &context, Scalar Sw_mob)
    {
        return std::min(Scalar(1),
                        std::max(Scalar(0),
                                 1 - Sw_mob));
    }
};
}

#endif
