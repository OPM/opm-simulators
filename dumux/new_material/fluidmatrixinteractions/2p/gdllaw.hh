/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
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
 * \file gdllaw.hh
 *
 * \brief Implements the experimental capillary-pressure saturation
 *        relations for gas diffusion layers in PEM fuel cells
 */
#ifndef DUMUX_GDL_LAW_HH
#define DUMUX_GDL_LAW_HH

#include "gdllawparams.hh"

namespace Dumux
{
/*!
 * \ingroup material
 *
 * \brief Implements the experimental capillary-pressure saturation
 *        relations for gas diffusion layers in PEM fuel cells
 *
 * \sa GdlLawParams
 */
template <class ScalarT, class ParamsT = GdlLawParams<ScalarT> >
class GdlLaw
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * \param Sw Absolute saturation of of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar pC(const Params &params, Scalar Sw)
    {
        int idx = findIdxSw_(params, Sw);

        // linear interpolation
        Scalar alpha = (Sw - params.Sw(idx))/(params.Sw(idx + 1) - params.Sw(idx));
        return alpha*params.pC(idx + 1) + (1-alpha)*params.pC(idx);
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        int idx = findIdxPc_(params, pC);

        // linear interpolation
        Scalar alpha = (pC - params.pC(idx))/(params.pC(idx + 1) - params.pC(idx));
        return alpha*params.Sw(idx + 1) + (1 - alpha)*params.Sw(idx);
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the absolute saturation.
    */
    static Scalar dpC_dSw(const Params &params, Scalar Sw)
    {
        int idx = findIdxSw_(params, Sw);

        return 
            (params.pC(idx + 1) - params.pC(idx)) /
            (params.Sw(idx + 1) - params.Sw(idx));
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        int idx = findIdxPc_(params, pC);

        return 
            (params.Sw(idx + 1) - params.Sw(idx)) /
            (params.pC(idx + 1) - params.pC(idx));
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param Sw The absolute saturation of the wetting phase.
     */
    static Scalar krw(const Params &params, Scalar Sw)
    {
        return Sw;
    };

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param Sw The absolute saturation of the wetting phase.
     */
    static Scalar krn(const Params &params, Scalar Sw)
    {
        Scalar Sn = 1 - Sw;
        return Sn;
    }

private:
    // returns the index of a sampling point depending on the
    // capillary pressure
    static int findIdxPc_(const Params &params, Scalar pC)
    {
        // do an interval bisection. we assume that the capillary
        // pressure decreases monotonously.
        int lowIdx = 0;
        int highIdx = params.numSamples() - 1;
        while (lowIdx < highIdx - 1) {
            int idx = lowIdx + (highIdx - lowIdx)/2;
            if (pC < params.pC(idx))
                lowIdx = idx;
            else
                highIdx = idx;
        };
        return lowIdx;
    };

    // returns the index of a sampling point depending on the
    // absolute wetting saturation
    static int findIdxSw_(const Params &params, Scalar Sw)
    {
        // do an interval bisection. we assume that the absolute
        // wetting saturation increases monotonously
        int lowIdx = 0;
        int highIdx = params.numSamples() - 1;
        while (lowIdx < highIdx - 1) {
            int idx = lowIdx + (highIdx - lowIdx)/2;
            if (params.Sw(idx) < Sw)
                lowIdx = idx;
            else
                highIdx = idx;
        };
        return lowIdx;
    };
};
}

#endif
