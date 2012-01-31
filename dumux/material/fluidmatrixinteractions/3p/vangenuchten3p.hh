// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Jochen Fritz, Benjamin Faigle                     *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \file VanGenuchten.hh Implementation of van Genuchten's capillary
 *                       pressure <-> saturation relation
 */
#ifndef VAN_GENUCHTEN_3P_HH
#define VAN_GENUCHTEN_3P_HH

#include "vangenuchtenparams3p.hh"

#include <algorithm>


namespace Dumux
{
/*!
 * \ingroup material
 *
 * \brief Implementation of van Genuchten's capillary pressure <->
 *        saturation relation. This class bundles the "raw" curves
 *        as static members and doesn't concern itself converting
 *        absolute to effective saturations and vince versa.
 *
 * \sa VanGenuchten, VanGenuchtenThreephase
 */
template <class ScalarT, class ParamsT = VanGenuchtenParams3p<ScalarT> >
class VanGenuchten3P
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     */
    static Scalar pC(const Params &params, Scalar Sw)
    {
        DUNE_THROW(Dune::NotImplemented, "Capillary pressures for three phases not implemented! Do it yourself!");
   }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        DUNE_THROW(Dune::NotImplemented, "Capillary pressures for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     *
    */
    static Scalar dpC_dSw(const Params &params, Scalar Sw)
    {
        DUNE_THROW(Dune::NotImplemented, "Capillary pressures for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        DUNE_THROW(Dune::NotImplemented, "Capillary pressures for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by van Genuchten's
     *        parameterization.
     *
     * The permeability of water in a 3p system equals the standard 2p description.
     * (see p61. in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.)
     *
     * \param The mobile saturation of all phases. (Sw used)
     */
    static Scalar krw(const Params &params,  const Dune::FieldVector<Scalar, 3> saturation)
    {
        assert(0 <= saturation[0] && saturation[0] <= 1);

        //transformation to effective saturation
        Scalar Se = (saturation[0] - params.Swr()) / (1-params.Swr());

        /* regularization */
        if(Se > 1.) return 1.;
        if(Se < machineEps_) return 0.;

        Scalar r = 1. - std::pow(1 - std::pow(Se, 1/params.vgM()), params.vgM());
        return std::sqrt(Se)*r*r;
    };

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        after the Model of Parker et al. (1987).
     *
     * See model 7 in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.
     * or more comprehensive in
     * "Estimation of primary drainage three-phase relative permeability for organic
     * liquid transport in the vadose zone", Leonardo I. Oliveira, Avery H. Demond,
     * Journal of Contaminant Hydrology 66 (2003), 261-285
     *
     *
     * \param The mobile saturation of all phases. (Sw and Sn used)
     */
    static Scalar krn(const Params &params, const Dune::FieldVector<Scalar, 3> saturation)
    {
        assert(0 <= saturation[0] && saturation[0] <= 1);   //check Sw
        assert(0 <= saturation[1] && saturation[1] <= 1);   //check Sn


        Scalar Swe = std::min((saturation[0] - params.Swr()) / (1 - params.Swr()), 1.);
        Scalar Ste = std::min((saturation[0] +  saturation[1] - params.Swr()) / (1 - params.Swr()), 1.);

        // regularization
        if(Swe < machineEps_) Swe = 0.;
        if(Ste < machineEps_) Ste = 0.;
        if(Ste - Swe < machineEps_) return 0.;

        Scalar krn_;
        krn_ = std::pow(1 - std::pow(Swe, 1/params.vgM()), params.vgM());
        krn_ -= std::pow(1 - std::pow(Ste, 1/params.vgM()), params.vgM());
        krn_ *= krn_;

        if (params.krRegardsSnr())
        {
            // regard Snr in the permeability of the n-phase, see Helmig1997
            double resIncluded = std::max(std::min((saturation[1] - params.Snr()/ (1-params.Swr())), 1.), 0.);
            krn_ *= std::sqrt(resIncluded );
        }
        else
            krn_ *= std::sqrt(saturation[1] / (1 - params.Swr()));   // Hint: (Ste - Swe) = Sn / (1-Srw)


        return krn_;
    };


    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * The permeability of gas in a 3p system equals the standard 2p description.
     * (see p61. in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.)
     *
     * \param The mobile saturation of all phases. (Sg used)
     */
    static Scalar krg(const Params &params, const Dune::FieldVector<Scalar, 3> saturation)
    {
        assert(0 <= saturation[2] && saturation[2] <= 1);

        // Se = (Sw+Sn - Sgr)/(1-Sgr)
        Scalar Se = std::min(((1-saturation[2]) - params.Sgr()) / (1 - params.Sgr()), 1.);

        /* regularization */
        if(Se > 1.) return 0.;
        if(Se < machineEps_) return 1.;

        return
            std::pow(1 - Se, 1.0/3) *
            std::pow(1 - std::pow(Se, 1/params.vgM()), 2*params.vgM());
    };

    /*!
     * \brief The relative permeability for a phase.
     *
     * \param Phase indicator, The saturation of all phases.
     */
    static Scalar kr(const int phase, const Params &params, const Dune::FieldVector<Scalar, 3> saturation)
    {
        switch (phase)
        {
        case 0:
            return krw(params, saturation);
            break;
        case 1:
            return krn(params, saturation);
            break;
        case 2:
            return krg(params, saturation);
            break;
        }
        return 0;
    };

protected:
    static const double machineEps_ = 1e-15;
};
}

#endif
