// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Holger Class                                      *
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
 * \file parkerVanGen3pparams.hh Implementation of van Genuchten's capillary
 *                       pressure <-> saturation relation
 */
#ifndef PARKERVANGEN_3P_HH
#define PARKERVANGEN_3P_HH


#include "parkerVanGen3pparams.hh"

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
template <class ScalarT, class ParamsT = ParkerVanGen3PParams<ScalarT> >
class ParkerVanGen3P
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
        DUNE_THROW(Dune::NotImplemented, "Capillary pressures for three phases is not so simple! Use pCGN, pCNW, and pcGW");
    }

    static Scalar pCGW(const Params &params, Scalar Sw)
    {
    /*
         Sw = wetting phase saturation, or,
              sum of wetting phase saturations
         alpha : VanGenuchten-alpha
    this function is just copied from MUFTE/pml/constrel3p3cni.c
    that is why variable names do not yet fulfill Dumux rules, TODO Change */

    Scalar r,Se,x,vg_m;
    Scalar pc,pc_prime,Se_regu;
    Scalar PC_VG_REG = 0.01;

    Se   = (Sw-params.Swr())/(1.-params.Sgr());

    /* Snr  = 0.0;   test version   */

    /* regularization */
    if (Se<0.0) Se=0.0;
    if (Se>1.0) Se=1.0;
    vg_m = 1.-1./params.vgN();

        if (Se>PC_VG_REG && Se<1-PC_VG_REG)
        {
            r = std::pow(Se,-1/vg_m);
            x = r-1;
            vg_m = 1-vg_m;
            x = std::pow(x,vg_m);
            r = x/params.vgAlpha();
            return(r);
        }
        else
        {
            /* value and derivative at regularization point */
            if (Se<=PC_VG_REG) Se_regu = PC_VG_REG; else Se_regu = 1-PC_VG_REG;
            pc       = std::pow(std::pow(Se_regu,-1/vg_m)-1,1/params.vgN())/params.vgAlpha();
            pc_prime = std::pow(std::pow(Se_regu,-1/vg_m)-1,1/params.vgN()-1)*std::pow(Se_regu,-1/vg_m-1)*(-1/vg_m)/params.vgAlpha()/(1-params.Sgr()-params.Swr())/params.vgN();

            /* evaluate tangential */
            r        = (Se-Se_regu)*pc_prime+pc;
            return(r);
        }
    }

    static Scalar pCNW(const Params &params, Scalar Sw)
    {
    /*
         Sw = wetting phase saturation, or,
              sum of wetting phase saturations
         alpha : VanGenuchten-alpha
    this function is just copied from MUFTE/pml/constrel3p3cni.c
    that is why variable names do not yet fulfill Dumux rules, TODO Change */

    Scalar r,Se,x,vg_m;
    Scalar pc,pc_prime,Se_regu;
    Scalar PC_VG_REG = 0.01;

    Se   = (Sw-params.Swr())/(1.-params.Snr());

    /* Snr  = 0.0;   test version   */

    /* regularization */
    if (Se<0.0) Se=0.0;
    if (Se>1.0) Se=1.0;
    vg_m = 1.-1./params.vgN();

        if (Se>PC_VG_REG && Se<1-PC_VG_REG)
        {
            r = std::pow(Se,-1/vg_m);
            x = r-1;
            vg_m = 1-vg_m;
            x = std::pow(x,vg_m);
            r = x/params.vgAlpha();
            return(r);
        }
        else
        {
            /* value and derivative at regularization point */
            if (Se<=PC_VG_REG) Se_regu = PC_VG_REG; else Se_regu = 1-PC_VG_REG;
            pc       = std::pow(std::pow(Se_regu,-1/vg_m)-1,1/params.vgN())/params.vgAlpha();
            pc_prime = std::pow(std::pow(Se_regu,-1/vg_m)-1,1/params.vgN()-1)*std::pow(Se_regu,-1/vg_m-1)*(-1/vg_m)/params.vgAlpha()/(1-params.Snr()-params.Swr())/params.vgN();

            /* evaluate tangential */
            r        = (Se-Se_regu)*pc_prime+pc;
            return(r);
        }
    }

    static Scalar pCGN(const Params &params, Scalar St)
    {
    /*
         St = sum of wetting (liquid) phase saturations
         alpha : VanGenuchten-alpha
    this function is just copied from MUFTE/pml/constrel3p3cni.c
    that is why variable names do not yet fulfill Dumux rules, TODO Change */

    Scalar r,Se,x,vg_m;
    Scalar pc,pc_prime,Se_regu;
    Scalar PC_VG_REG = 0.01;

    Se   = (St-params.Swrx())/(1.-params.Swrx());

    /* Snr  = 0.0;   test version   */

    /* regularization */
    if (Se<0.0) Se=0.0;
    if (Se>1.0) Se=1.0;
    vg_m = 1.-1./params.vgN();

        if (Se>PC_VG_REG && Se<1-PC_VG_REG)
        {
            r = std::pow(Se,-1/vg_m);
            x = r-1;
            vg_m = 1-vg_m;
            x = std::pow(x,vg_m);
            r = x/params.vgAlpha();
            return(r);
        }
        else
        {
            /* value and derivative at regularization point */
            if (Se<=PC_VG_REG) Se_regu = PC_VG_REG; else Se_regu = 1-PC_VG_REG;
            pc       = std::pow(std::pow(Se_regu,-1/vg_m)-1,1/params.vgN())/params.vgAlpha();
            pc_prime = std::pow(std::pow(Se_regu,-1/vg_m)-1,1/params.vgN()-1)*std::pow(Se_regu,-1/vg_m-1)*(-1/vg_m)/params.vgAlpha()/(1-params.Sgr()-params.Swrx())/params.vgN();

            /* evaluate tangential */
            r        = (Se-Se_regu)*pc_prime+pc;
            return(r);
        }
    }

    static Scalar pCAlpha(const Params &params, Scalar Sn)
    {
        /* continuous transition to zero */
        Scalar alpha,Sne;

        Sne=Sn;
        /* regularization */
        if (Sne<=0.001) Sne=0.0;
        if (Sne>=1.0) Sne=1.0;

        if (Sne>params.Snr()) alpha = 1.0;
        else
        {
         if (params.Snr()>=0.001) alpha = Sne/params.Snr();
         else          alpha = 0.0;
        }
        return(alpha);
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        DUNE_THROW(Dune::NotImplemented, "Sw(pc) for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     *
    */
    static Scalar dpC_dSw(const Params &params, Scalar Sw)
    {
        DUNE_THROW(Dune::NotImplemented, "dpC/dSw for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        DUNE_THROW(Dune::NotImplemented, "dSw/dpC for three phases not implemented! Do it yourself!");
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
     * \param Sn saturation of the NAPL phase.
     * \param Sg saturation of the gas phase.
     * \param saturation saturation of the water phase.
     * \param params Array of parameters.
     */
    static Scalar krw(const Params &params,  Scalar saturation, Scalar Sn, Scalar Sg)
    {

        //transformation to effective saturation
        Scalar Se = (saturation - params.Swr()) / (1-params.Swr());

        /* regularization */
        if(Se > 1.0) return 1.;
        if(Se < 0.0) return 0.;

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
     * \param Sw saturation of the water phase.
     * \param Sg saturation of the gas phase.
     * \param saturation saturation of the NAPL phase.
     * \param params Array of parameters.
     */
    static Scalar krn(const Params &params, Scalar Sw, Scalar saturation, Scalar Sg)
    {

        Scalar Swe = std::min((Sw - params.Swr()) / (1 - params.Swr()), 1.);
        Scalar Ste = std::min((Sw +  saturation - params.Swr()) / (1 - params.Swr()), 1.);

        // regularization
        if(Swe <= 0.0) Swe = 0.;
        if(Ste <= 0.0) Ste = 0.;
        if(Ste - Swe <= 0.0) return 0.;

        Scalar krn_;
        krn_ = std::pow(1 - std::pow(Swe, 1/params.vgM()), params.vgM());
        krn_ -= std::pow(1 - std::pow(Ste, 1/params.vgM()), params.vgM());
        krn_ *= krn_;

        if (params.krRegardsSnr())
        {
            // regard Snr in the permeability of the n-phase, see Helmig1997
            Scalar resIncluded = std::max(std::min((saturation - params.Snr()/ (1-params.Swr())), 1.), 0.);
            krn_ *= std::sqrt(resIncluded );
        }
        else
            krn_ *= std::sqrt(saturation / (1 - params.Swr()));   // Hint: (Ste - Swe) = Sn / (1-Srw)


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
     * \param Sw saturation of the water phase.
     * \param Sn saturation of the NAPL phase.
     * \param saturation saturation of the gas phase.
     * \param params Array of parameters.
     */
    static Scalar krg(const Params &params, Scalar Sw, Scalar Sn, Scalar saturation)
    {

        // Se = (Sw+Sn - Sgr)/(1-Sgr)
        Scalar Se = std::min(((1-saturation) - params.Sgr()) / (1 - params.Sgr()), 1.);


        /* regularization */
        if(Se > 1.0) return 0.0;
        if(Se < 0.0) return 1.0;
        Scalar scalFact = 1.;
        if (saturation<=0.1)
        {
          scalFact = (saturation - params.Sgr())/(0.1 - params.Sgr());
          if (scalFact < 0.) scalFact = 0.;
        }

        Scalar result = scalFact * std::pow(1 - Se, 1.0/3.) * std::pow(1 - std::pow(Se, 1/params.vgM()), 2*params.vgM());

        return result;
    };

    /*!
     * \brief The relative permeability for a phase.
     * \param Sw saturation of the water phase.
     * \param Sg saturation of the gas phase.
     * \param Sn saturation of the NAPL phase.
     * \param params Array of parameters.
     * \param phase indicator, The saturation of all phases.
     */
    static Scalar kr(const Params &params, const int phase, const Scalar Sw, const Scalar Sn, const Scalar Sg)
    {
        switch (phase)
        {
        case 0:
            return krw(params, Sw, Sn, Sg);
            break;
        case 1:
            return krn(params, Sw, Sn, Sg);
            break;
        case 2:
            return krg(params, Sw, Sn, Sg);
            break;
        }
        return 0;
    };

   /*
    * \brief the basis for calculating adsorbed NAPL in storage term
    * \param bulk density of porous medium, adsorption coefficient
    */
   static Scalar bulkDensTimesAdsorpCoeff (const Params &params)
   {
      return params.rhoBulk() * params.KdNAPL();
   }
};
}

#endif
