/*
  Copyright (C) 2008-2013 by Andreas Lauser
  Copyright (C) 2012 by Holger Class
  Copyright (C) 2012 by Vishal Jambhekar

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
 * \copydoc Opm::ThreePhaseParkerVanGenuchten
 */
#ifndef OPM_THREE_PHASE_PARKER_VAN_GENUCHTEN_HPP
#define OPM_THREE_PHASE_PARKER_VAN_GENUCHTEN_HPP

#include "ThreePhaseParkerVanGenuchtenParams.hpp"

#include <algorithm>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Implementation of three-phase capillary pressure and
 *        relative permeability relations proposed by Parker and van
 *        Genuchten.
 *
 * Reference: J.B. Kool, J.C. Parker, M.Th. van Genuchten: Parameter
 * Estimation for Unsaturated Flow and Transport Models -- A Review;
 * Journal of Hydrology, 91 (1987) 255-293
 */
template <class TraitsT,
          class ParamsT = ThreePhaseParkerVanGenuchtenParams<TraitsT> >
class ThreePhaseParkerVanGenuchten
{
public:
    static_assert(TraitsT::numPhases == 3,
                  "The number of phases considered by this capillary pressure "
                  "law is always three!");

    typedef TraitsT Traits;
    typedef ParamsT Params;
    typedef typename Traits::Scalar Scalar;

    static const int numPhases = 3;
    static const int wPhaseIdx = Traits::wPhaseIdx;
    static const int nPhaseIdx = Traits::nPhaseIdx;
    static const int gPhaseIdx = Traits::gPhaseIdx;

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = false;

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = false;

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
     * \brief Implements the three phase capillary pressure law
     *        proposed by Parker and van Genuchten.
     *
     * This material law is valid for three fluid phases and only
     * depends on the saturations.
     *
     * \param values Container for the return values
     * \param params Parameters
     * \param state The fluid state
     */
    template <class ContainerT, class FluidState>
    static void capillaryPressures(ContainerT &values,
                                   const Params &params,
                                   const FluidState &fluidState)
    {
        values[gPhaseIdx] = pcgn(params, fluidState);
        values[nPhaseIdx] = 0;
        values[wPhaseIdx] = - pcnw(params, fluidState);
    }

    /*!
     * \brief Capillary pressure between the gas and the non-wetting
     *        liquid (i.e., oil) phase.
     *
     * This is defined as
     * \f[
     * p_{c,gn} = p_g - p_n
     * \f]
     */
    template <class FluidState>
    static Scalar pcgn(const Params &params,
                       const FluidState &fluidState)
    {
        Scalar PC_VG_REG = 0.01;

        // sum of liquid saturations
        Scalar St =
            fluidState.saturation(wPhaseIdx)
            + fluidState.saturation(nPhaseIdx);

        Scalar Se = (St - params.Swrx())/(1. - params.Swrx());

        // regularization
        if (Se < 0.0)
            Se=0.0;
        if (Se > 1.0)
            Se=1.0;

        if (Se>PC_VG_REG && Se<1-PC_VG_REG)
        {
            Scalar x = std::pow(Se,-1/params.vgM()) - 1;
            return std::pow(x, 1 - params.vgM())/params.vgAlpha();
        }

        // value and derivative at regularization point
        Scalar Se_regu;
        if (Se<=PC_VG_REG)
            Se_regu = PC_VG_REG;
        else
            Se_regu = 1-PC_VG_REG;
        Scalar x = std::pow(Se_regu,-1/params.vgM())-1;
        Scalar pc = std::pow(x, 1/params.vgN())/params.vgAlpha();
        Scalar pc_prime =
            std::pow(x, 1/params.vgN()-1)
            * std::pow(Se_regu,-1/params.vgM()-1)
            / (-params.vgM())
            / params.vgAlpha()
            / (1 - params.Sgr() - params.Swrx())
            / params.vgN();

        // evaluate tangential
        return ((Se-Se_regu)*pc_prime + pc)/params.betaGN();
    }

    /*!
     * \brief Capillary pressure between the non-wetting liquid (i.e.,
     *        oil) and the wetting liquid (i.e., water) phase.
     *
     * This is defined as
     * \f[
     * p_{c,nw} = p_n - p_w
     * \f]
     */
    template <class FluidState>
    static Scalar pcnw(const Params &params,
                       const FluidState &fluidState)
    {
        Scalar Sw = fluidState.saturation(wPhaseIdx);
        Scalar Se = (Sw-params.Swr())/(1.-params.Snr());

        Scalar PC_VG_REG = 0.01;

        // regularization
        if (Se<0.0)
            Se=0.0;
        if (Se>1.0)
            Se=1.0;

        if (Se>PC_VG_REG && Se<1-PC_VG_REG) {
            Scalar x = std::pow(Se,-1/params.vgM()) - 1.0;
            x = std::pow(x, 1 - params.vgM());
            return x/params.vgAlpha();
        }

        // value and derivative at regularization point
        Scalar Se_regu;
        if (Se<=PC_VG_REG)
            Se_regu = PC_VG_REG;
        else
            Se_regu = 1.0 - PC_VG_REG;

        Scalar x = std::pow(Se_regu,-1/params.vgM())-1;
        Scalar pc = std::pow(x, 1/params.vgN())/params.vgAlpha();
        Scalar pc_prime =
            std::pow(x,1/params.vgN()-1)
            * std::pow(Se_regu, -1.0/params.vgM() - 1)
            / (-params.vgM())
            / params.vgAlpha()
            / (1-params.Snr()-params.Swr())
            / params.vgN();

        // evaluate tangential
        return ((Se-Se_regu)*pc_prime + pc)/params.betaNW();
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     */
    template <class ContainerT, class FluidState>
    static void saturations(ContainerT &values,
                            const Params &params,
                            const FluidState &fluidState)
    { OPM_THROW(std::logic_error, "Not implemented: inverse capillary pressures"); }

    /*!
     * \brief The saturation of the gas phase.
     */
    template <class FluidState>
    static Scalar Sg(const Params &params,
                     const FluidState &fluidState)
    { OPM_THROW(std::logic_error, "Not implemented: Sg()"); }

    /*!
     * \brief The saturation of the non-wetting (i.e., oil) phase.
     */
    template <class FluidState>
    static Scalar Sn(const Params &params,
                     const FluidState &fluidState)
    { OPM_THROW(std::logic_error, "Not implemented: Sn()"); }

    /*!
     * \brief The saturation of the wetting (i.e., water) phase.
     */
    template <class FluidState>
    static Scalar Sw(const Params &params,
                     const FluidState &fluidState)
    { OPM_THROW(std::logic_error, "Not implemented: Sw()"); }

    /*!
     * \brief The relative permeability of all phases.
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT &values,
                                       const Params &params,
                                       const FluidState &fluidState)
    {
        values[wPhaseIdx] = krw(params, fluidState);
        values[nPhaseIdx] = krn(params, fluidState);
        values[gPhaseIdx] = krg(params, fluidState);
    }

    /*!
     * \brief The relative permeability for the wetting phase of the
     *        medium implied by van Genuchten's parameterization.
     *
     * The permeability of water in a 3p system equals the standard 2p description.
     * (see p61. in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.)
     */
    template <class FluidState>
    static Scalar krw(const Params &params,
                      const FluidState &fluidState)
    {

        // transformation to effective saturation
        Scalar Se = (fluidState.saturation(wPhaseIdx) - params.Swr()) / (1-params.Swr());

        // regularization
        if(Se > 1.0) return 1.;
        if(Se < 0.0) return 0.;

        Scalar r = 1. - std::pow(1 - std::pow(Se, 1/params.vgM()), params.vgM());
        return std::sqrt(Se)*r*r;
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        due to the model of Parker et al. (1987).
     *
     * See model 7 of "Comparison of the Three-Phase Oil Relative
     * Permeability Models" M. Delshad and G. A. Pope, Transport in
     * Porous Media 4 (1989), 59-83; or -- more comprehensively --
     * "Estimation of primary drainage three-phase relative
     * permeability for organic liquid transport in the vadose zone",
     * L. I. Oliveira, A. H. Demond, Journal of Contaminant Hydrology
     * 66 (2003), 261-285
     */
    template <class FluidState>
    static Scalar krn(const Params &params,
                      const FluidState &fluidState)
    {
        Scalar Sn = fluidState.saturation(nPhaseIdx);
        Scalar Sw = fluidState.saturation(wPhaseIdx);
        Scalar Swe = std::min((Sw - params.Swr()) / (1 - params.Swr()), 1.);
        Scalar Ste = std::min((Sw + Sn - params.Swr()) / (1 - params.Swr()), 1.);

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
            // regard Snr in the permeability of the non-wetting
            // phase, see Helmig1997
            Scalar resIncluded =
                std::max(std::min(Sw - params.Snr() / (1-params.Swr()), 1.0),
                         0.0);
            krn_ *= std::sqrt(resIncluded );
        }
        else
            krn_ *= std::sqrt(Sn / (1 - params.Swr()));

        return krn_;
    }


    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the medium implied by van Genuchten's
     *        parameterization.
     *
     * The permeability of gas in a three-phase system equals the
     * standard two-phase description. (see p61. of "Comparison of the
     * Three-Phase Oil Relative Permeability Models" M.  Delshad and
     * G. A. Pope, Transport in Porous Media 4 (1989), 59-83.)
     */
    template <class FluidState>
    static Scalar krg(const Params &params,
                      const FluidState &fluidState)
    {
        Scalar Sg = fluidState.saturation(gPhaseIdx);
        Scalar Se = std::min(((1-Sg) - params.Sgr()) / (1 - params.Sgr()), 1.);

        // regularization
        if(Se > 1.0)
            return 0.0;
        if(Se < 0.0)
            return 1.0;

        Scalar scaleFactor = 1.;
        if (Sg<=0.1) {
            scaleFactor = (Sg - params.Sgr())/(0.1 - params.Sgr());
            if (scaleFactor < 0.)
                scaleFactor = 0.;
        }

        return scaleFactor
            * std::pow(1 - Se, 1.0/3.)
            * std::pow(1 - std::pow(Se, 1/params.vgM()), 2*params.vgM());
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
        OPM_THROW(std::logic_error,
                  "Not implemented: dCapillaryPressures_dSaturation()");
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
        OPM_THROW(std::logic_error,
                  "Not implemented: dRelativePermeabilities_dSaturation()");
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
};
} // namespace Opm

#endif
