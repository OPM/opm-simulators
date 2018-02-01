// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::ThreePhaseParkerVanGenuchten
 */
#ifndef OPM_THREE_PHASE_PARKER_VAN_GENUCHTEN_HPP
#define OPM_THREE_PHASE_PARKER_VAN_GENUCHTEN_HPP

#include "ThreePhaseParkerVanGenuchtenParams.hpp"

#include <opm/material/common/MathToolbox.hpp>

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
    static const int wettingPhaseIdx = Traits::wettingPhaseIdx;
    static const int nonWettingPhaseIdx = Traits::nonWettingPhaseIdx;
    static const int gasPhaseIdx = Traits::gasPhaseIdx;

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
    static void capillaryPressures(ContainerT& values,
                                   const Params& params,
                                   const FluidState& fluidState)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[gasPhaseIdx] = pcgn<FluidState, Evaluation>(params, fluidState);
        values[nonWettingPhaseIdx] = 0;
        values[wettingPhaseIdx] = - pcnw<FluidState, Evaluation>(params, fluidState);
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
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcgn(const Params& params, const FluidState& fluidState)
    {
        Scalar PC_VG_REG = 0.01;

        // sum of liquid saturations
        const auto& St =
            Opm::decay<Evaluation>(fluidState.saturation(wettingPhaseIdx))
            + Opm::decay<Evaluation>(fluidState.saturation(nonWettingPhaseIdx));

        Evaluation Se = (St - params.Swrx())/(1. - params.Swrx());

        // regularization
        if (Se < 0.0)
            Se=0.0;
        if (Se > 1.0)
            Se=1.0;

        if (Se>PC_VG_REG && Se<1-PC_VG_REG)
        {
            const Evaluation& x = Opm::pow(Se,-1/params.vgM()) - 1;
            return Opm::pow(x, 1.0 - params.vgM())/params.vgAlpha();
        }

        // value and derivative at regularization point
        Scalar Se_regu;
        if (Se<=PC_VG_REG)
            Se_regu = PC_VG_REG;
        else
            Se_regu = 1-PC_VG_REG;
        const Evaluation& x = std::pow(Se_regu,-1/params.vgM())-1;
        const Evaluation& pc = Opm::pow(x, 1.0/params.vgN())/params.vgAlpha();
        const Evaluation& pc_prime =
            Opm::pow(x, 1/params.vgN()-1)
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
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params& params, const FluidState& fluidState)
    {
        const Evaluation& Sw =
            Opm::decay<Evaluation>(fluidState.saturation(wettingPhaseIdx));
        Evaluation Se = (Sw-params.Swr())/(1.-params.Snr());

        Scalar PC_VG_REG = 0.01;

        // regularization
        if (Se<0.0)
            Se=0.0;
        if (Se>1.0)
            Se=1.0;

        if (Se>PC_VG_REG && Se<1-PC_VG_REG) {
            Evaluation x = Opm::pow(Se,-1/params.vgM()) - 1.0;
            x = Opm::pow(x, 1 - params.vgM());
            return x/params.vgAlpha();
        }

        // value and derivative at regularization point
        Scalar Se_regu;
        if (Se<=PC_VG_REG)
            Se_regu = PC_VG_REG;
        else
            Se_regu = 1.0 - PC_VG_REG;

        const Evaluation& x = std::pow(Se_regu,-1/params.vgM())-1;
        const Evaluation& pc = Opm::pow(x, 1/params.vgN())/params.vgAlpha();
        const Evaluation& pc_prime =
            Opm::pow(x,1/params.vgN()-1)
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
    static void saturations(ContainerT& /*values*/,
                            const Params& /*params*/,
                            const FluidState& /*fluidState*/)
    { throw std::logic_error("Not implemented: inverse capillary pressures"); }

    /*!
     * \brief The saturation of the gas phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sg(const Params& /*params*/, const FluidState& /*fluidState*/)
    { throw std::logic_error("Not implemented: Sg()"); }

    /*!
     * \brief The saturation of the non-wetting (i.e., oil) phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params& /*params*/, const FluidState& /*fluidState*/)
    { throw std::logic_error("Not implemented: Sn()"); }

    /*!
     * \brief The saturation of the wetting (i.e., water) phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params& /*params*/, const FluidState& /*fluidState*/)
    { throw std::logic_error("Not implemented: Sw()"); }

    /*!
     * \brief The relative permeability of all phases.
     */
    template <class ContainerT, class FluidState>
    static void relativePermeabilities(ContainerT& values,
                                       const Params& params,
                                       const FluidState& fluidState)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[wettingPhaseIdx] = krw<FluidState, Evaluation>(params, fluidState);
        values[nonWettingPhaseIdx] = krn<FluidState, Evaluation>(params, fluidState);
        values[gasPhaseIdx] = krg<FluidState, Evaluation>(params, fluidState);
    }

    /*!
     * \brief The relative permeability for the wetting phase of the
     *        medium implied by van Genuchten's parameterization.
     *
     * The permeability of water in a 3p system equals the standard 2p description.
     * (see p61. in "Comparison of the Three-Phase Oil Relative Permeability Models"
     * MOJDEH  DELSHAD and GARY A. POPE, Transport in Porous Media 4 (1989), 59-83.)
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params& params, const FluidState& fluidState)
    {
        const Evaluation& Sw =
            Opm::decay<Evaluation>(fluidState.saturation(wettingPhaseIdx));
        // transformation to effective saturation
        const Evaluation& Se = (Sw - params.Swr()) / (1-params.Swr());

        // regularization
        if(Se > 1.0) return 1.;
        if(Se < 0.0) return 0.;

        const Evaluation& r = 1. - Opm::pow(1 - Opm::pow(Se, 1/params.vgM()), params.vgM());
        return Opm::sqrt(Se)*r*r;
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
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params& params, const FluidState& fluidState)
    {
        const Evaluation& Sn =
            Opm::decay<Evaluation>(fluidState.saturation(nonWettingPhaseIdx));
        const Evaluation& Sw =
            Opm::decay<Evaluation>(fluidState.saturation(wettingPhaseIdx));
        Evaluation Swe = Opm::min((Sw - params.Swr()) / (1 - params.Swr()), 1.);
        Evaluation Ste = Opm::min((Sw + Sn - params.Swr()) / (1 - params.Swr()), 1.);

        // regularization
        if(Swe <= 0.0) Swe = 0.;
        if(Ste <= 0.0) Ste = 0.;
        if(Ste - Swe <= 0.0) return 0.;

        Evaluation krn_;
        krn_ = Opm::pow(1 - Opm::pow(Swe, 1/params.vgM()), params.vgM());
        krn_ -= Opm::pow(1 - Opm::pow(Ste, 1/params.vgM()), params.vgM());
        krn_ *= krn_;

        if (params.krRegardsSnr())
        {
            // regard Snr in the permeability of the non-wetting
            // phase, see Helmig1997
            const Evaluation& resIncluded =
                Opm::max(Opm::min(Sw - params.Snr() / (1-params.Swr()), 1.0), 0.0);
            krn_ *= Opm::sqrt(resIncluded );
        }
        else
            krn_ *= Opm::sqrt(Sn / (1 - params.Swr()));

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
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krg(const Params& params, const FluidState& fluidState)
    {
        const Evaluation& Sg =
            Opm::decay<Evaluation>(fluidState.saturation(gasPhaseIdx));
        const Evaluation& Se = Opm::min(((1-Sg) - params.Sgr()) / (1 - params.Sgr()), 1.);

        // regularization
        if(Se > 1.0)
            return 0.0;
        if(Se < 0.0)
            return 1.0;

        Evaluation scaleFactor = 1.;
        if (Sg<=0.1) {
            scaleFactor = (Sg - params.Sgr())/(0.1 - params.Sgr());
            if (scaleFactor < 0.)
                return 0.0;
        }

        return scaleFactor
            * Opm::pow(1 - Se, 1.0/3.)
            * Opm::pow(1 - Opm::pow(Se, 1/params.vgM()), 2*params.vgM());
    }
};
} // namespace Opm

#endif
