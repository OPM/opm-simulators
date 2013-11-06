// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2011 by Philipp Nuske                                *
 *   Copyright (C) 2009-2012 by Andreas Lauser                               *
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
 * \file
 * \copydoc Opm::EffToAbsLaw
 */
#ifndef OPM_EFF_TO_ABS_LAW_HH
#define OPM_EFF_TO_ABS_LAW_HH

#include "EffToAbsLawParams.hpp"

#include <opm/material/fluidstates/SaturationOverlayFluidState.hpp>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief This material law takes a material law defined for effective
 *        saturations and converts it to a material law defined on absolute
 *        saturations.
 *
 * The idea: "material laws" (like VanGenuchten or BrooksCorey) are
 * defined for effective saturations.  The numeric calculations
 * however are performed with absolute saturations. The EffToAbsLaw
 * class gets the "material laws" actually used as well as the
 * corresponding parameter container as template arguments.
 *
 * Subsequently, the desired function (pc, Sw... ) of the actually
 * used "material laws" are called but with the saturations already
 * converted from absolute to effective.
 *
 * This approach makes sure that in the "material laws" only effective
 * saturations are considered, which makes sense, as these laws only
 * deal with effective saturations. This also allows for changing the
 * calculation of the effective saturations easily, as this is subject
 * of discussion / may be problem specific.
 *
 * Additionally, handing over effective saturations to the "material
 * laws" in stead of them calculating effective saturations prevents
 * accidently "converting twice".
 *
 * This boils down to:
 * - the actual material laws (linear, VanGenuchten...) do not need to
 *   deal with any kind of conversion
 * - the definition of the material law in the spatial parameters is
 *   not really intuitive, but using it is: Hand in values, get back
 *   values, do not deal with conversion.
 */
template <class EffLawT, class ParamsT = EffToAbsLawParams<typename EffLawT::Params, EffLawT::numPhases> >
class EffToAbsLaw : public EffLawT::Traits
{
    typedef EffLawT EffLaw;

public:
    typedef typename EffLaw::Traits Traits;
    typedef ParamsT Params;
    typedef typename EffLaw::Scalar Scalar;

    //! The number of fluid phases
    static const int numPhases = EffLaw::numPhases;

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = EffLaw::implementsTwoPhaseApi;

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = EffLaw::implementsTwoPhaseSatApi;

    //! Specify whether the quantities defined by this material law
    //! are saturation dependent
    static const bool isSaturationDependent = EffLaw::isSaturationDependent;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the absolute pressure
    static const bool isPressureDependent = EffLaw::isPressureDependent;

    //! Specify whether the quantities defined by this material law
    //! are temperature dependent
    static const bool isTemperatureDependent = EffLaw::isTemperatureDependent;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the phase composition
    static const bool isCompositionDependent = EffLaw::isCompositionDependent;

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
        typedef Opm::SaturationOverlayFluidState<Scalar, FluidState> OverlayFluidState;

        OverlayFluidState overlayFs(fs);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        EffLaw::capillaryPressures(values, params, overlayFs);
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
        typedef Opm::SaturationOverlayFluidState<Scalar, FluidState> OverlayFluidState;

        OverlayFluidState overlayFs(fs);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        EffLaw::relativePermeabilities(values, params, overlayFs);
    }

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     *
     * \param Sw Absolute saturation of the wetting phase
     *           \f$\overline{S}_w\f$. It is converted to effective
     *           saturation and then handed over to the material law
     *           actually used for calculation.
     * \param params A object that stores the appropriate coefficients
     *                for the respective law.
     *
     * \return Capillary pressure [Pa] calculated by specific
     *         constitutive relation (e.g. Brooks & Corey, van
     *         Genuchten, linear...)
     */
    template <class FluidState>
    static Scalar pcwn(const Params &params, const FluidState &fs)
    {
        typedef Opm::SaturationOverlayFluidState<Scalar, FluidState> OverlayFluidState;

        static_assert(FluidState::numPhases == numPhases,
                      "The fluid state and the material law must exhibit the same "
                      "number of phases!");

        OverlayFluidState overlayFs(fs);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        return EffLaw::pcwn(params, overlayFs);
    }

    template <class ScalarT = Scalar>
    static typename std::enable_if<implementsTwoPhaseSatApi, ScalarT>::type
    twoPhaseSatPcwn(const Params &params, Scalar SwAbs)
    {
        Scalar SwEff = effectiveSaturation(params, SwAbs, Traits::wPhaseIdx);

        return EffLaw::twoPhaseSatPcwn(params, SwEff);
    }

    /*!
     * \brief The saturation-capillary pressure curves.
     */
    template <class Container, class FluidState>
    static void saturations(Container &values, const Params &params, const FluidState &fs)
    {
        EffLaw::saturations(values, params, fs);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            values[phaseIdx] = absoluteSaturation(params, values[phaseIdx], phaseIdx);
        }
    }

    /*!
     * \brief Calculate wetting liquid phase saturation given that
     *        the rest of the fluid state has been initialized
     */
    template <class FluidState>
    static Scalar Sw(const Params &params, const FluidState &fs)
    { return absoluteSaturation(params, EffLaw::Sw(params, fs), Traits::wPhaseIdx); }

    template <class ScalarT = Scalar>
    static typename std::enable_if<implementsTwoPhaseSatApi, ScalarT>::type
    twoPhaseSatSw(const Params &params, Scalar Sw)
    { return absoluteSaturation(params, EffLaw::twoPhaseSatSw(params, Sw), Traits::wPhaseIdx); }

    /*!
     * \brief Calculate non-wetting liquid phase saturation given that
     *        the rest of the fluid state has been initialized
     */
    template <class FluidState>
    static Scalar Sn(const Params &params, const FluidState &fs)
    { return absoluteSaturation(params, EffLaw::Sn(params, fs), Traits::nPhaseIdx); }

    template <class ScalarT = Scalar>
    static typename std::enable_if<implementsTwoPhaseSatApi, ScalarT>::type
    twoPhaseSatSn(const Params &params, Scalar Sw)
    { return absoluteSaturation(params, EffLaw::twoPhaseSatSn(params, Sw), Traits::nPhaseIdx); }

    /*!
     * \brief Calculate gas phase saturation given that the rest of
     *        the fluid state has been initialized
     *
     * This method is only available for at least three fluid phases
     */
    template <class FluidState, class ScalarT = Scalar>
    static typename std::enable_if< (Traits::numPhases > 2), ScalarT>::type
    Sg(const Params &params, const FluidState &fs)
    { return absoluteSaturation(params, EffLaw::Sg(params, fs), Traits::gPhaseIdx); }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure w.r.t the absolute saturation.
     *
     * In this case the chain rule needs to be applied:
     \f[
             p_c = p_c( \overline S_w (S_w))
             \rightarrow p_c ^\prime = \frac{\partial  p_c}{\partial \overline S_w} \frac{\partial \overline S_w}{\partial S_w}
     \f]
     * \param Sw        Absolute saturation of the wetting phase \f$\overline{S}_w\f$.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of \f$p_c\f$ w.r.t. effective saturation according to EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
    */
    static Scalar dpC_dSw(const Params &params, Scalar Sw)
    {
        return EffLaw::dpC_dSw(params, SabsToSeff(params, Sw) )*dSwe_dSw_(params);
    }

    /*!
     * \brief Returns the partial derivative of the absolute
     *        saturation w.r.t. the capillary pressure.
     *
     * In this case the chain rule needs to be applied:
     \f[
            S_w = S_w(\overline{S}_w (p_c) )
            \rightarrow S_w^\prime = \frac{\partial S_w}{\partial \overline S_w} \frac{\partial \overline S_w}{\partial p_c}
     \f]
     *
     *
     * \param pC        Capillary pressure \f$p_C\f$:
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Partial derivative of effective saturation w.r.t. \f$p_c\f$ according to EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        return EffLaw::dSw_dpC(params, pC)*dSw_dSwe_(params);
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param Sw        Absolute saturation of the wetting phase \f$\overline{S}_w\f$. It is converted to effective saturation
     *                  and then handed over to the material law actually used for calculation.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the wetting phase calculated as implied by EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     *
     */
    template <class FluidState>
    static Scalar krw(const Params &params, const FluidState &fs)
    {
        typedef Opm::SaturationOverlayFluidState<Scalar, FluidState> OverlayFluidState;

        static_assert(FluidState::numPhases == numPhases,
                      "The fluid state and the material law must exhibit the same "
                      "number of phases!");

        OverlayFluidState overlayFs(fs);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        return EffLaw::krw(params, overlayFs);
    }

    template <class ScalarT = Scalar>
    static typename std::enable_if<implementsTwoPhaseSatApi, ScalarT>::type
    twoPhaseSatKrw(const Params &params, Scalar Sw)
    { return EffLaw::twoPhaseSatKrw(params, effectiveSaturation(params, Sw, Traits::nPhaseIdx)); }

    /*!
     * \brief The relative permeability of the non-wetting phase.
     */
    template <class FluidState>
    static Scalar krn(const Params &params, const FluidState &fs)
    {
        typedef Opm::SaturationOverlayFluidState<Scalar, FluidState> OverlayFluidState;

        static_assert(FluidState::numPhases == numPhases,
                      "The fluid state and the material law must exhibit the same "
                      "number of phases!");

        OverlayFluidState overlayFs(fs);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        return EffLaw::krn(params, overlayFs);
    }

    template <class ScalarT = Scalar>
    static typename std::enable_if<implementsTwoPhaseSatApi, ScalarT>::type
    twoPhaseSatKrn(const Params &params, Scalar Sw)
    { return EffLaw::twoPhaseSatKrn(params, effectiveSaturation(params, Sw, Traits::nPhaseIdx)); }

    /*!
     * \brief The relative permability of the gas phase
     *
     * This method is only available for at least three fluid phases
     */
    template <class FluidState, class ScalarT=Scalar>
    static typename std::enable_if< (Traits::numPhases > 2), ScalarT>::type
    krg(const Params &params, const FluidState &fs)
    {
        typedef Opm::SaturationOverlayFluidState<Scalar, FluidState> OverlayFluidState;

        static_assert(FluidState::numPhases == numPhases,
                      "The fluid state and the material law must exhibit the same "
                      "number of phases!");

        OverlayFluidState overlayFs(fs);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        return EffLaw::krg(params, overlayFs);
    }

    /*!
     * \brief Convert an absolute saturation to an effective one.
     */
    static Scalar effectiveSaturation(const Params &params, Scalar S, int phaseIdx)
    { return (S - params.residualSaturation(phaseIdx))/(1 - params.sumResidualSaturations()); }

    /*!
     * \brief Convert an effective saturation to an absolute one.
     */
    static Scalar absoluteSaturation(const Params &params, Scalar S, int phaseIdx)
    { return S*(1 - params.sumResidualSaturations()) + params.residualSaturation(phaseIdx); }

private:
    /*!
     * \brief Convert an effective wetting saturation to an absolute one.
     *
     * \param Swe       Effective saturation of the non-wetting phase \f$\overline{S}_n\f$.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Absolute saturation of the non-wetting phase.
     */
    static Scalar SweToSw_(const Params &params, Scalar Swe)
    {
        return Swe*(1 - params.Swr() - params.Snr()) + params.Swr();
    }

    /*!
     * \brief           Derivative of the effective saturation w.r.t. the absolute saturation.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Derivative of the effective saturation w.r.t. the absolute saturation.
     */
    static Scalar dSwe_dSw_(const Params &params)
    { return 1.0/(1 - params.Swr() - params.Snr()); }

    /*!
     * \brief           Derivative of the absolute saturation w.r.t. the effective saturation.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Derivative of the absolute saturation w.r.t. the effective saturation.
     */
    static Scalar dSw_dSwe_(const Params &params)
    { return 1 - params.Swr() - params.Snr(); }
};
} // namespace Opm

#endif
