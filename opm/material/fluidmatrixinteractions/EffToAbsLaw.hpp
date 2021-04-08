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
 * \copydoc Opm::EffToAbsLaw
 */
#ifndef OPM_EFF_TO_ABS_LAW_HPP
#define OPM_EFF_TO_ABS_LAW_HPP

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
    static void capillaryPressures(Container& values, const Params& params, const FluidState& fs)
    {
        typedef Opm::SaturationOverlayFluidState<FluidState> OverlayFluidState;

        OverlayFluidState overlayFs(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        EffLaw::template capillaryPressures<Container, OverlayFluidState>(values, params, overlayFs);
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
    static void relativePermeabilities(Container& values, const Params& params, const FluidState& fs)
    {
        typedef Opm::SaturationOverlayFluidState<FluidState> OverlayFluidState;

        OverlayFluidState overlayFs(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        EffLaw::template relativePermeabilities<Container, OverlayFluidState>(values, params, overlayFs);
    }

    /*!
     * \brief The capillary pressure-saturation curve.
     *
     *
     * \param params A object that stores the appropriate coefficients
     *                for the respective law.
     *
     * \return Capillary pressure [Pa] calculated by specific
     *         constitutive relation (e.g. Brooks & Corey, van
     *         Genuchten, linear...)
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params& params, const FluidState& fs)
    {
        typedef Opm::SaturationOverlayFluidState<FluidState> OverlayFluidState;

        static_assert(FluidState::numPhases == numPhases,
                      "The fluid state and the material law must exhibit the same "
                      "number of phases!");

        OverlayFluidState overlayFs(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        return EffLaw::template pcnw<OverlayFluidState, Evaluation>(params, overlayFs);
    }

    template <class Evaluation>
    static typename std::enable_if<implementsTwoPhaseSatApi, Evaluation>::type
    twoPhaseSatPcnw(const Params& params, const Evaluation& SwAbs)
    {
        const Evaluation& SwEff = effectiveSaturation(params, SwAbs, Traits::wettingPhaseIdx);

        return EffLaw::twoPhaseSatPcnw(params, SwEff);
    }

    /*!
     * \brief The saturation-capillary pressure curves.
     */
    template <class Container, class FluidState>
    static void saturations(Container& values, const Params& params, const FluidState& fs)
    {
        EffLaw::template saturations<Container, FluidState>(values, params, fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            values[phaseIdx] = absoluteSaturation(params, values[phaseIdx], phaseIdx);
        }
    }

    /*!
     * \brief Calculate wetting liquid phase saturation given that
     *        the rest of the fluid state has been initialized
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params& params, const FluidState& fs)
    {
        return absoluteSaturation(params,
                                  EffLaw::template Sw<FluidState, Evaluation>(params, fs),
                                  Traits::wettingPhaseIdx);
    }

    template <class Evaluation>
    static typename std::enable_if<implementsTwoPhaseSatApi, Evaluation>::type
    twoPhaseSatSw(const Params& params, const Evaluation& Sw)
    { return absoluteSaturation(params,
                                EffLaw::twoPhaseSatSw(params, Sw),
                                Traits::wettingPhaseIdx); }

    /*!
     * \brief Calculate non-wetting liquid phase saturation given that
     *        the rest of the fluid state has been initialized
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params& params, const FluidState& fs)
    {
        return absoluteSaturation(params,
                                  EffLaw::template Sn<FluidState, Evaluation>(params, fs),
                                  Traits::nonWettingPhaseIdx);
    }

    template <class Evaluation>
    static typename std::enable_if<implementsTwoPhaseSatApi, Evaluation>::type
    twoPhaseSatSn(const Params& params, const Evaluation& Sw)
    {
        return absoluteSaturation(params,
                                  EffLaw::twoPhaseSatSn(params, Sw),
                                  Traits::nonWettingPhaseIdx);
    }

    /*!
     * \brief Calculate gas phase saturation given that the rest of
     *        the fluid state has been initialized
     *
     * This method is only available for at least three fluid phases
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static typename std::enable_if< (Traits::numPhases > 2), Evaluation>::type
    Sg(const Params& params, const FluidState& fs)
    {
        return absoluteSaturation(params,
                                  EffLaw::template Sg<FluidState, Evaluation>(params, fs),
                                  Traits::gasPhaseIdx);
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the wetting phase calculated as implied by EffLaw e.g. Brooks & Corey, van Genuchten, linear... .
     *
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params& params, const FluidState& fs)
    {
        typedef Opm::SaturationOverlayFluidState<FluidState> OverlayFluidState;

        static_assert(FluidState::numPhases == numPhases,
                      "The fluid state and the material law must exhibit the same "
                      "number of phases!");

        OverlayFluidState overlayFs(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        return EffLaw::template krw<OverlayFluidState, Evaluation>(params, overlayFs);
    }

    template <class Evaluation>
    static typename std::enable_if<implementsTwoPhaseSatApi, Evaluation>::type
    twoPhaseSatKrw(const Params& params, const Evaluation& Sw)
    { return EffLaw::twoPhaseSatKrw(params, effectiveSaturation(params, Sw, Traits::wettingPhaseIdx)); }

    /*!
     * \brief The relative permeability of the non-wetting phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params& params, const FluidState& fs)
    {
        typedef Opm::SaturationOverlayFluidState<FluidState> OverlayFluidState;

        static_assert(FluidState::numPhases == numPhases,
                      "The fluid state and the material law must exhibit the same "
                      "number of phases!");

        OverlayFluidState overlayFs(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        return EffLaw::template krn<OverlayFluidState, Evaluation>(params, overlayFs);
    }

    template <class Evaluation>
    static typename std::enable_if<implementsTwoPhaseSatApi, Evaluation>::type
    twoPhaseSatKrn(const Params& params, const Evaluation& Sw)
    { return EffLaw::twoPhaseSatKrn(params, effectiveSaturation(params, Sw, Traits::wettingPhaseIdx)); }

    /*!
     * \brief The relative permability of the gas phase
     *
     * This method is only available for at least three fluid phases
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static typename std::enable_if< (Traits::numPhases > 2), Evaluation>::type
    krg(const Params& params, const FluidState& fs)
    {
        typedef Opm::SaturationOverlayFluidState<FluidState> OverlayFluidState;

        static_assert(FluidState::numPhases == numPhases,
                      "The fluid state and the material law must exhibit the same "
                      "number of phases!");

        OverlayFluidState overlayFs(fs);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            overlayFs.setSaturation(phaseIdx,
                                    effectiveSaturation(params,
                                                        fs.saturation(phaseIdx),
                                                        phaseIdx));
        }

        return EffLaw::template krg<OverlayFluidState, Evaluation>(params, overlayFs);
    }

    /*!
     * \brief Convert an absolute saturation to an effective one.
     */
    template <class Evaluation>
    static Evaluation effectiveSaturation(const Params& params, const Evaluation& S, unsigned phaseIdx)
    { return (S - params.residualSaturation(phaseIdx))/(1.0 - params.sumResidualSaturations()); }

    /*!
     * \brief Convert an effective saturation to an absolute one.
     */
    template <class Evaluation>
    static Evaluation absoluteSaturation(const Params& params, const Evaluation& S, unsigned phaseIdx)
    { return S*(1.0 - params.sumResidualSaturations()) + params.residualSaturation(phaseIdx); }

private:
    /*!
     * \brief           Derivative of the effective saturation w.r.t. the absolute saturation.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Derivative of the effective saturation w.r.t. the absolute saturation.
     */
    static Scalar dSeff_dSabs_(const Params& params, int /*phaseIdx*/)
    { return 1.0/(1 - params.sumResidualSaturations()); }

    /*!
     * \brief           Derivative of the absolute saturation w.r.t. the effective saturation.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Derivative of the absolute saturation w.r.t. the effective saturation.
     */
    static Scalar dSabs_dSeff_(const Params& params, int /*phaseIdx*/)
    { return 1 - params.sumResidualSaturations(); }
};
} // namespace Opm

#endif
