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
 * \copydoc Opm::EclHysteresisTwoPhaseLaw
 */
#ifndef OPM_ECL_HYSTERESIS_TWO_PHASE_LAW_HPP
#define OPM_ECL_HYSTERESIS_TWO_PHASE_LAW_HPP

#include "EclHysteresisTwoPhaseLawParams.hpp"

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief This material law implements the hysteresis model of the ECL file format
 */
template <class EffectiveLawT,
          class ParamsT = EclHysteresisTwoPhaseLawParams<EffectiveLawT> >
class EclHysteresisTwoPhaseLaw : public EffectiveLawT::Traits
{
public:
    typedef EffectiveLawT EffectiveLaw;
    typedef typename EffectiveLaw::Params EffectiveLawParams;

    typedef typename EffectiveLaw::Traits Traits;
    typedef ParamsT Params;
    typedef typename EffectiveLaw::Scalar Scalar;

    enum { wettingPhaseIdx = Traits::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = Traits::nonWettingPhaseIdx };

    //! The number of fluid phases
    static const int numPhases = EffectiveLaw::numPhases;
    static_assert(numPhases == 2,
                  "The endpoint scaling applies to the nested twophase laws, not to "
                  "the threephase one!");

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = true;

    static_assert(EffectiveLaw::implementsTwoPhaseApi,
                  "The material laws put into EclEpsTwoPhaseLaw must implement the "
                  "two-phase material law API!");

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = true;

    static_assert(EffectiveLaw::implementsTwoPhaseSatApi,
                  "The material laws put into EclEpsTwoPhaseLaw must implement the "
                  "two-phase material law saturation API!");

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
    static void capillaryPressures(Container& /* values */,
                                   const Params& /* params */,
                                   const FluidState& /* fs */)
    {
        throw std::invalid_argument("The capillaryPressures(fs) method is not yet implemented");
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
    static void relativePermeabilities(Container& /* values */,
                                       const Params& /* params */,
                                       const FluidState& /* fs */)
    {
        throw std::invalid_argument("The pcnw(fs) method is not yet implemented");
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
    static Evaluation pcnw(const Params& /* params */,
                           const FluidState& /* fs */)
    {
        throw std::invalid_argument("The pcnw(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params& params, const Evaluation& Sw)
    {
        // TODO: capillary pressure hysteresis
        return EffectiveLaw::twoPhaseSatPcnw(params.drainageParams(), Sw);
/*
        if (!params.config().enableHysteresis() || params.config().pcHysteresisModel() < 0)
            return EffectiveLaw::twoPhaseSatPcnw(params.drainageParams(), Sw);

        if (Sw < params.SwMdc())
            return EffectiveLaw::twoPhaseSatPcnw(params.drainageParams(), Sw);

        const Evaluation& SwEff = Sw;

        //return EffectiveLaw::twoPhaseSatPcnw(params.imbibitionParams(), SwEff);
        return EffectiveLaw::twoPhaseSatPcnw(params.drainageParams(), SwEff);
*/
    }

    /*!
     * \brief The saturation-capillary pressure curves.
     */
    template <class Container, class FluidState>
    static void saturations(Container& /* values */,
                            const Params& /* params */,
                            const FluidState& /* fs */)
    {
        throw std::invalid_argument("The saturations(fs) method is not yet implemented");
    }

    /*!
     * \brief Calculate wetting liquid phase saturation given that
     *        the rest of the fluid state has been initialized
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params& /* params */,
                         const FluidState& /* fs */)
    {
        throw std::invalid_argument("The Sw(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatSw(const Params& /* params */,
                                    const Evaluation& /* pc */)
    {
        throw std::invalid_argument("The twoPhaseSatSw(pc) method is not yet implemented");
    }

    /*!
     * \brief Calculate non-wetting liquid phase saturation given that
     *        the rest of the fluid state has been initialized
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params& /* params */,
                         const FluidState& /* fs */)
    {
        throw std::invalid_argument("The Sn(pc) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatSn(const Params& /* params */,
                                    const Evaluation& /* pc */)
    {
        throw std::invalid_argument("The twoPhaseSatSn(pc) method is not yet implemented");
    }

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Relative permeability of the wetting phase calculated as implied by EffectiveLaw e.g. Brooks & Corey, van Genuchten, linear... .
     *
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params& /* params */,
                          const FluidState& /* fs */)
    {
        throw std::invalid_argument("The krw(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params& params, const Evaluation& Sw)
    {

        // if no relperm hysteresis is enabled, use the drainage curve
        if (!params.config().enableHysteresis() || params.config().krHysteresisModel() < 0)
            return EffectiveLaw::twoPhaseSatKrw(params.drainageParams(), Sw);

        // if it is enabled, use either the drainage or the imbibition curve. if the
        // imbibition curve is used, the saturation must be shifted.
        if (Sw <= params.krwSwMdc())
            return EffectiveLaw::twoPhaseSatKrw(params.drainageParams(), Sw);

        return EffectiveLaw::twoPhaseSatKrw(params.imbibitionParams(),
                                            Sw + params.deltaSwImbKrw());
    }

    /*!
     * \brief The relative permeability of the non-wetting phase.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params& /* params */,
                          const FluidState& /* fs */)
    {
        throw std::invalid_argument("The krn(fs) method is not yet implemented");
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrn(const Params& params, const Evaluation& Sw)
    {
        // if no relperm hysteresis is enabled, use the drainage curve
        if (!params.config().enableHysteresis() || params.config().krHysteresisModel() < 0)
            return EffectiveLaw::twoPhaseSatKrn(params.drainageParams(), Sw);

        // if it is enabled, use either the drainage or the imbibition curve. if the
        // imbibition curve is used, the saturation must be shifted.
        if (Sw <= params.krnSwMdc())
            return EffectiveLaw::twoPhaseSatKrn(params.drainageParams(), Sw);

        return EffectiveLaw::twoPhaseSatKrn(params.imbibitionParams(),
                                            Sw + params.deltaSwImbKrn());
    }
};
} // namespace Opm

#endif
