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
 * \copydoc Opm::ParkerLenhard
 */
#ifndef OPM_PARKER_LENHARD_HPP
#define OPM_PARKER_LENHARD_HPP

#include "ParkerLenhardParams.hpp"

#include <opm/material/fluidmatrixinteractions/VanGenuchten.hpp>

#include <algorithm>
#include <cassert>

namespace Opm {

/*!
 * \brief Represents a scanning curve in the Parker-Lenhard hysteresis model.
 *
 * The class has pointers to the scanning curves
 * with higher and lower loop number, this saving
 * the history of the imbibitions and drainages.
 */
template <class ScalarT>
class PLScanningCurve
{
public:
    typedef ScalarT Scalar;

    /*!
     * \brief Constructs main imbibition curve.
     *
     * Further scanning curves can be added with
     * setNext.
     */
    PLScanningCurve(Scalar Swr)
    {
        loopNum_ = 0;
        prev_ = new PLScanningCurve(NULL, // prev
                                    this, // next
                                    -1, // loop number
                                    Swr, // Sw
                                    1e12, // pcnw
                                    Swr, // SwMic
                                    Swr); // SwMdc
        next_ = NULL;

        Sw_ = 1.0;
        pcnw_ = 0.0;
        SwMic_ = 1.0;
        SwMdc_ = 1.0;
    }

protected:
    PLScanningCurve(PLScanningCurve* prevSC,
                    PLScanningCurve* nextSC,
                    int loopN,
                    Scalar SwReversal,
                    Scalar pcnwReversal,
                    Scalar SwMiCurve,
                    Scalar SwMdCurve)
    {
        prev_ = prevSC;
        next_ = nextSC;
        loopNum_ = loopN;
        Sw_ = SwReversal;
        pcnw_ = pcnwReversal;
        SwMic_ = SwMiCurve;
        SwMdc_ = SwMdCurve;
    }

public:
    /*!
     * \brief Destructor. After it was called
     *        all references to the next() curve are
     *        invalid!
     */
    ~PLScanningCurve()
    {
        if (loopNum_ == 0)
            delete prev_;
        if (loopNum_ >= 0)
            delete next_;
    }

    /*!
     * \brief Return the previous scanning curve, i.e. the curve
     *        with one less reversal than the current one.
     */
    PLScanningCurve* prev() const
    { return prev_; }

    /*!
     * \brief Return the next scanning curve, i.e. the curve
     *        with one more reversal than the current one.
     */
    PLScanningCurve* next() const
    { return next_; }

    /*!
     * \brief Set the next scanning curve.
     *
     * Next in the sense of the number of reversals
     * from imbibition to drainage or vince versa. If this
     * curve already has a list of next curves, it is
     * deleted and thus forgotten.
     */
    void setNext(Scalar SwReversal,
                 Scalar pcnwReversal,
                 Scalar SwMiCurve,
                 Scalar SwMdCurve)
    {
        // if next_ is NULL, delete does nothing, so
        // this is valid!!
        delete next_;

        next_ = new PLScanningCurve(this, // prev
                                    NULL, // next
                                    loopNum() + 1,
                                    SwReversal,
                                    pcnwReversal,
                                    SwMiCurve,
                                    SwMdCurve);
    }

    /*!
     * \brief Returns true iff the given effective saturation
     *        Swei is within the scope of the curve, i.e.
     *        whether Swei is part of the curve's
     *        domain and the curve thus applies to Swi.
     */
    bool isValidAt_Sw(Scalar SwReversal)
    {
        if (isImbib())
            // for inbibition the given saturation
            // must be between the start of the
            // current imbibition and the the start
            // of the last drainage
            return this->Sw() < SwReversal && SwReversal < prev_->Sw();
        else
            // for drainage the given saturation
            // must be between the start of the
            // last imbibition and the start
            // of the current drainage
            return prev_->Sw() < SwReversal && SwReversal < this->Sw();
    }

    /*!
     * \brief Returns true iff the scanning curve is a
     *        imbibition curve.
     */
    bool isImbib()
    { return loopNum()%2 == 1; }

    /*!
     * \brief Returns true iff the scanning curve is a
     *        drainage curve.
     */
    bool isDrain()
    { return !isImbib(); }

    /*!
     * \brief The loop number of the scanning curve.
     *
     * The MDC is 0, PISC is 1, PDSC is 2, ...
     */
    int loopNum()
    { return loopNum_; }

    /*!
     * \brief Absolute wetting-phase saturation at the
     *        scanning curve's reversal point.
     */
    Scalar Sw() const
    { return Sw_; }

    /*!
     * \brief Capillary pressure at the last reversal point.
     */
    Scalar pcnw() const
    { return pcnw_; }

    /*!
     * \brief Apparent saturation of the last reversal point on
     *        the pressure MIC.
     */
    Scalar SwMic()
    { return SwMic_; }

    /*!
     * \brief Apparent saturation of the last reversal point on
     *        the pressure MDC.
     */
    Scalar SwMdc()
    { return SwMdc_; }

private:
    PLScanningCurve* prev_;
    PLScanningCurve* next_;

    int loopNum_;

    Scalar Sw_;
    Scalar pcnw_;

    Scalar SwMdc_;
    Scalar SwMic_;
};

/*!
 * \ingroup material
 * \brief Implements the Parker-Lenhard twophase
 *        p_c-Sw hysteresis model. This class adheres to the twophase
 *        capillary pressure API.
 */
template <class TraitsT, class ParamsT = ParkerLenhardParams<TraitsT> >
class ParkerLenhard : public TraitsT
{
public:
    typedef TraitsT Traits;
    typedef ParamsT Params;
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases
    static const int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The Parker-Lenhard capillary pressure law only "
                  "applies to the case of two fluid phases");

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = true;

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = true;

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

    static_assert(Traits::numPhases == 2,
                  "The number of fluid phases must be two if you want to use "
                  "this material law!");

private:
    typedef typename ParamsT::VanGenuchten VanGenuchten;
    typedef Opm::PLScanningCurve<Scalar> ScanningCurve;

public:
    /*!
     * \brief Resets the hysteresis model to the
     *        initial parameters on the main drainage curve
     */
    static void reset(Params& params)
    {
        delete params.mdc(); // this will work even if mdc_ == NULL!
        params.setMdc(new ScanningCurve(params.SwrPc()));
        params.setCsc(params.mdc());
        params.setPisc(NULL);
        params.setCurrentSnr(0.0);
    }

    /*!
     * \brief Set the current absolute saturation for the
     *        current timestep
     */
    template <class FluidState>
    static void update(Params& params, const FluidState& fs)
    {
        Scalar Sw = Opm::scalarValue(fs.saturation(Traits::wettingPhaseIdx));

        if (Sw > 1 - 1e-5) {
            // if the absolute saturation is almost 1,
            // it means that we're back to the beginning
            reset(params);
            return;
        }

        // find the loop number which corrosponds to the
        // given effective saturation
        ScanningCurve* curve = findScanningCurve_(params, Sw);

        // calculate the apparent saturation on the MIC and MDC
        // which yield the same capillary pressure as the
        // Sw at the current scanning curve
        Scalar pc = pcnw<FluidState, Scalar>(params, fs);
        Scalar Sw_mic = VanGenuchten::twoPhaseSatSw(params.micParams(), pc);
        Scalar Sw_mdc = VanGenuchten::twoPhaseSatSw(params.mdcParams(), pc);

        curve->setNext(Sw, pc, Sw_mic, Sw_mdc);
        if (!curve->next())
            return;

        params.setCsc(curve);

        // if we're back on the MDC, we also have a new PISC!
        if (params.csc() == params.mdc()) {
            params.setPisc(params.mdc()->next());
            params.setCurrentSnr(computeCurrentSnr_(params, Sw));
        }
    }

    /*!
     * \brief Returns the capillary pressure dependening on
     *        the phase saturations.
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = 0.0; // reference phase
        values[Traits::nonWettingPhaseIdx] = pcnw<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief Returns the capillary pressure dependening on
     *        the phase saturations.
     */
    template <class Container, class FluidState>
    static void saturations(Container& /*values*/, const Params& /*params*/, const FluidState& /*fs*/)
    { throw std::logic_error("Not implemented: ParkerLenhard::saturations()"); }

    /*!
     * \brief Returns the relative permeabilities of the phases
     *        dependening on the phase saturations.
     */
    template <class Container, class FluidState>
    static void relativePermeabilities(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = krw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = krn<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief Returns the capillary pressure dependend on
     *        the phase saturations.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params& params, const FluidState& fs)
    {
        const Evaluation& Sw =
            Opm::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));

        return twoPhaseSatPcnw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params& params, const Evaluation& Sw)
    {
        // calculate the current apparent saturation
        ScanningCurve* sc = findScanningCurve_(params, Opm::scalarValue(Sw));

        // calculate the apparant saturation
        const Evaluation& Sw_app = absoluteToApparentSw_(params, Sw);

        // if the apparent saturation exceeds the 'legal' limits,
        // we also the underlying material law decide what to do.
        if (Sw_app > 1) {
            return 0.0; // VanGenuchten::pcnw(params.mdcParams(), Sw_app);
        }

        // put the apparent saturation into the main imbibition or
        // drainage curve
        Scalar SwAppCurSC = absoluteToApparentSw_(params, sc->Sw());
        Scalar SwAppPrevSC = absoluteToApparentSw_(params, sc->prev()->Sw());
        const Evaluation& pos = (Sw_app - SwAppCurSC)/(SwAppPrevSC - SwAppCurSC);
        if (sc->isImbib()) {
            const Evaluation& SwMic =
                pos * (sc->prev()->SwMic() - sc->SwMic()) + sc->SwMic();

            return VanGenuchten::twoPhaseSatPcnw(params.micParams(), SwMic);
        }
        else { // sc->isDrain()
            const Evaluation& SwMdc =
                pos*(sc->prev()->SwMdc() - sc->SwMdc()) + sc->SwMdc();

            return VanGenuchten::twoPhaseSatPcnw(params.mdcParams(), SwMdc);
        }
    }

    /*!
     * \brief Calculate the wetting phase saturations depending on
     *        the phase pressures.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sw(const Params& /*params*/, const FluidState& /*fs*/)
    { throw std::logic_error("Not implemented: ParkerLenhard::Sw()"); }

    template <class Evaluation>
    static Evaluation twoPhaseSatSw(const Params& /*params*/, const Evaluation& /*pc*/)
    { throw std::logic_error("Not implemented: ParkerLenhard::twoPhaseSatSw()"); }

    /*!
     * \brief Calculate the non-wetting phase saturations depending on
     *        the phase pressures.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation Sn(const Params& params, const FluidState& fs)
    { return 1 - Sw<FluidState, Evaluation>(params, fs); }

    template <class Evaluation>
    static Evaluation twoPhaseSatSn(const Params& /*params*/, const Evaluation& /*pc*/)
    { throw std::logic_error("Not implemented: ParkerLenhard::twoPhaseSatSn()"); }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params& params, const FluidState& fs)
    {
        const Evaluation& Sw =
            Opm::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));

        return twoPhaseSatKrw(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params& params, const Evaluation& Sw)
    {
        // for the effective permeability we only use Land's law and
        // the relative permeability of the main drainage curve.
        const Evaluation& Sw_app = absoluteToApparentSw_(params, Sw);
        return VanGenuchten::twoPhaseSatKrw(params.mdcParams(), Sw_app);
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the params.
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params& params, const FluidState& fs)
    {
        const Evaluation& Sw =
            Opm::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));

        return twoPhaseSatKrn(params, Sw);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrn(const Params& params, const Evaluation& Sw)
    {
        // for the effective permeability we only use Land's law and
        // the relative permeability of the main drainage curve.
        const Evaluation& Sw_app = absoluteToApparentSw_(params, Sw);
        return VanGenuchten::twoPhaseSatKrn(params.mdcParams(), Sw_app);
    }

    /*!
     * \brief Convert an absolute wetting saturation to an apparent one.
     */
    template <class Evaluation>
    static Evaluation absoluteToApparentSw_(const Params& params, const Evaluation& Sw)
    {
        return effectiveToApparentSw_(params, absoluteToEffectiveSw_(params, Sw));
    }

private:
    /*!
     * \brief Convert an absolute wetting saturation to an effective one.
     *
     * \param Sw        Absolute saturation of the wetting phase \f${S}_w\f$.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Effective saturation of the wetting phase.
     */
    template <class Evaluation>
    static Evaluation absoluteToEffectiveSw_(const Params& params, const Evaluation& Sw)
    { return (Sw - params.SwrPc())/(1 - params.SwrPc()); }

    /*!
     * \brief Convert an effective wetting saturation to an absolute one.
     *
     * \param Swe       Effective saturation of the non-wetting phase \f$\overline{S}_n\f$.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Absolute saturation of the non-wetting phase.
     */
    template <class Evaluation>
    static Evaluation effectiveToAbsoluteSw_(const Params& params, const Evaluation& Swe)
    { return Swe*(1 - params.SwrPc()) + params.SwrPc(); }

    // return the effctive residual non-wetting saturation, given an
    // effective wetting saturation
    template <class Evaluation>
    static Evaluation computeCurrentSnr_(const Params& params, const Evaluation& Sw)
    {
        // regularize
        if (Sw > 1 - params.Snr())
            return 0.0;
        if (Sw < params.SwrPc())
            return params.Snr();

        if (params.Snr() == 0.0)
            return 0.0;

        // use Land's law
        Scalar R = 1.0/params.Snr() - 1;
        const Evaluation& curSnr = (1 - Sw)/(1 + R*(1 - Sw));

        // the current effective residual non-wetting saturation must
        // be smaller than the residual non-wetting saturation
        assert(curSnr <= params.Snr());

        return curSnr;
    }

    // returns the trapped effective non-wetting saturation for a
    // given wetting phase saturation
    template <class Evaluation>
    static Evaluation trappedEffectiveSn_(const Params& params, const Evaluation& Sw)
    {
        const Evaluation& Swe = absoluteToEffectiveSw_(params, Sw);
        Scalar SwePisc = absoluteToEffectiveSw_(params, params.pisc()->Sw());

        Scalar Snre = absoluteToEffectiveSw_(params, params.currentSnr());
        return Snre*(Swe - SwePisc) / (1 - Snre - SwePisc);
    }

    // returns the apparent saturation of the wetting phase depending
    // on the effective saturation
    template <class Evaluation>
    static Evaluation effectiveToApparentSw_(const Params& params, const Evaluation& Swe)
    {
        if (params.pisc() == NULL ||
            Swe <= absoluteToEffectiveSw_(params, params.pisc()->Sw()))
        {
            // we are on the main drainage curve, i.e.  no non-wetting
            // fluid is trapped -> apparent saturation == effective
            // saturation
            return Swe;
        }

        // we are on a imbibition or drainage curve which is not the
        // main drainage curve -> apparent saturation == effective
        // saturation + trapped effective saturation
        return Swe + trappedEffectiveSn_(params, Swe);
    }

    // Returns the effective saturation to a given apparent one
    template <class Evaluation>
    static Evaluation apparentToEffectiveSw_(const Params& params, const Evaluation& Swapp)
    {
        Scalar SwePisc = absoluteToEffectiveSw_(params, params.pisc()->Sw());
        if (params.pisc() == NULL || Swapp <= SwePisc) {
            // we are on the main drainage curve, i.e.
            // no non-wetting fluid is trapped
            // -> apparent saturation == effective saturation
            return Swapp;
        }

        Scalar Snre = absoluteToEffectiveSw_(params.currentSnr());
        return
            (Swapp*(1 - Snre - SwePisc) + Snre*SwePisc)
            /(1 - SwePisc);
    }


    // find the loop on which the an effective
    // saturation has to be
    static ScanningCurve* findScanningCurve_(const Params& params, Scalar Sw)
    {
        if (params.pisc() == NULL || Sw <= params.pisc()->Sw()) {
            // we don't have a PISC yet, or the effective
            // saturation is smaller than the saturation where the
            // PISC begins. In this case are on the MDC
            return params.mdc();
        }

        // if we have a primary imbibition curve, and our current
        // effective saturation is higher than the beginning of
        // the secondary drainage curve. this means we are on the
        // PISC again.
        if (params.pisc()->next() == NULL ||
            params.pisc()->next()->Sw() < Sw)
        {
            return params.pisc();
        }

        ScanningCurve* curve = params.csc()->next();
        while (true) {
            assert(curve != params.mdc()->prev());
            if (curve->isValidAt_Sw(Sw)) {
                return curve;
            }
            curve = curve->prev();
        }
    }
};

} // namespace Opm

#endif // PARKER_LENHARD_HPP
