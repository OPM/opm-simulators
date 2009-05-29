/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser                                    *
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
 * \file ParkerLenhard.hh
 * \brief Template implementing Parker-Lenhard capillary pressure hystersis.
 */
#ifndef PARKER_LENHARD_HH
#define PARKER_LENHARD_HH

#include <dumux/auxiliary/apis.hh>

#include <dumux/new_material/twophasesat.hh>
#include <dumux/new_material/vangenuchten.hh>
#include <dumux/new_material/parkerlenhardstate.hh>

#include <dumux/auxiliary/spline.hh>
#include <dumux/auxiliary/expspline.hh>

#include <math.h>
#include <assert.h>

#include <algorithm>

#include <iostream>
#include <boost/format.hpp>

//#define PL_STRICT_ASSERTS

#ifdef PL_STRICT_ASSERTS
// we use a macro here, since if the assertation fails at run time,
// we want to see at which line it happened
#define ASSERT_RANGE(VAL, LEFT, RIGHT) assert(LEFT <= VAL && VAL <= RIGHT);
#else // !PL_STRICT_ASSERTS
template<class Scalar, class Scalar1, class Scalar2>
void ASSERT_RANGE(Scalar &var, Scalar1 left, Scalar2 right)
{ var = std::max(Scalar(left), std::min(Scalar(right), var)); }
#endif

namespace Dune
{
/*!
 * \internal
 * \brief Represents a scanning curve.
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
    typedef Dune::ExpSpline<Scalar> Spline;
    //        typedef Dune::Spline<Scalar> Spline;

    /*!
     * \brief Constructs main imbibtion curve.
     *
     * Further scanning curves can be added with
     * setNext.
     */
    PLScanningCurve()
    {
        loopNum_ = 0;
        prev_ = new PLScanningCurve(NULL, // prev
                                    this, // next
                                    -1,   // loop number
                                    0,    // Sw
                                    1e12, // pC
                                    0,    // Sw_app
                                    0,    // SwMic
                                    0);   // SwMdc
        next_ = NULL;

        Swe_ = 1.0;
        pC_  = 0.0;
        Sw_app_ = 1.0;
        SwMic_ = 1.0;
        SwMdc_ = 1.0;

        spline_ = NULL;
    }

protected:
    PLScanningCurve(PLScanningCurve *prev,
                    PLScanningCurve *next,
                    int    loopN,
                    Scalar Swe,
                    Scalar pC,
                    Scalar Sw_app,
                    Scalar SwMic,
                    Scalar SwMdc)
    {
        prev_ = prev;
        next_ = next;
        loopNum_ = loopN;
        Swe_ = Swe;
        pC_  = pC;
        Sw_app_ = Sw_app;
        SwMic_ = SwMic;
        SwMdc_ = SwMdc;

        spline_ = NULL;
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
        delete spline_;
    }

    /*!
     * \brief Return the previous scanning curve, i.e. the curve
     *        with one less reversal than the current one.
     */
    PLScanningCurve *prev() const
    { return prev_; }

    /*!
     * \brief Return the next scanning curve, i.e. the curve
     *        with one more reversal than the current one.
     */
    PLScanningCurve *next() const
    { return next_; }

    /*!
     * \brief Set the next scanning curve.
     *
     * Next in the sense of the number of reversals
     * from imbibition to drainage or vince versa. If this
     * curve already has a list of next curves, it is
     * deleted and thus forgotten.
     */
    void setNext(Scalar Swe,
                 Scalar pC,
                 Scalar Sw_app,
                 Scalar SwMic,
                 Scalar SwMdc)
    {
        // if next_ is NULL, delete does nothing, so
        // this is valid!!
        delete next_;

        next_ = new PLScanningCurve(this, // prev
                                    NULL, // next
                                    loopNum() + 1,
                                    Swe,
                                    pC,
                                    Sw_app,
                                    SwMic,
                                    SwMdc);
    }

    /*!
     * \brief Set the spline used for regularization.
     */
    void setSpline(Spline *spline)
    {
        delete spline_;
        spline_ = spline;
    };

    /*!
     * \brief Returns true if the spline shall be used for an
     *        absolute saturation.
     */
    bool useSpline(Scalar Sw)
    { return spline_ && spline_->applies(Sw); }

    /*!
     * \brief Returns the spline to be used for regularization
     *        of the current reversal point.
     */
    Spline *spline()
    { return spline_; }

    /*!
     * \brief Returns true iff the given effective saturation
     *        Swei is within the scope of the curve, i.e.
     *        whether Swei is part of the curve's
     *        domain and the curve thus applies to Swi.
     */
    bool isValidAt_Swe(Scalar Swei)
    {
        if (isImbib())
            // for inbibition the given saturation
            // must be between the start of the
            // current imbibition and the the start
            // of the last drainage
            return Swe() < Swei && Swei < prev_->Swe();
        else
            // for drainage the given saturation
            // must be between the start of the
            // last imbibition and the start
            // of the current drainage
            return prev_->Swe() < Swei && Swei < Swe();
    };

    /*!
     * \brief Returns true iff a given capillary pressure
     *        pC is within the scope of the curve, i.e.
     *        whether pC is part of the curve's
     *        image and the curve thus applies to pC.
     */
    bool isValidAt_pC(Scalar pCi)
    {
        if (isImbib())
            return prev_->pC() < pCi  && pCi < pC();
        else
            return pC() < pCi && pCi < prev_->pC();
    };

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
     * \brief Effective saturation at the
     *        last reversal point.
     */
    Scalar Swe() const
    { return Swe_; }

    /*!
     * \brief Capillary pressure at the last reversal point.
     */
    Scalar pC() const
    { return pC_; }

    /*!
     * \brief Apparent saturation at the
     *        last reversal point.
     */
    Scalar Sw_app() const
    { return Sw_app_; }

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
    PLScanningCurve *prev_;
    PLScanningCurve *next_;

    int    loopNum_;

    Scalar Swe_;
    Scalar pC_;
    Scalar Sw_app_;

    Scalar SwMdc_;
    Scalar SwMic_;

    Spline *spline_;
};

/*!
 * \ingroup material
 * \brief Implements the Parker-Lenhard twophase
 *        p_c-Sw hysteresis model. This class adheres to the twophase
 *        capillary pressure API.
 */
template <class StateT, class CapPressureT, bool useSplines=true>
class ParkerLenhard
{
public:
    typedef CapPressureT CapPressure;
    typedef StateT State;

    typedef typename State::Scalar Scalar;
    typedef Dune::PLScanningCurve<Scalar> ScanningCurve;
    typedef typename ScanningCurve::Spline Spline;

    typedef Dune::TwophaseSat<State> TwophaseSat;

    /*!
     * \brief Resets the hysteresis model to the
     *        initial state on the main drainage curve
     */
    static void reset(State &state)
    {
        Api::require<Api::ParkerLenhardState>(state);
        Api::require<Api::TwophaseSatParams>(state);

        delete state.mdc(); // this will work even if mdc_ == NULL!
        state.setMdc(new ScanningCurve());
        state.setCsc(state.mdc());
        state.setPisc(NULL);
        state.setSnrei(0.0);

        //            addSpline_(state, state.mdc()->prev());
        //            addSpline_(state, state.mdc());
    };

    /*!
     * \brief Returns true iff the given absolute saturation
     *        is viable (i.e. larger than the residual
     *        saturation of the wetting phase but smaller
     *        than the maximum saturation of the current
     *        PISC)
     */
    static bool isValidSw(const State &state, Scalar Sw)
    {
        Api::require<Api::ParkerLenhardParams>(state);
        Api::require<Api::TwophaseSatParams>(state);

        Scalar Swe = TwophaseSat::Swe(state, Sw);

        return 0 <= Swe && Swe <= 1. - state.Snrei();
    }

    /*!
     * \brief Set the current absolute saturation for the
     *        current timestep
     */
    static void updateState(State &state, Scalar Sw)
    {
        Api::require<Api::ParkerLenhardState>(state);
        Api::require<Api::TwophaseSatParams>(state);

        if (Sw > 1 - 1e-5) {
            // if the absolute saturation is almost 1,
            // it means that we're back to the beginning
            reset(state);
            return;
        }

        Scalar Swe = TwophaseSat::Swe(state, Sw);
        ASSERT_RANGE(Swe, 0, 1);

        // find the loop number which corrosponds to the
        // given effective saturation
        ScanningCurve *curve = findScanningCurve_Swe_(state, Swe);

        // check whether the new reversal point is still in the
        // range where we use a spline instead of the actual
        // material law.
        if (curve->useSpline(Sw))
            return;

        Scalar Sw_app = Swapp_(state, Swe);

        // calculate the apparent saturation on the MIC and MDC
        // which yield the same capillary pressure as the
        // Sw at the current scanning curve
        Scalar pc = pC(state, Sw);
        Scalar Sw_mic = CapPressure::Sw(state.micParams(), pc);
        Scalar Sw_mdc = CapPressure::Sw(state.mdcParams(), pc);

        curve->setNext(Swe, pc, Sw_app,
                       Sw_mic, Sw_mdc);
        if (!curve->next())
            return;

        state.setCsc(curve);

        // if we're back on the MDC, we also have a new PISC!
        if (state.csc() == state.mdc()) {
            state.setPisc(state.mdc()->next());
            state.setSnrei(Snrei_(state, Swe));
        }

        // add a spline to the newly created reversal point
        if (useSplines)
            addSpline_(state, curve->next());
    }


    /*!
     * \brief Returns the capillary pressure dependend on
     *        the absolute__ saturation of the wetting phase.
     */
    static Scalar pC(const State &state, Scalar Sw)
    {
        Api::require<Api::ParkerLenhardParams>(state);
        Api::require<Api::TwophaseSatParams>(state);

        Scalar Swe = TwophaseSat::Swe(state, Sw);

        // if the effective saturation is smaller than zero we let
        // the underlying material law decide what to do.
        if (Swe < 0) {
            return CapPressure::pC(state.mdcParams(), Swe);
        }

        // calculate the current apparent saturation
        ScanningCurve *sc = findScanningCurve_Swe_(state, Swe);

        // check whether we are too close to the reveersal point
        // that we a spline to get a C1 continous function
        if (sc == state.mdc() && sc->prev()->useSpline(Sw))
            return sc->prev()->spline()->eval(Sw);
        if (sc->useSpline(Sw))
            return sc->spline()->eval(Sw);

        // calculate the apparant saturation
        Scalar Sw_app = Swapp_(state, Swe);

        // if the apparent saturation exceeds the 'legal' limits,
        // we also the underlying material law decide what to do.
        if (Sw_app > 1) {
            return CapPressure::pC(state.mdcParams(), Sw_app);
        }

        // put the effective saturation into the capillary pressure model
        Scalar pos = (Sw_app - sc->Sw_app())/(sc->prev()->Sw_app() - sc->Sw_app());
        if (sc->isImbib()) {
            Scalar SwMic =
                pos * (sc->prev()->SwMic() - sc->SwMic())
                + sc->SwMic();

            return CapPressure::pC(state.micParams(), SwMic);
        }
        else { // sc->isDrain()
            Scalar SwMdc =
                pos*(sc->prev()->SwMdc() - sc->SwMdc())
                + sc->SwMdc();

            return CapPressure::pC(state.mdcParams(), SwMdc);
        }
    }

    /*!
     * \brief Returns the absolute saturation to a given on the
     *        capillary pressure.
     */
    static Scalar Sw(const State &state, Scalar pC)
    {
        Api::require<Api::ParkerLenhardParams>(state);
        Api::require<Api::TwophaseSatParams>(state);

        assert(pC >= 0);

        // find the relevant scanning curve based on
        // the given capillary pressure
        ScanningCurve *sc = findScanningCurve_pC_(state, pC);
        // if we're on the MDC, we're already done...
        if (sc == state.mdc()) {
            Scalar Swe = CapPressure::Sw(state.mdcParams(), pC);
            return TwophaseSat::Sw(state, Swe);
        }

        // like in the forward direction, we first calculate
        // the relative position on the scanning curve sc
        Scalar pos;
        if (sc->isImbib()) {
            Scalar SwMic = CapPressure::Sw(state.micParams(), pC);
            pos = (SwMic - sc->SwMic())/(sc->prev()->SwMic() - sc->SwMic());
        }
        else {
            Scalar SwMdc = CapPressure::Sw(state.mdcParams(), pC);
            pos = (SwMdc - sc->SwMdc())/(sc->prev()->SwMdc() - sc->SwMdc());
        }

        // now we can calculate the apparent saturation
        Scalar Sw_app = pos*(sc->prev()->Sw_app() - sc->Sw_app()) + sc->Sw_app();

        // which can finally converted into an absolute saturation
        return TwophaseSat::Sw(state, SweFromSwapp_(state, Sw_app));
    }

    /*!
     * \brief The derivative of the capillary pressure regarding
     *        the absolute saturation.
     */
    static Scalar dpC_dSw(const State &state, Scalar Sw)
    {
        Api::require<Api::ParkerLenhardParams>(state);
        Api::require<Api::TwophaseSatParams>(state);

        Scalar Swe = TwophaseSat::Swe(state, Sw);

        // if the effective saturation exceeds 1, we have a
        // problem, but we let the underlying capillary pressure
        // model deal with it. it might apply a regularization,
        // or it might crash ;)
        if (Swe > 1)
            return CapPressure::dpC_dSw(state.mdcParams(), Swe);

        // calculate the current apparent saturation
        ScanningCurve *sc = findScanningCurve_Swe_(state, Swe);
        if (sc->useSpline(Sw))
            return sc->spline()->evalDerivative(Sw);

        Scalar Sw_app = Swapp_(state, Swe);
        ASSERT_RANGE(Sw_app, 0, 1);

        // calculate the derivative of the linear interpolation parameter
        // in regard to the apparent saturation
        Scalar pos = (Sw_app - sc->Sw_app())/(sc->prev()->Sw_app() - sc->Sw_app());
        Scalar dpos_dSwapp = 1/(sc->prev()->Sw_app() - sc->Sw_app());

        ASSERT_RANGE(pos, 0, 1);

        if (sc->isImbib()) {
            Scalar SwMic =
                pos * (sc->prev()->SwMic() - sc->SwMic())
                + sc->SwMic();
            // the factor behind the pos variable is a constant
            Scalar dSwMic_dSwapp = dpos_dSwapp *
                (sc->prev()->SwMic() - sc->SwMic());
            ASSERT_RANGE(SwMic, 0, 1);

            // inner times outer derivative (-> chain rule)
            return CapPressure::dpC_dSw(state.micParams(), SwMic)*
                dSwMic_dSwapp*
                dSwapp_dSwe_(state, sc) *
                TwophaseSat::dSwe_dSw(state);
        }
        else { // sc->isDrain()
            Scalar SwMdc =
                pos*(sc->prev()->SwMdc() - sc->SwMdc())
                + sc->SwMdc();
            // the factor behind the pos variable is a constant
            Scalar dSwMdc_dSwapp = dpos_dSwapp *
                (sc->prev()->SwMdc() - sc->SwMdc());
            ASSERT_RANGE(SwMdc, 0, 1);

            // inner times outer derivative (-> chain rule)
            return CapPressure::dpC_dSw(state.mdcParams(), SwMdc)*
                dSwMdc_dSwapp *
                dSwapp_dSwe_(state, sc) *
                TwophaseSat::dSwe_dSw(state);
        }
    }

    /*!
     * \brief The derivative of the absolute saturation regarding
     *        the capillary pressure.
     */
    static Scalar dSw_dpC (const State &state, Scalar pC)
    {
        Api::require<Api::ParkerLenhardParams>(state);
        Api::require<Api::TwophaseSatParams>(state);

#ifdef PL_STRICT_ASSERTS
        assert(pC >= 0);
#else
        if (pC < 0) pC = 0;
#endif

        // find the relevant scanning curve based on
        // the given capillary pressure
        ScanningCurve *sc = findScanningCurve_pC_(state, pC);
        // if we're on the MDC, we're already done...
        if (sc == state.mdc()) {
            return CapPressure::dSw_dpC(state.mdcParams(), pC)*TwophaseSat::dSw_dSwe(state);
        }

        // like in the forward direction, we first calculate
        // the relative position on the scanning curve sc
        //            Scalar pos;
        Scalar dpos_dSwMc;
        Scalar dSwMc_dpC;
        if (sc->isImbib()) {
            //                Scalar SwMic = state.mic().Sw(pC);
            //                pos = (SwMic - sc->SwMic())/(sc->prev()->SwMic() - sc->SwMic());
            dpos_dSwMc = 1/(sc->prev()->SwMic() - sc->SwMic());
            dSwMc_dpC = CapPressure::dSw_dpC(state.micParams(), pC);
        }
        else {
            //                Scalar SwMdc = state.mdc().Sw(pC);
            //                pos = (SwMdc - sc->SwMdc())/(sc->prev()->SwMdc() - sc->SwMdc());
            dpos_dSwMc = 1/(sc->prev()->SwMdc() - sc->SwMdc());
            dSwMc_dpC = CapPressure::dSw_dpC(state.mdcParams(), pC);
        }

        // now we can calculate the apparent saturation
        //            Scalar Sw_app = pos*(sc->prev()->Sw_app() - sc->Sw_app()) + sc->Sw_app();
        Scalar dSwapp_dpos = sc->prev()->Sw_app() - sc->Sw_app();
        Scalar dSwapp_dSwMc = dSwapp_dpos*dpos_dSwMc;

        return dSwMc_dpC*dSwapp_dSwMc*
            dSwe_dSwapp_(state, sc)*TwophaseSat::dSw_dSwe(state);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium.
     */
    static Scalar krw(const State &state, Scalar Sw)
    {
        Api::require<Api::ParkerLenhardParams>(state);
        Api::require<Api::TwophaseSatParams>(state);

        ASSERT_RANGE(Sw, 0, 1);

        // for the effective permeability we use play-type
        // hystersis. (That's because it's basically impossible
        // to invert the krw-Sw realtion effectivly.)
        Scalar Swe = TwophaseSat::Swe(state, Sw);

        Scalar Sw_app = Swapp_(state, Swe);
        ASSERT_RANGE(Sw_app, 0, 1);
        return CapPressure::krw(state.mdcParams(), Sw_app);

#if 0 // TODO: saturation-permebility hysteresis
        if (state.pisc() && Swe > state.csc()->Swe()) {
            return CapPressure::krw(state.micParams(), Sw_app);
        }
        else { // sc->isDrain()
            return CapPressure::krw(state.mdcParams(), Sw_app);
        }
#endif
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the params.
     */
    static Scalar krn(const State &state, Scalar Sw)
    {
        Api::require<Api::ParkerLenhardParams>(state);
        Api::require<Api::TwophaseSatParams>(state);

        ASSERT_RANGE(Sw, 0, 1);

        // for the effective permeability we use play-type
        // hystersis. (That's because it's basically impossible
        // to invert the krw-Sw realtion effectivly.)
        Scalar Swe = TwophaseSat::Swe(state, Sw);

        Scalar Sw_app = Swapp_(state, Swe);
        ASSERT_RANGE(Sw_app, 0, 1);
        return CapPressure::krn(state.mdcParams(), Sw_app);

#if 0 // TODO: saturation-permebility hysteresis
        if (state.pisc() && Swe > state.csc()->Swe()) {
            return CapPressure::krn(state.micParams(), Sw_app);
        }
        else { // sc->isDrain()
            return CapPressure::krn(state.mdcParams(), Sw_app);
        }
#endif
    }

    static Scalar SwToSwapp(const State &state, Scalar Sw)
    {
        return Swapp_(state, TwophaseSat::Swe(state, Sw));
    }
private:
    // find the loop on which the an effective
    // saturation has to be
    static ScanningCurve *findScanningCurve_Swe_(const State &state, Scalar Swe)
    {
        if (state.pisc() == NULL || Swe <= state.pisc()->Swe()) {
            // we don't have a PISC yet, or the effective
            // saturation is smaller than the saturation where the
            // PISC begins. In this case are on the MDC
            return state.mdc();
        }

        // if we have a primary imbibition curve, and our current
        // effective saturation is higher than the beginning of
        // the secondary drainage curve. this means we are on the
        // PISC again.
        if (state.pisc()->next() == NULL ||
            state.pisc()->next()->Swe() < Swe)
        {
            return state.pisc();
        }

        ScanningCurve *curve = state.csc()->next();
        while (true) {
            assert(curve != state.mdc()->prev());
            if (curve->isValidAt_Swe(Swe)) {
                return curve;
            }
            curve = curve->prev();
        }
    }

    // find the loop on which an capillary pressure belongs to
    static ScanningCurve *findScanningCurve_pC_(const State &state, Scalar pC)
    {
        if (state.mdc()->next() == NULL) {
            // we don't have a PISC yet,
            // so we must be on the MDC
            // (i.e. csc_ == mdc_)
            return state.mdc();
        }

        ScanningCurve *curve = state.csc()->next();
        while (true) {
            assert(curve != state.mdc()->prev());
            if (curve->isValidAt_pC(pC))
                return curve;
            curve = curve->prev();
        }
    }

    // calculate and save Snrei_
    static Scalar Snrei_(State &state, Scalar Swei)
    {
        if (Swei > 1 || Swei < 0)
            return state.Snrei();

        Scalar Snrei;
        if (state.Snre() == 0.0) {
            return 0.0;
        }
        else {
            // use Land's law
            Scalar R = 1.0/state.Snre() - 1;
            Snrei = (1 - Swei)/(1 + R*(1 - Swei));
        }

        // if the effective saturation is smaller or equal 100%,
        // the current trapped saturation must be smaller than the
        // residual saturation
        assert(state.Snrei() <= state.Snre());

        // we need to make sure that there is sufficent "distance"
        // between Swei and 1-Snrei in order not to get very steep
        // slopes which cause terrible nummeric headaches
        Snrei = std::min(Snrei, (Scalar) 1 - (Swei + 5e-2));
        Snrei = std::max(Snrei, (Scalar) 0.0);

        return Snrei;
    };

    // returns the trapped effective saturation at j
    static Scalar Snrij_(const State &state, Scalar Swej)
    {
        return state.Snrei()*(Swej - state.pisc()->Swe()) /
            (1 - state.Snrei() - state.pisc()->Swe());
    };

    // returns the apparent saturation of the
    // wetting phase depending on the effective saturation
    static Scalar Swapp_(const State &state, Scalar Swe)
    {
        if (state.pisc() == NULL || Swe <= state.pisc()->Swe()) {
            // we are on the main drainage curve, i.e.
            // no non-wetting fluid is trapped
            // -> apparent saturation == effective saturation
            return Swe;
        }


        // we are on a imbibition or drainage curve
        // which is not the main drainage curve
        // -> apparent saturation ==
        //    effective saturation + trapped effective saturation
        return Swe + Snrij_(state, Swe);
    };

    // Returns the effective saturation to a given apparent one
    static Scalar SweFromSwapp_(const State &state, Scalar Sw_app)
    {
        if (state.pisc() == NULL || Sw_app <= state.pisc()->Swe()) {
            // we are on the main drainage curve, i.e.
            // no non-wetting fluid is trapped
            // -> apparent saturation == effective saturation
            return Sw_app;
        }

        return (Sw_app*(1 - state.Snrei()
                        - state.pisc()->Swe())
                + state.pisc()->Swe()*state.Snrei())
            /(1 - state.pisc()->Swe());
    };

    // returns the derivative of the apparent saturation in
    // regard to the effective one
    static Scalar dSwapp_dSwe_(const State &state, ScanningCurve *sc)
    {
        if (sc == state.mdc())
            return 1.0;

        return 1 + state.Snrei()/(1 - state.Snrei()
                                  - state.pisc()->Swe());
    }

    // returns the derivative of the apparent saturation in
    // regard to the effective one
    static Scalar dSwe_dSwapp_(const State &state, ScanningCurve *sc)
    {
        if (sc == state.mdc())
            return 1.0;

        return (1 - state.Snrei() - state.pisc()->Swe())
            / (1 - state.pisc()->Swe());
    }

    // calculate a spline for the reversal point of the curve
    // in order to get a C1 steady function
    static void addSpline_(State &state, ScanningCurve *curve)
    {
        Scalar Sw = TwophaseSat::Sw(state, curve->Swe());
        if (Sw < 0.20)
            return;
        Scalar x1;
        Scalar x2;
        if (!curve->prev()) {
            x1 = 0.0;
            x2 = 8e-2;
        }
        else {
            Scalar deltaSwe = std::min(8e-2,
                                       fabs(curve->Swe() - curve->prev()->Swe())/3);
            if (curve->isImbib()) {
                Scalar SweTmp = TwophaseSat::Swe(state, Sw + deltaSwe);
                //                        Scalar curveRange = (1.0 - state.Snrei()) - curve->Swe();
                //                        SweTmp = std::min(SweTmp, curve->Swe() + curveRange/2);
                if (curve == state.pisc())
                    SweTmp = std::min(SweTmp, 1.0 - state.Snrei());
                else
                    SweTmp = std::min(SweTmp, curve->prev()->Swe());

                x1 = Sw;
                x2 = TwophaseSat::Sw(state, SweTmp);
            }
            else {
                Scalar SweTmp = TwophaseSat::Swe(state, Sw - deltaSwe);
                SweTmp = std::max(SweTmp, curve->prev()->Swe());

                x1 = TwophaseSat::Sw(state, SweTmp);
                x2 = Sw;
            }
        }

        if (fabs(x2 - x1) < 2e-2 ||
            curve->prev()->useSpline(x1) ||
            curve->prev()->useSpline(x2))
            return;


        Scalar y1 = pC(state, x1);
        Scalar y2 = pC(state, x2);
        Scalar m1 = dpC_dSw(state, x1);
        Scalar m2 = dpC_dSw(state, x2);
        Spline *spline = new Spline(x1, x2, y1, y2, m1, m2);
        curve->setSpline(spline);

        /*
          if (spline->p() > 2000) {
          std::cerr << boost::format("BLA: p: %f, x1: %f, x2: %f, y1: %f, y2: %f, m1: %f, m2: %f\n")
          %spline->p()%x1%x2%y1%y2%m1%m2;

          while (curve) {
          std::cerr << "loopnum:" << curve->loopNum() << ", Sw:" << TwophaseSat::Sw(state, curve->Swe()) << "\n";
          curve = curve->prev();
          }
          //                    exit(1);
          };
        */
    };
};


#undef ASSERT_RANGE
}; // namespace Dune

#endif // PARKER_LENHARD_HH
