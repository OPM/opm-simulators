// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
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
 *
 * \brief An implementation of the Parker-Lenhard capillary pressure
 *        hystersis model.
 */
#ifndef DUMUX_PARKER_LENHARD_HH
#define DUMUX_PARKER_LENHARD_HH

#include "parkerlenhardparams.hh"

#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

#include <algorithm>
#include <iostream>
#include <cassert>

namespace Dumux
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
                                    -1, // loop number
                                    0, // Sw
                                    1e12, // pC
                                    0, // Sw_app
                                    0, // SwMic
                                    0);   // SwMdc
        next_ = NULL;

        Swe_ = 1.0;
        pC_ = 0.0;
        Sw_app_ = 1.0;
        SwMic_ = 1.0;
        SwMdc_ = 1.0;
    }

protected:
    PLScanningCurve(PLScanningCurve *prev,
                    PLScanningCurve *next,
                    int loopN,
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
        pC_ = pC;
        Sw_app_ = Sw_app;
        SwMic_ = SwMic;
        SwMdc_ = SwMdc;
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
    }

    /*!
     * \brief Returns true iff a given capillary pressure
     *        pC is within the scope of the curve, i.e.
     *        whether pC is part of the curve's
     *        image and the curve thus applies to pC.
     */
    bool isValidAt_pC(Scalar pCi)
    {
        if (isImbib())
            return prev_->pC() < pCi && pCi < pC();
        else
            return pC() < pCi && pCi < prev_->pC();
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

    int loopNum_;

    Scalar Swe_;
    Scalar pC_;
    Scalar Sw_app_;

    Scalar SwMdc_;
    Scalar SwMic_;
};

/*!
 * \ingroup material
 * \brief Implements the Parker-Lenhard twophase
 *        p_c-Sw hysteresis model. This class adheres to the twophase
 *        capillary pressure API.
 */
template <class ScalarT, class ParamsT = ParkerLenhardParams<ScalarT> >
class ParkerLenhard
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

private:
    typedef typename ParamsT::VanGenuchten VanGenuchten;
    typedef Dumux::PLScanningCurve<Scalar> ScanningCurve;

public:
    /*!
     * \brief Resets the hysteresis model to the
     *        initial parameters on the main drainage curve
     */
    static void reset(Params &params)
    {
        delete params.mdc(); // this will work even if mdc_ == NULL!
        params.setMdc(new ScanningCurve());
        params.setCsc(params.mdc());
        params.setPisc(NULL);
        params.setSnrei(0.0);
    }

    /*!
     * \brief Set the current absolute saturation for the
     *        current timestep
     */
    static void update(Params &params, Scalar Sw)
    {
        if (Sw > 1 - 1e-5) {
            // if the absolute saturation is almost 1,
            // it means that we're back to the beginning
            reset(params);
            return;
        }

        Scalar Swe = SwToSwe_(params, Sw);

        // find the loop number which corrosponds to the
        // given effective saturation
        ScanningCurve *curve = findScanningCurve_Swe_(params, Swe);

        Scalar Sw_app = Swapp_(params, Swe);

        // calculate the apparent saturation on the MIC and MDC
        // which yield the same capillary pressure as the
        // Sw at the current scanning curve
        Scalar pc = pC(params, Sw);
        Scalar Sw_mic = VanGenuchten::Sw(params.micParams(), pc);
        Scalar Sw_mdc = VanGenuchten::Sw(params.mdcParams(), pc);

        curve->setNext(Swe, pc, Sw_app,
                       Sw_mic, Sw_mdc);
        if (!curve->next())
            return;

        params.setCsc(curve);

        // if we're back on the MDC, we also have a new PISC!
        if (params.csc() == params.mdc()) {
            params.setPisc(params.mdc()->next());
            params.setSnrei(Snrei_(params, Swe));
        }
    }

    /*!
     * \brief Returns true iff the given absolute saturation
     *        is viable (i.e. larger than the residual
     *        saturation of the wetting phase but smaller
     *        than the maximum saturation of the current
     *        PISC)
     */
    static bool isValidSw(const Params &params, Scalar Sw)
    {
        Scalar Swe = SwToSwe_(params, Sw);

        return 0 <= Swe && Swe <= 1. - params.Snrei();
    }

    /*!
     * \brief Returns the capillary pressure dependend on
     *        the absolute__ saturation of the wetting phase.
     */
    static Scalar pC(const Params &params, Scalar Sw)
    {
        Scalar Swe = SwToSwe_(params, Sw);

        // if the effective saturation is smaller than zero we let
        // the underlying material law decide what to do.
        if (Swe < 0) {
            return VanGenuchten::pC(params.mdcParams(), Swe);
        }

        // calculate the current apparent saturation
        ScanningCurve *sc = findScanningCurve_Swe_(params, Swe);

        // calculate the apparant saturation
        Scalar Sw_app = Swapp_(params, Swe);

        // if the apparent saturation exceeds the 'legal' limits,
        // we also the underlying material law decide what to do.
        if (Sw_app > 1) {
            return 0.0; // VanGenuchten::pC(params.mdcParams(), Sw_app);
        }

        // put the effective saturation into the capillary pressure model
        Scalar pos = (Sw_app - sc->Sw_app())/(sc->prev()->Sw_app() - sc->Sw_app());
        if (sc->isImbib()) {
            Scalar SwMic =
                pos * (sc->prev()->SwMic() - sc->SwMic())
                + sc->SwMic();

            return VanGenuchten::pC(params.micParams(), SwMic);
        }
        else { // sc->isDrain()
            Scalar SwMdc =
                pos*(sc->prev()->SwMdc() - sc->SwMdc())
                + sc->SwMdc();

            return VanGenuchten::pC(params.mdcParams(), SwMdc);
        }
    }

    /*!
     * \brief Returns the absolute saturation to a given on the
     *        capillary pressure.
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        assert(pC >= 0);

        // find the relevant scanning curve based on
        // the given capillary pressure
        ScanningCurve *sc = findScanningCurve_pC_(params, pC);
        // if we're on the MDC, we're already done...
        if (sc == params.mdc()) {
            Scalar Swe = VanGenuchten::Sw(params.mdcParams(), pC);
            return SweToSw_(params, Swe);
        }

        // like in the forward direction, we first calculate
        // the relative position on the scanning curve sc
        Scalar pos;
        if (sc->isImbib()) {
            Scalar SwMic = VanGenuchten::Sw(params.micParams(), pC);
            pos = (SwMic - sc->SwMic())/(sc->prev()->SwMic() - sc->SwMic());
        }
        else {
            Scalar SwMdc = VanGenuchten::Sw(params.mdcParams(), pC);
            pos = (SwMdc - sc->SwMdc())/(sc->prev()->SwMdc() - sc->SwMdc());
        }

        // now we can calculate the apparent saturation
        Scalar Sw_app = pos*(sc->prev()->Sw_app() - sc->Sw_app()) + sc->Sw_app();

        // which can finally converted into an absolute saturation
        return SweToSw_(params, SweFromSwapp_(params, Sw_app));
    }

    /*!
     * \brief The derivative of the capillary pressure regarding
     *        the absolute saturation.
     */
    static Scalar dpC_dSw(const Params &params, Scalar Sw)
    {
        Scalar Swe = SwToSwe_(params, Sw);

        // if the effective saturation exceeds 1, we have a
        // problem, but we let the underlying capillary pressure
        // model deal with it. it might apply a regularization,
        // or it might crash ;)
        if (Swe > 1)
            return VanGenuchten::dpC_dSw(params.mdcParams(), Swe);

        // calculate the current apparent saturation
        ScanningCurve *sc = findScanningCurve_Swe_(params, Swe);
        Scalar Sw_app = Swapp_(params, Swe);

        // calculate the derivative of the linear interpolation parameter
        // in regard to the apparent saturation
        Scalar pos = (Sw_app - sc->Sw_app())/(sc->prev()->Sw_app() - sc->Sw_app());
        Scalar dpos_dSwapp = 1/(sc->prev()->Sw_app() - sc->Sw_app());

        if (sc->isImbib()) {
            Scalar SwMic =
                pos * (sc->prev()->SwMic() - sc->SwMic())
                + sc->SwMic();
            // the factor behind the pos variable is a constant
            Scalar dSwMic_dSwapp = dpos_dSwapp *
                (sc->prev()->SwMic() - sc->SwMic());

            // inner times outer derivative (-> chain rule)
            return VanGenuchten::dpC_dSw(params.micParams(), SwMic)*
                dSwMic_dSwapp*
                dSwapp_dSwe_(params, sc) *
                dSwe_dSw_(params);
        }
        else { // sc->isDrain()
            Scalar SwMdc =
                pos*(sc->prev()->SwMdc() - sc->SwMdc())
                + sc->SwMdc();
            // the factor behind the pos variable is a constant
            Scalar dSwMdc_dSwapp = dpos_dSwapp *
                (sc->prev()->SwMdc() - sc->SwMdc());

            // inner times outer derivative (-> chain rule)
            return VanGenuchten::dpC_dSw(params.mdcParams(), SwMdc)*
                dSwMdc_dSwapp *
                dSwapp_dSwe_(params, sc) *
                dSwe_dSw_(params);
        }
    }

    /*!
     * \brief The derivative of the absolute saturation regarding
     *        the capillary pressure.
     */
    static Scalar dSw_dpC (const Params &params, Scalar pC)
    {
        if (pC < 0) pC = 0;

        // find the relevant scanning curve based on
        // the given capillary pressure
        ScanningCurve *sc = findScanningCurve_pC_(params, pC);
        // if we're on the MDC, we're already done...
        if (sc == params.mdc()) {
            return VanGenuchten::dSw_dpC(params.mdcParams(), pC)*dSw_dSwe_(params);
        }

        // like in the forward direction, we first calculate
        // the relative position on the scanning curve sc
        //            Scalar pos;
        Scalar dpos_dSwMc;
        Scalar dSwMc_dpC;
        if (sc->isImbib()) {
            //                Scalar SwMic = params.mic().Sw(pC);
            //                pos = (SwMic - sc->SwMic())/(sc->prev()->SwMic() - sc->SwMic());
            dpos_dSwMc = 1/(sc->prev()->SwMic() - sc->SwMic());
            dSwMc_dpC = VanGenuchten::dSw_dpC(params.micParams(), pC);
        }
        else {
            //                Scalar SwMdc = params.mdc().Sw(pC);
            //                pos = (SwMdc - sc->SwMdc())/(sc->prev()->SwMdc() - sc->SwMdc());
            dpos_dSwMc = 1/(sc->prev()->SwMdc() - sc->SwMdc());
            dSwMc_dpC = VanGenuchten::dSw_dpC(params.mdcParams(), pC);
        }

        // now we can calculate the apparent saturation
        //            Scalar Sw_app = pos*(sc->prev()->Sw_app() - sc->Sw_app()) + sc->Sw_app();
        Scalar dSwapp_dpos = sc->prev()->Sw_app() - sc->Sw_app();
        Scalar dSwapp_dSwMc = dSwapp_dpos*dpos_dSwMc;

        return dSwMc_dpC*dSwapp_dSwMc*
            dSwe_dSwapp_(params, sc)*dSw_dSwe_(params);
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium.
     */
    static Scalar krw(const Params &params, Scalar Sw)
    {
        // for the effective permeability we use play-type
        // hystersis. (That's because it's basically impossible
        // to invert the krw-Sw realtion effectivly.)
        Scalar Swe = SwToSwe_(params, Sw);

        Scalar Sw_app = Swapp_(params, Swe);
        return VanGenuchten::krw(params.mdcParams(), Sw_app);

#if 0 // TODO: saturation-permebility hysteresis
        if (params.pisc() && Swe > params.csc()->Swe()) {
            return VanGenuchten::krw(params.micParams(), Sw_app);
        }
        else { // sc->isDrain()
            return VanGenuchten::krw(params.mdcParams(), Sw_app);
        }
#endif
    }

    /*!
     * \brief The relative permeability for the non-wetting phase
     *        of the params.
     */
    static Scalar krn(const Params &params, Scalar Sw)
    {
        // for the effective permeability we use play-type
        // hystersis. (That's because it's basically impossible
        // to invert the krw-Sw realtion effectivly.)
        Scalar Swe = SwToSwe_(params, Sw);

        Scalar Sw_app = Swapp_(params, Swe);
        return VanGenuchten::krn(params.mdcParams(), Sw_app);

#if 0 // TODO: saturation-permebility hysteresis
        if (params.pisc() && Swe > params.csc()->Swe()) {
            return VanGenuchten::krn(params.micParams(), Sw_app);
        }
        else { // sc->isDrain()
            return VanGenuchten::krn(params.mdcParams(), Sw_app);
        }
#endif
    }

private:
    /*!
     * \brief Convert an absolute wetting saturation to an apparent one.
     */
    static Scalar SwToSwapp(const Params &params, Scalar Sw)
    {
        return Swapp_(params, SwToSwe_(params, Sw));
    }

    /*!
     * \brief Convert an absolute wetting saturation to an effective one.
     *
     * \param Sw        Absolute saturation of the wetting phase \f${S}_w\f$.
     * \param params    A container object that is populated with the appropriate coefficients for the respective law.
     *                  Therefore, in the (problem specific) spatialParameters  first, the material law is chosen, and then the params container
     *                  is constructed accordingly. Afterwards the values are set there, too.
     * \return          Effective saturation of the wetting phase.
     */
    static Scalar SwToSwe_(const Params &params, Scalar Sw)
    { return (Sw - params.Swr())/(1 - params.Swr() - params.Snr()); }

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



    // find the loop on which the an effective
    // saturation has to be
    static ScanningCurve *findScanningCurve_Swe_(const Params &params, Scalar Swe)
    {
        if (params.pisc() == NULL || Swe <= params.pisc()->Swe()) {
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
            params.pisc()->next()->Swe() < Swe)
        {
            return params.pisc();
        }

        ScanningCurve *curve = params.csc()->next();
        while (true) {
            assert(curve != params.mdc()->prev());
            if (curve->isValidAt_Swe(Swe)) {
                return curve;
            }
            curve = curve->prev();
        }
    }

    // find the loop on which an capillary pressure belongs to
    static ScanningCurve *findScanningCurve_pC_(const Params &params, Scalar pC)
    {
        if (params.mdc()->next() == NULL) {
            // we don't have a PISC yet,
            // so we must be on the MDC
            // (i.e. csc_ == mdc_)
            return params.mdc();
        }

        ScanningCurve *curve = params.csc()->next();
        while (true) {
            assert(curve != params.mdc()->prev());
            if (curve->isValidAt_pC(pC))
                return curve;
            curve = curve->prev();
        }
    }

    // calculate and save Snrei_
    static Scalar Snrei_(Params &params, Scalar Swei)
    {
        if (Swei > 1 || Swei < 0)
            return params.Snrei();

        Scalar Snrei;
        if (params.Snre() == 0.0) {
            return 0.0;
        }
        else {
            // use Land's law
            Scalar R = 1.0/params.Snre() - 1;
            Snrei = (1 - Swei)/(1 + R*(1 - Swei));
        }

        // if the effective saturation is smaller or equal 100%,
        // the current trapped saturation must be smaller than the
        // residual saturation
        assert(params.Snrei() <= params.Snre());

        // we need to make sure that there is sufficent "distance"
        // between Swei and 1-Snrei in order not to get very steep
        // slopes which cause terrible nummeric headaches
        Snrei = std::min(Snrei, (Scalar) 1 - (Swei + 5e-2));
        Snrei = std::max(Snrei, (Scalar) 0.0);

        return Snrei;
    }

    // returns the trapped effective saturation at j
    static Scalar Snrij_(const Params &params, Scalar Swej)
    {
        return params.Snrei()*(Swej - params.pisc()->Swe()) /
            (1 - params.Snrei() - params.pisc()->Swe());
    }

    // returns the apparent saturation of the
    // wetting phase depending on the effective saturation
    static Scalar Swapp_(const Params &params, Scalar Swe)
    {
        if (params.pisc() == NULL || Swe <= params.pisc()->Swe()) {
            // we are on the main drainage curve, i.e.
            // no non-wetting fluid is trapped
            // -> apparent saturation == effective saturation
            return Swe;
        }


        // we are on a imbibition or drainage curve
        // which is not the main drainage curve
        // -> apparent saturation ==
        //    effective saturation + trapped effective saturation
        return Swe + Snrij_(params, Swe);
    }

    // Returns the effective saturation to a given apparent one
    static Scalar SweFromSwapp_(const Params &params, Scalar Sw_app)
    {
        if (params.pisc() == NULL || Sw_app <= params.pisc()->Swe()) {
            // we are on the main drainage curve, i.e.
            // no non-wetting fluid is trapped
            // -> apparent saturation == effective saturation
            return Sw_app;
        }

        return (Sw_app*(1 - params.Snrei()
                        - params.pisc()->Swe())
                + params.pisc()->Swe()*params.Snrei())
            /(1 - params.pisc()->Swe());
    }

    // returns the derivative of the apparent saturation in
    // regard to the effective one
    static Scalar dSwapp_dSwe_(const Params &params, ScanningCurve *sc)
    {
        if (sc == params.mdc())
            return 1.0;

        return 1 + params.Snrei()/(1 - params.Snrei()
                                  - params.pisc()->Swe());
    }

    // returns the derivative of the apparent saturation in
    // regard to the effective one
    static Scalar dSwe_dSwapp_(const Params &params, ScanningCurve *sc)
    {
        if (sc == params.mdc())
            return 1.0;

        return (1 - params.Snrei() - params.pisc()->Swe())
            / (1 - params.pisc()->Swe());
    }
};


}; // namespace Dumux

#endif // PARKER_LENHARD_HH
