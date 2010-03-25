/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
 * \file ParkerLenhardState.hh Specification of the state API for the
 *                           Parker-Lenhard hysteresis model.
 */
#ifndef PARKER_LENHARD_STATE_HH
#define PARKER_LENHARD_STATE_HH

#include <dumux/common/apis.hh>

#include "statehelpermacros.hh"
#include "twophasesatstate.hh"

#include <stdlib.h>

namespace Dumux
{
// forward declaration
template <class ScalarT>
class PLScanningCurve;


namespace Api
{
BEGIN_API_DEF(ParkerLenhardParams)
{
    typedef typename Implementation::Scalar Scalar;
    typedef typename Implementation::CapPressureParams PCParams;
    typedef PLScanningCurve<Scalar> ScanningCurve;

    require<TwophaseSatParams>(const_impl);

    Scalar tmp = 0.5;
    tmp = impl.Snrei();

    const PCParams *tmpPc;
    tmpPc = &const_impl.micParams();
    tmpPc = &const_impl.mdcParams();

    const ScanningCurve *tmpSC;
    tmpSC = const_impl.mdc();
    tmpSC = const_impl.pisc();
    tmpSC = const_impl.csc();
}
END_API_DEF;

BEGIN_API_DEF(ParkerLenhardState)
{
    typedef typename Implementation::Scalar Scalar;
    typedef typename Implementation::CapPressureParams PCParams;
    typedef PLScanningCurve<Scalar> ScanningCurve;

    require<ParkerLenhardParams>(impl);

    Scalar tmp = 0.5;
    impl.setSnrei(tmp);

    ScanningCurve *tmpSC = NULL;
    tmpSC = impl.mdc();   impl.setMdc(tmpSC);
    tmpSC = impl.pisc();  impl.setPisc(tmpSC);
    tmpSC = impl.csc();   impl.setCsc(tmpSC);
}
END_API_DEF;
}; // namespace Api


/*!
 * \brief A reference implementation of the state API for the
 *        Parker-Lenhard hysteresis model.
 */
template <class CapPressureParamsT>
class ParkerLenhardState : public TwophaseSatState<typename CapPressureParamsT::Scalar>
{
public:
    typedef CapPressureParamsT CapPressureParams;
    typedef typename CapPressureParams::Scalar Scalar;
    typedef Dumux::TwophaseSatState<Scalar> TwophaseSatState;
    typedef Dune::PLScanningCurve<Scalar>  ScanningCurve;

    ParkerLenhardState(Scalar Swr,
                       Scalar Snr,
                       const CapPressureParams &micParams,
                       const CapPressureParams &mdcParams)
        : TwophaseSatState(Swr, Snr),
          micParams_(micParams),
          mdcParams_(mdcParams)
    {
        Api::require<Api::ParkerLenhardState>(*this);

        Snrei_ = 0;
        mdc_ = new ScanningCurve();
        pisc_ = csc_ = NULL;
    }

    ~ParkerLenhardState()
    { delete mdc_; }

    /*!
     * \brief The parameters of the MIC for the capillary pressure
     *        model.
     */
    PROPERTY(CapPressureParams, micParams, setMicParams);

    /*!
     * \brief The parameters of the MDC for the capillary pressure
     *        model.
     */
    PROPERTY(CapPressureParams, mdcParams, setMdcParams);

    // current effective residual saturation
    PROPERTY(Scalar, Snrei, setSnrei);

    /*!
     * \brief The scanning curves
     */
    MUTABLE_PTR_PROPERTY(ScanningCurve, mdc, setMdc); // main drainage
    MUTABLE_PTR_PROPERTY(ScanningCurve, pisc, setPisc); // primary imbibition
    MUTABLE_PTR_PROPERTY(ScanningCurve, csc, setCsc); // current
};
}; // namespace Dumux

#endif
