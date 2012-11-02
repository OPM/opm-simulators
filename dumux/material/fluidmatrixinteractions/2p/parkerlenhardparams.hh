// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
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
 * \copydoc Dumux::ParkerLenhardParams
 */
#ifndef DUMUX_PARKER_LENHARD_PARAMS_HH
#define DUMUX_PARKER_LENHARD_PARAMS_HH

#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>

namespace Dumux
{
// forward declaration
template <class ScalarT>
class PLScanningCurve;

/*!
 * \brief Default parameter class for the Parker-Lenhard hysteresis
 *        model.
 */
template <class ScalarT>
class ParkerLenhardParams
{
public:
    typedef ScalarT Scalar;
    typedef Dumux::RegularizedVanGenuchten<Scalar> VanGenuchten;
    typedef typename VanGenuchten::Params VanGenuchtenParams;
    typedef PLScanningCurve<Scalar> ScanningCurve;

    ParkerLenhardParams()
    {
        Snrei_ = 0;
        mdc_ = new ScanningCurve();
        pisc_ = csc_ = NULL;
    }

    ParkerLenhardParams(const ParkerLenhardParams &)
    {
        Snrei_ = 0;
        mdc_ = new ScanningCurve();
        pisc_ = csc_ = NULL;
    }

    ~ParkerLenhardParams()
    { delete mdc_; }

    /*!
     * \brief Returns the parameters of the main imbibition curve (which uses
     *        the van Genuchten capillary pressure model).
     */
    const VanGenuchtenParams &micParams() const
    { return *micParams_; }

    /*!
     * \brief Sets the parameters of the main imbibition curve (which uses
     *        the van Genuchten capillary pressure model).
     */
    void setMicParams(const VanGenuchtenParams *val)
    { micParams_ = val; }

    /*!
     * \brief Returns the parameters of the main drainage curve (which uses
     *        the van Genuchten capillary pressure model).
     */
    const VanGenuchtenParams &mdcParams() const
    { return *mdcParams_; }

    /*!
     * \brief Sets the parameters of the main drainage curve (which uses
     *        the van Genuchten capillary pressure model).
     */
    void setMdcParams(const VanGenuchtenParams *val)
    { mdcParams_ = val; }

    /*!
     * \brief Returns non-wetting phase residual saturation.
     */
    Scalar Snr() const
    { return Snr_; }

    /*!
     * \brief Set the  non-wetting phase residual saturation.
     */
    void setSnr(Scalar val)
    { Snr_ = val; }

    /*!
     * \brief Returns wetting phase residual saturation.
     */
    Scalar Swr() const
    { return Swr_; }

    /*!
     * \brief Set the  wetting phase residual saturation.
     */
    void setSwr(Scalar val)
    { Swr_ = val; }

    /*!
     * \brief Returns the current effective residual saturation.
     */
    Scalar Snrei() const
    { return Snrei_; }

    /*!
     * \brief Set the current effective residual saturation.
     */
    void setSnrei(Scalar val)
    { Snrei_ = val; }

    /*!
     * \brief Returns the effective residual saturation of the non-wetting phase.
     */
    Scalar Snre() const
    { return Snre_; }

    /*!
     * \brief Set the effective residual saturation of the non-wetting phase.
     */
    void setSnre(Scalar val)
    { Snre_ = val; }

    /*!
     * \brief Returns the main drainage curve
     */
    ScanningCurve *mdc() const
    { return mdc_; }

    /*!
     * \brief Set the main drainage curve.
     */
    void setMdc(ScanningCurve *val)
    { mdc_ = val; }

    /*!
     * \brief Returns the primary imbibition scanning curve
     */
    ScanningCurve *pisc() const
    { return pisc_; }

    /*!
     * \brief Set the primary imbibition scanning curve.
     */
    void setPisc(ScanningCurve *val)
    { pisc_ = val; }

    /*!
     * \brief Returns the current scanning curve
     */
    ScanningCurve *csc() const
    { return csc_; }

    /*!
     * \brief Set the current scanning curve.
     */
    void setCsc(ScanningCurve *val)
    { csc_ = val; }


private:
    const VanGenuchtenParams *micParams_;
    const VanGenuchtenParams *mdcParams_;
    Scalar Swr_;
    Scalar Snr_;
    Scalar Snre_;
    Scalar Snrei_;
    mutable ScanningCurve *mdc_;
    mutable ScanningCurve *pisc_;
    mutable ScanningCurve *csc_;
};
} // namespace Dumux

#endif
