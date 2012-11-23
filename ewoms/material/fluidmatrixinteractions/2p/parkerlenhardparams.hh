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
 * \copydoc Ewoms::ParkerLenhardParams
 */
#ifndef EWOMS_PARKER_LENHARD_PARAMS_HH
#define EWOMS_PARKER_LENHARD_PARAMS_HH

#include <ewoms/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>

namespace Ewoms
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
    typedef Ewoms::RegularizedVanGenuchten<Scalar> VanGenuchten;
    typedef typename VanGenuchten::Params VanGenuchtenParams;
    typedef PLScanningCurve<Scalar> ScanningCurve;

    ParkerLenhardParams()
    {
        currentSnr_ = 0;
        mdc_ = new ScanningCurve(/*Swr=*/0);
        pisc_ = csc_ = NULL;
    }

    ParkerLenhardParams(const ParkerLenhardParams &p)
    {
        currentSnr_ = 0;
        SwrPc_ = p.SwrPc();
        mdc_ = new ScanningCurve(p.SwrPc());
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
     * \brief Returns wetting phase residual saturation for the capillary pressure curve.
     */
    Scalar SwrPc() const
    { return SwrPc_; }

    /*!
     * \brief Returns wetting phase residual saturation for the residual saturation curves.
     */
    Scalar SwrKr() const
    { return SwrKr_; }

    /*!
     * \brief Set the wetting phase residual saturation for the
     *        capillary pressure and the relative permeabilities.
     */
    void setSwr(Scalar pcSwr, Scalar krSwr = -1)
    {
        SwrPc_ = pcSwr;
        if (krSwr < 0)
            SwrKr_ = pcSwr;
        else
            SwrKr_ = krSwr;
    }

    /*!
     * \brief Returns the current effective residual saturation.
     */
    Scalar currentSnr() const
    { return currentSnr_; }

    /*!
     * \brief Set the current effective residual saturation.
     */
    void setCurrentSnr(Scalar val)
    { currentSnr_ = val; }

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
    Scalar SwrPc_;
    Scalar SwrKr_;
    Scalar Snr_;
    Scalar currentSnr_;
    mutable ScanningCurve *mdc_;
    mutable ScanningCurve *pisc_;
    mutable ScanningCurve *csc_;
};
} // namespace Ewoms

#endif
