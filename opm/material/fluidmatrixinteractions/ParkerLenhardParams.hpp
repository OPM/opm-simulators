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
 * \copydoc Opm::ParkerLenhardParams
 */
#ifndef OPM_PARKER_LENHARD_PARAMS_HPP
#define OPM_PARKER_LENHARD_PARAMS_HPP

#include <opm/material/fluidmatrixinteractions/RegularizedVanGenuchten.hpp>

namespace Opm
{
// forward declaration
template <class ScalarT>
class PLScanningCurve;

/*!
 * \brief Default parameter class for the Parker-Lenhard hysteresis
 *        model.
 */
template <class TraitsT>
class ParkerLenhardParams
{
public:
    typedef typename TraitsT::Scalar Scalar;
    typedef Opm::RegularizedVanGenuchten<TraitsT> VanGenuchten;
    typedef typename VanGenuchten::Params VanGenuchtenParams;
    typedef PLScanningCurve<Scalar> ScanningCurve;

    ParkerLenhardParams()
    {
        currentSnr_ = 0;
        mdc_ = new ScanningCurve(/*Swr=*/0);
        pisc_ = csc_ = NULL;

#ifndef NDEBUG
        finalized_ = false;
#endif
    }

    ParkerLenhardParams(const ParkerLenhardParams &p)
    {
        currentSnr_ = 0;
        SwrPc_ = p.SwrPc();
        mdc_ = new ScanningCurve(p.SwrPc());
        pisc_ = csc_ = NULL;

#ifndef NDEBUG
        finalized_ = p.finalized_;
#endif
    }

    ~ParkerLenhardParams()
    { delete mdc_; }

    /*!
     * \brief Calculate all dependent quantities once the independent
     *        quantities of the parameter object have been set.
     */
    void finalize()
    {
#ifndef NDEBUG
        finalized_ = true;
#endif
    }

    /*!
     * \brief Returns the parameters of the main imbibition curve (which uses
     *        the van Genuchten capillary pressure model).
     */
    const VanGenuchtenParams &micParams() const
    { assertFinalized_(); return *micParams_; }

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
    { assertFinalized_(); return *mdcParams_; }

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
    { assertFinalized_(); return Snr_; }

    /*!
     * \brief Set the  non-wetting phase residual saturation.
     */
    void setSnr(Scalar val)
    { Snr_ = val; }

    /*!
     * \brief Returns wetting phase residual saturation for the capillary pressure curve.
     */
    Scalar SwrPc() const
    { assertFinalized_(); return SwrPc_; }

    /*!
     * \brief Returns wetting phase residual saturation for the residual saturation curves.
     */
    Scalar SwrKr() const
    { assertFinalized_(); return SwrKr_; }

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
    { assertFinalized_(); return currentSnr_; }

    /*!
     * \brief Set the current effective residual saturation.
     */
    void setCurrentSnr(Scalar val)
    { currentSnr_ = val; }

    /*!
     * \brief Returns the main drainage curve
     */
    ScanningCurve *mdc() const
    { assertFinalized_(); return mdc_; }

    /*!
     * \brief Set the main drainage curve.
     */
    void setMdc(ScanningCurve *val)
    { mdc_ = val; }

    /*!
     * \brief Returns the primary imbibition scanning curve
     */
    ScanningCurve *pisc() const
    { assertFinalized_(); return pisc_; }

    /*!
     * \brief Set the primary imbibition scanning curve.
     */
    void setPisc(ScanningCurve *val)
    { pisc_ = val; }

    /*!
     * \brief Returns the current scanning curve
     */
    ScanningCurve *csc() const
    { assertFinalized_(); return csc_; }

    /*!
     * \brief Set the current scanning curve.
     */
    void setCsc(ScanningCurve *val)
    { csc_ = val; }

private:
#ifndef NDEBUG
    void assertFinalized_() const
    { assert(finalized_); }

    bool finalized_;
#else
    void assertFinalized_() const
    { }
#endif

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
} // namespace Opm

#endif
