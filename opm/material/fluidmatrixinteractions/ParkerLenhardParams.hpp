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
 * \copydoc Opm::ParkerLenhardParams
 */
#ifndef OPM_PARKER_LENHARD_PARAMS_HPP
#define OPM_PARKER_LENHARD_PARAMS_HPP

#include <opm/material/fluidmatrixinteractions/RegularizedVanGenuchten.hpp>
#include <opm/material/common/EnsureFinalized.hpp>

#include <cassert>

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
class ParkerLenhardParams : public EnsureFinalized
{
public:
    using EnsureFinalized :: finalize;

    typedef typename TraitsT::Scalar Scalar;
    typedef Opm::RegularizedVanGenuchten<TraitsT> VanGenuchten;
    typedef typename VanGenuchten::Params VanGenuchtenParams;
    typedef PLScanningCurve<Scalar> ScanningCurve;

    ParkerLenhardParams()
    {
        currentSnr_ = 0;
        mdc_ = new ScanningCurve(/*Swr=*/0);
        pisc_ = csc_ = NULL;
    }

    ParkerLenhardParams(const ParkerLenhardParams& p)
        : EnsureFinalized( p )
    {
        currentSnr_ = 0;
        SwrPc_ = p.SwrPc_;
        mdc_ = new ScanningCurve(SwrPc_);
        pisc_ = csc_ = NULL;
    }

    ~ParkerLenhardParams()
    { delete mdc_; }

    /*!
     * \brief Returns the parameters of the main imbibition curve (which uses
     *        the van Genuchten capillary pressure model).
     */
    const VanGenuchtenParams& micParams() const
    { EnsureFinalized::check(); return *micParams_; }

    /*!
     * \brief Sets the parameters of the main imbibition curve (which uses
     *        the van Genuchten capillary pressure model).
     */
    void setMicParams(const VanGenuchtenParams* val)
    { micParams_ = val; }

    /*!
     * \brief Returns the parameters of the main drainage curve (which uses
     *        the van Genuchten capillary pressure model).
     */
    const VanGenuchtenParams& mdcParams() const
    { EnsureFinalized::check(); return *mdcParams_; }

    /*!
     * \brief Sets the parameters of the main drainage curve (which uses
     *        the van Genuchten capillary pressure model).
     */
    void setMdcParams(const VanGenuchtenParams* val)
    { mdcParams_ = val; }

    /*!
     * \brief Returns non-wetting phase residual saturation.
     */
    Scalar Snr() const
    { EnsureFinalized::check(); return Snr_; }

    /*!
     * \brief Set the  non-wetting phase residual saturation.
     */
    void setSnr(Scalar val)
    { Snr_ = val; }

    /*!
     * \brief Returns wetting phase residual saturation for the capillary pressure curve.
     */
    Scalar SwrPc() const
    { EnsureFinalized::check(); return SwrPc_; }

    /*!
     * \brief Returns wetting phase residual saturation for the residual saturation curves.
     */
    Scalar SwrKr() const
    { EnsureFinalized::check(); return SwrKr_; }

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
    { EnsureFinalized::check(); return currentSnr_; }

    /*!
     * \brief Set the current effective residual saturation.
     */
    void setCurrentSnr(Scalar val)
    { currentSnr_ = val; }

    /*!
     * \brief Returns the main drainage curve
     */
    ScanningCurve* mdc() const
    { EnsureFinalized::check(); return mdc_; }

    /*!
     * \brief Set the main drainage curve.
     */
    void setMdc(ScanningCurve* val)
    { mdc_ = val; }

    /*!
     * \brief Returns the primary imbibition scanning curve
     */
    ScanningCurve* pisc() const
    { EnsureFinalized::check(); return pisc_; }

    /*!
     * \brief Set the primary imbibition scanning curve.
     */
    void setPisc(ScanningCurve* val)
    { pisc_ = val; }

    /*!
     * \brief Returns the current scanning curve
     */
    ScanningCurve* csc() const
    { EnsureFinalized::check(); return csc_; }

    /*!
     * \brief Set the current scanning curve.
     */
    void setCsc(ScanningCurve* val)
    { csc_ = val; }

private:
    const VanGenuchtenParams* micParams_;
    const VanGenuchtenParams* mdcParams_;
    Scalar SwrPc_;
    Scalar SwrKr_;
    Scalar Snr_;
    Scalar currentSnr_;
    mutable ScanningCurve* mdc_;
    mutable ScanningCurve* pisc_;
    mutable ScanningCurve* csc_;
};
} // namespace Opm

#endif
