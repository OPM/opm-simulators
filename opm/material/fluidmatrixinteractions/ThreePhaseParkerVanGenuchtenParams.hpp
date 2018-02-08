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
 * \copydoc Opm::ThreePhaseParkerVanGenuchtenParams
 */
#ifndef OPM_THREE_PHASE_PARKER_VAN_GENUCHTEN_PARAMS_HPP
#define OPM_THREE_PHASE_PARKER_VAN_GENUCHTEN_PARAMS_HPP

#include <dune/common/fvector.hh>

#include <opm/material/common/EnsureFinalized.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <cassert>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specification of the material params for the three-phase van
 *        Genuchten capillary pressure model.
 *
 * In comparison to the two-phase version, this parameter object also
 * includes the residual saturations, as their handling is very
 * model-specific.
 */
template<class TraitsT>
class ThreePhaseParkerVanGenuchtenParams : public EnsureFinalized
{
public:
    using EnsureFinalized :: finalize;

    typedef TraitsT Traits;
    typedef typename Traits::Scalar Scalar;

    ThreePhaseParkerVanGenuchtenParams()
    {
        betaNW_ = 1.0;
        betaGN_ = 1.0;
    }

    /*!
     * \brief Return the \f$\alpha\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgAlpha() const
    { EnsureFinalized::check(); return vgAlpha_; }

    /*!
     * \brief Set the \f$\alpha\f$ shape parameter of van Genuchten's
     *        curve.
     */
    void setVgAlpha(Scalar v)
    { vgAlpha_ = v; }

    /*!
     * \brief Return the \f$m\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgM() const
    { EnsureFinalized::check(); return vgM_; }

    /*!
     * \brief Set the \f$m\f$ shape parameter of van Genuchten's
     *        curve.
     *
     * The \f$n\f$ shape parameter is set to \f$n = \frac{1}{1 - m}\f$
     */
    void setVgM(Scalar m)
    { vgM_ = m; vgN_ = 1/(1 - vgM_); }

    /*!
     * \brief Return the \f$n\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgN() const
    { EnsureFinalized::check(); return vgN_; }

    /*!
     * \brief Set the \f$n\f$ shape parameter of van Genuchten's
     *        curve.
     *
     * The \f$n\f$ shape parameter is set to \f$m = 1 - \frac{1}{n}\f$
     */
    void setVgN(Scalar n)
    { vgN_ = n; vgM_ = 1 - 1/vgN_; }

    /*!
     * \brief Return the residual wetting saturation.
     */
    Scalar Swr() const
    { EnsureFinalized::check(); return Swr_; }

    /*!
     * \brief Set the residual wetting saturation.
     */
    void setSwr(Scalar input)
    { Swr_ = input; }

    /*!
     * \brief Return the residual non-wetting saturation.
     */
    Scalar Snr() const
    { EnsureFinalized::check(); return Snr_; }

    /*!
     * \brief Set the residual non-wetting saturation.
     */
    void setSnr(Scalar input)
    { Snr_ = input; }

    /*!
     * \brief Return the residual gas saturation.
     */
    Scalar Sgr() const
    { EnsureFinalized::check(); return Sgr_; }

    /*!
     * \brief Set the residual gas saturation.
     */
    void setSgr(Scalar input)
    { Sgr_ = input; }

    Scalar Swrx() const
    { EnsureFinalized::check(); return Swrx_; }

    /*!
     * \brief Set the residual gas saturation.
     */
    void setSwrx(Scalar input)
    { Swrx_ = input; }

    /*!
     * \brief defines the scaling parameters of capillary pressure between the phases (=1 for Gas-Water)
     */
    void setBetaNW(Scalar input)
    { betaNW_ = input; }

    void setBetaGN(Scalar input)
    { betaGN_ = input; }

    /*!
     * \brief Return the values for the beta scaling parameters of capillary pressure between the phases
     */
    Scalar betaNW() const
    { EnsureFinalized::check(); return betaNW_; }

    Scalar betaGN() const
    { EnsureFinalized::check(); return betaGN_; }

    /*!
     * \brief defines if residual n-phase saturation should be regarded in its relative permeability.
     */
    void setkrRegardsSnr(bool input)
    { krRegardsSnr_ = input; }
    /*!
     * \brief Calls if residual n-phase saturation should be regarded in its relative permeability.
     */
    bool krRegardsSnr() const
    { EnsureFinalized::check(); return krRegardsSnr_; }

    void checkDefined() const
    {
        Valgrind::CheckDefined(vgAlpha_);
        Valgrind::CheckDefined(vgM_);
        Valgrind::CheckDefined(vgN_);
        Valgrind::CheckDefined(Swr_);
        Valgrind::CheckDefined(Snr_);
        Valgrind::CheckDefined(Sgr_);
        Valgrind::CheckDefined(Swrx_);
        Valgrind::CheckDefined(betaNW_);
        Valgrind::CheckDefined(betaGN_);
        Valgrind::CheckDefined(krRegardsSnr_);
    }

private:
    Scalar vgAlpha_;
    Scalar vgM_;
    Scalar vgN_;
    Scalar Swr_;
    Scalar Snr_;
    Scalar Sgr_;
    Scalar Swrx_; // Swr + Snr

    Scalar betaNW_;
    Scalar betaGN_;

    bool krRegardsSnr_ ;
};
} // namespace Opm

#endif
