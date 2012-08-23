// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Holger Class                                      *
 *   Copyright (C) 2012 by Bernd Flemisch                                    *
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
 * \ingroup fluidmatrixinteractionsparams
 *
 * \brief Specification of the material params for the three-phase van
 *        Genuchten capillary pressure model.
 */
#ifndef DUMUX_3P_PARKER_VAN_GENUCHTEN_PARAMS_HH
#define DUMUX_3P_PARKER_VAN_GENUCHTEN_PARAMS_HH

#include <dune/common/fvector.hh>

#include <dumux/common/valgrind.hh>

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionsparams
 *
 * \brief Specification of the material params for the three-phase van
 *        Genuchten capillary pressure model.
 *
 * In comparison to the two-phase version, this parameter object also
 * includes the residual saturations, as their handling is very
 * model-specific.
 */
template<class ScalarT>
class ParkerVanGen3PParams
{
public:
    typedef ScalarT Scalar;

    ParkerVanGen3PParams()
    {betaGW_ = betaNW_ = betaGN_ = 1.;}

    ParkerVanGen3PParams(Scalar vgAlpha, Scalar vgN, Scalar KdNAPL, Scalar rhoBulk, Dune::FieldVector<Scalar, 4> residualSaturation, Scalar betaNW = 1., Scalar betaGN = 1., Scalar betaGW = 1., bool regardSnr=false)
    {
        setVgAlpha(vgAlpha);
        setVgN(vgN);
        setSwr(residualSaturation[0]);
        setSnr(residualSaturation[1]);
        setSgr(residualSaturation[2]);
        setSwrx(residualSaturation[3]);
        setkrRegardsSnr(regardSnr);
        setKdNAPL(KdNAPL);
        setBetaNW(betaNW);
        setBetaGN(betaGN);
        setBetaGW(betaGW);
        setRhoBulk(rhoBulk);
    }

    /*!
     * \brief Return the \f$\alpha\f$ shape parameter of van Genuchten's
     *        curve.
     */
    Scalar vgAlpha() const
    { return vgAlpha_; }

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
    { return vgM_; }

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
    { return vgN_; }

    /*!
     * \brief Set the \f$n\f$ shape parameter of van Genuchten's
     *        curve.
     *
     * The \f$n\f$ shape parameter is set to \f$m = 1 - \frac{1}{n}\f$
     */
    void setVgN(Scalar n)
    { vgN_ = n; vgM_ = 1 - 1/vgN_; }

    /*!
     * \brief Return the residual saturation.
     */
    Scalar satResidual(int phaseIdx) const
    {
        switch (phaseIdx)
        {
        case 0:
            return Swr_;
            break;
        case 1:
            return Snr_;
            break;
        case 2:
            return Sgr_;
            break;
        };
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Set all residual saturations.
     */
    void setResiduals(Dune::FieldVector<Scalar, 3> residualSaturation)
    {
        setSwr(residualSaturation[0]);
        setSnr(residualSaturation[1]);
        setSgr(residualSaturation[2]);
    }


    /*!
     * \brief Return the residual wetting saturation.
     */
    Scalar Swr() const
    { return Swr_; }

    /*!
     * \brief Set the residual wetting saturation.
     */
    void setSwr(Scalar input)
    { Swr_ = input; }

    /*!
     * \brief Return the residual non-wetting saturation.
     */
    Scalar Snr() const
    { return Snr_; }

    /*!
     * \brief Set the residual non-wetting saturation.
     */
    void setSnr(Scalar input)
    { Snr_ = input; }

    /*!
     * \brief Return the residual gas saturation.
     */
    Scalar Sgr() const
    { return Sgr_; }

    /*!
     * \brief Set the residual gas saturation.
     */
    void setSgr(Scalar input)
    { Sgr_ = input; }

    Scalar Swrx() const
    { return Swrx_; }

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

    void setBetaGW(Scalar input)
    { betaGW_ = input; }

    /*!
     * \brief Return the values for the beta scaling parameters of capillary pressure between the phases
     */
    Scalar betaNW() const
    { return betaNW_; }

    Scalar betaGN() const
    { return betaGN_; }

    Scalar betaGW() const
    { return betaGW_; }

    /*!
     * \brief defines if residual n-phase saturation should be regarded in its relative permeability.
     */
    void setkrRegardsSnr(bool input)
    { krRegardsSnr_ = input; }
    /*!
     * \brief Calls if residual n-phase saturation should be regarded in its relative permeability.
     */
    bool krRegardsSnr() const
    { return krRegardsSnr_; }


    /*!
     * \brief Return the bulk density of the porous medium
     */
    Scalar rhoBulk() const
    { return rhoBulk_; }

    /*!
     * \brief Set the bulk density of the porous medium
     */
    void setRhoBulk(Scalar input)
    { rhoBulk_ = input; }

    /*!
     * \brief Return the adsorption coefficient
     */
    Scalar KdNAPL() const
    { return KdNAPL_; }

    /*!
     * \brief Set the adsorption coefficient
     */
    void setKdNAPL(Scalar input)
    { KdNAPL_ = input; }

    void checkDefined() const
    {
        Valgrind::CheckDefined(vgAlpha_);
        Valgrind::CheckDefined(vgM_);
        Valgrind::CheckDefined(vgN_);
        Valgrind::CheckDefined(Swr_);
        Valgrind::CheckDefined(Snr_);
        Valgrind::CheckDefined(Sgr_);
        Valgrind::CheckDefined(Swrx_);
        Valgrind::CheckDefined(KdNAPL_);
        Valgrind::CheckDefined(rhoBulk_);
        Valgrind::CheckDefined(krRegardsSnr_);
    }

private:
    Scalar vgAlpha_;
    Scalar vgM_;
    Scalar vgN_;
    Scalar Swr_;
    Scalar Snr_;
    Scalar Sgr_;
    Scalar Swrx_;     /* (Sw+Sn)_r */

    Scalar KdNAPL_;
    Scalar rhoBulk_;

    Scalar betaNW_;
    Scalar betaGN_;
    Scalar betaGW_;

    bool krRegardsSnr_ ;
};
} // namespace Dumux

#endif
