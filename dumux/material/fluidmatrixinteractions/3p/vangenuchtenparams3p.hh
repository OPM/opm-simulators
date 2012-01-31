// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Jochen Fritz, Benjamin Faigle                     *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
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
 * \file vangenuchtenparams.hh Specification of the material params
 *       for the van Genuchten capillary pressure model.
 *
 * In comparison to the 2p version, this parameter container also includes
 * the residual saturations, as their inclusion is very model-specific.
 */
#ifndef VAN_GENUCHTEN_PARAMS_3P_HH
#define VAN_GENUCHTEN_PARAMS_3P_HH

namespace Dumux
{
/*!
 * \brief Reference implementation of a van Genuchten params
 */
template<class ScalarT>
class VanGenuchtenParams3p
{
public:
    typedef ScalarT Scalar;

    VanGenuchtenParams3p()
    {}

    VanGenuchtenParams3p(Scalar vgAlpha, Scalar vgN, Dune::FieldVector<Scalar, 3> residualSaturation, bool regardSnr=false)
    {
        setVgAlpha(vgAlpha);
        setVgN(vgN);
        setSwr(residualSaturation[0]);
        setSnr(residualSaturation[1]);
        setSgr(residualSaturation[2]);
        setkrRegardsSnr(regardSnr);
    };

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

private:
    Scalar vgAlpha_;
    Scalar vgM_;
    Scalar vgN_;
    Scalar Swr_;
    Scalar Snr_;
    Scalar Sgr_;

    bool krRegardsSnr_ ;
};
} // namespace Dumux

#endif
