/*****************************************************************************
 *   Copyright (C) 2010 by Jochen Fritz, Benjamin Faigle                     *
 *   Institute of Hydraulic Engineering                                      *
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
 * \file linearmaterialparams.hh Reference implementation of params for
 *                           the linear material material.
 */
#ifndef LINEAR_MATERIAL_PARAMS_3P_HH
#define LINEAR_MATERIAL_PARAMS_3P_HH

namespace Dumux
{
/*!
 * \brief Reference implementation of params for the linear material
 *        material.
 */
template<class ScalarT>
class LinearMaterialParams3P
{
public:
    typedef ScalarT Scalar;

    LinearMaterialParams3P()
    {}

    LinearMaterialParams3P(Scalar Swr, Scalar Snr, Scalar Sgr)
    {
        setSwr(Swr);
        setSnr(Snr);
        setSgr(Sgr);
    };

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

private:
    Scalar Swr_;
    Scalar Snr_;
    Scalar Sgr_;
};
} // namespace Dumux

#endif
