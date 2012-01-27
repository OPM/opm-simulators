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
 * \file linearmaterial.hh Implements a linear saturation-capillary
 *                    pressure relation
 */
#ifndef LINEAR_MATERIAL_3P_HH
#define LINEAR_MATERIAL_3P_HH

#include "linearmaterialparams3p.hh"

#include <algorithm>

namespace Dumux
{
/*!
 * \ingroup material
 *
 * \brief Implements a linear saturation-capillary pressure relation
 *
 *
 * \sa LinearMaterialParams3p
 */
template <class ScalarT, class ParamsT = LinearMaterialParams3P<ScalarT> >
class LinearMaterial3P
{
public:
    typedef ParamsT Params;
    typedef typename Params::Scalar Scalar;

    /*!
     * \brief The linear capillary pressure-saturation curve.
     *
     * \param Swe Effective saturation of of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar pC(const Params &params, Scalar Swe)
    {
        DUNE_THROW(Dune::NotImplemented, "Capillary pressures for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief The saturation-capillary pressure curve.
     *
     *
     * \param pC Capillary pressure \f$\p_C\f$
     * \return The effective saturaion of the wetting phase \f$\overline{S}_w\f$
     */
    static Scalar Sw(const Params &params, Scalar pC)
    {
        DUNE_THROW(Dune::NotImplemented, "Capillary pressures for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief Returns the partial derivative of the capillary
     *        pressure to the effective saturation.
     *
    */
    static Scalar dpC_dSw(const Params &params, Scalar Swe)
    {
        DUNE_THROW(Dune::NotImplemented, "Capillary pressures for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief Returns the partial derivative of the effective
     *        saturation to the capillary pressure.
     */
    static Scalar dSw_dpC(const Params &params, Scalar pC)
    {
        DUNE_THROW(Dune::NotImplemented, "Capillary pressures for three phases not implemented! Do it yourself!");
    }

    /*!
     * \brief The relative permeability for a phase.
     *
     * \param Phase indicator, The saturation of all phases.
     */
    static Scalar kr(const int phase, const Params &params, const Dune::FieldVector<Scalar, 3> saturation)
    {
        switch (phase)
        {
        case 0:
            return krw(params, saturation);
            break;
        case 1:
            return krn(params, saturation);
            break;
        case 2:
            return krg(params, saturation);
            break;
        }
        return 0;
    };

    /*!
     * \brief The relative permeability for the wetting phase.
     *
     * \param The saturation of all phases.
     */
    static Scalar krw(const Params &params, const Dune::FieldVector<Scalar, 3> saturation)
    {
        Scalar krw_ = std::max(std::min((saturation[0] - params.Swr())/(1- params.Swr()), 1.), 0.);
        return krw_;
    };

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param The saturation of all phases.
     */
    static Scalar krn(const Params &params, const Dune::FieldVector<Scalar, 3> saturation)
    {
        Scalar krn_ = std::max(std::min((saturation[1] - params.Snr())/(1- params.Snr()), 1.), 0.);
        return krn_;
    }

    /*!
     * \brief The relative permeability for the non-wetting phase.
     *
     * \param The saturation of all phases.
     */
    static Scalar krg(const Params &params, const Dune::FieldVector<Scalar, 3> saturation)
    {
        Scalar krg_ = std::max(std::min((saturation[2] - params.Sgr())/(1- params.Sgr()), 1.), 0.);
        return krg_;
    }
};
}

#endif
