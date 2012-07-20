// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
 * \file
 * \ingroup fluidmatrixinteractionsparams
 *
 * \brief Reference implementation of the parameter object for the
 *        linear material law.
 */
#ifndef LINEAR_MATERIAL_PARAMS_HH
#define LINEAR_MATERIAL_PARAMS_HH

namespace Dumux
{
/*!
 * \ingroup fluidmatrixinteractionsparams
 *
 * \brief Reference implementation of the parameter object for the
 *        linear material law.
 */
template<class ScalarT>
class LinearMaterialParams
{
public:
    typedef ScalarT Scalar;

    LinearMaterialParams()
    {
        setEntryPC(0);
        setMaxPC(0);
    }

    LinearMaterialParams(Scalar entryPC, Scalar maxPC)
    {
        setEntryPC(entryPC);
        setMaxPC(maxPC);
    }


    /*!
     * \brief Return the entry pressure for the linear material law.
     *
     * The entry pressure is reached at \f$\overline S_w = 1\f$
     */
    Scalar entryPC() const
    { return entryPC_; }

    /*!
     * \brief Set the entry pressure for the linear material law.
     *
     * The entry pressure is reached at \f$ \overline S_w = 1\f$
     */
    void setEntryPC(Scalar v)
    { entryPC_ = v; }

    /*!
     * \brief Return the maximum capillary pressure for the linear material law.
     *
     * The maximum capillary pressure is reached at \f$ \overline S_w = 0\f$
     */
    Scalar maxPC() const
    { return maxPC_; }

    /*!
     * \brief Set the maximum capillary pressure for the linear material law.
     *
     * The maximum capillary pressure is reached at \f$ \overline S_w = 0\f$
     */
    void setMaxPC(Scalar v)
    { maxPC_ = v; }


private:
    Scalar entryPC_;
    Scalar maxPC_;
};
} // namespace Dumux

#endif
