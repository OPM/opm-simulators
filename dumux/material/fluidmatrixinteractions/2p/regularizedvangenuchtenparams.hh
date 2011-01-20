// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Andreas Lauser, Bernd Flemisch                    *
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
 * \file
 *
 * \brief   Parameters that are necessary for the \em regularization of
 *          VanGenuchten "material law".
 */

#ifndef REGULARIZED_VAN_GENUCHTEN_PARAMS_HH
#define REGULARIZED_VAN_GENUCHTEN_PARAMS_HH

#include "vangenuchtenparams.hh"

namespace Dumux
{
/*!
 *
 *
 * \brief   Parameters that are necessary for the \em regularization of
 *          VanGenuchten "material law".
 *
 * \ingroup fluidmatrixinteractionsparams
 */
template<class ScalarT>
class RegularizedVanGenuchtenParams : public VanGenuchtenParams<ScalarT>
{
public:
    typedef ScalarT Scalar;
    typedef VanGenuchtenParams<Scalar> Parent;
    typedef RegularizedVanGenuchtenParams<Scalar> Self;

    RegularizedVanGenuchtenParams()
    {}

    RegularizedVanGenuchtenParams(Scalar vgAlpha,
                                   Scalar vgN)
        : Parent(vgAlpha, vgN)
    {};

    /*!
     * \brief Threshold saturation below which the capillary pressure
     *        is regularized.
     *
     * This is just 1%. If you need a different value, overload this
     * class.
     */
    Scalar pCLowSw() const
    {
        // Some problems are very sensitive to this value
        // (e.g. makeing it smaller might result in negative
        // pressures), if you change it here, you will almost
        // certainly break someone's code!
        //
        // If you want to use a different regularization threshold,
        // overload this class and supply the new class as second
        // template parameter for the RegularizedVanGenuchten law!
        return /* PLEASE DO _NOT_ */ 1e-2; /* CHANGE THIS VALUE. READ
                                            * COMMENT ABOVE! */
   }

    /*!
     * \brief Threshold saturation above which the capillary pressure
     *        is regularized.
     *
     * This is just 99%. If you need a different value, overload this
     * class.
     */
    Scalar pCHighSw() const
    {
        // Some problems are very sensitive to this value
        // (e.g. makeing it smaller might result in negative
        // pressures), if you change it here, you will almost
        // certainly break someone's code!
        //
        // If you want to use a different regularization threshold,
        // overload this class and supply the new class as second
        // template parameter for the RegularizedVanGenuchten law!
        return /* PLEASE DO _NOT_ */ 99e-2; /* CHANGE THIS VALUE. READ
                                             * COMMENT ABOVE! */
    }

    /*!
     * \brief Threshold saturation below which the relative
     *        permeability of the non-wetting phase gets regulatized.
     *
     * This is just 10%. If you need a different value, overload this
     * class.
     */
    Scalar krnLowSw() const
    { return 0.10; }

    /*!
     * \brief Threshold saturation above which the relative
     *        permeability of the wetting phase gets regulatized.
     *
     * This is just 90%. If you need a different value, overload this
     * class.
     */
    Scalar krwHighSw() const
    { return 0.90; }

};
}; // namespace Dumux

#endif
