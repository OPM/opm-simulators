/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file 
 *
 * \brief Binary coefficients for nitrogen and oxygen.
 */
#ifndef DUMUX_BINARY_COEFF_N2_O2_HH
#define DUMUX_BINARY_COEFF_N2_O2_HH

#include "henryiapws.hh"
#include "fullermethod.hh"

#include <dumux/new_material/components/o2.hh>
#include <dumux/new_material/components/n2.hh>

namespace Dumux
{
namespace BinaryCoeff
{

/*!
 * \brief Binary coefficients for nitrogen and oxygen.
 */
class N2_O2
{
public:
    /*!
     * \brief Henry coefficent \f$[N/m^2]\f$  for molecular oxygen in liquid nitrogen.
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        DUNE_THROW(Dune::NotImplemented, "henry coefficient for oxygen in liquid nitrogen");
    };

    /*!
     * \brief Binary diffusion coefficent [m^2/s] for molecular nitrogen and oxygen.
     *
     * \copybody fullerMethod()
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        typedef Dumux::N2<Scalar> N2;
        typedef Dumux::O2<Scalar> O2;
        
        // atomic diffusion volumes
        const Scalar SigmaNu[2] = { 18.1 /* N2 */,  16.3 /* O2 */ };
        // molar masses [g/mol]
        const Scalar M[2] = { N2::molarMass()*1e3, O2::molarMass()*1e3 };
        return fullerMethod(M, SigmaNu, temperature, pressure);
    };

    /*!
     * \brief Diffusion coefficent [m^2/s] for molecular oxygen in liquid nitrogen.
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "diffusion coefficient for liquid oxygen and nitrogen");
    };
};

}
} // end namepace

#endif
