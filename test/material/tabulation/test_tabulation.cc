// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
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
 * \brief This is a program to test the tabulation class for of
 *        individual components.
 *
 * It either prints "success" or "error", it does not do anything
 * else.
 */
#include "config.h"

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

template <class Scalar>
void isSame(Scalar v, Scalar vRef, Scalar tol=5e-4)
{
    if (std::abs( (v - vRef)/vRef ) > tol) {
        std::cout << "\nerror: " << (v - vRef)/vRef << "\n";
        exit(1);
    }
}

int main()
{
    typedef double Scalar;
    typedef Dumux::H2O<Scalar> IapwsH2O;
    typedef Dumux::TabulatedComponent<Scalar, IapwsH2O> TabulatedH2O;

    Scalar tempMin = 274.15;
    Scalar tempMax = 622.15;
    int nTemp = (int) (tempMax - tempMin)*3/2;

    Scalar pMin = 10.00;
    Scalar pMax = IapwsH2O::vaporPressure(tempMax*1.1);
    int nPress = 600;
    
    std::cout << "Creating tabulation with " << nTemp*nPress << " entries per quantity\n";
    TabulatedH2O::init(tempMin, tempMax, nTemp,
                       pMin, pMax, nPress);
    
    std::cout << "Checking tabulation\n";

    int m = nTemp*3;
    int n = nPress*3;   
    for (int i = 0; i < m; ++i) {
        Scalar T = tempMin + (tempMax - tempMin)*Scalar(i)/m;

        if (i % std::max(1, m/1000) == 0) {
            std::cout << Scalar(i)/m*100 << "% done        \r";
            std::cout.flush();
        }
        
        isSame(TabulatedH2O::vaporPressure(T), 
               IapwsH2O::vaporPressure(T), 
               1e-3);
        for (int j = 0; j < n; ++j) {
            Scalar p = pMin + (pMax - pMin)*Scalar(j)/n;
            if (p < IapwsH2O::vaporPressure(T) * 1.03) {
                Scalar tol = 5e-4;
                if (p > IapwsH2O::vaporPressure(T))
                    tol = 5e-2;
                isSame(TabulatedH2O::gasEnthalpy(T,p), IapwsH2O::gasEnthalpy(T,p), tol);
                isSame(TabulatedH2O::gasInternalEnergy(T,p), IapwsH2O::gasInternalEnergy(T,p), tol);
                isSame(TabulatedH2O::gasDensity(T,p), IapwsH2O::gasDensity(T,p), tol);
                isSame(TabulatedH2O::gasViscosity(T,p), IapwsH2O::gasViscosity(T,p), tol);
            }
            
            if (p > IapwsH2O::vaporPressure(T) / 1.03) {
                Scalar tol = 5e-4;
                if (p < IapwsH2O::vaporPressure(T))
                    tol = 5e-2;
                isSame(TabulatedH2O::liquidEnthalpy(T,p), IapwsH2O::liquidEnthalpy(T,p), tol);
                isSame(TabulatedH2O::liquidInternalEnergy(T,p), IapwsH2O::liquidInternalEnergy(T,p), tol);
                isSame(TabulatedH2O::liquidDensity(T,p), IapwsH2O::liquidDensity(T,p), tol);
                isSame(TabulatedH2O::liquidViscosity(T,p), IapwsH2O::liquidViscosity(T,p), tol);
            }
        }
    }
    std::cout << "\nsuccess\n";
    return 0;
}
