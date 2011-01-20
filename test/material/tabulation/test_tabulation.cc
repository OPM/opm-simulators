// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
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
 * \brief This is a program to test the tabulation class for of
 *        individual components.
 *
 * It either prints "success" or "error", it does not do anything
 * else.
 */
#include "config.h"

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

bool success;

template <class Scalar>
void isSame(const char *str, Scalar v, Scalar vRef, Scalar tol=5e-4)
{
    if (std::abs( (v - vRef)/vRef ) > tol) {
        std::cout << "error for \"" << str << "\": "  << (v - vRef)/vRef*100 << "% difference\n";
        success = false;
        //exit(1);
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
    int nPress = 200;

    std::cout << "Creating tabulation with " << nTemp*nPress << " entries per quantity\n";
    TabulatedH2O::init(tempMin, tempMax, nTemp,
                       pMin, pMax, nPress);

    std::cout << "Checking tabulation\n";
    success = true;
    int m = nTemp*3;
    int n = nPress*3;
    for (int i = 0; i < m; ++i) {
        Scalar T = tempMin + (tempMax - tempMin)*Scalar(i)/m;

        if (i % std::max(1, m/1000) == 0) {
            std::cout << Scalar(i)/m*100 << "% done        \r";
            std::cout.flush();
        }

        isSame("vaporPressure",
               TabulatedH2O::vaporPressure(T),
               IapwsH2O::vaporPressure(T),
               1e-3);
        for (int j = 0; j < n; ++j) {
            Scalar p = pMin + (pMax - pMin)*Scalar(j)/n;
            if (p < IapwsH2O::vaporPressure(T) * 1.03) {
                Scalar tol = 5e-4;
                if (p > IapwsH2O::vaporPressure(T))
                    tol = 5e-2;
                Scalar rho = IapwsH2O::gasDensity(T,p);
                //std::cerr << T << " " << p << " " << IapwsH2O::gasPressure(T,rho) << " " << TabulatedH2O::gasPressure(T,rho) << "\n";
                isSame("Iapws::gasPressure", IapwsH2O::gasPressure(T,rho), p, 1e-6);
                isSame("gasPressure", TabulatedH2O::gasPressure(T,rho), p, 2e-2);
                isSame("gasEnthalpy", TabulatedH2O::gasEnthalpy(T,p), IapwsH2O::gasEnthalpy(T,p), tol);
                isSame("gasInternalEnergy", TabulatedH2O::gasInternalEnergy(T,p), IapwsH2O::gasInternalEnergy(T,p), tol);
                isSame("gasDensity", TabulatedH2O::gasDensity(T,p), rho, tol);
                isSame("gasViscosity", TabulatedH2O::gasViscosity(T,p), IapwsH2O::gasViscosity(T,p), tol);
            }

            if (p > IapwsH2O::vaporPressure(T) / 1.03) {
                Scalar tol = 5e-4;
                if (p < IapwsH2O::vaporPressure(T))
                    tol = 5e-2;
                Scalar rho = IapwsH2O::liquidDensity(T,p);
                //std::cerr << T << " " << p << " " << IapwsH2O::liquidPressure(T,rho) << " " << TabulatedH2O::liquidPressure(T,rho) << "\n";
                isSame("Iapws::gasPressure", IapwsH2O::liquidPressure(T,rho), p, 1e-6);
                isSame("liquidPressure", TabulatedH2O::liquidPressure(T,rho), p, 2e-2);
                isSame("liquidEnthalpy", TabulatedH2O::liquidEnthalpy(T,p), IapwsH2O::liquidEnthalpy(T,p), tol);
                isSame("liquidInternalEnergy", TabulatedH2O::liquidInternalEnergy(T,p), IapwsH2O::liquidInternalEnergy(T,p), tol);
                isSame("liquidDensity", TabulatedH2O::liquidDensity(T,p), rho, tol);
                isSame("liquidViscosity", TabulatedH2O::liquidViscosity(T,p), IapwsH2O::liquidViscosity(T,p), tol);
            }
        }
        //std::cerr << "\n";
    }

    if (success)
        std::cout << "\nsuccess\n";
    return 0;
}
