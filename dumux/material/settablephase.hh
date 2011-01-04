// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser
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
 * \brief A fluid state for a single phase which an be arbitrarily modified.
 */
#ifndef DUMUX_SETTABLE_PHASE_HH
#define DUMUX_SETTABLE_PHASE_HH

#include <dune/common/fmatrix.hh>

namespace Dumux
{
/*!
 * \brief A fluid state for a single phase which an be arbitrarily modified.
 */
template <class Scalar, class FluidSystem>
class SettablePhase
{
    enum { numComponents = FluidSystem::numComponents };

public:
    SettablePhase()
    { }

    /*!
     * \brief Return the mole fraction of a component within a phase.
     */
    Scalar moleFrac(int compIdx) const
    { return moleFrac_[compIdx]; }

    /*!
     * \brief Return the mass fraction of a component within a phase.
     */
    Scalar massFrac(int compIdx) const
    { return massFrac_[compIdx]; }

    /*!
     * \brief Return the sum of the concentrations of all components
     *        for a phase.
     *
     * Unit: \f$\mathrm{[mol/m^3]}\f$
     */
    Scalar phaseConcentration() const
    { return density_/meanMolarMass_; }

    /*!
     * \brief Return the concentrations of an individual component in
     *        a phase.
     *
     * Unit: \f$\mathrm{[mol/m^3]}\f$
     */
    Scalar concentration(int compIdx) const
    { return moleFrac_[compIdx]*phaseConcentration(); }

    /*!
     * \brief Return the density of a phase.
     *
     * Unit: \f$\mathrm{[kg/m^3]}\f$
     */
    Scalar density() const
    { return density_; }

    /*!
     * \brief Return the average molar mass of a phase.
     *
     * This is the sum of all molar masses times their respective mole
     * fractions in the phase.
     *
     * Unit: \f$\mathrm{[kg/m^3]}\f$
     */
    Scalar meanMolarMass() const
    { return meanMolarMass_; }

    /*!
     * \brief Return the total pressure of the phase.
     *
     * Unit: \f$\mathrm{[Pa] = [N/m^2]}\f$
     */
    Scalar pressure() const
    { return pressure_; };

    /*!
     * \brief Calculates the mole fractions from the mass fractions.
     *
     * Also updates the mean molar mass.
     */
    void XTox()
    {
        typedef Dune::FieldMatrix<Scalar, numComponents, numComponents> Matrix;
        typedef Dune::FieldVector<Scalar, numComponents> Vector;

        Matrix M;
        Vector b;

        Scalar sumX = 0;
        for (int j = 0; j < numComponents; ++j)
            sumX += massFrac_[j];

        // Calculate the linear system of equations which determines
        // the conversion from mass fractions to mole fractions
        for (int i = 0; i < numComponents - 1; ++i) {
            for (int j = 0; j < numComponents; ++j) {
                M[i][j] = FluidSystem::molarMass(j)*massFrac_[i];
            }

            // main diagonal is different from the remaining matrix
            M[i][i] -= FluidSystem::molarMass(i)*sumX;

            // right hand side
            b[i] = 0;
        }

        // set up last row
        for (int j = 0; j < numComponents; ++j)
            M[numComponents-1][j] = 1;
        b[numComponents-1] = sumX;

        // solve the system of equations
        M.solve(moleFrac, massFrac);

        // calculate mean molar mass
        meanMolarMass_ = 0;
        for (int i = 0; i < numComponents; ++i) {
            meanMolarMass_ += moleFrac_[i] * FluidSystem::molarMass(i);
        }
        meanMolarMass_ /= sumX;
    }

    /*!
     * \brief Calculates the mass fractions from the mole fractions.
     *
     * Also updates the mean molar mass.
     */
    void xToX()
    {
        // calculate mean molar mass
        meanMolarMass_ = 0;
        Scalar sumx = 0;
        for (int i = 0; i < numComponents; ++i) {
            sumx += moleFrac_[i];
            meanMolarMass_ += moleFrac_[i] * FluidSystem::molarMass(i);
        }
        meanMolarMass_ /= sumx;

        // calculate mass fractions
        for (int i = 0; i < numComponents; ++i)
            massFrac_[i] = moleFrac_[i]*FluidSystem::molarMass(i)/meanMolarMass_;
    }

    Scalar pressure_;
    Scalar density_;
    Scalar meanMolarMass_;
    Scalar moleFrac_[numComponents];
    Scalar massFrac_[numComponents];
};

} // namespace Dumux

#endif
