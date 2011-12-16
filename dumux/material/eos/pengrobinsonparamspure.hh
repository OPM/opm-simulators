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
 * \brief The Peng-Robinson parameters for a pure component
 *
 * See:
 *
 * R. Reid, et al.: The Properties of Gases and Liquids, 4th edition,
 * McGraw-Hill, 1987, pp. 43-44
 */
#ifndef DUMUX_PENG_ROBINSON_PARAMS_PURE_HH
#define DUMUX_PENG_ROBINSON_PARAMS_PURE_HH

#include "pengrobinsonparams.hh"

#include <dumux/material/constants.hh>

namespace Dumux
{
/*!
 * \brief Stores, provides access to and calculates the Peng-Robinson
 *        parameters of a pure component.
 *
 * See:
 *
 * R. Reid, et al.: The Properties of Gases and Liquids, 4th edition,
 * McGraw-Hill, 1987, pp. 43-44
 */
template <class Scalar, class ComponentT>
class PengRobinsonParamsPure : public PengRobinsonParams<Scalar>
{
    typedef PengRobinsonParams<Scalar> ParentType;

    // the ideal gas constant
    static const Scalar R = Dumux::Constants<Scalar>::R;

public:
    typedef ComponentT Component;

    /*!
     * \brief Calculates the "a" and "b" Peng-Robinson parameters for
     *        the component.
     *
     * See:
     *
     * R. Reid, et al.: The Properties of Gases and Liquids, 4th edition,
     * McGraw-Hill, 1987, pp. 43-44
     */
    void update(Scalar T, Scalar p)
    {
        temperature_ = T;
        pressure_ = p;

        Scalar pc = Component::criticalPressure();
        Scalar omega = Component::acentricFactor();
        Scalar Tr = T/Component::criticalTemperature();
        Scalar RTc = R*Component::criticalTemperature();
        Scalar f_omega = 0.37464 + omega*(1.54226 - omega*0.26992);
        Scalar tmp = 1.0 + f_omega*(1.0 - std::sqrt(Tr));
        tmp = tmp*tmp;

        Scalar a = 0.4572355 * RTc*RTc/pc * tmp;
        Scalar b =  0.0777961 * RTc/pc;

        this->setA(a);
        this->setB(b);
    }

    /*!
     * \brief Sets the molar volume [m^3/mol] of the substance.
     *
     * The phaseIdx parameter is there to adhere to the common
     * interface with the multi-phase stuff and is just ignored.
     */
    void setMolarVolume(int phaseIdx, Scalar Vm)
    { setMolarVolume(Vm); }

    /*!
     * \brief Sets the molar volume [m^3/mol] of the substance.
     */
    void setMolarVolume(Scalar Vm)
    { molarVolume_ = Vm; }

    /*!
     * \brief Returns the temperature [K] of the system.
     *
     * The phaseIdx parameter is there to adhere to the common
     * interface with the multi-phase stuff and is just ignored.
     */
    Scalar temperature(int phaseIdx = 0) const
    { return temperature_; }

    /*!
     * \brief Returns the pressure [Pa] of the system.
     *
     * The phaseIdx parameter is there to adhere to the common
     * interface with the multi-phase stuff and is just ignored.
     */
    Scalar pressure(int phaseIdx = 0) const
    { return pressure_; }

    /*!
     * \brief Returns the molar volume [m^3/mol] of the substance.
     *
     * The phaseIdx parameter is there to adhere to the common
     * interface with the multi-phase stuff and is just ignored.
     */
    Scalar molarVolume(int phaseIdx = 0) const
    { return molarVolume_; }

    /*!
     * \brief Returns the attractive parameter "a" [Pa (m^3/mol)^2] for the cubic EOS.
     *
     * The phaseIdx parameter is there to adhere to the common
     * interface with the multi-phase stuff and is just ignored.
     */
    Scalar a(int phaseIdx = 0) const
    { return ParentType::a(); }

    /*!
     * \brief Returns the covolume of the substance [m^3/mol]
     *
     * The phaseIdx parameter is there to adhere to the common
     * interface with the multi-phase stuff and is just ignored.
     */
    Scalar b(int phaseIdx = 0) const
    { return ParentType::b(); }


private:
    Scalar temperature_;
    Scalar pressure_;
    Scalar molarVolume_;
};
} // end namepace

#endif
