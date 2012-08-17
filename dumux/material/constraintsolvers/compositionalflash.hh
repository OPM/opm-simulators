// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Benjamin Faigle                                   *
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
 *
 * \brief Determines the pressures and saturations of all fluid phases
 *        given the total mass of all components.
 */
#ifndef DUMUX_COMPOSITIONAL_FLASH_HH
#define DUMUX_COMPOSITIONAL_FLASH_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/decoupled/2p2c/2p2cproperties.hh>
#include <dumux/decoupled/2p2c/2p2cfluidstate.hh>
#include <dumux/decoupled/2p2c/pseudo1p2cfluidstate.hh>

namespace Dumux
{
/*!
 * \brief Flash calculation routines for compositional decoupled models
 *
 *        Routines for isothermal and isobaric 2p2c and 1p2c flash.
 *  \tparam TypeTag The property Type Tag
 */
template <class TypeTag>
class CompositionalFlash
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar)      Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;
    typedef PseudoOnePTwoCFluidState<TypeTag> FluidState1p2c;

    enum {  numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
            numComponents = GET_PROP_VALUE(TypeTag, NumComponents)};

    enum{
        wPhaseIdx = 0,
        nPhaseIdx = 1,
        wCompIdx = 0,
        nCompIdx = 1
    };

public:
    typedef Dune::FieldVector<Scalar, numComponents> ComponentVector;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;
/*!
 * \name Concentration flash for a given feed fraction
 */
//@{
    //! 2p2c Flash for constant p & t if concentrations (feed mass fraction) is given.
    /*!
     * Routine goes as follows:
     * - determination of the equilibrium constants from the fluid system
     * - determination of maximum solubilities (mole fractions) according to phase pressures
     * - comparison with Z1 to determine phase presence => phase mass fractions
     * - round off fluid properties
     * \param fluidState The decoupled fluid State
     * \param Z1 Feed mass fraction: Mass of comp1 per total mass \f$\mathrm{[-]}\f$
     * \param phasePressure Vector holding the pressure \f$\mathrm{[Pa]}\f$
     * \param poro Porosity \f$\mathrm{[-]}\f$
     * \param temperature Temperature \f$\mathrm{[K]}\f$
     */
    static void concentrationFlash2p2c(FluidState &fluidState,
                             const Scalar &Z1,
                             const PhaseVector &phasePressure,
                             const Scalar &porosity,
                             const Scalar &temperature)
    {
        // set the temperature, pressure
        fluidState.setTemperature(temperature);
        fluidState.setPressure(wPhaseIdx, phasePressure[wPhaseIdx]);
        fluidState.setPressure(nPhaseIdx, phasePressure[nPhaseIdx]);

        //mole equilibrium ratios K for in case wPhase is reference phase
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);
        double k1 = FluidSystem::fugacityCoefficient(fluidState, paramCache, wPhaseIdx, wCompIdx);    // = p^wComp_vap
        double k2 = FluidSystem::fugacityCoefficient(fluidState, paramCache, wPhaseIdx, nCompIdx);    // = H^nComp_w

        // get mole fraction from equilibrium konstants
        fluidState.setMoleFraction(wPhaseIdx,wCompIdx, ((1. - k2) / (k1 -k2)));
        fluidState.setMoleFraction(nPhaseIdx,wCompIdx, (fluidState.moleFraction(wPhaseIdx,wCompIdx) * k1));

        // transform mole to mass fractions
        fluidState.setMassFraction(wPhaseIdx, wCompIdx,
                (fluidState.moleFraction(wPhaseIdx,wCompIdx) * FluidSystem::molarMass(wCompIdx)
                / ( fluidState.moleFraction(wPhaseIdx,wCompIdx) * FluidSystem::molarMass(wCompIdx)
                    + (1.-fluidState.moleFraction(wPhaseIdx,wCompIdx)) * FluidSystem::molarMass(nCompIdx) )));
        fluidState.setMassFraction(nPhaseIdx,wCompIdx,
                (fluidState.moleFraction(nPhaseIdx,wCompIdx) * FluidSystem::molarMass(wCompIdx)
                / ( fluidState.moleFraction(nPhaseIdx,wCompIdx) * FluidSystem::molarMass(wCompIdx)
                    + (1.-fluidState.moleFraction(nPhaseIdx,wCompIdx)) * FluidSystem::molarMass(nCompIdx) )));

        //mass equilibrium ratios
        Scalar equilRatio_[numPhases][numComponents];
        equilRatio_[nPhaseIdx][wCompIdx] = fluidState.massFraction(nPhaseIdx,wCompIdx)
                / fluidState.massFraction(wPhaseIdx, wCompIdx);     // = Xn1 / Xw1 = K1
        equilRatio_[nPhaseIdx][nCompIdx] = (1.-fluidState.massFraction(nPhaseIdx, wCompIdx))
                / (1.-fluidState.massFraction(wPhaseIdx, wCompIdx)); // =(1.-Xn1) / (1.-Xw1)     = K2
        equilRatio_[wPhaseIdx][nCompIdx] = equilRatio_[wPhaseIdx][wCompIdx] = 1.;

        // phase fraction of nPhase [mass/totalmass]
        fluidState.setNu(nPhaseIdx, 0.);

        // check if there is enough of component 1 to form a phase
        if (Z1 > fluidState.massFraction(nPhaseIdx,wCompIdx)
                             && Z1 < fluidState.massFraction(wPhaseIdx,wCompIdx))
            fluidState.setNu(nPhaseIdx, -((equilRatio_[nPhaseIdx][wCompIdx]-1)*Z1 + (equilRatio_[nPhaseIdx][nCompIdx]-1)*(1-Z1))
                                / (equilRatio_[nPhaseIdx][wCompIdx]-1) / (equilRatio_[nPhaseIdx][nCompIdx] -1));
        else if (Z1 <= fluidState.massFraction(nPhaseIdx,wCompIdx)) // too little wComp to form a phase
        {
            fluidState.setNu(nPhaseIdx, 1.); // only nPhase
            fluidState.setMassFraction(nPhaseIdx,wCompIdx, Z1); // hence, assign complete mass dissolved into nPhase

            // store as moleFractions
            Scalar xw_n = Z1 /*=Xw_n*/ / FluidSystem::molarMass(wCompIdx);  // = moles of compIdx
            xw_n /= ( Z1 /*=Xw_n*/ / FluidSystem::molarMass(wCompIdx)
                    +(1- Z1 /*=Xn_n*/) / FluidSystem::molarMass(nCompIdx) ); // /= total moles in phase

            fluidState.setMoleFraction(nPhaseIdx,wCompIdx, xw_n);

//            // opposing non-present phase is already set to equilibrium mass fraction
//            fluidState.setMassFraction(wPhaseIdx,wCompIdx, 1.); // non present phase is set to be pure
//            fluidState.setMoleFraction(wPhaseIdx,wCompIdx, 1.); // non present phase is set to be pure
        }
        else    // (Z1 >= Xw1) => no nPhase
        {
            fluidState.setNu(nPhaseIdx, 0.); // no second phase
            fluidState.setMassFraction(wPhaseIdx, wCompIdx, Z1);

            // store as moleFractions
            Scalar xw_w = Z1 /*=Xw_w*/ / FluidSystem::molarMass(wCompIdx);  // = moles of compIdx
            xw_w /= ( Z1 /*=Xw_w*/ / FluidSystem::molarMass(wCompIdx)
                    +(1- Z1 /*=Xn_w*/) / FluidSystem::molarMass(nCompIdx) ); // /= total moles in phase
            fluidState.setMoleFraction(wPhaseIdx, wCompIdx, xw_w);

//            // opposing non-present phase is already set to equilibrium mass fraction
//            fluidState.setMassFraction(nPhaseIdx,wCompIdx, 0.); // non present phase is set to be pure
//            fluidState.setMoleFraction(nPhaseIdx,wCompIdx, 0.); // non present phase is set to be pure
        }

        // complete array of mass fractions
        fluidState.setMassFraction(wPhaseIdx, nCompIdx, 1. - fluidState.massFraction(wPhaseIdx,wCompIdx));
        fluidState.setMassFraction(nPhaseIdx, nCompIdx, 1. - fluidState.massFraction(nPhaseIdx,wCompIdx));
        // complete array of mole fractions
        fluidState.setMoleFraction(wPhaseIdx, nCompIdx, 1. - fluidState.moleFraction(wPhaseIdx,wCompIdx));
        fluidState.setMoleFraction(nPhaseIdx, nCompIdx, 1. - fluidState.moleFraction(nPhaseIdx,wCompIdx));

        // complete phase mass fractions
        fluidState.setNu(wPhaseIdx, 1. - fluidState.phaseMassFraction(nPhaseIdx));

        // get densities with correct composition
        fluidState.setDensity(wPhaseIdx, FluidSystem::density(fluidState, paramCache, wPhaseIdx));
        fluidState.setDensity(nPhaseIdx, FluidSystem::density(fluidState, paramCache, nPhaseIdx));

        Scalar Sw = fluidState.phaseMassFraction(wPhaseIdx) / fluidState.density(wPhaseIdx);
        Sw /= (fluidState.phaseMassFraction(wPhaseIdx)/fluidState.density(wPhaseIdx)
                    + fluidState.phaseMassFraction(nPhaseIdx)/fluidState.density(nPhaseIdx));
        fluidState.setSaturation(wPhaseIdx, Sw);
    };

    //! The simplest possible update routine for 1p2c "flash" calculations
    /*!
     * Routine goes as follows:
     * - Check if we are in single phase condition
     * - Assign total concentration to the present phase
     *
     * \param Z1 Feed mass fraction \f$\mathrm{[-]}\f$
     * \param phasePressure Vector holding the pressure \f$\mathrm{[Pa]}\f$
     * \param presentPhaseIdx Subdomain Index = Indication which phase is present
     * \param temperature Temperature \f$\mathrm{[K]}\f$
     */
    static void concentrationFlash1p2c(FluidState1p2c& fluidState, const Scalar& Z1,const Dune::FieldVector<Scalar, numPhases> phasePressure,const int presentPhaseIdx, const Scalar& temperature)
    {
        // set the temperature, pressure
        fluidState.setTemperature(temperature);
        fluidState.setPressure(wPhaseIdx, phasePressure[wPhaseIdx]);
        fluidState.setPressure(nPhaseIdx, phasePressure[nPhaseIdx]);

        fluidState.setPresentPhaseIdx(presentPhaseIdx);
        fluidState.setMassFraction(presentPhaseIdx,wCompIdx, Z1);

        // calculate mole fraction and average molar mass
        Scalar xw_alpha= Z1 / FluidSystem::molarMass(wCompIdx);
        xw_alpha /= ( Z1 / FluidSystem::molarMass(wCompIdx)
                + (1.-Z1) / FluidSystem::molarMass(nCompIdx));
        fluidState.setMoleFraction(presentPhaseIdx, wCompIdx, xw_alpha);

//        if (presentPhaseIdx == wPhaseIdx)
//        {
//
////            fluidState.setMassFraction(wPhaseIdx,wCompIdx, 0.;
//
//
//
//
//
////            fluidState.moleFractionWater_[nPhaseIdx] = 0.;
//
//            fluidState.setPresentPhaseIdx(presentPhaseIdx);
//        }
//        else if (presentPhaseIdx == nPhaseIdx)
//        {
//            fluidState.setMassFraction[wPhaseIdx] = 0.;
//            fluidState.setMassFraction[nPhaseIdx] = Z1;
//
//            // interested in nComp => 1-X1
//            fluidState.moleFractionWater_[nPhaseIdx] = ( Z1 / FluidSystem::molarMass(0) );   // = moles of compIdx
//            fluidState.moleFractionWater_[nPhaseIdx] /= (Z1/ FluidSystem::molarMass(0)
//                           + (1.-Z1) / FluidSystem::molarMass(1) );    // /= total moles in phase
//            fluidState.moleFractionWater_[nPhaseIdx] = 0.;
//
//            fluidState.presentPhaseIdx_ = nPhaseIdx;
//        }
//        else
//            Dune::dgrave << __FILE__ <<": Twophase conditions in single-phase flash!"
//                << " Z1 is "<<Z1<< std::endl;

        fluidState.setAverageMolarMass(presentPhaseIdx,
                fluidState.massFraction(presentPhaseIdx, wCompIdx)*FluidSystem::molarMass(wCompIdx)
                +fluidState.massFraction(presentPhaseIdx, nCompIdx)*FluidSystem::molarMass(nCompIdx));

        typename FluidSystem::ParameterCache paramCache;
        paramCache.updatePhase(fluidState, presentPhaseIdx);
        fluidState.setDensity(presentPhaseIdx, FluidSystem::density(fluidState, paramCache, presentPhaseIdx));

        return;
    }
//@}

/*!
 * \name Saturation flash for a given saturation (e.g. at boundary)
 */
//@{
    //! a flash routine for 2p2c systems if the saturation instead of total concentration is known.
    /*!
     * Routine goes as follows:
     * - determination of the equilibrium constants from the fluid system
     * - determination of maximum solubilities (mole fractions) according to phase pressures
     * - round off fluid properties
     * \param fluidState The decoupled fluid state
     * \param sat Saturation of phase 1 \f$\mathrm{[-]}\f$
     * \param phasePressure Vector holding the pressure \f$\mathrm{[Pa]}\f$
     * \param poro Porosity \f$\mathrm{[-]}\f$
     * \param temperature Temperature \f$\mathrm{[K]}\f$
     */
    static void saturationFlash2p2c(FluidState &fluidState,
            const Scalar &saturation,
            const PhaseVector &phasePressure,
            const Scalar &porosity,
            const Scalar &temperature)
    {
        if (saturation == 0. || saturation == 1.)
                    Dune::dinfo << "saturation initial and boundary conditions set to zero or one!"
                        << " assuming fully saturated compositional conditions" << std::endl;

        // set the temperature, pressure
        fluidState.setTemperature(temperature);
        fluidState.setPressure(wPhaseIdx, phasePressure[wPhaseIdx]);
        fluidState.setPressure(nPhaseIdx, phasePressure[nPhaseIdx]);

        //in contrast to the standard update() method, satflash() does not calculate nu.
        fluidState.setNu(wPhaseIdx, NAN);
        fluidState.setNu(nPhaseIdx, NAN);

        //mole equilibrium ratios K for in case wPhase is reference phase
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);
        double k1 = FluidSystem::fugacityCoefficient(fluidState, paramCache, wPhaseIdx, wCompIdx);    // = p^wComp_vap
        double k2 = FluidSystem::fugacityCoefficient(fluidState, paramCache, wPhaseIdx, nCompIdx);    // = H^nComp_w

        // get mole fraction from equilibrium konstants
        fluidState.setMoleFraction(wPhaseIdx,wCompIdx, ((1. - k2) / (k1 -k2)));
        fluidState.setMoleFraction(nPhaseIdx,wCompIdx, (fluidState.moleFraction(wPhaseIdx,wCompIdx) * k1));

        // transform mole to mass fractions
        fluidState.setMassFraction(wPhaseIdx, wCompIdx,
                (fluidState.moleFraction(wPhaseIdx,wCompIdx) * FluidSystem::molarMass(wCompIdx)
                / ( fluidState.moleFraction(wPhaseIdx,wCompIdx) * FluidSystem::molarMass(wCompIdx)
                    + (1.-fluidState.moleFraction(wPhaseIdx,wCompIdx)) * FluidSystem::molarMass(nCompIdx) )));
        fluidState.setMassFraction(nPhaseIdx,wCompIdx,
                (fluidState.moleFraction(nPhaseIdx,wCompIdx) * FluidSystem::molarMass(wCompIdx)
                / ( fluidState.moleFraction(nPhaseIdx,wCompIdx) * FluidSystem::molarMass(wCompIdx)
                    + (1.-fluidState.moleFraction(nPhaseIdx,wCompIdx)) * FluidSystem::molarMass(nCompIdx) )));

        // complete array of mass fractions
        fluidState.setMassFraction(wPhaseIdx, nCompIdx, 1. - fluidState.massFraction(wPhaseIdx,wCompIdx));
        fluidState.setMassFraction(nPhaseIdx, nCompIdx, 1. - fluidState.massFraction(nPhaseIdx,wCompIdx));
        // complete array of mole fractions
        fluidState.setMoleFraction(wPhaseIdx, nCompIdx, 1. - fluidState.moleFraction(wPhaseIdx,wCompIdx));
        fluidState.setMoleFraction(nPhaseIdx, nCompIdx, 1. - fluidState.moleFraction(nPhaseIdx,wCompIdx));

        // get densities with correct composition
        fluidState.setDensity(wPhaseIdx, FluidSystem::density(fluidState, paramCache, wPhaseIdx));
        fluidState.setDensity(nPhaseIdx, FluidSystem::density(fluidState, paramCache, nPhaseIdx));

        // set saturation
        fluidState.setSaturation(wPhaseIdx, saturation);
    }
//@}
};

} // end namespace Dumux

#endif
