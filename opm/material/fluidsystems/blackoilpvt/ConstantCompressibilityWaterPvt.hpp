/*
  Copyright (C) 2015 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \copydoc Opm::ConstantCompressibilityWaterPvt
 */
#ifndef OPM_CONSTANT_COMPRESSIBILITY_WATER_HPP
#define OPM_CONSTANT_COMPRESSIBILITY_WATER_HPP

#include "WaterPvtInterface.hpp"

#include <opm/material/OpmFinal.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

#include <vector>

namespace Opm {

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phase
 *        without vaporized oil.
 */
template <class Scalar>
class ConstantCompressibilityWaterPvt : public WaterPvtInterface<Scalar>
{
    typedef FluidSystems::BlackOil<Scalar> BlackOilFluidSystem;

    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

    static const int oilPhaseIdx = BlackOilFluidSystem::oilPhaseIdx;
    static const int gasPhaseIdx = BlackOilFluidSystem::gasPhaseIdx;
    static const int waterPhaseIdx = BlackOilFluidSystem::waterPhaseIdx;

    static const int oilCompIdx = BlackOilFluidSystem::oilCompIdx;
    static const int gasCompIdx = BlackOilFluidSystem::gasCompIdx;
    static const int waterCompIdx = BlackOilFluidSystem::waterCompIdx;

public:
    void setNumRegions(int numRegions)
    {
        waterReferencePressure_.resize(numRegions);
        waterReferenceFormationVolumeFactor_.resize(numRegions);
        waterCompressibility_.resize(numRegions);
        waterViscosity_.resize(numRegions);
        waterViscosibility_.resize(numRegions);

        for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            setReferenceFormationVolumeFactor(regionIdx, 1.0);
            setReferencePressure(regionIdx, BlackOilFluidSystem::surfacePressure);
        }
    }

#if HAVE_OPM_PARSER
    /*!
     * \brief Sets the pressure-dependent water viscosity and density
     *        using a table stemming from the Eclipse PVTW keyword.
     */
    void setPvtw(int regionIdx, DeckKeywordConstPtr pvtwKeyword)
    {
        assert(static_cast<int>(pvtwKeyword->size()) >= regionIdx);

        auto pvtwRecord = pvtwKeyword->getRecord(regionIdx);
        waterReferencePressure_[regionIdx] =
            pvtwRecord->getItem("P_REF")->getSIDouble(0);
        waterReferenceFormationVolumeFactor_[regionIdx] =
            pvtwRecord->getItem("WATER_VOL_FACTOR")->getSIDouble(0);
        waterCompressibility_[regionIdx] =
            pvtwRecord->getItem("WATER_COMPRESSIBILITY")->getSIDouble(0);
        waterViscosity_[regionIdx] =
            pvtwRecord->getItem("WATER_VISCOSITY")->getSIDouble(0);
        waterViscosibility_[regionIdx] =
            pvtwRecord->getItem("WATER_VISCOSIBILITY")->getSIDouble(0);
    }
#endif

    /*!
     * \brief Set the viscosity and "viscosibility" of the water phase.
     */
    void setViscosity(int regionIdx, Scalar muw, Scalar waterViscosibility = 0.0)
    {
        waterViscosity_[regionIdx] = muw;
        waterViscosibility_[regionIdx] = waterViscosibility;
    }

    /*!
     * \brief Set the compressibility of the water phase.
     */
    void setCompressibility(int regionIdx, Scalar waterCompressibility)
    { waterCompressibility_[regionIdx] = waterCompressibility; }

    /*!
     * \brief Set the water reference pressure [Pa]
     */
    void setReferencePressure(int regionIdx, Scalar p)
    { waterReferencePressure_[regionIdx] = p; }

    /*!
     * \brief Set the water reference formation volume factor [-]
     */
    void setReferenceFormationVolumeFactor(int regionIdx, Scalar BwRef)
    { waterReferenceFormationVolumeFactor_[regionIdx] = BwRef; }

    /*!
     * \brief Set the water "viscosibility" [1/ (Pa s)]
     */
    void setViscosibility(int regionIdx, Scalar muComp)
    { waterViscosibility_[regionIdx] = muComp; }

    /*!
     * \brief Finish initializing the water phase PVT properties.
     */
    void initEnd()
    { }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    Scalar viscosity(int regionIdx,
                     Scalar temperature,
                     Scalar pressure) const OPM_FINAL
    {
        // Eclipse calculates the viscosity in a weird way: it
        // calcultes the product of B_w and mu_w and then divides the
        // result by B_w...
        Scalar BwMuwRef = waterViscosity_[regionIdx]*waterReferenceFormationVolumeFactor_[regionIdx];
        Scalar Bw = formationVolumeFactor(regionIdx, temperature, pressure);

        Scalar pRef = waterReferencePressure_[regionIdx];
        Scalar Y =
            (waterCompressibility_[regionIdx] - waterViscosibility_[regionIdx])
            * (pressure - pRef);
        return BwMuwRef/((1 + Y*(1 + Y/2))*Bw);
    }

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    Scalar density(int regionIdx,
                   Scalar temperature,
                   Scalar pressure) const OPM_FINAL
    {
        Scalar Bw = formationVolumeFactor(regionIdx, temperature, pressure);
        Scalar rhowRef = BlackOilFluidSystem::referenceDensity(waterPhaseIdx, regionIdx);
        return rhowRef/Bw;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    Scalar formationVolumeFactor(int regionIdx,
                                 Scalar temperature,
                                 Scalar pressure) const OPM_FINAL
    {
        // cf. ECLiPSE 2011 technical description, p. 116
        Scalar pRef = waterReferencePressure_[regionIdx];
        Scalar X = waterCompressibility_[regionIdx]*(pressure - pRef);

        Scalar BwRef = waterReferenceFormationVolumeFactor_[regionIdx];

        // TODO (?): consider the salt concentration of the brine
        return BwRef/(1 + X*(1 + X/2));
    }

    /*!
     * \brief Returns the fugacity coefficient [Pa] of a component in the fluid phase given
     *        a set of parameters.
     */
    Scalar fugacityCoefficient(int regionIdx,
                               Scalar temperature,
                               Scalar pressure,
                               int compIdx) const OPM_FINAL
    {
        // set the affinity of the gas and oil components to the water phase to be 10
        // orders of magnitute smaller than that of the water component. for this we use
        // a pseudo-realistic vapor pressure of water as a starting point. (we just set
        // it to 30 kPa to ease interpreting the results.)
        const Scalar pvWater = 30e3;
        if (compIdx == BlackOilFluidSystem::oilCompIdx)
            return 1e10*pvWater / pressure;
        else if (compIdx == BlackOilFluidSystem::gasCompIdx)
            return 1.01e10*pvWater / pressure;

        return pvWater / pressure;
    }

private:
    std::vector<Scalar> waterReferencePressure_;
    std::vector<Scalar> waterReferenceFormationVolumeFactor_;
    std::vector<Scalar> waterCompressibility_;
    std::vector<Scalar> waterViscosity_;
    std::vector<Scalar> waterViscosibility_;
};

} // namespace Opm

#endif
