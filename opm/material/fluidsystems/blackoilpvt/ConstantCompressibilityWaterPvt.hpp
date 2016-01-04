// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#ifndef OPM_CONSTANT_COMPRESSIBILITY_WATER_PVT_HPP
#define OPM_CONSTANT_COMPRESSIBILITY_WATER_PVT_HPP

#include <opm/material/common/Tabulated1DFunction.hpp>

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
class ConstantCompressibilityWaterPvt
{
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

public:
#if HAVE_OPM_PARSER
    /*!
     * \brief Sets the pressure-dependent water viscosity and density
     *        using a table stemming from the Eclipse PVTW keyword.
     */
    void initFromDeck(DeckConstPtr deck, EclipseStateConstPtr /*eclState*/)
    {
        DeckKeywordConstPtr pvtwKeyword = deck->getKeyword("PVTW");
        DeckKeywordConstPtr densityKeyword = deck->getKeyword("DENSITY");

        assert(pvtwKeyword->size() == densityKeyword->size());

        size_t numRegions = pvtwKeyword->size();
        setNumRegions(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            auto pvtwRecord = pvtwKeyword->getRecord(regionIdx);
            auto densityRecord = densityKeyword->getRecord(regionIdx);

            waterReferenceDensity_[regionIdx] =
                densityRecord->getItem("WATER")->getSIDouble(0);

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
    }
#endif

    void setNumRegions(size_t numRegions)
    {
        waterReferenceDensity_.resize(numRegions);
        waterReferencePressure_.resize(numRegions);
        waterReferenceFormationVolumeFactor_.resize(numRegions);
        waterCompressibility_.resize(numRegions);
        waterViscosity_.resize(numRegions);
        waterViscosibility_.resize(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            setReferenceDensities(regionIdx, 650.0, 1.0, 1000.0);
            setReferenceFormationVolumeFactor(regionIdx, 1.0);
            setReferencePressure(regionIdx, 1e5);
        }
    }

    /*!
     * \brief Set the water reference density [kg / m^3]
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar /*rhoRefOil*/,
                               Scalar /*rhoRefGas*/,
                               Scalar rhoRefWater)
    { waterReferenceDensity_[regionIdx] = rhoRefWater; }

    /*!
     * \brief Set the water reference pressure [Pa]
     */
    void setReferencePressure(unsigned regionIdx, Scalar p)
    { waterReferencePressure_[regionIdx] = p; }

    /*!
     * \brief Set the viscosity and "viscosibility" of the water phase.
     */
    void setViscosity(unsigned regionIdx, Scalar muw, Scalar waterViscosibility = 0.0)
    {
        waterViscosity_[regionIdx] = muw;
        waterViscosibility_[regionIdx] = waterViscosibility;
    }

    /*!
     * \brief Set the compressibility of the water phase.
     */
    void setCompressibility(unsigned regionIdx, Scalar waterCompressibility)
    { waterCompressibility_[regionIdx] = waterCompressibility; }

    /*!
     * \brief Set the water reference formation volume factor [-]
     */
    void setReferenceFormationVolumeFactor(unsigned regionIdx, Scalar BwRef)
    { waterReferenceFormationVolumeFactor_[regionIdx] = BwRef; }

    /*!
     * \brief Set the water "viscosibility" [1/ (Pa s)]
     */
    void setViscosibility(unsigned regionIdx, Scalar muComp)
    { waterViscosibility_[regionIdx] = muComp; }

    /*!
     * \brief Finish initializing the water phase PVT properties.
     */
    void initEnd()
    { }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure) const
    {
        // Eclipse calculates the viscosity in a weird way: it
        // calcultes the product of B_w and mu_w and then divides the
        // result by B_w...
        Scalar muwRef = waterViscosity_[regionIdx];

        // note: this is NOT equivalent to the equation given by the ECL RM. It is
        // equivalent to the code which was used by opm-core at the time when this was
        // written.
        Scalar pRef = waterReferencePressure_[regionIdx];
        const Evaluation& x = (-waterViscosibility_[regionIdx])*(pressure - pRef);
        const Evaluation& d = 1.0 + x*(1.0 + x/2.0);
        return muwRef/d;
    }

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation density(unsigned regionIdx,
                       const Evaluation& temperature,
                       const Evaluation& pressure) const
    {
        const Evaluation& Bw = formationVolumeFactor(regionIdx, temperature, pressure);
        return waterReferenceDensity_[regionIdx]/Bw;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation formationVolumeFactor(unsigned regionIdx,
                                     const Evaluation& /*temperature*/,
                                     const Evaluation& pressure) const
    {
        // cf. ECLiPSE 2011 technical description, p. 116
        Scalar pRef = waterReferencePressure_[regionIdx];
        const Evaluation& X = waterCompressibility_[regionIdx]*(pressure - pRef);

        Scalar BwRef = waterReferenceFormationVolumeFactor_[regionIdx];

        // TODO (?): consider the salt concentration of the brine
        return BwRef/(1 + X*(1 + X/2));
    }

private:
    std::vector<Scalar> waterReferenceDensity_;
    std::vector<Scalar> waterReferencePressure_;
    std::vector<Scalar> waterReferenceFormationVolumeFactor_;
    std::vector<Scalar> waterCompressibility_;
    std::vector<Scalar> waterViscosity_;
    std::vector<Scalar> waterViscosibility_;
};

} // namespace Opm

#endif
