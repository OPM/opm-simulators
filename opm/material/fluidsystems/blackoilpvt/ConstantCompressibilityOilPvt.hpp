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
 * \copydoc Opm::ConstantCompressibilityOilPvt
 */
#ifndef OPM_CONSTANT_COMPRESSIBILITY_OIL_PVT_HPP
#define OPM_CONSTANT_COMPRESSIBILITY_OIL_PVT_HPP

#include <opm/material/common/OpmFinal.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/Spline.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Deck/DeckRecord.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil phase
 *        without dissolved gas and constant compressibility/"viscosibility".
 */
template <class Scalar>
class ConstantCompressibilityOilPvt
{
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

public:
#if HAVE_OPM_PARSER
    /*!
     * \brief Sets the pressure-dependent oil viscosity and density
     *        using the Eclipse PVCDO keyword.
     */
    /*!
     * \brief Initialize the oil parameters via the data specified by the PVTO ECL keyword.
     */
    void initFromDeck(DeckConstPtr deck, EclipseStateConstPtr /*eclState*/)
    {
        DeckKeywordConstPtr pvcdoKeyword = deck->getKeyword("PVCDO");
        DeckKeywordConstPtr densityKeyword = deck->getKeyword("DENSITY");

        assert(pvcdoKeyword->size() == densityKeyword->size());

        size_t numRegions = pvcdoKeyword->size();
        setNumRegions(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            Scalar rhoRefO = densityKeyword->getRecord(regionIdx)->getItem("OIL")->getSIDouble(0);
            Scalar rhoRefG = densityKeyword->getRecord(regionIdx)->getItem("GAS")->getSIDouble(0);
            Scalar rhoRefW = densityKeyword->getRecord(regionIdx)->getItem("WATER")->getSIDouble(0);

            setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);

            auto pvcdoRecord = pvcdoKeyword->getRecord(regionIdx);
            oilReferencePressure_[regionIdx] =
                pvcdoRecord->getItem("P_REF")->getSIDouble(0);
            oilReferenceFormationVolumeFactor_[regionIdx] =
                pvcdoRecord->getItem("OIL_VOL_FACTOR")->getSIDouble(0);
            oilCompressibility_[regionIdx] =
                pvcdoRecord->getItem("OIL_COMPRESSIBILITY")->getSIDouble(0);
            oilViscosity_[regionIdx] =
                pvcdoRecord->getItem("OIL_VISCOSITY")->getSIDouble(0);
            oilViscosibility_[regionIdx] =
                pvcdoRecord->getItem("OIL_VISCOSIBILITY")->getSIDouble(0);
        }
    }
#endif

    void setNumRegions(size_t numRegions)
    {
        oilReferenceDensity_.resize(numRegions);
        oilReferencePressure_.resize(numRegions);
        oilReferenceFormationVolumeFactor_.resize(numRegions);
        oilCompressibility_.resize(numRegions);
        oilViscosity_.resize(numRegions);
        oilViscosibility_.resize(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            setReferenceFormationVolumeFactor(regionIdx, 1.0);
            setReferencePressure(regionIdx, 1.03125);
        }
    }

    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar rhoRefOil,
                               Scalar /*rhoRefGas*/,
                               Scalar /*rhoRefWater*/)
    { oilReferenceDensity_[regionIdx] = rhoRefOil; }

    /*!
     * \brief Set the viscosity and "viscosibility" of the oil phase.
     */
    void setViscosity(unsigned regionIdx, Scalar muo, Scalar oilViscosibility = 0.0)
    {
        oilViscosity_[regionIdx] = muo;
        oilViscosibility_[regionIdx] = oilViscosibility;
    }

    /*!
     * \brief Set the compressibility of the oil phase.
     */
    void setCompressibility(unsigned regionIdx, Scalar oilCompressibility)
    { oilCompressibility_[regionIdx] = oilCompressibility; }

    /*!
     * \brief Set the oil reference pressure [Pa]
     */
    void setReferencePressure(unsigned regionIdx, Scalar p)
    { oilReferencePressure_[regionIdx] = p; }

    /*!
     * \brief Set the oil reference formation volume factor [-]
     */
    void setReferenceFormationVolumeFactor(unsigned regionIdx, Scalar BoRef)
    { oilReferenceFormationVolumeFactor_[regionIdx] = BoRef; }

    /*!
     * \brief Set the oil "viscosibility" [1/ (Pa s)]
     */
    void setViscosibility(unsigned regionIdx, Scalar muComp)
    { oilViscosibility_[regionIdx] = muComp; }

    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd()
    { }

    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& temperature,
                         const Evaluation& pressure,
                         const Evaluation& /*Rs*/) const
    { return saturatedViscosity(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of gas saturated oil given a pressure.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& pressure) const
    {
        // Eclipse calculates the viscosity in a weird way: it
        // calcultes the product of B_w and mu_w and then divides the
        // result by B_w...
        Scalar BoMuoRef = oilViscosity_[regionIdx]*oilReferenceFormationVolumeFactor_[regionIdx];
        const Evaluation& Bo = saturatedFormationVolumeFactor(regionIdx, temperature, pressure);

        Scalar pRef = oilReferencePressure_[regionIdx];
        const Evaluation& Y =
            (oilCompressibility_[regionIdx] - oilViscosibility_[regionIdx])
            * (pressure - pRef);
        return BoMuoRef/((1 + Y*(1 + Y/2))*Bo);
    }

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation density(unsigned regionIdx,
                       const Evaluation& temperature,
                       const Evaluation& pressure,
                       const Evaluation& /*Rs*/) const
    { return saturatedDensity(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the density [kg/m^3] of gas saturated oil given a pressure.
     */
    template <class Evaluation>
    Evaluation saturatedDensity(unsigned regionIdx,
                                const Evaluation& temperature,
                                const Evaluation& pressure) const
    {
        Scalar rhooRef = oilReferenceDensity_[regionIdx];

        const Evaluation& Bo = saturatedFormationVolumeFactor(regionIdx, temperature, pressure);
        return rhooRef/Bo;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation formationVolumeFactor(unsigned regionIdx,
                                     const Evaluation& temperature,
                                     const Evaluation& pressure,
                                     const Evaluation& /*Rs*/) const
    { return saturatedFormationVolumeFactor(regionIdx, temperature, pressure); }

    /*!
     * \brief Returns the formation volume factor [-] of gas saturated oil.
     *
     * Note that constant compressibility oil is a special case of dead oil and dead oil
     * is always gas saturated by by definition.
     */
    template <class Evaluation>
    Evaluation saturatedFormationVolumeFactor(unsigned regionIdx,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& pressure) const
    {
        // cf. ECLiPSE 2011 technical description, p. 116
        Scalar pRef = oilReferencePressure_[regionIdx];
        const Evaluation& X = oilCompressibility_[regionIdx]*(pressure - pRef);

        Scalar BoRef = oilReferenceFormationVolumeFactor_[regionIdx];
        return BoRef/(1 + X*(1 + X/2));
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned /*regionIdx*/,
                                             const Evaluation& /*temperature*/,
                                             const Evaluation& /*pressure*/) const
    { return 0.0; /* this is dead oil! */ }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param Rs The surface volume of gas component dissolved in what will yield one cubic meter of oil at the surface [-]
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned /*regionIdx*/,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& /*Rs*/) const
    { return 0.0; /* this is dead oil, so there isn't any meaningful saturation pressure! */ }

private:
    std::vector<Scalar> oilReferenceDensity_;
    std::vector<Scalar> oilReferencePressure_;
    std::vector<Scalar> oilReferenceFormationVolumeFactor_;
    std::vector<Scalar> oilCompressibility_;
    std::vector<Scalar> oilViscosity_;
    std::vector<Scalar> oilViscosibility_;
};

} // namespace Opm

#endif
