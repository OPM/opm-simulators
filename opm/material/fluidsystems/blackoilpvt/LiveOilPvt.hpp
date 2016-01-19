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
 * \copydoc Opm::LiveOilPvt
 */
#ifndef OPM_LIVE_OIL_PVT_HPP
#define OPM_LIVE_OIL_PVT_HPP

#include <opm/material/Constants.hpp>

#include <opm/material/common/OpmFinal.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/Spline.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil phas
 *        with dissolved gas.
 */
template <class Scalar>
class LiveOilPvt
{
    typedef Opm::UniformXTabulated2DFunction<Scalar> TabulatedTwoDFunction;
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef Opm::Spline<Scalar> Spline;
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

public:
#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize the oil parameters via the data specified by the PVTO ECL keyword.
     */
    void initFromDeck(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        const auto& pvtoTables = eclState->getTableManager()->getPvtoTables();
        DeckKeywordConstPtr densityKeyword = deck->getKeyword("DENSITY");

        assert(pvtoTables.size() == densityKeyword->size());

        size_t numRegions = pvtoTables.size();
        setNumRegions(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            Scalar rhoRefO = densityKeyword->getRecord(regionIdx)->getItem("OIL")->getSIDouble(0);
            Scalar rhoRefG = densityKeyword->getRecord(regionIdx)->getItem("GAS")->getSIDouble(0);
            Scalar rhoRefW = densityKeyword->getRecord(regionIdx)->getItem("WATER")->getSIDouble(0);

            setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);
        }

        // initialize the internal table objects
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            const auto& pvtoTable = pvtoTables[regionIdx];

            const auto saturatedTable = pvtoTable.getSaturatedTable();
            assert(saturatedTable.numRows() > 1);

            auto& oilMu = oilMuTable_[regionIdx];
            auto& satOilMu = saturatedOilMuTable_[regionIdx];
            auto& invOilB = inverseOilBTable_[regionIdx];
            auto& invSatOilB = inverseSaturatedOilBTable_[regionIdx];
            auto& gasDissolutionFac = saturatedGasDissolutionFactorTable_[regionIdx];
            std::vector<Scalar> invSatOilBArray;
            std::vector<Scalar> satOilMuArray;

            // extract the table for the gas dissolution and the oil formation volume factors
            for (unsigned outerIdx = 0; outerIdx < saturatedTable.numRows(); ++ outerIdx) {
                Scalar Rs    = saturatedTable.get("RS", outerIdx);
                Scalar BoSat = saturatedTable.get("BO", outerIdx);
                Scalar muoSat = saturatedTable.get("MU", outerIdx);

                satOilMuArray.push_back(muoSat);
                invSatOilBArray.push_back(1.0/BoSat);

                invOilB.appendXPos(Rs);
                oilMu.appendXPos(Rs);

                assert(invOilB.numX() == outerIdx + 1);
                assert(oilMu.numX() == outerIdx + 1);

                const auto& underSaturatedTable = pvtoTable.getUnderSaturatedTable(outerIdx);
                size_t numRows = underSaturatedTable.numRows();
                for (unsigned innerIdx = 0; innerIdx < numRows; ++ innerIdx) {
                    Scalar po = underSaturatedTable.get("P", innerIdx);
                    Scalar Bo = underSaturatedTable.get("BO", innerIdx);
                    Scalar muo = underSaturatedTable.get("MU", innerIdx);

                    invOilB.appendSamplePoint(outerIdx, po, 1.0/Bo);
                    oilMu.appendSamplePoint(outerIdx, po, muo);
                }
            }

            // update the tables for the formation volume factor and for the gas
            // dissolution factor of saturated oil
            {
                std::vector<double> tmpPressureColumn = saturatedTable.getColumn("P").vectorCopy();
                std::vector<double> tmpGasSolubilityColumn = saturatedTable.getColumn("RS").vectorCopy();
                std::vector<double> tmpMuColumn = saturatedTable.getColumn("MU").vectorCopy();

                invSatOilB.setXYContainers(tmpPressureColumn, invSatOilBArray);
                satOilMu.setXYContainers(tmpPressureColumn, satOilMuArray);
                gasDissolutionFac.setXYContainers(tmpPressureColumn, tmpGasSolubilityColumn);
            }

            updateSaturationPressureSpline_(regionIdx);
            // make sure to have at least two sample points per Rs value
            for (unsigned xIdx = 0; xIdx < invOilB.numX(); ++xIdx) {
                // a single sample point is definitely needed
                assert(invOilB.numY(xIdx) > 0);

                // everything is fine if the current table has two or more sampling points
                // for a given mole fraction
                if (invOilB.numY(xIdx) > 1)
                    continue;

                // find the master table which will be used as a template to extend the
                // current line. We define master table as the first table which has values
                // for undersaturated oil...
                size_t masterTableIdx = xIdx + 1;
                for (; masterTableIdx < saturatedTable.numRows(); ++masterTableIdx)
                {
                    if (pvtoTable.getUnderSaturatedTable(masterTableIdx).numRows() > 1)
                        break;
                }

                if (masterTableIdx >= saturatedTable.numRows())
                    OPM_THROW(std::runtime_error,
                              "PVTO tables are invalid: The last table must exhibit at least one "
                              "entry for undersaturated oil!");

                // extend the current table using the master table.
                extendPvtoTable_(regionIdx,
                                 xIdx,
                                 pvtoTable.getUnderSaturatedTable(xIdx),
                                 pvtoTable.getUnderSaturatedTable(masterTableIdx));
            }
        }
    }

private:
    void extendPvtoTable_(unsigned regionIdx,
                          unsigned xIdx,
                          const SimpleTable& curTable,
                          const SimpleTable& masterTable)
    {
        std::vector<double> pressuresArray = curTable.getColumn("P").vectorCopy();
        std::vector<double> oilBArray = curTable.getColumn("BO").vectorCopy();
        std::vector<double> oilMuArray = curTable.getColumn("MU").vectorCopy();

        auto& invOilB = inverseOilBTable_[regionIdx];
        auto& oilMu = oilMuTable_[regionIdx];

        for (unsigned newRowIdx = 1; newRowIdx < masterTable.numRows(); ++ newRowIdx) {
            const auto& pressureColumn = masterTable.getColumn("P");
            const auto& BOColumn = masterTable.getColumn("BO");
            const auto& viscosityColumn = masterTable.getColumn("MU");

            // compute the oil pressure for the new entry
            Scalar diffPo = pressureColumn[newRowIdx] - pressureColumn[newRowIdx - 1];
            Scalar newPo = pressuresArray.back() + diffPo;

            // calculate the compressibility of the master table
            Scalar B1 = BOColumn[newRowIdx];
            Scalar B2 = BOColumn[newRowIdx - 1];
            Scalar x = (B1 - B2)/( (B1 + B2)/2.0 );

            // calculate the oil formation volume factor which exhibits the same
            // compressibility for the new pressure
            Scalar newBo = oilBArray.back()*(1.0 + x/2.0)/(1.0 - x/2.0);

            // calculate the "viscosibility" of the master table
            Scalar mu1 = viscosityColumn[newRowIdx];
            Scalar mu2 = viscosityColumn[newRowIdx - 1];
            Scalar xMu = (mu1 - mu2)/( (mu1 + mu2)/2.0 );

            // calculate the oil formation volume factor which exhibits the same
            // compressibility for the new pressure
            Scalar newMuo = oilMuArray.back()*(1.0 + xMu/2)/(1.0 - xMu/2.0);

            // append the new values to the arrays which we use to compute the additional
            // values ...
            pressuresArray.push_back(newPo);
            oilBArray.push_back(newBo);
            oilMuArray.push_back(newMuo);

            // ... and register them with the internal table objects
            invOilB.appendSamplePoint(xIdx, newPo, 1.0/newBo);
            oilMu.appendSamplePoint(xIdx, newPo, newMuo);
        }
    }

public:
#endif // HAVE_OPM_PARSER

    void setNumRegions(size_t numRegions)
    {
        oilReferenceDensity_.resize(numRegions);
        gasReferenceDensity_.resize(numRegions);
        inverseOilBTable_.resize(numRegions);
        inverseOilBMuTable_.resize(numRegions);
        inverseSaturatedOilBTable_.resize(numRegions);
        inverseSaturatedOilBMuTable_.resize(numRegions);
        oilMuTable_.resize(numRegions);
        saturatedOilMuTable_.resize(numRegions);
        saturatedGasDissolutionFactorTable_.resize(numRegions);
        saturationPressureSpline_.resize(numRegions);
    }

    /*!
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setReferenceDensities(unsigned regionIdx,
                               Scalar rhoRefOil,
                               Scalar rhoRefGas,
                               Scalar /*rhoRefWater*/)
    {
        oilReferenceDensity_[regionIdx] = rhoRefOil;
        gasReferenceDensity_[regionIdx] = rhoRefGas;
    }

    /*!
     * \brief Initialize the function for the gas dissolution factor \f$R_s\f$
     *
     * \param samplePoints A container of (x,y) values.
     */
    void setSaturatedOilGasDissolutionFactor(unsigned regionIdx, const SamplingPoints &samplePoints)
    { saturatedGasDissolutionFactorTable_[regionIdx].setContainerOfTuples(samplePoints); }

    /*!
     * \brief Initialize the function for the oil formation volume factor
     *
     * The oil formation volume factor \f$B_o\f$ is a function of \f$(p_o, X_o^G)\f$ and
     * represents the partial density of the oil component in the oil phase at a given
     * pressure. This method only requires the volume factor of gas-saturated oil (which
     * only depends on pressure) while the dependence on the gas mass fraction is
     * guesstimated.
     */
    void setSaturatedOilFormationVolumeFactor(unsigned regionIdx, const SamplingPoints &samplePoints)
    {
        Scalar T = 273.15 + 15.56; // [K]
        auto& invOilB = inverseOilBTable_[regionIdx];

        updateSaturationPressureSpline_(regionIdx);

        // calculate a table of estimated densities of undersatured gas
        for (size_t pIdx = 0; pIdx < samplePoints.size(); ++pIdx) {
            Scalar p1 = std::get<0>(samplePoints[pIdx]);
            Scalar p2 = p1 * 2.0;

            Scalar Bo1 = std::get<1>(samplePoints[pIdx]);
            Scalar drhoo_dp = (1.1200 - 1.1189)/((5000 - 4000)*6894.76);
            Scalar Bo2 = Bo1/(1.0 + (p2 - p1)*drhoo_dp);

            Scalar Rs = saturatedGasDissolutionFactor(regionIdx, T, p1);

            invOilB.appendXPos(Rs);
            invOilB.appendSamplePoint(pIdx, p1, 1.0/Bo1);
            invOilB.appendSamplePoint(pIdx, p2, 1.0/Bo2);
        }
    }

    /*!
     * \brief Initialize the spline for the oil formation volume factor
     *
     * The oil formation volume factor \f$B_o\f$ is a function of \f$(p_o, X_o^G)\f$ and
     * represents the partial density of the oil component in the oil phase at a given
     * pressure.
     *
     * This method sets \f$1/B_o(R_s, p_o)\f$. Note that instead of the mass fraction of
     * the gas component in the oil phase, this function depends on the gas dissolution
     * factor. Also note, that the order of the arguments needs to be \f$(R_s, p_o)\f$
     * and not the other way around.
     */
    void setInverseOilFormationVolumeFactor(unsigned regionIdx, const TabulatedTwoDFunction& invBo)
    { inverseOilBTable_[regionIdx] = invBo; }

    /*!
     * \brief Initialize the viscosity of the oil phase.
     *
     * This is a function of \f$(R_s, p_o)\f$...
     */
    void setOilViscosity(unsigned regionIdx, const TabulatedTwoDFunction& muo)
    { oilMuTable_[regionIdx] = muo; }

    /*!
     * \brief Initialize the phase viscosity for gas saturated oil
     *
     * The oil viscosity is a function of \f$(p_o, X_o^G)\f$, but this method only
     * requires the viscosity of gas-saturated oil (which only depends on pressure) while
     * there is assumed to be no dependence on the gas mass fraction...
     */
    void setSaturatedOilViscosity(unsigned regionIdx, const SamplingPoints &samplePoints)
    {
        Scalar T = 273.15 + 15.56; // [K]

        // update the table for the saturated oil
        saturatedOilMuTable_[regionIdx].setContainerOfTuples(samplePoints);

        // calculate a table of estimated viscosities depending on pressure and gas mass
        // fraction for untersaturated oil to make the other code happy
        for (size_t pIdx = 0; pIdx < samplePoints.size(); ++pIdx) {
            Scalar p1 = std::get<0>(samplePoints[pIdx]);
            Scalar p2 = p1 * 2.0;

            // no pressure dependence of the viscosity
            Scalar mu1 = std::get<1>(samplePoints[pIdx]);
            Scalar mu2 = mu1;

            Scalar Rs = saturatedGasDissolutionFactor(regionIdx, T, p1);

            oilMuTable_[regionIdx].appendXPos(Rs);
            oilMuTable_[regionIdx].appendSamplePoint(pIdx, p1, mu1);
            oilMuTable_[regionIdx].appendSamplePoint(pIdx, p2, mu2);
        }
    }

    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd()
    {
        // calculate the final 2D functions which are used for interpolation.
        size_t numRegions = oilMuTable_.size();
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            // calculate the table which stores the inverse of the product of the oil
            // formation volume factor and the oil viscosity
            const auto& oilMu = oilMuTable_[regionIdx];
            const auto& satOilMu = saturatedOilMuTable_[regionIdx];
            const auto& invOilB = inverseOilBTable_[regionIdx];
            assert(oilMu.numX() == invOilB.numX());

            auto& invOilBMu = inverseOilBMuTable_[regionIdx];
            auto& invSatOilB = inverseSaturatedOilBTable_[regionIdx];
            auto& invSatOilBMu = inverseSaturatedOilBMuTable_[regionIdx];

            std::vector<Scalar> satPressuresArray;
            std::vector<Scalar> invSatOilBArray;
            std::vector<Scalar> invSatOilBMuArray;
            for (unsigned rsIdx = 0; rsIdx < oilMu.numX(); ++rsIdx) {
                invOilBMu.appendXPos(oilMu.xAt(rsIdx));

                assert(oilMu.numY(rsIdx) == invOilB.numY(rsIdx));

                size_t numPressures = oilMu.numY(rsIdx);
                for (unsigned pIdx = 0; pIdx < numPressures; ++pIdx)
                    invOilBMu.appendSamplePoint(rsIdx,
                                                oilMu.yAt(rsIdx, pIdx),
                                                invOilB.valueAt(rsIdx, pIdx)
                                                / oilMu.valueAt(rsIdx, pIdx));

                // the sampling points in UniformXTabulated2DFunction are always sorted
                // in ascending order. Thus, the value for saturated oil is the first one
                // (i.e., the one for the lowest pressure value)
                satPressuresArray.push_back(oilMu.yAt(rsIdx, 0));
                invSatOilBArray.push_back(invOilB.valueAt(rsIdx, 0));
                invSatOilBMuArray.push_back(invSatOilBArray.back()/satOilMu.valueAt(rsIdx));
            }

            invSatOilB.setXYContainers(satPressuresArray, invSatOilBArray);
            invSatOilBMu.setXYContainers(satPressuresArray, invSatOilBMuArray);

            updateSaturationPressureSpline_(regionIdx);
        }
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation viscosity(unsigned regionIdx,
                         const Evaluation& /*temperature*/,
                         const Evaluation& pressure,
                         const Evaluation& Rs) const
    {
        // ATTENTION: Rs is the first axis!
        const Evaluation& invBo = inverseOilBTable_[regionIdx].eval(Rs, pressure, /*extrapolate=*/true);
        const Evaluation& invMuoBo = inverseOilBMuTable_[regionIdx].eval(Rs, pressure, /*extrapolate=*/true);

        return invBo/invMuoBo;
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& pressure) const
    {
        // ATTENTION: Rs is the first axis!
        const Evaluation& invBo = inverseSaturatedOilBTable_[regionIdx].eval(pressure, /*extrapolate=*/true);
        const Evaluation& invMuoBo = inverseSaturatedOilBMuTable_[regionIdx].eval(pressure, /*extrapolate=*/true);

        return invBo/invMuoBo;
    }

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation density(unsigned regionIdx,
                       const Evaluation& temperature,
                       const Evaluation& pressure,
                       const Evaluation& Rs) const
    {
        Scalar rhooRef = oilReferenceDensity_[regionIdx];
        Scalar rhogRef = gasReferenceDensity_[regionIdx];
        Valgrind::CheckDefined(rhooRef);
        Valgrind::CheckDefined(rhogRef);

        const Evaluation& Bo = formationVolumeFactor(regionIdx, temperature, pressure, Rs);
        Valgrind::CheckDefined(Bo);

        Evaluation rhoo = rhooRef/Bo;

        // the oil formation volume factor just represents the partial density of the oil
        // component in the oil phase. to get the total density of the phase, we have to
        // add the partial density of the gas component.
        rhoo += rhogRef*Rs/Bo;

        return rhoo;
    }

    /*!
     * \brief Returns the density [kg/m^3] of gas saturated oil for a given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedDensity(unsigned regionIdx,
                                const Evaluation& temperature,
                                const Evaluation& pressure) const
    {
        Scalar rhooRef = oilReferenceDensity_[regionIdx];
        Scalar rhogRef = gasReferenceDensity_[regionIdx];
        Valgrind::CheckDefined(rhooRef);
        Valgrind::CheckDefined(rhogRef);

        const Evaluation& Bo = saturatedFormationVolumeFactor(regionIdx, temperature, pressure);
        Valgrind::CheckDefined(Bo);

        Evaluation rhoo = rhooRef/Bo;

        // the oil formation volume factor just represents the partial density of the oil
        // component in the oil phase. to get the total density of the phase, we have to
        // add the partial density of the gas component.
        const Evaluation& RsSat = saturatedGasDissolutionFactor(regionIdx, temperature, pressure);
        rhoo += rhogRef*RsSat/Bo;

        return rhoo;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation formationVolumeFactor(unsigned regionIdx,
                                     const Evaluation& /*temperature*/,
                                     const Evaluation& pressure,
                                     const Evaluation& Rs) const
    {
        // ATTENTION: Rs is represented by the _first_ axis!
        return 1.0 / inverseOilBTable_[regionIdx].eval(Rs, pressure, /*extrapolate=*/true);
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation saturatedFormationVolumeFactor(unsigned regionIdx,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& pressure) const
    {
        // ATTENTION: Rs is represented by the _first_ axis!
        return 1.0 / inverseSaturatedOilBTable_[regionIdx].eval(pressure, /*extrapolate=*/true);
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& /*temperature*/,
                                             const Evaluation& pressure) const
    { return saturatedGasDissolutionFactorTable_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param Rs The surface volume of gas component dissolved in what will yield one cubic meter of oil at the surface [-]
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature,
                                  const Evaluation& Rs) const
    {
        typedef Opm::MathToolbox<Evaluation> Toolbox;

        // use the saturation pressure spline to get a pretty good initial value
        Evaluation pSat = saturationPressureSpline_[regionIdx].eval(Rs, /*extrapolate=*/true);
        Evaluation eps = pSat*1e-11;

        // Newton method to do the remaining work. If the initial
        // value is good, this should only take two to three
        // iterations...
        for (int i = 0; i < 20; ++i) {
            const Evaluation& f = saturatedGasDissolutionFactor(regionIdx, temperature, pSat) - Rs;
            const Evaluation& fPrime = ((saturatedGasDissolutionFactor(regionIdx, temperature, pSat + eps) - Rs) - f)/eps;

            const Evaluation& delta = f/fPrime;
            pSat -= delta;

            Scalar absDelta = std::abs(Toolbox::value(delta));
            if (absDelta < Toolbox::value(pSat) * 1e-10 || absDelta < 1e-4)
                return pSat;
        }

        OPM_THROW(NumericalProblem, "Could find the oil saturation pressure for X_o^G = " << Rs);
    }

private:
    void updateSaturationPressureSpline_(unsigned regionIdx)
    {
        auto& gasDissolutionFac = saturatedGasDissolutionFactorTable_[regionIdx];

        // create the spline representing saturation pressure
        // depending of the mass fraction in gas
        size_t n = gasDissolutionFac.numSamples()*5;
        Scalar delta = (gasDissolutionFac.xMax() - gasDissolutionFac.xMin())/(n + 1);

        SamplingPoints pSatSamplePoints;
        Scalar Rs = 0;
        for (size_t i=0; i <= n; ++ i) {
            Scalar pSat = gasDissolutionFac.xMin() + i*delta;
            Rs = saturatedGasDissolutionFactor(regionIdx,
                                               /*temperature=*/Scalar(1e100),
                                               pSat);

            std::pair<Scalar, Scalar> val(Rs, pSat);
            pSatSamplePoints.push_back(val);
        }
        saturationPressureSpline_[regionIdx].setContainerOfTuples(pSatSamplePoints,
                                                                  /*type=*/Spline::Monotonic);
    }

    std::vector<Scalar> gasReferenceDensity_;
    std::vector<Scalar> oilReferenceDensity_;
    std::vector<TabulatedTwoDFunction> inverseOilBTable_;
    std::vector<TabulatedTwoDFunction> oilMuTable_;
    std::vector<TabulatedTwoDFunction> inverseOilBMuTable_;
    std::vector<TabulatedOneDFunction> saturatedOilMuTable_;
    std::vector<TabulatedOneDFunction> inverseSaturatedOilBTable_;
    std::vector<TabulatedOneDFunction> inverseSaturatedOilBMuTable_;
    std::vector<TabulatedOneDFunction> saturatedGasDissolutionFactorTable_;
    std::vector<Spline> saturationPressureSpline_;
};

} // namespace Opm

#endif
