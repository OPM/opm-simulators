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
 * \copydoc Opm::WetGasPvt
 */
#ifndef OPM_WET_GAS_PVT_HPP
#define OPM_WET_GAS_PVT_HPP

#include <opm/material/Constants.hpp>

#include <opm/material/common/OpmFinal.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/Spline.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

namespace Opm {

template <class Scalar>
class OilPvtMultiplexer;

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phas
 *        with vaporized oil.
 */
template <class Scalar>
class WetGasPvt
{
    typedef Opm::OilPvtMultiplexer<Scalar> OilPvtMultiplexer;

    typedef Opm::UniformXTabulated2DFunction<Scalar> TabulatedTwoDFunction;
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef Opm::Spline<Scalar> Spline;
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

public:
#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize the parameters for wet gas using an ECL deck.
     *
     * This method assumes that the deck features valid DENSITY and PVTG keywords.
     */
    void initFromDeck(DeckConstPtr deck, EclipseStateConstPtr eclState)
    {
        const auto& pvtgTables = eclState->getTableManager()->getPvtgTables();
        DeckKeywordConstPtr densityKeyword = deck->getKeyword("DENSITY");

        assert(pvtgTables.size() == densityKeyword->size());

        size_t numRegions = pvtgTables.size();
        setNumRegions(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            Scalar rhoRefO = densityKeyword->getRecord(regionIdx)->getItem("OIL")->getSIDouble(0);
            Scalar rhoRefG = densityKeyword->getRecord(regionIdx)->getItem("GAS")->getSIDouble(0);
            Scalar rhoRefW = densityKeyword->getRecord(regionIdx)->getItem("WATER")->getSIDouble(0);

            setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);

            // determine the molar masses of the components
            Scalar p = 1.01325e5; // surface pressure, [Pa]
            Scalar T = 273.15 + 15.56; // surface temperature, [K]
            Scalar MO = 175e-3; // [kg/mol]
            Scalar MG = Opm::Constants<Scalar>::R*T*rhoRefG / p; // [kg/mol], consequence of the ideal gas law
            Scalar MW = 18.0e-3; // [kg/mol]
            // TODO (?): the molar mass of the components can possibly specified
            // explicitly in the deck.
            setMolarMasses(regionIdx, MO, MG, MW);
        }

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            const auto& pvtgTable = pvtgTables[regionIdx];

            const auto saturatedTable = pvtgTable.getOuterTable();
            assert(saturatedTable->numRows() > 1);

            auto& gasMu = gasMu_[regionIdx];
            auto& invGasB = inverseGasB_[regionIdx];
            auto& invSatGasB = inverseSaturatedGasB_[regionIdx];
            auto& invSatGasBMu = inverseSaturatedGasBMu_[regionIdx];
            auto& oilVaporizationFac = oilVaporizationFactorTable_[regionIdx];

            oilVaporizationFac.setXYArrays(saturatedTable->numRows(),
                                           saturatedTable->getPressureColumn(),
                                           saturatedTable->getOilSolubilityColumn());

            std::vector<Scalar> invSatGasBArray;
            std::vector<Scalar> invSatGasBMuArray;

            // extract the table for the gas dissolution and the oil formation volume factors
            for (size_t outerIdx = 0; outerIdx < saturatedTable->numRows(); ++ outerIdx) {
                Scalar pg = saturatedTable->getPressureColumn()[outerIdx];
                Scalar B = saturatedTable->getGasFormationFactorColumn()[outerIdx];
                Scalar mu = saturatedTable->getGasViscosityColumn()[outerIdx];

                invGasB.appendXPos(pg);
                gasMu.appendXPos(pg);

                invSatGasBArray.push_back(1.0/B);
                invSatGasBMuArray.push_back(1.0/(mu*B));

                assert(invGasB.numX() == outerIdx + 1);
                assert(gasMu.numX() == outerIdx + 1);

                const auto underSaturatedTable = pvtgTable.getInnerTable(outerIdx);
                size_t numRows = underSaturatedTable->numRows();
                for (size_t innerIdx = 0; innerIdx < numRows; ++ innerIdx) {
                    Scalar Rv = underSaturatedTable->getOilSolubilityColumn()[innerIdx];
                    Scalar Bg = underSaturatedTable->getGasFormationFactorColumn()[innerIdx];
                    Scalar mug = underSaturatedTable->getGasViscosityColumn()[innerIdx];

                    invGasB.appendSamplePoint(outerIdx, Rv, 1.0/Bg);
                    gasMu.appendSamplePoint(outerIdx, Rv, mug);
                }
            }

            invSatGasB.setXYContainers(saturatedTable->getPressureColumn(), invSatGasBArray);
            invSatGasBMu.setXYContainers(saturatedTable->getPressureColumn(), invSatGasBMuArray);

            // make sure to have at least two sample points per gas pressure value
            for (size_t xIdx = 0; xIdx < invGasB.numX(); ++xIdx) {
               // a single sample point is definitely needed
                assert(invGasB.numY(xIdx) > 0);

                // everything is fine if the current table has two or more sampling points
                // for a given mole fraction
                if (invGasB.numY(xIdx) > 1)
                    continue;

                // find the master table which will be used as a template to extend the
                // current line. We define master table as the first table which has values
                // for undersaturated gas...
                size_t masterTableIdx = xIdx + 1;
                for (; masterTableIdx < pvtgTable.getOuterTable()->numRows(); ++masterTableIdx)
                {
                    if (pvtgTable.getInnerTable(masterTableIdx)->numRows() > 1)
                        break;
                }

                if (masterTableIdx >= pvtgTable.getOuterTable()->numRows())
                    OPM_THROW(std::runtime_error,
                              "PVTG tables are invalid: The last table must exhibit at least one "
                              "entry for undersaturated gas!");

                // extend the current table using the master table.
                extendPvtgTable_(regionIdx,
                                 xIdx,
                                 *pvtgTable.getInnerTable(xIdx),
                                 *pvtgTable.getInnerTable(masterTableIdx));
            }
        }
    }

private:
    void extendPvtgTable_(int regionIdx,
                          int xIdx,
                          const PvtgInnerTable& curTable,
                          const PvtgInnerTable& masterTable)
    {
        std::vector<Scalar> RvArray(curTable.getOilSolubilityColumn());
        std::vector<Scalar> gasBArray(curTable.getGasFormationFactorColumn());
        std::vector<Scalar> gasMuArray(curTable.getGasViscosityColumn());

        auto& invGasB = inverseGasB_[regionIdx];
        auto& gasMu = gasMu_[regionIdx];

        for (unsigned newRowIdx = 1; newRowIdx < masterTable.numRows(); ++ newRowIdx) {
            // compute the gas pressure for the new entry
            Scalar diffRv =
                masterTable.getOilSolubilityColumn()[newRowIdx]
                - masterTable.getOilSolubilityColumn()[newRowIdx - 1];
            Scalar newRv = RvArray.back() + diffRv;

            // calculate the compressibility of the master table
            Scalar B1 = masterTable.getGasFormationFactorColumn()[newRowIdx];
            Scalar B2 = masterTable.getGasFormationFactorColumn()[newRowIdx - 1];
            Scalar x = (B1 - B2)/( (B1 + B2)/2.0 );

            // calculate the gas formation volume factor which exhibits the same
            // "compressibility" for the new value of Rv
            Scalar newBg = gasBArray.back()*(1.0 + x/2.0)/(1.0 - x/2.0);

            // calculate the "viscosibility" of the master table
            Scalar mu1 = masterTable.getGasViscosityColumn()[newRowIdx];
            Scalar mu2 = masterTable.getGasViscosityColumn()[newRowIdx - 1];
            Scalar xMu = (mu1 - mu2)/( (mu1 + mu2)/2.0 );

            // calculate the gas formation volume factor which exhibits the same
            // compressibility for the new pressure
            Scalar newMug = gasMuArray.back()*(1.0 + xMu/2)/(1.0 - xMu/2.0);

            // append the new values to the arrays which we use to compute the additional
            // values ...
            RvArray.push_back(newRv);
            gasBArray.push_back(newBg);
            gasMuArray.push_back(newMug);

            // ... and register them with the internal table objects
            invGasB.appendSamplePoint(xIdx, newRv, 1.0/newBg);
            gasMu.appendSamplePoint(xIdx, newRv, newMug);
        }
    }

public:
#endif // HAVE_OPM_PARSER

    void setNumRegions(size_t numRegions)
    {
        oilMolarMass_.resize(numRegions);
        gasMolarMass_.resize(numRegions);
        oilReferenceDensity_.resize(numRegions);
        gasReferenceDensity_.resize(numRegions);
        inverseGasB_.resize(numRegions);
        inverseGasBMu_.resize(numRegions);
        inverseSaturatedGasB_.resize(numRegions);
        inverseSaturatedGasBMu_.resize(numRegions);
        gasMu_.resize(numRegions);
        oilVaporizationFactorTable_.resize(numRegions);
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
     * \brief Initialize the reference densities of all fluids for a given PVT region
     */
    void setMolarMasses(unsigned regionIdx,
                        Scalar MOil,
                        Scalar MGas,
                        Scalar /*MWater*/)
    {
        oilMolarMass_[regionIdx] = MOil;
        gasMolarMass_[regionIdx] = MGas;
    }

    /*!
     * \brief Initialize the function for the oil vaporization factor \f$R_v\f$
     *
     * \param samplePoints A container of (x,y) values.
     */
    void setSaturatedGasOilVaporizationFactor(unsigned regionIdx, const SamplingPoints &samplePoints)
    { oilVaporizationFactorTable_[regionIdx].setContainerOfTuples(samplePoints); }

    /*!
     * \brief Initialize the function for the gas formation volume factor
     *
     * The gas formation volume factor \f$B_g\f$ is a function of \f$(p_g, X_g^O)\f$ and
     * represents the partial density of the oil component in the gas phase at a given
     * pressure. This method only requires the volume factor of oil-saturated gas (which
     * only depends on pressure) while the dependence on the oil mass fraction is
     * guesstimated...
     */
    void setSaturatedGasFormationVolumeFactor(unsigned regionIdx, const SamplingPoints &samplePoints)
    {
        auto& invGasB = inverseGasB_[regionIdx];

        auto &RvTable = oilVaporizationFactorTable_[regionIdx];

        Scalar T = 273.15 + 15.56; // [K]

        Scalar RvMin = 0.0;
        Scalar RvMax = RvTable.eval(oilVaporizationFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

        Scalar poMin = samplePoints.front().first;
        Scalar poMax = samplePoints.back().first;

        size_t nRv = 20;
        size_t nP = samplePoints.size()*2;

        Scalar rhogRef = gasReferenceDensity_[regionIdx];
        Scalar rhooRef = oilReferenceDensity_[regionIdx];

        Spline gasFormationVolumeFactorSpline;
        gasFormationVolumeFactorSpline.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);

        updateSaturationPressureSpline_(regionIdx);

        // calculate a table of estimated densities depending on pressure and gas mass
        // fraction. note that this assumes oil of constant compressibility. (having said
        // that, if only the saturated gas densities are available, there's not much
        // choice.)
        for (size_t RvIdx = 0; RvIdx < nRv; ++RvIdx) {
            Scalar Rv = RvMin + (RvMax - RvMin)*RvIdx/nRv;

            invGasB.appendXPos(Rv);

            for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
                Scalar pg = poMin + (poMax - poMin)*pIdx/nP;

                Scalar poSat = gasSaturationPressure(regionIdx, T, Rv);
                Scalar BgSat = gasFormationVolumeFactorSpline.eval(poSat, /*extrapolate=*/true);
                Scalar drhoo_dp = (1.1200 - 1.1189)/((5000 - 4000)*6894.76);
                Scalar rhoo = rhooRef/BgSat*(1 + drhoo_dp*(pg - poSat));

                Scalar Bg = rhooRef/rhoo;

                invGasB.appendSamplePoint(RvIdx, pg, 1.0/Bg);
            }
        }
    }

    /*!
     * \brief Initialize the function for the gas formation volume factor
     *
     * The gas formation volume factor \f$B_g\f$ is a function of \f$(p_g, X_g^O)\f$ and
     * represents the partial density of the oil component in the gas phase at a given
     * pressure.
     *
     * This method sets \f$1/B_g(R_v, p_g)\f$. Note that instead of the mass fraction of
     * the oil component in the gas phase, this function depends on the gas dissolution
     * factor. Also note, that the order of the arguments needs to be \f$(R_s, p_o)\f$
     * and not the other way around.
     */
    void setInverseGasFormationVolumeFactor(unsigned regionIdx, const TabulatedTwoDFunction& invBg)
    { inverseGasB_[regionIdx] = invBg; }

    /*!
     * \brief Initialize the viscosity of the gas phase.
     *
     * This is a function of \f$(R_s, p_o)\f$...
     */
    void setGasViscosity(unsigned regionIdx, const TabulatedTwoDFunction& mug)
    { gasMu_[regionIdx] = mug; }

    /*!
     * \brief Initialize the phase viscosity for oil saturated gas
     *
     * The gas viscosity is a function of \f$(p_g, X_g^O)\f$, but this method only
     * requires the viscosity of oil-saturated gas (which only depends on pressure) while
     * there is assumed to be no dependence on the gas mass fraction...
     */
    void setSaturatedGasViscosity(unsigned regionIdx, const SamplingPoints &samplePoints  )
    {
        auto& oilVaporizationFac = oilVaporizationFactorTable_[regionIdx];

        Scalar RvMin = 0.0;
        Scalar RvMax = oilVaporizationFac.eval(oilVaporizationFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

        Scalar poMin = samplePoints.front().first;
        Scalar poMax = samplePoints.back().first;

        size_t nRv = 20;
        size_t nP = samplePoints.size()*2;

        Spline mugSpline;
        mugSpline.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);

        // calculate a table of estimated densities depending on pressure and gas mass
        // fraction
        for (size_t RvIdx = 0; RvIdx < nRv; ++RvIdx) {
            Scalar Rv = RvMin + (RvMax - RvMin)*RvIdx/nRv;

            gasMu_[regionIdx].appendXPos(Rv);

            for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
                Scalar pg = poMin + (poMax - poMin)*pIdx/nP;
                Scalar mug = mugSpline.eval(pg, /*extrapolate=*/true);

                gasMu_[regionIdx].appendSamplePoint(RvIdx, pg, mug);
            }
        }
    }

    /*!
     * \brief Finish initializing the gas phase PVT properties.
     */
    void initEnd(const OilPvtMultiplexer *oilPvt)
    {
        oilPvt_ = oilPvt;

        // calculate the final 2D functions which are used for interpolation.
        size_t numRegions = gasMu_.size();
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            // calculate the table which stores the inverse of the product of the gas
            // formation volume factor and the gas viscosity
            const auto& gasMu = gasMu_[regionIdx];
            const auto& invGasB = inverseGasB_[regionIdx];
            assert(gasMu.numX() == invGasB.numX());

            auto& invGasBMu = inverseGasBMu_[regionIdx];
            auto& invSatGasB = inverseSaturatedGasB_[regionIdx];
            auto& invSatGasBMu = inverseSaturatedGasBMu_[regionIdx];

            std::vector<Scalar> satPressuresArray;
            std::vector<Scalar> invSatGasBArray;
            std::vector<Scalar> invSatGasBMuArray;
            for (size_t pIdx = 0; pIdx < gasMu.numX(); ++pIdx) {
                invGasBMu.appendXPos(gasMu.xAt(pIdx));

                assert(gasMu.numY(pIdx) == invGasB.numY(pIdx));

                size_t numRv = gasMu.numY(pIdx);
                for (size_t rvIdx = 0; rvIdx < numRv; ++rvIdx)
                    invGasBMu.appendSamplePoint(pIdx,
                                                gasMu.yAt(pIdx, rvIdx),
                                                invGasB.valueAt(pIdx, rvIdx)
                                                / gasMu.valueAt(pIdx, rvIdx));

                // the sampling points in UniformXTabulated2DFunction are always sorted
                // in ascending order. Thus, the value for saturated gas is the last one
                // (i.e., the one with the largest Rv value)
                satPressuresArray.push_back(gasMu.xAt(pIdx));
                invSatGasBArray.push_back(invGasB.valueAt(pIdx, numRv - 1));
                invSatGasBMuArray.push_back(invGasBMu.valueAt(pIdx, numRv - 1));
            }

            invSatGasB.setXYContainers(satPressuresArray, invSatGasBArray);
            invSatGasBMu.setXYContainers(satPressuresArray, invSatGasBMuArray);

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
                         const Evaluation& Rv) const
    {
        const Evaluation& invBg = inverseGasB_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);
        const Evaluation& invMugBg = inverseGasBMu_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);

        return invBg/invMugBg;
    }

    /*!
     * \brief Returns the dynamic viscosity [Pa s] of oil saturated gas at a given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedViscosity(unsigned regionIdx,
                                  const Evaluation& /*temperature*/,
                                  const Evaluation& pressure) const
    {
        const Evaluation& invBg = inverseSaturatedGasB_[regionIdx].eval(pressure, /*extrapolate=*/true);
        const Evaluation& invMugBg = inverseSaturatedGasBMu_[regionIdx].eval(pressure, /*extrapolate=*/true);

        return invBg/invMugBg;
    }

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    template <class Evaluation>
    Evaluation density(unsigned regionIdx,
                       const Evaluation& temperature,
                       const Evaluation& pressure,
                       const Evaluation& Rv) const
    {
        const Evaluation& Bg = formationVolumeFactor(regionIdx, temperature, pressure, Rv);

        Scalar rhooRef = oilReferenceDensity_[regionIdx];
        Scalar rhogRef = gasReferenceDensity_[regionIdx];
        Evaluation rhog = rhogRef/Bg;

        // the oil formation volume factor just represents the partial density of the gas
        // component in the gas phase. to get the total density of the phase, we have to
        // add the partial density of the oil component.
        rhog += (rhooRef*Rv)/Bg;

        return rhog;
    }

    /*!
     * \brief Returns the density [kg/m^3] of oil saturated gas at a given pressure.
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

        const Evaluation& Bg = saturatedFormationVolumeFactor(regionIdx, temperature, pressure);

        Evaluation rhog = rhogRef/Bg;
        const Evaluation& RvSat = oilVaporizationFactor(regionIdx, temperature, pressure);
        rhog += (rhooRef*RvSat)/Bg;

        return rhog;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation formationVolumeFactor(unsigned regionIdx,
                                     const Evaluation& /*temperature*/,
                                     const Evaluation& pressure,
                                     const Evaluation& Rv) const
    { return 1.0 / inverseGasB_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true); }

    /*!
     * \brief Returns the formation volume factor [-] of oil saturated gas at a given pressure.
     */
    template <class Evaluation>
    Evaluation saturatedFormationVolumeFactor(unsigned regionIdx,
                                              const Evaluation& /*temperature*/,
                                              const Evaluation& pressure) const
    { return 1.0 / inverseSaturatedGasB_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the fugacity coefficient [Pa] of a component in the fluid phase given
     *        a set of parameters.
     */
    template <class Evaluation>
    Evaluation fugacityCoefficientGas(unsigned /*regionIdx*/,
                                      const Evaluation& /*temperature*/,
                                      const Evaluation& /*pressure*/) const
    {
        // the fugacity coefficient of the gas component in the gas phase is assumed to
        // be that of an ideal gas.
        return 1.0;
    }

    template <class Evaluation>
    Evaluation fugacityCoefficientOil(unsigned regionIdx,
                                      const Evaluation& temperature,
                                      const Evaluation& pressure) const
    {
        // the fugacity coefficient of the oil component in the wet gas phase:
        //
        // first, retrieve the mole fraction of gas a saturated oil would exhibit at the
        // given pressure
        const Evaluation& x_gOSat = saturatedOilMoleFraction(regionIdx, temperature, pressure);

        // then, scale the oil component's gas phase fugacity coefficient, so that the
        // gas phase ends up at the right composition if we were doing a flash experiment
        const Evaluation& phi_oO = oilPvt_->fugacityCoefficientOil(regionIdx, temperature, pressure);

        return phi_oO / x_gOSat;
    }

    template <class Evaluation>
    Evaluation fugacityCoefficientWater(unsigned regionIdx,
                                        const Evaluation& temperature,
                                        const Evaluation& pressure) const
    {
        // assume that the affinity of the water component to the gas phase is much
        // smaller than that of the gas component
        return 1e8*fugacityCoefficientWater(regionIdx, temperature, pressure);
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation oilVaporizationFactor(unsigned regionIdx,
                                     const Evaluation& /*temperature*/,
                                     const Evaluation& pressure) const
    { return oilVaporizationFactorTable_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *        depending on its mass fraction of the oil component
     *
     * \param Rv The surface volume of oil component dissolved in what will yield one cubic meter of gas at the surface [-]
     */
    template <class Evaluation>
    Evaluation gasSaturationPressure(unsigned regionIdx,
                                     const Evaluation& temperature,
                                     const Evaluation& Rv) const
    {
        typedef Opm::MathToolbox<Evaluation> Toolbox;

        // use the saturation pressure spline to get a pretty good initial value
        Evaluation pSat = saturationPressureSpline_[regionIdx].eval(Rv, /*extrapolate=*/true);
        const Evaluation& eps = pSat*1e-11;

        // Newton method to do the remaining work. If the initial
        // value is good, this should only take two to three
        // iterations...
        for (unsigned i = 0; i < 20; ++i) {
            const Evaluation& f = oilVaporizationFactor(regionIdx, temperature, pSat) - Rv;
            const Evaluation& fPrime = ((oilVaporizationFactor(regionIdx, temperature, pSat + eps) - Rv) - f)/eps;

            const Evaluation& delta = f/fPrime;
            pSat -= delta;

            if (std::abs(Toolbox::value(delta)) < std::abs(Toolbox::value(pSat)) * 1e-10)
                return pSat;
        }

        OPM_THROW(NumericalProblem, "Could find the gas saturation pressure for X_g^O = " << Rv);
    }

    template <class Evaluation>
    Evaluation saturatedOilMassFraction(unsigned regionIdx,
                                        const Evaluation& temperature,
                                        const Evaluation& pressure) const
    {
        Scalar rho_gRef = gasReferenceDensity_[regionIdx];
        Scalar rho_oRef = oilReferenceDensity_[regionIdx];

        // calculate the mass of the oil component [kg/m^3] in the gas phase. This is
        // equivalent to the oil vaporization factor [m^3/m^3] at current pressure times
        // the oil density [kg/m^3] at standard pressure
        const Evaluation& Rv = oilVaporizationFactor(regionIdx, temperature, pressure);
        const Evaluation& rho_gO = Rv * rho_oRef;

        // we now have the total density of saturated oil and the partial density of the
        // oil component within it. The gas mass fraction is the ratio of these two.
        return rho_gO/(rho_gRef + rho_gO);
    }

    template <class Evaluation>
    Evaluation saturatedOilMoleFraction(unsigned regionIdx,
                                        const Evaluation& temperature,
                                        const Evaluation& pressure) const
    {
        // calculate the mass fractions of gas and oil
        const Evaluation& Rv = saturatedOilMassFraction(regionIdx, temperature, pressure);

        // which can be converted to mole fractions, given the
        // components' molar masses
        Scalar MG = gasMolarMass_[regionIdx];
        Scalar MO = oilMolarMass_[regionIdx];

        const Evaluation& avgMolarMass = MO/(1 + (1 - Rv)*(MO/MG - 1));
        return Rv*avgMolarMass/MO;
    }

private:
    void updateSaturationPressureSpline_(unsigned regionIdx)
    {
        auto& oilVaporizationFac = oilVaporizationFactorTable_[regionIdx];

        // create the spline representing saturation pressure
        // depending of the mass fraction in gas
        size_t n = oilVaporizationFac.numSamples()*5;
        Scalar delta = (oilVaporizationFac.xMax() - oilVaporizationFac.xMin())/(n + 1);

        SamplingPoints pSatSamplePoints;
        Scalar Rv = 0;
        for (size_t i = 0; i <= n; ++ i) {
            Scalar pSat = oilVaporizationFac.xMin() + i*delta;
            Rv = oilVaporizationFactor(regionIdx, /*temperature=*/Scalar(1e100), pSat);

            std::pair<Scalar, Scalar> val(Rv, pSat);
            pSatSamplePoints.push_back(val);
        }
        saturationPressureSpline_[regionIdx].setContainerOfTuples(pSatSamplePoints,
                                                                  /*type=*/Spline::Monotonic);
    }

    const OilPvtMultiplexer *oilPvt_;

    std::vector<Scalar> gasMolarMass_;
    std::vector<Scalar> oilMolarMass_;
    std::vector<Scalar> gasReferenceDensity_;
    std::vector<Scalar> oilReferenceDensity_;
    std::vector<TabulatedTwoDFunction> inverseGasB_;
    std::vector<TabulatedOneDFunction> inverseSaturatedGasB_;
    std::vector<TabulatedTwoDFunction> gasMu_;
    std::vector<TabulatedTwoDFunction> inverseGasBMu_;
    std::vector<TabulatedOneDFunction> inverseSaturatedGasBMu_;
    std::vector<TabulatedOneDFunction> oilVaporizationFactorTable_;
    std::vector<Spline> saturationPressureSpline_;
};

} // namespace Opm

#endif
