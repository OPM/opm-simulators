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
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

namespace Opm {
template <class Scalar>
class GasPvtMultiplexer;

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil phas
 *        with dissolved gas.
 */
template <class Scalar>
class LiveOilPvt
{
    typedef Opm::GasPvtMultiplexer<Scalar> GasPvtMultiplexer;

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

            // determine the molar masses of the components
            Scalar p = 1.01325e5; // surface pressure, [Pa]
            Scalar T = 273.15 + 15.56; // surface temperature, [K]
            Scalar MO = 175e-3; // [kg/mol]
            Scalar MG = Opm::Constants<Scalar>::R*T*rhoRefG / p; // [kg/mol], consequence of the ideal gas law
            Scalar MW = 18.0e-3; // [kg/mol]
            // TODO (?): the molar mass of the components can possibly specified
            // explicitly in the deck.
            setMolarMasses(regionIdx, MO, MG, MW);

            const auto& pvtoTable = pvtoTables[regionIdx];

            const auto saturatedTable = pvtoTable.getOuterTable();
            assert(saturatedTable->numRows() > 1);

            auto& oilMu = oilMuTable_[regionIdx];
            auto& invOilB = inverseOilBTable_[regionIdx];
            auto& gasDissolutionFac = gasDissolutionFactorTable_[regionIdx];

            gasDissolutionFac.setXYArrays(saturatedTable->numRows(),
                                             saturatedTable->getPressureColumn(),
                                             saturatedTable->getGasSolubilityColumn());

            // extract the table for the gas dissolution and the oil formation volume factors
            for (unsigned outerIdx = 0; outerIdx < saturatedTable->numRows(); ++ outerIdx) {
                Scalar Rs = saturatedTable->getGasSolubilityColumn()[outerIdx];

                invOilB.appendXPos(Rs);
                oilMu.appendXPos(Rs);

                assert(invOilB.numX() == outerIdx + 1);
                assert(oilMu.numX() == outerIdx + 1);

                const auto underSaturatedTable = pvtoTable.getInnerTable(outerIdx);
                size_t numRows = underSaturatedTable->numRows();
                for (unsigned innerIdx = 0; innerIdx < numRows; ++ innerIdx) {
                    Scalar po = underSaturatedTable->getPressureColumn()[innerIdx];
                    Scalar Bo = underSaturatedTable->getOilFormationFactorColumn()[innerIdx];
                    Scalar muo = underSaturatedTable->getOilViscosityColumn()[innerIdx];

                    invOilB.appendSamplePoint(outerIdx, po, 1.0/Bo);
                    oilMu.appendSamplePoint(outerIdx, po, muo);
                }
            }

            // make sure to have at least two sample points per mole fraction
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
                for (; masterTableIdx < pvtoTable.getOuterTable()->numRows(); ++masterTableIdx)
                {
                    if (pvtoTable.getInnerTable(masterTableIdx)->numRows() > 1)
                        break;
                }

                if (masterTableIdx >= pvtoTable.getOuterTable()->numRows())
                    OPM_THROW(std::runtime_error,
                              "PVTO tables are invalid: The last table must exhibit at least one "
                              "entry for undersaturated oil!");

                // extend the current table using the master table. this is done by assuming
                // that the current table exhibits the same ratios of the oil formation
                // volume factors and viscosities for identical pressure rations as in the
                // master table.
                const auto masterTable = pvtoTable.getInnerTable(masterTableIdx);
                const auto curTable = pvtoTable.getInnerTable(xIdx);
                for (unsigned newRowIdx = 1; newRowIdx < masterTable->numRows(); ++ newRowIdx) {
                    Scalar alphaPo =
                        masterTable->getPressureColumn()[newRowIdx]
                        / masterTable->getPressureColumn()[0];

                    Scalar alphaBo =
                        masterTable->getOilFormationFactorColumn()[newRowIdx]
                        / masterTable->getOilFormationFactorColumn()[0];

                    Scalar alphaMuo =
                        masterTable->getOilViscosityColumn()[newRowIdx]
                        / masterTable->getOilViscosityColumn()[0];

                    Scalar newPo = curTable->getPressureColumn()[0]*alphaPo;
                    Scalar newBo = curTable->getOilFormationFactorColumn()[0]*alphaBo;
                    Scalar newMuo = curTable->getOilViscosityColumn()[0]*alphaMuo;

                    invOilB.appendSamplePoint(xIdx, newPo, 1.0/newBo);
                    oilMu.appendSamplePoint(xIdx, newPo, newMuo);
                }
            }
        }
    }
#endif // HAVE_OPM_PARSER

    void setNumRegions(size_t numRegions)
    {
        oilMolarMass_.resize(numRegions);
        gasMolarMass_.resize(numRegions);
        oilReferenceDensity_.resize(numRegions);
        gasReferenceDensity_.resize(numRegions);
        inverseOilBTable_.resize(numRegions);
        inverseOilBMuTable_.resize(numRegions);
        oilMuTable_.resize(numRegions);
        gasDissolutionFactorTable_.resize(numRegions);
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
     * \brief Initialize the function for the gas dissolution factor \f$R_s\f$
     *
     * \param samplePoints A container of (x,y) values.
     */
    void setSaturatedOilGasDissolutionFactor(unsigned regionIdx, const SamplingPoints &samplePoints)
    { gasDissolutionFactorTable_[regionIdx].setContainerOfTuples(samplePoints); }

    /*!
     * \brief Initialize the function for the oil formation volume factor
     *
     * The oil formation volume factor \f$B_o\f$ is a function of \f$(p_o, X_o^G)\f$ and
     * represents the partial density of the oil component in the oil phase at a given
     * pressure. This method only requires the volume factor of gas-saturated oil (which
     * only depends on pressure) while the dependence on the gas mass fraction is
     * guesstimated...
     */
    void setSaturatedOilFormationVolumeFactor(unsigned regionIdx, const SamplingPoints &samplePoints)
    {
        auto& invOilB = inverseOilBTable_[regionIdx];

        auto &RsTable = gasDissolutionFactorTable_[regionIdx];

        Scalar T = 273.15 + 15.56; // [K]

        Scalar RsMin = 0.0;
        Scalar RsMax = RsTable.eval(gasDissolutionFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

        Scalar poMin = samplePoints.front().first;
        Scalar poMax = samplePoints.back().first;

        size_t nRs = 20;
        size_t nP = samplePoints.size()*2;

        Spline oilFormationVolumeFactorSpline;
        oilFormationVolumeFactorSpline.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);

        updateSaturationPressureSpline_(regionIdx);

        // calculate a table of estimated densities depending on pressure and gas mass
        // fraction
        for (size_t RsIdx = 0; RsIdx < nRs; ++RsIdx) {
            Scalar Rs = RsMin + (RsMax - RsMin)*RsIdx/nRs;

            invOilB.appendXPos(Rs);

            for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
                Scalar po = poMin + (poMax - poMin)*pIdx/nP;

                Scalar poSat = oilSaturationPressure(regionIdx, T, Rs);
                Scalar BoSat = oilFormationVolumeFactorSpline.eval(poSat, /*extrapolate=*/true);
                Scalar drhoo_dp = (1.1200 - 1.1189)/((5000 - 4000)*6894.76);
                Scalar rhoo = oilReferenceDensity_[regionIdx]/BoSat*(1 + drhoo_dp*(po - poSat));

                Scalar Bo = oilReferenceDensity_[regionIdx]/rhoo;

                invOilB.appendSamplePoint(RsIdx, po, 1.0/Bo);
            }
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
    void setSaturatedOilViscosity(unsigned regionIdx, const SamplingPoints &samplePoints  )
    {
        auto& gasDissolutionFac = gasDissolutionFactorTable_[regionIdx];

        Scalar RsMin = 0.0;
        Scalar RsMax = gasDissolutionFac.eval(gasDissolutionFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

        Scalar poMin = samplePoints.front().first;
        Scalar poMax = samplePoints.back().first;

        size_t nRs = 20;
        size_t nP = samplePoints.size()*2;

        Spline muoSpline;
        muoSpline.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);

        // calculate a table of estimated densities depending on pressure and gas mass
        // fraction
        for (size_t RsIdx = 0; RsIdx < nRs; ++RsIdx) {
            Scalar Rs = RsMin + (RsMax - RsMin)*RsIdx/nRs;

            oilMuTable_[regionIdx].appendXPos(Rs);

            for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
                Scalar po = poMin + (poMax - poMin)*pIdx/nP;
                Scalar muo = muoSpline.eval(po, /*extrapolate=*/true);

                oilMuTable_[regionIdx].appendSamplePoint(RsIdx, po, muo);
            }
        }
    }

    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd(const GasPvtMultiplexer *gasPvt)
    {
        gasPvt_ = gasPvt;

        // calculate the final 2D functions which are used for interpolation.
        size_t numRegions = oilMuTable_.size();
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            // calculate the table which stores the inverse of the product of the oil
            // formation volume factor and the oil viscosity
            const auto& oilMu = oilMuTable_[regionIdx];
            const auto& invOilB = inverseOilBTable_[regionIdx];
            assert(oilMu.numX() == invOilB.numX());

            auto& invOilBMu = inverseOilBMuTable_[regionIdx];

            for (unsigned rsIdx = 0; rsIdx < oilMu.numX(); ++rsIdx) {
                invOilBMu.appendXPos(oilMu.xAt(rsIdx));

                assert(oilMu.numY(rsIdx) == invOilB.numY(rsIdx));

                size_t numPressures = oilMu.numY(rsIdx);
                for (unsigned pIdx = 0; pIdx < numPressures; ++pIdx)
                    invOilBMu.appendSamplePoint(rsIdx,
                                                oilMu.yAt(rsIdx, pIdx),
                                                invOilB.valueAt(rsIdx, pIdx)*
                                                1/oilMu.valueAt(rsIdx, pIdx));
            }

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
     * \brief Returns the fugacity coefficient [Pa] of a component in the fluid phase given
     *        a set of parameters.
     */
    template <class Evaluation>
    Evaluation fugacityCoefficientOil(unsigned /*regionIdx*/,
                                      const Evaluation& /*temperature*/,
                                      const Evaluation& pressure) const
    {
        // set the oil component fugacity coefficient in oil phase
        // arbitrarily. we use some pseudo-realistic value for the vapor
        // pressure to ease physical interpretation of the results
        return 20e3/pressure;
    }

    template <class Evaluation>
    Evaluation fugacityCoefficientWater(unsigned regionIdx,
                                        const Evaluation& temperature,
                                        const Evaluation& pressure) const
    {
        // assume that the affinity of the water component to the
        // oil phase is one million times smaller than that of the
        // oil component
        return 1e8*fugacityCoefficientOil(regionIdx, temperature, pressure);
    }


    template <class Evaluation>
    Evaluation fugacityCoefficientGas(unsigned regionIdx,
                                      const Evaluation& temperature,
                                      const Evaluation& pressure) const
    {
        /////////////
        // the fugacity coefficient of the gas component:
        //
        // first, retrieve the mole fraction of gas a saturated oil
        // would exhibit at the given pressure
        const Evaluation& x_oGSat = saturatedOilGasMoleFraction(regionIdx, temperature, pressure);

        // then, scale the gas component's gas phase fugacity
        // coefficient, so that the oil phase ends up at the right
        // composition if we were doing a flash experiment
        const Evaluation& phi_gG = gasPvt_->fugacityCoefficientGas(regionIdx, temperature, pressure);

        return phi_gG / x_oGSat;
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    template <class Evaluation>
    Evaluation gasDissolutionFactor(unsigned regionIdx,
                                    const Evaluation& /*temperature*/,
                                    const Evaluation& pressure) const
    { return gasDissolutionFactorTable_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param Rs The surface volume of gas component dissolved in what will yield one cubic meter of oil at the surface [-]
     */
    template <class Evaluation>
    Evaluation oilSaturationPressure(unsigned regionIdx,
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
            const Evaluation& f = gasDissolutionFactor(regionIdx, temperature, pSat) - Rs;
            const Evaluation& fPrime = ((gasDissolutionFactor(regionIdx, temperature, pSat + eps) - Rs) - f)/eps;

            const Evaluation& delta = f/fPrime;
            pSat -= delta;

            Scalar absDelta = std::abs(Toolbox::value(delta));
            if (absDelta < Toolbox::value(pSat) * 1e-10 || absDelta < 1e-4)
                return pSat;
        }

        OPM_THROW(NumericalProblem, "Could find the oil saturation pressure for X_o^G = " << Rs);
    }

    template <class Evaluation>
    Evaluation saturatedOilGasMassFraction(unsigned regionIdx,
                                           const Evaluation& temperature,
                                           const Evaluation& pressure) const
    {
        Scalar rho_gRef = gasReferenceDensity_[regionIdx];
        Scalar rho_oRef = oilReferenceDensity_[regionIdx];

        // calculate the mass of the gas component [kg/m^3] in the oil phase. This is
        // equivalent to the gas dissolution factor [m^3/m^3] at current pressure times
        // the gas density [kg/m^3] at standard pressure
        const Evaluation& rho_oG = gasDissolutionFactor(regionIdx, temperature, pressure) * rho_gRef;

        // we now have the total density of saturated oil and the partial density of the
        // gas component within it. The gas mass fraction is the ratio of these two.
        return rho_oG/(rho_oRef + rho_oG);
    }

    template <class Evaluation>
    Evaluation saturatedOilGasMoleFraction(unsigned regionIdx,
                                           const Evaluation& temperature,
                                           const Evaluation& pressure) const
    {
        // calculate the mass fractions of gas and oil
        const Evaluation& XoG = saturatedOilGasMassFraction(regionIdx, temperature, pressure);
        Valgrind::CheckDefined(XoG);

        // which can be converted to mole fractions, given the
        // components' molar masses
        Scalar MG = gasMolarMass_[regionIdx];
        Scalar MO = oilMolarMass_[regionIdx];

        Evaluation avgMolarMass = MO/(1 + XoG*(MO/MG - 1));
        return XoG*avgMolarMass/MG;
    }

private:
    void updateSaturationPressureSpline_(unsigned regionIdx)
    {
        auto& gasDissolutionFac = gasDissolutionFactorTable_[regionIdx];

        // create the spline representing saturation pressure
        // depending of the mass fraction in gas
        size_t n = gasDissolutionFac.numSamples()*5;
        Scalar delta = (gasDissolutionFac.xMax() - gasDissolutionFac.xMin())/(n + 1);

        SamplingPoints pSatSamplePoints;
        Scalar Rs = 0;
        for (size_t i=0; i <= n; ++ i) {
            Scalar pSat = gasDissolutionFac.xMin() + i*delta;
            Rs = gasDissolutionFactor(regionIdx,
                                      /*temperature=*/Scalar(1e100),
                                      pSat);

            std::pair<Scalar, Scalar> val(Rs, pSat);
            pSatSamplePoints.push_back(val);
        }
        saturationPressureSpline_[regionIdx].setContainerOfTuples(pSatSamplePoints,
                                                                  /*type=*/Spline::Monotonic);
    }

    const GasPvtMultiplexer *gasPvt_;

    std::vector<Scalar> gasMolarMass_;
    std::vector<Scalar> oilMolarMass_;
    std::vector<Scalar> gasReferenceDensity_;
    std::vector<Scalar> oilReferenceDensity_;
    std::vector<TabulatedTwoDFunction> inverseOilBTable_;
    std::vector<TabulatedTwoDFunction> oilMuTable_;
    std::vector<TabulatedTwoDFunction> inverseOilBMuTable_;
    std::vector<TabulatedOneDFunction> gasDissolutionFactorTable_;
    std::vector<Spline> saturationPressureSpline_;
};

} // namespace Opm

#endif
