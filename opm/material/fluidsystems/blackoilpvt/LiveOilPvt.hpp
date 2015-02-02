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

#include "OilPvtInterface.hpp"

#include <opm/material/OpmFinal.hpp>
#include <opm/material/UniformXTabulated2DFunction.hpp>
#include <opm/material/Tabulated1DFunction.hpp>
#include <opm/core/utility/Spline.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

namespace Opm {

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil phas
 *        with dissolved gas.
 */
template <class Scalar>
class LiveOilPvt : public OilPvtInterface<Scalar>
{
    typedef FluidSystems::BlackOil<Scalar> BlackOilFluidSystem;

    typedef Opm::UniformXTabulated2DFunction<Scalar> TabulatedTwoDFunction;
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef Opm::Spline<Scalar> Spline;
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
        if (static_cast<int>(inverseOilB_.size()) < numRegions) {
            inverseOilB_.resize(numRegions);
            inverseOilBMu_.resize(numRegions);
            oilMu_.resize(numRegions);
            gasDissolutionFactor_.resize(numRegions);
            saturationPressureSpline_.resize(numRegions);
        }
    }

    /*!
     * \brief Initialize the function for the gas dissolution factor \f$R_s\f$
     *
     * \param samplePoints A container of (x,y) values.
     */
    void setSaturatedOilGasDissolutionFactor(int regionIdx, const SamplingPoints &samplePoints)
    { gasDissolutionFactor_[regionIdx].setContainerOfTuples(samplePoints); }

    /*!
     * \brief Initialize the function for the oil formation volume factor
     *
     * The oil formation volume factor \f$B_o\f$ is a function of \f$(p_o, X_o^G)\f$ and
     * represents the partial density of the oil component in the oil phase at a given
     * pressure. This method only requires the volume factor of gas-saturated oil (which
     * only depends on pressure) while the dependence on the gas mass fraction is
     * guesstimated...
     */
    void setSaturatedOilFormationVolumeFactor(int regionIdx, const SamplingPoints &samplePoints)
    {
        auto& invOilB = inverseOilB_[regionIdx];

        auto &Rs = gasDissolutionFactor_[regionIdx];

        Scalar T = BlackOilFluidSystem::surfaceTemperature;

        Scalar RsMin = 0.0;
        Scalar RsMax = Rs.eval(gasDissolutionFactor_[regionIdx].xMax(), /*extrapolate=*/true);

        Scalar poMin = samplePoints.front().first;
        Scalar poMax = samplePoints.back().first;

        size_t nRs = 20;
        size_t nP = samplePoints.size()*2;

        Scalar rhogRef = BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx);
        Scalar rhooRef = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx);

        Spline oilFormationVolumeFactorSpline;
        oilFormationVolumeFactorSpline.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);

        updateSaturationPressureSpline_(regionIdx);

        // calculate a table of estimated densities depending on pressure and gas mass
        // fraction
        for (size_t RsIdx = 0; RsIdx < nRs; ++RsIdx) {
            Scalar Rs = RsMin + (RsMax - RsMin)*RsIdx/nRs;
            Scalar XoG = Rs/(rhooRef/rhogRef + Rs);

            invOilB.appendXPos(Rs);

            for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
                Scalar po = poMin + (poMax - poMin)*pIdx/nP;

                Scalar poSat = oilSaturationPressure(regionIdx, T, XoG);
                Scalar BoSat = oilFormationVolumeFactorSpline.eval(poSat, /*extrapolate=*/true);
                Scalar drhoo_dp = (1.1200 - 1.1189)/((5000 - 4000)*6894.76);
                Scalar rhoo = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx)/BoSat*(1 + drhoo_dp*(po - poSat));

                Scalar Bo = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx)/rhoo;

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
    void setInverseOilFormationVolumeFactor(int regionIdx, const TabulatedTwoDFunction& invBo)
    { inverseOilB_[regionIdx] = invBo; }

    /*!
     * \brief Initialize the viscosity of the oil phase.
     *
     * This is a function of \f$(R_s, p_o)\f$...
     */
    void setOilViscosity(int regionIdx, const TabulatedTwoDFunction& muo)
    { oilMu_[regionIdx] = muo; }

    /*!
     * \brief Initialize the phase viscosity for gas saturated oil
     *
     * The oil viscosity is a function of \f$(p_o, X_o^G)\f$, but this method only
     * requires the viscosity of gas-saturated oil (which only depends on pressure) while
     * there is assumed to be no dependence on the gas mass fraction...
     */
    void setSaturatedOilViscosity(int regionIdx, const SamplingPoints &samplePoints  )
    {
        auto& gasDissolutionFactor = gasDissolutionFactor_[regionIdx];

        Scalar RsMin = 0.0;
        Scalar RsMax = gasDissolutionFactor.eval(gasDissolutionFactor_[regionIdx].xMax(), /*extrapolate=*/true);

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

            oilMu_[regionIdx].appendXPos(Rs);

            for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
                Scalar po = poMin + (poMax - poMin)*pIdx/nP;
                Scalar muo = muoSpline.eval(po, /*extrapolate=*/true);

                oilMu_[regionIdx].appendSamplePoint(RsIdx, po, muo);
            }
        }
    }

#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize the oil parameters via the data specified by the PVTO ECL keyword.
     */
    void setPvtoTable(int regionIdx, const PvtoTable &pvtoTable)
    {
        const auto saturatedTable = pvtoTable.getOuterTable();
        assert(saturatedTable->numRows() > 1);

        auto& oilMu = oilMu_[regionIdx];
        auto& invOilB = inverseOilB_[regionIdx];
        auto& gasDissolutionFactor = gasDissolutionFactor_[regionIdx];

        gasDissolutionFactor.setXYArrays(saturatedTable->numRows(),
                                         saturatedTable->getPressureColumn(),
                                         saturatedTable->getGasSolubilityColumn());

        // extract the table for the gas dissolution and the oil formation volume factors
        for (int outerIdx = 0; outerIdx < static_cast<int>(saturatedTable->numRows()); ++ outerIdx) {
            Scalar Rs = saturatedTable->getGasSolubilityColumn()[outerIdx];

            invOilB.appendXPos(Rs);
            oilMu.appendXPos(Rs);

            assert(invOilB.numX() == outerIdx + 1);
            assert(oilMu.numX() == outerIdx + 1);

            const auto underSaturatedTable = pvtoTable.getInnerTable(outerIdx);
            int numRows = underSaturatedTable->numRows();
            for (int innerIdx = 0; innerIdx < numRows; ++ innerIdx) {
                Scalar po = underSaturatedTable->getPressureColumn()[innerIdx];
                Scalar Bo = underSaturatedTable->getOilFormationFactorColumn()[innerIdx];
                Scalar muo = underSaturatedTable->getOilViscosityColumn()[innerIdx];

                invOilB.appendSamplePoint(outerIdx, po, 1.0/Bo);
                oilMu.appendSamplePoint(outerIdx, po, muo);
            }
        }

        // make sure to have at least two sample points per mole fraction
        for (int xIdx = 0; xIdx < invOilB.numX(); ++xIdx) {
            // a single sample point is definitely needed
            assert(invOilB.numY(xIdx) > 0);

            // everything is fine if the current table has two or more sampling points
            // for a given mole fraction
            if (invOilB.numY(xIdx) > 1)
                continue;

            // find the master table which will be used as a template to extend the
            // current line. We define master table as the first table which has values
            // for undersaturated oil...
            int masterTableIdx = xIdx + 1;
            for (; masterTableIdx < static_cast<int>(pvtoTable.getOuterTable()->numRows());
                 ++masterTableIdx)
            {
                if (pvtoTable.getInnerTable(masterTableIdx)->numRows() > 1)
                    break;
            }

            if (masterTableIdx >= static_cast<int>(pvtoTable.getOuterTable()->numRows()))
                OPM_THROW(std::runtime_error,
                          "PVTO tables are invalid: The last table must exhibit at least one "
                          "entry for undersaturated oil!");

            // extend the current table using the master table. this is done by assuming
            // that the current table exhibits the same ratios of the oil formation
            // volume factors and viscosities for identical pressure rations as in the
            // master table.
            const auto masterTable = pvtoTable.getInnerTable(masterTableIdx);
            const auto curTable = pvtoTable.getInnerTable(xIdx);
            for (int newRowIdx = 1;
                 newRowIdx < static_cast<int>(masterTable->numRows());
                 ++ newRowIdx)
            {
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
#endif // HAVE_OPM_PARSER

    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd()
    {
        // calculate the final 2D functions which are used for interpolation.
        int numRegions = oilMu_.size();
        for (int regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            // calculate the table which stores the inverse of the product of the oil
            // formation volume factor and the oil viscosity
            const auto& oilMu = oilMu_[regionIdx];
            const auto& invOilB = inverseOilB_[regionIdx];
            assert(oilMu.numX() == invOilB.numX());

            auto& invOilBMu = inverseOilBMu_[regionIdx];

            for (int rsIdx = 0; rsIdx < oilMu.numX(); ++rsIdx) {
                invOilBMu.appendXPos(oilMu.xAt(rsIdx));

                assert(oilMu.numY(rsIdx) == invOilB.numY(rsIdx));

                int numPressures = oilMu.numY(rsIdx);
                for (int pIdx = 0; pIdx < numPressures; ++pIdx)
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
    Scalar viscosity(int regionIdx,
                     Scalar temperature,
                     Scalar pressure,
                     Scalar XoG) const OPM_FINAL
    {
        Scalar Rs =
            XoG/(1 - XoG)
            * BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx)
            / BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx);

        // ATTENTION: Rs is the first axis!
        Scalar invBo = inverseOilB_[regionIdx].eval(Rs, pressure, /*extrapolate=*/true);
        Scalar invMuoBo = inverseOilBMu_[regionIdx].eval(Rs, pressure, /*extrapolate=*/true);

        return invBo/invMuoBo;
    }

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    Scalar density(int regionIdx,
                   Scalar temperature,
                   Scalar pressure,
                   Scalar XoG) const OPM_FINAL
    {
        Scalar rhooRef = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx);
        Scalar rhogRef = BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx);

        Scalar Bo = formationVolumeFactor(regionIdx, temperature, pressure, XoG);
        Scalar rhoo = rhooRef/Bo;

        // the oil formation volume factor just represents the partial density of the oil
        // component in the oil phase. to get the total density of the phase, we have to
        // add the partial density of the gas component.
        Scalar Rs = XoG/(1 - XoG) * rhooRef/rhogRef;
        rhoo += rhogRef*Rs/Bo;

        return rhoo;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    Scalar formationVolumeFactor(int regionIdx,
                                 Scalar temperature,
                                 Scalar pressure,
                                 Scalar XoG) const OPM_FINAL
    {
        Scalar Rs =
            XoG/(1-XoG)
            *BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx)
            /BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx);

        // ATTENTION: Rs is represented by the _first_ axis!
        return 1.0 / inverseOilB_[regionIdx].eval(Rs, pressure, /*extrapolate=*/true);
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
        // set the oil component fugacity coefficient in oil phase
        // arbitrarily. we use some pseudo-realistic value for the vapor
        // pressure to ease physical interpretation of the results
        Scalar phi_oO = 20e3/pressure;

        if (compIdx == BlackOilFluidSystem::oilCompIdx)
            return phi_oO;
        else if (compIdx == BlackOilFluidSystem::waterCompIdx)
            // assume that the affinity of the water component to the
            // oil phase is one million times smaller than that of the
            // oil component
            return 1e8*phi_oO;

        /////////////
        // the rest of this method determines the fugacity coefficient
        // of the gas component:
        //
        // first, retrieve the mole fraction of gas a saturated oil
        // would exhibit at the given pressure
        Scalar x_oGSat = saturatedOilGasMoleFraction(regionIdx, temperature, pressure);

        // then, scale the gas component's gas phase fugacity
        // coefficient, so that the oil phase ends up at the right
        // composition if we were doing a flash experiment
        Scalar phi_gG = BlackOilFluidSystem::fugCoefficientInGas(gasCompIdx,
                                                                 temperature,
                                                                 pressure,
                                                                 regionIdx);


        return phi_gG / x_oGSat;
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    Scalar gasDissolutionFactor(int regionIdx,
                                Scalar temperature,
                                Scalar pressure) const OPM_FINAL
    { return gasDissolutionFactor_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param XoG The mass fraction of the gas component in the oil phase [-]
     */
    Scalar oilSaturationPressure(int regionIdx,
                                 Scalar temperature,
                                 Scalar XoG) const OPM_FINAL
    {
        // use the saturation pressure spline to get a pretty good initial value
        Scalar pSat = saturationPressureSpline_[regionIdx].eval(XoG, /*extrapolate=*/true);

        // Newton method to do the remaining work. If the initial
        // value is good, this should only take two to three
        // iterations...
        for (int i = 0; i < 20; ++i) {
            Scalar f = saturatedOilGasMassFraction(regionIdx, temperature, pSat) - XoG;
            Scalar eps = pSat*1e-11;
            Scalar fPrime = ((saturatedOilGasMassFraction(regionIdx, temperature, pSat + eps) - XoG) - f)/eps;

            Scalar delta = f/fPrime;
            pSat -= delta;

            if (std::abs(delta) < pSat * 1e-10)
                return pSat;
        }

        OPM_THROW(NumericalProblem, "Could find the oil saturation pressure for X_o^g = " << XoG);
    }

    Scalar saturatedOilGasMassFraction(int regionIdx,
                                       Scalar temperature,
                                       Scalar pressure) const OPM_FINAL
    {
        Scalar rho_gRef = BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx);
        Scalar rho_oRef = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx);

        // calculate the mass of the gas component [kg/m^3] in the oil phase. This is
        // equivalent to the gas dissolution factor [m^3/m^3] at current pressure times
        // the gas density [kg/m^3] at standard pressure
        Scalar rho_oG = gasDissolutionFactor(regionIdx, temperature, pressure) * rho_gRef;

        // we now have the total density of saturated oil and the partial density of the
        // gas component within it. The gas mass fraction is the ratio of these two.
        return rho_oG/(rho_oRef + rho_oG);
    }

    Scalar saturatedOilGasMoleFraction(int regionIdx,
                                       Scalar temperature,
                                       Scalar pressure) const OPM_FINAL
    {
        // calculate the mass fractions of gas and oil
        Scalar XoG = saturatedOilGasMassFraction(regionIdx, temperature, pressure);

        // which can be converted to mole fractions, given the
        // components' molar masses
        Scalar MG = BlackOilFluidSystem::molarMass(gasCompIdx, regionIdx);
        Scalar MO = BlackOilFluidSystem::molarMass(oilCompIdx, regionIdx);

        Scalar avgMolarMass = MO/(1 + XoG*(MO/MG - 1));
        Scalar xoG = XoG*avgMolarMass/MG;

        return xoG;
    }

private:
    void updateSaturationPressureSpline_(int regionIdx)
    {
        auto& gasDissolutionFactor = gasDissolutionFactor_[regionIdx];

        // create the spline representing saturation pressure
        // depending of the mass fraction in gas
        int n = gasDissolutionFactor.numSamples()*5;
        int delta = (gasDissolutionFactor.xMax() - gasDissolutionFactor.xMin())/(n + 1);

        SamplingPoints pSatSamplePoints;
        Scalar XoG = 0;
        for (int i=0; i <= n; ++ i) {
            Scalar pSat = gasDissolutionFactor.xMin() + i*delta;
            XoG = saturatedOilGasMassFraction(regionIdx, /*temperature=*/1e100, pSat);

            std::pair<Scalar, Scalar> val(XoG, pSat);
            pSatSamplePoints.push_back(val);
        }
        saturationPressureSpline_[regionIdx].setContainerOfTuples(pSatSamplePoints,
                                                                  /*type=*/Spline::Monotonic);
    }

    std::vector<TabulatedTwoDFunction> inverseOilB_;
    std::vector<TabulatedTwoDFunction> oilMu_;
    std::vector<TabulatedTwoDFunction> inverseOilBMu_;
    std::vector<TabulatedOneDFunction> gasDissolutionFactor_;
    std::vector<Spline> saturationPressureSpline_;
};

} // namespace Opm

#endif
