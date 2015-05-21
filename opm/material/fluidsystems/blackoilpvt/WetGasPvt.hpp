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

#include "GasPvtInterface.hpp"

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/material/common/OpmFinal.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/Spline.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phas
 *        with vaporized oil.
 */
template <class Scalar, class Evaluation = Scalar>
class WetGasPvt
    : public GasPvtInterfaceTemplateWrapper<Scalar, Evaluation, WetGasPvt<Scalar, Evaluation> >
{
    friend class GasPvtInterfaceTemplateWrapper<Scalar, Evaluation, WetGasPvt<Scalar, Evaluation> >;

    typedef FluidSystems::BlackOil<Scalar, Evaluation> BlackOilFluidSystem;

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
        inverseGasB_.resize(numRegions);
        inverseGasBMu_.resize(numRegions);
        gasMu_.resize(numRegions);
        oilVaporizationFactorTable_.resize(numRegions);
        saturationPressureSpline_.resize(numRegions);
    }

    /*!
     * \brief Initialize the function for the oil vaporization factor \f$R_v\f$
     *
     * \param samplePoints A container of (x,y) values.
     */
    void setSaturatedGasOilVaporizationFactor(int regionIdx, const SamplingPoints &samplePoints)
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
    void setSaturatedGasFormationVolumeFactor(int regionIdx, const SamplingPoints &samplePoints)
    {
        auto& invGasB = inverseGasB_[regionIdx];

        auto &Rv = oilVaporizationFactorTable_[regionIdx];

        Scalar T = BlackOilFluidSystem::surfaceTemperature;

        Scalar RvMin = 0.0;
        Scalar RvMax = Rv.eval(oilVaporizationFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

        Scalar poMin = samplePoints.front().first;
        Scalar poMax = samplePoints.back().first;

        size_t nRv = 20;
        size_t nP = samplePoints.size()*2;

        Scalar rhogRef = BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx);
        Scalar rhooRef = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx);

        Spline gasFormationVolumeFactorSpline;
        gasFormationVolumeFactorSpline.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);

        updateSaturationPressureSpline_(regionIdx);

        // calculate a table of estimated densities depending on pressure and gas mass
        // fraction
        for (size_t RvIdx = 0; RvIdx < nRv; ++RvIdx) {
            Scalar Rv = RvMin + (RvMax - RvMin)*RvIdx/nRv;
            Scalar XgO = Rv/(rhooRef/rhogRef + Rv);

            invGasB.appendXPos(Rv);

            for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
                Scalar pg = poMin + (poMax - poMin)*pIdx/nP;

                Scalar poSat = gasSaturationPressure(regionIdx, T, XgO);
                Scalar BgSat = gasFormationVolumeFactorSpline.eval(poSat, /*extrapolate=*/true);
                Scalar drhoo_dp = (1.1200 - 1.1189)/((5000 - 4000)*6894.76);
                Scalar rhoo = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx)/BgSat*(1 + drhoo_dp*(pg - poSat));

                Scalar Bg = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx)/rhoo;

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
    void setInverseGasFormationVolumeFactor(int regionIdx, const TabulatedTwoDFunction& invBg)
    { inverseGasB_[regionIdx] = invBg; }

    /*!
     * \brief Initialize the viscosity of the gas phase.
     *
     * This is a function of \f$(R_s, p_o)\f$...
     */
    void setGasViscosity(int regionIdx, const TabulatedTwoDFunction& mug)
    { gasMu_[regionIdx] = mug; }

    /*!
     * \brief Initialize the phase viscosity for oil saturated gas
     *
     * The gas viscosity is a function of \f$(p_g, X_g^O)\f$, but this method only
     * requires the viscosity of oil-saturated gas (which only depends on pressure) while
     * there is assumed to be no dependence on the gas mass fraction...
     */
    void setSaturatedGasViscosity(int regionIdx, const SamplingPoints &samplePoints  )
    {
        auto& oilVaporizationFactor = oilVaporizationFactorTable_[regionIdx];

        Scalar RvMin = 0.0;
        Scalar RvMax = oilVaporizationFactor.eval(oilVaporizationFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

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

#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize the oil parameters via the data specified by the PVTO ECL keyword.
     */
    void setPvtgTable(int regionIdx, const PvtgTable &pvtgTable)
    {
        const auto saturatedTable = pvtgTable.getOuterTable();
        assert(saturatedTable->numRows() > 1);

        auto& gasMu = gasMu_[regionIdx];
        auto& invGasB = inverseGasB_[regionIdx];
        auto& oilVaporizationFactor = oilVaporizationFactorTable_[regionIdx];

        oilVaporizationFactor.setXYArrays(saturatedTable->numRows(),
                                          saturatedTable->getPressureColumn(),
                                          saturatedTable->getOilSolubilityColumn());

        // extract the table for the gas dissolution and the oil formation volume factors
        for (int outerIdx = 0; outerIdx < static_cast<int>(saturatedTable->numRows()); ++ outerIdx) {
            Scalar pg = saturatedTable->getPressureColumn()[outerIdx];

            invGasB.appendXPos(pg);
            gasMu.appendXPos(pg);

            assert(invGasB.numX() == outerIdx + 1);
            assert(gasMu.numX() == outerIdx + 1);

            const auto underSaturatedTable = pvtgTable.getInnerTable(outerIdx);
            int numRows = underSaturatedTable->numRows();
            for (int innerIdx = 0; innerIdx < numRows; ++ innerIdx) {
                Scalar Rv = underSaturatedTable->getOilSolubilityColumn()[innerIdx];
                Scalar Bg = underSaturatedTable->getGasFormationFactorColumn()[innerIdx];
                Scalar mug = underSaturatedTable->getGasViscosityColumn()[innerIdx];

                invGasB.appendSamplePoint(outerIdx, Rv, 1.0/Bg);
                gasMu.appendSamplePoint(outerIdx, Rv, mug);
            }
        }

        // make sure to have at least two sample points per mole fraction
        for (int xIdx = 0; xIdx < invGasB.numX(); ++xIdx) {
            // a single sample point is definitely needed
            assert(invGasB.numY(xIdx) > 0);

            // everything is fine if the current table has two or more sampling points
            // for a given mole fraction
            if (invGasB.numY(xIdx) > 1)
                continue;

            // find the master table which will be used as a template to extend the
            // current line. We define master table as the first table which has values
            // for undersaturated gas...
            int masterTableIdx = xIdx + 1;
            for (; masterTableIdx < static_cast<int>(pvtgTable.getOuterTable()->numRows());
                 ++masterTableIdx)
            {
                if (pvtgTable.getInnerTable(masterTableIdx)->numRows() > 1)
                    break;
            }

            if (masterTableIdx >= static_cast<int>(pvtgTable.getOuterTable()->numRows()))
                OPM_THROW(std::runtime_error,
                          "PVTG tables are invalid: The last table must exhibit at least one "
                          "entry for undersaturated gas!");

            // extend the current table using the master table. this is done by assuming
            // that the current table exhibits the same ratios of the gas formation
            // volume factors and viscosities for identical pressure rations as in the
            // master table.
            const auto masterTable = pvtgTable.getInnerTable(masterTableIdx);
            const auto curTable = pvtgTable.getInnerTable(xIdx);
            for (int newRowIdx = 1;
                 newRowIdx < static_cast<int>(masterTable->numRows());
                 ++ newRowIdx)
            {
                Scalar alphaRv =
                    masterTable->getOilSolubilityColumn()[newRowIdx]
                    / masterTable->getOilSolubilityColumn()[0];

                Scalar alphaBg =
                    masterTable->getGasFormationFactorColumn()[newRowIdx]
                    / masterTable->getGasFormationFactorColumn()[0];

                Scalar alphaMug =
                    masterTable->getGasViscosityColumn()[newRowIdx]
                    / masterTable->getGasViscosityColumn()[0];

                Scalar newRv = curTable->getOilSolubilityColumn()[0]*alphaRv;
                Scalar newBg = curTable->getGasFormationFactorColumn()[0]*alphaBg;
                Scalar newMug = curTable->getGasViscosityColumn()[0]*alphaMug;

                invGasB.appendSamplePoint(xIdx, newRv, 1.0/newBg);
                gasMu.appendSamplePoint(xIdx, newRv, newMug);
            }
        }
    }
#endif // HAVE_OPM_PARSER

    /*!
     * \brief Finish initializing the gas phase PVT properties.
     */
    void initEnd()
    {
        // calculate the final 2D functions which are used for interpolation.
        int numRegions = gasMu_.size();
        for (int regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            // calculate the table which stores the inverse of the product of the gas
            // formation volume factor and the gas viscosity
            const auto& gasMu = gasMu_[regionIdx];
            const auto& invGasB = inverseGasB_[regionIdx];
            assert(gasMu.numX() == invGasB.numX());

            auto& invGasBMu = inverseGasBMu_[regionIdx];

            for (int pIdx = 0; pIdx < gasMu.numX(); ++pIdx) {
                invGasBMu.appendXPos(gasMu.xAt(pIdx));

                assert(gasMu.numY(pIdx) == invGasB.numY(pIdx));

                int numPressures = gasMu.numY(pIdx);
                for (int rvIdx = 0; rvIdx < numPressures; ++rvIdx)
                    invGasBMu.appendSamplePoint(pIdx,
                                                gasMu.yAt(pIdx, rvIdx),
                                                invGasB.valueAt(pIdx, rvIdx)*
                                                1/gasMu.valueAt(pIdx, rvIdx));
            }

            updateSaturationPressureSpline_(regionIdx);
        }
    }

private:
    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class LhsEval>
    LhsEval viscosity_(int regionIdx,
                       const LhsEval& temperature,
                       const LhsEval& pressure,
                       const LhsEval& XgO) const
    {
        const LhsEval& Rv =
            XgO/(1 - XgO)
            * (BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx)
               / BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx));

        const LhsEval& invBg = inverseGasB_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);
        const LhsEval& invMugBg = inverseGasBMu_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);

        return invBg/invMugBg;
    }

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    template <class LhsEval>
    LhsEval density_(int regionIdx,
                     const LhsEval& temperature,
                     const LhsEval& pressure,
                     const LhsEval& XgO) const
    {
        Scalar rhooRef = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx);
        Scalar rhogRef = BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx);

        const LhsEval& Bg = formationVolumeFactor_(regionIdx, temperature, pressure, XgO);
        LhsEval rhoo = rhooRef/Bg;

        // the oil formation volume factor just represents the partial density of the gas
        // component in the gas phase. to get the total density of the phase, we have to
        // add the partial density of the oil component.
        const LhsEval& Rv = XgO/(1 - XgO) * (rhooRef/rhogRef);
        rhoo += (rhogRef*Rv)/Bg;

        return rhoo;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class LhsEval>
    LhsEval formationVolumeFactor_(int regionIdx,
                                   const LhsEval& temperature,
                                   const LhsEval& pressure,
                                   const LhsEval& XgO) const
    {
        const LhsEval& Rv =
            XgO/(1-XgO)
            *BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx)
            /BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx);

        return 1.0 / inverseGasB_[regionIdx].eval(pressure, Rv, /*extrapolate=*/true);
    }

    /*!
     * \brief Returns the fugacity coefficient [Pa] of a component in the fluid phase given
     *        a set of parameters.
     */
    template <class LhsEval>
    LhsEval fugacityCoefficient_(int regionIdx,
                                 const LhsEval& temperature,
                                 const LhsEval& pressure,
                                 int compIdx) const
    {
        typedef Opm::MathToolbox<LhsEval> Toolbox;

        // set the oil component fugacity coefficient in oil phase
        // arbitrarily. we use some pseudo-realistic value for the vapor
        // pressure to ease physical interpretation of the results
        LhsEval phi_gG = Toolbox::createConstant(1.0);

        if (compIdx == BlackOilFluidSystem::gasCompIdx)
            return phi_gG;
        else if (compIdx == BlackOilFluidSystem::waterCompIdx)
            // assume that the affinity of the water component to the gas phase is much
            // smaller than that of the gas component
            return 1e8*phi_gG;

        /////////////
        // the rest of this method determines the fugacity coefficient
        // of the oil component:
        //
        // first, retrieve the mole fraction of gas a saturated oil
        // would exhibit at the given pressure
        const LhsEval& x_gOSat = saturatedGasOilMoleFraction_(regionIdx, temperature, pressure);

        // then, scale the oil component's gas phase fugacity
        // coefficient, so that the oil phase ends up at the right
        // composition if we were doing a flash experiment
        const LhsEval& phi_oO = BlackOilFluidSystem::fugCoefficientInOil(oilCompIdx,
                                                                         temperature,
                                                                         pressure,
                                                                         regionIdx);

        return phi_oO / x_gOSat;
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    template <class LhsEval>
    LhsEval oilVaporizationFactor_(int regionIdx,
                                   const LhsEval& temperature,
                                   const LhsEval& pressure) const
    { return oilVaporizationFactorTable_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the saturation pressure of the gas phase [Pa]
     *        depending on its mass fraction of the oil component
     *
     * \param XgO The mass fraction of the oil component in the gas phase [-]
     */
    template <class LhsEval>
    LhsEval gasSaturationPressure_(int regionIdx,
                                   const LhsEval& temperature,
                                   const LhsEval& XgO) const
    {
        typedef Opm::MathToolbox<LhsEval> Toolbox;

        // use the saturation pressure spline to get a pretty good initial value
        LhsEval pSat = saturationPressureSpline_[regionIdx].eval(XgO, /*extrapolate=*/true);
        const LhsEval& eps = pSat*1e-11;

        // Newton method to do the remaining work. If the initial
        // value is good, this should only take two to three
        // iterations...
        for (int i = 0; i < 20; ++i) {
            const LhsEval& f = saturatedGasOilMassFraction_(regionIdx, temperature, pSat) - XgO;
            const LhsEval& fPrime = ((saturatedGasOilMassFraction_(regionIdx, temperature, pSat + eps) - XgO) - f)/eps;

            const LhsEval& delta = f/fPrime;
            pSat -= delta;

            if (std::abs(Toolbox::value(delta)) < std::abs(Toolbox::value(pSat)) * 1e-10)
                return pSat;
        }

        OPM_THROW(NumericalIssue, "Could find the gas saturation pressure for X_g^O = " << XgO);
    }

    template <class LhsEval>
    LhsEval saturatedGasOilMassFraction_(int regionIdx,
                                         const LhsEval& temperature,
                                         const LhsEval& pressure) const
    {
        Scalar rho_gRef = BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx);
        Scalar rho_oRef = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx);

        // calculate the mass of the oil component [kg/m^3] in the oil phase. This is
        // equivalent to the gas dissolution factor [m^3/m^3] at current pressure times
        // the gas density [kg/m^3] at standard pressure
        const LhsEval& rho_oG = oilVaporizationFactor_(regionIdx, temperature, pressure) * rho_gRef;

        // we now have the total density of saturated oil and the partial density of the
        // oil component within it. The gas mass fraction is the ratio of these two.
        return rho_oG/(rho_oRef + rho_oG);
    }

    template <class LhsEval>
    LhsEval saturatedGasOilMoleFraction_(int regionIdx,
                                         const LhsEval& temperature,
                                         const LhsEval& pressure) const
    {
        // calculate the mass fractions of gas and oil
        const LhsEval& XgO = saturatedGasOilMassFraction_(regionIdx, temperature, pressure);

        // which can be converted to mole fractions, given the
        // components' molar masses
        Scalar MG = BlackOilFluidSystem::molarMass(gasCompIdx, regionIdx);
        Scalar MO = BlackOilFluidSystem::molarMass(oilCompIdx, regionIdx);

        const LhsEval& avgMolarMass = MO/(1 + (1 - XgO)*(MO/MG - 1));
        const LhsEval& xgO = XgO*avgMolarMass/MG;

        return xgO;
    }

private:
    void updateSaturationPressureSpline_(int regionIdx)
    {
        auto& oilVaporizationFactor = oilVaporizationFactorTable_[regionIdx];

        // create the spline representing saturation pressure
        // depending of the mass fraction in gas
        int n = oilVaporizationFactor.numSamples()*5;
        int delta = (oilVaporizationFactor.xMax() - oilVaporizationFactor.xMin())/(n + 1);

        SamplingPoints pSatSamplePoints;
        Scalar XgO = 0;
        for (int i=0; i <= n; ++ i) {
            Scalar pSat = oilVaporizationFactor.xMin() + i*delta;
            XgO = saturatedGasOilMassFraction_(regionIdx, /*temperature=*/Scalar(1e100), pSat);

            std::pair<Scalar, Scalar> val(XgO, pSat);
            pSatSamplePoints.push_back(val);
        }
        saturationPressureSpline_[regionIdx].setContainerOfTuples(pSatSamplePoints,
                                                                  /*type=*/Spline::Monotonic);
    }

    std::vector<TabulatedTwoDFunction> inverseGasB_;
    std::vector<TabulatedTwoDFunction> gasMu_;
    std::vector<TabulatedTwoDFunction> inverseGasBMu_;
    std::vector<TabulatedOneDFunction> oilVaporizationFactorTable_;
    std::vector<Spline> saturationPressureSpline_;
};

} // namespace Opm

#endif
