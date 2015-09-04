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

#include "OilPvtInterface.hpp"

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
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil phas
 *        with dissolved gas.
 */
template <class Scalar, class Evaluation = Scalar>
class LiveOilPvt : public OilPvtInterfaceTemplateWrapper<Scalar, Evaluation, LiveOilPvt<Scalar, Evaluation> >
{
    friend class OilPvtInterfaceTemplateWrapper<Scalar, Evaluation, LiveOilPvt<Scalar, Evaluation> >;

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
        if (static_cast<int>(inverseOilBTable_.size()) < numRegions) {
            inverseOilBTable_.resize(numRegions);
            inverseOilBMuTable_.resize(numRegions);
            oilMuTable_.resize(numRegions);
            gasDissolutionFactorTable_.resize(numRegions);
            saturationPressureSpline_.resize(numRegions);
        }
    }

    /*!
     * \brief Initialize the function for the gas dissolution factor \f$R_s\f$
     *
     * \param samplePoints A container of (x,y) values.
     */
    void setSaturatedOilGasDissolutionFactor(int regionIdx, const SamplingPoints &samplePoints)
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
    void setSaturatedOilFormationVolumeFactor(int regionIdx, const SamplingPoints &samplePoints)
    {
        auto& invOilB = inverseOilBTable_[regionIdx];

        auto &Rs = gasDissolutionFactorTable_[regionIdx];

        Scalar T = BlackOilFluidSystem::surfaceTemperature;

        Scalar RsMin = 0.0;
        Scalar RsMax = Rs.eval(gasDissolutionFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

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

                Scalar poSat = oilSaturationPressure_(regionIdx, T, XoG);
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
    { inverseOilBTable_[regionIdx] = invBo; }

    /*!
     * \brief Initialize the viscosity of the oil phase.
     *
     * This is a function of \f$(R_s, p_o)\f$...
     */
    void setOilViscosity(int regionIdx, const TabulatedTwoDFunction& muo)
    { oilMuTable_[regionIdx] = muo; }

    /*!
     * \brief Initialize the phase viscosity for gas saturated oil
     *
     * The oil viscosity is a function of \f$(p_o, X_o^G)\f$, but this method only
     * requires the viscosity of gas-saturated oil (which only depends on pressure) while
     * there is assumed to be no dependence on the gas mass fraction...
     */
    void setSaturatedOilViscosity(int regionIdx, const SamplingPoints &samplePoints  )
    {
        auto& gasDissolutionFactor = gasDissolutionFactorTable_[regionIdx];

        Scalar RsMin = 0.0;
        Scalar RsMax = gasDissolutionFactor.eval(gasDissolutionFactorTable_[regionIdx].xMax(), /*extrapolate=*/true);

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

#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize the oil parameters via the data specified by the PVTO ECL keyword.
     */
    void setPvtoTable(int regionIdx, const PvtoTable &pvtoTable)
    {
        const auto saturatedTable = pvtoTable.getOuterTable();
        assert(saturatedTable->numRows() > 1);

        auto& oilMu = oilMuTable_[regionIdx];
        auto& invOilB = inverseOilBTable_[regionIdx];
        auto& gasDissolutionFactor = gasDissolutionFactorTable_[regionIdx];

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
        int numRegions = oilMuTable_.size();
        for (int regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            // calculate the table which stores the inverse of the product of the oil
            // formation volume factor and the oil viscosity
            const auto& oilMu = oilMuTable_[regionIdx];
            const auto& invOilB = inverseOilBTable_[regionIdx];
            assert(oilMu.numX() == invOilB.numX());

            auto& invOilBMu = inverseOilBMuTable_[regionIdx];

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

private:
    /*!
     * \brief Returns the dynamic viscosity [Pa s] of the fluid phase given a set of parameters.
     */
    template <class LhsEval>
    LhsEval viscosity_(int regionIdx,
                       const LhsEval& temperature,
                       const LhsEval& pressure,
                       const LhsEval& XoG) const
    {
        const LhsEval& Rs =
            XoG/(1 - XoG)
            * (BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx)
               / BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx));

        // ATTENTION: Rs is the first axis!
        const LhsEval& invBo = inverseOilBTable_[regionIdx].eval(Rs, pressure, /*extrapolate=*/true);
        const LhsEval& invMuoBo = inverseOilBMuTable_[regionIdx].eval(Rs, pressure, /*extrapolate=*/true);

        return invBo/invMuoBo;
    }

    /*!
     * \brief Returns the density [kg/m^3] of the fluid phase given a set of parameters.
     */
    template <class LhsEval>
    LhsEval density_(int regionIdx,
                     const LhsEval& temperature,
                     const LhsEval& pressure,
                     const LhsEval& XoG) const
    {
        Scalar rhooRef = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx);
        Scalar rhogRef = BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx);
        Valgrind::CheckDefined(rhooRef);
        Valgrind::CheckDefined(rhogRef);

        const LhsEval& Bo = formationVolumeFactor_(regionIdx, temperature, pressure, XoG);
        Valgrind::CheckDefined(Bo);

        LhsEval rhoo = rhooRef/Bo;

        // the oil formation volume factor just represents the partial density of the oil
        // component in the oil phase. to get the total density of the phase, we have to
        // add the partial density of the gas component.
        const LhsEval Rs = XoG/(1 - XoG) * rhooRef/rhogRef;
        rhoo += rhogRef*Rs/Bo;

        return rhoo;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class LhsEval>
    LhsEval formationVolumeFactor_(int regionIdx,
                                   const LhsEval& temperature,
                                   const LhsEval& pressure,
                                   const LhsEval& XoG) const
    {
        const LhsEval& Rs =
            XoG/(1-XoG)
            *(BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx)
              /BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx));
        Valgrind::CheckDefined(Rs);
        Valgrind::CheckDefined(XoG);

        // ATTENTION: Rs is represented by the _first_ axis!
        return 1.0 / inverseOilBTable_[regionIdx].eval(Rs, pressure, /*extrapolate=*/true);
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
        // set the oil component fugacity coefficient in oil phase
        // arbitrarily. we use some pseudo-realistic value for the vapor
        // pressure to ease physical interpretation of the results
        const LhsEval& phi_oO = 20e3/pressure;

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
        const LhsEval& x_oGSat = saturatedOilGasMoleFraction_(regionIdx, temperature, pressure);

        // then, scale the gas component's gas phase fugacity
        // coefficient, so that the oil phase ends up at the right
        // composition if we were doing a flash experiment
        const LhsEval& phi_gG = BlackOilFluidSystem::fugCoefficientInGas(gasCompIdx,
                                                                         temperature,
                                                                         pressure,
                                                                         regionIdx);

        return phi_gG / x_oGSat;
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    template <class LhsEval>
    LhsEval gasDissolutionFactor_(int regionIdx,
                                  const LhsEval& temperature,
                                  const LhsEval& pressure) const
    { return gasDissolutionFactorTable_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param XoG The mass fraction of the gas component in the oil phase [-]
     */
    template <class LhsEval>
    LhsEval oilSaturationPressure_(int regionIdx,
                                   const LhsEval& temperature,
                                   const LhsEval& XoG) const
    {
        typedef Opm::MathToolbox<LhsEval> Toolbox;

        // use the saturation pressure spline to get a pretty good initial value
        LhsEval pSat = saturationPressureSpline_[regionIdx].eval(XoG, /*extrapolate=*/true);
        LhsEval eps = pSat*1e-11;

        // Newton method to do the remaining work. If the initial
        // value is good, this should only take two to three
        // iterations...
        for (int i = 0; i < 20; ++i) {
            const LhsEval& f = saturatedOilGasMassFraction_(regionIdx, temperature, pSat) - XoG;
            const LhsEval& fPrime = ((saturatedOilGasMassFraction_(regionIdx, temperature, pSat + eps) - XoG) - f)/eps;

            const LhsEval& delta = f/fPrime;
            pSat -= delta;

            Scalar absDelta = std::abs(Toolbox::value(delta));
            if (absDelta < Toolbox::value(pSat) * 1e-10 || absDelta < 1e-4)
                return pSat;
        }

        OPM_THROW(NumericalIssue, "Could find the oil saturation pressure for X_o^G = " << XoG);
    }

    template <class LhsEval>
    LhsEval saturatedOilGasMassFraction_(int regionIdx,
                                         const LhsEval& temperature,
                                         const LhsEval& pressure) const
    {
        Scalar rho_gRef = BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx);
        Scalar rho_oRef = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx);

        // calculate the mass of the gas component [kg/m^3] in the oil phase. This is
        // equivalent to the gas dissolution factor [m^3/m^3] at current pressure times
        // the gas density [kg/m^3] at standard pressure
        const LhsEval& rho_oG = gasDissolutionFactor_(regionIdx, temperature, pressure) * rho_gRef;

        // we now have the total density of saturated oil and the partial density of the
        // gas component within it. The gas mass fraction is the ratio of these two.
        return rho_oG/(rho_oRef + rho_oG);
    }

    template <class LhsEval>
    LhsEval saturatedOilGasMoleFraction_(int regionIdx,
                                         const LhsEval& temperature,
                                         const LhsEval& pressure) const
    {
        // calculate the mass fractions of gas and oil
        const LhsEval& XoG = saturatedOilGasMassFraction_(regionIdx, temperature, pressure);

        // which can be converted to mole fractions, given the
        // components' molar masses
        Scalar MG = BlackOilFluidSystem::molarMass(gasCompIdx, regionIdx);
        Scalar MO = BlackOilFluidSystem::molarMass(oilCompIdx, regionIdx);

        LhsEval avgMolarMass = MO/(1 + XoG*(MO/MG - 1));
        return XoG*avgMolarMass/MG;
    }

    void updateSaturationPressureSpline_(int regionIdx)
    {
        auto& gasDissolutionFactor = gasDissolutionFactorTable_[regionIdx];

        // create the spline representing saturation pressure
        // depending of the mass fraction in gas
        int n = gasDissolutionFactor.numSamples()*5;
        int delta = (gasDissolutionFactor.xMax() - gasDissolutionFactor.xMin())/(n + 1);

        SamplingPoints pSatSamplePoints;
        Scalar XoG = 0;
        for (int i=0; i <= n; ++ i) {
            Scalar pSat = gasDissolutionFactor.xMin() + i*delta;
            XoG = saturatedOilGasMassFraction_(regionIdx,
                                               /*temperature=*/Scalar(1e100),
                                               pSat);

            std::pair<Scalar, Scalar> val(XoG, pSat);
            pSatSamplePoints.push_back(val);
        }
        saturationPressureSpline_[regionIdx].setContainerOfTuples(pSatSamplePoints,
                                                                  /*type=*/Spline::Monotonic);
    }

    std::vector<TabulatedTwoDFunction> inverseOilBTable_;
    std::vector<TabulatedTwoDFunction> oilMuTable_;
    std::vector<TabulatedTwoDFunction> inverseOilBMuTable_;
    std::vector<TabulatedOneDFunction> gasDissolutionFactorTable_;
    std::vector<Spline> saturationPressureSpline_;
};

} // namespace Opm

#endif
