// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::LiveOilPvt
 */
#ifndef OPM_LIVE_OIL_PVT_HPP
#define OPM_LIVE_OIL_PVT_HPP

#include <opm/material/Constants.hpp>
#include <opm/material/common/OpmFinal.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>

#if HAVE_ECL_INPUT
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#endif

#if HAVE_OPM_COMMON
#include <opm/common/OpmLog/OpmLog.hpp>
#endif

namespace Opm {
/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil phas
 *        with dissolved gas.
 */
template <class Scalar>
class LiveOilPvt
{
    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

public:
    typedef Opm::UniformXTabulated2DFunction<Scalar> TabulatedTwoDFunction;
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;

    LiveOilPvt()
    {
        vapPar2_ = 0.0;
    }

    LiveOilPvt(const std::vector<Scalar>& gasReferenceDensity,
               const std::vector<Scalar>& oilReferenceDensity,
               const std::vector<TabulatedTwoDFunction>& inverseOilBTable,
               const std::vector<TabulatedTwoDFunction>& oilMuTable,
               const std::vector<TabulatedTwoDFunction>& inverseOilBMuTable,
               const std::vector<TabulatedOneDFunction>& saturatedOilMuTable,
               const std::vector<TabulatedOneDFunction>& inverseSaturatedOilBTable,
               const std::vector<TabulatedOneDFunction>& inverseSaturatedOilBMuTable,
               const std::vector<TabulatedOneDFunction>& saturatedGasDissolutionFactorTable,
               const std::vector<TabulatedOneDFunction>& saturationPressure,
               Scalar vapPar2)
        : gasReferenceDensity_(gasReferenceDensity)
        , oilReferenceDensity_(oilReferenceDensity)
        , inverseOilBTable_(inverseOilBTable)
        , oilMuTable_(oilMuTable)
        , inverseOilBMuTable_(inverseOilBMuTable)
        , saturatedOilMuTable_(saturatedOilMuTable)
        , inverseSaturatedOilBTable_(inverseSaturatedOilBTable)
        , inverseSaturatedOilBMuTable_(inverseSaturatedOilBMuTable)
        , saturatedGasDissolutionFactorTable_(saturatedGasDissolutionFactorTable)
        , saturationPressure_(saturationPressure)
        , vapPar2_(vapPar2)
    { }

#if HAVE_ECL_INPUT
    /*!
     * \brief Initialize the oil parameters via the data specified by the PVTO ECL keyword.
     */
    void initFromDeck(const Deck& deck, const EclipseState& eclState)
    {
        const auto& pvtoTables = eclState.getTableManager().getPvtoTables();
        const auto& densityKeyword = deck.getKeyword("DENSITY");

        assert(pvtoTables.size() == densityKeyword.size());

        size_t numRegions = pvtoTables.size();
        setNumRegions(numRegions);

        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            Scalar rhoRefO = densityKeyword.getRecord(regionIdx).getItem("OIL").getSIDouble(0);
            Scalar rhoRefG = densityKeyword.getRecord(regionIdx).getItem("GAS").getSIDouble(0);
            Scalar rhoRefW = densityKeyword.getRecord(regionIdx).getItem("WATER").getSIDouble(0);

            setReferenceDensities(regionIdx, rhoRefO, rhoRefG, rhoRefW);
        }

        // initialize the internal table objects
        for (unsigned regionIdx = 0; regionIdx < numRegions; ++ regionIdx) {
            const auto& pvtoTable = pvtoTables[regionIdx];

            const auto& saturatedTable = pvtoTable.getSaturatedTable();
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
                const auto& tmpPressureColumn = saturatedTable.getColumn("P");
                const auto& tmpGasSolubilityColumn = saturatedTable.getColumn("RS");

                invSatOilB.setXYContainers(tmpPressureColumn, invSatOilBArray);
                satOilMu.setXYContainers(tmpPressureColumn, satOilMuArray);
                gasDissolutionFac.setXYContainers(tmpPressureColumn, tmpGasSolubilityColumn);
            }

            updateSaturationPressure_(regionIdx);
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
                    throw std::runtime_error("PVTO tables are invalid: The last table must exhibit at least one "
                                             "entry for undersaturated oil!");

                // extend the current table using the master table.
                extendPvtoTable_(regionIdx,
                                 xIdx,
                                 pvtoTable.getUnderSaturatedTable(xIdx),
                                 pvtoTable.getUnderSaturatedTable(masterTableIdx));
            }
        }

        vapPar2_ = 0.0;
        if (deck.hasKeyword("VAPPARS")) {
            const auto& vapParsKeyword = deck.getKeyword("VAPPARS");
            vapPar2_ = vapParsKeyword.getRecord(0).getItem("OIL_DENSITY_PROPENSITY").template get<double>(0);
        }

        initEnd();
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
#endif // HAVE_ECL_INPUT

    void setNumRegions(size_t numRegions)
    {
        oilReferenceDensity_.resize(numRegions);
        gasReferenceDensity_.resize(numRegions);
        inverseOilBTable_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::LeftExtreme});
        inverseOilBMuTable_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::LeftExtreme});
        inverseSaturatedOilBTable_.resize(numRegions);
        inverseSaturatedOilBMuTable_.resize(numRegions);
        oilMuTable_.resize(numRegions, TabulatedTwoDFunction{TabulatedTwoDFunction::InterpolationPolicy::LeftExtreme});
        saturatedOilMuTable_.resize(numRegions);
        saturatedGasDissolutionFactorTable_.resize(numRegions);
        saturationPressure_.resize(numRegions);
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
    void setSaturatedOilGasDissolutionFactor(unsigned regionIdx, const SamplingPoints& samplePoints)
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
    void setSaturatedOilFormationVolumeFactor(unsigned regionIdx, const SamplingPoints& samplePoints)
    {
        Scalar T = 273.15 + 15.56; // [K]
        auto& invOilB = inverseOilBTable_[regionIdx];

        updateSaturationPressure_(regionIdx);

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
     * \brief Initialize the function for the oil formation volume factor
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
    void setSaturatedOilViscosity(unsigned regionIdx, const SamplingPoints& samplePoints)
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

            updateSaturationPressure_(regionIdx);
        }
    }

    /*!
     * \brief Return the number of PVT regions which are considered by this PVT-object.
     */
    unsigned numRegions() const
    { return inverseOilBMuTable_.size(); }

    /*!
     * \brief Returns the specific enthalpy [J/kg] of oil given a set of parameters.
     */
    template <class Evaluation>
    Evaluation internalEnergy(unsigned regionIdx OPM_UNUSED,
                        const Evaluation& temperature OPM_UNUSED,
                        const Evaluation& pressure OPM_UNUSED,
                        const Evaluation& Rs OPM_UNUSED) const
    {
        throw std::runtime_error("Requested the enthalpy of oil but the thermal option is not enabled");
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
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation inverseFormationVolumeFactor(unsigned regionIdx,
                                            const Evaluation& /*temperature*/,
                                            const Evaluation& pressure,
                                            const Evaluation& Rs) const
    {
        // ATTENTION: Rs is represented by the _first_ axis!
        return inverseOilBTable_[regionIdx].eval(Rs, pressure, /*extrapolate=*/true);
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class Evaluation>
    Evaluation saturatedInverseFormationVolumeFactor(unsigned regionIdx,
                                                     const Evaluation& /*temperature*/,
                                                     const Evaluation& pressure) const
    {
        // ATTENTION: Rs is represented by the _first_ axis!
        return inverseSaturatedOilBTable_[regionIdx].eval(pressure, /*extrapolate=*/true);
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
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     *
     * This variant of the method prevents all the oil to be vaporized even if the gas
     * phase is still not saturated. This is physically quite dubious but it corresponds
     * to how the Eclipse 100 simulator handles this. (cf the VAPPARS keyword.)
     */
    template <class Evaluation>
    Evaluation saturatedGasDissolutionFactor(unsigned regionIdx,
                                             const Evaluation& /*temperature*/,
                                             const Evaluation& pressure,
                                             const Evaluation& oilSaturation,
                                             Evaluation maxOilSaturation) const
    {
        Evaluation tmp =
            saturatedGasDissolutionFactorTable_[regionIdx].eval(pressure, /*extrapolate=*/true);

        // apply the vaporization parameters for the gas phase (cf. the Eclipse VAPPARS
        // keyword)
        maxOilSaturation = Opm::min(maxOilSaturation, Scalar(1.0));
        if (vapPar2_ > 0.0 && maxOilSaturation > 0.01 && oilSaturation < maxOilSaturation) {
            static const Scalar eps = 0.001;
            const Evaluation& So = Opm::max(oilSaturation, eps);
            tmp *= Opm::max(1e-3, Opm::pow(So/maxOilSaturation, vapPar2_));
        }

        return tmp;
    }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param Rs The surface volume of gas component dissolved in what will yield one cubic meter of oil at the surface [-]
     */
    template <class Evaluation>
    Evaluation saturationPressure(unsigned regionIdx,
                                  const Evaluation& temperature OPM_UNUSED,
                                  const Evaluation& Rs) const
    {
        typedef Opm::MathToolbox<Evaluation> Toolbox;

        const auto& RsTable = saturatedGasDissolutionFactorTable_[regionIdx];
        const Scalar eps = std::numeric_limits<typename Toolbox::Scalar>::epsilon()*1e6;

        // use the saturation pressure function to get a pretty good initial value
        Evaluation pSat = saturationPressure_[regionIdx].eval(Rs, /*extrapolate=*/true);

        // Newton method to do the remaining work. If the initial
        // value is good, this should only take two to three
        // iterations...
        bool onProbation = false;
        for (int i = 0; i < 20; ++i) {
            const Evaluation& f = RsTable.eval(pSat, /*extrapolate=*/true) - Rs;
            const Evaluation& fPrime = RsTable.evalDerivative(pSat, /*extrapolate=*/true);

            // If the derivative is "zero" Newton will not converge,
            // so simply return our initial guess.
            if (std::abs(Opm::scalarValue(fPrime)) < 1.0e-30) {
                return pSat;
            }

            const Evaluation& delta = f/fPrime;

            pSat -= delta;

            if (pSat < 0.0) {
                // if the pressure is lower than 0 Pascals, we set it back to 0. if this
                // happens twice, we give up and just return 0 Pa...
                if (onProbation)
                    return 0.0;

                onProbation = true;
                pSat = 0.0;
            }

            if (std::abs(Opm::scalarValue(delta)) < std::abs(Opm::scalarValue(pSat))*eps)
                return pSat;
        }

        std::stringstream errlog;
        errlog << "Finding saturation pressure did not converge:"
               << " pSat = " << pSat
               << ", Rs = " << Rs;
#if HAVE_OPM_COMMON
        OpmLog::debug("Live oil saturation pressure", errlog.str());
#endif
        throw NumericalIssue(errlog.str());
    }

    const std::vector<Scalar>& gasReferenceDensity() const
    { return gasReferenceDensity_; }

    const std::vector<Scalar>& oilReferenceDensity() const
    { return oilReferenceDensity_; }

    const std::vector<TabulatedTwoDFunction>& inverseOilBTable() const
    { return inverseOilBTable_; }

    const std::vector<TabulatedTwoDFunction>& oilMuTable() const
    { return oilMuTable_; }

    const std::vector<TabulatedTwoDFunction>& inverseOilBMuTable() const
    { return inverseOilBMuTable_; }

    const std::vector<TabulatedOneDFunction>& saturatedOilMuTable() const
    { return saturatedOilMuTable_; }

    const std::vector<TabulatedOneDFunction>& inverseSaturatedOilBTable() const
    { return inverseSaturatedOilBTable_; }

    const std::vector<TabulatedOneDFunction>& inverseSaturatedOilBMuTable() const
    { return inverseSaturatedOilBMuTable_; }

    const std::vector<TabulatedOneDFunction>& saturatedGasDissolutionFactorTable() const
    { return saturatedGasDissolutionFactorTable_; }

    const std::vector<TabulatedOneDFunction>& saturationPressure() const
    { return saturationPressure_; }

    Scalar vapPar2() const
    { return vapPar2_; }

    bool operator==(const LiveOilPvt<Scalar>& data) const
    {
        return this->gasReferenceDensity() == data.gasReferenceDensity() &&
               this->oilReferenceDensity() == data.oilReferenceDensity() &&
               this->inverseOilBTable() == data.inverseOilBTable() &&
               this->oilMuTable() == data.oilMuTable() &&
               this->inverseOilBMuTable() == data.inverseOilBMuTable() &&
               this->saturatedOilMuTable() == data.saturatedOilMuTable() &&
               this->inverseSaturatedOilBTable() == data.inverseSaturatedOilBTable() &&
               this->inverseSaturatedOilBMuTable() == data.inverseSaturatedOilBMuTable() &&
               this->saturatedGasDissolutionFactorTable() == data.saturatedGasDissolutionFactorTable() &&
               this->vapPar2() == data.vapPar2();
    }

private:
    void updateSaturationPressure_(unsigned regionIdx)
    {
        typedef std::pair<Scalar, Scalar> Pair;
        const auto& gasDissolutionFac = saturatedGasDissolutionFactorTable_[regionIdx];

        // create the function representing saturation pressure depending of the mass
        // fraction in gas
        size_t n = gasDissolutionFac.numSamples();
        Scalar delta = (gasDissolutionFac.xMax() - gasDissolutionFac.xMin())/Scalar(n + 1);

        SamplingPoints pSatSamplePoints;
        Scalar Rs = 0;
        for (size_t i=0; i <= n; ++ i) {
            Scalar pSat = gasDissolutionFac.xMin() + Scalar(i)*delta;
            Rs = saturatedGasDissolutionFactor(regionIdx,
                                               /*temperature=*/Scalar(1e30),
                                               pSat);

            Pair val(Rs, pSat);
            pSatSamplePoints.push_back(val);
        }

        //Prune duplicate Rs values (can occur, and will cause problems in further interpolation)
        auto x_coord_comparator = [](const Pair& a, const Pair& b) { return a.first == b.first; };
        auto last = std::unique(pSatSamplePoints.begin(), pSatSamplePoints.end(), x_coord_comparator);
        pSatSamplePoints.erase(last, pSatSamplePoints.end());

        saturationPressure_[regionIdx].setContainerOfTuples(pSatSamplePoints);
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
    std::vector<TabulatedOneDFunction> saturationPressure_;

    Scalar vapPar2_;
};

} // namespace Opm

#endif
