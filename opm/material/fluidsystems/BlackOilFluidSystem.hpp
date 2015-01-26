/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 * \copydoc Opm::FluidSystems::BlackOil
 */
#ifndef OPM_BLACK_OIL_FLUID_SYSTEM_HPP
#define OPM_BLACK_OIL_FLUID_SYSTEM_HPP

#include <opm/material/IdealGas.hpp>
#include <opm/material/Constants.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/fluidsystems/NullParameterCache.hpp>
#include <opm/material/UniformXTabulated2DFunction.hpp>
#include <opm/material/Tabulated1DFunction.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif // HAVE_OPM_PARSER

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Spline.hpp>

#include <array>
#include <vector>

namespace Opm {
namespace FluidSystems {

/*!
 * \brief A fluid system which uses the black-oil parameters
 *        to calculate termodynamically meaningful quantities.
 */
template <class Scalar>
class BlackOil
    : public BaseFluidSystem<Scalar, BlackOil<Scalar> >
{
    typedef Opm::UniformXTabulated2DFunction<Scalar> TabulatedTwoDFunction;
    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
    typedef Opm::Spline<Scalar> Spline;

    typedef std::vector<std::pair<Scalar, Scalar> > SamplingPoints;

public:
    //! \copydoc BaseFluidSystem::ParameterCache
    class ParameterCache : public Opm::NullParameterCache
    {
    public:
        ParameterCache(int regionIdx=0)
        { regionIdx_ = 0; }

        /*!
         * \brief Return the index of the region which should be used to determine the
         *        thermodynamic properties
         */
        int regionIndex() const
        { return regionIdx_; }

        /*!
         * \brief Set the index of the region which should be used to determine the
         *        thermodynamic properties
         */
        void setRegionIndex(int val)
        { regionIdx_ = val; }

    private:
        int regionIdx_;
    };

    /****************************************
     * Fluid phase parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numPhases
    static const int numPhases = 3;

    //! Index of the oil phase
    static const int oilPhaseIdx = 0;
    //! Index of the water phase
    static const int waterPhaseIdx = 1;
    //! Index of the gas phase
    static const int gasPhaseIdx = 2;

    //! The pressure at the surface
    static const Scalar surfacePressure;

    /*!
     * \copydoc BaseFluidSystem::init
     *
     * \attention For this fluid system, this method just throws a
     *            <tt>std::logic_error</tt> as there is no
     *            way to generically calculate the required black oil
     *            parameters. Instead of this method, use
     * \code
     * FluidSystem::initBegin();
     * // set the black oil parameters
     * FluidSystem::initEnd();
     * \endcode
     */
    static void init()
    {
        OPM_THROW(std::logic_error,
                  "There is no generic init() method for this fluid system. The "
                  << "black-oil fluid system must be initialized using:\n"
                  << "    FluidSystem::initBegin()\n"
                  << "    // set black oil parameters\n"
                  << "    FluidSystem::initEnd()\n");
    }

    /*!
     * \brief Begin the initialization of the black oil fluid system.
     *
     * After calling this method the reference densities, all dissolution and formation
     * volume factors, the oil bubble pressure, all viscosities and the water
     * compressibility must be set. Before the fluid system can be used, initEnd() must
     * be called to finalize the initialization.
     */
    static void initBegin()
    {
        setWaterReferenceFormationVolumeFactor(1.0);
        setWaterViscosibility(0.0);
        setWaterReferencePressure(surfacePressure);
    }

#if HAVE_OPM_PARSER
    /*!
     * \brief Sets the pressure-dependent oil viscosity, density and
     *        gas content using a table stemming from the Eclipse PVTO
     *        keyword.
     *
     * When calling this method, the reference densities of the fluids _must_ have already
     * been set!
     */
    static void setPvtoTable(const PvtoTable &pvtoTable, int regionIdx=0)
    {
        const auto saturatedTable = pvtoTable.getOuterTable();
        assert(saturatedTable->numRows() > 1);

        resizeArrays_(regionIdx);

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

    /*!
     * \brief Sets the pressure-dependent water viscosity and density
     *        using a table stemming from the Eclipse PVTW keyword.
     *
     * This function also sets the reference viscosity and the reference
     * density of water, but these can be overwritten using setReference*().
     */
    static void setPvtw(DeckKeywordConstPtr pvtwKeyword, int regionIdx=0)
    {
        assert(static_cast<int>(pvtwKeyword->size()) >= regionIdx);

        resizeArrays_(regionIdx);

        auto pvtwRecord = pvtwKeyword->getRecord(regionIdx);
        waterReferencePressure_[regionIdx] =
            pvtwRecord->getItem("P_REF")->getSIDouble(0);
        waterReferenceFormationVolumeFactor_[regionIdx] =
            pvtwRecord->getItem("WATER_VOL_FACTOR")->getSIDouble(0);
        waterCompressibility_[regionIdx] =
            pvtwRecord->getItem("WATER_COMPRESSIBILITY")->getSIDouble(0);
        waterViscosity__[regionIdx] =
            pvtwRecord->getItem("WATER_VISCOSITY")->getSIDouble(0);
        waterViscosibility_[regionIdx] =
            pvtwRecord->getItem("WATER_VISCOSIBILITY")->getSIDouble(0);
    }

    /*!
     * \brief Sets the pressure-dependent viscosity and density of dry
     *        gas using a table stemming from the Eclipse PVDG
     *        keyword.
     *
     * This function also sets the reference viscosity and the reference
     * density of gas, but these can be overwritten using
     * setReference*().
     */
    static void setPvdgTable(const PvdgTable &pvdgTable, int regionIdx=0)
    {
        int numSamples = pvdgTable.numRows();
        assert(numSamples > 1);

        resizeArrays_(regionIdx);

        // say 99.97% of all time: "premature optimization is the root of all
        // evil". Eclipse does it this way for no good reason!
        std::vector<Scalar> invB(pvdgTable.numRows());
        const auto& Bg = pvdgTable.getFormationFactorColumn();
        for (unsigned i = 0; i < Bg.size(); ++ i) {
            invB[i] = 1.0/Bg[i];
        }

        inverseGasB_[regionIdx].setXYArrays(numSamples, pvdgTable.getPressureColumn(), invB);
        gasMu_[regionIdx].setXYArrays(numSamples, pvdgTable.getPressureColumn(), pvdgTable.getViscosityColumn());
    }
#endif

    /*!
     * \brief Initialize the values of the reference densities
     *
     * \param rhoOil The reference density of (gas saturated) oil phase.
     * \param rhoWater The reference density of the water phase.
     * \param rhoGas The reference density of the gas phase.
     */
    static void setReferenceDensities(Scalar rhoOil,
                                      Scalar rhoWater,
                                      Scalar rhoGas,
                                      int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        referenceDensity_[regionIdx][oilPhaseIdx] = rhoOil;
        referenceDensity_[regionIdx][waterPhaseIdx] = rhoWater;
        referenceDensity_[regionIdx][gasPhaseIdx] = rhoGas;
    }

    /*!
     * \brief Initialize the function for the gas dissolution factor \f$R_s\f$
     *
     * \param samplePoints A container of (x,y) values.
     */
    static void setSaturatedOilGasDissolutionFactor(const SamplingPoints &samplePoints,
                                                    int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        gasDissolutionFactor_[regionIdx].setContainerOfTuples(samplePoints);
    }

    /*!
     * \brief Initialize the function for the oil formation volume factor
     *
     * The oil formation volume factor \f$B_o\f$ is a function of \f$(p_o, X_o^G)\f$ and
     * represents the partial density of the oil component in the oil phase at a given
     * pressure. This method only requires the volume factor of gas-saturated oil (which
     * only depends on pressure) while the dependence on the gas mass fraction is
     * guesstimated...
     */
    static void setSaturatedOilFormationVolumeFactor(const SamplingPoints &samplePoints,
                                                     int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        auto& invOilB = inverseOilB_[regionIdx];

        auto &Rs = gasDissolutionFactor_[regionIdx];

        Scalar RsMin = 0.0;
        Scalar RsMax = Rs.eval(gasDissolutionFactor_[regionIdx].xMax(), /*extrapolate=*/true);

        Scalar poMin = samplePoints.front().first;
        Scalar poMax = samplePoints.back().first;

        size_t nRs = 20;
        size_t nP = samplePoints.size()*2;

        Scalar rhogRef = referenceDensity(gasPhaseIdx, regionIdx);
        Scalar rhooRef = referenceDensity(oilPhaseIdx, regionIdx);

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

                Scalar poSat = oilSaturationPressure(XoG, regionIdx);
                Scalar BoSat = oilFormationVolumeFactorSpline.eval(poSat, /*extrapolate=*/true);
                Scalar drhoo_dp = (1.1200 - 1.1189)/((5000 - 4000)*6894.76);
                Scalar rhoo = referenceDensity(oilPhaseIdx, regionIdx)/BoSat*(1 + drhoo_dp*(po - poSat));

                Scalar Bo = referenceDensity(oilPhaseIdx, regionIdx)/rhoo;

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
    static void setInverseOilFormationVolumeFactor(const TabulatedTwoDFunction &invBo, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        inverseOilB_[regionIdx] = invBo;
    }

    /*!
     * \brief Initialize the spline for the viscosity of the oil phase.
     *
     * This is a function of \f$(R_s, p_o)\f$...
     */
    static void setOilViscosity(const TabulatedTwoDFunction &muo, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        oilMu_[regionIdx] = muo;
    }

    /*!
     * \brief Initialize the phase viscosity for gas saturated oil
     *
     * The oil viscosity is a function of \f$(p_o, X_o^G)\f$, but this method only
     * requires the viscosity of gas-saturated oil (which only depends on pressure) while
     * there is assumed to be no dependence on the gas mass fraction...
     */
    static void setSaturatedOilViscosity(const SamplingPoints &samplePoints, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

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

    /*!
     * \brief Initialize the function for the formation volume factor of dry gas
     *
     * \param samplePoints A container of \f$(p_g, B_g)\f$ values
     */
    static void setGasFormationVolumeFactor(const SamplingPoints &samplePoints, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        SamplingPoints tmp(samplePoints);
        auto it = tmp.begin();
        const auto& endIt = tmp.end();
        for (; it != endIt; ++ it)
            std::get<1>(*it) = 1.0/std::get<1>(*it);

        inverseGasB_[regionIdx].setContainerOfTuples(tmp);
        assert(inverseGasB_[regionIdx].monotonic());
    }

    /*!
     * \brief Initialize the function for the viscosity of dry gas
     *
     * \param samplePoints A container of \f$(p_g, \mu_g)\f$ values
     */
    static void setGasViscosity(const SamplingPoints &samplePoints, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        gasMu_[regionIdx].setContainerOfTuples(samplePoints);
    }

    /*!
     * \brief Set the water reference pressure [Pa]
     */
    static void setWaterReferencePressure(Scalar p, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        waterReferencePressure_[regionIdx] = p;
    }

    /*!
     * \brief Set the water reference formation volume factor [-]
     */
    static void setWaterReferenceFormationVolumeFactor(Scalar BwRef, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        waterReferenceFormationVolumeFactor_[regionIdx] = BwRef;
    }

    /*!
     * \brief Set the water reference viscosity [Pa s]
     */
    static void setWaterReferenceViscosity(Scalar muWater, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        waterViscosity__[regionIdx] = muWater;
    }

    /*!
     * \brief Set the water "viscosibility" [1/ (Pa s)]
     */
    static void setWaterViscosibility(Scalar muCompWater, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        waterViscosibility_[regionIdx] = muCompWater;
    }

    /*!
     * \brief Set the water compressibility [1 / Pa]
     */
    static void setWaterCompressibility(Scalar cWater, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        waterCompressibility_[regionIdx] = cWater;
    }

    /*!
     * \brief Finish initializing the black oil fluid system.
     */
    static void initEnd()
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

            // calculate the table which stores the inverse of the product of the gas
            // formation volume factor and the gas viscosity
            const auto& gasMu = gasMu_[regionIdx];
            const auto& invGasB = inverseGasB_[regionIdx];
            assert(gasMu.numSamples() == invGasB.numSamples());

            std::vector<Scalar> pressureValues(gasMu.numSamples());
            std::vector<Scalar> invGasBMuValues(gasMu.numSamples());
            for (int pIdx = 0; pIdx < gasMu.numSamples(); ++pIdx) {
                pressureValues[pIdx] = invGasB.xAt(pIdx);
                invGasBMuValues[pIdx] = invGasB.valueAt(pIdx) * (1.0/gasMu.valueAt(pIdx));
            }

            inverseGasBMu_[regionIdx].setXYContainers(pressureValues, invGasBMuValues);

            updateSaturationPressureSpline_(regionIdx);
        }

        // calculate molar masses

        // water is simple: 18 g/mol
        molarMass_[waterCompIdx] = 18e-3;

        // for gas, we take the density at standard conditions and assume it to be ideal
        Scalar p = surfacePressure;
        Scalar rho_g = referenceDensity_[/*regionIdx=*/0][gasPhaseIdx];
        Scalar T = 297.15;
        molarMass_[gasCompIdx] = Opm::Constants<Scalar>::R*T*rho_g / p;

        // finally, for oil phase, we take the molar mass from the
        // spe9 paper
        molarMass_[oilCompIdx] = 175e-3; // kg/mol
    }

    //! \copydoc BaseFluidSystem::phaseName
    static const char *phaseName(const int phaseIdx)
    {
        static const char *name[] = { "oil", "water", "gas" };

        assert(0 <= phaseIdx && phaseIdx < numPhases + 1);
        return name[phaseIdx];
    }

    //! \copydoc BaseFluidSystem::isLiquid
    static bool isLiquid(const int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return phaseIdx != gasPhaseIdx;
    }

    /****************************************
     * Component related parameters
     ****************************************/

    //! \copydoc BaseFluidSystem::numComponents
    static const int numComponents = 3;

    //! Index of the oil component
    static const int oilCompIdx = 0;
    //! Index of the water component
    static const int waterCompIdx = 1;
    //! Index of the gas component
    static const int gasCompIdx = 2;

    //! \copydoc BaseFluidSystem::componentName
    static const char *componentName(const int compIdx)
    {
        static const char *name[] = { "O", "W", "G" };

        assert(0 <= compIdx && compIdx < numComponents);
        return name[compIdx];
    }

    //! \copydoc BaseFluidSystem::molarMass
    static Scalar molarMass(const int compIdx)
    { return molarMass_[compIdx]; }

    //! \copydoc BaseFluidSystem::isIdealMixture
    static bool isIdealMixture(int phaseIdx)
    {
        // fugacity coefficients are only pressure dependent -> we
        // have an ideal mixture
        return true;
    }

    //! \copydoc BaseFluidSystem::isCompressible
    static bool isCompressible(int phaseIdx)
    { return true; /* all phases are compressible */ }

    //! \copydoc BaseFluidSystem::isIdealGas
    static bool isIdealGas(int phaseIdx)
    { return false; }

    /****************************************
     * thermodynamic relations
     ****************************************/
    //! \copydoc BaseFluidSystem::density
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          ParameterCache &paramCache,
                          const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        Scalar p = fluidState.pressure(phaseIdx);
        int regionIdx = paramCache.regionIndex();

        switch (phaseIdx) {
        case waterPhaseIdx: return waterDensity(p, regionIdx);
        case gasPhaseIdx: return gasDensity(p, regionIdx);
        case oilPhaseIdx: return oilDensity(p,
                                            fluidState.massFraction(oilPhaseIdx, gasCompIdx),
                                            regionIdx);
        }

        OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
    }

    //! \copydoc BaseFluidSystem::fugacityCoefficient
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      const ParameterCache &paramCache,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);
        assert(0 <= compIdx  && compIdx <= numComponents);

        Scalar p = fluidState.pressure(phaseIdx);
        int regionIdx = paramCache.regionIndex();

        switch (phaseIdx) {
        case waterPhaseIdx: return fugCoefficientInWater(compIdx, p, regionIdx);
        case gasPhaseIdx: return fugCoefficientInGas(compIdx, p, regionIdx);
        case oilPhaseIdx: return fugCoefficientInOil(compIdx, p, regionIdx);
        }

        OPM_THROW(std::logic_error, "Unhandled phase or component index");
    }

    //! \copydoc BaseFluidSystem::viscosity
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            const ParameterCache &paramCache,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx <= numPhases);

        Scalar p = fluidState.pressure(phaseIdx);
        int regionIdx = paramCache.regionIndex();

        switch (phaseIdx) {
        case oilPhaseIdx: {
            Scalar XoG = fluidState.massFraction(oilPhaseIdx, gasCompIdx);
            Scalar Rs =
                XoG/(1 - XoG)
                * referenceDensity(oilPhaseIdx)
                / referenceDensity(gasPhaseIdx);

            return oilViscosityRs_(p, Rs, regionIdx);
        }
        case waterPhaseIdx: return waterViscosity_(p, regionIdx);
        case gasPhaseIdx: return gasViscosity_(p, regionIdx);
        }

        OPM_THROW(std::logic_error, "Unhandled phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the density of a fluid phase at surface pressure [kg/m^3]
     *
     * \copydoc Doxygen::phaseIdxParam
     */
    static Scalar referenceDensity(int phaseIdx, int regionIdx=0)
    { return referenceDensity_[regionIdx][phaseIdx]; }

    /*!
     * \brief Returns the oil formation volume factor \f$B_o\f$ of saturated oil for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar saturatedOilFormationVolumeFactor(Scalar pressure, int regionIdx=0)
    {
        Valgrind::CheckDefined(pressure);

        // calculate the mass fractions of gas and oil
        Scalar XoG = saturatedOilGasMassFraction(pressure, regionIdx);

        // ATTENTION: XoG is represented by the _first_ axis!
        return inverseOilB_[regionIdx].eval(XoG, pressure, /*extrapolate=*/true);
    }

    /*!
     * \brief Return the formation volume factor of water.
     */
    static Scalar waterFormationVolumeFactor(Scalar pressure, int regionIdx=0)
    {
        // cf. ECLiPSE 2011 technical description, p. 116
        Scalar pRef = waterReferencePressure_[regionIdx];
        Scalar X = waterCompressibility_[regionIdx]*(pressure - pRef);

        Scalar BwRef = waterReferenceFormationVolumeFactor_[regionIdx];

        // TODO (?): consider the salt concentration of the brine
        return BwRef/(1 + X*(1 + X/2));
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar gasDissolutionFactor(Scalar pressure, int regionIdx=0)
    { return gasDissolutionFactor_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the water phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar fugCoefficientInWater(int compIdx, Scalar pressure, int regionIdx=0)
    {
        // set the affinity of the gas and oil components to the water
        // phase to be 6 orders of magnitute smaller than that of the
        // water component. for this we use a pseudo-realistic vapor
        // pressure of water as a starting point. (we just set it to
        // 30 kPa to ease interpreting the results.)
        const Scalar pvWater = 30e3;
        if (compIdx == oilCompIdx)
            return 1e6*pvWater / pressure;
        else if (compIdx == gasCompIdx)
            return 1.01e6*pvWater / pressure;

        return pvWater / pressure;
    }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the gas phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar fugCoefficientInGas(int compIdx, Scalar pressure, int regionIdx=0)
    {
        // make the gas component more affine to the gas phase than the other components
        if (compIdx == gasCompIdx)
            return 1e-3;
        return 1.0;
    }

    /*!
     * \brief Returns the fugacity coefficient of a given component in the oil phase
     *
     * \param compIdx The index of the component of interest
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar fugCoefficientInOil(int compIdx, Scalar pressure, int regionIdx=0)
    {
        // set the oil component fugacity coefficient in oil phase
        // arbitrarily. we use some pseudo-realistic value for the vapor
        // pressure to ease physical interpretation of the results
        Scalar phi_oO = 20e3/pressure;

        if (compIdx == oilCompIdx)
            return phi_oO;
        else if (compIdx == waterCompIdx)
            // assume that the affinity of the water component to the
            // oil phase is one million times smaller than that of the
            // oil component
            return 1e6*phi_oO;

        /////////////
        // the rest of this method determines the fugacity coefficient
        // of the gas component:
        //
        // first, retrieve the mole fraction of gas a saturated oil
        // would exhibit at the given pressure
        Scalar x_oGf = saturatedOilGasMoleFraction(pressure, regionIdx);

        // then, scale the gas component's gas phase fugacity
        // coefficient, so that the oil phase ends up at the right
        // composition if we were doing a flash experiment
        Scalar phi_gG = fugCoefficientInGas(gasCompIdx, pressure, regionIdx);
        return phi_gG / x_oGf;
    }

    /*!
     * \brief Returns the saturation pressure of the oil phase [Pa]
     *        depending on its mass fraction of the gas component
     *
     * \param X_oG The mass fraction of the gas component in the oil phase [-]
     */
    static Scalar oilSaturationPressure(Scalar X_oG, int regionIdx=0)
    {
        // use the saturation pressure spline to get a pretty good initial value
        Scalar pSat = saturationPressureSpline_[regionIdx].eval(X_oG, /*extrapolate=*/true);

        // Newton method to do the remaining work. If the initial
        // value is good, this should only take two to three
        // iterations...
        for (int i = 0; i < 20; ++i) {
            Scalar f = saturatedOilGasMassFraction(pSat, regionIdx) - X_oG;
            Scalar eps = pSat*1e-11;
            Scalar fPrime = ((saturatedOilGasMassFraction(pSat + eps, regionIdx) - X_oG) - f)/eps;

            Scalar delta = f/fPrime;
            pSat -= delta;

            if (std::abs(delta) < pSat * 1e-10)
                return pSat;
        }

        OPM_THROW(NumericalProblem, "Could find the oil saturation pressure for X_o^g = " << X_oG);
    }

    // the mass fraction of the gas component in the oil phase in a
    // flash experiment
    static Scalar saturatedOilGasMassFraction(Scalar pressure, int regionIdx=0)
    {
        Scalar rho_gRef = referenceDensity(gasPhaseIdx, regionIdx);

        // calculate the mass of the gas component [kg/m^3] in the oil phase. This is
        // equivalent to the gas dissolution factor [m^3/m^3] at current pressure times
        // the gas density [kg/m^3] at standard pressure
        Scalar rho_oG = gasDissolutionFactor(pressure, regionIdx) * rho_gRef;

        // we now have the total density of saturated oil and the partial density of the
        // gas component within it. The gas mass fraction is the ratio of these two.
        return rho_oG/(referenceDensity(oilPhaseIdx, regionIdx) + rho_oG);
    }

    // the mole fraction of the gas component of a gas-saturated oil phase
    static Scalar saturatedOilGasMoleFraction(Scalar pressure, int regionIdx=0)
    {
        // calculate the mass fractions of gas and oil
        Scalar XoG = saturatedOilGasMassFraction(pressure, regionIdx);

        // which can be converted to mole fractions, given the
        // components' molar masses
        Scalar MG = molarMass(gasCompIdx);
        Scalar MO = molarMass(oilCompIdx);

        Scalar avgMolarMass = MO/(1 + XoG*(MO/MG - 1));
        Scalar xoG = XoG*avgMolarMass/MG;

        return xoG;
    }

    /*!
     * \brief Return the normalized formation volume factor of (potentially)
     *        under-saturated oil.
     */
    static Scalar oilFormationVolumeFactor(Scalar oilPressure, Scalar XoG, int regionIdx=0)
    {
        Scalar Rs = XoG/(1-XoG)*referenceDensity(oilPhaseIdx)/referenceDensity(gasPhaseIdx);

        // ATTENTION: Rs is represented by the _first_ axis!
        return 1.0 / inverseOilB_[regionIdx].eval(Rs, oilPressure, /*extrapolate=*/true);
    }

    /*!
     * \brief Return the normalized formation volume factor of (potentially)
     *        under-saturated oil.
     *
     * This is the version of the method which takes the gas dissolution factor instead
     * of the gas mass fraction in oil as an argument.
     */
    static Scalar oilFormationVolumeFactorRs(Scalar oilPressure, Scalar Rs, int regionIdx=0)
    {
        // ATTENTION: Rs is represented by the _first_ axis!
        return 1.0 / inverseOilB_[regionIdx].eval(Rs, oilPressure, /*extrapolate=*/true);
    }

    /*!
     * \brief Return the density of (potentially) under-saturated oil.
     */
    static Scalar oilDensity(Scalar oilPressure, Scalar XoG, int regionIdx=0)
    {
        Scalar rhooRef = referenceDensity_[regionIdx][oilPhaseIdx];
        Scalar rhogRef = referenceDensity_[regionIdx][gasPhaseIdx];

        Scalar Bo = oilFormationVolumeFactor(oilPressure, XoG, regionIdx);
        Scalar rhoo = rhooRef/Bo;

        // the oil formation volume factor just represents the partial density of the oil
        // component in the oil phase. to get the total density of the phase, we have to
        // add the partial density of the gas component.
        Scalar Rs = XoG/(1 - XoG) * rhooRef/rhogRef;
        rhoo += rhogRef*Rs/Bo;

        return rhoo;
    }

    /*!
     * \brief Return the density of gas-saturated oil.
     */
    static Scalar saturatedOilDensity(Scalar pressure, int regionIdx=0)
    {
        // mass fraction of gas-saturated oil
        Scalar XoG = saturatedOilGasMassFraction(pressure, regionIdx);
        return oilDensity(pressure, XoG, regionIdx);
    }

    /*!
     * \brief Return the formation volume factor of gas.
     */
    static Scalar gasFormationVolumeFactor(Scalar pressure, int regionIdx=0)
    { return 1.0/inverseGasB_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    /*!
     * \brief Return the density of dry gas.
     */
    static Scalar gasDensity(Scalar gasPressure, int regionIdx)
    {
        // gas formation volume factor at reservoir pressure
        Scalar Bg = gasFormationVolumeFactor(gasPressure, regionIdx);
        return referenceDensity_[regionIdx][gasPhaseIdx]/Bg;
    }

    /*!
     * \brief Return the density of water.
     */
    static Scalar waterDensity(Scalar pressure, int regionIdx)
    {
        Scalar Bw = waterFormationVolumeFactor(pressure, regionIdx);
        Scalar rhowRef = referenceDensity(waterPhaseIdx, regionIdx);
        return rhowRef/Bw;
    }

private:
    static void resizeArrays_(int regionIdx)
    {
        if (static_cast<int>(inverseOilB_.size()) <= regionIdx) {
            int numRegions = regionIdx + 1;
            inverseOilB_.resize(numRegions);
            inverseGasB_.resize(numRegions);
            inverseGasBMu_.resize(numRegions);
            inverseOilBMu_.resize(numRegions);
            oilMu_.resize(numRegions);
            gasMu_.resize(numRegions);
            gasDissolutionFactor_.resize(numRegions);
            saturationPressureSpline_.resize(numRegions);
            waterReferencePressure_.resize(numRegions);
            waterReferenceFormationVolumeFactor_.resize(numRegions);
            waterCompressibility_.resize(numRegions);
            waterViscosity__.resize(numRegions);
            waterViscosibility_.resize(numRegions);

            referenceDensity_.resize(numRegions);
        }
    }

    static void updateSaturationPressureSpline_(int regionIdx)
    {
        resizeArrays_(regionIdx);

        auto& gasDissolutionFactor = gasDissolutionFactor_[regionIdx];

        // create the spline representing saturation pressure
        // depending of the mass fraction in gas
        int n = gasDissolutionFactor.numSamples()*5;
        int delta = (gasDissolutionFactor.xMax() - gasDissolutionFactor.xMin())/(n + 1);

        SamplingPoints pSatSamplePoints;
        Scalar X_oG = 0;
        for (int i=0; i <= n; ++ i) {
            Scalar pSat = gasDissolutionFactor.xMin() + i*delta;
            X_oG = saturatedOilGasMassFraction(pSat, regionIdx);

            std::pair<Scalar, Scalar> val(X_oG, pSat);
            pSatSamplePoints.push_back(val);
        }
        saturationPressureSpline_[regionIdx].setContainerOfTuples(pSatSamplePoints,
                                                                  /*type=*/Spline::Monotonic);
    }

    static Scalar oilViscosityRs_(Scalar oilPressure, Scalar Rs, int regionIdx)
    {
        // ATTENTION: Rs is the first axis!
        Scalar invBo = inverseOilB_[regionIdx].eval(Rs, oilPressure, /*extrapolate=*/true);
        Scalar invMuoBo = inverseOilBMu_[regionIdx].eval(Rs, oilPressure, /*extrapolate=*/true);

        return invBo/invMuoBo;
    }

    static Scalar gasViscosity_(Scalar gasPressure, int regionIdx)
    {
        Scalar invBg = inverseGasB_[regionIdx].eval(gasPressure, /*extrapolate=*/true);
        Scalar invMugBg = inverseGasBMu_[regionIdx].eval(gasPressure, /*extrapolate=*/true);

        return invBg/invMugBg;
    }

    static Scalar waterViscosity_(Scalar pressure, int regionIdx)
    {
        // Eclipse calculates the viscosity in a weird way: it
        // calcultes the product of B_w and mu_w and then divides the
        // result by B_w...
        Scalar BwMuwRef = waterViscosity__[regionIdx]*waterReferenceFormationVolumeFactor_[regionIdx];
        Scalar Bw = waterFormationVolumeFactor(pressure, regionIdx);

        Scalar pRef = waterReferencePressure_[regionIdx];
        Scalar Y =
            (waterCompressibility_[regionIdx] - waterViscosibility_[regionIdx])
            * (pressure - pRef);
        return BwMuwRef/((1 + Y*(1 + Y/2))*Bw);
    }

    static std::vector<TabulatedTwoDFunction> inverseOilB_;
    static std::vector<TabulatedTwoDFunction> oilMu_;
    static std::vector<TabulatedTwoDFunction> inverseOilBMu_;
    static std::vector<TabulatedOneDFunction> gasDissolutionFactor_;
    static std::vector<Spline> saturationPressureSpline_;

    static std::vector<TabulatedOneDFunction> inverseGasB_;
    static std::vector<TabulatedOneDFunction> gasMu_;
    static std::vector<TabulatedOneDFunction> inverseGasBMu_;

    static std::vector<Scalar> waterReferencePressure_;
    static std::vector<Scalar> waterReferenceFormationVolumeFactor_;
    static std::vector<Scalar> waterCompressibility_;
    static std::vector<Scalar> waterViscosity__;
    static std::vector<Scalar> waterViscosibility_;

    // HACK for GCC 4.4: the array size has to be specified using the literal value '3'
    // here, because GCC 4.4 seems to be unable to determine the number of phases from
    // the BlackOil fluid system in the attribute declaration below...
    static std::vector<std::array<Scalar, /*numPhases=*/3> > referenceDensity_;

    static Scalar molarMass_[numComponents];
};

template <class Scalar>
const Scalar
BlackOil<Scalar>::surfacePressure = 101325.0; // [Pa]

template <class Scalar>
std::vector<typename BlackOil<Scalar>::TabulatedTwoDFunction>
BlackOil<Scalar>::inverseOilB_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::TabulatedTwoDFunction>
BlackOil<Scalar>::inverseOilBMu_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::TabulatedOneDFunction>
BlackOil<Scalar>::gasDissolutionFactor_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::Spline>
BlackOil<Scalar>::saturationPressureSpline_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::TabulatedOneDFunction>
BlackOil<Scalar>::inverseGasB_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::TabulatedOneDFunction>
BlackOil<Scalar>::inverseGasBMu_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::TabulatedOneDFunction>
BlackOil<Scalar>::gasMu_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::TabulatedTwoDFunction>
BlackOil<Scalar>::oilMu_;

template <class Scalar>
std::vector<Scalar>
BlackOil<Scalar>::waterReferencePressure_;

template <class Scalar>
std::vector<Scalar>
BlackOil<Scalar>::waterReferenceFormationVolumeFactor_;

template <class Scalar>
std::vector<Scalar>
BlackOil<Scalar>::waterCompressibility_;

template <class Scalar>
std::vector<Scalar>
BlackOil<Scalar>::waterViscosity__;

template <class Scalar>
std::vector<Scalar>
BlackOil<Scalar>::waterViscosibility_;

template <class Scalar>
std::vector<std::array<Scalar, 3> >
BlackOil<Scalar>::referenceDensity_;

template <class Scalar>
Scalar
BlackOil<Scalar>::molarMass_[BlackOil<Scalar>::numComponents];
}} // namespace Opm, FluidSystems

#endif
