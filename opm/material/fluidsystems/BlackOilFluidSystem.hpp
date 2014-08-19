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

#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Spline.hpp>

#include <opm/material/IdealGas.hpp>
#include <opm/material/Constants.hpp>
#include <opm/material/components/H2O.hpp>
#include <opm/material/fluidsystems/BaseFluidSystem.hpp>
#include <opm/material/fluidsystems/NullParameterCache.hpp>
#include <opm/material/UniformXTabulated2DFunction.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/Utility/PvtoTable.hpp>
#include <opm/parser/eclipse/Utility/PvtwTable.hpp>
#include <opm/parser/eclipse/Utility/PvdgTable.hpp>
#endif // HAVE_OPM_PARSER

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
    typedef Opm::Spline<Scalar> Spline;
    typedef std::vector<std::pair<Scalar, Scalar> > SplineSamplingPoints;

    typedef Opm::UniformXTabulated2DFunction<Scalar> TabulatedFunction;

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
     * After calling this method the reference densities, all formation
     * and volume factors, the oil bubble pressure, all viscosities
     * and the water compressibility must be set. Before the fluid
     * system can be used, initEnd() must be called to finalize the
     * initialization.
     */
    static void initBegin()
    {
        setWaterReferenceFormationFactor(1.0);
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

        auto& oilViscosity = oilViscosity_[regionIdx];
        auto& oilFormationVolumeFactor = oilFormationVolumeFactor_[regionIdx];
        auto& gasDissolutionFactorSpline = gasDissolutionFactorSpline_[regionIdx];

        gasDissolutionFactorSpline.setXYArrays(saturatedTable->numRows(),
                                               saturatedTable->getPressureColumn(),
                                               saturatedTable->getGasSolubilityColumn(),
                                               /*type=*/Spline::Monotonic);
        updateSaturationPressureSpline_(regionIdx);

        Scalar rhooRef = referenceDensity(oilPhaseIdx, regionIdx);
        Scalar rhogRef = referenceDensity(gasPhaseIdx, regionIdx);

        // extract the table for the oil formation factor
        for (int outerIdx = 0; outerIdx < saturatedTable->numRows(); ++ outerIdx) {
            Scalar Rs = saturatedTable->getGasSolubilityColumn()[outerIdx];

            Scalar XoG = Rs/(rhooRef/rhogRef + Rs);

            oilFormationVolumeFactor.appendXPos(XoG);
            oilViscosity.appendXPos(XoG);

            assert(oilFormationVolumeFactor.numX() == outerIdx + 1);
            assert(oilViscosity.numX() == outerIdx + 1);

            const auto underSaturatedTable = pvtoTable.getInnerTable(outerIdx);
            for (int innerIdx = underSaturatedTable->numRows() - 1; innerIdx >= 0; -- innerIdx) {
                Scalar po = underSaturatedTable->getPressureColumn()[innerIdx];
                Scalar Bo = underSaturatedTable->getOilFormationFactorColumn()[innerIdx];
                Scalar muo = underSaturatedTable->getOilViscosityColumn()[innerIdx];

                oilFormationVolumeFactor.appendSamplePoint(outerIdx, po, Bo);
                oilViscosity.appendSamplePoint(outerIdx, po, muo);
            }
        }

        // make sure to have at least two sample points per mole fraction
        for (int xIdx = 0; xIdx < oilFormationVolumeFactor.numX(); ++xIdx) {
            // a single sample point is definitely needed
            assert(oilFormationVolumeFactor.numY(xIdx) > 0);

            // everything is fine if the current table has two or more sampling points
            // for a given mole fraction
            if (oilFormationVolumeFactor.numY(xIdx) > 1)
                continue;

            // find the master table which will be used as a template to extend the
            // current line. We define master table as the first table which has values
            // for undersaturated oil...
            int masterTableIdx = xIdx + 1;
            for (; masterTableIdx < pvtoTable.getOuterTable()->numRows(); ++masterTableIdx) {
                if (pvtoTable.getInnerTable(masterTableIdx)->numRows() > 1)
                    break;
            }

            if (masterTableIdx >= pvtoTable.getOuterTable()->numRows())
                OPM_THROW(std::runtime_error,
                          "PVTO tables are invalid: The last table must exhibit at least one "
                          "entry for undersaturated oil!");

            // extend the current table using the master table. this is done by assuming
            // that the current table exhibits the same ratios of the oil formation
            // factors and viscosities for identical pressure rations as in the master
            // table.
            const auto masterTable = pvtoTable.getInnerTable(masterTableIdx);
            const auto curTable = pvtoTable.getInnerTable(xIdx);
            for (int newRowIdx = 1; newRowIdx < masterTable->numRows(); ++ newRowIdx) {
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

                oilFormationVolumeFactor.appendSamplePoint(xIdx, newPo, newBo);
                oilViscosity.appendSamplePoint(xIdx, newPo, newMuo);
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
    static void setPvtwTable(const PvtwTable &pvtwTable, int regionIdx=0)
    {
        assert(pvtwTable.numRows() > 0);

        resizeArrays_(regionIdx);

        // actually the PVTW does not specify a table, but for now we use a table wrapper
        // anyway because we don't want to break the opm-parser API at this point...
        std::vector<double> pressureCol = pvtwTable.getPressureColumn();
        std::vector<double> refBwCol = pvtwTable.getFormationFactorColumn();
        std::vector<double> compressCol = pvtwTable.getCompressibilityColumn();
        std::vector<double> viscosityCol = pvtwTable.getViscosityColumn();
        std::vector<double> viscosibilityCol = pvtwTable.getViscosibilityColumn();

        assert(pressureCol.size() == 1);

        waterReferencePressureScalar_[regionIdx] = pressureCol[0];
        waterReferenceFormationFactorScalar_[regionIdx] = refBwCol[0];
        waterCompressibilityScalar_[regionIdx] = compressCol[0];
        waterViscosityScalar_[regionIdx] = viscosityCol[0];
        waterViscosibilityScalar_[regionIdx] = viscosibilityCol[0];
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

        gasFormationVolumeFactorSpline_[regionIdx].setXYArrays(numSamples,
                                                               pvdgTable.getPressureColumn(),
                                                               pvdgTable.getFormationFactorColumn(),
                                                               /*type=*/Spline::Monotonic);

        gasViscositySpline_[regionIdx].setXYArrays(numSamples,
                                                   pvdgTable.getPressureColumn(),
                                                   pvdgTable.getViscosityColumn(),
                                                   /*type=*/Spline::Monotonic);
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
     * \brief Initialize the spline for the gas dissolution factor \f$R_s\f$
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setSaturatedOilGasDissolutionFactor(const SplineSamplingPoints &samplePoints,
                                                    int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        gasDissolutionFactorSpline_[regionIdx].setContainerOfTuples(samplePoints,
                                                                    /*type=*/Spline::Monotonic);
        assert(gasDissolutionFactorSpline_[regionIdx].monotonic());

        updateSaturationPressureSpline_(regionIdx);
    }

    /*!
     * \brief Initialize the spline for the oil formation volume factor
     *
     * The oil density/volume factor is a function of (p_o, X_o^G), but this method here
     * only requires the volume factor of gas-saturated oil (which only depends on
     * pressure) while the dependence on the gas mass fraction is estimated...
     */
    static void setSaturatedOilFormationVolumeFactor(const SplineSamplingPoints &samplePoints,
                                                     int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        auto& oilFormationVolumeFactor = oilFormationVolumeFactor_[regionIdx];

        auto &RsSpline = gasDissolutionFactorSpline_[regionIdx];

        Scalar XoGMin = 0.0;
        Scalar RsMax = RsSpline.eval(gasDissolutionFactorSpline_[regionIdx].xMax());
        Scalar XoGMax =
            RsMax*referenceDensity(gasPhaseIdx, regionIdx)
            / (RsMax*referenceDensity(gasPhaseIdx, regionIdx)
               + referenceDensity(oilPhaseIdx, regionIdx));

        Scalar poMin = samplePoints.front().first;
        Scalar poMax = samplePoints.back().first;

        size_t nX = 20;
        size_t nP = samplePoints.size()*2;

        Spline oilFormationVolumeFactorSpline;
        oilFormationVolumeFactorSpline.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);

        // calculate a table of estimated densities depending on pressure and gas mass
        // fraction
        for (size_t XIdx = 0; XIdx < nX; ++XIdx) {
            Scalar XoG = XoGMin + (XoGMax - XoGMin)*XIdx/nX;

            oilFormationVolumeFactor.appendXPos(XoG);

            for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
                Scalar po = poMin + (poMax - poMin)*pIdx/nP;

                Scalar poSat = oilSaturationPressure(XoG, regionIdx);
                Scalar BoSat = oilFormationVolumeFactorSpline.eval(poSat, /*extrapolate=*/true);
                Scalar drhoo_dp = (1.1200 - 1.1189)/((5000 - 4000)*6894.76);
                Scalar rhoo = referenceDensity(oilPhaseIdx, regionIdx)/BoSat*(1 + drhoo_dp*(po - poSat));

                Scalar Bo = referenceDensity(oilPhaseIdx, regionIdx)/rhoo;

                oilFormationVolumeFactor.appendSamplePoint(XIdx, po, Bo);
            }
        }
    }

    /*!
     * \brief Initialize the spline for the oil formation volume factor
     *
     * This is a function of (Rs, po)...
     */
    static void setOilFormationVolumeFactor(const TabulatedFunction &Bo, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        oilFormationVolumeFactor_[regionIdx] = Bo;
    }

    /*!
     * \brief Initialize the spline for the viscosity of gas-saturated
     *        oil.
     *
     * This is a function of (Rs, po)...
     */
    static void setOilViscosity(const TabulatedFunction &muo, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        oilViscosity_[regionIdx] = muo;
    }

    /*!
     * \brief Initialize the oil viscosity
     *
     * The oil viscosity is a function of (p_o, X_o^G), but this method here only requires
     * the viscosity of gas-saturated oil (which only depends on pressure) while the
     * dependence on the gas mass fraction is assumed to be zero...
     */
    static void setSaturatedOilViscosity(const SplineSamplingPoints &samplePoints, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        auto& gasDissolutionFactorSpline = gasDissolutionFactorSpline_[regionIdx];

        Scalar XoGMin = 0.0;
        Scalar RsMax = gasDissolutionFactorSpline.eval(gasDissolutionFactorSpline_[regionIdx].xMax());
        Scalar XoGMax =
            RsMax*referenceDensity(gasPhaseIdx, regionIdx)
            /(RsMax*referenceDensity(gasPhaseIdx, regionIdx)
              + referenceDensity(oilPhaseIdx, regionIdx));

        Scalar poMin = samplePoints.front().first;
        Scalar poMax = samplePoints.back().first;

        size_t nX = 20;
        size_t nP = samplePoints.size()*2;

        auto& oilViscosity = oilViscosity_[regionIdx];

        Spline muoSpline;
        muoSpline.setContainerOfTuples(samplePoints, /*type=*/Spline::Monotonic);

        // calculate a table of estimated densities depending on pressure and gas mass
        // fraction
        for (size_t XIdx = 0; XIdx < nX; ++XIdx) {
            Scalar XoG = XoGMin + (XoGMax - XoGMin)*XIdx/nX;

            oilViscosity.appendXPos(XoG);

            for (size_t pIdx = 0; pIdx < nP; ++pIdx) {
                Scalar po = poMin + (poMax - poMin)*pIdx/nP;
                Scalar muo = muoSpline.eval(po);

                oilViscosity.appendSamplePoint(XIdx, po, muo);
            }

        }
    }

    /*!
     * \brief Initialize the spline for the formation volume factor of dry gas
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setGasFormationVolumeFactor(const SplineSamplingPoints &samplePoints,
                                            int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        gasFormationVolumeFactorSpline_[regionIdx].setContainerOfTuples(samplePoints,
                                                                        /*type=*/Spline::Monotonic);
        assert(gasFormationVolumeFactorSpline_[regionIdx].monotonic());
    }

    /*!
     * \brief Initialize the spline for the viscosity of dry gas
     *
     * \param samplePoints A container of (x,y) values which is suitable to be passed to a spline.
     */
    static void setGasViscosity(const SplineSamplingPoints &samplePoints, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        gasViscositySpline_[regionIdx].setContainerOfTuples(samplePoints,
                                                            /*type=*/Spline::Monotonic);
        assert(gasViscositySpline_[regionIdx].monotonic());
    }

    /*!
     * \brief Set the water reference pressure [Pa]
     */
    static void setWaterReferencePressure(Scalar p, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        waterReferencePressureScalar_[regionIdx] = p;
    }

    /*!
     * \brief Set the water reference formation volume factor [-]
     */
    static void setWaterReferenceFormationFactor(Scalar BwRef, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        waterReferenceFormationFactorScalar_[regionIdx] = BwRef;
    }

    /*!
     * \brief Set the water viscosity [Pa s]
     */
    static void setWaterViscosity(Scalar muWater, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        waterViscosityScalar_[regionIdx] = muWater;
    }

    /*!
     * \brief Set the water "viscosibility" [1/ (Pa s)]
     */
    static void setWaterViscosibility(Scalar muCompWater, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        waterViscosibilityScalar_[regionIdx] = muCompWater;
    }

    /*!
     * \brief Set the water compressibility [1 / Pa]
     */
    static void setWaterCompressibility(Scalar cWater, int regionIdx=0)
    {
        resizeArrays_(regionIdx);

        waterCompressibilityScalar_[regionIdx] = cWater;
    }

    /*!
     * \brief Finish initializing the black oil fluid system.
     */
    static void initEnd()
    {
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
            // ATTENTION: XoG is represented by the _first_ axis!
            return oilViscosity_[regionIdx].eval(XoG, p);
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
        return oilFormationVolumeFactor_[regionIdx].eval(XoG, pressure);
    }

    /*!
     * \brief Return the formation volume factor of gas.
     */
    static Scalar gasFormationVolumeFactor(Scalar pressure, int regionIdx=0)
    {
        Scalar BgRaw = gasFormationVolumeFactorSpline_[regionIdx].eval(pressure, /*extrapolate=*/true);
        return BgRaw;
    }

    /*!
     * \brief Return the formation volume factor of water.
     */
    static Scalar waterFormationVolumeFactor(Scalar pressure, int regionIdx=0)
    {
        // cf. ECLiPSE 2011 technical description, p. 116
        Scalar pRef = waterReferencePressureScalar_[regionIdx];
        Scalar X = waterCompressibilityScalar_[regionIdx]*(pressure - pRef);

        Scalar BwRef = waterReferenceFormationFactorScalar_[regionIdx];

        // TODO (?): consider the salt concentration of the brine
        return BwRef/(1 + X*(1 + X/2));
    }

    /*!
     * \brief Returns the gas formation factor \f$R_s\f$ for a given pressure
     *
     * \param pressure The pressure of interest [Pa]
     */
    static Scalar gasDissolutionFactor(Scalar pressure, int regionIdx=0)
    { return gasDissolutionFactorSpline_[regionIdx].eval(pressure, /*extrapolate=*/true); }

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
            return 1e3*pvWater / pressure;
        else if (compIdx == gasCompIdx)
            return 1e6*pvWater / pressure;

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
        // assume an ideal gas
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
        // use the saturation pressure spline to get a pretty good
        // initial value
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
        // equivalent to the gas formation factor [m^3/m^3] at current pressure times the
        // gas density [kg/m^3] at standard pressure
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
     *
     * "Normalized" means "1 at surface pressure".
     */
    static Scalar oilFormationVolumeFactor(Scalar oilPressure, Scalar XoG, int regionIdx=0)
    {
        // ATTENTION: XoG is represented by the _first_ axis!
        Scalar BoRaw = oilFormationVolumeFactor_[regionIdx].eval(XoG, oilPressure);
        return BoRaw;
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

        // for some reason, the gas concentration is not considered in the oil formation
        // volume factor. WTF? While I have no idea why this should be there, it is
        // analogous to what's done in opm-autodiff which seems to be in agreement what
        // the commercial simulator does...
        Scalar Rs = XoG*rhooRef / ((1-XoG)*rhogRef);
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
     * \brief Return the density of dry gas.
     */
    static Scalar gasDensity(Scalar pressure, int regionIdx)
    {
        // gas formation volume factor at reservoir pressure
        Scalar Bg = gasFormationVolumeFactor(pressure, regionIdx);
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
        if (static_cast<int>(oilViscosity_.size()) <= regionIdx) {
            int numRegions = regionIdx + 1;
            oilFormationVolumeFactor_.resize(numRegions);
            gasFormationVolumeFactorSpline_.resize(numRegions);
            oilViscosity_.resize(numRegions);
            gasViscositySpline_.resize(numRegions);
            gasDissolutionFactorSpline_.resize(numRegions);
            saturationPressureSpline_.resize(numRegions);
            waterReferencePressureScalar_.resize(numRegions);
            waterReferenceFormationFactorScalar_.resize(numRegions);
            waterCompressibilityScalar_.resize(numRegions);
            waterViscosityScalar_.resize(numRegions);
            waterViscosibilityScalar_.resize(numRegions);

            referenceDensity_.resize(numRegions);
        }
    }

    static void updateSaturationPressureSpline_(int regionIdx)
    {
        resizeArrays_(regionIdx);

        auto& gasDissolutionFactorSpline = gasDissolutionFactorSpline_[regionIdx];

        // create the spline representing saturation pressure
        // depending of the mass fraction in gas
        int n = gasDissolutionFactorSpline.numSamples()*5;
        int delta =
            (gasDissolutionFactorSpline.xMax() - gasDissolutionFactorSpline.xMin())/(n + 1);

        SplineSamplingPoints pSatSamplePoints;
        Scalar X_oG = 0;
        for (int i=0; i <= n; ++ i) {
            Scalar pSat = gasDissolutionFactorSpline.xMin() + i*delta;
            X_oG = saturatedOilGasMassFraction(pSat, regionIdx);

            std::pair<Scalar, Scalar> val(X_oG, pSat);
            pSatSamplePoints.push_back(val);
        }
        saturationPressureSpline_[regionIdx].setContainerOfTuples(pSatSamplePoints,
                                                                  /*type=*/Spline::Monotonic);
    }

    static Scalar gasViscosity_(Scalar pressure, int regionIdx)
    { return gasViscositySpline_[regionIdx].eval(pressure, /*extrapolate=*/true); }

    static Scalar waterViscosity_(Scalar pressure, int regionIdx)
    {
        Scalar muRef = waterViscosityScalar_[regionIdx];
        Scalar Cnu = waterViscosibilityScalar_[regionIdx];
        Scalar pRef = waterReferencePressureScalar_[regionIdx];
        Scalar deltamu = (pressure - pRef)*Cnu*muRef;
        return muRef + deltamu;
    }

    static std::vector<TabulatedFunction> oilFormationVolumeFactor_;
    static std::vector<TabulatedFunction> oilViscosity_;
    static std::vector<Spline> gasDissolutionFactorSpline_;
    static std::vector<Spline> saturationPressureSpline_;

    static std::vector<Spline> gasFormationVolumeFactorSpline_;
    static std::vector<Spline> gasViscositySpline_;

    static std::vector<Scalar> waterReferencePressureScalar_;
    static std::vector<Scalar> waterReferenceFormationFactorScalar_;
    static std::vector<Scalar> waterCompressibilityScalar_;
    static std::vector<Scalar> waterViscosityScalar_;
    static std::vector<Scalar> waterViscosibilityScalar_;

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
std::vector<typename BlackOil<Scalar>::TabulatedFunction>
BlackOil<Scalar>::oilFormationVolumeFactor_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::TabulatedFunction>
BlackOil<Scalar>::oilViscosity_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::Spline>
BlackOil<Scalar>::gasDissolutionFactorSpline_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::Spline>
BlackOil<Scalar>::gasFormationVolumeFactorSpline_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::Spline>
BlackOil<Scalar>::saturationPressureSpline_;

template <class Scalar>
std::vector<typename BlackOil<Scalar>::Spline>
BlackOil<Scalar>::gasViscositySpline_;

template <class Scalar>
std::vector<Scalar>
BlackOil<Scalar>::waterReferencePressureScalar_;

template <class Scalar>
std::vector<Scalar>
BlackOil<Scalar>::waterReferenceFormationFactorScalar_;

template <class Scalar>
std::vector<Scalar>
BlackOil<Scalar>::waterCompressibilityScalar_;

template <class Scalar>
std::vector<Scalar>
BlackOil<Scalar>::waterViscosityScalar_;

template <class Scalar>
std::vector<Scalar>
BlackOil<Scalar>::waterViscosibilityScalar_;

template <class Scalar>
std::vector<std::array<Scalar, 3> >
BlackOil<Scalar>::referenceDensity_;

template <class Scalar>
Scalar
BlackOil<Scalar>::molarMass_[BlackOil<Scalar>::numComponents];
}} // namespace Opm, FluidSystems

#endif
