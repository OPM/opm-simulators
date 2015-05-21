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
 * \copydoc Opm::DryGasPvt
 */
#ifndef OPM_DRY_GAS_PVT_HPP
#define OPM_DRY_GAS_PVT_HPP

#include "GasPvtInterface.hpp"

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/material/common/OpmFinal.hpp>
#include <opm/material/common/Tabulated1DFunction.hpp>

#if HAVE_OPM_PARSER
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#endif

#include <vector>

namespace Opm {

/*!
 * \brief This class represents the Pressure-Volume-Temperature relations of the gas phase
 *        without vaporized oil.
 */
template <class Scalar, class Evaluation = Scalar>
class DryGasPvt
    : public GasPvtInterfaceTemplateWrapper<Scalar, Evaluation, DryGasPvt<Scalar, Evaluation> >
{
    friend class GasPvtInterfaceTemplateWrapper<Scalar, Evaluation, DryGasPvt<Scalar, Evaluation> >;

    typedef FluidSystems::BlackOil<Scalar, Evaluation> BlackOilFluidSystem;

    typedef Opm::Tabulated1DFunction<Scalar> TabulatedOneDFunction;
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
    }

#if HAVE_OPM_PARSER
    void setPvdgTable(int regionIdx, const PvdgTable& pvdgTable)
    {
        int numSamples = pvdgTable.numRows();
        assert(numSamples > 1);

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
     * \brief Initialize the viscosity of the gas phase.
     *
     * This is a function of \f$(p_g)\f$...
     */
    void setGasViscosity(int regionIdx, const TabulatedOneDFunction& mug)
    { gasMu_[regionIdx] = mug; }

    /*!
     * \brief Initialize the function for the formation volume factor of dry gas
     *
     * \param samplePoints A container of \f$(p_g, B_g)\f$ values
     */
    void setGasFormationVolumeFactor(int regionIdx, const SamplingPoints &samplePoints)
    {
        SamplingPoints tmp(samplePoints);
        auto it = tmp.begin();
        const auto& endIt = tmp.end();
        for (; it != endIt; ++ it)
            std::get<1>(*it) = 1.0/std::get<1>(*it);

        inverseGasB_[regionIdx].setContainerOfTuples(tmp);
        assert(inverseGasB_[regionIdx].monotonic());
    }

    /*!
     * \brief Finish initializing the oil phase PVT properties.
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
            assert(gasMu.numSamples() == invGasB.numSamples());

            std::vector<Scalar> pressureValues(gasMu.numSamples());
            std::vector<Scalar> invGasBMuValues(gasMu.numSamples());
            for (int pIdx = 0; pIdx < gasMu.numSamples(); ++pIdx) {
                pressureValues[pIdx] = invGasB.xAt(pIdx);
                invGasBMuValues[pIdx] = invGasB.valueAt(pIdx) * (1.0/gasMu.valueAt(pIdx));
            }

            inverseGasBMu_[regionIdx].setXYContainers(pressureValues, invGasBMuValues);
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
        const LhsEval& invBg = inverseGasB_[regionIdx].eval(pressure, /*extrapolate=*/true);
        const LhsEval& invMugBg = inverseGasBMu_[regionIdx].eval(pressure, /*extrapolate=*/true);

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
        // gas formation volume factor at reservoir pressure
        const LhsEval& Bg = formationVolumeFactor_(regionIdx, temperature, pressure, XgO);
        return BlackOilFluidSystem::referenceDensity(gasPhaseIdx, regionIdx)/Bg;
    }

    /*!
     * \brief Returns the formation volume factor [-] of the fluid phase.
     */
    template <class LhsEval>
    LhsEval formationVolumeFactor_(int regionIdx,
                                   const LhsEval& temperature,
                                   const LhsEval& pressure,
                                   const LhsEval& XgO) const
    { return 1.0/inverseGasB_[regionIdx].eval(pressure, /*extrapolate=*/true); }

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

        // make the gas component more affine to the gas phase than the other components
        if (compIdx == BlackOilFluidSystem::gasCompIdx)
            return Toolbox::createConstant(1.0);
        return Toolbox::createConstant(1e6);
    }

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
        return Toolbox::createConstant(0.0);  // this is dry gas!
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    template <class LhsEval>
    LhsEval oilVaporizationFactor_(int regionIdx,
                                   const LhsEval& temperature,
                                   const LhsEval& pressure) const
    {
        typedef Opm::MathToolbox<LhsEval> Toolbox;
        return Toolbox::createConstant(0.0);  // this is dry gas!
    }

    template <class LhsEval>
    LhsEval saturatedGasOilMassFraction_(int regionIdx,
                                         const LhsEval& temperature,
                                         const LhsEval& pressure) const
    {
        typedef Opm::MathToolbox<LhsEval> Toolbox;
        return Toolbox::createConstant(0.0);  // this is dry gas!
    }

    template <class LhsEval>
    LhsEval saturatedGasOilMoleFraction_(int regionIdx,
                                         const LhsEval& temperature,
                                         const LhsEval& pressure) const
    {
        typedef Opm::MathToolbox<LhsEval> Toolbox;
        return Toolbox::createConstant(0.0);  // this is dry gas!
    }

private:
    std::vector<TabulatedOneDFunction> inverseGasB_;
    std::vector<TabulatedOneDFunction> gasMu_;
    std::vector<TabulatedOneDFunction> inverseGasBMu_;
};

} // namespace Opm

#endif
