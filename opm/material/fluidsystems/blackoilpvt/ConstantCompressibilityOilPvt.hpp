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
 * \copydoc Opm::ConstantCompressibilityOilPvt
 */
#ifndef OPM_CONSTANT_COMPRESSIBILITY_OIL_PVT_HPP
#define OPM_CONSTANT_COMPRESSIBILITY_OIL_PVT_HPP

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
 * \brief This class represents the Pressure-Volume-Temperature relations of the oil phase
 *        without dissolved gas and constant compressibility/"viscosibility".
 */
template <class Scalar, class Evaluation = Scalar>
class ConstantCompressibilityOilPvt  : public OilPvtInterfaceTemplateWrapper<Scalar,
                                                                             Evaluation,
                                                                             ConstantCompressibilityOilPvt<Scalar, Evaluation> >
{
    friend class OilPvtInterfaceTemplateWrapper<Scalar,
                                                Evaluation,
                                                ConstantCompressibilityOilPvt<Scalar, Evaluation> >;

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
        oilReferencePressure_.resize(numRegions);
        oilReferenceFormationVolumeFactor_.resize(numRegions);
        oilCompressibility_.resize(numRegions);
        oilViscosity_.resize(numRegions);
        oilViscosibility_.resize(numRegions);

        for (int regionIdx = 0; regionIdx < numRegions; ++regionIdx) {
            setReferenceFormationVolumeFactor(regionIdx, 1.0);
            setReferencePressure(regionIdx, BlackOilFluidSystem::surfacePressure);
        }
    }

#if HAVE_OPM_PARSER
    /*!
     * \brief Sets the pressure-dependent oil viscosity and density
     *        using the Eclipse PVCDO keyword.
     */
    void setPvcdo(int regionIdx, DeckKeywordConstPtr pvcdoKeyword)
    {
        assert(static_cast<int>(pvcdoKeyword->size()) >= regionIdx);

        auto pvcdoRecord = pvcdoKeyword->getRecord(regionIdx);
        oilReferencePressure_[regionIdx] =
            pvcdoRecord->getItem("P_REF")->getSIDouble(0);
        oilReferenceFormationVolumeFactor_[regionIdx] =
            pvcdoRecord->getItem("OIL_VOL_FACTOR")->getSIDouble(0);
        oilCompressibility_[regionIdx] =
            pvcdoRecord->getItem("OIL_COMPRESSIBILITY")->getSIDouble(0);
        oilViscosity_[regionIdx] =
            pvcdoRecord->getItem("OIL_VISCOSITY")->getSIDouble(0);
        oilViscosibility_[regionIdx] =
            pvcdoRecord->getItem("OIL_VISCOSIBILITY")->getSIDouble(0);
    }
#endif

    /*!
     * \brief Set the viscosity and "viscosibility" of the oil phase.
     */
    void setViscosity(int regionIdx, Scalar muo, Scalar oilViscosibility = 0.0)
    {
        oilViscosity_[regionIdx] = muo;
        oilViscosibility_[regionIdx] = oilViscosibility;
    }

    /*!
     * \brief Set the compressibility of the oil phase.
     */
    void setCompressibility(int regionIdx, Scalar oilCompressibility)
    { oilCompressibility_[regionIdx] = oilCompressibility; }

    /*!
     * \brief Set the oil reference pressure [Pa]
     */
    void setReferencePressure(int regionIdx, Scalar p)
    { oilReferencePressure_[regionIdx] = p; }

    /*!
     * \brief Set the oil reference formation volume factor [-]
     */
    void setReferenceFormationVolumeFactor(int regionIdx, Scalar BoRef)
    { oilReferenceFormationVolumeFactor_[regionIdx] = BoRef; }

    /*!
     * \brief Set the oil "viscosibility" [1/ (Pa s)]
     */
    void setViscosibility(int regionIdx, Scalar muComp)
    { oilViscosibility_[regionIdx] = muComp; }

    /*!
     * \brief Finish initializing the oil phase PVT properties.
     */
    void initEnd()
    { }

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
        // Eclipse calculates the viscosity in a weird way: it
        // calcultes the product of B_w and mu_w and then divides the
        // result by B_w...
        Scalar BoMuoRef = oilViscosity_[regionIdx]*oilReferenceFormationVolumeFactor_[regionIdx];
        const LhsEval& Bo = formationVolumeFactor_(regionIdx, temperature, pressure, XoG);

        Scalar pRef = oilReferencePressure_[regionIdx];
        const LhsEval& Y =
            (oilCompressibility_[regionIdx] - oilViscosibility_[regionIdx])
            * (pressure - pRef);
        return BoMuoRef/((1 + Y*(1 + Y/2))*Bo);
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
        const LhsEval& Bo = formationVolumeFactor_(regionIdx, temperature, pressure, XoG);
        Scalar rhooRef = BlackOilFluidSystem::referenceDensity(oilPhaseIdx, regionIdx);
        return rhooRef/Bo;
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
        // cf. ECLiPSE 2011 technical description, p. 116
        Scalar pRef = oilReferencePressure_[regionIdx];
        const LhsEval& X = oilCompressibility_[regionIdx]*(pressure - pRef);

        Scalar BoRef = oilReferenceFormationVolumeFactor_[regionIdx];
        return BoRef/(1 + X*(1 + X/2));
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
        else if (compIdx == BlackOilFluidSystem::oilCompIdx)
            // assume that the affinity of the oil component to the
            // oil phase is one million times smaller than that of the
            // oil component
            return 1e8*phi_oO;

        assert(compIdx == BlackOilFluidSystem::gasCompIdx);
        // gas is immiscible with dead oil as well...
        return 1.01e8*phi_oO;
    }

    /*!
     * \brief Returns the gas dissolution factor \f$R_s\f$ [m^3/m^3] of the oil phase.
     */
    template <class LhsEval>
    LhsEval gasDissolutionFactor_(int regionIdx,
                                  const LhsEval& temperature,
                                  const LhsEval& pressure) const
    {
        typedef Opm::MathToolbox<LhsEval> Toolbox;

        return Toolbox::createConstant(0.0); /* this is dead oil! */
    }

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

        return Toolbox::createConstant(0.0); /* this is dead oil, so there isn't any meaningful saturation pressure! */
    }

    template <class LhsEval>
    LhsEval saturatedOilGasMassFraction_(int regionIdx,
                                         const LhsEval& temperature,
                                         const LhsEval& pressure) const
    {
        typedef Opm::MathToolbox<LhsEval> Toolbox;

        return Toolbox::createConstant(0.0); /* this is dead oil! */
    }

    template <class LhsEval>
    LhsEval saturatedOilGasMoleFraction_(int regionIdx,
                                         const LhsEval& temperature,
                                         const LhsEval& pressure) const
    {
        typedef Opm::MathToolbox<LhsEval> Toolbox;

        return Toolbox::createConstant(0.0); /* this is dead oil! */
    }

private:
    std::vector<Scalar> oilReferencePressure_;
    std::vector<Scalar> oilReferenceFormationVolumeFactor_;
    std::vector<Scalar> oilCompressibility_;
    std::vector<Scalar> oilViscosity_;
    std::vector<Scalar> oilViscosibility_;
};

} // namespace Opm

#endif
