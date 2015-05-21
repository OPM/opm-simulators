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
 * \copydoc Opm::DeadOilPvt
 */
#ifndef OPM_DEAD_OIL_PVT_HPP
#define OPM_DEAD_OIL_PVT_HPP

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
 *        without dissolved gas.
 */
template <class Scalar, class Evaluation = Scalar>
class DeadOilPvt : public OilPvtInterfaceTemplateWrapper<Scalar,
                                                         Evaluation,
                                                         DeadOilPvt<Scalar, Evaluation> >
{
    friend class OilPvtInterfaceTemplateWrapper<Scalar,
                                                Evaluation,
                                                DeadOilPvt<Scalar, Evaluation> >;

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
        if (static_cast<int>(inverseOilB_.size()) < numRegions) {
            inverseOilB_.resize(numRegions);
            inverseOilBMu_.resize(numRegions);
            oilMu_.resize(numRegions);
        }
    }

    /*!
     * \brief Initialize the function for the oil formation volume factor
     *
     * The oil formation volume factor \f$B_o\f$ is a function of \f$(p_o, X_o^G)\f$ and
     * represents the partial density of the oil component in the oil phase at a given
     * pressure.
     *
     * This method sets \f$1/B_o(p_o)\f$. Note that the mass fraction of the gas
     * component in the oil phase is missing when assuming dead oil.
     */
    void setInverseOilFormationVolumeFactor(int regionIdx, const TabulatedOneDFunction& invBo)
    { inverseOilB_[regionIdx] = invBo; }

    /*!
     * \brief Initialize the viscosity of the oil phase.
     *
     * This is a function of \f$(R_s, p_o)\f$...
     */
    void setOilViscosity(int regionIdx, const TabulatedOneDFunction& muo)
    { oilMu_[regionIdx] = muo; }


#if HAVE_OPM_PARSER
    /*!
     * \brief Initialize the oil parameters via the data specified by the PVDO ECL keyword.
     */
    void setPvdoTable(int regionIdx, const PvdoTable &pvdoTable)
    {
        assert(pvdoTable.numRows() > 1);

        const auto& BColumn(pvdoTable.getFormationFactorColumn());
        std::vector<Scalar> invBColumn(BColumn.size());
        for (unsigned i = 0; i < invBColumn.size(); ++i)
            invBColumn[i] = 1/BColumn[i];

        inverseOilB_[regionIdx].setXYArrays(pvdoTable.numRows(),
                                            pvdoTable.getPressureColumn(),
                                            invBColumn);
        oilMu_[regionIdx].setXYArrays(pvdoTable.numRows(),
                                      pvdoTable.getPressureColumn(),
                                      pvdoTable.getViscosityColumn());

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
            assert(oilMu.numSamples() == invOilB.numSamples());

            std::vector<Scalar> invBMuColumn;
            std::vector<Scalar> pressureColumn;
            invBMuColumn.resize(oilMu.numSamples());
            pressureColumn.resize(oilMu.numSamples());

            for (int pIdx = 0; pIdx < oilMu.numSamples(); ++pIdx) {
                pressureColumn[pIdx] = invOilB.xAt(pIdx);
                invBMuColumn[pIdx] = invOilB.valueAt(pIdx)*1/oilMu.valueAt(pIdx);
            }

            inverseOilBMu_[regionIdx].setXYArrays(pressureColumn.size(),
                                                  pressureColumn,
                                                  invBMuColumn);
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
        const LhsEval& invBo = inverseOilB_[regionIdx].eval(pressure, /*extrapolate=*/true);
        const LhsEval& invMuoBo = inverseOilBMu_[regionIdx].eval(pressure, /*extrapolate=*/true);

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

        const LhsEval& Bo = formationVolumeFactor_(regionIdx, temperature, pressure, XoG);
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
    { return 1.0 / inverseOilB_[regionIdx].eval(pressure, /*extrapolate=*/true); }

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
    std::vector<TabulatedOneDFunction> inverseOilB_;
    std::vector<TabulatedOneDFunction> oilMu_;
    std::vector<TabulatedOneDFunction> inverseOilBMu_;
};

} // namespace Opm

#endif
