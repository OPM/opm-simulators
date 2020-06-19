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
 *
 * \copydoc Opm::RichardsPrimaryVariables
 */
#ifndef EWOMS_RICHARDS_PRIMARY_VARIABLES_HH
#define EWOMS_RICHARDS_PRIMARY_VARIABLES_HH

#include "richardsproperties.hh"

#include <opm/models/discretization/common/fvbaseprimaryvariables.hh>

#include <opm/material/constraintsolvers/ImmiscibleFlash.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/common/fvector.hh>

namespace Opm {

/*!
 * \ingroup RichardsModel
 *
 * \brief Represents the primary variables used in the Richards model.
 *
 * This class is basically a Dune::FieldVector which can retrieve its
 * contents from an aribitatry fluid state.
 */
template <class TypeTag>
class RichardsPrimaryVariables : public FvBasePrimaryVariables<TypeTag>
{
    using ParentType = FvBasePrimaryVariables<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using EnergyModule = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    // primary variable indices
    enum { pressureWIdx = Indices::pressureWIdx };

    enum { liquidPhaseIdx = getPropValue<TypeTag, Properties::LiquidPhaseIndex>() };
    enum { gasPhaseIdx = getPropValue<TypeTag, Properties::GasPhaseIndex>() };

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

    using ComponentVector = Dune::FieldVector<Scalar, numComponents>;
    using PhaseVector = Dune::FieldVector<Scalar, numPhases>;
    using Toolbox = typename Opm::MathToolbox<Evaluation>;
    using ImmiscibleFlash = Opm::ImmiscibleFlash<Scalar, FluidSystem>;

public:
    RichardsPrimaryVariables() : ParentType()
    { Opm::Valgrind::SetUndefined(*this); }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    RichardsPrimaryVariables(Scalar value) : ParentType(value)
    {}

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const
     * ImmisciblePrimaryVariables& )
     */
    RichardsPrimaryVariables(const RichardsPrimaryVariables& value) = default;
    RichardsPrimaryVariables& operator=(const RichardsPrimaryVariables& value) = default;

    /*!
     * \brief Set the primary variables with the wetting phase
     *        pressure, saturation and temperature.
     *
     * \param T The temperature [K]
     * \param pw The pressure of the wetting phase [Pa]
     * \param Sw The saturation of the wetting phase []
     * \param matParams The capillary pressure law parameters
     */
    void assignImmiscibleFromWetting(Scalar T, Scalar pw, Scalar Sw,
                                     const MaterialLawParams& matParams)
    {
        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;

        fs.setTemperature(T);
        fs.setSaturation(liquidPhaseIdx, Sw);
        fs.setSaturation(gasPhaseIdx, 1 - Sw);

        // set phase pressures
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(liquidPhaseIdx, pw);
        fs.setPressure(gasPhaseIdx, pw + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));

        assignNaive(fs);
    }

    /*!
     * \brief Set the primary variables with the non-wetting phase
     *        pressure, saturation and temperature.
     *
     * \param T The temperature [K]
     * \param pn The pressure of the non-wetting phase [Pa]
     * \param Sn The saturation of the non-wetting phase []
     * \param matParams The capillary pressure law parameters
     */
    void assignImmiscibleFromNonWetting(Scalar T, Scalar pn, Scalar Sn,
                                        const MaterialLawParams& matParams)
    {
        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fs;

        fs.setTemperature(T);
        fs.setSaturation(liquidPhaseIdx, 1 - Sn);
        fs.setSaturation(gasPhaseIdx, Sn);

        // set phase pressures
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(gasPhaseIdx, pn);
        fs.setPressure(gasPhaseIdx, pn + (pC[liquidPhaseIdx] - pC[gasPhaseIdx]));

        assignNaive(fs);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignMassConservative
     */
    template <class FluidState>
    void assignMassConservative(const FluidState& fluidState,
                                const MaterialLawParams& matParams,
                                bool isInEquilibrium OPM_UNUSED= false)
    {
        ComponentVector globalMolarities(0.0);
        for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                globalMolarities[compIdx] +=
                    fluidState.molarity(phaseIdx, compIdx) * fluidState.saturation(phaseIdx);
            }
        }

        Opm::ImmiscibleFluidState<Scalar, FluidSystem> fsFlash;
        fsFlash.assign(fluidState);
        typename FluidSystem::ParameterCache paramCache;
        ImmiscibleFlash::template solve<MaterialLaw>(fsFlash, paramCache,
                                                     matParams,
                                                     globalMolarities);

        assignNaive(fsFlash);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState)
    {
        // assign the phase temperatures. this is out-sourced to
        // the energy module
        EnergyModule::setPriVarTemperatures(*this, fluidState);

        (*this)[pressureWIdx] = fluidState.pressure(liquidPhaseIdx);
    }
};

} // namespace Opm

#endif
