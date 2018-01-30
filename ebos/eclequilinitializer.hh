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
/**
 * \file
 *
 * \copydoc Ewoms::EclEquilInitializer
 */
#ifndef EWOMS_ECL_EQUIL_INITIALIZER_HH
#define EWOMS_ECL_EQUIL_INITIALIZER_HH

#include "equil/initstateequil.hh"

#include <ewoms/common/propertysystem.hh>

#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

#include <vector>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(MaterialLaw);
NEW_PROP_TAG(EnableEnergy);
}

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Computes the initial condition based on the EQUIL keyword from ECL.
 *
 * So far, it uses the "initStateEquil()" function from opm-core. Since this method is
 * very much glued into the opm-core data structures, it should be reimplemented in the
 * medium to long term for some significant memory savings and less significant
 * performance improvements.
 */
template <class TypeTag>
class EclEquilInitializer
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

    enum { numPhases = FluidSystem::numPhases };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    enum { numComponents = FluidSystem::numComponents };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };

    enum { dimWorld = GridView::dimensionworld };
    enum { enableTemperature = GET_PROP_VALUE(TypeTag, EnableTemperature) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    // NB: setting the enableEnergy argument to true enables storage of enthalpy and
    // internal energy!
    typedef Opm::BlackOilFluidState<Scalar,
                                    FluidSystem,
                                    enableTemperature,
                                    enableEnergy> ScalarFluidState;

public:
    template <class EclMaterialLawManager>
    EclEquilInitializer(const Simulator& simulator,
                        EclMaterialLawManager& materialLawManager)
        : simulator_(simulator)
    {
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        unsigned numElems = vanguard.grid().size(0);
        unsigned numCartesianElems = vanguard.cartesianSize();
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

        EQUIL::DeckDependent::InitialStateComputer<TypeTag> initialState(materialLawManager,
                                                                         eclState,
                                                                         vanguard.grid(),
                                                                         simulator.problem().gravity()[dimWorld - 1]);

        // copy the result into the array of initial fluid states
        initialFluidStates_.resize(numCartesianElems);
        for (unsigned int elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            unsigned cartesianElemIdx = vanguard.cartesianIndex(elemIdx);
            auto& fluidState = initialFluidStates_[cartesianElemIdx];

            // get the PVT region index of the current element
            unsigned regionIdx = simulator_.problem().pvtRegionIndex(elemIdx);
            fluidState.setPvtRegionIndex(regionIdx);

            // set the phase saturations
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (FluidSystem::phaseIsActive(phaseIdx))
                    fluidState.setSaturation(phaseIdx, initialState.saturation()[phaseIdx][elemIdx]);
                else
                    fluidState.setSaturation(phaseIdx, 0.0);
            }

            if (FluidSystem::enableDissolvedGas())
                fluidState.setRs(initialState.rs()[elemIdx]);
            else
                fluidState.setRs(0.0);

            if (FluidSystem::enableVaporizedOil())
                fluidState.setRv(initialState.rv()[elemIdx]);
            else
                fluidState.setRv(0.0);


            // set the temperature.
            if (enableTemperature || enableEnergy)
                fluidState.setTemperature(initialState.temperature()[elemIdx]);

            // set the phase pressures.
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fluidState.setPressure(phaseIdx, initialState.press()[phaseIdx][elemIdx]);
        }
    }

    /*!
     * \brief Return the initial thermodynamic state which should be used as the initial
     *        condition.
     *
     * This is supposed to correspond to hydrostatic conditions.
     */
    const ScalarFluidState& initialFluidState(unsigned elemIdx) const
    {
        const auto& vanguard = simulator_.vanguard();

        unsigned cartesianElemIdx = vanguard.cartesianIndex(elemIdx);
        return initialFluidStates_[cartesianElemIdx];
    }

protected:
    const Simulator& simulator_;

    std::vector<ScalarFluidState> initialFluidStates_;
};
} // namespace Ewoms

#endif
