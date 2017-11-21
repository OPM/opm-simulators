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

#include <ewoms/common/propertysystem.hh>

#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>

// the ordering of these includes matters. do not touch it if you're not prepared to deal
// with some trouble!
#include <dune/grid/cpgrid/GridHelpers.hpp>
#include <opm/core/simulator/initStateEquil.hpp>


#include <vector>

namespace Ewoms {
namespace Properties {
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(GridView);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(MaterialLaw);
NEW_PROP_TAG(EnableSwatinit);
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

    typedef Opm::CompositionalFluidState<Scalar, FluidSystem> ScalarFluidState;

    enum { numPhases = FluidSystem::numPhases };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    enum { numComponents = FluidSystem::numComponents };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };

    enum { dimWorld = GridView::dimensionworld };

public:
    template <class EclMaterialLawManager>
    EclEquilInitializer(const Simulator& simulator,
                        EclMaterialLawManager& materialLawManager,
                        bool enableSwatinit)
        : simulator_(simulator)
    {
        const auto& gridManager = simulator.gridManager();
        const auto& deck = gridManager.deck();
        const auto& eclState = gridManager.eclState();
        const auto& equilGrid = gridManager.equilGrid();

        unsigned numElems = gridManager.grid().size(0);
        unsigned numEquilElems = gridManager.equilGrid().size(0);
        unsigned numCartesianElems = gridManager.cartesianSize();
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typedef Opm::ThreePhaseMaterialTraits<double,
                                              /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                              /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                              /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx> EquilTraits;

        // create a separate instance of the material law manager just because opm-core
        // only supports double as the type for scalars (but ebos may use float or quad)
        std::vector<int> compressedToCartesianEquilElemIdx(numEquilElems);
        std::vector<int> equilCartesianToCompressed( gridManager.equilCartesianSize(), -1 );

        for (unsigned equilElemIdx = 0; equilElemIdx < numEquilElems; ++equilElemIdx)
        {
            unsigned int equilCartesianIdx = gridManager.equilCartesianIndex(equilElemIdx);
            compressedToCartesianEquilElemIdx[equilElemIdx] = equilCartesianIdx;
            equilCartesianToCompressed[ equilCartesianIdx ] = equilElemIdx;
        }

        Opm::EclMaterialLawManager<EquilTraits> equilMaterialLawManager =
            Opm::EclMaterialLawManager<EquilTraits>();
        equilMaterialLawManager.initFromDeck(deck, eclState, compressedToCartesianEquilElemIdx);

        Opm::EQUIL::DeckDependent::InitialStateComputer<FluidSystem> initialState(equilMaterialLawManager,
                                                                                  gridManager.eclState(),
                                                                                  equilGrid,
                                                                                  simulator.problem().gravity()[dimWorld - 1],
                                                                                  enableSwatinit);


        std::vector<int> localToEquilIndex( numElems, -1 );
        for( unsigned int elemIdx = 0; elemIdx < numElems; ++elemIdx )
        {
            const int cartesianIndex = gridManager.cartesianIndex( elemIdx );
            assert( equilCartesianToCompressed[ cartesianIndex ] >= 0 );
            localToEquilIndex[ elemIdx ] = equilCartesianToCompressed[ cartesianIndex ];
        }

        // copy the result into the array of initial fluid states
        initialFluidStates_.resize(numCartesianElems);
        for (unsigned int elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            unsigned cartesianElemIdx = gridManager.cartesianIndex(elemIdx);
            auto& fluidState = initialFluidStates_[cartesianElemIdx];

            const unsigned int equilElemIdx = localToEquilIndex[ elemIdx ];

            // get the PVT region index of the current element
            unsigned regionIdx = simulator_.problem().pvtRegionIndex(elemIdx);

            // set the phase saturations
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                Scalar S;
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    S = 0.0;
                else {
                    S = initialState.saturation()[phaseIdx][equilElemIdx];
                }
                fluidState.setSaturation(phaseIdx, S);
            }

            // set the temperature
            Scalar T = FluidSystem::surfaceTemperature;
            fluidState.setTemperature(T);

            // set the phase pressures.
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                fluidState.setPressure(phaseIdx, initialState.press()[phaseIdx][equilElemIdx]);

            // reset the phase compositions
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                for (unsigned compIdx = 0; compIdx < numComponents; ++compIdx)
                    fluidState.setMoleFraction(phaseIdx, compIdx, 0.0);

            // the composition of the water phase is simple: it only consists of the
            // water component.
            fluidState.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);

            if (FluidSystem::enableDissolvedGas()) {
                // for gas and oil we have to translate surface volumes to mole fractions
                // before we can set the composition in the fluid state
                Scalar Rs = initialState.rs()[equilElemIdx];
                Scalar RsSat = FluidSystem::saturatedDissolutionFactor(fluidState, oilPhaseIdx, regionIdx);

                if (Rs > RsSat)
                    Rs = RsSat;

                // convert the Rs factor to mole fraction dissolved gas in oil
                Scalar XoG = FluidSystem::convertRsToXoG(Rs, regionIdx);
                Scalar xoG = FluidSystem::convertXoGToxoG(XoG, regionIdx);

                fluidState.setMoleFraction(oilPhaseIdx, oilCompIdx, 1 - xoG);
                fluidState.setMoleFraction(oilPhaseIdx, gasCompIdx, xoG);
            }

            // retrieve the surface volume of vaporized gas
            if (FluidSystem::enableVaporizedOil()) {
                Scalar Rv = initialState.rv()[equilElemIdx];
                Scalar RvSat = FluidSystem::saturatedDissolutionFactor(fluidState, gasPhaseIdx, regionIdx);

                if (Rv > RvSat)
                    Rv = RvSat;

                // convert the Rs factor to mole fraction dissolved gas in oil
                Scalar XgO = FluidSystem::convertRvToXgO(Rv, regionIdx);
                Scalar xgO = FluidSystem::convertXgOToxgO(XgO, regionIdx);

                fluidState.setMoleFraction(gasPhaseIdx, oilCompIdx, xgO);
                fluidState.setMoleFraction(gasPhaseIdx, gasCompIdx, 1 - xgO);
            }

            // deal with the changed pressure scaling due to SWATINIT if SWATINIT is
            // requested to be applied. this is quite hacky but hey it works!
            if (enableSwatinit) {
                const auto& equilScalingPoints =
                    equilMaterialLawManager.oilWaterScaledEpsPointsDrainage(equilElemIdx);
                auto& scalingPoints =
                    materialLawManager.oilWaterScaledEpsPointsDrainage(elemIdx);
                scalingPoints.setMaxPcnw(equilScalingPoints.maxPcnw());
            }
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
        const auto& gridManager = simulator_.gridManager();

        unsigned cartesianElemIdx = gridManager.cartesianIndex(elemIdx);
        return initialFluidStates_[cartesianElemIdx];
    }

protected:
    const Simulator& simulator_;

    std::vector<ScalarFluidState> initialFluidStates_;
};
} // namespace Ewoms

#endif
