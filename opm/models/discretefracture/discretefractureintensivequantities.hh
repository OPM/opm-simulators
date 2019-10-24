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
 * \copydoc Opm::DiscreteFractureIntensiveQuantities
 */
#ifndef EWOMS_DISCRETE_FRACTURE_INTENSIVE_QUANTITIES_HH
#define EWOMS_DISCRETE_FRACTURE_INTENSIVE_QUANTITIES_HH

#include "discretefractureproperties.hh"

#include <opm/models/immiscible/immiscibleintensivequantities.hh>

#include <opm/material/common/Valgrind.hpp>

namespace Opm {

/*!
 * \ingroup DiscreteFractureModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the discret fracture immiscible multi-phase
 *        model.
 */
template <class TypeTag>
class DiscreteFractureIntensiveQuantities : public ImmiscibleIntensiveQuantities<TypeTag>
{
    typedef ImmiscibleIntensiveQuantities<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    enum { numPhases = FluidSystem::numPhases };
    enum { dimWorld = GridView::dimensionworld };

    static_assert(dimWorld == 2, "The fracture module currently is only "
                                 "implemented for the 2D case!");
    static_assert(numPhases == 2, "The fracture module currently is only "
                                  "implemented for two fluid phases!");

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { wettingPhaseIdx = MaterialLaw::wettingPhaseIdx };
    enum { nonWettingPhaseIdx = MaterialLaw::nonWettingPhaseIdx };
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef Opm::ImmiscibleFluidState<Scalar, FluidSystem,
                                      /*storeEnthalpy=*/enableEnergy> FluidState;

public:
    DiscreteFractureIntensiveQuantities()
    { }

    DiscreteFractureIntensiveQuantities(const DiscreteFractureIntensiveQuantities& other) = default;

    DiscreteFractureIntensiveQuantities& operator=(const DiscreteFractureIntensiveQuantities& other) = default;

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned vertexIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, vertexIdx, timeIdx);

        const auto& problem = elemCtx.problem();
        const auto& fractureMapper = problem.fractureMapper();
        unsigned globalVertexIdx = elemCtx.globalSpaceIndex(vertexIdx, timeIdx);

        Opm::Valgrind::SetUndefined(fractureFluidState_);
        Opm::Valgrind::SetUndefined(fractureVolume_);
        Opm::Valgrind::SetUndefined(fracturePorosity_);
        Opm::Valgrind::SetUndefined(fractureIntrinsicPermeability_);
        Opm::Valgrind::SetUndefined(fractureRelativePermeabilities_);

        // do nothing if there is no fracture within the current degree of freedom
        if (!fractureMapper.isFractureVertex(globalVertexIdx)) {
            fractureVolume_ = 0;
            return;
        }

        // Make sure that the wetting saturation in the matrix fluid
        // state does not get larger than 1
        Scalar SwMatrix =
            std::min<Scalar>(1.0, this->fluidState_.saturation(wettingPhaseIdx));
        this->fluidState_.setSaturation(wettingPhaseIdx, SwMatrix);
        this->fluidState_.setSaturation(nonWettingPhaseIdx, 1 - SwMatrix);

        // retrieve the facture width and intrinsic permeability from
        // the problem
        fracturePorosity_ =
            problem.fracturePorosity(elemCtx, vertexIdx, timeIdx);
        fractureIntrinsicPermeability_ =
            problem.fractureIntrinsicPermeability(elemCtx, vertexIdx, timeIdx);

        // compute the fracture volume for the current sub-control
        // volume. note, that we don't take overlaps of fractures into
        // account for this.
        fractureVolume_ = 0;
        const auto& vertexPos = elemCtx.pos(vertexIdx, timeIdx);
        for (unsigned vertex2Idx = 0; vertex2Idx < elemCtx.numDof(/*timeIdx=*/0); ++ vertex2Idx) {
            unsigned globalVertex2Idx = elemCtx.globalSpaceIndex(vertex2Idx, timeIdx);

            if (vertexIdx == vertex2Idx ||
                !fractureMapper.isFractureEdge(globalVertexIdx, globalVertex2Idx))
                continue;

            Scalar fractureWidth =
                problem.fractureWidth(elemCtx, vertexIdx, vertex2Idx, timeIdx);

            auto distVec = elemCtx.pos(vertex2Idx, timeIdx);
            distVec -= vertexPos;

            Scalar edgeLength = distVec.two_norm();

            // the fracture is always adjacent to two sub-control
            // volumes of the control volume, so when calculating the
            // volume of the fracture which gets attributed to one
            // SCV, the fracture width needs to divided by 2. Also,
            // only half of the edge is located in the current control
            // volume, so its length also needs to divided by 2.
            fractureVolume_ += (fractureWidth / 2) * (edgeLength / 2);
        }

        //////////
        // set the fluid state for the fracture.
        //////////

        // start with the same fluid state as in the matrix. This
        // implies equal saturations, pressures, temperatures,
        // enthalpies, etc.
        fractureFluidState_.assign(this->fluidState_);

        // ask the problem for the material law parameters of the
        // fracture.
        const auto& fractureMatParams =
            problem.fractureMaterialLawParams(elemCtx, vertexIdx, timeIdx);

        // calculate the fracture saturations which would be required
        // to be consistent with the pressures
        Scalar saturations[numPhases];
        MaterialLaw::saturations(saturations, fractureMatParams,
                                 fractureFluidState_);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            fractureFluidState_.setSaturation(phaseIdx, saturations[phaseIdx]);

        // Make sure that the wetting saturation in the fracture does
        // not get negative
        Scalar SwFracture =
            std::max<Scalar>(0.0, fractureFluidState_.saturation(wettingPhaseIdx));
        fractureFluidState_.setSaturation(wettingPhaseIdx, SwFracture);
        fractureFluidState_.setSaturation(nonWettingPhaseIdx, 1 - SwFracture);

        // calculate the relative permeabilities of the fracture
        MaterialLaw::relativePermeabilities(fractureRelativePermeabilities_,
                                            fractureMatParams,
                                            fractureFluidState_);

        // make sure that valgrind complains if the fluid state is not
        // fully defined.
        fractureFluidState_.checkDefined();
    }

public:
    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar fractureRelativePermeability(unsigned phaseIdx) const
    { return fractureRelativePermeabilities_[phaseIdx]; }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar fractureMobility(unsigned phaseIdx) const
    {
        return fractureRelativePermeabilities_[phaseIdx]
               / fractureFluidState_.viscosity(phaseIdx);
    }

    /*!
     * \brief Returns the average porosity within the fracture.
     */
    Scalar fracturePorosity() const
    { return fracturePorosity_; }

    /*!
     * \brief Returns the average intrinsic permeability within the
     *        fracture.
     */
    const DimMatrix& fractureIntrinsicPermeability() const
    { return fractureIntrinsicPermeability_; }

    /*!
     * \brief Return the volume [m^2] occupied by fractures within the
     *        given sub-control volume.
     */
    Scalar fractureVolume() const
    { return fractureVolume_; }

    /*!
     * \brief Returns a fluid state object which represents the
     *        thermodynamic state of the fluids within the fracture.
     */
    const FluidState& fractureFluidState() const
    { return fractureFluidState_; }

protected:
    FluidState fractureFluidState_;
    Scalar fractureVolume_;
    Scalar fracturePorosity_;
    DimMatrix fractureIntrinsicPermeability_;
    Scalar fractureRelativePermeabilities_[numPhases];
};

} // namespace Opm

#endif
