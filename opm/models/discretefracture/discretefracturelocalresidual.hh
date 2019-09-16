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
 * \copydoc Opm::DiscreteFractureLocalResidual
 */
#ifndef EWOMS_DISCRETE_FRACTURE_LOCAL_RESIDUAL_BASE_HH
#define EWOMS_DISCRETE_FRACTURE_LOCAL_RESIDUAL_BASE_HH

#include <ewoms/models/immiscible/immisciblelocalresidual.hh>

namespace Opm {

/*!
 * \ingroup DiscreteFractureModel
 *
 * \brief Calculates the local residual of the discrete fracture
 *        immiscible multi-phase model.
 */
template <class TypeTag>
class DiscreteFractureLocalResidual : public ImmiscibleLocalResidual<TypeTag>
{
    typedef ImmiscibleLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };

    typedef Opm::EnergyModule<TypeTag, enableEnergy> EnergyModule;

public:
    /*!
     * \brief Adds the amount all conservation quantities (e.g. phase
     *        mass) within a single fluid phase
     *
     * \copydetails Doxygen::storageParam
     * \copydetails Doxygen::dofCtxParams
     * \copydetails Doxygen::phaseIdxParam
     */
    void addPhaseStorage(EqVector& storage,
                         const ElementContext& elemCtx,
                         unsigned dofIdx,
                         unsigned timeIdx,
                         unsigned phaseIdx) const
    {
        EqVector phaseStorage(0.0);
        ParentType::addPhaseStorage(phaseStorage, elemCtx, dofIdx, timeIdx, phaseIdx);

        const auto& problem = elemCtx.problem();
        const auto& fractureMapper = problem.fractureMapper();
        unsigned globalIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);

        if (!fractureMapper.isFractureVertex(globalIdx)) {
            // don't do anything in addition to the immiscible model for degrees of
            // freedom that do not feature fractures
            storage += phaseStorage;
            return;
        }

        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& scv = elemCtx.stencil(timeIdx).subControlVolume(dofIdx);

        // reduce the matrix storage by the fracture volume
        phaseStorage *= 1 - intQuants.fractureVolume()/scv.volume();

        // add the storage term inside the fractures
        const auto& fsFracture = intQuants.fractureFluidState();

        phaseStorage[conti0EqIdx + phaseIdx] +=
            intQuants.fracturePorosity()*
            fsFracture.saturation(phaseIdx) *
            fsFracture.density(phaseIdx) *
            intQuants.fractureVolume()/scv.volume();

        EnergyModule::addFracturePhaseStorage(phaseStorage, intQuants, scv,
                                              phaseIdx);

        // add the result to the overall storage term
        storage += phaseStorage;
    }

    /*!
     * \copydoc FvBaseLocalResidual::computeFlux
     */
    void computeFlux(RateVector& flux,
                     const ElementContext& elemCtx,
                     unsigned scvfIdx,
                     unsigned timeIdx) const
    {
        ParentType::computeFlux(flux, elemCtx, scvfIdx, timeIdx);

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        unsigned i = extQuants.interiorIndex();
        unsigned j = extQuants.exteriorIndex();
        unsigned I = elemCtx.globalSpaceIndex(i, timeIdx);
        unsigned J = elemCtx.globalSpaceIndex(j, timeIdx);
        const auto& fractureMapper = elemCtx.problem().fractureMapper();
        if (!fractureMapper.isFractureEdge(I, J))
            // do nothing if the edge from i to j is not part of a
            // fracture
            return;

        const auto& scvf = elemCtx.stencil(timeIdx).interiorFace(scvfIdx);
        Scalar scvfArea = scvf.area();

        // advective mass fluxes of all phases
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            // reduce the matrix mass flux by the width of the scv
            // face that is occupied by the fracture. As usual, the
            // fracture is shared between two SCVs, so the its width
            // needs to be divided by two.
            flux[conti0EqIdx + phaseIdx] *=
                1 - extQuants.fractureWidth() / (2 * scvfArea);

            // intensive quantities of the upstream and the downstream DOFs
            unsigned upIdx = static_cast<unsigned>(extQuants.upstreamIndex(phaseIdx));
            const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
            flux[conti0EqIdx + phaseIdx] +=
                extQuants.fractureVolumeFlux(phaseIdx) * up.fractureFluidState().density(phaseIdx);
        }

        EnergyModule::handleFractureFlux(flux, elemCtx, scvfIdx, timeIdx);
    }
};

} // namespace Opm

#endif
