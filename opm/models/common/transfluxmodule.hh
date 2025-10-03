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
 * \brief This file contains the flux module that uses transmissibilities
 *
 * The transmissibility approach to fluxes used here is limited
 * to the two-point flux approximation
 */
#ifndef EWOMS_TRANS_FLUX_MODULE_HH
#define EWOMS_TRANS_FLUX_MODULE_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/common/multiphasebaseproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/discretization/ecfv/ecfvdiscretization.hh>
#include <opm/models/utils/signum.hh>

#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace Opm {

template <class TypeTag>
class TransIntensiveQuantities;

template <class TypeTag>
class TransExtensiveQuantities;

template <class TypeTag>
class TransBaseProblem;

/*!
 * \brief Specifies a flux module which uses transmissibilities.
 */
template <class TypeTag>
struct TransFluxModule
{
    using FluxIntensiveQuantities = TransIntensiveQuantities<TypeTag>;
    using FluxExtensiveQuantities = TransExtensiveQuantities<TypeTag>;
    using FluxBaseProblem = TransBaseProblem<TypeTag>;

    /*!
     * \brief Register all run-time parameters for the flux module.
     */
    static void registerParameters()
    {}
};

/*!
 * \brief Provides the defaults for the parameters required by the
 *        transmissibility based volume flux calculation.
 */
template <class TypeTag>
class TransBaseProblem
{};

/*!
 * \brief Provides the intensive quantities for the transmissibility based flux module
 */
template <class TypeTag>
class TransIntensiveQuantities
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

protected:
    void update_(const ElementContext&, unsigned, unsigned)
    {}
};

/*!
 * \brief Provides the transmissibility based flux module
 */
template <class TypeTag>
class TransExtensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Discretization = GetPropType<TypeTag, Properties::Discretization>;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = FluidSystem::numPhases };

    typedef MathToolbox<Evaluation> Toolbox;
    typedef Dune::FieldVector<Evaluation, dimWorld> EvalDimVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

public:
    /*!
     * \brief Return the intrinsic permeability tensor at a face [m^2]
     */
    const DimMatrix& intrinsicPermeability() const
    {
        throw std::logic_error("The ECL transmissibility module does not "
                               "provide an explicit intrinsic permeability");
    }

    /*!
     * \brief Return the pressure potential gradient of a fluid phase at the
     *        face's integration point [Pa/m]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const EvalDimVector& potentialGrad(unsigned) const
    {
        throw std::logic_error("The ECL transmissibility module does not "
                               "provide explicit potential gradients");
    }

    /*!
     * \brief Return the gravity corrected pressure difference between the interior and
     *        the exterior of a face.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Evaluation& pressureDifference(unsigned phaseIdx) const
    { return pressureDifference_[phaseIdx]; }

    /*!
     * \brief Return the filter velocity of a fluid phase at the face's integration point
     *        [m/s]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const EvalDimVector& filterVelocity(unsigned) const
    {
        throw std::logic_error("The ECL transmissibility module does not "
                               "provide explicit filter velocities");
    }

    /*!
     * \brief Return the volume flux of a fluid phase at the face's integration point
     *        \f$[m^3/s / m^2]\f$
     *
     * This is the fluid volume of a phase per second and per square meter of face
     * area.
     *
     * \param phaseIdx The index of the fluid phase
     */
    const Evaluation& volumeFlux(unsigned phaseIdx) const
    { return volumeFlux_[phaseIdx]; }

protected:
    /*!
     * \brief Returns the local index of the degree of freedom in which is
     *        in upstream direction.
     *
     * i.e., the DOF which exhibits a higher effective pressure for
     * the given phase.
     */
    unsigned upstreamIndex_(unsigned phaseIdx) const
    {
        assert(phaseIdx < numPhases);

        return upIdx_[phaseIdx];
    }

    /*!
     * \brief Returns the local index of the degree of freedom in which is
     *        in downstream direction.
     *
     * i.e., the DOF which exhibits a lower effective pressure for the
     * given phase.
     */
    unsigned downstreamIndex_(unsigned phaseIdx) const
    {
        assert(phaseIdx < numPhases);

        return dnIdx_[phaseIdx];
    }

    void updateSolvent(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    { asImp_().updateVolumeFluxTrans(elemCtx, scvfIdx, timeIdx); }

    void updatePolymer(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    { asImp_().updateShearMultipliers(elemCtx, scvfIdx, timeIdx); }

    /*!
     * \brief Update the required gradients for interior faces
     */
    void calculateGradients_(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    {
        Valgrind::SetUndefined(*this);

        // only valid for element center finite volume discretization
        static_assert(std::is_same_v<Discretization, EcfvDiscretization<TypeTag>>);

        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        interiorDofIdx_ = scvf.interiorIndex();
        exteriorDofIdx_ = scvf.exteriorIndex();
        assert(interiorDofIdx_ != exteriorDofIdx_);

        const unsigned I = stencil.globalSpaceIndex(interiorDofIdx_);
        const unsigned J = stencil.globalSpaceIndex(exteriorDofIdx_);

        const Scalar trans = transmissibility_(elemCtx, scvfIdx, timeIdx);

        // estimate the gravity correction: for performance reasons we use a simplified
        // approach for this flux module that assumes that gravity is constant and always
        // acts into the downwards direction. (i.e., no centrifuge experiments, sorry.)
        const Scalar g = elemCtx.problem().gravity()[dimWorld - 1];

        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx_, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx_, timeIdx);

        const Scalar zIn = dofCenterDepth_(elemCtx, interiorDofIdx_, timeIdx);
        const Scalar zEx = dofCenterDepth_(elemCtx, exteriorDofIdx_, timeIdx);
        // the distances from the DOF's depths. (i.e., the additional depth of the
        // exterior DOF)
        const Scalar distZ = zIn - zEx;

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            // check shortcut: if the mobility of the phase is zero in the interior as
            // well as the exterior DOF, we can skip looking at the phase.
            if (intQuantsIn.mobility(phaseIdx) <= 0.0 &&
                intQuantsEx.mobility(phaseIdx) <= 0.0)
            {
                upIdx_[phaseIdx] = interiorDofIdx_;
                dnIdx_[phaseIdx] = exteriorDofIdx_;
                pressureDifference_[phaseIdx] = 0.0;
                volumeFlux_[phaseIdx] = 0.0;
                continue;
            }

            // do the gravity correction: compute the hydrostatic pressure for the
            // external at the depth of the internal one
            const Evaluation& rhoIn = intQuantsIn.fluidState().density(phaseIdx);
            const Scalar rhoEx = Toolbox::value(intQuantsEx.fluidState().density(phaseIdx));
            const Evaluation rhoAvg = (rhoIn + rhoEx) / 2;

            const Evaluation& pressureInterior = intQuantsIn.fluidState().pressure(phaseIdx);
            Evaluation pressureExterior = Toolbox::value(intQuantsEx.fluidState().pressure(phaseIdx));

            pressureExterior += rhoAvg * (distZ * g);

            pressureDifference_[phaseIdx] = pressureExterior - pressureInterior;

            // decide the upstream index for the phase. for this we make sure that the
            // degree of freedom which is regarded upstream if both pressures are equal
            // is always the same: if the pressure is equal, the DOF with the lower
            // global index is regarded to be the upstream one.
            if (pressureDifference_[phaseIdx] > 0.0) {
                upIdx_[phaseIdx] = exteriorDofIdx_;
                dnIdx_[phaseIdx] = interiorDofIdx_;
            }
            else if (pressureDifference_[phaseIdx] < 0.0) {
                upIdx_[phaseIdx] = interiorDofIdx_;
                dnIdx_[phaseIdx] = exteriorDofIdx_;
            }
            else {
                // if the pressure difference is zero, we chose the DOF which has the
                // larger volume associated to it as upstream DOF
                const Scalar Vin = elemCtx.dofVolume(interiorDofIdx_, /*timeIdx=*/0);
                const Scalar Vex = elemCtx.dofVolume(exteriorDofIdx_, /*timeIdx=*/0);
                if (Vin > Vex) {
                    upIdx_[phaseIdx] = interiorDofIdx_;
                    dnIdx_[phaseIdx] = exteriorDofIdx_;
                }
                else if (Vin < Vex) {
                    upIdx_[phaseIdx] = exteriorDofIdx_;
                    dnIdx_[phaseIdx] = interiorDofIdx_;
                }
                else {
                    assert(Vin == Vex);
                    // if the volumes are also equal, we pick the DOF which exhibits the
                    // smaller global index
                    if (I < J) {
                        upIdx_[phaseIdx] = interiorDofIdx_;
                        dnIdx_[phaseIdx] = exteriorDofIdx_;
                    }
                    else {
                        upIdx_[phaseIdx] = exteriorDofIdx_;
                        dnIdx_[phaseIdx] = interiorDofIdx_;
                    }
                }
            }

            // this is slightly hacky because in the automatic differentiation case, it
            // only works for the element centered finite volume method. for ebos this
            // does not matter, though.
            const unsigned upstreamIdx = upstreamIndex_(phaseIdx);
            const auto& up = elemCtx.intensiveQuantities(upstreamIdx, timeIdx);

            if (upstreamIdx == interiorDofIdx_) {
                volumeFlux_[phaseIdx] =
                    pressureDifference_[phaseIdx] * up.mobility(phaseIdx) * (-trans);
            }
            else {
                volumeFlux_[phaseIdx] =
                    pressureDifference_[phaseIdx] * (Toolbox::value(up.mobility(phaseIdx)) * (-trans));
            }
        }
    }

    /*!
     * \brief Update the required gradients for boundary faces
     */
    template <class FluidState>
    void calculateBoundaryGradients_(const ElementContext& elemCtx,
                                     unsigned scvfIdx,
                                     unsigned timeIdx,
                                     const FluidState& exFluidState)
    {
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.boundaryFace(scvfIdx);

        interiorDofIdx_ = scvf.interiorIndex();

        const Scalar trans = transmissibilityBoundary_(elemCtx, scvfIdx, timeIdx);

        // estimate the gravity correction: for performance reasons we use a simplified
        // approach for this flux module that assumes that gravity is constant and always
        // acts into the downwards direction. (i.e., no centrifuge experiments, sorry.)
        const Scalar g = elemCtx.problem().gravity()[dimWorld - 1];

        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx_, timeIdx);

        // this is quite hacky because the dune grid interface does not provide a
        // cellCenterDepth() method (so we ask the problem to provide it). The "good"
        // solution would be to take the Z coordinate of the element centroids, but since
        // ECL seems to like to be inconsistent on that front, it needs to be done like
        // here...
        const Scalar zIn = dofCenterDepth_(elemCtx, interiorDofIdx_, timeIdx);
        const Scalar zEx = scvf.integrationPos()[dimWorld - 1];

        // the distances from the DOF's depths. (i.e., the additional depth of the
        // exterior DOF)
        const Scalar distZ = zIn - zEx;

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            // do the gravity correction: compute the hydrostatic pressure for the
            // integration position
            const Evaluation& rhoIn = intQuantsIn.fluidState().density(phaseIdx);
            const auto& rhoEx = exFluidState.density(phaseIdx);
            const Evaluation rhoAvg = (rhoIn + rhoEx) / 2;

            const Evaluation& pressureInterior = intQuantsIn.fluidState().pressure(phaseIdx);
            Evaluation pressureExterior = exFluidState.pressure(phaseIdx);
            pressureExterior += rhoAvg * (distZ * g);

            pressureDifference_[phaseIdx] = pressureExterior - pressureInterior;

            // decide the upstream index for the phase. for this we make sure that the
            // degree of freedom which is regarded upstream if both pressures are equal
            // is always the same: if the pressure is equal, the DOF with the lower
            // global index is regarded to be the upstream one.
            if (pressureDifference_[phaseIdx] > 0.0) {
                upIdx_[phaseIdx] = -1;
                dnIdx_[phaseIdx] = interiorDofIdx_;
            }
            else {
                upIdx_[phaseIdx] = interiorDofIdx_;
                dnIdx_[phaseIdx] = -1;
            }

            const short upstreamIdx = upstreamIndex_(phaseIdx);
            if (upstreamIdx == interiorDofIdx_) {
                // this is slightly hacky because in the automatic differentiation case, it
                // only works for the element centered finite volume method. for ebos this
                // does not matter, though.
                const auto& up = elemCtx.intensiveQuantities(upstreamIdx, timeIdx);

                volumeFlux_[phaseIdx] =
                    pressureDifference_[phaseIdx] * up.mobility(phaseIdx) * (-trans);
            }
            else {
                // compute the phase mobility using the material law parameters of the
                // interior element. \todo {this could probably be done more efficiently}
                const auto& matParams =
                    elemCtx.problem().materialLawParams(elemCtx,
                                                        interiorDofIdx_,
                                                        /*timeIdx=*/0);
                std::array<typename FluidState::Scalar,numPhases> kr;
                MaterialLaw::relativePermeabilities(kr, matParams, exFluidState);

                const auto& mob = kr[phaseIdx] / exFluidState.viscosity(phaseIdx);
                volumeFlux_[phaseIdx] =
                    pressureDifference_[phaseIdx] * mob * (-trans);
            }
        }
    }

    /*!
     * \brief Update the volumetric fluxes for all fluid phases on the interior faces of the context
     */
    void calculateFluxes_(const ElementContext&, unsigned, unsigned)
    {}

    void calculateBoundaryFluxes_(const ElementContext&, unsigned, unsigned)
    {}

private:
    Scalar transmissibility_(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx) const
    {
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& face = stencil.interiorFace(scvfIdx);
        const auto& interiorPos = stencil.subControlVolume(face.interiorIndex()).globalPos();
        const auto& exteriorPos = stencil.subControlVolume(face.exteriorIndex()).globalPos();
        const auto distVec0 = face.integrationPos() - interiorPos;
        const auto distVec1 = face.integrationPos() - exteriorPos;
        const Scalar ndotDistIn = std::abs(face.normal() * distVec0);
        const Scalar ndotDistExt = std::abs(face.normal() * distVec1);

        const Scalar distSquaredIn = distVec0 * distVec0;
        const Scalar distSquaredExt = distVec1 * distVec1;
        const auto& K0mat = elemCtx.problem().intrinsicPermeability(elemCtx, face.interiorIndex(), timeIdx);
        const auto& K1mat = elemCtx.problem().intrinsicPermeability(elemCtx, face.exteriorIndex(), timeIdx);

        // the permeability per definition aligns with the grid
        // we only support diagonal permeability tensor
        // and can therefore neglect off-diagonal values
        int idx = 0;
        Scalar val = 0.0;
        for (unsigned i = 0; i < dimWorld; ++i) {
            if (std::abs(face.normal()[i]) > val) {
                val = std::abs(face.normal()[i]);
                idx = i;
            }
        }
        const Scalar& K0 = K0mat[idx][idx];
        const Scalar& K1 = K1mat[idx][idx];
        const Scalar T0 = K0 * ndotDistIn / distSquaredIn;
        const Scalar T1 = K1 * ndotDistExt / distSquaredExt;
        return T0 * T1 / (T0 + T1);
    }
    Scalar transmissibilityBoundary_(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx) const
    {
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& face = stencil.interiorFace(scvfIdx);
        const auto& interiorPos = stencil.subControlVolume(face.interiorIndex()).globalPos();
        const auto distVec0 = face.integrationPos() - interiorPos;
        const Scalar ndotDistIn = face.normal() * distVec0;
        const Scalar distSquaredIn = distVec0 * distVec0;
        const auto& K0mat = elemCtx.problem().intrinsicPermeability(elemCtx, face.interiorIndex(), timeIdx);

        // the permeability per definition aligns with the grid
        // we only support diagonal permeability tensor
        // and can therefore neglect off-diagonal values
        int idx = 0;
        Scalar val = 0.0;
        for (unsigned i = 0; i < dimWorld; ++i) {
            if (std::abs(face.normal()[i]) > val) {
                val = std::abs(face.normal()[i]);
                idx = i;
            }
        }
        const Scalar& K0 = K0mat[idx][idx];
        const Scalar T0 = K0 * ndotDistIn / distSquaredIn;
        return T0;
    }

    template <class Context>
    Scalar dofCenterDepth_(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const auto& pos = context.pos(spaceIdx, timeIdx);
        return pos[dimWorld-1];
    }

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // the volumetric flux of all phases [m^3/s]
    std::array<Evaluation, numPhases> volumeFlux_;

    // the difference in effective pressure between the exterior and the interior degree
    // of freedom [Pa]
    std::array<Evaluation, numPhases> pressureDifference_;

    // the local indices of the interior and exterior degrees of freedom
    unsigned short interiorDofIdx_{};
    unsigned short exteriorDofIdx_{};
    std::array<short, numPhases> upIdx_{};
    std::array<short, numPhases> dnIdx_{};
};

} // namespace Opm

#endif
