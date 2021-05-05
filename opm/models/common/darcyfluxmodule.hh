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
 * \brief This file contains the necessary classes to calculate the
 *        volumetric fluxes out of a pressure potential gradient using the
 *        Darcy relation.
 */
#ifndef EWOMS_DARCY_FLUX_MODULE_HH
#define EWOMS_DARCY_FLUX_MODULE_HH

#include "multiphasebaseproperties.hh"
#include <opm/models/common/quantitycallbacks.hh>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <cmath>

namespace Opm {

template <class TypeTag>
class DarcyIntensiveQuantities;

template <class TypeTag>
class DarcyExtensiveQuantities;

template <class TypeTag>
class DarcyBaseProblem;

/*!
 * \ingroup FluxModules
 * \brief Specifies a flux module which uses the Darcy relation.
 */
template <class TypeTag>
struct DarcyFluxModule
{
    using FluxIntensiveQuantities = DarcyIntensiveQuantities<TypeTag>;
    using FluxExtensiveQuantities = DarcyExtensiveQuantities<TypeTag>;
    using FluxBaseProblem = DarcyBaseProblem<TypeTag>;

    /*!
     * \brief Register all run-time parameters for the flux module.
     */
    static void registerParameters()
    { }
};

/*!
 * \ingroup FluxModules
 * \brief Provides the defaults for the parameters required by the
 *        Darcy velocity approach.
 */
template <class TypeTag>
class DarcyBaseProblem
{ };

/*!
 * \ingroup FluxModules
 * \brief Provides the intensive quantities for the Darcy flux module
 */
template <class TypeTag>
class DarcyIntensiveQuantities
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
protected:
    void update_(const ElementContext& elemCtx OPM_UNUSED,
                 unsigned dofIdx OPM_UNUSED,
                 unsigned timeIdx OPM_UNUSED)
    { }
};

/*!
 * \ingroup FluxModules
 * \brief Provides the Darcy flux module
 *
 * The commonly used Darcy relation looses its validity for Reynolds numbers \f$ Re <
 * 1\f$.  If one encounters flow velocities in porous media above this threshold, the
 * Forchheimer approach can be used.
 *
 * The Darcy equation is given by the following relation:
 *
 * \f[
  \vec{v}_\alpha =
  \left( \nabla p_\alpha - \rho_\alpha \vec{g}\right)
  \frac{\mu_\alpha}{k_{r,\alpha} K}
 \f]
 */
template <class TypeTag>
class DarcyExtensiveQuantities
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };

    using Toolbox = MathToolbox<Evaluation>;
    using ParameterCache = typename FluidSystem::template ParameterCache<Evaluation>;
    using EvalDimVector = Dune::FieldVector<Evaluation, dimWorld>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \brief Returns the intrinsic permeability tensor for a given
     *        sub-control volume face.
     */
    const DimMatrix& intrinsicPermability() const
    { return K_; }

    /*!
     * \brief Return the pressure potential gradient of a fluid phase
     *        at the face's integration point [Pa/m]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const EvalDimVector& potentialGrad(unsigned phaseIdx) const
    { return potentialGrad_[phaseIdx]; }

    /*!
     * \brief Return the filter velocity of a fluid phase at the
     *        face's integration point [m/s]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const EvalDimVector& filterVelocity(unsigned phaseIdx) const
    { return filterVelocity_[phaseIdx]; }

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
    short upstreamIndex_(unsigned phaseIdx) const
    { return upstreamDofIdx_[phaseIdx]; }

    short downstreamIndex_(unsigned phaseIdx) const
    { return downstreamDofIdx_[phaseIdx]; }

    /*!
     * \brief Calculate the gradients which are required to determine the volumetric fluxes
     *
     * The the upwind directions is also determined by method.
     */
    void calculateGradients_(const ElementContext& elemCtx,
                             unsigned faceIdx,
                             unsigned timeIdx)
    {
        const auto& gradCalc = elemCtx.gradientCalculator();
        PressureCallback<TypeTag> pressureCallback(elemCtx);

        const auto& scvf = elemCtx.stencil(timeIdx).interiorFace(faceIdx);
        const auto& faceNormal = scvf.normal();

        unsigned i = scvf.interiorIndex();
        unsigned j = scvf.exteriorIndex();
        interiorDofIdx_ = static_cast<short>(i);
        exteriorDofIdx_ = static_cast<short>(j);
        unsigned focusDofIdx = elemCtx.focusDofIndex();

        // calculate the "raw" pressure gradient
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                Valgrind::SetUndefined(potentialGrad_[phaseIdx]);
                continue;
            }

            pressureCallback.setPhaseIndex(phaseIdx);
            gradCalc.calculateGradient(potentialGrad_[phaseIdx],
                                       elemCtx,
                                       faceIdx,
                                       pressureCallback);
            Valgrind::CheckDefined(potentialGrad_[phaseIdx]);
        }

        // correct the pressure gradients by the gravitational acceleration
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity)) {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            const auto& gIn = elemCtx.problem().gravity(elemCtx, i, timeIdx);
            const auto& gEx = elemCtx.problem().gravity(elemCtx, j, timeIdx);

            const auto& intQuantsIn = elemCtx.intensiveQuantities(i, timeIdx);
            const auto& intQuantsEx = elemCtx.intensiveQuantities(j, timeIdx);

            const auto& posIn = elemCtx.pos(i, timeIdx);
            const auto& posEx = elemCtx.pos(j, timeIdx);
            const auto& posFace = scvf.integrationPos();

            // the distance between the centers of the control volumes
            DimVector distVecIn(posIn);
            DimVector distVecEx(posEx);
            DimVector distVecTotal(posEx);

            distVecIn -= posFace;
            distVecEx -= posFace;
            distVecTotal -= posIn;
            Scalar absDistTotalSquared = distVecTotal.two_norm2();
            for (unsigned phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
                if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                    continue;

                // calculate the hydrostatic pressure at the integration point of the face
                Evaluation pStatIn;

                if (std::is_same<Scalar, Evaluation>::value ||
                    interiorDofIdx_ == static_cast<int>(focusDofIdx))
                {
                    const Evaluation& rhoIn = intQuantsIn.fluidState().density(phaseIdx);
                    pStatIn = - rhoIn*(gIn*distVecIn);
                }
                else {
                    Scalar rhoIn = Toolbox::value(intQuantsIn.fluidState().density(phaseIdx));
                    pStatIn = - rhoIn*(gIn*distVecIn);
                }

                // the quantities on the exterior side of the face do not influence the
                // result for the TPFA scheme, so they can be treated as scalar values.
                Evaluation pStatEx;

                if (std::is_same<Scalar, Evaluation>::value ||
                    exteriorDofIdx_ == static_cast<int>(focusDofIdx))
                {
                    const Evaluation& rhoEx = intQuantsEx.fluidState().density(phaseIdx);
                    pStatEx = - rhoEx*(gEx*distVecEx);
                }
                else {
                    Scalar rhoEx = Toolbox::value(intQuantsEx.fluidState().density(phaseIdx));
                    pStatEx = - rhoEx*(gEx*distVecEx);
                }

                // compute the hydrostatic gradient between the two control volumes (this
                // gradient exhibitis the same direction as the vector between the two
                // control volume centers and the length (pStaticExterior -
                // pStaticInterior)/distanceInteriorToExterior
                Dune::FieldVector<Evaluation, dimWorld> f(distVecTotal);
                f *= (pStatEx - pStatIn)/absDistTotalSquared;

                // calculate the final potential gradient
                for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                    potentialGrad_[phaseIdx][dimIdx] += f[dimIdx];

                for (unsigned dimIdx = 0; dimIdx < potentialGrad_[phaseIdx].size(); ++dimIdx) {
                    if (!isfinite(potentialGrad_[phaseIdx][dimIdx])) {
                        throw NumericalIssue("Non-finite potential gradient for phase '"
                                             +std::string(FluidSystem::phaseName(phaseIdx))+"'");
                    }
                }
            }
        }

        Valgrind::SetUndefined(K_);
        elemCtx.problem().intersectionIntrinsicPermeability(K_, elemCtx, faceIdx, timeIdx);
        Valgrind::CheckDefined(K_);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                Valgrind::SetUndefined(potentialGrad_[phaseIdx]);
                continue;
            }

            // determine the upstream and downstream DOFs
            Evaluation tmp = 0.0;
            for (unsigned dimIdx = 0; dimIdx < faceNormal.size(); ++dimIdx)
                tmp += potentialGrad_[phaseIdx][dimIdx]*faceNormal[dimIdx];

            if (tmp > 0) {
                upstreamDofIdx_[phaseIdx] = exteriorDofIdx_;
                downstreamDofIdx_[phaseIdx] = interiorDofIdx_;
            }
            else {
                upstreamDofIdx_[phaseIdx] = interiorDofIdx_;
                downstreamDofIdx_[phaseIdx] = exteriorDofIdx_;
            }

            // we only carry the derivatives along if the upstream DOF is the one which
            // we currently focus on
            const auto& up = elemCtx.intensiveQuantities(upstreamDofIdx_[phaseIdx], timeIdx);
            if (upstreamDofIdx_[phaseIdx] == static_cast<int>(focusDofIdx))
                mobility_[phaseIdx] = up.mobility(phaseIdx);
            else
                mobility_[phaseIdx] = Toolbox::value(up.mobility(phaseIdx));
        }
    }

    /*!
     * \brief Calculate the gradients at the grid boundary which are required to
     *        determine the volumetric fluxes
     *
     * The the upwind directions is also determined by method.
     */
    template <class FluidState>
    void calculateBoundaryGradients_(const ElementContext& elemCtx,
                                     unsigned boundaryFaceIdx,
                                     unsigned timeIdx,
                                     const FluidState& fluidState)
    {
        const auto& gradCalc = elemCtx.gradientCalculator();
        BoundaryPressureCallback<TypeTag, FluidState> pressureCallback(elemCtx, fluidState);

        // calculate the pressure gradient
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                Valgrind::SetUndefined(potentialGrad_[phaseIdx]);
                continue;
            }

            pressureCallback.setPhaseIndex(phaseIdx);
            gradCalc.calculateBoundaryGradient(potentialGrad_[phaseIdx],
                                               elemCtx,
                                               boundaryFaceIdx,
                                               pressureCallback);
            Valgrind::CheckDefined(potentialGrad_[phaseIdx]);
        }

        const auto& scvf = elemCtx.stencil(timeIdx).boundaryFace(boundaryFaceIdx);
        auto i = scvf.interiorIndex();
        interiorDofIdx_ = static_cast<short>(i);
        exteriorDofIdx_ = -1;
        int focusDofIdx = elemCtx.focusDofIndex();

        // calculate the intrinsic permeability
        const auto& intQuantsIn = elemCtx.intensiveQuantities(i, timeIdx);
        K_ = intQuantsIn.intrinsicPermeability();

        // correct the pressure gradients by the gravitational acceleration
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity)) {
            // estimate the gravitational acceleration at a given SCV face
            // using the arithmetic mean
            const auto& gIn = elemCtx.problem().gravity(elemCtx, i, timeIdx);
            const auto& posIn = elemCtx.pos(i, timeIdx);
            const auto& posFace = scvf.integrationPos();

            // the distance between the face center and the center of the control volume
            DimVector distVecIn(posIn);
            distVecIn -= posFace;
            Scalar absDistSquared = distVecIn.two_norm2();
            Scalar gTimesDist = gIn*distVecIn;

            for (unsigned phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
                if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                    continue;

                // calculate the hydrostatic pressure at the integration point of the face
                Evaluation rhoIn = intQuantsIn.fluidState().density(phaseIdx);
                Evaluation pStatIn = - gTimesDist*rhoIn;

                Valgrind::CheckDefined(pStatIn);
                // compute the hydrostatic gradient between the control volume and face integration
                // point. This gradient exhibitis the same direction as the vector between the
                // control volume center and face integration point (-distVecIn / absDist) and the
                // length of the vector is -pStaticIn / absDist. Note that the two negatives become
                // + and the 1 / (absDist * absDist) -> absDistSquared.
                EvalDimVector f(distVecIn);
                f *= pStatIn / absDistSquared;

                // calculate the final potential gradient
                for (unsigned dimIdx = 0; dimIdx < dimWorld; ++dimIdx)
                    potentialGrad_[phaseIdx][dimIdx] += f[dimIdx];

                Valgrind::CheckDefined(potentialGrad_[phaseIdx]);
                for (unsigned dimIdx = 0; dimIdx < potentialGrad_[phaseIdx].size(); ++dimIdx) {
                    if (!isfinite(potentialGrad_[phaseIdx][dimIdx])) {
                        throw NumericalIssue("Non finite potential gradient for phase '"
                                             +std::string(FluidSystem::phaseName(phaseIdx))+"'");
                    }
                }
            }
        }

        // determine the upstream and downstream DOFs
        const auto& faceNormal = scvf.normal();

        const auto& matParams = elemCtx.problem().materialLawParams(elemCtx, i, timeIdx);

        Scalar kr[numPhases];
        MaterialLaw::relativePermeabilities(kr, matParams, fluidState);
        Valgrind::CheckDefined(kr);

        for (unsigned phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            Evaluation tmp = 0.0;
            for (unsigned dimIdx = 0; dimIdx < faceNormal.size(); ++dimIdx)
                tmp += potentialGrad_[phaseIdx][dimIdx]*faceNormal[dimIdx];

            if (tmp > 0) {
                upstreamDofIdx_[phaseIdx] = exteriorDofIdx_;
                downstreamDofIdx_[phaseIdx] = interiorDofIdx_;
            }
            else {
                upstreamDofIdx_[phaseIdx] = interiorDofIdx_;
                downstreamDofIdx_[phaseIdx] = exteriorDofIdx_;
            }

            // take the phase mobility from the DOF in upstream direction
            if (upstreamDofIdx_[phaseIdx] < 0) {
                if (interiorDofIdx_ == focusDofIdx)
                    mobility_[phaseIdx] =
                        kr[phaseIdx] / fluidState.viscosity(phaseIdx);
                else
                    mobility_[phaseIdx] =
                        Toolbox::value(kr[phaseIdx])
                        / Toolbox::value(fluidState.viscosity(phaseIdx));
            }
            else if (upstreamDofIdx_[phaseIdx] != focusDofIdx)
                mobility_[phaseIdx] = Toolbox::value(intQuantsIn.mobility(phaseIdx));
            else
                mobility_[phaseIdx] = intQuantsIn.mobility(phaseIdx);
            Valgrind::CheckDefined(mobility_[phaseIdx]);
        }
    }

    /*!
     * \brief Calculate the volumetric fluxes of all phases
     *
     * The pressure potentials and upwind directions must already be
     * determined before calling this method!
     */
    void calculateFluxes_(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    {
        const auto& scvf = elemCtx.stencil(timeIdx).interiorFace(scvfIdx);
        const DimVector& normal = scvf.normal();
        Valgrind::CheckDefined(normal);

        for (unsigned phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            filterVelocity_[phaseIdx] = 0.0;
            volumeFlux_[phaseIdx] = 0.0;
            if (!elemCtx.model().phaseIsConsidered(phaseIdx))
                continue;

            asImp_().calculateFilterVelocity_(phaseIdx);
            Valgrind::CheckDefined(filterVelocity_[phaseIdx]);

            volumeFlux_[phaseIdx] = 0.0;
            for (unsigned i = 0; i < normal.size(); ++i)
                volumeFlux_[phaseIdx] += filterVelocity_[phaseIdx][i] * normal[i];
        }
    }

    /*!
     * \brief Calculate the volumetric fluxes at a boundary face of all fluid phases
     *
     * The pressure potentials and upwind directions must already be determined before
     * calling this method!
     */
    void calculateBoundaryFluxes_(const ElementContext& elemCtx,
                                  unsigned boundaryFaceIdx,
                                  unsigned timeIdx)
    {
        const auto& scvf = elemCtx.stencil(timeIdx).boundaryFace(boundaryFaceIdx);
        const DimVector& normal = scvf.normal();
        Valgrind::CheckDefined(normal);

        for (unsigned phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!elemCtx.model().phaseIsConsidered(phaseIdx)) {
                filterVelocity_[phaseIdx] = 0.0;
                volumeFlux_[phaseIdx] = 0.0;
                continue;
            }

            asImp_().calculateFilterVelocity_(phaseIdx);
            Valgrind::CheckDefined(filterVelocity_[phaseIdx]);
            volumeFlux_[phaseIdx] = 0.0;
            for (unsigned i = 0; i < normal.size(); ++i)
                volumeFlux_[phaseIdx] += filterVelocity_[phaseIdx][i] * normal[i];
        }
    }

    void calculateFilterVelocity_(unsigned phaseIdx)
    {
#ifndef NDEBUG
        assert(isfinite(mobility_[phaseIdx]));
        for (unsigned i = 0; i < K_.M(); ++ i)
            for (unsigned j = 0; j < K_.N(); ++ j)
                assert(std::isfinite(K_[i][j]));
#endif

        K_.mv(potentialGrad_[phaseIdx], filterVelocity_[phaseIdx]);
        filterVelocity_[phaseIdx] *= - mobility_[phaseIdx];

#ifndef NDEBUG
        for (unsigned i = 0; i < filterVelocity_[phaseIdx].size(); ++ i)
            assert(isfinite(filterVelocity_[phaseIdx][i]));
#endif
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

protected:
    // intrinsic permeability tensor and its square root
    DimMatrix K_;

    // mobilities of all fluid phases [1 / (Pa s)]
    Evaluation mobility_[numPhases];

    // filter velocities of all phases [m/s]
    EvalDimVector filterVelocity_[numPhases];

    // the volumetric flux of all fluid phases over the control
    // volume's face [m^3/s / m^2]
    Evaluation volumeFlux_[numPhases];

    // pressure potential gradients of all phases [Pa / m]
    EvalDimVector potentialGrad_[numPhases];

    // upstream, downstream, interior and exterior DOFs
    short upstreamDofIdx_[numPhases];
    short downstreamDofIdx_[numPhases];
    short interiorDofIdx_;
    short exteriorDofIdx_;
};

} // namespace Opm

#endif
