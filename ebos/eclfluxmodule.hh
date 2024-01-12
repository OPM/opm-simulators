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
 * \brief This file contains the flux module which is used for ECL problems
 *
 * This approach to fluxes is very specific to two-point flux approximation and applies
 * what the Eclipse Technical Description calls the "NEWTRAN" transmissibility approach.
 */
#ifndef EWOMS_ECL_FLUX_MODULE_HH
#define EWOMS_ECL_FLUX_MODULE_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/utils/signum.hh>

#include <array>

namespace Opm {

template <class TypeTag>
class EclTransIntensiveQuantities;

template <class TypeTag>
class EclTransExtensiveQuantities;

template <class TypeTag>
class EclTransBaseProblem;

/*!
 * \ingroup EclBlackOilSimulator
 * \brief Specifies a flux module which uses ECL transmissibilities.
 */
template <class TypeTag>
struct EclTransFluxModule
{
    using FluxIntensiveQuantities = EclTransIntensiveQuantities<TypeTag>;
    using FluxExtensiveQuantities = EclTransExtensiveQuantities<TypeTag>;
    using FluxBaseProblem = EclTransBaseProblem<TypeTag>;

    /*!
     * \brief Register all run-time parameters for the flux module.
     */
    static void registerParameters()
    { }
};

/*!
 * \ingroup EclBlackOilSimulator
 * \brief Provides the defaults for the parameters required by the
 *        transmissibility based volume flux calculation.
 */
template <class TypeTag>
class EclTransBaseProblem
{ };

/*!
 * \ingroup EclBlackOilSimulator
 * \brief Provides the intensive quantities for the ECL flux module
 */
template <class TypeTag>
class EclTransIntensiveQuantities
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
protected:
    void update_(const ElementContext&, unsigned, unsigned)
    { }
};

/*!
 * \ingroup EclBlackOilSimulator
 * \brief Provides the ECL flux module
 */
template <class TypeTag>
class EclTransExtensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;

    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;

    enum { dimWorld = GridView::dimensionworld };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { numPhases = FluidSystem::numPhases };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };
    enum { enableExtbo = getPropValue<TypeTag, Properties::EnableExtbo>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };

    using Toolbox = MathToolbox<Evaluation>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using EvalDimVector = Dune::FieldVector<Evaluation, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \brief Return the intrinsic permeability tensor at a face [m^2]
     */
    const DimMatrix& intrinsicPermeability() const
    {
        throw std::invalid_argument("The ECL transmissibility module does not provide an explicit intrinsic permeability");
    }

    /*!
     * \brief Return the pressure potential gradient of a fluid phase at the
     *        face's integration point [Pa/m]
     *
     * \param phaseIdx The index of the fluid phase
     */
    const EvalDimVector& potentialGrad(unsigned) const
    {
        throw std::invalid_argument("The ECL transmissibility module does not provide explicit potential gradients");
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
        throw std::invalid_argument("The ECL transmissibility module does not provide explicit filter velocities");
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

public:

    static void volumeAndPhasePressureDifferences(std::array<short, numPhases>& upIdx,
                                                  std::array<short, numPhases>& dnIdx,
                                                  Evaluation (&volumeFlux)[numPhases],
                                                  Evaluation (&pressureDifferences)[numPhases],
                                                  const ElementContext& elemCtx,
                                                  unsigned scvfIdx,
                                                  unsigned timeIdx)
    {
        const auto& problem = elemCtx.problem();
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);
        unsigned interiorDofIdx = scvf.interiorIndex();
        unsigned exteriorDofIdx = scvf.exteriorIndex();

        assert(interiorDofIdx != exteriorDofIdx);

        unsigned I = stencil.globalSpaceIndex(interiorDofIdx);
        unsigned J = stencil.globalSpaceIndex(exteriorDofIdx);
        Scalar trans = problem.transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
        Scalar faceArea = scvf.area();
        Scalar thpres = problem.thresholdPressure(I, J);

        // estimate the gravity correction: for performance reasons we use a simplified
        // approach for this flux module that assumes that gravity is constant and always
        // acts into the downwards direction. (i.e., no centrifuge experiments, sorry.)
        Scalar g = elemCtx.problem().gravity()[dimWorld - 1];

        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx, timeIdx);

        // this is quite hacky because the dune grid interface does not provide a
        // cellCenterDepth() method (so we ask the problem to provide it). The "good"
        // solution would be to take the Z coordinate of the element centroids, but since
        // ECL seems to like to be inconsistent on that front, it needs to be done like
        // here...
        Scalar zIn = problem.dofCenterDepth(elemCtx, interiorDofIdx, timeIdx);
        Scalar zEx = problem.dofCenterDepth(elemCtx, exteriorDofIdx, timeIdx);

        // the distances from the DOF's depths. (i.e., the additional depth of the
        // exterior DOF)
        Scalar distZ = zIn - zEx;

        Scalar Vin = elemCtx.dofVolume(interiorDofIdx, /*timeIdx=*/0);
        Scalar Vex = elemCtx.dofVolume(exteriorDofIdx, /*timeIdx=*/0);

        for (unsigned phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;
            calculatePhasePressureDiff_(upIdx[phaseIdx],
                                        dnIdx[phaseIdx],
                                        pressureDifferences[phaseIdx],
                                        intQuantsIn,
                                        intQuantsEx,
                                        phaseIdx,//input
                                        interiorDofIdx,//input
                                        exteriorDofIdx,//input
                                        Vin,
                                        Vex,
                                        I,
                                        J,
                                        distZ*g,
                                        thpres);
            if (pressureDifferences[phaseIdx] == 0) {
                volumeFlux[phaseIdx] = 0.0;
                continue;
            }

            const bool upwindIsInterior = (static_cast<unsigned>(upIdx[phaseIdx]) == interiorDofIdx);
            const IntensiveQuantities& up = upwindIsInterior ? intQuantsIn : intQuantsEx;
            // Use arithmetic average (more accurate with harmonic, but that requires recomputing the transmissbility)
            const Evaluation transMult = (intQuantsIn.rockCompTransMultiplier() + Toolbox::value(intQuantsEx.rockCompTransMultiplier()))/2;

            const auto& materialLawManager = problem.materialLawManager();
            FaceDir::DirEnum facedir = FaceDir::DirEnum::Unknown;
            if (materialLawManager->hasDirectionalRelperms()) {
                facedir = scvf.faceDirFromDirId();  // direction (X, Y, or Z) of the face
            }
            if (upwindIsInterior)
                volumeFlux[phaseIdx] =
                    pressureDifferences[phaseIdx]*up.mobility(phaseIdx, facedir)*transMult*(-trans/faceArea);
            else
                volumeFlux[phaseIdx] =
                    pressureDifferences[phaseIdx]*
                        (Toolbox::value(up.mobility(phaseIdx, facedir))*transMult*(-trans/faceArea));
        }
    }

    template<class EvalType>
    static void calculatePhasePressureDiff_(short& upIdx,
                                            short& dnIdx,
                                            EvalType& pressureDifference,
                                            const IntensiveQuantities& intQuantsIn,
                                            const IntensiveQuantities& intQuantsEx,
                                            const unsigned phaseIdx,
                                            const unsigned interiorDofIdx,
                                            const unsigned exteriorDofIdx,
                                            const Scalar Vin,
                                            const Scalar Vex,
                                            const unsigned globalIndexIn,
                                            const unsigned globalIndexEx,
                                            const Scalar distZg,
                                            const Scalar thpres
        )
    {

        // check shortcut: if the mobility of the phase is zero in the interior as
        // well as the exterior DOF, we can skip looking at the phase.
        if (intQuantsIn.mobility(phaseIdx) <= 0.0 &&
            intQuantsEx.mobility(phaseIdx) <= 0.0)
        {
            upIdx = interiorDofIdx;
            dnIdx = exteriorDofIdx;
            pressureDifference = 0.0;
            return;
        }

        // do the gravity correction: compute the hydrostatic pressure for the
        // external at the depth of the internal one
        const Evaluation& rhoIn = intQuantsIn.fluidState().density(phaseIdx);
        Scalar rhoEx = Toolbox::value(intQuantsEx.fluidState().density(phaseIdx));
        Evaluation rhoAvg = (rhoIn + rhoEx)/2;

        const Evaluation& pressureInterior = intQuantsIn.fluidState().pressure(phaseIdx);
        Evaluation pressureExterior = Toolbox::value(intQuantsEx.fluidState().pressure(phaseIdx));
        if (enableExtbo) // added stability; particulary useful for solvent migrating in pure water
                         // where the solvent fraction displays a 0/1 behaviour ...
            pressureExterior += Toolbox::value(rhoAvg)*(distZg);
        else
            pressureExterior += rhoAvg*(distZg);

        pressureDifference = pressureExterior - pressureInterior;

        // decide the upstream index for the phase. for this we make sure that the
        // degree of freedom which is regarded upstream if both pressures are equal
        // is always the same: if the pressure is equal, the DOF with the lower
        // global index is regarded to be the upstream one.
        if (pressureDifference > 0.0) {
            upIdx = exteriorDofIdx;
            dnIdx = interiorDofIdx;
        }
        else if (pressureDifference < 0.0) {
            upIdx = interiorDofIdx;
            dnIdx = exteriorDofIdx;
        }
        else {
            // if the pressure difference is zero, we chose the DOF which has the
            // larger volume associated to it as upstream DOF
            if (Vin > Vex) {
                upIdx = interiorDofIdx;
                dnIdx = exteriorDofIdx;
            }
            else if (Vin < Vex) {
                upIdx = exteriorDofIdx;
                dnIdx = interiorDofIdx;
            }
            else {
                assert(Vin == Vex);
                // if the volumes are also equal, we pick the DOF which exhibits the
                // smaller global index
                if (globalIndexIn < globalIndexEx) {
                    upIdx = interiorDofIdx;
                    dnIdx = exteriorDofIdx;
                }
                else {
                    upIdx  = exteriorDofIdx;
                    dnIdx = interiorDofIdx;
                }
            }
        }

        // apply the threshold pressure for the intersection. note that the concept
        // of threshold pressure is a quite big hack that only makes sense for ECL
        // datasets. (and even there, its physical justification is quite
        // questionable IMO.)
        if (thpres > 0.0) {
            if (std::abs(Toolbox::value(pressureDifference)) > thpres) {
                if (pressureDifference < 0.0)
                    pressureDifference += thpres;
                else
                    pressureDifference -= thpres;
            }
            else {
                pressureDifference = 0.0;
            }
        }
    }

protected:
    /*!
     * \brief Update the required gradients for interior faces
     */
    void calculateGradients_(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    {
        Valgrind::SetUndefined(*this);

        volumeAndPhasePressureDifferences(upIdx_ , dnIdx_, volumeFlux_, pressureDifference_, elemCtx, scvfIdx, timeIdx);
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
        const auto& scvf = elemCtx.stencil(timeIdx).boundaryFace(scvfIdx);
        const Scalar faceArea = scvf.area();
        const Scalar zEx = scvf.integrationPos()[dimWorld - 1];
        const auto& problem = elemCtx.problem();
        const unsigned globalSpaceIdx = elemCtx.globalSpaceIndex(0, timeIdx);
        const auto& intQuantsIn = elemCtx.intensiveQuantities(0, timeIdx);

        calculateBoundaryGradients_(problem,
                                    globalSpaceIdx,
                                    intQuantsIn,
                                    scvfIdx,
                                    faceArea,
                                    zEx,
                                    exFluidState,
                                    upIdx_,
                                    dnIdx_,
                                    volumeFlux_,
                                    pressureDifference_);

        // Treating solvent here and not in the static method, since that would require more
        // extensive refactoring. It means that the TpfaLinearizer will not support bcs for solvent until this is
        // addressed.
        if constexpr (enableSolvent) {
            if (upIdx_[gasPhaseIdx] == 0) {
                const Scalar trans = problem.transmissibilityBoundary(globalSpaceIdx, scvfIdx);
                const Scalar transModified = trans * Toolbox::value(intQuantsIn.rockCompTransMultiplier());
                const auto solventFlux = pressureDifference_[gasPhaseIdx] * intQuantsIn.mobility(gasPhaseIdx) * (-transModified/faceArea);
                asImp_().setSolventVolumeFlux(solventFlux);
            } else {
                asImp_().setSolventVolumeFlux(0.0);
            }
        }
    }

public:
    /*!
     * \brief Update the required gradients for boundary faces
     */
    template <class Problem, class FluidState, class EvaluationContainer>
    static void calculateBoundaryGradients_(const Problem& problem,
                                            const unsigned globalSpaceIdx,
                                            const IntensiveQuantities& intQuantsIn,
                                            const unsigned bfIdx,
                                            const double faceArea,
                                            const double zEx,
                                            const FluidState& exFluidState,
                                            std::array<short, numPhases>& upIdx,
                                            std::array<short, numPhases>& dnIdx,
                                            EvaluationContainer& volumeFlux,
                                            EvaluationContainer& pressureDifference)
    {

        bool enableBoundaryMassFlux = problem.nonTrivialBoundaryConditions();
        if (!enableBoundaryMassFlux)
            return;

        Scalar trans = problem.transmissibilityBoundary(globalSpaceIdx, bfIdx);

        // estimate the gravity correction: for performance reasons we use a simplified
        // approach for this flux module that assumes that gravity is constant and always
        // acts into the downwards direction. (i.e., no centrifuge experiments, sorry.)
        Scalar g = problem.gravity()[dimWorld - 1];

        // this is quite hacky because the dune grid interface does not provide a
        // cellCenterDepth() method (so we ask the problem to provide it). The "good"
        // solution would be to take the Z coordinate of the element centroids, but since
        // ECL seems to like to be inconsistent on that front, it needs to be done like
        // here...
        Scalar zIn = problem.dofCenterDepth(globalSpaceIdx);

        // the distances from the DOF's depths. (i.e., the additional depth of the
        // exterior DOF)
        Scalar distZ = zIn - zEx;

        for (unsigned phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            // do the gravity correction: compute the hydrostatic pressure for the
            // integration position
            const Evaluation& rhoIn = intQuantsIn.fluidState().density(phaseIdx);
            const auto& rhoEx = exFluidState.density(phaseIdx);
            Evaluation rhoAvg = (rhoIn + rhoEx)/2;

            const Evaluation& pressureInterior = intQuantsIn.fluidState().pressure(phaseIdx);
            Evaluation pressureExterior = exFluidState.pressure(phaseIdx);
            pressureExterior += rhoAvg*(distZ*g);

            pressureDifference[phaseIdx] = pressureExterior - pressureInterior;

            // decide the upstream index for the phase. for this we make sure that the
            // degree of freedom which is regarded upstream if both pressures are equal
            // is always the same: if the pressure is equal, the DOF with the lower
            // global index is regarded to be the upstream one.
            const unsigned interiorDofIdx = 0; // Valid only for cell-centered FV.
            if (pressureDifference[phaseIdx] > 0.0) {
                upIdx[phaseIdx] = -1;
                dnIdx[phaseIdx] = interiorDofIdx;
            }
            else {
                upIdx[phaseIdx] = interiorDofIdx;
                dnIdx[phaseIdx] = -1;
            }

            Evaluation transModified = trans;

            if (upIdx[phaseIdx] == interiorDofIdx) {

                // this is slightly hacky because in the automatic differentiation case, it
                // only works for the element centered finite volume method. for ebos this
                // does not matter, though.
                const auto& up = intQuantsIn;

                // deal with water induced rock compaction
                const Scalar transMult = Toolbox::value(up.rockCompTransMultiplier());
                transModified *= transMult;

                volumeFlux[phaseIdx] =
                    pressureDifference[phaseIdx]*up.mobility(phaseIdx)*(-transModified/faceArea);
            }
            else {
                // compute the phase mobility using the material law parameters of the
                // interior element. TODO: this could probably be done more efficiently
                const auto& matParams = problem.materialLawParams(globalSpaceIdx);
                std::array<typename FluidState::Scalar,numPhases> kr;
                MaterialLaw::relativePermeabilities(kr, matParams, exFluidState);

                const auto& mob = kr[phaseIdx]/exFluidState.viscosity(phaseIdx);
                volumeFlux[phaseIdx] =
                    pressureDifference[phaseIdx]*mob*(-transModified/faceArea);
            }
        }
    }

protected:

    /*!
     * \brief Update the volumetric fluxes for all fluid phases on the interior faces of the context
     */
    void calculateFluxes_(const ElementContext&, unsigned, unsigned)
    { }

    void calculateBoundaryFluxes_(const ElementContext&, unsigned, unsigned)
    {}

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    // the volumetric flux of all phases [m^3/s]
    Evaluation volumeFlux_[numPhases];

    // the difference in effective pressure between the exterior and the interior degree
    // of freedom [Pa]
    Evaluation pressureDifference_[numPhases];

    // the local indices of the interior and exterior degrees of freedom
    std::array<short, numPhases> upIdx_;
    std::array<short, numPhases> dnIdx_;
 };

} // namespace Opm

#endif
