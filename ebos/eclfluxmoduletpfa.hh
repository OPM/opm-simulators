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
#ifndef EWOMS_ECL_FLUX_TPFA_MODULE_HH
#define EWOMS_ECL_FLUX_TPFA_MODULE_HH

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/utils/signum.hh>

#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Opm {

template <class TypeTag>
class EclTransIntensiveQuantities;

template <class TypeTag>
class EclTransExtensiveQuantities;

template <class TypeTag>
class EclTransBaseProblem;

// /*!
//  * \ingroup EclBlackOilSimulator
//  * \brief Specifies a flux module which uses ECL transmissibilities.
//  */
// template <class TypeTag>
// struct EclTransFluxModule
// {
//     typedef EclTransIntensiveQuantities<TypeTag> FluxIntensiveQuantities;
//     typedef EclTransExtensiveQuantities<TypeTag> FluxExtensiveQuantities;
//     typedef EclTransBaseProblem<TypeTag> FluxBaseProblem;

//     /*!
//      * \brief Register all run-time parameters for the flux module.
//      */
//     static void registerParameters()
//     { }
// };

// /*!
//  * \ingroup EclBlackOilSimulator
//  * \brief Provides the defaults for the parameters required by the
//  *        transmissibility based volume flux calculation.
//  */
// template <class TypeTag>
// class EclTransBaseProblem
// { };

// /*!
//  * \ingroup EclBlackOilSimulator
//  * \brief Provides the intensive quantities for the ECL flux module
//  */
// template <class TypeTag>
// class EclTransIntensiveQuantities
// {
//     using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
// protected:
//     void update_(const ElementContext&, unsigned, unsigned)
//     { }
// };

/*!
 * \ingroup EclBlackOilSimulator
 * \brief Provides the ECL flux module
 */
template <class TypeTag>
class EclTransExtensiveQuantitiesTPFA
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

    typedef MathToolbox<Evaluation> Toolbox;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> EvalDimVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;

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
    {
        //throw std::invalid_argument("Should not acces volume flux for eclmoduletpfa");
        return volumeFlux_[phaseIdx];
    }

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
        //throw std::invalid_argument("No upstreamIndex");
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
        //throw std::invalid_argument("No downstream index");
        assert(phaseIdx < numPhases);

        return dnIdx_[phaseIdx];
    }

    void updateSolvent(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    { asImp_().updateVolumeFluxTrans(elemCtx, scvfIdx, timeIdx); }

    void updatePolymer(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    { asImp_().updateShearMultipliers(elemCtx, scvfIdx, timeIdx); }


public:
    template<class EvalType>
    static void calculatePhasePressureDiff_(short& upIdx,
                                            short& dnIdx,
                                            EvalType& pressureDifference,
                                            const IntensiveQuantities& intQuantsIn,
                                            const IntensiveQuantities& intQuantsEx,
                                            const unsigned timeIdx,
                                            const unsigned phaseIdx,
                                            const unsigned interiorDofIdx,
                                            const unsigned exteriorDofIdx,
                                            const Scalar& Vin,
                                            const Scalar& Vex,
                                            const unsigned& globalIndexIn,
                                            const unsigned& globalIndexEx,
                                            const Scalar& distZg,
                                            const Scalar& thpres
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


    static void volumeAndPhasePressureDifferences(short (&upIdx)[numPhases] ,
                                           Evaluation (&volumeFlux)[numPhases],
                                           Evaluation (&pressureDifferences)[numPhases],
                                           const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    {

        // Valgrind::SetUndefined(*this);

        const auto& problem = elemCtx.problem();
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);
        unsigned interiorDofIdx = scvf.interiorIndex();
        unsigned exteriorDofIdx = scvf.exteriorIndex();
        const auto& globalIndexIn = stencil.globalSpaceIndex(interiorDofIdx);
        const auto& globalIndexEx = stencil.globalSpaceIndex(exteriorDofIdx);

        assert(interiorDofIdx != exteriorDofIdx);

        // unsigned I = stencil.globalSpaceIndex(interiorDofIdx_);
        // unsigned J = stencil.globalSpaceIndex(exteriorDofIdx_);
        Scalar Vin = elemCtx.dofVolume(interiorDofIdx, /*timeIdx=*/0);
        Scalar Vex = elemCtx.dofVolume(exteriorDofIdx, /*timeIdx=*/0);

        Scalar trans = problem.transmissibility(elemCtx, interiorDofIdx, exteriorDofIdx);
        Scalar faceArea = scvf.area();
        Scalar thpres = problem.thresholdPressure(globalIndexIn, globalIndexEx);

        // estimate the gravity correction: for performance reasons we use a simplified
        // approach for this flux module that assumes that gravity is constant and always
        // acts into the downwards direction. (i.e., no centrifuge experiments, sorry.)
        constexpr Scalar g = 9.8;

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

        for (unsigned phaseIdx=0; phaseIdx < numPhases; phaseIdx++) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;
            short dnIdx;
            calculatePhasePressureDiff_(upIdx[phaseIdx],
                                        dnIdx,
                                        pressureDifferences[phaseIdx],
                                        intQuantsIn,
                                        intQuantsEx,
                                        timeIdx,//input
                                        phaseIdx,//input
                                        interiorDofIdx,//input
                                        exteriorDofIdx,//intput
                                        Vin,
                                        Vex,
                                        globalIndexIn,
                                        globalIndexEx,
                                        distZ*g,
                                        thpres);
            if (pressureDifferences[phaseIdx] == 0) {
                volumeFlux[phaseIdx] = 0.0;
                continue;
            }

            const IntensiveQuantities& up = (upIdx[phaseIdx] == interiorDofIdx) ? intQuantsIn : intQuantsEx;
            // TODO: should the rock compaction transmissibility multiplier be upstreamed
            // or averaged? all fluids should see the same compaction?!
            const Evaluation& transMult = up.rockCompTransMultiplier();

            if (upIdx[phaseIdx] == interiorDofIdx)
                volumeFlux[phaseIdx] =
                    pressureDifferences[phaseIdx]*up.mobility(phaseIdx)*transMult*(-trans/faceArea);
            else
                volumeFlux[phaseIdx] =
                    pressureDifferences[phaseIdx]*(Toolbox::value(up.mobility(phaseIdx))*Toolbox::value(transMult)*(-trans/faceArea));
        }
    }

protected:
    /*!
     * \brief Update the required gradients for interior faces
     */
    void calculateGradients_(const ElementContext& elemCtx, unsigned scvfIdx, unsigned timeIdx)
    {
        //throw std::invalid_argument("No calculateGradients_");

        Valgrind::SetUndefined(*this);

        const auto& problem = elemCtx.problem();
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        interiorDofIdx_ = scvf.interiorIndex();
        exteriorDofIdx_ = scvf.exteriorIndex();
        assert(interiorDofIdx_ != exteriorDofIdx_);

        unsigned I = stencil.globalSpaceIndex(interiorDofIdx_);
        unsigned J = stencil.globalSpaceIndex(exteriorDofIdx_);

        Scalar trans = problem.transmissibility(elemCtx, interiorDofIdx_, exteriorDofIdx_);
        Scalar faceArea = scvf.area();
        Scalar thpres = problem.thresholdPressure(I, J);

        // estimate the gravity correction: for performance reasons we use a simplified
        // approach for this flux module that assumes that gravity is constant and always
        // acts into the downwards direction. (i.e., no centrifuge experiments, sorry.)
        Scalar g = elemCtx.problem().gravity()[dimWorld - 1];

        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx_, timeIdx);
        const auto& intQuantsEx = elemCtx.intensiveQuantities(exteriorDofIdx_, timeIdx);

        // this is quite hacky because the dune grid interface does not provide a
        // cellCenterDepth() method (so we ask the problem to provide it). The "good"
        // solution would be to take the Z coordinate of the element centroids, but since
        // ECL seems to like to be inconsistent on that front, it needs to be done like
        // here...
        Scalar zIn = problem.dofCenterDepth(elemCtx, interiorDofIdx_, timeIdx);
        Scalar zEx = problem.dofCenterDepth(elemCtx, exteriorDofIdx_, timeIdx);

        // the distances from the DOF's depths. (i.e., the additional depth of the
        // exterior DOF)
        Scalar distZ = zIn - zEx;

        Scalar Vin = elemCtx.dofVolume(interiorDofIdx_, /*timeIdx=*/0);
        Scalar Vex = elemCtx.dofVolume(exteriorDofIdx_, /*timeIdx=*/0);

        for (unsigned phaseIdx = 0; phaseIdx < numPhases; phaseIdx++) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;
            calculatePhasePressureDiff_(upIdx_[phaseIdx],
                                        dnIdx_[phaseIdx],
                                        pressureDifference_[phaseIdx],
                                        intQuantsIn,
                                        intQuantsEx,
                                        timeIdx, // input
                                        phaseIdx, // input
                                        interiorDofIdx_, // input
                                        exteriorDofIdx_, // intput
                                        Vin,
                                        Vex,
                                        I,
                                        J,
                                        distZ * g,
                                        thpres);
            if (pressureDifference_[phaseIdx] == 0) {
                volumeFlux_[phaseIdx] = 0.0;
                continue;
            }
            const IntensiveQuantities& up = (upIdx_[phaseIdx] == interiorDofIdx_) ? intQuantsIn : intQuantsEx;
            // TODO: should the rock compaction transmissibility multiplier be upstreamed
            // or averaged? all fluids should see the same compaction?!
            const Evaluation& transMult = up.rockCompTransMultiplier();

            if (upIdx_[phaseIdx] == interiorDofIdx_)
                volumeFlux_[phaseIdx] =
                    pressureDifference_[phaseIdx]*up.mobility(phaseIdx)*transMult*(-trans/faceArea);
            else
                volumeFlux_[phaseIdx] =
                    pressureDifference_[phaseIdx]*(Toolbox::value(up.mobility(phaseIdx))*Toolbox::value(transMult)*(-trans/faceArea));
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
        // throw std::invalid_argument("No calculateGradients for boundary");
        const auto& problem = elemCtx.problem();

        bool enableBoundaryMassFlux = problem.nonTrivialBoundaryConditions();
        if (!enableBoundaryMassFlux)
            return;

        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.boundaryFace(scvfIdx);

        interiorDofIdx_ = scvf.interiorIndex();

        Scalar trans = problem.transmissibilityBoundary(elemCtx, scvfIdx);
        Scalar faceArea = scvf.area();

        // estimate the gravity correction: for performance reasons we use a simplified
        // approach for this flux module that assumes that gravity is constant and always
        // acts into the downwards direction. (i.e., no centrifuge experiments, sorry.)
        Scalar g = elemCtx.problem().gravity()[dimWorld - 1];

        const auto& intQuantsIn = elemCtx.intensiveQuantities(interiorDofIdx_, timeIdx);

        // this is quite hacky because the dune grid interface does not provide a
        // cellCenterDepth() method (so we ask the problem to provide it). The "good"
        // solution would be to take the Z coordinate of the element centroids, but since
        // ECL seems to like to be inconsistent on that front, it needs to be done like
        // here...
        Scalar zIn = problem.dofCenterDepth(elemCtx, interiorDofIdx_, timeIdx);
        Scalar zEx = scvf.integrationPos()[dimWorld - 1];

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

            Evaluation transModified = trans;

            short upstreamIdx = upstreamIndex_(phaseIdx);
            if (upstreamIdx == interiorDofIdx_) {

                // this is slightly hacky because in the automatic differentiation case, it
                // only works for the element centered finite volume method. for ebos this
                // does not matter, though.
                const auto& up = elemCtx.intensiveQuantities(upstreamIdx, timeIdx);

                // deal with water induced rock compaction
                const double transMult = Toolbox::value(up.rockCompTransMultiplier());
                transModified *= transMult;

                volumeFlux_[phaseIdx] =
                    pressureDifference_[phaseIdx]*up.mobility(phaseIdx)*(-transModified/faceArea);

                if (enableSolvent && phaseIdx == gasPhaseIdx)
                    asImp_().setSolventVolumeFlux( pressureDifference_[phaseIdx]*up.solventMobility()*(-transModified/faceArea));
            }
            else {
                // compute the phase mobility using the material law parameters of the
                // interior element. TODO: this could probably be done more efficiently
                const auto& matParams =
                    elemCtx.problem().materialLawParams(elemCtx,
                                                        interiorDofIdx_,
                                                        /*timeIdx=*/0);
                typename FluidState::Scalar kr[numPhases];
                MaterialLaw::relativePermeabilities(kr, matParams, exFluidState);

                const auto& mob = kr[phaseIdx]/exFluidState.viscosity(phaseIdx);
                volumeFlux_[phaseIdx] =
                    pressureDifference_[phaseIdx]*mob*(-transModified/faceArea);

                // Solvent inflow is not yet supported
                if (enableSolvent && phaseIdx == gasPhaseIdx)
                    asImp_().setSolventVolumeFlux(0.0);
            }
        }
    }

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
    unsigned short interiorDofIdx_;
    unsigned short exteriorDofIdx_;
    short upIdx_[numPhases];
    short dnIdx_[numPhases];
};

}// namespace Opm

#endif
