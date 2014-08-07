/*
  Copyright (C) 2014 by Andreas Lauser

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
/**
 * \file
 *
 * \copydoc Ewoms::EclPeacemanWell
 */
#ifndef EWOMS_ECL_PEACEMAN_WELL_HH
#define EWOMS_ECL_PEACEMAN_WELL_HH

#include <opm/core/utility/PropertySystem.hpp>
#include <opm/core/utility/Average.hpp>

#include <dune/common/version.hh>

namespace Opm {
namespace Properties {
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Discretization);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(ElementContext);
NEW_PROP_TAG(RateVector);
NEW_PROP_TAG(GridView);
}}

namespace Ewoms {
/*!
 * \brief The well model of Peaceman.
 *
 * This class is tailored for the element centered finite volume
 * discretization, assumes a vertical borehole and is intended to be
 * used by the EclWellManager.
 *
 * See:
 *
 * Z. Chen, G. Huan, Y. Ma: Computational Methods for Multiphase
 * Flows in Porous Media, 1st edition, SIAM, 2006, pp. 445-446
 *
 * and
 *
 * D. W. Peaceman: Interpretation of well-block pressures in numerical
 * reservoir simulation, The 52nd Annual SPE Fall Technical Conference
 * and Exhibition, Denver, CO., 1977
 */
template <class TypeTag>
class EclPeacemanWell
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Discretization) Discretization;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    // the dimension of the simulator's world
    static const int dimWorld = GridView::dimensionworld;

    // convenient access to the number of phases and the number of
    // components
    static const int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);
    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);

    // convenient access to the phase and component indices. If the
    // compiler bails out here, you're probably using an incompatible
    // fluid system. Try to use the Ewoms::BlackOilFluidSystem...
    static const int gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static const int oilPhaseIdx = FluidSystem::oilPhaseIdx;
    static const int waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static const int oilCompIdx = FluidSystem::oilCompIdx;
    static const int waterCompIdx = FluidSystem::waterCompIdx;
    static const int gasCompIdx = FluidSystem::gasCompIdx;

    typedef Opm::CompositionalFluidState<Scalar, FluidSystem, /*storeEnthalpy=*/false> FluidState;

    // all quantities that need to be stored per degree of freedom that intersects the
    // well.
    struct DofVariables {
        // The connection transmissibility factor to be used for a given
        // DOF.  This object doubles by defining which DOFs are part of
        // the well.
        Scalar connectionTransmissibilityFactor;

        // the volumetric reservoir rates for each fluid phase and each
        // degree of freedom. This is calculated at the beginning of each
        // iteration and used to impose rate limits.
        std::array<Scalar, numPhases> unconstraintRates;

        // the volumetric surface rates for each fluid phase and each
        // degree of freedom. This is calculated at the beginning of each
        // iteration and used to impose rate limits.
        std::array<Scalar, numPhases> unconstraintSurfaceRates;

        // the effective size of an element in each direction. This is defined as the
        // distance of the face centers along the respective axis.
        std::array<Scalar, dimWorld> effectiveSize;

        // the radius of the well for the given degree of freedom
        Scalar boreholeRadius;

        // The skin factor of the well at the given degree of freedom
        Scalar skinFactor;
    };

    // some safety checks/caveats
    static_assert(std::is_same<Discretization, EcfvDiscretization<TypeTag> >::value,
                  "The Peaceman well model is only implemented for the "
                  "element-centered finite volume discretization!");
    static_assert(dimWorld == 3,
                  "The Peaceman well model is only implemented for 3D grids!");

public:
    enum ControlMode {
        BottomHolePressure,
        TopHolePressure,
        VolumetricSurfaceRate,
        VolumetricReservoirRate
    };

    enum WellType {
        Undefined,
        Injector,
        Producer
    };

    EclPeacemanWell(const Simulator &simulator)
        : simulator_(simulator)
    {
        // set the composition of the injected fluids based. If
        // somebody is stupid enough to inject oil, we assume he wants
        // to loose his fortune on dry oil...
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            for (int compIdx = 0; compIdx < numComponents; ++ compIdx)
                injectionFluidState_.setMoleFraction(phaseIdx, compIdx, 0.0);
        injectionFluidState_.setMoleFraction(gasPhaseIdx, gasCompIdx, 1.0);
        injectionFluidState_.setMoleFraction(waterPhaseIdx, waterCompIdx, 1.0);
        injectionFluidState_.setMoleFraction(oilPhaseIdx, oilCompIdx, 1.0);

        // set the temperature to 25 deg C, just so that it is set
        injectionFluidState_.setTemperature(273.15 + 25);
    }

    /*!
     * \brief Begin the specification of the well.
     *
     * The specification process is the following:
     *
     * beginSpec()
     * setName("FOO");
     * // add degrees of freedom to the well
     * for (dof in wellDofs)
     *    addDof(dof);
     *
     *    // set the radius of the well at the dof [m].
     *    // optional, if not specified, it is assumed to be 0.1524m
     *    setRadius(dof, someRadius);
     *
     *    // set the skin factor of the well.
     *    // optional, if not specified, it is assumed to be 0
     *    setSkinFactor(dof, someSkinFactor);
     * endSpec()
     *
     * // specify the phase which is supposed to be injected. (Optional,
     * // if unspecified, the well will throw an
     * // exception if it would inject something.)
     * setInjectedPhaseIndex(phaseIdx);
     *
     * // set maximum production rate at reservoir conditions
     * // (kg/s, optional, if not specified, the well is assumed to be
     * // shut for production)
     * setMaximumReservoirRate(someMassRate);
     *
     * // set maximum injection rate at reservoir conditions
     * // (kg/s, optional, if not specified, the well is assumed to be
     * // shut for injection)
     * setMinmumReservoirRate(someMassRate);
     *
     * // set the relative weight of the mass rate of a fluid phase.
     * // (Optional, if unspecified each phase exhibits a weight of 1)
     * setPhaseWeight(phaseIdx, someWeight);
     *
     * // set maximum production rate at surface conditions
     * // (kg/s, optional, if not specified, the well is assumed to be
     * // not limited by the surface rate)
     * setMaximumSurfaceRate(someMassRate);
     *
     * // set maximum production rate at surface conditions
     * // (kg/s, optional, if not specified, the well is assumed to be
     * // not limited by the surface rate)
     * setMinimumSurfaceRate(someMassRate);
     *
     * // set the minimum pressure at the bottom of the well (Pa,
     * // optional, if not specified, the well is assumes it estimates
     * // the bottom hole pressure based on the top hole pressure
     * // assuming hydrostatic conditions.)
     * setMinimumBottomHolePressure(somePressure);
     *
     * // set the pressure at the top of the well (Pa,
     * // optional, if not specified, the top hole pressure is
     * // assumed to be 1 bar)
     * setTopHolePressure(somePressure);
     *
     * // set the control mode of the well [m].
     * // optional, if not specified, it is assumed to be "BottomHolePressure"
     * setControlMode(Well::TopHolePressure);
     *
     * // set the top hole pressure of the well [Pa]
     * // only required if the control mode is "TopHolePressure"
     * setTopHolePressure(1e5);
     */
    void beginSpec()
    {
        // this is going to be increased by any realistic grid. Shall we bet?
        bottomDepth_ = -1e100;
        bottomDofGlobalIdx_ = -1;

        // By default, take the bottom hole pressure as a given
        controlMode_ = ControlMode::BottomHolePressure;

        // use one bar for the default bottom and top hole
        // pressures. For the bottom hole pressure, this is probably
        // off by at least one magnitude...
        targetBhp_ = 1e5;
        targetThp_ = 1e5;

        // By default, all fluids exhibit the weight 1.0
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            volumetricWeight_[phaseIdx] = 1.0;

        wellType_ = Undefined;

        wellTotalVolume_ = 0.0;
    }

    /*!
     * \brief Set the relative weight of the volumetric phase rates.
     */
    void setVolumetricPhaseWeights(Scalar oilWeight, Scalar gasWeight, Scalar waterWeight)
    {
        volumetricWeight_[oilPhaseIdx] = oilWeight;
        volumetricWeight_[gasPhaseIdx] = gasWeight;
        volumetricWeight_[waterPhaseIdx] = waterWeight;
    }

    /*!
     * \brief Return the human-readable name of the well
     *
     * Well, let's say "readable by some humans".
     */
    const std::string &name() const
    { return name_; }

    /*!
     * \brief Set the human-readable name of the well
     */
    void setName(const std::string &newName)
    { name_ = newName; }

    /*!
     * \brief Add a degree of freedom to the well.
     */
    template <class Context>
    void addDof(const Context &context,
                int dofIdx,
                Scalar connectionTransmissibilityFactor = 1.0)
    {
        int globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        if (applies(globalDofIdx))
            // we already have this DOF in the well!
            return;

        const auto &dofPos = context.pos(dofIdx, /*timeIdx=*/0);

        DofVariables &dofVars = dofVariables_[globalDofIdx];
        dofVars.connectionTransmissibilityFactor = connectionTransmissibilityFactor;

        wellTotalVolume_ += context.model().dofTotalVolume(globalDofIdx);

        // determine the size of the element
        dofVars.effectiveSize.fill(0.0);

        // we assume all elements to be hexahedrons!
        assert(context.element().template count</*codim=*/dimWorld>() == 8);

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2,3)
        const auto &refElem = Dune::ReferenceElements<Scalar, /*dim=*/3>::cube();
#else
        const auto &refElem = Dune::GenericReferenceElements<Scalar, /*dim=*/3>::cube();
#endif

        // determine the current element's effective size
        const auto &elem = context.element();
        int faceIdx = 0;
        int numFaces = refElem.size(/*codim=*/1);
        for (; faceIdx < numFaces; ++faceIdx) {
            const auto &faceCenterLocal = refElem.position(faceIdx, /*codim=*/1);
            const auto &faceCenter = elem.geometry().global(faceCenterLocal);

            switch (faceIdx) {
            case 0:
                dofVars.effectiveSize[0] -= faceCenter[0];
                break;
            case 1:
                dofVars.effectiveSize[0] += faceCenter[0];
                break;
            case 2:
                dofVars.effectiveSize[1] -= faceCenter[1];
                break;
            case 3:
                dofVars.effectiveSize[1] += faceCenter[1];
                break;
            case 4:
                dofVars.effectiveSize[2] -= faceCenter[2];
                break;
            case 5:
                dofVars.effectiveSize[2] += faceCenter[2];
                break;
            }
        }

        // default borehole radius: 1/2 foot
        dofVars.boreholeRadius = 0.3048/2;

        // default skin factor: 0
        dofVars.skinFactor = 0;

        // we assume that the z-coordinate represents depth (and not
        // height) here...
        if (dofPos[2] > bottomDepth_) {
            bottomDofGlobalIdx_ = globalDofIdx;
            bottomDepth_ = dofPos[2];
        }
    }

    /*!
     * \brief Finalize the specification of the borehole.
     */
    void endSpec()
    {
        const auto& comm = simulator_.gridView().comm();

        // determine the maximum depth of the well over all processes
        bottomDepth_ = comm.max(bottomDepth_);

        // the total volume of the well must also be summed over all processes
        wellTotalVolume_ = comm.sum(wellTotalVolume_);
    }

    /*!
     * \brief Set the control mode of the well.
     *
     * This specifies which quantities are assumed to be externally
     * given and which must be calculated based on those.
     */
    void setControlMode(ControlMode controlMode)
    { controlMode_ = controlMode; }

    /*!
     * \brief Set the type of the well (i.e., injector or producer).
     */
    void setWellType(WellType wellType)
    { wellType_ = wellType; }

    /*!
     * \brief Returns the type of the well (i.e., injector or producer).
     */
    WellType wellType() const
    { return wellType_; }

    /*!
     * \brief Set the index of fluid phase to be injected.
     *
     * This is only relevant if the well type is an injector.
     */
    void setInjectedPhaseIndex(int injPhaseIdx)
    { injectedPhaseIdx_ = injPhaseIdx; }

    /*!
     * \brief The Z-coordinate of the well's deepest degree of freedom
     */
    Scalar bottomDepth() const
    { return bottomDepth_; }

    /*!
     * \brief Set whether the well should be closed or not
     */
    void setOpen(bool yesno)
    { isOpen_ = yesno; }

    /*!
     * \brief Return whether the well is closed or not
     */
    bool isOpen() const
    { return isOpen_; }

    /*!
     * \brief Return true iff a degree of freedom is directly affected
     *        by the well
     */
    bool applies(int globalDofIdx) const
    { return dofVariables_.count(globalDofIdx) > 0; }

    /*!
     * \brief Set the maximum bottom hole pressure [Pa] of the well.
     */
    void setTargetBottomHolePressure(Scalar val)
    { targetBhp_ = val; }

    /*!
     * \brief Set the top hole pressure [Pa] of the well.
     */
    void setTargetTopHolePressure(Scalar val)
    { targetThp_ = val; }

    /*!
     * \brief Set the maximum combined rate of the fluids at the surface.
     */
    void setMaximumSurfaceRate(Scalar value)
    { maximumSurfaceRate_ = value; }

    /*!
     * \brief Set the maximum combined rate of the fluids at the surface.
     */
    void setMaximumReservoirRate(Scalar value)
    { maximumReservoirRate_ = value; }

    /*!
     * \brief Set the skin factor of the well
     */
    template <class Context>
    void setSkinFactor(const Context &context, int dofIdx, Scalar value)
    {
        int globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        dofVariables_[globalDofIdx].skinFactor = value;
    }

    /*!
     * \brief Set the borehole radius of the well
     */
    template <class Context>
    void setRadius(const Context &context, int dofIdx, Scalar value)
    {
        int globalDofIdx = context.globalSpaceIndex(dofIdx, /*timeIdx=*/0);
        dofVariables_[globalDofIdx].boreholeRadius = value;
    }

    /*!
     * \brief Informs the well that a time step has just begun.
     */
    void beginTimeStep()
    {
        // nothing to do, yet
    }

    /*!
     * \brief Informs the well that an iteration has just begun.
     *
     * The beginIteration*() methods, the well calculates the bottom
     * and top hole pressures, the actual unconstraint production and
     * injection rates, etc. The callback is split into three parts as
     * this arrangement avoids iterating over the whole grid and to
     * re-calculate the volume variables for each well.
     *
     * This is supposed to prepare the well object to do the
     * computations which are required to do the DOF specific
     * things.
     */
    void beginIterationPreProcess()
    {
        // calculate the bottom hole pressure to be actually used
        if (controlMode_ == ControlMode::TopHolePressure) {
            // assume a density of 650 kg/m^3 for the bottom hole pressure
            // calculation
            Scalar rho = 650.0;
            effectiveBottomHolePressure_ = targetThp_ + rho*bottomDepth_;

            // set the maximum rates to unlimited
            maximumReservoirRate_ = 1e100;
            maximumSurfaceRate_ = 1e100;
        }
        else if (controlMode_ == ControlMode::BottomHolePressure) {
            effectiveBottomHolePressure_ = targetBhp_;

            // set the maximum rates to unlimited
            maximumReservoirRate_ = 1e100;
            maximumSurfaceRate_ = 1e100;
        }
        else {
            // calculate the bottom hole pressure limit from the top-hole pressure
            // limit. this is a HACK since the effective density must be given and is
            // assumed to be constant...
            Scalar rhoEff = 650; // kg/m^3
            Scalar bhpFromThp = targetThp_ + rhoEff*bottomDepth_;

            if (wellType_ == WellType::Injector)
                effectiveBottomHolePressure_ = std::max(bhpFromThp, targetBhp_);
            else if (wellType_ == WellType::Producer)
                effectiveBottomHolePressure_ = std::min(bhpFromThp, targetBhp_);
        }

        if (wellType_ == WellType::Injector)
            observedBhp_ = - 1e100;
        else if (wellType_ == WellType::Producer)
            observedBhp_ = 1e100;

        // make it very likely that we screw up if we control for {surface,reservoir}
        // rate, but depend on the {reservoir,surface} rate somewhere...
        if (controlMode_ == ControlMode::VolumetricSurfaceRate)
            maximumReservoirRate_ = std::numeric_limits<Scalar>::quiet_NaN();
        else if (controlMode_ == ControlMode::VolumetricReservoirRate)
            maximumSurfaceRate_ = std::numeric_limits<Scalar>::quiet_NaN();

        // reset the unconstraint rates for the complete well. ("unconstraint" ==
        // unaffected by rate limits.)
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            unconstraintReservoirRates_[phaseIdx] = 0.0;
            unconstraintSurfaceRates_[phaseIdx] = 0.0;
            currentSurfaceRates_[phaseIdx] = 0.0;
        }
    }

    /*!
     * \brief Do the DOF specific part at the beginning of each iteration
     */
    template <class Context>
    void beginIterationAccumulate(Context &context, int timeIdx)
    {
        std::array<Scalar, numPhases> reservoirVolRates;
        RateVector massRate;
        for (int dofIdx = 0; dofIdx < context.numPrimaryDof(timeIdx); ++dofIdx) {
            int globalDofIdx = context.globalSpaceIndex(dofIdx, timeIdx);
            if (!applies(globalDofIdx))
                continue;

            const DofVariables &dofVars = dofVariables_.at(globalDofIdx);

            computeUnconstraintVolumetricDofRates_(reservoirVolRates, dofVars, context, dofIdx, timeIdx);

            dofVariables_[globalDofIdx].unconstraintRates = reservoirVolRates;
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                massRate.setVolumetricRate(context.intensiveQuantities(dofIdx, timeIdx).fluidState(),
                                           phaseIdx,
                                           reservoirVolRates[phaseIdx]);

                unconstraintReservoirRates_[phaseIdx] += reservoirVolRates[phaseIdx];
            }

            std::array<Scalar, numPhases> dofSurfaceRate;
            const auto& intQuants = context.intensiveQuantities(dofIdx, timeIdx);
            computeSurfaceRates_(dofSurfaceRate,
                                 reservoirVolRates,
                                 intQuants.fluidState());

            dofVariables_[globalDofIdx].unconstraintSurfaceRates = dofSurfaceRate;
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
                unconstraintSurfaceRates_[phaseIdx] += dofSurfaceRate[phaseIdx];

            if (globalDofIdx == bottomDofGlobalIdx_)
                observedBhp_ = intQuants.fluidState().pressure(oilPhaseIdx);
        }
    }

    /*!
     * \brief Informs the well that an iteration has just begun.
     *
     * This is the post-processing part which uses the results of the
     * accumulation callback.
     */
    void beginIterationPostProcess()
    {
        const auto& comm = simulator_.gridView().comm();
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            unconstraintReservoirRates_[phaseIdx] = comm.sum(unconstraintReservoirRates_[phaseIdx]);
            unconstraintSurfaceRates_[phaseIdx] = comm.sum(unconstraintSurfaceRates_[phaseIdx]);
        }

        // determine the grid-global observed bottom hole pressure
        if (wellType_ == Producer)
            observedBhp_ = comm.min(observedBhp_);
        else if (wellType_ == Injector)
            observedBhp_ = comm.max(observedBhp_);

        // determine the rate-limited surface rates
        Scalar alpha = 1.0;
        if (controlMode_ == VolumetricSurfaceRate) {
            Scalar weightedSurfRate = computeWeightedRate_(unconstraintSurfaceRates_);
            if (std::abs(weightedSurfRate) > maximumSurfaceRate_)
                alpha = std::abs(maximumSurfaceRate_/weightedSurfRate);
        }
        else if (controlMode_ == VolumetricReservoirRate) {
            Scalar weightedResvRate = computeWeightedRate_(unconstraintReservoirRates_);
            if (std::abs(weightedResvRate) > maximumReservoirRate_)
                alpha = std::abs(maximumReservoirRate_/weightedResvRate);
        }
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            currentSurfaceRates_[phaseIdx] = unconstraintSurfaceRates_[phaseIdx]*alpha;
        }
    }

    /*!
     * \brief Called by the simulator after each Newton-Raphson iteration.
     */
    void endIteration()
    { }

    /*!
     * \brief Called by the simulator after each time step.
     */
    void endTimeStep()
    {
        if (simulator_.gridView().comm().rank() == 0) {
            Scalar weightedLimitedSurfaceRate = computeWeightedRate_(currentSurfaceRates_);

            std::cout << "Well '" << name() << "':\n";
            std::cout << " Control mode: " << controlMode_ << "\n";
            std::cout << " Target BHP: " << targetBhp_ << "\n";
            std::cout << " Observed BHP: " << observedBhp_ << "\n";
            std::cout << " Unconstraint phase-specific surface rates:\n";
            std::cout << "   oil=" << unconstraintSurfaceRates_[oilPhaseIdx]
                      << " m^3/s (=" << 543439.65*unconstraintSurfaceRates_[oilPhaseIdx] << " STB/day)\n";
            std::cout << "   gas=" << unconstraintSurfaceRates_[gasPhaseIdx]
                      << " m^3/s (=" << 3051.1872*unconstraintSurfaceRates_[gasPhaseIdx] << " MCF/day)\n";
            std::cout << "   water=" << unconstraintSurfaceRates_[waterPhaseIdx]
                      << " m^3/s (=" << 543439.65*unconstraintSurfaceRates_[waterPhaseIdx] << " STB/day)\n";
            std::cout << " Rate-limited phase-specific surface rates:\n";
            std::cout << "   oil=" << currentSurfaceRates_[oilPhaseIdx]
                      << " m^3/s (=" << 543439.65*currentSurfaceRates_[oilPhaseIdx] << " STB/day)\n";
            std::cout << "   gas=" << currentSurfaceRates_[gasPhaseIdx]
                      << " m^3/s (=" << 3051.1872*currentSurfaceRates_[gasPhaseIdx] << " MCF/day)\n";
            std::cout << "   water=" << currentSurfaceRates_[waterPhaseIdx]
                      << " m^3/s (=" << 543439.65*currentSurfaceRates_[waterPhaseIdx] << " STB/day)\n";
            std::cout << " Rate-limited weighted limited rate: " << weightedLimitedSurfaceRate << "\n";
            std::cout << " Maximum weighted surface rate: " << maximumSurfaceRate_ << "\n";
        }
    }

    /*!
     * \brief Computes the source term for a degree of freedom.
     */
    template <class Context>
    void computeTotalRatesForDof(RateVector &q,
                                 const Context &context,
                                 int dofIdx,
                                 int timeIdx) const
    {
        q = 0.0;

        int globalDofIdx = context.globalSpaceIndex(dofIdx, timeIdx);
        auto dofVarsIt = dofVariables_.find(globalDofIdx);
        if (!isOpen_ || dofVarsIt == dofVariables_.end())
            return;

        std::array<Scalar, numPhases> volumetricRates;
        computeUnconstraintVolumetricDofRates_(volumetricRates,
                                               dofVarsIt->second,
                                               context,
                                               dofIdx,
                                               timeIdx);

        limitVolumetricReservoirRates_(volumetricRates,
                                       dofVarsIt->second,
                                       context,
                                       dofIdx,
                                       timeIdx);

        // convert to mass rates
        RateVector phaseRate;
        const auto &volVars = context.intensiveQuantities(dofIdx, timeIdx);
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            phaseRate.setVolumetricRate(volVars.fluidState(), phaseIdx, volumetricRates[phaseIdx]);
            q += phaseRate;
        }
        Valgrind::CheckDefined(q);
    }

    /*!
     * \brief This method writes the complete state of the well
     *        to the harddisk.
     */
    template <class Restarter>
    void serialize(Restarter &res)
    {
        res.serializeSectionBegin("PeacemanWell");

        res.serializeStream()
            << targetThp_ << " "
            << targetBhp_ << " "
            << controlMode_ << " "
            << wellType_ << " "
            << maximumSurfaceRate_ << " "
            << maximumReservoirRate_ << " "
            << isOpen_ << " "
            << injectedPhaseIdx_ << " ";

        // fluid state
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            res.serializeStream()
                << volumetricWeight_[phaseIdx] << " ";
        res.serializeSectionEnd();
    }

    /*!
     * \brief This method restores the complete state of the well
     *        from disk.
     *
     * It is the inverse of the serialize() method.
     *
     * \tparam Restarter The deserializer type
     *
     * \param res The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter &res)
    {
        res.deserializeSectionBegin("PeacemanWell");
        res.deserializeStream()
            >> targetThp_
            >> targetBhp_
            >> controlMode_
            >> wellType_
            >> maximumSurfaceRate_
            >> maximumReservoirRate_
            >> isOpen_
            >> injectedPhaseIdx_;

        // fluid state
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx)
            res.serializeStream()
                >> volumetricWeight_[phaseIdx];
        res.deserializeSectionEnd();
    }

protected:
    template <class Context>
    void computeUnconstraintVolumetricDofRates_(std::array<Scalar, numPhases> &volRates,
                                                const DofVariables &dofVars,
                                                const Context &context,
                                                int dofIdx,
                                                int timeIdx) const
    {
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            volRates[phaseIdx] = 0.0;

        const auto &intQuants = context.intensiveQuantities(dofIdx, timeIdx);

        RateVector phaseRate;

        Scalar depth = context.pos(dofIdx, timeIdx)[2];

        // connection transmissibility factor for the current DOF.
        //Scalar connTransFac = dofVars.connectionTransmissibilityFactor;

        // Intrinsic permeability. E100 uses the geometric average of the X and the Y
        // permability as the effective one...
        const auto &K = intQuants.intrinsicPermeability();
        Scalar Kvertical = Opm::utils::geometricAverage(K[0][0], K[1][1]);

        // calculate the equivalence radius of the well inside the cell. This seems to be
        // E100 vodoo...
        Scalar Dx = dofVars.effectiveSize[0];
        Scalar Dy = dofVars.effectiveSize[1];
        Scalar Dz = dofVars.effectiveSize[2];

        Scalar tmp = std::sqrt(K[1][1]/K[0][0]);
        Scalar tmp2 = std::sqrt(tmp);
        Scalar rEquiv = 0.28*std::sqrt(Dx*Dx*tmp + Dy*Dy/tmp)/(tmp2 + 1/tmp2);

        // the well borehole radius for the cell
        Scalar rWell = dofVars.boreholeRadius;

        // the skin factor of the well at the cell
        Scalar skinFactor = dofVars.skinFactor;

        // gravity constant
        Scalar g = 9.81;

        typename FluidSystem::ParameterCache paramCache;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // well model due to Peaceman; see Chen et al., p. 449

            // bottom hole pressure
            Scalar pbh = effectiveBottomHolePressure_;

            // phase pressure in grid cell
            Scalar p = intQuants.fluidState().pressure(phaseIdx);

            // density of fluid phase
            Scalar rho;

            Scalar lambda;
            if (wellType_ == Producer) {
                //assert(p < pbh);
                rho = intQuants.fluidState().density(phaseIdx);
                lambda = intQuants.mobility(phaseIdx);
            }
            else if (wellType_ == Injector) {
                //assert(p > pbh);

                if (phaseIdx != injectedPhaseIdx_)
                    continue;

                injectionFluidState_.setPressure(phaseIdx, p);

                typename FluidSystem::ParameterCache paramCache;
                paramCache.updateAll(injectionFluidState_);

                rho = FluidSystem::density(injectionFluidState_, paramCache, phaseIdx);
                lambda = 1.0/FluidSystem::viscosity(injectionFluidState_, paramCache, phaseIdx);
            }
            else
                OPM_THROW(std::logic_error,
                          "Type of well \"" << name() << "\" is undefined");

            Valgrind::CheckDefined(pbh);
            Valgrind::CheckDefined(p);
            Valgrind::CheckDefined(rho);
            Valgrind::CheckDefined(lambda);

            // pressure in the borehole ("hole pressure") at the given location
            Scalar ph = pbh + rho*g*(bottomDepth_ - depth);

            // volumetric flux of the phase from the well to the reservoir
            volRates[phaseIdx] =
                lambda*(ph - p)*Kvertical*Dz*2*M_PI
                / (std::log(rEquiv/rWell) + skinFactor);

            // make sure that injector wells only inject (-> positive rates) and
            // producers only produce (-> negative rates). TODO: this is not what happens
            // in the physical world, as cross-flow may occur...
            if (wellType_ == Injector)
                volRates[phaseIdx] = std::max(volRates[phaseIdx], 0.0);
            else if (wellType_ == Producer)
                volRates[phaseIdx] = std::min(volRates[phaseIdx], 0.0);

            Valgrind::CheckDefined(g);
            Valgrind::CheckDefined(ph);
            Valgrind::CheckDefined(volRates[phaseIdx]);
        }
    }

    /*!
     * \brief Given the volumetric rates for all phases, return the
     *        corresponding weighted rate
     *
     * The weights are user-specified and can be set using
     * setVolumetricPhaseWeights()
     */
    Scalar computeWeightedRate_(const std::array<Scalar, numPhases> &volRates) const
    {
        Scalar result = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
            result += volRates[phaseIdx]*volumetricWeight_[phaseIdx];
        return result;
    }

    /*!
     * \brief Convert volumetric reservoir rates into volumetric volume rates.
     *
     * This requires the density and composition of the phases and
     * thus the applicable fluid state.
     */
    template <class FluidState>
    void computeSurfaceRates_(std::array<Scalar, numPhases> &surfaceRates,
                              const std::array<Scalar, numPhases> &reservoirRate,
                              const FluidState &fluidState) const
    {
        // If your compiler bails out here, you have not
        // chosen the correct fluid system. Currently,
        // only Opm::FluidSystems::BlackOil is supported,
        // sorry...
        Scalar rhoOilSurface = FluidSystem::surfaceDensity(oilPhaseIdx, /*regionIdx=*/0);
        Scalar rhoGasSurface = FluidSystem::surfaceDensity(gasPhaseIdx, /*regionIdx=*/0);
        Scalar rhoWaterSurface = FluidSystem::surfaceDensity(waterPhaseIdx, /*regionIdx=*/0);

        // oil
        surfaceRates[oilPhaseIdx] =
            reservoirRate[oilPhaseIdx]
            * fluidState.density(oilPhaseIdx)
            * fluidState.massFraction(oilPhaseIdx, oilCompIdx)
            / rhoOilSurface;

        // gas
        surfaceRates[gasPhaseIdx] =
            // gas in gas phase
            reservoirRate[gasPhaseIdx]
            * fluidState.density(gasPhaseIdx)
            / rhoGasSurface
            +
            // gas in oil phase
            reservoirRate[oilPhaseIdx]
            * fluidState.density(oilPhaseIdx)
            * fluidState.massFraction(oilPhaseIdx, gasCompIdx)
            / rhoGasSurface;

        // water
        surfaceRates[waterPhaseIdx] =
            reservoirRate[waterPhaseIdx]
            * fluidState.density(waterPhaseIdx)
            / rhoWaterSurface;
    }

    /*!
     * \brief Calculate the final mass rate which ought to be used
     *        after the user specified rate limits have been applied.
     *
     * The input rates are the volumetric reservoir phase rates which
     * emerge if the user-defined rate limits are not considered.
     */
    template <class Context>
    void limitVolumetricReservoirRates_(std::array<Scalar, numPhases> &reservoirDofVolRates,
                                        const DofVariables &dofVars,
                                        const Context &context,
                                        int dofIdx,
                                        int timeIdx) const
    {
        // we don't look at the rates if we control for one of the pressures...
        if (controlMode_ == ControlMode::BottomHolePressure ||
            controlMode_ == ControlMode::TopHolePressure)
            return;

        //int globalDofIdx = context.globalSpaceIndex(dofIdx, timeIdx);
        if (controlMode_ == ControlMode::VolumetricSurfaceRate) {
            // convert the volumetric reservoir rates to to volumetric surface rates.
            Scalar weightedSurfaceRate;

            std::array<Scalar, numPhases> surfaceDofVolRates;
            computeSurfaceRates_(surfaceDofVolRates,
                                 reservoirDofVolRates,
                                 context.intensiveQuantities(dofIdx, timeIdx).fluidState());

            // subtract the effect of the unmodified degree of freedom to the total well
            // rate and add the effect of the potentially modified one. (i.e., add
            // the difference due to the modified primary variables at the DOF.)
            weightedSurfaceRate = computeWeightedRate_(unconstraintSurfaceRates_);
            weightedSurfaceRate -= computeWeightedRate_(dofVars.unconstraintSurfaceRates);
            weightedSurfaceRate += computeWeightedRate_(surfaceDofVolRates);

            // if we're below the limit, we're gold
            if (std::abs(weightedSurfaceRate) <= maximumSurfaceRate_)
                return;

            // if not, reduce the well's rate. so far, we reduce the reservoir rate
            // proportionally. that is slightly wrong, but we don't care...
            Scalar alpha = maximumSurfaceRate_ / std::abs(weightedSurfaceRate);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                reservoirDofVolRates[phaseIdx] *= alpha;
        }
        else if (controlMode_ == ControlMode::VolumetricReservoirRate) {
            // calulate the current total rate of the well: first subtract the rate of the
            // DOF from the prestine well rates, then add the just calculated rate to it.
            Scalar weightedReservoirRate = computeWeightedRate_(unconstraintReservoirRates_);
            weightedReservoirRate -= computeWeightedRate_(dofVars.unconstraintRates);
            weightedReservoirRate += computeWeightedRate_(reservoirDofVolRates);

            // if we're below the limit, we're gold
            if (std::abs(weightedReservoirRate) <= maximumReservoirRate_)
                return;

            // if not, we have to reduce the total rate. We do this proportionally to the
            // volume produced of each fluid.
            Scalar alpha = maximumReservoirRate_ / std::abs(weightedReservoirRate);
            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                reservoirDofVolRates[phaseIdx] *= alpha;
        }
    }

    const Simulator &simulator_;

    std::string name_;

    std::unordered_map<int, DofVariables> dofVariables_;

    // the sum of the total volumes of all the degrees of freedoms that interact with the well
    Scalar wellTotalVolume_;

    // The assumed bottom and top hole pressures as specified by the user
    Scalar targetBhp_;
    Scalar targetThp_;

    // real pressure seen at the bottom of the borehole
    Scalar observedBhp_;

    // The sum of the unconstraint volumetric reservoir rates of all
    // degrees of freedom in the well for all fluid phases. This is
    // calculated at the beginning of each iteration and used to
    // impose rate limits. (basically, this can be calculated from the
    // above structure but it would be quite slow because this number
    // is required for each DOF...)
    std::array<Scalar, numPhases> unconstraintReservoirRates_;

    // the same as the above but as surface rate
    std::array<Scalar, numPhases> unconstraintSurfaceRates_;

    // the total rate of the well with limits applied
    std::array<Scalar, numPhases> currentSurfaceRates_;

    // specifies the quantities which are controlled for (i.e., which
    // should be assumed to be externally specified and which should
    // be computed based on those)
    ControlMode controlMode_;

    // the type of the well (injector, producer or undefined)
    WellType wellType_;

    // The bottom hole pressure to be used by the well model. This may
    // be computed from the top hole pressure (if the control mode is
    // TopHolePressure), or it may be just the user-specified bottom
    // hole pressure if the control mode is BottomHolePressure.
    Scalar effectiveBottomHolePressure_;

    // The maximum weighted volumetric surface rates specified by the
    // user. This is used to apply rate limits and it is to be read as
    // the maximum absolute value of the rate, i.e., the well can
    // produce or inject the given amount.
    Scalar maximumSurfaceRate_;

    // The maximum weighted volumetric reservoir rates specified by
    // the user. This is used to apply rate limits and it is to be
    // read as the maximum absolute value of the rate, i.e., the well
    // can produce or inject the given amount.
    Scalar maximumReservoirRate_;

    // Specifies whether the well is currently shut or not. If true,
    // this has the same effect as setting the minimum and maximum
    // well rates to zero, but with this the well can be shut and
    // opened without remembering the well rates
    bool isOpen_;

    // The relative weight of the volumetric rate of each fluid
    Scalar volumetricWeight_[numPhases];

    // The thermodynamic state of the fluid which gets injected
    //
    // The fact that this attribute is mutable is kind of an hack
    // which can be avoided using a PressureOverlayFluidState, but
    // then performance would be slightly worse...
    mutable FluidState injectionFluidState_;

    int injectedPhaseIdx_;

    // the depth of the deepest DOF. (actually, the center of this
    // DOF, but the difference should be minimal.)
    Scalar bottomDepth_;

    // global index of the DOF at the bottom of the well
    int bottomDofGlobalIdx_;
};
} // namespace Ewoms

#endif
