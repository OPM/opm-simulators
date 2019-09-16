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
 * \brief Contains the classes required to extend the black-oil model by energy.
 */
#ifndef EWOMS_BLACK_OIL_ENERGY_MODULE_HH
#define EWOMS_BLACK_OIL_ENERGY_MODULE_HH

#include "blackoilproperties.hh"
#include <opm/models/io/vtkblackoilenergymodule.hh>
#include <opm/models/common/quantitycallbacks.hh>

#include <opm/material/common/Tabulated1DFunction.hpp>

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/Unused.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/fvector.hh>

#include <string>

namespace Opm {
/*!
 * \ingroup BlackOil
 * \brief Contains the high level supplements required to extend the black oil
 *        model by energy.
 */
template <class TypeTag, bool enableEnergyV = GET_PROP_VALUE(TypeTag, EnableEnergy)>
class BlackOilEnergyModule
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    static constexpr unsigned temperatureIdx = Indices::temperatureIdx;
    static constexpr unsigned contiEnergyEqIdx = Indices::contiEnergyEqIdx;

    static constexpr unsigned enableEnergy = enableEnergyV;
    static constexpr unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr unsigned numPhases = FluidSystem::numPhases;

public:
    /*!
     * \brief Register all run-time parameters for the black-oil energy module.
     */
    static void registerParameters()
    {
        if (!enableEnergy)
            // energys have been disabled at compile time
            return;

        Opm::VtkBlackOilEnergyModule<TypeTag>::registerParameters();
    }

    /*!
     * \brief Register all energy specific VTK and ECL output modules.
     */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enableEnergy)
            // energys have been disabled at compile time
            return;

        model.addOutputModule(new Opm::VtkBlackOilEnergyModule<TypeTag>(simulator));
    }

    static bool primaryVarApplies(unsigned pvIdx)
    {
        if (!enableEnergy)
            // energys have been disabled at compile time
            return false;

        return pvIdx == temperatureIdx;
    }

    static std::string primaryVarName(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        return "temperature";
    }

    static Scalar primaryVarWeight(unsigned pvIdx OPM_OPTIM_UNUSED)
    {
        assert(primaryVarApplies(pvIdx));

        // TODO: it may be beneficial to chose this differently.
        return static_cast<Scalar>(1.0);
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableEnergy)
            return false;

        return eqIdx == contiEnergyEqIdx;
    }

    static std::string eqName(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        return "conti^energy";
    }

    static Scalar eqWeight(unsigned eqIdx OPM_OPTIM_UNUSED)
    {
        assert(eqApplies(eqIdx));

        return 1.0;
    }

    // must be called after water storage is computed
    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enableEnergy)
            return;

        const auto& poro = Opm::decay<LhsEval>(intQuants.porosity());

        // accumulate the internal energy of the fluids
        const auto& fs = intQuants.fluidState();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            const auto& u = Opm::decay<LhsEval>(fs.internalEnergy(phaseIdx));
            const auto& S = Opm::decay<LhsEval>(fs.saturation(phaseIdx));
            const auto& rho = Opm::decay<LhsEval>(fs.density(phaseIdx));

            storage[contiEnergyEqIdx] += poro*S*u*rho;
        }

        // add the internal energy of the rock
        Scalar refPoro = intQuants.referencePorosity();
        const auto& uRock = Opm::decay<LhsEval>(intQuants.rockInternalEnergy());
        storage[contiEnergyEqIdx] += (1.0 - refPoro)*uRock;
        storage[contiEnergyEqIdx] *= GET_PROP_VALUE(TypeTag, BlackOilEnergyScalingFactor);
    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)
    {
        if (!enableEnergy)
            return;

        flux[contiEnergyEqIdx] = 0.0;

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned focusIdx = elemCtx.focusDofIndex();
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            unsigned upIdx = extQuants.upstreamIndex(phaseIdx);
            if (upIdx == focusIdx)
                addPhaseEnthalpyFlux_<Evaluation>(flux, phaseIdx, elemCtx, scvfIdx, timeIdx);
            else
                addPhaseEnthalpyFlux_<Scalar>(flux, phaseIdx, elemCtx, scvfIdx, timeIdx);
        }

        // diffusive energy flux
        flux[contiEnergyEqIdx] += extQuants.energyFlux();
        flux[contiEnergyEqIdx] *= GET_PROP_VALUE(TypeTag, BlackOilEnergyScalingFactor);
    }

    template <class UpstreamEval>
    static void addPhaseEnthalpyFlux_(RateVector& flux,
                                      unsigned phaseIdx,
                                      const ElementContext& elemCtx,
                                      unsigned scvfIdx,
                                      unsigned timeIdx)
    {
        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        unsigned upIdx = extQuants.upstreamIndex(phaseIdx);
        const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& fs = up.fluidState();

        const auto& volFlux = extQuants.volumeFlux(phaseIdx);
        flux[contiEnergyEqIdx] +=
            Opm::decay<UpstreamEval>(fs.enthalpy(phaseIdx))
            * Opm::decay<UpstreamEval>(fs.density(phaseIdx))
            * volFlux;
    }

    static void addToEnthalpyRate(RateVector& flux,
                                  const Evaluation& hRate)
    {
        if (!enableEnergy)
            return;

        flux[contiEnergyEqIdx] += hRate;
    }

    /*!
     * \brief Assign the energy specific primary variables to a PrimaryVariables object
     */
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  Scalar temperature OPM_UNUSED)
    {
        if (!enableEnergy)
            return;

        priVars[temperatureIdx] = temperatureIdx;
    }

    /*!
     * \brief Assign the energy specific primary variables to a PrimaryVariables object
     */
    template <class FluidState>
    static void assignPrimaryVars(PrimaryVariables& priVars,
                                  const FluidState& fluidState)
    {
        if (!enableEnergy)
            return;

        priVars[temperatureIdx] = fluidState.temperature(/*phaseIdx=*/0);
    }

    /*!
     * \brief Do a Newton-Raphson update the primary variables of the energys.
     */
    static void updatePrimaryVars(PrimaryVariables& newPv,
                                  const PrimaryVariables& oldPv,
                                  const EqVector& delta)
    {
        if (!enableEnergy)
            return;

        // do a plain unchopped Newton update
        newPv[temperatureIdx] = oldPv[temperatureIdx] - delta[temperatureIdx];
    }

    /*!
     * \brief Return how much a Newton-Raphson update is considered an error
     */
    static Scalar computeUpdateError(const PrimaryVariables& oldPv OPM_UNUSED,
                                     const EqVector& delta OPM_UNUSED)
    {
        // do not consider consider the cange of energy primary variables for
        // convergence
        // TODO: maybe this should be changed
        return static_cast<Scalar>(0.0);
    }

    /*!
     * \brief Return how much a residual is considered an error
     */
    static Scalar computeResidualError(const EqVector& resid)
    {
        // do not weight the residual of energy when it comes to convergence
        return std::abs(Opm::scalarValue(resid[contiEnergyEqIdx]));
    }

    template <class DofEntity>
    static void serializeEntity(const Model& model, std::ostream& outstream, const DofEntity& dof)
    {
        if (!enableEnergy)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        const PrimaryVariables& priVars = model.solution(/*timeIdx=*/0)[dofIdx];
        outstream << priVars[temperatureIdx];
    }

    template <class DofEntity>
    static void deserializeEntity(Model& model, std::istream& instream, const DofEntity& dof)
    {
        if (!enableEnergy)
            return;

        unsigned dofIdx = model.dofMapper().index(dof);
        PrimaryVariables& priVars0 = model.solution(/*timeIdx=*/0)[dofIdx];
        PrimaryVariables& priVars1 = model.solution(/*timeIdx=*/1)[dofIdx];

        instream >> priVars0[temperatureIdx];

        // set the primary variables for the beginning of the current time step.
        priVars1 = priVars0[temperatureIdx];
    }
};

/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilEnergyIntensiveQuantities
 *
 * \brief Provides the volumetric quantities required for the equations needed by the
 *        energys extension of the black-oil model.
 */
template <class TypeTag, bool enableEnergyV = GET_PROP_VALUE(TypeTag, EnableEnergy)>
class BlackOilEnergyIntensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, SolidEnergyLaw) SolidEnergyLaw;
    typedef typename GET_PROP_TYPE(TypeTag, ThermalConductionLaw) ThermalConductionLaw;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef BlackOilEnergyModule<TypeTag> EnergyModule;

    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    static constexpr int temperatureIdx = Indices::temperatureIdx;
    static constexpr int waterPhaseIdx = FluidSystem::waterPhaseIdx;


public:
    /*!
     * \brief Update the temperature of the intensive quantity's fluid state
     *
     */
    void updateTemperature_(const ElementContext& elemCtx,
                            unsigned dofIdx,
                            unsigned timeIdx)
    {
        auto& fs = asImp_().fluidState_;
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        // set temperature
        fs.setTemperature(priVars.makeEvaluation(temperatureIdx, timeIdx));
    }

    /*!
     * \brief Compute the intensive quantities needed to handle energy conservation
     *
     */
    void updateEnergyQuantities_(const ElementContext& elemCtx,
                                 unsigned dofIdx,
                                 unsigned timeIdx,
                                 const typename FluidSystem::template ParameterCache<Evaluation>& paramCache)
    {
        auto& fs = asImp_().fluidState_;

        // compute the specific enthalpy of the fluids, the specific enthalpy of the rock
        // and the thermal condictivity coefficients
        for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                fs.setEnthalpy(phaseIdx, 0.0);
                continue;
            }

            const auto& h = FluidSystem::enthalpy(fs, paramCache, phaseIdx);
            fs.setEnthalpy(phaseIdx, h);
        }

        const auto& solidEnergyLawParams = elemCtx.problem().solidEnergyLawParams(elemCtx, dofIdx, timeIdx);
        rockInternalEnergy_ = SolidEnergyLaw::solidInternalEnergy(solidEnergyLawParams, fs);

        const auto& thermalConductionLawParams = elemCtx.problem().thermalConductionLawParams(elemCtx, dofIdx, timeIdx);
        totalThermalConductivity_ = ThermalConductionLaw::thermalConductivity(thermalConductionLawParams, fs);
    }

    const Evaluation& rockInternalEnergy() const
    { return rockInternalEnergy_; }

    const Evaluation& totalThermalConductivity() const
    { return totalThermalConductivity_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation rockInternalEnergy_;
    Evaluation totalThermalConductivity_;
};

template <class TypeTag>
class BlackOilEnergyIntensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    static constexpr bool enableTemperature = GET_PROP_VALUE(TypeTag, EnableTemperature);

public:
    void updateTemperature_(const ElementContext& elemCtx,
                            unsigned dofIdx,
                            unsigned timeIdx)
    {
        if (enableTemperature) {
            // even if energy is conserved, the temperature can vary over the spatial
            // domain if the EnableTemperature property is set to true
            auto& fs = asImp_().fluidState_;
            Scalar T = elemCtx.problem().temperature(elemCtx, dofIdx, timeIdx);
            fs.setTemperature(T);
        }
    }

    void updateEnergyQuantities_(const ElementContext& elemCtx OPM_UNUSED,
                                 unsigned dofIdx OPM_UNUSED,
                                 unsigned timeIdx OPM_UNUSED,
                                 const typename FluidSystem::template ParameterCache<Evaluation>& paramCache OPM_UNUSED)
    { }

    const Evaluation& rockInternalEnergy() const
    { throw std::logic_error("Requested the rock internal energy, which is "
                             "unavailable because energy is not conserved"); }

    const Evaluation& totalThermalConductivity() const
    { throw std::logic_error("Requested the total thermal conductivity, which is "
                             "unavailable because energy is not conserved"); }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }
};


/*!
 * \ingroup BlackOil
 * \class Opm::BlackOilEnergyExtensiveQuantities
 *
 * \brief Provides the energy specific extensive quantities to the generic black-oil
 *        module's extensive quantities.
 */
template <class TypeTag, bool enableEnergyV = GET_PROP_VALUE(TypeTag, EnableEnergy)>
class BlackOilEnergyExtensiveQuantities
{
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, ExtensiveQuantities) ExtensiveQuantities;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Opm::MathToolbox<Evaluation> Toolbox;

    typedef BlackOilEnergyModule<TypeTag> EnergyModule;

    static const int dimWorld = GridView::dimensionworld;
    typedef Dune::FieldVector<Scalar, dimWorld> DimVector;
    typedef Dune::FieldVector<Evaluation, dimWorld> DimEvalVector;

public:
    void updateEnergy(const ElementContext& elemCtx,
                      unsigned scvfIdx,
                      unsigned timeIdx)
    {
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        Scalar faceArea = scvf.area();
        unsigned inIdx = scvf.interiorIndex();
        unsigned exIdx = scvf.exteriorIndex();
        const auto& inIq = elemCtx.intensiveQuantities(inIdx, timeIdx);
        const auto& exIq = elemCtx.intensiveQuantities(exIdx, timeIdx);
        const auto& inFs = inIq.fluidState();
        const auto& exFs = exIq.fluidState();

        Evaluation deltaT;
        if (elemCtx.focusDofIndex() == inIdx)
            deltaT =
                Opm::decay<Scalar>(exFs.temperature(/*phaseIdx=*/0))
                - inFs.temperature(/*phaseIdx=*/0);
        else if (elemCtx.focusDofIndex() == exIdx)
            deltaT =
                exFs.temperature(/*phaseIdx=*/0)
                - Opm::decay<Scalar>(inFs.temperature(/*phaseIdx=*/0));
        else
            deltaT =
                Opm::decay<Scalar>(exFs.temperature(/*phaseIdx=*/0))
                - Opm::decay<Scalar>(inFs.temperature(/*phaseIdx=*/0));

        Evaluation inLambda;
        if (elemCtx.focusDofIndex() == inIdx)
            inLambda = inIq.totalThermalConductivity();
        else
            inLambda = Opm::decay<Scalar>(inIq.totalThermalConductivity());

        Evaluation exLambda;
        if (elemCtx.focusDofIndex() == exIdx)
            exLambda = exIq.totalThermalConductivity();
        else
            exLambda = Opm::decay<Scalar>(exIq.totalThermalConductivity());

        auto distVec = elemCtx.pos(exIdx, timeIdx);
        distVec -= elemCtx.pos(inIdx, timeIdx);

        Evaluation H;
        if (inLambda > 0.0 && exLambda > 0.0) {
            // compute the "thermal transmissibility". In contrast to the normal
            // transmissibility this cannot be done as a preprocessing step because the
            // average thermal thermal conductivity is analogous to the permeability but
            // depends on the solution.
            Scalar alpha = elemCtx.problem().thermalHalfTransmissibility(elemCtx, scvfIdx, timeIdx);
            const Evaluation& inH = inLambda*alpha;
            const Evaluation& exH = exLambda*alpha;
            H = 1.0/(1.0/inH + 1.0/exH);
        }
        else
            H = 0.0;

        energyFlux_ = deltaT * (-H/faceArea);
    }

    template <class Context, class BoundaryFluidState>
    void updateEnergyBoundary(const Context& ctx,
                              unsigned scvfIdx,
                              unsigned timeIdx,
                              const BoundaryFluidState& boundaryFs)
    {
        const auto& stencil = ctx.stencil(timeIdx);
        const auto& scvf = stencil.boundaryFace(scvfIdx);

        unsigned inIdx = scvf.interiorIndex();
        const auto& inIq = ctx.intensiveQuantities(inIdx, timeIdx);
        const auto& inFs = inIq.fluidState();

        Evaluation deltaT;
        if (ctx.focusDofIndex() == inIdx)
            deltaT =
                boundaryFs.temperature(/*phaseIdx=*/0)
                - inFs.temperature(/*phaseIdx=*/0);
        else
            deltaT =
                Opm::decay<Scalar>(boundaryFs.temperature(/*phaseIdx=*/0))
                - Opm::decay<Scalar>(inFs.temperature(/*phaseIdx=*/0));

        Evaluation lambda;
        if (ctx.focusDofIndex() == inIdx)
            lambda = inIq.totalThermalConductivity();
        else
            lambda = Opm::decay<Scalar>(inIq.totalThermalConductivity());

        auto distVec = scvf.integrationPos();
        distVec -= ctx.pos(inIdx, timeIdx);

        if (lambda > 0.0) {
            // compute the "thermal transmissibility". In contrast to the normal
            // transmissibility this cannot be done as a preprocessing step because the
            // average thermal conductivity is analogous to the permeability but depends
            // on the solution.
            Scalar alpha = ctx.problem().thermalHalfTransmissibilityBoundary(ctx, scvfIdx);
            energyFlux_ = deltaT*lambda*(-alpha);
        }
        else
            energyFlux_ = 0.0;
    }

    const Evaluation& energyFlux()  const
    { return energyFlux_; }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation energyFlux_;
};

template <class TypeTag>
class BlackOilEnergyExtensiveQuantities<TypeTag, false>
{
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;

public:
    void updateEnergy(const ElementContext& elemCtx OPM_UNUSED,
                      unsigned scvfIdx OPM_UNUSED,
                      unsigned timeIdx OPM_UNUSED)
    {}

    template <class Context, class BoundaryFluidState>
    void updateEnergyBoundary(const Context& ctx OPM_UNUSED,
                              unsigned scvfIdx OPM_UNUSED,
                              unsigned timeIdx OPM_UNUSED,
                              const BoundaryFluidState& boundaryFs OPM_UNUSED)
    {}

    const Evaluation& energyFlux()  const
    { throw std::logic_error("Requested the energy flux, but energy is not conserved"); }
};


} // namespace Opm

#endif
