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
 * \copydoc Opm::BlackOilIntensiveQuantities
 */
#ifndef EWOMS_ECL_BLACK_OIL_INTENSIVE_QUANTITIES_HH
#define EWOMS_ECL_BLACK_OIL_INTENSIVE_QUANTITIES_HH
#include <opm/models/blackoil/blackoilintensivequantities.hh>
//#include "eclblackoilprimaryvariabels.hh"
namespace Opm {
/*!
 * \ingroup BlackOilModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the black-oil model.
 */
template <class TypeTag>
class EclBlackOilIntensiveQuantities 
    : public GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities)
    , public GET_PROP_TYPE(TypeTag, FluxModule)::FluxIntensiveQuantities
    , public BlackOilSolventIntensiveQuantities<TypeTag>
    , public BlackOilPolymerIntensiveQuantities<TypeTag>
    , public BlackOilFoamIntensiveQuantities<TypeTag>
    , public BlackOilBrineIntensiveQuantities<TypeTag>
    , public BlackOilEnergyIntensiveQuantities<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, DiscIntensiveQuantities) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluxModule) FluxModule;

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };
    enum { enablePolymer = GET_PROP_VALUE(TypeTag, EnablePolymer) };
    enum { enableFoam = GET_PROP_VALUE(TypeTag, EnableFoam) };
    enum { enableBrine = GET_PROP_VALUE(TypeTag, EnableBrine) };
    enum { enableTemperature = GET_PROP_VALUE(TypeTag, EnableTemperature) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents = GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { dimWorld = GridView::dimensionworld };
    enum { compositionSwitchIdx = Indices::compositionSwitchIdx };

    static const bool compositionSwitchEnabled = Indices::gasEnabled;
    static const bool waterEnabled = Indices::waterEnabled;

    typedef Opm::MathToolbox<Evaluation> Toolbox;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef typename FluxModule::FluxIntensiveQuantities FluxIntensiveQuantities;
    typedef Opm::BlackOilFluidState<Evaluation, FluidSystem, enableTemperature, enableEnergy, compositionSwitchEnabled,  enableBrine, Indices::numPhases > FluidState;

public:
    EclBlackOilIntensiveQuantities()
    {
        if (compositionSwitchEnabled) {
            fluidState_.setRs(0.0);
            fluidState_.setRv(0.0);
        }
    }

    EclBlackOilIntensiveQuantities(const EclBlackOilIntensiveQuantities& other) = default;

    EclBlackOilIntensiveQuantities& operator=(const EclBlackOilIntensiveQuantities& other) = default;

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);
        const auto& problem = elemCtx.problem();
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        
        const auto linearizationType = elemCtx.linearizationType();
        if(not(linearizationType.type == LinearizationType::seqtransport)){
            assert(priVars.simulationType() == PrimaryVariables::SimulationType::Implicit);
        }else{
            assert(priVars.simulationType() == PrimaryVariables::SimulationType::Seq);
        }
        

        asImp_().updateTemperature_(elemCtx, dofIdx, timeIdx);

        unsigned globalSpaceIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
        unsigned pvtRegionIdx = priVars.pvtRegionIndex();
        fluidState_.setPvtRegionIndex(pvtRegionIdx);

        asImp_().updateSaltConcentration_(elemCtx, dofIdx, timeIdx);

        // extract the water and the gas saturations for convenience
        Evaluation Sw = 0.0;
        if (waterEnabled) {
            if (priVars.primaryVarsMeaning() == PrimaryVariables::OnePhase_p) {
                Sw = 1.0;
            } else {
                Sw = priVars.makeEvaluation(Indices::waterSaturationIdx, timeIdx, linearizationType);
            }
        }
        Evaluation Sg = 0.0;
        if (compositionSwitchEnabled)
        {
            if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
                // -> threephase case
                assert( priVars.primaryVarsMeaning() != PrimaryVariables::OnePhase_p );
                Sg = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            } else if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv) {
                // -> gas-water case
                Sg = 1.0 - Sw;

                // deal with solvent
                if (enableSolvent)
                    Sg -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx, linearizationType);
            }
            else
            {
                assert(priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Rs);
                // -> oil-water case
                Sg = 0.0;
            }
        }

        Opm::Valgrind::CheckDefined(Sg);
        Opm::Valgrind::CheckDefined(Sw);

        Evaluation So = 1.0 - Sw - Sg;

        // deal with solvent
        if (enableSolvent)
            So -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx, linearizationType);

        if (FluidSystem::phaseIsActive(waterPhaseIdx)){
            fluidState_.setSaturation(waterPhaseIdx, Sw);
        }
        if (FluidSystem::phaseIsActive(gasPhaseIdx)){
            fluidState_.setSaturation(gasPhaseIdx, Sg);
        }
        if (FluidSystem::phaseIsActive(oilPhaseIdx)){
            fluidState_.setSaturation(oilPhaseIdx, So);
        }
        asImp_().solventPreSatFuncUpdate_(elemCtx, dofIdx, timeIdx);

        // now we compute all phase pressures
        Evaluation pC[numPhases];
        const auto& materialParams = problem.materialLawParams(elemCtx, dofIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);

        // NB which index is correct!!
        if (FluidSystem::phaseIsActive(waterPhaseIdx)){
            fluidState_.setPc(waterPhaseIdx, pC[waterPhaseIdx]);
        }
        if (FluidSystem::phaseIsActive(gasPhaseIdx)){
            fluidState_.setPc(gasPhaseIdx, pC[gasPhaseIdx]);
        }
        if (FluidSystem::phaseIsActive(oilPhaseIdx)){
            fluidState_.setPc(oilPhaseIdx, pC[oilPhaseIdx]);
        }

        // oil is the reference phase for pressure
        Evaluation totalSaturation = 1.0;
        if (linearizationType.type == Opm::LinearizationType::seqtransport) {
            totalSaturation = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx, linearizationType);
            //assert(totalSaturation<2);// debug to be sure no pressure is pressent
            //assert(timeIdx==0 ||  totalSaturation == 1);// for pure sequential with multiple steps this is always true.
            if(timeIdx == 1){
                totalSaturation = elemCtx.problem().totalSaturation(globalSpaceIdx);
            }
            
        }else if (linearizationType.type == Opm::LinearizationType::pressure) {
            // we are doing sequation pressure where totalSaturation may not be 1 any more
            if(timeIdx == 1){//pressure solves try to make new value i.e. timeIdx==0 so that totalSaturation==1;
                totalSaturation = elemCtx.problem().totalSaturation(globalSpaceIdx);
            }
        }
        fluidState_.setTotalSaturation(totalSaturation);
        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv) {
            Evaluation pg;
            if (linearizationType.type == Opm::LinearizationType::seqtransport) {
                pg = Toolbox::value(
                    elemCtx.problem().pressure(globalSpaceIdx,
                                                  gasPhaseIdx,
                                                  timeIdx)); // get current value assuem this is not updated
            } else {
                pg = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx, linearizationType);
                //assert(pg>2.0);//enure not totalSaturation is present
            }
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                if (FluidSystem::phaseIsActive(phaseIdx))
                    fluidState_.setPressure(phaseIdx, pg + (pC[phaseIdx] - pC[gasPhaseIdx]));
        } else {
            Evaluation po;
            if (linearizationType.type == Opm::LinearizationType::seqtransport) {
                // get current value assume this is never updated
                po = Toolbox::value(elemCtx.problem().pressure(globalSpaceIdx, oilPhaseIdx, timeIdx));
            } else {
                po = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx, linearizationType);
                //assert(po>2.0);//enure not totalSaturation is present
            }
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                if (FluidSystem::phaseIsActive(phaseIdx))
                    fluidState_.setPressure(phaseIdx, po + (pC[phaseIdx] - pC[oilPhaseIdx]));
        }

        // calculate relative permeabilities. note that we store the result into the
        // mobility_ class attribute. the division by the phase viscosity happens later.
        MaterialLaw::relativePermeabilities(mobility_, materialParams, fluidState_);
        Opm::Valgrind::CheckDefined(mobility_);

        // update the Saturation functions for the blackoil solvent module.
        asImp_().solventPostSatFuncUpdate_(elemCtx, dofIdx, timeIdx);

        Evaluation SoMax = 0.0;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            SoMax = Opm::max(fluidState_.saturation(oilPhaseIdx),
                             elemCtx.problem().maxOilSaturation(globalSpaceIdx));
        }

        // take the meaning of the switiching primary variable into account for the gas
        // and oil phase compositions
        if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
            // in the threephase case, gas and oil phases are potentially present, i.e.,
            // we use the compositions of the gas-saturated oil and oil-saturated gas.
            if (FluidSystem::enableDissolvedGas()) {
                Scalar RsMax = elemCtx.problem().maxGasDissolutionFactor(timeIdx, globalSpaceIdx);
                const Evaluation& RsSat =
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            oilPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);
                fluidState_.setRs(Opm::min(RsMax, RsSat));
            }
            else if (compositionSwitchEnabled)
                fluidState_.setRs(0.0);

            if (FluidSystem::enableVaporizedOil()) {
                Scalar RvMax = elemCtx.problem().maxOilVaporizationFactor(timeIdx, globalSpaceIdx);
                const Evaluation& RvSat =
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);
                fluidState_.setRv(Opm::min(RvMax, RvSat));
            }
            else if (compositionSwitchEnabled)
                fluidState_.setRv(0.0);
        }
        else if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_po_Rs) {
            // if the switching variable is the mole fraction of the gas component in the
            Scalar RsMax = elemCtx.problem().maxGasDissolutionFactor(timeIdx, globalSpaceIdx);
            if(linearizationType.type == Opm::LinearizationType::pressure){
                const Evaluation& RsSat =
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            oilPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);
                RsMax = Opm::min(RsMax,Toolbox::value(RsSat));
                auto priVarsCopy = priVars;
                priVarsCopy[Indices::compositionSwitchIdx] = Opm::min(RsMax,priVarsCopy[Indices::compositionSwitchIdx]);
                const auto& Rs = priVarsCopy.makeEvaluation(Indices::compositionSwitchIdx, timeIdx, linearizationType);                
                fluidState_.setRs(Rs);
            }else{
                // oil phase, we can directly set the composition of the oil phase
                const auto& Rs = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx, linearizationType);
                fluidState_.setRs(Opm::min(RsMax, Rs));
            }    
            if (FluidSystem::enableVaporizedOil()) {
                // the gas phase is not present, but we need to compute its "composition"
                // for the gravity correction anyway
                Scalar RvMax = elemCtx.problem().maxOilVaporizationFactor(timeIdx, globalSpaceIdx);
                const auto& RvSat =
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);

                fluidState_.setRv(Opm::min(RvMax, RvSat));
            }
            else
                fluidState_.setRv(0.0);
        }
        else if (priVars.primaryVarsMeaning() == PrimaryVariables::Sw_pg_Rv) {
            
            if(linearizationType.type == Opm::LinearizationType::pressure){
                const auto& RvSat =
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);
                auto priVarsCopy = priVars;
                priVarsCopy[Indices::compositionSwitchIdx] = Opm::min(Toolbox::value(RvSat),priVarsCopy[Indices::compositionSwitchIdx]);
                const auto& Rv = priVarsCopy.makeEvaluation(Indices::compositionSwitchIdx, timeIdx, linearizationType);
                fluidState_.setRv(Rv);
            }else{
                const auto& Rv = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx, linearizationType);
                fluidState_.setRv(Rv);
            }

            if (FluidSystem::enableDissolvedGas()) {
                // the oil phase is not present, but we need to compute its "composition" for
                // the gravity correction anyway
                Scalar RsMax = elemCtx.problem().maxGasDissolutionFactor(timeIdx, globalSpaceIdx);
                const auto& RsSat =
                    FluidSystem::saturatedDissolutionFactor(fluidState_,
                                                            oilPhaseIdx,
                                                            pvtRegionIdx,
                                                            SoMax);

                fluidState_.setRs(Opm::min(RsMax, RsSat));
            } else {
                fluidState_.setRs(0.0);
            }
        } else {
            assert(priVars.primaryVarsMeaning() == PrimaryVariables::OnePhase_p);
        }

        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.setRegionIndex(pvtRegionIdx);
        if(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)){
            paramCache.setMaxOilSat(SoMax);
        }
        paramCache.updateAll(fluidState_);

        // compute the phase densities and transform the phase permeabilities into mobilities
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            const auto& b = FluidSystem::inverseFormationVolumeFactor(fluidState_, phaseIdx, pvtRegionIdx);
            fluidState_.setInvB(phaseIdx, b);

            const auto& mu = FluidSystem::viscosity(fluidState_, paramCache, phaseIdx);
            mobility_[phaseIdx] /= mu;
        }
        Opm::Valgrind::CheckDefined(mobility_);

        // calculate the phase densities
        Evaluation rho;
        if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
            rho = fluidState_.invB(waterPhaseIdx);
            rho *= FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
            fluidState_.setDensity(waterPhaseIdx, rho);
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            rho = fluidState_.invB(gasPhaseIdx);
            rho *= FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableVaporizedOil()) {
                rho +=
                    fluidState_.invB(gasPhaseIdx) *
                    fluidState_.Rv() *
                    FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            }
            fluidState_.setDensity(gasPhaseIdx, rho);
        }

        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            rho = fluidState_.invB(oilPhaseIdx);
            rho *= FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableDissolvedGas()) {
                rho +=
                    fluidState_.invB(oilPhaseIdx) *
                    fluidState_.Rs() *
                    FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            }
            fluidState_.setDensity(oilPhaseIdx, rho);
        }

        // retrieve the porosity from the problem
        referencePorosity_ = problem.porosity(elemCtx, dofIdx, timeIdx);
        porosity_ = referencePorosity_;

        // the porosity must be modified by the compressibility of the
        // rock...
        Scalar rockCompressibility = problem.rockCompressibility(elemCtx, dofIdx, timeIdx);
        if (rockCompressibility > 0.0) {
            Scalar rockRefPressure = problem.rockReferencePressure(elemCtx, dofIdx, timeIdx);
            Evaluation x;
            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                x = rockCompressibility*(fluidState_.pressure(oilPhaseIdx) - rockRefPressure);
            } else if (FluidSystem::phaseIsActive(waterPhaseIdx)){
            } else {
                x = rockCompressibility*(fluidState_.pressure(gasPhaseIdx) - rockRefPressure);
            }
            porosity_ *= 1.0 + x + 0.5*x*x;
        }

        // deal with water induced rock compaction
        porosity_ *= problem.template rockCompPoroMultiplier<Evaluation>(*this, globalSpaceIdx);

        // this line is doing things for sequential simulation to avoid modifying storage term.
        porosity_ *= totalSaturation;

        asImp_().solventPvtUpdate_(elemCtx, dofIdx, timeIdx);
        asImp_().polymerPropertiesUpdate_(elemCtx, dofIdx, timeIdx);
        asImp_().updateEnergyQuantities_(elemCtx, dofIdx, timeIdx, paramCache);
        asImp_().foamPropertiesUpdate_(elemCtx, dofIdx, timeIdx);

        // update the quantities which are required by the chosen
        // velocity model
        // TODO: use the focus time here
        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

#ifndef NDEBUG
        // some safety checks in debug mode
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            assert(Opm::isfinite(fluidState_.density(phaseIdx)));
            assert(Opm::isfinite(fluidState_.saturation(phaseIdx)));
            assert(Opm::isfinite(fluidState_.temperature(phaseIdx)));
            assert(Opm::isfinite(fluidState_.pressure(phaseIdx)));
            assert(Opm::isfinite(fluidState_.invB(phaseIdx)));
        }
        assert(Opm::isfinite(fluidState_.Rs()));
        assert(Opm::isfinite(fluidState_.Rv()));
#endif
    }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::mobility
     */
    const Evaluation& mobility(unsigned phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the index of the PVT region used to calculate the thermodynamic
     *        quantities.
     *
     * This allows to specify different Pressure-Volume-Temperature (PVT) relations in
     * different parts of the spatial domain. Note that this concept should be seen as a
     * work-around of the fact that the black-oil model does not capture the
     * thermodynamics well enough. (Because there is, err, only a single real world with
     * in which all substances follow the same physical laws and hence the same
     * thermodynamics.) Anyway: Since the ECL file format uses multiple PVT regions, we
     * support it as well in our black-oil model. (Note that, if it is not explicitly
     * specified, the PVT region index is 0.)
     */
    auto pvtRegionIndex() const
        -> decltype(std::declval<FluidState>().pvtRegionIndex())
    { return fluidState_.pvtRegionIndex(); }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
     */
    Evaluation relativePermeability(unsigned phaseIdx) const
    {
        // warning: slow
        return fluidState_.viscosity(phaseIdx)*mobility(phaseIdx);
    }

    /*!
     * \brief Returns the porosity of the rock at reference conditions.
     *
     * I.e., the porosity of rock which is not perturbed by pressure and temperature
     * changes.
     */
    Scalar referencePorosity() const
    { return referencePorosity_; }

private:
    friend BlackOilSolventIntensiveQuantities<TypeTag>;
    friend BlackOilPolymerIntensiveQuantities<TypeTag>;
    friend BlackOilEnergyIntensiveQuantities<TypeTag>;
    friend BlackOilFoamIntensiveQuantities<TypeTag>;
    friend BlackOilBrineIntensiveQuantities<TypeTag>;


    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    FluidState fluidState_;
    Scalar referencePorosity_;
    Evaluation porosity_;
    Evaluation mobility_[numPhases];
};

} // namespace Opm

#endif
