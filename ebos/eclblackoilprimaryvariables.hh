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
 * \copydoc Opm::BlackOilPrimaryVariables
 */
#ifndef EWOMS_ECL_BLACK_OIL_PRIMARY_VARIABLES_HH
#define EWOMS_ECL_BLACK_OIL_PRIMARY_VARIABLES_HH

#include <opm/models/blackoil/blackoilprimaryvariables.hh>
namespace Opm{
/*!
 * \ingroup BlackOilModel
 *
 * \brief Represents the primary variables used by the black-oil model.
 */
template <class TypeTag>
class EclBlackOilPrimaryVariables : public BlackOilPrimaryVariables<TypeTag>
{
    using BaseType = BlackOilPrimaryVariables<TypeTag>;

    using typename BaseType::Scalar;
    using typename BaseType::Problem;
    using typename BaseType::Indices;
    using typename BaseType::FluidSystem;
    using typename BaseType::MaterialLawParams;
    using typename BaseType::EnergyModule;

    using BaseType::numPhases;
    using BaseType::oilPhaseIdx;
    using BaseType::gasPhaseIdx;
    using BaseType::waterPhaseIdx;
    using BaseType::compositionSwitchEnabled;


    using BaseType::waterEnabled;
    using BaseType::solventSaturation_;
    using BaseType::asImp_;

    using BaseType::primaryVarsMeaning_;
    using BaseType::computeCapillaryPressures_;


    using BaseType::pvtRegionIdx_;

    using BaseType::waterSaturationIdx;
    using BaseType::pressureSwitchIdx;
    using BaseType::compositionSwitchIdx;

public:
    using BaseType::Sw_pg_Rv;
    using BaseType::Sw_po_Rs;
    using BaseType::Sw_po_Sg;

    using BaseType::checkDefined;
    using BaseType::OnePhase_p;
    using BaseType::primaryVarsMeaning;
    using BaseType::setPrimaryVarsMeaning;

    using BaseType::operator=;

    enum SimulationType {
        Implicit,
        Seq
    };
    
    EclBlackOilPrimaryVariables()
        : BaseType()
    {
        Opm::Valgrind::SetUndefined(*this);
        simulationType_ = SimulationType::Implicit;
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(Scalar)
     */
    EclBlackOilPrimaryVariables(Scalar value)
        : BaseType(value)
    {
        Opm::Valgrind::SetUndefined(simulationType_);
    }

    /*!
     * \copydoc ImmisciblePrimaryVariables::ImmisciblePrimaryVariables(const ImmisciblePrimaryVariables& )
     */
    EclBlackOilPrimaryVariables(const EclBlackOilPrimaryVariables& value) = default;


    SimulationType simulationType() const
    { return simulationType_; }
    

    /*
    void setSimulaiontType(SimulationType simulationType)
    { simulaionType_ = simulationType; }
    */

    /*!
     * \copydoc ImmisciblePrimaryVariables::assignNaive
     */
    template <class FluidState>
    void assignNaive(const FluidState& fluidState)
    {
        using ConstEvaluation = typename std::remove_reference<typename FluidState::Scalar>::type;
        using FsEvaluation = typename std::remove_const<ConstEvaluation>::type;
        using FsToolbox = typename Opm::MathToolbox<FsEvaluation>;

        bool gasPresent = FluidSystem::phaseIsActive(gasPhaseIdx)?(fluidState.saturation(gasPhaseIdx) > 0.0):false;
        bool oilPresent = FluidSystem::phaseIsActive(oilPhaseIdx)?(fluidState.saturation(oilPhaseIdx) > 0.0):false;
        static const Scalar thresholdWaterFilledCell = 1.0 - 1e-6;
        bool onlyWater = FluidSystem::phaseIsActive(waterPhaseIdx)?(fluidState.saturation(waterPhaseIdx) > thresholdWaterFilledCell):false;

        // deal with the primary variables for the energy extension
        EnergyModule::assignPrimaryVars(*this, fluidState);

        // determine the meaning of the primary variables
        if (FluidSystem::numActivePhases() == 1) {
            primaryVarsMeaning_ = OnePhase_p;
        }
        else if ((gasPresent && oilPresent) || (onlyWater && FluidSystem::phaseIsActive(oilPhaseIdx))) {
            // gas and oil: both hydrocarbon phases are in equilibrium (i.e., saturated
            // with the "protagonist" component of the other phase.)
            primaryVarsMeaning_ = Sw_po_Sg;
        }
        else if (oilPresent) {
            // only oil: if dissolved gas is enabled, we need to consider the oil phase
            // composition, if it is disabled, the gas component must stick to its phase
            if (FluidSystem::enableDissolvedGas())
                primaryVarsMeaning_ = Sw_po_Rs;
            else
                primaryVarsMeaning_ = Sw_po_Sg;
        }
        else {
            assert(gasPresent);
            // only gas: if vaporized oil is enabled, we need to consider the gas phase
            // composition, if it is disabled, the oil component must stick to its phase
            if (FluidSystem::enableVaporizedOil())
                primaryVarsMeaning_ = Sw_pg_Rv;
            else
                primaryVarsMeaning_ = Sw_po_Sg;
        }

        // assign the actual primary variables
        if (primaryVarsMeaning() == OnePhase_p) {
            if (waterEnabled) {
                (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));
                (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(waterPhaseIdx));
            } else {
                throw std::logic_error("For single-phase runs, only pure water is presently allowed.");
            }

        }
        else if (primaryVarsMeaning() == Sw_po_Sg) {
            if (waterEnabled)
                (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));
            (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(oilPhaseIdx));
            if( compositionSwitchEnabled )
                (*this)[compositionSwitchIdx] = FsToolbox::value(fluidState.saturation(gasPhaseIdx));
        }
        else if (primaryVarsMeaning() == Sw_po_Rs) {
            const auto& Rs = Opm::BlackOil::getRs_<FluidSystem, FluidState, Scalar>(fluidState, pvtRegionIdx_);

            if (waterEnabled)
                (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));
            (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(oilPhaseIdx));
            if( compositionSwitchEnabled )
                (*this)[compositionSwitchIdx] = Rs;
        }
        else {
            assert(primaryVarsMeaning() == Sw_pg_Rv);

            const auto& Rv = Opm::BlackOil::getRv_<FluidSystem, FluidState, Scalar>(fluidState, pvtRegionIdx_);
            if (waterEnabled)
                (*this)[waterSaturationIdx] = FsToolbox::value(fluidState.saturation(waterPhaseIdx));

            (*this)[pressureSwitchIdx] = FsToolbox::value(fluidState.pressure(gasPhaseIdx));
            if( compositionSwitchEnabled )
                (*this)[compositionSwitchIdx] = Rv;
        }

        checkDefined();
    }


    /*!
     * \brief Adapt the interpretation of the switching variables to be physically
     *        meaningful.
     *
     * If the meaning of the primary variables changes, their values are also adapted in a
     * meaningful manner. (e.g. if the gas phase appears and the composition switching
     * variable changes its meaning from the gas dissolution factor Rs to the gas
     * saturation Sg, the value for this variable is set to zero.)
     * A Scalar eps can be passed to make the switching condition more strict.
     * Useful for avoiding ocsilation in the primaryVarsMeaning.
     *
     * \return true Iff the interpretation of one of the switching variables was changed
     */
    bool adaptPrimaryVariables(const Problem& problem, unsigned globalDofIdx, Scalar eps = 0.0)
    {
        const LinearizationType linearizationType = problem.simulator().model().linearizer().getLinearizationType();
        static const Scalar thresholdWaterFilledCell = 1.0 - eps;

        // this function accesses quite a few black-oil specific low-level functions
        // directly for better performance (instead of going the canonical way through
        // the IntensiveQuantities). The reason is that most intensive quantities are not
        // required to be able to decide if the primary variables needs to be switched or
        // not, so it would be a waste to compute them.
        if (primaryVarsMeaning() == OnePhase_p){
            return false;
        }
        Scalar pressure=-1.0;
        Scalar totalsaturation=-1.0;
        if(simulationType_ == SimulationType::Implicit){
            if(linearizationType.type == LinearizationType::pressure){
                return false; //No switchin during pressure solve;
            }
        }
        
        if(simulationType_ == SimulationType::Implicit){
            pressure = (*this)[Indices::pressureSwitchIdx];
            totalsaturation = problem.totalSaturation(globalDofIdx,/*timeidx*/ 0);//NB hack
        }else if (simulationType_ == SimulationType::Seq){
            if(primaryVarsMeaning() == Sw_po_Sg){
                pressure = problem.pressure(globalDofIdx, oilPhaseIdx,/*timeidx*/ 0);//NB hack
            }else if(primaryVarsMeaning() == Sw_po_Rs){
                pressure = problem.pressure(globalDofIdx, oilPhaseIdx,/*timeidx*/ 0);//NB hack
            }else if(primaryVarsMeaning() == Sw_pg_Rv){
                pressure = problem.pressure(globalDofIdx, gasPhaseIdx,/*timeidx*/ 0);//NB hack
            }else{
                assert(false);
            }
            totalsaturation = (*this)[Indices::pressureSwitchIdx];
        }else{
            assert(false);
        }
        totalsaturation = std::max(0.1,totalsaturation);
        totalsaturation = std::min(totalsaturation,5.0);
        
        if(linearizationType.type == LinearizationType::seqtransport){
            (*this)[Indices::pressureSwitchIdx] = totalsaturation;
            this->simulationType_ = SimulationType::Seq;
        }else{
            (*this)[Indices::pressureSwitchIdx] = pressure;
            simulationType_ = SimulationType::Implicit;
        }
        
        Scalar Sw = 0.0;
        if (waterEnabled)
            Sw = (*this)[Indices::waterSaturationIdx];

        if (primaryVarsMeaning() == Sw_po_Sg) {

            // special case for cells with almost only water
            if (Sw >= thresholdWaterFilledCell) {

                // make sure water saturations does not exceed 1.0
                if (waterEnabled)
                    (*this)[Indices::waterSaturationIdx] = 1.0;
                // the hydrocarbon gas saturation is set to 0.0
                if (compositionSwitchEnabled)
                    (*this)[Indices::compositionSwitchIdx] = 0.0;

                return false;
            }

            // phase equilibrium, i.e., both hydrocarbon phases are present.
            Scalar Sg = 0.0;
            if (compositionSwitchEnabled)
                Sg = (*this)[Indices::compositionSwitchIdx];

            Scalar So = 1.0 - Sw - Sg - solventSaturation_();

            Scalar So2 = 1.0 - Sw - solventSaturation_();
            if (Sg < -eps && So2 > 0.0 && FluidSystem::enableDissolvedGas()) {
                // the hydrocarbon gas phase disappeared and some oil phase is left,
                // i.e., switch the primary variables to { Sw, po, Rs }.
                //
                // by a "lucky" coincidence the pressure switching variable already
                // represents the oil phase pressure, so we do not need to change
                // this. For the gas dissolution factor, we use the low-level blackoil
                // PVT objects to calculate the mole fraction of gas saturated oil.
                Scalar po;
                if (simulationType_ == SimulationType::Implicit){
                    po = (*this)[Indices::pressureSwitchIdx];
                    assert(po>2.0);
                }else{
                    po = problem.pressure(globalDofIdx, oilPhaseIdx,/*timeidx*/ 0);
                }
                Scalar T = asImp_().temperature_();
                Scalar SoMax = problem.maxOilSaturation(globalDofIdx);
                Scalar RsMax = problem.maxGasDissolutionFactor(/*timeIdx=*/0, globalDofIdx);
                Scalar RsSat = FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_,
                                                                                   T,
                                                                                   po,
                                                                                   So2,
                                                                                   SoMax);

                setPrimaryVarsMeaning(Sw_po_Rs);
                if (compositionSwitchEnabled)
                    (*this)[Indices::compositionSwitchIdx] =
                        std::min(RsMax, RsSat);
                return true;
            }

            Scalar Sg2 = 1.0 - Sw - solventSaturation_();
            if (So < -eps && Sg2 > 0.0 && FluidSystem::enableVaporizedOil()) {
                // the oil phase disappeared and some hydrocarbon gas phase is still
                // present, i.e., switch the primary variables to { Sw, pg, Rv }.
                Scalar po;
                if (simulationType_ == SimulationType::Implicit){
                    po = (*this)[Indices::pressureSwitchIdx];
                    assert(po>2.0);
                }else{
                    po = problem.pressure(globalDofIdx, oilPhaseIdx,/*timeidx*/ 0);
                }
                // we only have the oil pressure readily available, but we need the gas
                // pressure, i.e. we must determine capillary pressure
                Scalar pC[numPhases] = { 0.0 };
                const MaterialLawParams& matParams = problem.materialLawParams(globalDofIdx);
                computeCapillaryPressures_(pC, /*So=*/0.0, Sg2 + solventSaturation_(), Sw, matParams);
                Scalar pg = po + (pC[gasPhaseIdx] - pC[oilPhaseIdx]);

                // we start at the Rv value that corresponds to that of oil-saturated
                // hydrocarbon gas
                Scalar T = asImp_().temperature_();
                Scalar SoMax = problem.maxOilSaturation(globalDofIdx);
                Scalar RvMax = problem.maxOilVaporizationFactor(/*timeIdx=*/0, globalDofIdx);
                Scalar RvSat =
                    FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_,
                                                                         T,
                                                                         pg,
                                                                         Scalar(0),
                                                                         SoMax);
                setPrimaryVarsMeaning(Sw_pg_Rv);
                if (not(linearizationType.type == Opm::LinearizationType::seqtransport)) {
                    (*this)[Indices::pressureSwitchIdx] = pg;
                }
                if (compositionSwitchEnabled)
                    (*this)[Indices::compositionSwitchIdx] = std::min(RvMax, RvSat);

                return true;
            }

            return false;
        }
        else if (primaryVarsMeaning() == Sw_po_Rs) {
            assert(compositionSwitchEnabled);

            // special case for cells with almost only water
            if (Sw >= thresholdWaterFilledCell) {
                // switch back to phase equilibrium mode if the oil phase vanishes (i.e.,
                // the water-only case)
                setPrimaryVarsMeaning(Sw_po_Sg);
                if (waterEnabled)
                    (*this)[Indices::waterSaturationIdx] = 1.0; // water saturation

                (*this)[Indices::compositionSwitchIdx] = 0.0; // hydrocarbon gas saturation

                return true;
            }

            // Only the oil and the water phases are present. The hydrocarbon gas phase
            // appears as soon as more of the gas component is present in the oil phase
            // than what saturated oil can hold.
            Scalar T = asImp_().temperature_();
            Scalar po;
            if (simulationType_ == SimulationType::Implicit) {
                po = (*this)[Indices::pressureSwitchIdx];
                assert(po>2.0);
            } else {
                po = problem.pressure(globalDofIdx, oilPhaseIdx,/*timeidx*/ 0);                
            }
            Scalar So = 1.0 - Sw - solventSaturation_();
            Scalar SoMax = std::max(So, problem.maxOilSaturation(globalDofIdx));
            Scalar RsMax = problem.maxGasDissolutionFactor(/*timeIdx=*/0, globalDofIdx);
            Scalar RsSat =
                FluidSystem::oilPvt().saturatedGasDissolutionFactor(pvtRegionIdx_,
                                                                    T,
                                                                    po,
                                                                    So,
                                                                    SoMax);

            Scalar Rs = (*this)[Indices::compositionSwitchIdx];
            if (Rs > std::min(RsMax, RsSat*(1.0 + eps))) {
                // the gas phase appears, i.e., switch the primary variables to { Sw, po,
                // Sg }.
                setPrimaryVarsMeaning(Sw_po_Sg);
                (*this)[Indices::compositionSwitchIdx] = 0.0; // hydrocarbon gas saturation

                return true;
            }

            return false;
        }
        else {
            assert(primaryVarsMeaning() == Sw_pg_Rv);
            assert(compositionSwitchEnabled);

            Scalar pg;
            if (simulationType_ == SimulationType::Implicit) {
                // This should be the gas pressure
                pg = (*this)[Indices::pressureSwitchIdx];
                assert(pg>2.0);
            } else {
                pg = problem.pressure(globalDofIdx, gasPhaseIdx,/*timeidx*/ 0);
            }
            Scalar Sg = 1.0 - Sw - solventSaturation_();

            // special case for cells with almost only water
            if (Sw >= thresholdWaterFilledCell) {
                // switch to phase equilibrium mode because the hydrocarbon gas phase
                // disappears. here we need the capillary pressures to calculate the oil
                // phase pressure using the gas phase pressure
                Scalar pC[numPhases] = { 0.0 };
                const MaterialLawParams& matParams = problem.materialLawParams(globalDofIdx);
                computeCapillaryPressures_(pC,
                                           /*So=*/0.0,
                                           /*Sg=*/Sg + solventSaturation_(),
                                           Sw,
                                           matParams);
                Scalar po = pg + (pC[oilPhaseIdx] - pC[gasPhaseIdx]);

                setPrimaryVarsMeaning(Sw_po_Sg);
                if (waterEnabled)
                    (*this)[Indices::waterSaturationIdx] = 1.0;
                if(not(linearizationType.type == Opm::LinearizationType::seqtransport)){
                    (*this)[Indices::pressureSwitchIdx] = po;
                }
                (*this)[Indices::compositionSwitchIdx] = 0.0; // hydrocarbon gas saturation

                return true;
            }

            // Only the gas and the water phases are present. The oil phase appears as
            // soon as more of the oil component is present in the hydrocarbon gas phase
            // than what saturated gas contains. Note that we use the blackoil specific
            // low-level PVT objects here for performance reasons.
            Scalar T = asImp_().temperature_();
            Scalar SoMax = problem.maxOilSaturation(globalDofIdx);
            Scalar RvMax = problem.maxOilVaporizationFactor(/*timeIdx=*/0, globalDofIdx);
            Scalar RvSat =
                FluidSystem::gasPvt().saturatedOilVaporizationFactor(pvtRegionIdx_,
                                                                     T,
                                                                     pg,
                                                                     /*So=*/Scalar(0.0),
                                                                     SoMax);

            Scalar Rv = (*this)[Indices::compositionSwitchIdx];
            if (Rv > std::min(RvMax, RvSat*(1.0 + eps))) {
                // switch to phase equilibrium mode because the oil phase appears. here
                // we also need the capillary pressures to calculate the oil phase
                // pressure using the gas phase pressure
                Scalar pC[numPhases] = { 0.0 };
                const MaterialLawParams& matParams = problem.materialLawParams(globalDofIdx);
                computeCapillaryPressures_(pC,
                                           /*So=*/0.0,
                                           /*Sg=*/Sg + solventSaturation_(),
                                           Sw,
                                           matParams);
                Scalar po = pg + (pC[oilPhaseIdx] - pC[gasPhaseIdx]);

                setPrimaryVarsMeaning(Sw_po_Sg);
                if (linearizationType.type != Opm::LinearizationType::seqtransport) {
                    (*this)[Indices::pressureSwitchIdx] = po;
                }
                (*this)[Indices::compositionSwitchIdx] = Sg; // hydrocarbon gas saturation

                return true;
            }

            return false;
        }

        assert(false);
        return false;
    }

protected:
    SimulationType simulationType_;
};


} // namespace Opm

#endif
