#ifndef ECL_PROBLEM_OIL_STANDARD_HPP
#define ECL_PROBLEM_OIL_STANDARD_HPP
namespace Opm{
    template<typename TypeTag>
    class EclProblemStandard: public EclProblem<TypeTag>{
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
        using DirectionalMobilityPtr = ::Opm::Utility::CopyablePtr<DirectionalMobility<TypeTag, Evaluation>>;
        using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using SolidEnergyLaw = GetPropType<TypeTag, Properties::SolidEnergyLaw>;
        using EclMaterialLawManager = typename GetProp<TypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
        using EclThermalLawManager = typename GetProp<TypeTag, Properties::SolidEnergyLaw>::EclThermalLawManager;
        using MaterialLawParams = typename EclMaterialLawManager::MaterialLawParams;
        using SolidEnergyLawParams = typename EclThermalLawManager::SolidEnergyLawParams;
        using ThermalConductionLawParams = typename EclThermalLawManager::ThermalConductionLawParams;
        enum { numPhases = FluidSystem::numPhases };
    public:
        EclProblemStandard(Simulator& simulator): EclProblem<TypeTag>(simulator){
        }

    public:
        /*
          Function which only update the intensive quantites if hystresis paramters in given cell is changed.
         */
        bool updateExplicitQuantities_(){
            OPM_TIMEBLOCK(updateExpliciteQuantitiesSimple);
            const auto& primaryVars = this->model().solution(/*timeIdx*/0);
            const auto& problem = this->simulator().problem();
            const auto& model = this->simulator().model();
            const auto& vanguard = this->simulator().vanguard();
            size_t numGridDof = primaryVars.size();
            std::string failureMsg("Error updating explicit quantities");
            bool any_changed = false;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (unsigned dofIdx = 0; dofIdx < numGridDof; ++dofIdx) {
                bool changed=false;
                // using the cached explitely need for updating
                const auto& iq = model.intensiveQuantities(dofIdx, /*timeIdx=*/ 0);
                if (!this->maxWaterSaturation_.empty()){
                    const bool invalidateFromMaxWaterSat = this->updateMaxWaterSaturation_(dofIdx,iq);
                    changed = changed || invalidateFromMaxWaterSat;
                }
                if (!this->minOilPressure_.empty()){
                    const bool invalidateFromMinPressure = this->updateMinPressure_(dofIdx,iq);
                    changed = changed || invalidateFromMinPressure;
                }

                // update hysteresis and max oil saturation used in vappars
                if (this->materialLawManager()->enableHysteresis()){
                    //need bools in opm-material
                    const bool invalidateFromHyst = this->materialLawManager()->updateHysteresis(iq.fluidState(),dofIdx);
                    changed = changed || invalidateFromHyst;
                }
                int episodeIdx = this->episodeIndex();
                if (this->vapparsActive(episodeIdx)) {
                    const bool invalidateFromMaxOilSat = this->updateMaxOilSaturation_(dofIdx,iq);
                    changed = changed || invalidateFromMaxOilSat;
                }

                // the derivatives may have change
                any_changed = changed || any_changed;
                if(changed){
                    OPM_TIMEBLOCK_LOCAL(updateIntensiveQuantityHyst);
                    const auto& primaryVar = primaryVars[dofIdx];
                    //auto& intquant = this->intensiveQuantityCache_[timeIdx][dofIdx];
                    // to get it working with and without cached normally &iq = iq
                    // auto* iql = const_cast<IntensiveQuantities&>(model.cachedIntensiveQuantities(dofIdx, /*timeIdx=*/ 0));
                    // if(iql){
                    //     iql->update(problem, primaryVar, dofIdx,/*timeIdx*/0);
                    // }
                    auto& iql = const_cast<IntensiveQuantities&>(iq);
                    iql.update(problem, primaryVar, dofIdx,/*timeIdx*/0);
                    model.updateCachedIntensiveQuantities(iql,dofIdx,/*timeIdx*/0);

                }
                if constexpr (getPropValue<TypeTag, Properties::EnablePolymer>()){
                    updateMaxPolymerAdsorption_(dofIdx,iq);
                }
            }
            return any_changed;
        }
        private:
        friend EclProblem<TypeTag>;

    };


}


#endif
