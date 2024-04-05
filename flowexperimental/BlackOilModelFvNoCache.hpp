#ifndef BLACK_OIL_MODEL_FV_NOCACHE_HPP
#define BLACK_OIL_MODEL_FV_NOCACHE_HPP
#include <opm/simulators/flow/FIBlackoilModel.hpp>
//#include "BlackOilModelFvFast.hpp"
namespace Opm{
    template<typename TypeTag>
    class BlackOilModelFvNoCache: public FIBlackOilModel<TypeTag>{
        using Parent = FIBlackOilModel<TypeTag>;
        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    public:
        BlackOilModelFvNoCache(Simulator& simulator): FIBlackOilModel<TypeTag>(simulator){
        }

        IntensiveQuantities intensiveQuantities(unsigned globalIdx, unsigned timeIdx) const{
            OPM_TIMEBLOCK_LOCAL(intensiveQuantitiesNoCache);
            const auto& primaryVar = this->solution(timeIdx)[globalIdx];
            const auto& problem = this->simulator_.problem();
            //IntensiveQuantities* intQuants = &(this->intensiveQuantityCache_[timeIdx][globalIdx]);
            if (!(this->enableIntensiveQuantityCache_) ||
                !(this->intensiveQuantityCacheUpToDate_[timeIdx][globalIdx])){
                IntensiveQuantities intQuants;
                intQuants.update(problem,primaryVar, globalIdx, timeIdx);
                return intQuants;// reqiored for updating extrution factor
            }else{
                IntensiveQuantities intQuants = (this->intensiveQuantityCache_[timeIdx][globalIdx]);
                return intQuants;
            }

        }

    };
}
#endif
