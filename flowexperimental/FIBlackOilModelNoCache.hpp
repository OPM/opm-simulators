#ifndef OPM_FI_BLACK_OIL_MODEL_NOCACHE_HPP
#define OPM_FI_BLACK_OIL_MODEL_NOCACHE_HPP

#include <opm/simulators/flow/FIBlackoilModel.hpp>

namespace Opm {

template<typename TypeTag>
class FIBlackOilModelNoCache: public FIBlackOilModel<TypeTag>
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;

public:
    explicit FIBlackOilModelNoCache(Simulator& simulator)
        : FIBlackOilModel<TypeTag>(simulator)
    {}

    IntensiveQuantities intensiveQuantities(unsigned globalIdx, unsigned timeIdx) const
    {
        OPM_TIMEBLOCK_LOCAL(intensiveQuantitiesNoCache, Subsystem::PvtProps | Subsystem::SatProps);
        const auto& primaryVar = this->solution(timeIdx)[globalIdx];
        const auto& problem = this->simulator_.problem();
        if (!(this->enableIntensiveQuantityCache_) ||
            !(this->intensiveQuantityCacheUpToDate_[timeIdx][globalIdx])) {
            IntensiveQuantities intQuants;
            intQuants.update(problem,primaryVar, globalIdx, timeIdx);
            return intQuants;// reqiored for updating extrution factor
        } else {
            IntensiveQuantities intQuants = (this->intensiveQuantityCache_[timeIdx][globalIdx]);
            return intQuants;
        }
    }
};

} // namespace Opm

#endif
