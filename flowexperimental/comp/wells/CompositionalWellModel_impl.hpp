
#include <opm/input/eclipse/Schedule/Schedule.hpp>

namespace Opm
{

template <typename TypeTag>
CompositionalWellModel<TypeTag>::
CompositionalWellModel(Simulator& simulator)
    : simulator_(simulator)
      //, schedule_(schedule)
{
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
beginTimeStep()
{
    const int reportStepIdx = simulator_.episodeIndex();
    createWellContainer(reportStepIdx);
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
createWellContainer(unsigned report_step)
{
    const auto& schedule = simulator_.vanguard().schedule();
    auto wells = schedule.getWells(report_step);
    well_container_.clear();
    for (const auto& well : wells) {
        well_container_.emplace_back(std::make_shared<CompWell<Scalar>>(well.name()));
    }
}

} // end of namespace Opm