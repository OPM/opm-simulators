
#include <opm/input/eclipse/Schedule/Schedule.hpp>

namespace Opm
{

template <typename TypeTag>
CompositionalWellModel<TypeTag>::
CompositionalWellModel(Simulator& simulator)
    : simulator_(simulator)
    , schedule_(simulator.vanguard().schedule())
    , eclState_(simulator.vanguard().eclState())
    , comm_(simulator.gridView().comm())
      //, schedule_(schedule)
{
}

template <typename TypeTag>
void
CompositionalWellModel<TypeTag>::
beginReportStep(unsigned report_step)
{
    // TODO: not considering the parallel running yet
    wells_ecl_ = schedule_.getWells(report_step);
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
beginTimeStep()
{
    createWellContainer();
}

template <typename TypeTag>
void CompositionalWellModel<TypeTag>::
createWellContainer()
{
    const auto& schedule = simulator_.vanguard().schedule();
    const auto nw = wells_ecl_.size(); // not considering the parallel running yet
    well_container_.clear();
    for (int w = 0; w < nw; ++w) {
        well_container_.emplace_back(std::make_shared<CompWell<TypeTag>>(wells_ecl_[w], w));
    }
}

} // end of namespace Opm