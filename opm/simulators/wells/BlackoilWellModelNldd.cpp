/*
  Copyright 2016 - 2019 SINTEF Digital, Mathematics & Cybernetics.
  Copyright 2016 - 2018 Equinor ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 Norce AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>
#include <opm/simulators/wells/BlackoilWellModelNldd.hpp>

namespace Opm {

template <typename Scalar>
std::vector<Scalar>
BlackoilWellModelNlddGeneric<Scalar>::
getPrimaryVarsDomain(const int domainIdx) const
{
    std::vector<Scalar> ret;
    for (const auto& well : genWellModel_.genericWells()) {
        if (this->well_domain_.at(well->name()) == domainIdx) {
            const auto& pv = well->getPrimaryVars();
            ret.insert(ret.end(), pv.begin(), pv.end());
        }
    }
    return ret;
}

template <typename Scalar>
void
BlackoilWellModelNlddGeneric<Scalar>::
setPrimaryVarsDomain(const int domainIdx, const std::vector<Scalar>& vars)
{
    std::size_t offset = 0;
    for (const auto& well : genWellModel_.genericWells()) {
        if (this->well_domain_.at(well->name()) == domainIdx) {
            int num_pri_vars = well->setPrimaryVars(vars.begin() + offset);
            offset += num_pri_vars;
        }
    }
    assert(offset == vars.size());
}

template <typename Scalar>
void
BlackoilWellModelNlddGeneric<Scalar>::
findWellDomains(const std::vector<const SubDomainIndices*>& domains)
{
    // TODO: This loop nest may be slow for very large numbers of
    // domains and wells, but that has not been observed on tests so
    // far.  Using the partition vector instead would be faster if we
    // need to change.
    for (const auto& wellPtr : genWellModel_.genericWells()) {
        const int first_well_cell = wellPtr->cells().front();
        for (const auto& domain : domains) {
            auto cell_present = [domain](const auto cell)
            {
                return std::binary_search(domain->cells.begin(),
                                          domain->cells.end(), cell);
            };

            if (cell_present(first_well_cell)) {
                // Assuming that if the first well cell is found in a domain,
                // then all of that well's cells are in that same domain.
                this->well_domain_[wellPtr->name()] = domain->index;

                // Verify that all of that well's cells are in that same domain.
                for (int well_cell : wellPtr->cells()) {
                    if (! cell_present(well_cell)) {
                        OPM_THROW(std::runtime_error,
                                  fmt::format("Well '{}' found on multiple domains.",
                                              wellPtr->name()));
                    }
                }
            }
        }
    }
}

template <typename Scalar>
void
BlackoilWellModelNlddGeneric<Scalar>::
logDomains() const
{
    // Write well/domain info to the DBG file.
    const int rank = genWellModel_.comm().rank();
    DeferredLogger local_log;
    if (!this->well_domain_.empty()) {
        std::ostringstream os;
        os << "Well name      Rank      Domain\n";
        for (const auto& [wname, domain] : this->well_domain_) {
            os << wname << std::setw(19 - wname.size()) << rank << std::setw(12) << domain << '\n';
        }
        local_log.debug(os.str());
    }
    auto global_log = gatherDeferredLogger(local_log, genWellModel_.comm());
    if (genWellModel_.terminalOutput()) {
        global_log.logMessages();
    }
}

template <typename Scalar>
void
BlackoilWellModelNlddGeneric<Scalar>::
calcLocalIndices(const std::vector<const SubDomainIndices*>& domains)
{
    well_local_cells_.clear();
    well_local_cells_.reserve(genWellModel_.genericWells().size(), 10);
    std::vector<int> local_cells;
    for (const auto& well : genWellModel_.genericWells()) {
        const auto& global_cells = well->cells();
        const int domain_index = this->well_domain_.at(well->name());
        const auto& domain_cells = domains[domain_index]->cells;
        local_cells.resize(global_cells.size());

        // find the local cell index for each well cell in the domain
        // assume domain_cells is sorted
        for (size_t i = 0; i < global_cells.size(); ++i) {
            auto it = std::lower_bound(domain_cells.begin(), domain_cells.end(), global_cells[i]);
            if (it != domain_cells.end() && *it == global_cells[i]) {
                local_cells[i] = std::distance(domain_cells.begin(), it);
            } else {
                OPM_THROW(std::runtime_error, fmt::format("Cell {} not found in domain {}",
                                                          global_cells[i], domain_index));
            }
        }
        well_local_cells_.appendRow(local_cells.begin(), local_cells.end());
    }
}

template <typename Scalar>
void
BlackoilWellModelNlddGeneric<Scalar>::
calcDomains(const std::vector<const SubDomainIndices*>& domains)
{
    const Opm::Parallel::Communication& comm = genWellModel_.comm();

    OPM_BEGIN_PARALLEL_TRY_CATCH();
    this->findWellDomains(domains);
    OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel::setupDomains(): "
                               "well found on multiple domains.", comm);

    // Write well/domain info to the DBG file.
    this->logDomains();

    // Pre-calculate the local cell indices for each well
    this->calcLocalIndices(domains);
}

template class BlackoilWellModelNlddGeneric<double>;

#if FLOW_INSTANTIATE_FLOAT
template class BlackoilWellModelNlddGeneric<float>;
#endif

} // namespace Opm
