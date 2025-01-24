/*
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2017 Statoil ASA.

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

#ifndef OPM_WELLCONNECTIONAUXILIARYMODULE_HEADER_INCLUDED
#define OPM_WELLCONNECTIONAUXILIARYMODULE_HEADER_INCLUDED

#include <opm/models/discretization/common/baseauxiliarymodule.hh>

#include <opm/simulators/flow/SubDomain.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#if HAVE_MPI
#include <opm/simulators/utils/MPISerializer.hpp>
#endif

namespace Opm {

template<class TypeTag, class Model>
class WellConnectionAuxiliaryModule : public BaseAuxiliaryModule<TypeTag>
{
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;

public:
    using NeighborSet = typename
        ::Opm::BaseAuxiliaryModule<TypeTag>::NeighborSet;

    using Domain = SubDomain<Grid>;

    WellConnectionAuxiliaryModule(Model& model, Parallel::Communication comm)
        : model_(model)
        , lin_comm_(std::move(comm))
    {
    }

    unsigned numDofs() const override
    {
        // No extra dofs are inserted for wells.
        return 0;
    }

    void addNeighbors(std::vector<NeighborSet>& neighbors) const override
    {
        if (!model_.addMatrixContributions()) {
            return;
        }

        // Create cartesian to compressed mapping
        const auto& schedule_wells = model_.schedule().getWellsatEnd();
        auto possibleFutureConnections = model_.schedule().getPossibleFutureConnections();

#if HAVE_MPI
        // Communicate Map to other processes, since it is only available on rank 0
        Parallel::MpiSerializer ser(lin_comm_);
        ser.broadcast(possibleFutureConnections);
#endif
        // initialize the additional cell connections introduced by wells.
        for (const auto& well : schedule_wells)
        {
            std::vector<int> wellCells = model_.getCellsForConnections(well);
            // Now add the cells of the possible future connections
            const auto possibleFutureConnectionSetIt = possibleFutureConnections.find(well.name());
            if (possibleFutureConnectionSetIt != possibleFutureConnections.end()) {
                for (const auto& global_index : possibleFutureConnectionSetIt->second) {
                    int compressed_idx = model_.compressedIndexForInterior(global_index);
                    if (compressed_idx >= 0) { // Ignore connections in inactive/remote cells.
                        wellCells.push_back(compressed_idx);
                    }
                }
            }
            for (int cellIdx : wellCells) {
                neighbors[cellIdx].insert(wellCells.begin(),
                                          wellCells.end());
            }
        }
    }

    void applyInitial() override
    {}

    void linearize(SparseMatrixAdapter& jacobian, GlobalEqVector& res) override
    {
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        for (const auto& well : model_) {
            this->linearizeSingleWell(jacobian, res, well);
        }
        OPM_END_PARALLEL_TRY_CATCH("BlackoilWellModel::linearize failed: ", lin_comm_);
    }

    void postSolve(GlobalEqVector& deltaX) override
    {
        model_.recoverWellSolutionAndUpdateWellState(deltaX);
    }

    void linearizeDomain(const Domain& domain,
                         SparseMatrixAdapter& jacobian,
                         GlobalEqVector& res)
    {
        // Note: no point in trying to do a parallel gathering
        // try/catch here, as this function is not called in
        // parallel but for each individual domain of each rank.
        for (const auto& well : model_) {
            if (model_.well_domain().at(well->name()) == domain.index) {
                this->linearizeSingleWell(jacobian, res, well);
            }
        }
    }

    void postSolveDomain(const GlobalEqVector& deltaX, const Domain& domain)
    {
        model_.recoverWellSolutionAndUpdateWellStateDomain(deltaX, domain.index);
    }

    template <class Restarter>
    void deserialize(Restarter& /* res */)
    {
        // TODO (?)
    }

    /*!
     * \brief This method writes the complete state of the well
     *        to the harddisk.
     */
    template <class Restarter>
    void serialize(Restarter& /* res*/)
    {
        // TODO (?)
    }

private:
    template<class WellType>
    void linearizeSingleWell(SparseMatrixAdapter& jacobian,
                             GlobalEqVector& res,
                             const WellType& well)
    {
        if (model_.addMatrixContributions()) {
            well->addWellContributions(jacobian);
        }

        const auto& cells = well->cells();
        linearize_res_local_.resize(cells.size());

        for (size_t i = 0; i < cells.size(); ++i) {
           linearize_res_local_[i] = res[cells[i]];
        }

        well->apply(linearize_res_local_);

        for (size_t i = 0; i < cells.size(); ++i) {
            res[cells[i]] = linearize_res_local_[i];
        }
    }

    Model& model_;
    GlobalEqVector linearize_res_local_{};
    Parallel::Communication lin_comm_;
};

} // end namespace OPM
#endif
