/*
  Copyright 2024, SINTEF Digital

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

#ifndef OPM_COMPOSITIONAL_WELL_MODEL_HPP
#define OPM_COMPOSITIONAL_WELL_MODEL_HPP

#include <opm/output/data/Wells.hpp>

#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <flowexperimental/comp/wells/CompWell.hpp>
#include <flowexperimental/comp/wells/CompWellState.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/simulators/wells/WellConnectionAuxiliaryModule.hpp>

#include <opm/simulators/wells/PerforationData.hpp>


#include <vector>

namespace Opm {

class Schedule;

template<typename TypeTag>
class CompWellModel : WellConnectionAuxiliaryModule<TypeTag, CompWellModel<TypeTag>>
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using CompConnectionData = PerforationData<Scalar>;

    // using NeighborSet = typename BaseAuxiliaryModule<TypeTag>::NeighborSet;

public:
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    using WellConnectionModule = WellConnectionAuxiliaryModule<TypeTag, CompWellModel<TypeTag>>;

    using Indices = GetPropType<TypeTag, Properties::Indices>;
    static const int numEq = Indices::numEq;
    using VectorBlockType = Dune::FieldVector<Scalar, numEq>;
    typedef Dune::BlockVector<VectorBlockType> BVector;


    // TODO: Scalar will probably to be TypeTag later
    using CompWellPtr = std::shared_ptr<CompWell<TypeTag> >;
    explicit CompWellModel(Simulator& /*simulator*/);

    // No extra dofs are inserted for wells. (we use a Schur complement.)
    [[nodiscard]] unsigned numDofs() const override
    { return 0; }

//    void addNeighbors(std::vector<NeighborSet>& /* neighbors */) const override
//    {}

    void applyInitial() override {}

    // void linearize(SparseMatrixAdapter& /*matrix*/, GlobalEqVector& /*residual*/) override;

    template <class Restarter>
    void serialize(Restarter& /*res*/)
    {}

    template <class Restarter>
    void deserialize(Restarter& /*res*/)
    {}


    void beginEpisode() { beginReportStep(simulator_.episodeIndex()); }
    void beginReportStep(unsigned report_step);
    void beginTimeStep();
    void beginIteration();

    void init();
    void endIteration() const {}
    void endTimeStep() {}
    void endEpisode() {}

    void computeTotalRatesForDof(RateVector& /*rate*/, unsigned /*globalIdx*/) const;
    //
    [[nodiscard]] data::Wells wellData() const;

    [[nodiscard]] data::WellBlockAveragePressures wellBlockAveragePressures() const {
         return data::WellBlockAveragePressures{};
    }
    [[nodiscard]] data::GroupAndNetworkValues groupAndNetworkData(const int&) const {
         return data::GroupAndNetworkValues{};
    }
    [[nodiscard]] WellTestState wellTestState() const {
         return WellTestState{};
    }

    // using the solution x to recover the solution xw for wells and applying
    // xw to update Well State
    void recoverWellSolutionAndUpdateWellState(const BVector& x);

    // some functions to compile
    bool addMatrixContributions() const { return false; }
    const Schedule& schedule() const { return schedule_; }
    auto begin() const { return well_container_.begin(); }
    auto end() const { return well_container_.end(); }

    bool getWellConvergence() const;

    // the following functions are not used while added to avoid modifying WellConnectionAuxiliaryModule.hpp
    std::size_t compressedIndexForInterior(std::size_t cartesian_cell_idx) const;

    std::vector<int> getCellsForConnections(const Well& well) const;

private:
     Simulator& simulator_;
     const Schedule& schedule_;
     const SummaryState& summary_state_;
     const EclipseState& ecl_state_;
     const Parallel::Communication& comm_;

     // we might need something lighter
     const CompositionalConfig& comp_config_;

     // we will need two to handle the changes between time stepping
     CompWellState<FluidSystem, Scalar> comp_well_states_;

     // this is needed for parallel running, not all the wells will be in the same process
     std::vector<Well> wells_ecl_;
     std::vector<std::vector<CompConnectionData> > well_connection_data_;
     // const Schedule& schedule_;
     std::vector<CompWellPtr> well_container_;

     mutable BVector x_local_;

     std::size_t local_num_cells_{0};

     void createWellContainer();
     void initWellContainer();

     void initWellConnectionData();
     void initWellState();

     void assemble(const int iterationIdx,
                   const double dt);

     void calculateExplicitQuantities();
};

} // end of namespace Opm

#include "CompWellModel_impl.hpp"

#endif // OPM_COMPOSITIONAL_WELL_MODEL_HPP
