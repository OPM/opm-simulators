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

#ifndef OPM_COMP_WELLINTERFACE_HPP
#define OPM_COMP_WELLINTERFACE_HPP

#include <opm/models/utils/propertysystem.hh>
#include <opm/simulators/wells/PerforationData.hpp>

#include <flowexperimental/comp/wells/SingleCompWellState.hpp>

#include <string>

namespace Opm {

template <typename TypeTag> // TODO: do we need to use TypeTag here?
class CompWellInterface
{
public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>; // FlashRateVector at the moment

    using VectorBlockType = Dune::FieldVector<Scalar, Indices::numEq>;
    using BVector = Dune::BlockVector<VectorBlockType>;
    using CompConnectionData = PerforationData<Scalar>;
    using SingleWellState = SingleCompWellState<FluidSystem, Scalar>;

    CompWellInterface(const Well& well,
                      const int index_of_well,
                      const std::vector<CompConnectionData>& well_connection_data);

    virtual ~CompWellInterface() = default;

    virtual void init();

    const std::string& name() const;

    virtual void calculateExplicitQuantities(const Simulator& simulator,
                                             const SingleWellState& well_state) = 0;

    virtual void updatePrimaryVariables(const Simulator& simulator,
                                        const SingleWellState& well_state) = 0;

    virtual bool iterateWellEq(const Simulator& simulator,
                               const Scalar dt,
                               SingleWellState& well_state) = 0;

    void addCellRates(RateVector& rates, unsigned cellIdx) const;

    const std::vector<std::size_t>& cells() const { return well_cells_; }

    virtual void apply(BVector& r) const = 0;

    /// using the solution x to recover the solution xw for wells and applying
    /// xw to update Well State
    virtual void recoverWellSolutionAndUpdateWellState(const BVector& x,
                                                       SingleWellState& well_state) = 0;
    virtual bool getConvergence() const = 0;

    virtual void addWellContributions(SparseMatrixAdapter&) const = 0;

protected:

    const Well& well_ecl_;
    int index_of_well_{-1};

    int number_of_connection_ {};
    Scalar reference_depth_ {}; // TODO: we might not need it since it in well_ecl_.

    std::vector<RateVector> connectionRates_;

    // cell index for each well connection
    // TODO: maybe it should be called connection_cells
    std::vector<std::size_t> well_cells_;
    // cell index for each well connection
    // TODO: should it called trans_index
    std::vector<Scalar> well_index_;
    std::vector<int> saturation_table_number_;

    // std::string name_;

};


} // end of namespace Opm

#include "CompWellInterface_impl.hpp"

#endif // OPM_COMP_WELLINTERFACE_HPP
