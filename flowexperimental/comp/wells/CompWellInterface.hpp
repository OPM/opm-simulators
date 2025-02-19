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

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <string>

#include "SingleCompWellState.hpp"

namespace Opm
{
template <typename Scalar>
class CompConnectionData;

template <typename TypeTag> // TODO: do we need to use TypeTag here?
class CompWellInterface
{
public:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>; // FlashRateVector at the moment

    CompWellInterface(const Well& well,
                      const int index_of_well,
                      const std::vector<CompConnectionData<Scalar>>& well_connection_data);

    const std::string& name() const;

    virtual void calculateExplicitQuantities(const Simulator& simulator,
                                             const SingleCompWellState<Scalar>& well_state) = 0;

    virtual void updatePrimaryVariables(const Simulator& simulator,
                                        const SingleCompWellState<Scalar>& well_state) = 0;

    // TODO: not sure this funciton will be used
    // but leaviing it here for protoyping purpsoe
    void solveWellEq(const Simulator& simulator,
                     SingleCompWellState<Scalar>& well_state);

    virtual bool iterateWellEq(const Simulator& simulator,
                               const Scalar dt,
                               SingleCompWellState<Scalar>& well_state) = 0;

    void addCellRates(RateVector& rates, unsigned cellIdx) const;

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
    virtual void init();

};


} // end of namespace Opm

#include "CompWellInterface_impl.hpp"

#endif // OPM_COMP_WELLINTERFACE_HPP
