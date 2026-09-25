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

#include <string>

namespace Opm {
template <typename TypeTag>
CompWellInterface<TypeTag>::
CompWellInterface(const Well& well,
                  const int index_of_well,
                  const std::vector<CompConnectionData>& well_connection_data)
    : well_ecl_(well)
    , index_of_well_(index_of_well)
    , number_of_connection_(well_connection_data.size())
    , reference_depth_(well.getRefDepth())
    , connectionRates_(number_of_connection_)
{
    {
        well_cells_.resize(number_of_connection_);
        well_index_.resize(number_of_connection_);
        saturation_table_number_.resize(number_of_connection_, 0);
        int connection_idx = 0;
        for (const auto& connection_data : well_connection_data) {
            well_cells_[connection_idx] = connection_data.cell_index;
            well_index_[connection_idx] = connection_data.connection_transmissibility_factor;
            saturation_table_number_[connection_idx] = connection_data.satnum_id;
            ++connection_idx;
        }
        // TODO: saturation_table_number
    }

}

template <typename TypeTag>
void
CompWellInterface<TypeTag>::
init()
{
    // more things to add here
}

template <typename TypeTag>
const std::string&
CompWellInterface<TypeTag>::
name() const
{
    return this->well_ecl_.name();
}

template <typename TypeTag>
void
CompWellInterface<TypeTag>::
addCellRates(RateVector& rates, unsigned cellIdx) const
{
    for (int con = 0; con < this->number_of_connection_; ++con) {
        if (this->well_cells_[con] == cellIdx) {
            for (int i = 0; i < RateVector::dimension; ++i) {
                rates[i] += connectionRates_[con][i];
            }
        }
    }
}


} // end of namespace Opm
