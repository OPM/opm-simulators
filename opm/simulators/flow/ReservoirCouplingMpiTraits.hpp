/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_RESERVOIR_COUPLING_MPI_TRAITS_HPP
#define OPM_RESERVOIR_COUPLING_MPI_TRAITS_HPP

#include <dune/common/parallel/mpitraits.hh>

#include <opm/simulators/flow/ReservoirCoupling.hpp>

#include <mpi.h>

namespace Dune {

template <>
struct MPITraits<::Opm::ReservoirCoupling::Potentials> {
    using Potentials = ::Opm::ReservoirCoupling::Potentials;
    constexpr static std::size_t num_fields = Potentials::num_fields;
public:
    static inline MPI_Datatype getType() {
        if (type == MPI_DATATYPE_NULL) { // Will be NULL only once (the first time)
            // Array of block lengths (number of elements per field)
            int block_lengths[num_fields] = {1, 1, 1};

            // Array of displacements for each field
            MPI_Aint displacements[num_fields];
            Potentials dummy;
            MPI_Aint base;
            MPI_Get_address(&dummy, &base);
            MPI_Get_address(&dummy.oil_rate, &displacements[Potentials::OIL_IDX]);
            MPI_Get_address(&dummy.gas_rate, &displacements[Potentials::GAS_IDX]);
            MPI_Get_address(&dummy.water_rate, &displacements[Potentials::WATER_IDX]);

            // Adjust displacements relative to the base
            for (std::size_t i = 0; i < num_fields; ++i) {
                displacements[i] -= base;
            }

            // Array of MPI data types for each field
            MPI_Datatype types[num_fields] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

            // Create the MPI datatype
            MPI_Datatype tmp_type;
            MPI_Type_create_struct(num_fields, block_lengths, displacements, types, &tmp_type);

            // Resize the datatype to account for possible padding issues
            MPI_Type_create_resized(tmp_type, 0, sizeof(Potentials), &type);

            MPI_Type_commit(&type);
            MPI_Type_free(&tmp_type);
        }
        return type;
    }

private:
    // Initial value of MPI_DATATYPE_NULL is used to indicate that the type
    //  has not been created yet
    static inline MPI_Datatype type = MPI_DATATYPE_NULL;
};

} // namespace Dune

#endif // OPM_RESERVOIR_COUPLING_MPI_TRAITS_HPP
