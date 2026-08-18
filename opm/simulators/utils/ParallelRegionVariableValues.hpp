/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_PARALLEL_REGION_VARIABLE_VALUES_HPP
#define OPM_PARALLEL_REGION_VARIABLE_VALUES_HPP

#include <opm/output/data/RegionVariableValues.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <memory>

/// \file Component to collect per-region values of region level summary
/// quantities.

namespace Opm {

    /// Numerical values for set of region level summary variables defined
    /// over a collection of region sets.
    ///
    /// This is the parallel implementation of the region variable values
    /// class.  It is intended to be used in conjunction with the parallel
    /// region variable descriptor class, and synchronises the region value
    /// "increment" values between MPI ranks.
    class ParallelRegionVariableValues : public data::RegionVariableValues
    {
    public:
        /// \param[in] comm MPI communication object. Copied into the values object;
        /// the underlying communicator must remain valid for the values object's lifetime.
        explicit ParallelRegionVariableValues(const Parallel::Communication& comm);

        /// Virtual destructor to support inheritance.
        virtual ~ParallelRegionVariableValues() = default;

        /// Create a deep copy of this object.
        ///
        /// Intended for when an object holds a pointer to a base class
        /// and needs to make a copy of the derived class object.
        ///
        /// \return A deep copy of this object.
        std::unique_ptr<data::RegionVariableValues> clone() const override;

    protected:
        /// Compute canonical increment_ values.
        ///
        /// Sums values across MPI ranks.
        void communicateIncrement() override;

    private:
        /// MPI communication object.
        Parallel::Communication comm_;
    };

} // namespace Opm

#endif // OPM_PARALLEL_REGION_VARIABLE_VALUES_HPP
