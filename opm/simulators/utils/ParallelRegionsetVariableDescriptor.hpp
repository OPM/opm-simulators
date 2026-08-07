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

#ifndef OPM_PARALLEL_REGIONSET_VARIABLE_DESCRIPTOR_HPP
#define OPM_PARALLEL_REGIONSET_VARIABLE_DESCRIPTOR_HPP

#include <opm/output/data/RegionsetVariableDescriptor.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <memory>

/// \file Component to describe a collection of region sets.

namespace Opm {

    /// Basic information about a collection of region sets.
    ///
    /// In particular, this class tracks the maximum region ID for each
    /// region set registered in the collection.  Common region sets in this
    /// context include the built-in FIPNUM set as well as user defined
    /// region sets named FIP*, such as FIPABC or FIP_UZ.  There is
    /// nevertheless nothing that will prevent including a different
    /// region set such as PVTNUM or KRNUMX in this collection.
    ///
    /// The primary use case for this class is assisting in allocating
    /// region level summary vectors over region sets.  That allocation
    /// requires knowing the maximum expected region ID for each pertinent
    /// region set.
    ///
    /// Constructing a descriptor object is a multi-step process.
    ///
    ///  -# Create the object through the default constructor
    ///  -# Call member function prepareDescriptorSet() to set the stage for
    ///     adding new region sets.
    ///  -# Incorporate one or more region sets through the addRegionSet()
    ///     overload set.
    ///  -# Call member function finaliseDescriptorSet() to analyse the data
    ///     and set up internal data structures.
    ///
    /// Once member function finaliseDescriptorSet() has been called, the
    /// object should be treated as read-only and only its \c const member
    /// functions should be invoked by client code.  Calling member function
    /// prepareDescriptorSet() after finaliseDescriptorSet() is allowed, but
    /// doing so will discard all information collected until that point and
    /// will require re-registering all pertinent region sets.
    ///
    /// This is the parallel implementation of the descriptor class.  It is
    /// intended to be used in conjunction with the parallel region variable
    /// values class, and synchronises the maximum region IDs between MPI ranks.
    class ParallelRegionsetVariableDescriptor : public data::RegionsetVariableDescriptor
    {
    public:
        /// Deleted default constructor.
        ParallelRegionsetVariableDescriptor() = delete;

        /// Constructor with MPI communication object.
        ///
        /// \param[in] comm MPI communication object.  Communication object
        /// must outlive the descriptor object.
        explicit ParallelRegionsetVariableDescriptor(const Parallel::Communication& comm);

        /// Destructor.
        virtual ~ParallelRegionsetVariableDescriptor() = default;

        /// Copy constructor.
        ParallelRegionsetVariableDescriptor(const ParallelRegionsetVariableDescriptor&) = default;

        /// Move constructor.
        ParallelRegionsetVariableDescriptor(ParallelRegionsetVariableDescriptor&&) = default;

        /// Assignment operator.
        ParallelRegionsetVariableDescriptor&
        operator=(const ParallelRegionsetVariableDescriptor&) = default;

        /// Move-assignment operator.
        ParallelRegionsetVariableDescriptor&
        operator=(ParallelRegionsetVariableDescriptor&&) = default;

        /// Create a deep copy of this object.
        ///
        /// Intended for when an object holds a pointer to a base class
        /// and needs to make a copy of the derived class object.
        ///
        /// \return A deep copy of this object.
        std::unique_ptr<data::RegionsetVariableDescriptor> clone() const override;

    protected:
        /// Exchange maximum region set IDs if needed.
        ///
        /// Parallel implementation that exchanges maximum region IDs across
        /// all ranks.
        void communicateGlobalRegsetMaxIDs() override;

    private:
        /// MPI communication object.
        Parallel::Communication comm_;
    };

} // namespace Opm

#endif // OPM_PARALLEL_REGIONSET_VARIABLE_DESCRIPTOR_HPP
