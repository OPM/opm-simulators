// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2022 Equinor AS

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

#ifndef OPM_INTERREG_FLOWS_MODULE_HPP
#define OPM_INTERREG_FLOWS_MODULE_HPP

#include <opm/output/data/InterRegFlowMap.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

/// \file
///
/// MPI-aware facility for converting collection of tuples of region ID
/// pairs and associate flow rates into a sparse (CSR) adjacency matrix
/// representation of a graph.  Supports O(nnz) compression.

namespace Opm {

    class InterRegFlowMap;

    /// Form CSR adjacency matrix representation of inter-region flow rate
    /// graph provided as a list of connections between regions on local MPI
    /// rank.  Pertains to a single FIP definition array (e.g., FIPNUM).
    class InterRegFlowMapSingleFIP
    {
    public:
        /// Minimal characteristics of a cell from a simulation grid.
        struct Cell {
            /// Cell's active index on local rank.
            int activeIndex{-1};

            /// Cell's global cell ID.
            int cartesianIndex{-1};

            /// Whether or not cell is interior to local rank.
            bool isInterior{true};
        };

        friend class InterRegFlowMap;

        /// Constructor
        ///
        /// \param[in] region Local rank's FIP region definition array.
        explicit InterRegFlowMapSingleFIP(const std::vector<int>& region);

        /// Add flow rate connection between regions.
        ///
        /// \param[in] source Cell from which the flow nominally originates.
        ///
        /// \param[in] destination Cell into which flow nominally goes.
        ///
        /// \param[in] rates Flow rates associated to single connection.
        ///
        /// If both cells are in the same region, or if neither cell is
        /// interior to this MPI rank, then this function does nothing.  If
        /// one cell is interior to this MPI rank and the other isn't, then
        /// this function will include the flow rate contribution if and
        /// only if the cell with the smallest associate region ID is
        /// interior to this MPI rank.
        void addConnection(const Cell& source,
                           const Cell& destination,
                           const data::InterRegFlowMap::FlowRates& rates);

        /// Form CSR adjacency matrix representation of input graph from
        /// connections established in previous calls to addConnection().
        ///
        /// Number of rows in the CSR representation is the maximum FIP
        /// region ID.
        void compress();

        /// Clear all internal buffers, but preserve allocated capacity.
        void clear();

        /// Get read-only access to the underlying CSR representation.
        ///
        /// Mostly intended for summary output purposes.
        const data::InterRegFlowMap& getInterRegFlows() const;

        /// Retrieve maximum FIP region ID on local MPI rank.
        std::size_t getLocalMaxRegionID() const;

        /// Assign maximum FIP region ID across all MPI ranks.
        ///
        /// Fails if global maximum is smaller than local maximum region ID.
        ///
        /// \param[in] regID Global maximum FIP region ID for this FIP
        ///    definition array across all MPI ranks.
        ///
        /// \return Whether or not assignment succeeded.
        bool assignGlobalMaxRegionID(const std::size_t regID);

        /// Serialise internal representation to MPI message buffer
        ///
        /// \tparam MessageBufferType Linear MPI message buffer.  API should
        ///    be similar to Dune::MessageBufferIF
        ///
        /// \param[in,out] buffer Linear MPI message buffer instance.
        ///    Function appends a partially linearised representation of
        ///    \code *this \endcode to the buffer contents.
        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const
        {
            // Note: this->region_ is *intentionally* omitted here.
            buffer.write(this->maxGlobalRegionID_);
            this->iregFlow_.write(buffer);
        }

        /// Reconstitute internal object representation from MPI message
        /// buffer
        ///
        /// This object (\code *this \endcode) is not usable in subsequent
        /// calls to \code addConnection() \endcode following a call to
        /// member function \code read() \endcode.
        ///
        /// \tparam MessageBufferType Linear MPI message buffer.  API should
        ///    be similar to Dune::MessageBufferIF
        ///
        /// \param[in,out] buffer Linear MPI message buffer instance.
        ///    Function reads a partially linearised representation of \code
        ///    *this \endcode from the buffer contents and advances the
        ///    buffer's read position.
        template <class MessageBufferType>
        void read(MessageBufferType& buffer)
        {
            auto otherMaxRegionID = 0 * this->maxGlobalRegionID_;
            buffer.read(otherMaxRegionID);
            this->maxGlobalRegionID_ = std::max(this->maxGlobalRegionID_, otherMaxRegionID);

            this->iregFlow_.read(buffer);
            this->isReadFromStream_ = true;
        }

    private:
        /// Zero-based FIP region IDs on local MPI rank.
        std::vector<int> region_{};

        /// Maximum one-based FIP region ID on local MPI rank.
        std::size_t maxLocalRegionID_{0};

        /// Maximum one-based FIP region ID for this FIP region definition
        /// array across all MPI ranks.
        std::size_t maxGlobalRegionID_{0};

        /// Rank-local inter-regional flow map.
        data::InterRegFlowMap iregFlow_{};

        /// Whether or not this object contains contributions deserialised
        /// from a stream.  For error detection.
        bool isReadFromStream_{false};

        /// Default constructor.
        InterRegFlowMapSingleFIP() = default;
    };

    /// Inter-region flow accumulation maps for all region definition arrays
    class InterRegFlowMap
    {
    public:
        /// Minimal representation of a single named region defintion.
        ///
        /// Name is typically "FIPNUM" or any additional "FIP*" array names.
        struct SingleRegion {
            /// Region definition array name
            std::string name;

            /// Region definition array.
            std::reference_wrapper<const std::vector<int>> definition;
        };

        /// Characteristics of a cell from a simulation grid.
        using Cell = InterRegFlowMapSingleFIP::Cell;

        /// Default constructor.
        InterRegFlowMap() = default;

        /// Special purpose constructor for global object being collected on
        /// the I/O rank.
        ///
        /// Only knows about the FIP region set names.
        ///
        /// \param[in] names Sorted sequence of FIP region names.
        static InterRegFlowMap
        createMapFromNames(std::vector<std::string> names);

        /// Constructor.
        ///
        /// \param[in] numCells Number of cells on local MPI rank, including
        ///    overlap cells if applicable.
        ///
        /// \param[in] regions All applicable region definition arrays.
        ///
        /// \param[in] declaredMaxRegID Declared maximum region ID in the
        ///    run-typically from the TABDIMS and/or REGDIMS keywords.  Used
        ///    for sizing internal data structures if greater than zero.
        explicit InterRegFlowMap(const std::size_t                numCells,
                                 const std::vector<SingleRegion>& regions,
                                 const std::size_t                declaredMaxRegID = 0);

        InterRegFlowMap(const InterRegFlowMap& rhs) = default;
        InterRegFlowMap(InterRegFlowMap&& rhs) noexcept = default;

        InterRegFlowMap& operator=(const InterRegFlowMap& rhs) = default;
        InterRegFlowMap& operator=(InterRegFlowMap&& rhs) noexcept = default;

        /// Add flow rate connection between regions for all region
        /// definitions.
        ///
        /// \param[in] source Cell from which the flow nominally originates.
        ///
        /// \param[in] destination Cell into which flow nominally goes.
        ///
        /// \param[in] rates Flow rates associated to single connection.
        ///
        /// If both cells are in the same region, or if neither cell is
        /// interior to this MPI rank, then this function does nothing.  If
        /// one cell is interior to this MPI rank and the other isn't, then
        /// this function will include the flow rate contribution if and
        /// only if the cell with the smallest associate region ID is
        /// interior to this MPI rank.
        void addConnection(const Cell& source,
                           const Cell& destination,
                           const data::InterRegFlowMap::FlowRates& rates);

        /// Form CSR adjacency matrix representation of input graph from
        /// connections established in previous calls to addConnection().
        ///
        /// Number of rows in the CSR representation is the maximum FIP
        /// region ID.
        void compress();

        /// Clear all internal buffers, but preserve allocated capacity.
        void clear();

        /// Names of all applicable region definition arrays.
        ///
        /// Mostly intended for summary output purposes.
        const std::vector<std::string>& names() const;

        /// Get read-only access to the underlying CSR representation.
        ///
        /// Mostly intended for summary output purposes.
        std::vector<data::InterRegFlowMap> getInterRegFlows() const;

        /// Retrieve maximum FIP region ID on local MPI rank.
        std::vector<std::size_t> getLocalMaxRegionID() const;

        /// Assign maximum FIP region ID across all MPI ranks.
        ///
        /// Fails if global maximum is smaller than local maximum region ID.
        ///
        /// \param[in] regID Global maximum FIP region ID for this FIP
        ///    definition array across all MPI ranks.
        ///
        /// \return Whether or not assignment succeeded.
        bool assignGlobalMaxRegionID(const std::vector<std::size_t>& regID);

        /// Whether or not previous read() operation succeeded.
        bool readIsConsistent() const;

        bool wantInterRegflowSummary() const
        {
            return !this->regionMaps_.empty();
        }

        /// Serialise internal representation to MPI message buffer
        ///
        /// \tparam MessageBufferType Linear MPI message buffer.  API should
        ///    be similar to Dune::MessageBufferIF
        ///
        /// \param[in,out] buffer Linear MPI message buffer instance.
        ///    Function appends a partially linearised representation of
        ///    \code *this \endcode to the buffer contents.
        template <class MessageBufferType>
        void write(MessageBufferType& buffer) const
        {
            buffer.write(this->names_.size());
            for (const auto& name : this->names_) {
                buffer.write(name);
            }

            for (const auto& map : this->regionMaps_) {
                map.write(buffer);
            }
        }

        /// Reconstitute internal object representation from MPI message
        /// buffer
        ///
        /// This object (\code *this \endcode) is not usable in subsequent
        /// calls to \code addConnection() \endcode following a call to
        /// member function \code read() \endcode.
        ///
        /// \tparam MessageBufferType Linear MPI message buffer.  API should
        ///    be similar to Dune::MessageBufferIF
        ///
        /// \param[in,out] buffer Linear MPI message buffer instance.
        ///    Function reads a partially linearised representation of \code
        ///    *this \endcode from the buffer contents and advances the
        ///    buffer's read position.
        template <class MessageBufferType>
        void read(MessageBufferType& buffer)
        {
            const auto names = this->readNames(buffer);

            if (names == this->names_) {
                // Input stream holds the maps in expected order (common
                // case).  Read the maps and merge with internal values.
                for (auto& map : this->regionMaps_) {
                    map.read(buffer);
                }
            }
            else {
                // Input stream does not hold the maps in expected order (or
                // different number of maps).  Unexpected.  Read the values
                // from the input stream, but do not merge with internal
                // values.
                auto map = InterRegFlowMapSingleFIP {
                    std::vector<int>(this->numCells_, 1)
                };

                const auto numMaps = names.size();
                for (auto n = 0*numMaps; n < numMaps; ++n) {
                    map.read(buffer);
                }

                this->readIsConsistent_ = false;
            }
        }

    private:
        /// Inter-region flow accumulators.  One accumulator map for each
        /// region definition array.
        std::vector<InterRegFlowMapSingleFIP> regionMaps_{};

        /// Names of region definition arrays.  Typically "FIPNUM" and other
        /// "FIPXYZ" array names.
        std::vector<std::string> names_;

        /// Number of cells, including overlap, reachable from this MPI
        /// rank.
        std::size_t numCells_{0};

        /// Wheter or not read() successfully updated this object from an
        /// input stream.
        bool readIsConsistent_{true};

        /// Retrieve array names from an input stream.
        ///
        /// Needed for consistency checks.
        template <class MessageBufferType>
        std::vector<std::string> readNames(MessageBufferType& buffer) const
        {
            auto numNames = 0 * this->names_.size();
            buffer.read(numNames);

            auto names = std::vector<std::string>(numNames);
            for (auto name = 0*numNames; name < numNames; ++name) {
                buffer.read(names[name]);
            }

            return names;
        }
    };
} // namespace Opm

#endif // OPM_INTERREG_FLOWS_MODULE_HPP
