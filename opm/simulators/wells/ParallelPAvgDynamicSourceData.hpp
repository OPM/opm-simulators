/*
  Copyright 2023 Equinor ASA.

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

#ifndef PARALLEL_PAVG_DYNAMIC_SOURCE_DATA_HPP
#define PARALLEL_PAVG_DYNAMIC_SOURCE_DATA_HPP

#include <opm/input/eclipse/Schedule/Well/PAvgDynamicSourceData.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <cstddef>
#include <functional>
#include <vector>

namespace Opm {

/// Dynamic source data for block-average pressure calculations.
/// Specialisation for parallel runs.
class ParallelPAvgDynamicSourceData : public PAvgDynamicSourceData
{
public:
    /// Translate globally unique, linearised Cartesian cell indices to
    /// local, on-rank, cell indices.  Assumed to return a negative value
    /// result if the input cell index is not owned by the current rank.
    using GlobalToLocal = std::function<int(const std::size_t)>;

    /// Collect source term contributions from local, on-rank, cell.
    ///
    /// Called as
    /// \code
    ///    eval(cellIndex, sourceTerm)
    /// \endcode
    /// in which \c cellIndex is the local, on-rank, cell index in the range
    /// 0 to #active cells on rank - 1.  Function \c eval is expected to
    /// fill in/assign all \c sourceTerm items for this cell.
    using Evaluator = std::function<void(int, SourceDataSpan<double>)>;

    /// Constructor
    ///
    /// \param[in] comm MPI communication object.  Typically \code
    ///   grid.comm() \endcode from the main simulation grid.
    ///
    /// \param[in] sourceLocations Known locations, typically linearised
    ///   global call IDs, for which to enable collecting/reporting dynamic
    ///   source data.  Typically \code allWBPCells() \endcode from a \c
    ///   PAvgCalculatorCollection.
    ///
    /// \param[in] localCellIdx Translation from global, Cartesian cell
    ///   indices to local, on-rank, cell indices.
    ParallelPAvgDynamicSourceData(const Parallel::Communication&  comm,
                                  const std::vector<std::size_t>& sourceLocations,
                                  GlobalToLocal                   localCellIdx);

    /// Clear contents of local source term contributions.
    ///
    /// Mostly useful when collecting source term contributions along the
    /// well bore.
    void setToZero();

    /// Reconstruct Source Data backing storage and internal mapping tables
    ///
    /// Effectively replaces the original object formed by the constructor.
    /// Mainly intended for updating objects as new wells and/or new
    /// reservoir connections are introduced.
    ///
    /// \param[in] sourceLocations Known locations, typically linearised
    ///   global call IDs, for which to enable collecting/reporting dynamic
    ///   source data.  Typically \code allWBPCells() \endcode from a \c
    ///   PAvgCalculatorCollection.
    ///
    /// \param[in] localCellIdx Translation from global, Cartesian cell
    ///   indices to local, on-rank, cell indices.
    void reconstruct(const std::vector<std::size_t>& sourceLocations,
                     GlobalToLocal                   localCellIdx);

    /// Compute local, on-rank, contributions to the collection of source
    /// terms.
    ///
    /// \param[in] eval Source term evaluator object.
    void collectLocalSources(Evaluator eval);

    /// Exchange local contributions to build full, global view of all
    /// source terms.
    void synchroniseSources();

private:
    /// Identifier for local source term element.
    struct LocalLocation
    {
        /// Source term element index.
        std::size_t ix{};

        /// Local cell index for this source term (0..#active on rank - 1).
        int cell{};
    };

    /// MPI communication object.
    std::reference_wrapper<const Parallel::Communication> comm_;

    /// Subset of source locations owned by current rank.
    std::vector<LocalLocation> locations_{};

    /// Source data values owned by current rank.
    std::vector<double> localSrc_{};

    /// Translation map from element index to storage index in
    /// PAvgDynamicSourceData::src_.
    std::vector<std::vector<double>::size_type> storageIndex_{};

    /// Receive size from all ranks (allgatherv()).
    std::vector<int> allSizes_{}; // Type int to meet API requirements.

    /// Receive displacements for all ranks (allgatherv()).
    std::vector<int> startPointers_{}; // Type int to meet API requirements.

    /// Translate element index into storage index.
    ///
    /// Customisation point.
    ///
    /// \param[in] elemIndex Source element index.
    ///
    /// \return Storage (starting) index in PAvgDynamicSourceData::src_.
    [[nodiscard]] std::vector<double>::size_type
    storageIndex(std::vector<double>::size_type elemIndex) const override;

    /// Identify local source term elements on rank and build communication
    /// pattern for all source terms.
    ///
    /// Assigns \c storageIndex_, \c allSizes_, and \c startPointers_.
    ///
    /// \param[in] sourceLocations Known locations, typically linearised
    ///   global call IDs, for which to enable collecting/reporting dynamic
    ///   source data.  Typically \code allWBPCells() \endcode from a \c
    ///   PAvgCalculatorCollection.
    ///
    /// \param[in] localCellIdx Translation from global, Cartesian cell
    ///   indices to local, on-rank, cell indices.
    void finaliseConstruction(const std::vector<std::size_t>& sourceLocations,
                              GlobalToLocal                   localCellIdx);

    /// Form mutable data span into non-default backing store.
    ///
    /// \param[in] localIx Logical element index into \c localSrc_.
    ///
    /// \return Mutable view into \c localSrc_.
    [[nodiscard]] SourceDataSpan<double>
    localSourceTerm(const std::size_t localIx);

    /// Build communication pattern for all source terms.
    ///
    /// Assigns \c storageIndex_, \c allSizes_, and \c startPointers_.
    void defineCommunication();
};

} // namespace Opm

#endif // PARALLEL_PAVG_DYNAMIC_SOURCE_DATA_HPP
