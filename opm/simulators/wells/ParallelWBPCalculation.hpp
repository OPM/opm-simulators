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

#ifndef PARALLEL_WBP_CALCULATION_HPP
#define PARALLEL_WBP_CALCULATION_HPP

#include <opm/simulators/wells/ParallelPAvgDynamicSourceData.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/simulators/wells/PerforationData.hpp>

#include <opm/input/eclipse/Schedule/Well/PAvgCalculator.hpp>
#include <opm/input/eclipse/Schedule/Well/PAvgCalculatorCollection.hpp>

#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

namespace Opm {
    class GridDims;
    class ParallelWellInfo;
    class PAvg;
    class Well;
}

namespace Opm {

/// Parallel facility for managing the on-rank collection and global
/// distribution of WBPn source values as well as local calculation and
/// distributed reduction of the inferred WBPn report values.
class ParallelWBPCalculation
{
public:
    /// Callback for inferring the source locations which are active on the
    /// current MPI rank.
    using GlobalToLocal = ParallelPAvgDynamicSourceData::GlobalToLocal;

    /// Callback for evaluating WBPn source terms on the current MPI rank.
    using Evaluator = ParallelPAvgDynamicSourceData::Evaluator;

    /// Callback for constructing a source term evaluation function on the
    /// current MPI rank.  Needed for deferred construction of per-well
    /// source term evaluation functions.
    using EvaluatorFactory = std::function<Evaluator()>;

    /// Constructor.
    ///
    /// \param[in] cellIndexMap Cell index triple map ((I,J,K) <-> global).
    ///
    /// \param[in] gridComm Main, grid level, global communicator.
    explicit ParallelWBPCalculation(const GridDims&                cellIndexMap,
                                    const Parallel::Communication& gridComm);

    /// Assign translation function for inferring the on-rank IDs of the
    /// known source locations.
    ///
    /// \param[in] localCellIdx Translation from global, Cartesian cell
    ///   indices to local, on-rank, cell indices.
    ///
    /// \return \code *this \endcode
    ParallelWBPCalculation& localCellIndex(GlobalToLocal localCellIdx);

    /// Assign evaluation function for computing the on-rank, cell level
    /// WBPn source terms.
    ///
    /// \param[in] evalCellSrc Source term evaluation function.
    ///
    /// \return \code *this \endcode
    ParallelWBPCalculation& evalCellSource(Evaluator evalCellSrc);

    /// Create, or reassign, a WBPn calculation object for a particular
    /// well.
    ///
    /// \param[in] well Well for which to create a WBPn calculation object.
    ///
    /// \param[in] parallelWellInfo Communicator object for the ranks
    ///   sharing this well.
    ///
    /// \param[in] localConnIdx Local (on-rank) connection index.  Sized
    ///   according to \code well.getConnections().size() \endcode, but with
    ///   a non-negative entry only for those active connections that
    ///   intersect the current rank.  If \code localConnIdx[i] == j
    ///   \endcode, then the \c i-th global connection is the \c j-th active
    ///   connection on the current rank.  Use a negative value to identify
    ///   a connection that is either not flowing or which does not
    ///   intersect the current MPI rank.
    ///
    /// \param[in] makeWellSourceEvaluator Factory function to support
    ///   deferred creation of an evaluation function for the per-connection
    ///   WBP source terms.
    ///
    /// \return Calculator object index.  Must be used in subsequent calls
    ///   to inferBlockAveragePressures() and averagePressures() for this \p
    ///   well.
    std::size_t
    createCalculator(const Well&             well,
                     const ParallelWellInfo& parallelWellInfo,
                     const std::vector<int>& localConnIdx,
                     EvaluatorFactory        makeWellSourceEvaluator);

    /// Set up communication patterns for both cell and connection level
    /// source terms and partial/intermediate WBPn results.
    ///
    /// Clients must call this function once all calculation objects have
    /// been created, and strictly before the first call to member function
    /// \c collectDynamicValues().
    void defineCommunication();

    /// Collect all on-rank source term value and distribute those on-rank
    /// values to all other MPI ranks.
    ///
    /// Will call the registered source term evaluation functions for all
    /// on-rank source locations.  Once this function returns, all ranks
    /// have a full view of all cell level source term values, and all ranks
    /// which share an individual well have a full view of the
    /// per-connection source term values for that well.
    ///
    /// Will \c throw an object of type \code std::logic_error \endcode if
    /// the communication patterns have not been fully defined through a
    /// prior call to member function \c defineCommunication().
    void collectDynamicValues();

    /// Compute WBPn report values for a single well.
    ///
    /// \param[in] calcIndex Calculator object index.  Return value from a
    ///   previous call to member function \c createCalculator().
    ///
    /// \param[in] controls Pressure averaging procedure controls for this
    ///   well.
    ///
    /// \param[in] gravity Strength of gravity acceleration.
    ///
    /// \param[in] refDepth WBPn reference depth.  Typically \code
    ///   Well::getWPaveRefDepth() \endcode.
    void inferBlockAveragePressures(const std::size_t calcIndex,
                                    const PAvg&       controls,
                                    const double      gravity,
                                    const double      refDepth);

    /// Retrieve results from most recent WBPn value calculation for
    /// specified well.
    ///
    /// \param[in] calcIndex Calculator object index.  Return value from a
    ///   previous call to member function \c createCalculator().
    ///
    /// \return Result set from most recent call to member function \c
    ///   inferBlockAveragePressures() for \c calcIndex.
    const PAvgCalculator::Result&
    averagePressures(const std::size_t calcIndex) const;

private:
    /// Callable wrapper for the local, per-well reservoir connections.
    class LocalConnSet
    {
    public:
        /// Constructor.
        ///
        /// \param[in] localConnIdx Local (on-rank) connection index.  Sized
        ///   according to total number of connections of a well, but with a
        ///   non-negative entry only for those active connections that
        ///   intersect the current rank.  If \code localConnIdx[i] == j
        ///   \endcode, then the \c i-th global connection is the \c j-th
        ///   active connection on the current rank.  Use a negative value
        ///   to identify a connection that is either not flowing or which
        ///   does not intersect the current MPI rank.
        explicit LocalConnSet(const std::vector<int>& localConnIdx);

        /// Retrieve local, on-rank connection index of a well's global
        /// connection index.
        ///
        /// \param[in] connIdx Well's global connection index
        ///
        /// \return Local, on-rank, connection index of well's global index.
        ///   Negative value if \p connIdx does not intersect the current
        ///   rank.
        int localIndex(const std::size_t connIdx) const;

    private:
        /// Local (on-rank) connection index.  Sized according to total
        /// number of connections of a well, but with a non-negative entry
        /// only for those active connections that intersect the current
        /// rank.  If \code localConnIdx[i] == j \endcode, then the \c i-th
        /// global connection is the \c j-th active connection on the
        /// current rank.  Use a negative value to identify a connection
        /// that is either not flowing or which does not intersect the
        /// current MPI rank.
        std::vector<int> localConnIdx_{};
    };

    /// Parallel collection of individual source terms
    class SourceData
    {
    public:
        /// Conversion operator.
        ///
        /// Enables transparent usage of this object as an argument to \code
        /// PAvgCalculator::inferBlockAveragePressures() \endcode.
        operator const ParallelPAvgDynamicSourceData&() const
        {
            return *this->srcData_;
        }

        /// Constructor.
        ///
        /// \param[in] comm MPI communicator.
        explicit SourceData(const Parallel::Communication& comm);

        /// Assign translation function for inferring the on-rank IDs of the
        /// known source locations.
        ///
        /// \param[in] localIdx Translation from global indices to local,
        ///   on-rank indices.
        ///
        /// \return \code *this \endcode
        SourceData& localIndex(GlobalToLocal localIdx);

        /// Assign evaluation function for computing the on-rank, WBPn
        /// source terms.
        ///
        /// \param[in] eval Source term evaluation function.
        ///
        /// \return \code *this \endcode
        SourceData& evaluator(Evaluator eval);

        /// Assign evaluation function factory for computing the on-rank,
        /// WBPn source terms.
        ///
        /// \param[in] evalFactory Factory function to support deferred
        ///   creation of an evaluation function for the per-connection WBP
        ///   source terms.
        ///
        /// \return \code *this \endcode
        SourceData& evaluatorFactory(EvaluatorFactory evalFactory);

        /// Construct or update wrapped source term collection to account
        /// for new set of source locations.
        ///
        /// This function is a prerequisite for subsequent calls to \c
        /// collectDynamicValues().
        ///
        /// \param[in] sourceLocations Known locations, possibly linearised
        ///   global cell IDs, for which to enable collecting/reporting
        ///   dynamic source data.
        void buildStructure(const std::vector<std::size_t>& sourceLocations);

        /// Collect all on-rank source term value and distribute those
        /// on-rank values to all other MPI ranks.
        ///
        /// Will call the registered source term evaluation functions for
        /// all on-rank source locations.
        ///
        /// Will \c throw an object of type \code std::logic_error \endcode
        /// if the communication pattern has not been fully defined through
        /// a prior call to member function \c buildStructure().
        void collectDynamicValues();

        /// Provide read-only access to the associate MPI communicator.
        ///
        /// \return Read-only MPI communicator.
        const Parallel::Communication& comm() const
        {
            return this->comm_.get();
        }

        /// Convert collection of global indices to local, on-rank, indices.
        ///
        /// Applies \c localIndex() to each element of the input array.
        ///
        /// \param[in] globalIndex Collection of global indices.
        ///
        /// \return Local, on-rank, indices for each index in \p
        ///   globalIndex.  Negative local index value for global indices
        ///   which are either not active or not on the current MPI rank.
        std::vector<int>
        getLocalIndex(const std::vector<std::size_t>& globalIndex) const;

    private:
        /// Type of wrapped object.
        using DataPtr = std::unique_ptr<ParallelPAvgDynamicSourceData>;

        /// MPI communicator for this source data object.
        std::reference_wrapper<const Parallel::Communication> comm_;

        /// Translation from global indices to local, on-rank, indices.
        GlobalToLocal localIdx_{};

        /// Source term evaluator object.  Empty if we need deferred
        /// initialisation, e.g., for well connections.
        Evaluator eval_{};

        /// Creation function for source term evaluation functions.  Empty
        /// if evaluation function has already been assigned, e.g., for
        /// reservoir cells.
        EvaluatorFactory evalFactory_{};

        /// Parallel WBPn source term object.
        DataPtr srcData_{};
    };

    /// Calculation object IDs.
    using WellID = std::vector<SourceData>::size_type;

    /// Cell index triple map ((I,J,K) <-> global).
    std::reference_wrapper<const GridDims> cellIndexMap_;

    /// Source term object for the reservoir cells.
    SourceData reservoirSrc_;

    /// Collection of WBPn calculation objects.  One object for each well on
    /// rank.
    PAvgCalculatorCollection calculators_{};

    /// Source term objects for each well on rank.
    std::vector<SourceData> wellConnSrc_{};

    /// Local connection indices for each well on rank.
    std::vector<LocalConnSet> localConnSet_{};

    /// Eliminate inactive cells from the source locations backing \c
    /// reservoirSrc_.
    void pruneInactiveWBPCells();

    /// Serial implementation of pruneInactiveWBPCells().
    void pruneInactiveWBPCellsSerial();

    /// Parallel implementation of pruneInactiveWBPCells().
    void pruneInactiveWBPCellsParallel();

    /// Define communication patterns for \c reservoirSrc_.
    void defineReservoirCommunication();

    /// Define communication patterns for source terms pertaining to the
    /// reservoir connections of a single well.
    ///
    /// \param[in] well Well for which to define communication pattern.
    void defineWellCommunication(const std::size_t well);

    /// Aggregate pertinent source terms for the WBPn calculation object of
    /// a single well.
    ///
    /// \param[in] well Well for which to aggregate the pertient source
    /// terms.
    ///
    /// \return WBPn source terms aggregated for \p well.
    PAvgCalculator::Sources makeEvaluationSources(const WellID well) const;
};

} // namespace Opm

#endif // PARALLEL_WBP_CALCULATION_HPP
