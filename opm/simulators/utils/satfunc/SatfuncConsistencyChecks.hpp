/*
  Copyright 2024 Equinor AS

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

#ifndef OPM_SATFUNC_CONSISTENCY_CHECK_MODULE_HPP
#define OPM_SATFUNC_CONSISTENCY_CHECK_MODULE_HPP

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <cstddef>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <string_view>
#include <vector>

namespace Opm {
    template <typename Scalar>
    struct EclEpsScalingPointsInfo;
} // namespace Opm

namespace Opm {

    /// Platform for running sets of consistency checks against collection
    /// of saturation function end-points
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class SatfuncConsistencyChecks
    {
    public:
        /// Call-back interface for an individual check.
        ///
        /// Specific checks are expected to inherit from this base class.
        class Check
        {
        public:
            /// Run specific check against a set of saturation function end-points.
            ///
            /// \param[in] endPoints Set of saturation function end-points.
            ///    Might for instance be the scaled end-points of the
            ///    drainage functions in a single grid block or the unscaled
            ///    end-points of the tabulated saturation functions in a
            ///    single saturation region.
            virtual void test(const EclEpsScalingPointsInfo<Scalar>& endPoints) = 0;

            /// Whether or not last set of end-points violated this particular check.
            virtual bool isViolated() const = 0;

            /// Whether or not this check is critical to the simulator's
            /// ability to run the case.
            ///
            /// Violating critical checks should typically stop the run.
            virtual bool isCritical() const = 0;

            /// Number of \c Scalar values involved in the check.
            virtual std::size_t numExportedCheckValues() const = 0;

            /// Get a linearised copy of the \c Sclar values involved in the check.
            ///
            /// \param[in,out] exportedCheckValues Pointer to contiguous
            ///    sequence of at least numExportedCheckValues() \c Scalars.
            ///    It is the responsibility of exportCheckValues() to
            ///    populate this sequence with sensible values.
            virtual void exportCheckValues(Scalar* exportedCheckValues) const = 0;

            /// Descriptive textual summary of this check.
            ///
            /// Might for instance be something along the lines of
            ///    Oil-phase scaled end-points
            virtual std::string description() const = 0;

            /// Textual representation of the consistency condition.
            virtual std::string condition() const = 0;

            /// Retrieve names of the exported check values.
            ///
            /// Order should match that of exportCheckValues().
            ///
            /// \param[in,out] headers Pointer to contiguous sequence of at
            ///    least numExportedCheckValues() strings.  It is the
            ///    responsibility of columnNames() to populate this sequence
            ///    with sensible values.
            virtual void columnNames(std::string* headers) const = 0;
        };

        /// Severity level for consistency condition violation.
        enum class ViolationLevel : std::size_t {
            /// Consistency condition violated, but we're able to continue
            /// the run.
            Standard,

            /// Consistency condition violated and we're not able to
            /// continue the run.
            Critical,

            /// Implementation helper.  Must be last enumerator.
            NumLevels,
        };

        /// Call-back function type for outputting a single record of a
        /// consistency condition violation report.
        using ReportRecordOutput = std::function<void(std::string_view)>;

        /// Call-back function type for formatting a numeric end-point ID.
        using PointIDFormatCallback = std::function<std::string(std::size_t)>;

        /// Constructor
        ///
        /// \param[in] pointName Name/category of the points in this set of
        ///    checks.  Might for instance be "Grid block" or "Saturation
        ///    region".  Will be used as a column header.
        ///
        /// \param[in] numSamplePoints Upper bound on the number of
        ///    end-point check violations to preserve for reporting
        ///    purposes.  Should normally be a small number like 5 or 10.
        explicit SatfuncConsistencyChecks(std::string_view  pointName,
                                          const std::size_t numSamplePoints);

        /// Destructor.
        ~SatfuncConsistencyChecks() = default;

        /// Deleted copy constructor
        SatfuncConsistencyChecks(const SatfuncConsistencyChecks& rhs) = delete;

        /// Move-constructor.
        ///
        /// \param[in,out] rhs Source object.  Left in a "valid but
        ///    unspecified" state on exit.
        SatfuncConsistencyChecks(SatfuncConsistencyChecks&& rhs);

        /// Deleted assignment operator.
        SatfuncConsistencyChecks&
        operator=(const SatfuncConsistencyChecks& rhs) = delete;

        /// Move-assignment operator.
        ///
        /// \param[in,out] rhs Source object.  Left in a "valid but
        ///    unspecified" state on exit.
        ///
        /// \return \code *this \endcode.
        SatfuncConsistencyChecks& operator=(SatfuncConsistencyChecks&& rhs);

        /// Replace formatting function for end-point IDs.
        ///
        /// The default formatting function is just the identity (\code
        /// std::to_string() \endcode) which is useful for testing, but
        /// which will for instance not capture Cartesian structure.
        ///
        /// \param[in] formatPointID Call-back function type for formatting
        ///    a numeric end-point ID.
        ///
        /// \return \code *this \endcode
        SatfuncConsistencyChecks&
        setPointIDFormatCallback(const PointIDFormatCallback& formatPointID)
        {
            this->formatPointID_ = formatPointID;
            return *this;
        }

        /// Clear current set of end-point checks.
        void resetCheckSet();

        /// Add specific check to in-progress check set.
        ///
        /// \param[in] check Particular end-point check.
        void addCheck(std::unique_ptr<Check> check);

        /// Commit current set of checks and build requisite internal
        /// support structures.
        void finaliseCheckSet();

        /// Run current set of checks against a specific set of end-points.
        ///
        /// \param[in] pointID Numeric identifier for this particular set of
        ///    end-points.  Typically a saturation region or a cell ID.
        ///
        /// \param[in] endPoints Set of saturation function end-points.
        ///    Might for instance be the scaled end-points of the drainage
        ///    functions in a single grid block or the unscaled end-points
        ///    of the tabulated saturation functions in a single saturation
        ///    region.  Will be passed directly on to \code Check::test()
        ///    \endcode for each check in the current set.
        void checkEndpoints(const std::size_t                      pointID,
                            const EclEpsScalingPointsInfo<Scalar>& endPoints);

        /// Collect consistency violations from all ranks in MPI communicator.
        ///
        /// Incorporates violation counts and sampled failure points into
        /// the internal structures on each rank.  Aggregate results useful
        /// for subsequent call to reportFailures() on root process.
        ///
        /// \param[in] root MPI root process.  This is the process onto
        ///    which the counts and samples will be collected.  Typically
        ///    the index of the IO rank.
        ///
        /// \param[in] comm MPI communication object.
        void collectFailures(int root, const Parallel::Communication& comm);

        /// Whether or not any checks failed at the \c Standard level.
        bool anyFailedChecks() const;

        /// Whether or not any checks failed at the \c Critical level.
        bool anyFailedCriticalChecks() const;

        /// Generate textual summary output of all failed consistency checks
        /// at specific level.
        ///
        /// Reports only those conditions/checks for which there is at least
        /// one violation.
        ///
        /// In a parallel run it is only safe to call this function on the
        /// MPI process to which the consistency check violations were
        /// collected in a previous call to collectFailures().
        ///
        /// \param[in] level Report's severity level.
        ///
        /// \param[in] emitReportRecord Call-back function for outputting a
        ///    single record/line of a violation report.  Typically a
        ///    wrapper of \code OpmLog::info() \endcode.  It is the
        ///    responsibility of emitReportRecord() to properly display the
        ///    text lines to end users.
        void reportFailures(const ViolationLevel      level,
                            const ReportRecordOutput& emitReportRecord) const;

    private:
        /// Uniform random bit generator for "reservoir sampling".
        ///
        /// We don't need mt19937 for this problem.  The linear congruential
        /// engine has a smaller internal state than mt19937 and that aspect
        /// is more important here.
        using RandomBitGenerator = std::minstd_rand;

        /// Sample of consistency check violations at single severity level.
        struct ViolationSample
        {
            /// Number of consistency check violations.
            ///
            /// Size equal to number of consistency checks.
            std::vector<std::size_t> count{};

            /// Sample of point IDs for violated consistency checks.
            ///
            /// \c numSamplePoints_ allocated for each consistency check.
            /// Number of valid entries for check \c i is minimum of \c
            /// numSamplePoints_ and \code count[i] \endcode.
            std::vector<std::size_t> pointID{};

            /// Scalar values for each sampled point.
            ///
            /// \c numSamplePoints_ allocated for each individual check, and
            /// the number of values per check determined by \code
            /// Check::numExportedCheckValues() \endcode.  Number of valid
            /// entries for check \c i is minimum of \c numSamplePoints_ and
            /// \code count[i] \endcode.
            std::vector<Scalar> checkValues{};

            /// Clear contents of all data members.
            void clear();
        };

        /// Collection of consistency check violations.
        ///
        /// One set of violations for each severity level.
        using ViolationCollection = std::array
            <ViolationSample, static_cast<std::size_t>(ViolationLevel::NumLevels)>;

        /// Name/category of the points in this set of checks.
        ///
        /// Set by constructor.  Might for instance be "Grid block" or
        /// "Saturation region".  Will be used as a column header.
        std::string pointName_{};

        /// Maximum number of points retained for reporting purposes.
        std::size_t numSamplePoints_;

        /// Formatting function for end-point IDs.
        ///
        /// The default formatting function is just the identity (\code
        /// std::to_string() \endcode) which is useful for testing, but
        /// which will for instance not capture Cartesian structure.
        PointIDFormatCallback formatPointID_{};

        /// Start offsets into \code ViolationSample::checkValues \endcode
        /// for each individual check.
        std::vector<typename std::vector<Scalar>::size_type> startCheckValues_{};

        /// Collection of sampled consistency check violations.
        ///
        /// One collection for each severity level.
        ViolationCollection violations_{};

        /// Collection of checks to run against the saturation function
        /// end-points.
        std::vector<std::unique_ptr<Check>> battery_{};

        /// Random bit generator for point sampling.
        ///
        /// Represented as a pointer in order to avoid allocation and
        /// initialisation when the facility is not needed.  This is useful
        /// when none of the consistency checks are violated which in turn
        /// is a common case in production runs.
        std::unique_ptr<RandomBitGenerator> urbg_{};

        /// Collect violations of single severity level from all ranks in
        /// MPI communicator.
        ///
        /// Incorporates violation counts and sampled failure points into
        /// the internal structures on each rank.  Aggregate results useful
        /// for subsequent call to reportFailures().
        ///
        /// \param[in] root MPI root process.  This is the process/rank onto
        ///    which the counts and samples will be collected.  Typically
        ///    the index of the IO rank.
        ///
        /// \param[in] comm MPI communication object.
        ///
        /// \param[in, out] violation Current rank's violation structure for
        ///    a single severity level.  Holds aggregate values across all
        ///    ranks, including updated sample points, on return.
        void collectFailures(int                            root,
                             const Parallel::Communication& comm,
                             ViolationSample&               violation);

        /// Allocate and initialise backing storage for a single set of
        /// sampled consistency check violations.
        ///
        /// \param[out] violation On return, zeroed or otherwise initialised
        ///    violation sample of proper size.
        void buildStructure(ViolationSample& violation);

        /// Internalise a single violation into internal data structures.
        ///
        /// Counts the violation and uses "reservoir sampling"
        /// (https://en.wikipedia.org/wiki/Reservoir_sampling) to determine
        /// whether or not to include the specific point into the reporting
        /// sample.
        ///
        /// \tparam PopulateCheckValues Call-back function type
        /// encapsulating block of code populate sequence of check values
        /// for a single, failed consistency check.  Expected to be a
        /// callable type with a function call operator of the form
        /// \code
        ///    void operator()(Scalar* checkValues) const
        /// \endcode
        /// in which the \c checkValues points the start of a sequence of
        /// values associated to particular check.  The call-back function
        /// is expected to know how many values are in a valid sequence and
        /// to fill in exactly this many values.
        ///
        /// \param[in, out] violation Current rank's violation sample at
        ///    particular severity level.
        ///
        /// \param[in] checkIx Numerical check index in the range
        ///    [0..battery_.size()).
        ///
        /// \param[in] pointID Numeric identifier for this particular set of
        ///    end-points.  Typically a saturation region or a cell ID.
        ///
        /// \param[in] populateCheckValues Call-back function to populate a
        ///    sequence of values pertaining to specified check.  Typically
        ///    \code Check::exportCheckValues() \endcode or a copy routine
        ///    to incorporate samples from multiple MPI ranks.
        template <typename PopulateCheckValues>
        void processViolation(ViolationSample&      violation,
                              const std::size_t     checkIx,
                              const std::size_t     pointID,
                              PopulateCheckValues&& populateCheckValues);

        /// Internalise a single violation into internal data structures.
        ///
        /// Counts the violation and uses "reservoir sampling"
        /// (https://en.wikipedia.org/wiki/Reservoir_sampling) to determine
        /// whether or not to include the specific point into the reporting
        /// sample.
        ///
        /// \param[in] level Violation severity level.
        ///
        /// \param[in] checkIx Numerical check index in the range
        ///    [0..battery_.size()).
        ///
        /// \param[in] pointID Numeric identifier for this particular set of
        ///    end-points.  Typically a saturation region or a cell ID.
        void processViolation(const ViolationLevel level,
                              const std::size_t    checkIx,
                              const std::size_t    pointID);

        /// Incorporate single severity level's set of violations from
        /// single MPI rank into current rank's internal data structures.
        ///
        /// \param[in] count Start of sequence of failure counts for all
        ///    checks from single MPI rank.
        ///
        /// \param[in] pointID Start of sequence of sampled point IDs for
        ///    all checks from a single MPI rank.
        ///
        /// \param[in] checkValues Start of sequence of sampled check values
        ///    for all checks from a single MPI rank.
        ///
        /// \param[in, out] violation
        void incorporateRankViolations(const std::size_t* count,
                                       const std::size_t* pointID,
                                       const Scalar*      checkValues,
                                       ViolationSample&   violation);

        /// Generate random index in the sample size.
        ///
        /// \param[in] sampleSize Total number of violations of a particular
        ///    check.
        ///
        /// \return Uniformly distributed random integer in the range
        ///    0..sampleSize-1, inclusive.
        std::size_t getSampleIndex(const std::size_t sampleSize);

        /// Ensure that random bit generator exists and is properly
        /// initialised.
        void ensureRandomBitGeneratorIsInitialised();

        /// Compute start offset into ViolationSample::pointID for
        /// particular check.
        ///
        /// \param[in] checkIx Numerical check index in the range
        ///    [0..battery_.size()).
        ///
        /// \return Start offset into ViolationSample::pointID.
        std::vector<std::size_t>::size_type
        violationPointIDStart(const std::size_t checkIx) const;

        /// Compute start offset into ViolationSample::checkValues for
        /// particular check and sample index.
        ///
        /// \param[in] checkIx Numerical check index in the range
        ///    [0..battery_.size()).
        ///
        /// \param[in] sampleIx Numerical index in the range
        ///    [0..min(count[checkIx], numSamplePoints_)).
        ///
        /// \return Start offset into ViolationSample::checkValues.
        typename std::vector<Scalar>::size_type
        violationValueStart(const std::size_t checkIx,
                            const std::size_t sampleIx) const;

        /// Emit violation report header.
        ///
        /// Includes such information as the \code Check::description()
        /// \endcode, the \code Check::condition() \endcode and the total
        /// number of violations of this specific check.
        ///
        /// \param[in] currentCheck Current check object.  Needed for the
        ///    description and condition values.
        ///
        /// \param[in] violationCount Total number of condition violations.
        ///
        /// \param[in] emitReportRecord Call-back function for outputting a
        ///    single record/line of a violation report.
        void writeReportHeader(const Check*              currentCheck,
                               const std::size_t         violationCount,
                               const ReportRecordOutput& emitReportRecord) const;

        /// Emit tabulated sample of check violations.
        ///
        /// \param[in] nValueChar Minimum number of characters (field width)
        ///    needed to display a single floating-point number ("checkValues").
        ///
        /// \param[in] currentCheck Current check object.  Needed for the
        ///    column headers in the tabulated output.
        ///
        /// \param[in] violation Sample of consistency check violations.
        ///
        /// \param[in] checkIx Numerical check index in the range
        ///    [0..battery_.size()).
        ///
        /// \param[in] emitReportRecord Call-back function for outputting a
        ///    single record/line of a violation report.
        void writeTabulatedReportSample(const std::size_t         nValueChar,
                                        const Check*              currentCheck,
                                        const ViolationSample&    violation,
                                        const std::size_t         checkIx,
                                        const ReportRecordOutput& emitReportRecord) const;

        /// Generate textual representation of all sampled point IDs for a
        /// single consistency check.
        ///
        /// \param[in] violation Sample of consistency check violations.
        ///
        /// \param[in] checkIx Numerical check index in the range
        ///    [0..battery_.size()).
        ///
        /// \return Sequence of formatted point IDs and the minimum number
        ///    of characters needed to represent each one of these.
        std::pair<std::vector<std::string>, std::string::size_type>
        formatPointIDs(const ViolationSample& violation,
                       const std::size_t      checkIx) const;

        /// Generate sequence of table column headers for a single
        /// consistency check.
        ///
        /// \param[in] currentCheck Current check object.  Needed for the
        ///    \code Check::columnNames() \endcode.
        ///
        /// \return Copy of the check's column names.
        std::vector<std::string>
        collectColumnHeaders(const Check* currentCheck) const;

        /// Generate sequence of sample indices, sorted ascendingly on the
        /// corresponding point ID.
        ///
        /// \param[in] currentCheck Current check object.  Needed for the
        ///    \code Check::columnNames() \endcode.
        ///
        /// \param[in] checkIx Numerical check index in the range
        ///    [0..battery_.size()).
        ///
        /// \return Ascendingly sorted sample indices.
        std::vector<std::size_t>
        sortedPointIndices(const ViolationSample& violation,
                           const std::size_t      checkIx) const;

        /// Compute number of sample points for a single check's violations.
        ///
        /// Effectively the minimum of the number of violations of that
        /// check and the maximum number of sample points (\code
        /// this->numSamplePoints_ \endcode).
        ///
        /// \param[in] violation Sample of consistency check violations.
        ///
        /// \param[in] checkIx Numerical check index in the range
        ///    [0..battery_.size()).
        ///
        /// \return Number of active sample points.
        std::size_t numPoints(const ViolationSample& violation,
                              const std::size_t      checkIx) const;

        /// Compute number of sample points for a single check's violations.
        ///
        /// Effectively the minimum of the number of violations of that
        /// check and the maximum number of sample points (\code
        /// this->numSamplePoints_ \endcode).
        ///
        /// \param[in] violationCount Total number of check violations.
        ///
        /// \return Number of active sample points.
        std::size_t numPoints(const std::size_t violationCount) const;

        /// Whether or not any checks failed at specified severity level.
        ///
        /// \param[in] level Violation severity level.
        ///
        /// \return Whether or not any checks failed at severity level \p
        /// level.
        bool anyFailedChecks(const ViolationLevel level) const;

        /// Convert severity level into collection index.
        ///
        /// \param[in] level Violation severity level.
        ///
        /// \return Corresponding index into \code this->violations_
        /// \endcode.
        auto index(const ViolationLevel level) const
        {
            return static_cast<typename ViolationCollection::size_type>(level);
        }

        /// Run block of code for each registered consistency check.
        ///
        /// Mutable version.
        ///
        /// \tparam Body Call-back function type encapsulating block of code
        /// to run for each registered consistency check.  Expected to be a
        /// callable type with a function call operator of the form
        /// \code
        ///    void Body::operator()(Check* checkPtr, std::size_t i)
        /// \endcode
        /// in which the \c checkPtr points to a (mutable) \c Check object
        /// and \c i is the sequence index in which the particular check was
        /// registered in a call to addCheck().
        ///
        /// \param[in,out] body Block of code to run for each registered
        /// consistency check.
        template <typename Body>
        void checkLoop(Body&& body);

        /// Run block of code for each registered consistency check.
        ///
        /// Immutable version.
        ///
        /// \tparam Body Call-back function type encapsulating block of code
        /// to run for each registered consistency check.  Expected to be a
        /// callable type with a function call operator of the form
        /// \code
        ///    void Body::operator()(const Check* checkPtr, std::size_t i) const
        /// \endcode
        /// in which the \c checkPtr points to an immutable \c Check object
        /// and \c i is the sequence index in which the particular check was
        /// registered in a call to addCheck().
        ///
        /// \param[in] body Block of code to run for each registered
        /// consistency check.
        template <typename Body>
        void checkLoop(Body&& body) const;
    };

} // namespace Opm

#endif // OPM_SATFUNC_CONSISTENCY_CHECK_MODULE_HPP
