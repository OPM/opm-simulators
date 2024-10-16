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

#ifndef SATFUNC_CONSISTENCY_CHECK_MANAGER_HPP_INCLUDED
#define SATFUNC_CONSISTENCY_CHECK_MANAGER_HPP_INCLUDED

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/input/eclipse/EclipseState/Grid/SatfuncPropertyInitializers.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsGridProperties.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/simulators/utils/satfunc/SatfuncCheckPointInterface.hpp>
#include <opm/simulators/utils/satfunc/SatfuncConsistencyChecks.hpp>

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace Opm::Satfunc::PhaseChecks {
    template <typename Scalar>
    class UnscaledSatfuncCheckPoint;
} // namespace Opm::Satfunc::PhaseChecks

namespace Opm::Satfunc::PhaseChecks {

    /// Define and execute saturation function consistency checks for all
    /// cells in model.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class SatfuncConsistencyCheckManager
    {
    public:
        /// Callback for translating active cell index to globally unique
        /// point ID.
        using LocalToGlobal = std::function<std::size_t(int)>;

        /// Call-back function type for outputting a single record of a
        /// consistency condition violation report.
        using ReportRecordOutput = typename
            SatfuncConsistencyChecks<Scalar>::ReportRecordOutput;

        /// Severity level for consistency condition violation.
        using ViolationLevel = typename
            SatfuncConsistencyChecks<Scalar>::ViolationLevel;

        /// Constructor
        ///
        /// Creates a collection of saturation function checks based on the
        /// characteristics of the simulation model, e.g., whether or not
        /// end-point scaling is active or whether or not the run uses the
        /// alternative (three-point) scaling method.
        ///
        /// \param[in] numSamplePoints Upper bound on the number of
        /// end-point check violations to preserve for reporting purposes.
        /// Should normally be a small number like 5 or 10.
        ///
        /// \param[in] eclipseState Container of static properties such as
        /// the scaled saturation function end-points.
        ///
        /// \param[in] localToGlobal Callback for translating active cell
        /// indices to globally unique point IDs.
        explicit SatfuncConsistencyCheckManager(const std::size_t    numSamplePoints,
                                                const EclipseState&  eclipseState,
                                                const LocalToGlobal& localToGlobal);

        /// Set rank to which failure reports should be collected
        ///
        /// \param[in] root Failure report destination rank.  Should
        /// normally be the run's I/O rank.
        ///
        /// \return \code *this \endcode
        SatfuncConsistencyCheckManager& collectFailuresTo(const int root)
        {
            this->root_ = root;
            return *this;
        }

        /// Execute collection of saturation function consistency checks for
        /// all cells in simulation model.
        ///
        /// \tparam GridView Dune grid view type.
        ///
        /// \tparam GetCellIndex Callback function type for translating an
        /// active cell object into a numeric index.  Assumed to support a
        /// function call operator of the form
        /// \code
        ///    int operator()(const Element& e)
        /// \endcode
        /// in which \c Element is the type representing a co-dimension zero
        /// entity in the grid view.
        ///
        /// \param[in] gv Grid view for which to analyse the saturation
        /// function consistency.  Each MPI rank will analyse its interior
        /// cells only, and any failure reports will be subsequently
        /// gathered on the root process defined by collectFailuresTo().
        ///
        /// \param[in] getCellIndex Callback function for computing a
        /// numeric lookup index associated to each interior element of the
        /// grid view.
        template <typename GridView, typename GetCellIndex>
        void run(const GridView& gv, GetCellIndex&& getCellIndex)
        {
            this->isRoot_ = gv.comm().rank() == this->root_;

            this->warnIfDirectionalOrIrreversibleEPS();

            for (const auto& elem : elements(gv, Dune::Partitions::interior)) {
                this->runCellChecks(getCellIndex(elem));
            }

            gv.comm().barrier();

            this->collectFailures(gv.comm());
        }

        /// Whether or not any checks failed at the \c Standard level.
        bool anyFailedStandardChecks() const;

        /// Whether or not any checks failed at the \c Critical level.
        bool anyFailedCriticalChecks() const;

        /// Generate textual summary output of all failed consistency checks
        /// at specific level.
        ///
        /// Reports only those conditions/checks for which there is at least
        /// one violation.
        ///
        /// In a parallel run it is only safe to call this function on the
        /// root process defined by collectFailuresTo().
        ///
        /// \param[in] level Report's severity level.
        ///
        /// \param[in] emitReportRecord Call-back function for outputting a
        /// single record/line of a violation report.  Typically a wrapper
        /// of \code OpmLog::info() \endcode.  It is the responsibility of
        /// emitReportRecord() to properly display the text lines to end
        /// users.
        void reportFailures(const ViolationLevel      level,
                            const ReportRecordOutput& emitReportRecord) const;

    private:
        /// Association between points on a specific saturation function
        /// curves and the saturation function consistency checks to run on
        /// those points.
        struct CurveCollection
        {
            /// Constructor
            ///
            /// Convenience only, as this enables constructing objects using
            /// vector<>::emplace_back().
            ///
            /// \param[in] point Callback protocol for defining and
            /// populating saturation function end-points on a single
            /// saturation function curve.  Typically represents either a
            /// collection of per-region, tabulated and unscaled saturation
            /// functions or a collection of per-cell scaled saturation
            /// functions.
            ///
            /// \param[in] pointName Name/category of the points in this set
            /// of checks.  Might for instance be "Grid block" or
            /// "Saturation region".  Will be forwarded as a constructor
            /// argument to \c SatfuncConsistencyChecks from whence it will
            /// be used as a column header.
            ///
            /// \param[in] numSamplePoints Upper bound on the number of
            /// end-point check violations to preserve for reporting
            /// purposes.  Will be forwarded as a constructor argument to \c
            /// SatfuncConsistencyChecks.
            explicit CurveCollection(std::unique_ptr<SatfuncCheckPointInterface<Scalar>> point,
                                     std::string_view  pointName,
                                     const std::size_t numSamplePoints);

            /// Callback protocol for defining and populating saturation
            /// function end-points on a single saturation function curve.
            /// Typically represents either a collection of per-region,
            /// tabulated and unscaled saturation functions or a collection
            /// of per-cell scaled saturation functions.
            std::unique_ptr<SatfuncCheckPointInterface<Scalar>> point;

            /// Set of consistency checks to run against \c point.
            SatfuncConsistencyChecks<Scalar> checks;
        };

        /// Container of static properties such as the scaled saturation
        /// function end-points.
        ///
        /// Also used to query if and which end-point scaling behaviour is
        /// active in the run.
        std::reference_wrapper<const EclipseState> eclipseState_;

        /// Callback for translating active cell indices to globally unique
        /// point IDs.
        ///
        /// Mostly stored for convenience.  Could arguably be forwarded to
        /// the per-cell checks instead.
        LocalToGlobal localToGlobal_;

        /// Raw table end-points.
        ///
        /// Minimum, critical, and maximum saturation points for each phase
        /// for all tabulated saturation functions.
        satfunc::RawTableEndPoints rtep_{};

        /// Raw saturation function values.
        ///
        /// Maximum function values for all saturation functions in addition
        /// to relative permeability values at critical saturation points.
        satfunc::RawFunctionValues rfunc_{};

        /// Access interface for scaled saturation function end-points.
        ///
        /// Represented as a vector in order to support expansion to
        /// hysteretic cases and/or directionally dependent end-point
        /// scaling.
        std::vector<EclEpsGridProperties> gridProps_{};

        /// All saturation function checks that will be run for all interior
        /// cells in a grid view.
        std::vector<CurveCollection> curves_{};

        /// Rank to which failure reports should be collected.
        int root_{0};

        /// Whether or not the current rank coincides with \c root_ in the
        /// grid view's communicator.
        bool isRoot_{false};

        /// Issue a warning on the \c root_ rank if the run uses directional
        /// or irreversible end-point scaling.
        ///
        /// Those scaled curves are currently not included in the saturation
        /// function consistency analysis.
        void warnIfDirectionalOrIrreversibleEPS() const;

        /// Run all configured saturation function checks for a single
        /// active cell.
        ///
        /// \param[in] cellIdx Numeric lookup index associated to an
        /// interior element/cell of a grid view.
        void runCellChecks(const int cellIdx);

        /// Configure all pertinent saturation function consistency checks.
        ///
        /// \param[in] numSamplePoints Upper bound on the number of
        /// end-point check violations to preserve for reporting purposes.
        void configureCurveChecks(const std::size_t numSamplePoints);

        /// Configure saturation function consistency checks for per-region,
        /// unscaled saturation functions.
        ///
        /// \param[in] regionName Region set for which to configure
        /// consistency checks.  Typically a well-known saturation function
        /// region property name like SATNUM or IMBNUM.
        ///
        /// \param[in] numSamplePoints Upper bound on the number of
        /// end-point check violations to preserve for reporting purposes.
        ///
        /// \return Callbacks for inferring the unscaled end-points of the
        /// saturation region \p regionName.  Nullptr if the region index
        /// property array does not exist.  This should, arguably, be an
        /// optional<> instead to better reflect the intended semantics, but
        /// then we would also need to include the header for class template
        /// UnscaledSatfuncCheckPoint<> here.
        std::unique_ptr<UnscaledSatfuncCheckPoint<Scalar>>
        configureUnscaledCurveChecks(const std::string& regionName,
                                     const std::size_t  numSamplePoints);

        /// Configure saturation function consistency checks for per-cell,
        /// scaled saturation functions.
        ///
        /// \param[in] unscaledChecks Callbacks for inferring the unscaled
        /// end-points of the underlying saturation region.  Typically the
        /// return value from a previous call to member function
        /// configureUnscaledCurveChecks().
        ///
        /// \param[in] useImbibition Whether or not to configure consistency
        /// checks for the imbibition curves.
        ///
        /// \param[in] numSamplePoints Upper bound on the number of
        /// end-point check violations to preserve for reporting purposes.
        void configureScaledCurveChecks(const UnscaledSatfuncCheckPoint<Scalar>& unscaledChecks,
                                        const bool useImbibition,
                                        const std::size_t numSamplePoints);

        /// Add set of particular end-point checks to each configured curve
        void addChecks();

        /// Collect consistency violations from all ranks in MPI communicator.
        ///
        /// Incorporates violation counts and sampled failure points into
        /// the internal structures on each rank.  Aggregate results useful
        /// for subsequent call to reportFailures() on root process.
        ///
        /// \param[in] comm MPI communication object.
        void collectFailures(const Parallel::Communication& comm);

        /// Run a function for each configured curve.
        ///
        /// Mutable version.
        ///
        /// \tparam Body Callback function type representing a block of code
        /// to run for each configured curve.  Typically a class generated
        /// by a lambda expression.
        ///
        /// \param[in] body Block of code to run for each configured curve.
        /// May mutate curves in place.
        template <typename Body>
        void curveLoop(Body&& body);

        /// Run a function for each configured curve.
        ///
        /// Immutable version.
        ///
        /// \tparam Body Callback function type representing a block of code
        /// to run for each configured curve.  Typically a class generated
        /// by a lambda expression.
        ///
        /// \param[in] body Block of code to run for each configured curve.
        /// May not mutate curves in place.
        template <typename Body>
        void curveLoop(Body&& body) const;
    };

} // namespace Opm::Satfunc::PhaseChecks

#endif // SATFUNC_CONSISTENCY_CHECK_MANAGER_HPP_INCLUDED
