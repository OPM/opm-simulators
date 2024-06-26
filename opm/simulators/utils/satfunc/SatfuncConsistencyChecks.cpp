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

#include <config.h>

#include <opm/simulators/utils/satfunc/SatfuncConsistencyChecks.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

#include <fmt/format.h>

// ===========================================================================
// Public member functions for SatfuncConsistencyChecks Template
// ===========================================================================

template <typename Scalar>
Opm::SatfuncConsistencyChecks<Scalar>::
SatfuncConsistencyChecks(std::string_view  pointName,
                         const std::size_t numSamplePoints)
    : pointName_       { pointName }
    , numSamplePoints_ { numSamplePoints }
    , formatPointID_   { [](const std::size_t i) { return fmt::format("{}", i); } }
{}

template <typename Scalar>
Opm::SatfuncConsistencyChecks<Scalar>::
SatfuncConsistencyChecks(SatfuncConsistencyChecks&& rhs)
    : pointName_        { std::move(rhs.pointName_) }
    , numSamplePoints_  { rhs.numSamplePoints_ }
    , formatPointID_    { std::move(rhs.formatPointID_) }
    , startCheckValues_ { std::move(rhs.startCheckValues_) }
    , violations_       { std::move(rhs.violations_) }
    , battery_          { std::move(rhs.battery_) }
{}

template <typename Scalar>
Opm::SatfuncConsistencyChecks<Scalar>&
Opm::SatfuncConsistencyChecks<Scalar>::operator=(SatfuncConsistencyChecks&& rhs)
{
    this->pointName_        = std::move(rhs.pointName_);
    this->numSamplePoints_  = rhs.numSamplePoints_;
    this->formatPointID_    = std::move(rhs.formatPointID_);
    this->startCheckValues_ = std::move(rhs.startCheckValues_);
    this->violations_       = std::move(rhs.violations_);
    this->battery_          = std::move(rhs.battery_);

    this->urbg_.reset();

    return *this;
}

template <typename Scalar>
void Opm::SatfuncConsistencyChecks<Scalar>::resetCheckSet()
{
    this->startCheckValues_.clear();
    this->startCheckValues_.push_back(0);

    for (auto& violation : this->violations_) {
        violation.clear();
    }

    this->battery_.clear();
    this->urbg_.reset();
}

template <typename Scalar>
void Opm::SatfuncConsistencyChecks<Scalar>::addCheck(std::unique_ptr<Check> check)
{
    this->battery_.push_back(std::move(check));

    const auto numCheckValues = this->battery_.back()->numExportedCheckValues();
    this->startCheckValues_.push_back(this->numSamplePoints_ * numCheckValues);
}

template <typename Scalar>
void Opm::SatfuncConsistencyChecks<Scalar>::finaliseCheckSet()
{
    std::partial_sum(this->startCheckValues_.begin(),
                     this->startCheckValues_.end(),
                     this->startCheckValues_.begin());

    for (auto& violation : this->violations_) {
        this->buildStructure(violation);
    }
}

template <typename Scalar>
void Opm::SatfuncConsistencyChecks<Scalar>::
checkEndpoints(const std::size_t                      pointID,
               const EclEpsScalingPointsInfo<Scalar>& endPoints)
{
    this->checkLoop([pointID, &endPoints, this]
                    (Check* currentCheck, const std::size_t checkIx)
    {
        currentCheck->test(endPoints);

        if (! currentCheck->isViolated()) {
            // Check holds for this set of end-points.  Nothing to do.
            return;
        }

        // If we get here then the check does not hold for this set of
        // end-points.  Process the violation at the prescribed level of
        // attention.  Critical violations typically end the run whereas
        // a standard level violation typically generates warnings only.

        const auto level = currentCheck->isCritical()
            ? ViolationLevel::Critical
            : ViolationLevel::Standard;

        this->processViolation(level, checkIx, pointID);
    });
}

template <typename Scalar>
bool Opm::SatfuncConsistencyChecks<Scalar>::anyFailedChecks() const
{
    return this->anyFailedChecks(ViolationLevel::Standard);
}

template <typename Scalar>
bool Opm::SatfuncConsistencyChecks<Scalar>::anyFailedCriticalChecks() const
{
    return this->anyFailedChecks(ViolationLevel::Critical);
}

template <typename Scalar>
void Opm::SatfuncConsistencyChecks<Scalar>::
reportFailures(const ViolationLevel      level,
               const ReportRecordOutput& emitReportRecord) const
{
    this->checkLoop([this,
                     &emitReportRecord,
                     nValueChar = fmt::formatted_size("{:> 8.6e}", 1.0),
                     &violation = this->violations_[this->index(level)]]
                    (const Check* currentCheck, const std::size_t checkIx)
    {
        if (violation.count[checkIx] == 0) {
            return;
        }

        this->writeReportHeader(currentCheck,
                                violation.count[checkIx],
                                emitReportRecord);

        this->writeTabulatedReportSample(nValueChar,
                                         currentCheck,
                                         violation,
                                         checkIx,
                                         emitReportRecord);
    });
}

// ===========================================================================
// Private member functions for SatfuncConsistencyChecks Template
// ===========================================================================

template <typename Scalar>
void Opm::SatfuncConsistencyChecks<Scalar>::ViolationSample::clear()
{
    this->count.clear();
    this->pointID.clear();
    this->checkValues.clear();
}

// ---------------------------------------------------------------------------

template <typename Scalar>
void
Opm::SatfuncConsistencyChecks<Scalar>::
buildStructure(ViolationSample& violation)
{
    violation.count.assign(this->battery_.size(), 0);
    violation.pointID.resize(this->battery_.size() * this->numSamplePoints_,
                             static_cast<std::size_t>(0xdeadc0deUL));
    violation.checkValues.resize(this->startCheckValues_.back());

    if constexpr (std::numeric_limits<Scalar>::has_quiet_NaN) {
        std::fill(violation.checkValues.begin(),
                  violation.checkValues.end(),
                  std::numeric_limits<Scalar>::quiet_NaN());
    }
}

template <typename Scalar>
void Opm::SatfuncConsistencyChecks<Scalar>::
processViolation(const ViolationLevel level,
                 const std::size_t    checkIx,
                 const std::size_t    pointID)
{
    auto& violation = this->violations_[this->index(level)];

    const auto nViol = ++violation.count[checkIx];

    // Special case handling for number of violations not exceeding number
    // of sample points.  Needed in order to guarantee that the full table
    // is populated before starting the random replacement stage.
    const auto sampleIx = (nViol <= this->numSamplePoints_)
        ? (nViol - 1)
        : this->getSampleIndex(nViol);

    if (sampleIx >= this->numSamplePoints_) {
        // Reservoir sampling algorithm
        // (https://en.wikipedia.org/wiki/Reservoir_sampling) says that this
        // particular set of end-points should *not* be included in the
        // reported violations.  No more work needed in this case.
        return;
    }

    // If we get here, then this set of end-points should be included in the
    // reported violations.  Record the pointID and the corresponding check
    // values in their appropriate locations.

    violation.pointID[checkIx*this->numSamplePoints_ + sampleIx] = pointID;

    auto* exportedCheckValues = violation.checkValues.data()
        + this->violationValueStart(checkIx, sampleIx);

    this->battery_[checkIx]->exportCheckValues(exportedCheckValues);
}

namespace {

    std::vector<std::string::size_type>
    computeFieldWidths(const std::vector<std::string>& columnHeaders,
                       const std::string::size_type    minColWidth)
    {
        auto fieldWidths = std::vector<std::size_t>(columnHeaders.size());

        std::transform(columnHeaders.begin(), columnHeaders.end(),
                       fieldWidths.begin(),
                       [minColWidth](const std::string& header)
                       { return std::max(minColWidth, header.size()); });

        return fieldWidths;
    }

    std::string
    createTableSeparator(const std::string::size_type               fwPointID,
                         const std::vector<std::string::size_type>& fieldWidths)
    {
        using namespace fmt::literals;

        // Note: "+2" for one blank space on each side of the string value.
        auto separator = fmt::format("+{name:-<{width}}",
                                     "name"_a = "",
                                     "width"_a = fwPointID + 2);

        for (const auto& fieldWidth : fieldWidths) {
            separator += fmt::format("+{name:-<{width}}",
                                     "name"_a = "",
                                     "width"_a = fieldWidth + 2);
        }

        separator += '+';

        return separator;
    }

    template <typename EmitRecord>
    void writeTableHeader(const std::string_view::size_type          fwPointID,
                          std::string_view                           pointName,
                          const std::vector<std::string::size_type>& fieldWidths,
                          const std::vector<std::string>&            columnHeaders,
                          EmitRecord&&                               emitRecord)
    {
        using namespace fmt::literals;

        auto tableHeader = fmt::format("| {name:<{width}} ",
                                       "name"_a = pointName,
                                       "width"_a = fwPointID);

        for (auto colIx = 0*columnHeaders.size(); colIx < columnHeaders.size(); ++colIx) {
            tableHeader += fmt::format("| {name:<{width}} ",
                                       "name"_a = columnHeaders[colIx],
                                       "width"_a = fieldWidths[colIx]);
        }

        emitRecord(tableHeader + '|');
    }

    template <typename Scalar, typename EmitRecord>
    void writeTableRecord(const std::string_view::size_type          fwPointID,
                          std::string_view                           pointID,
                          const std::vector<std::string::size_type>& fieldWidths,
                          const Scalar*                              checkValues,
                          EmitRecord&&                               emitRecord)
    {
        using namespace fmt::literals;

        auto record = fmt::format("| {pointID:<{width}} ",
                                  "width"_a = fwPointID,
                                  "pointID"_a = pointID);

        for (auto colIx = 0*fieldWidths.size(); colIx < fieldWidths.size(); ++colIx) {
            record += fmt::format("| {checkValue:>{width}.6e} ",
                                  "width"_a = fieldWidths[colIx],
                                  "checkValue"_a = checkValues[colIx]);
        }

        emitRecord(record + '|');
    }

} // Anonymous namespace

template <typename Scalar>
void Opm::SatfuncConsistencyChecks<Scalar>::
writeReportHeader(const Check*              currentCheck,
                  const std::size_t         violationCount,
                  const ReportRecordOutput& emitReportRecord) const
{
    const auto* sampleMsg = (violationCount > this->numSamplePoints_)
        ? "Sample Violations"
        : "List of Violations";

    emitReportRecord(fmt::format("Consistency Problem:\n"
                                 "  {}\n"
                                 "  {}\n"
                                 "  Total Violations: {}\n\n"
                                 "{}",
                                 currentCheck->description(),
                                 currentCheck->condition(),
                                 violationCount, sampleMsg));
}

template <typename Scalar>
void Opm::SatfuncConsistencyChecks<Scalar>::
writeTabulatedReportSample(const std::size_t         nValueChar,
                           const Check*              currentCheck,
                           const ViolationSample&    violation,
                           const std::size_t         checkIx,
                           const ReportRecordOutput& emitReportRecord) const
{
    const auto formattedPointIDs = this->formatPointIDs(violation, checkIx);
    const auto fieldWidthPointID =
        std::max(formattedPointIDs.second, this->pointName_.size());

    const auto columnHeaders = this->collectColumnHeaders(currentCheck);
    const auto fieldWidths   = computeFieldWidths(columnHeaders, nValueChar);

    const auto separator = createTableSeparator(fieldWidthPointID, fieldWidths);

    // Output separator to start table output.
    emitReportRecord(separator);

    // Output column headers.
    writeTableHeader(fieldWidthPointID, this->pointName_,
                     fieldWidths, columnHeaders,
                     emitReportRecord);

    // Output separator to start table value output.
    emitReportRecord(separator);

    // Emit sampled check violations in order sorted on the pointID.
    for (const auto& i : this->sortedPointIndices(violation, checkIx)) {
        const auto* checkValues = violation.checkValues.data()
            + this->violationValueStart(checkIx, i);

        writeTableRecord(fieldWidthPointID, formattedPointIDs.first[i],
                         fieldWidths, checkValues,
                         emitReportRecord);
    }

    // Output separator to end table output.
    //
    // Note: We emit two blank lines after final separator in order to
    // generate some vertical space for the case of multiple failing checks.
    emitReportRecord(fmt::format("{}\n\n", separator));
}

template <typename Scalar>
std::pair<std::vector<std::string>, std::string::size_type>
Opm::SatfuncConsistencyChecks<Scalar>::
formatPointIDs(const ViolationSample& violation,
               const std::size_t      checkIx) const
{
    auto formattedPointIDs = std::pair
        <std::vector<std::string>,
         std::string::size_type>
        {
            std::piecewise_construct,
            std::forward_as_tuple(),
            std::forward_as_tuple(std::string::size_type{0})
        };

    const auto nPoints = this->numPoints(violation, checkIx);

    formattedPointIDs.first.reserve(nPoints);

    const auto* pointIDs = violation.pointID.data()
        + (checkIx * this->numSamplePoints_);

    for (auto point = 0*nPoints; point < nPoints; ++point) {
        formattedPointIDs.first.push_back
            (this->formatPointID_(pointIDs[point]));

        formattedPointIDs.second =
            std::max(formattedPointIDs.second,
                     formattedPointIDs.first.back().size());
    }

    return formattedPointIDs;
}

template <typename Scalar>
std::vector<std::string>
Opm::SatfuncConsistencyChecks<Scalar>::
collectColumnHeaders(const Check* currentCheck) const
{
    auto headers = std::vector<std::string>
        (currentCheck->numExportedCheckValues());

    currentCheck->columnNames(headers.data());

    return headers;
}

template <typename Scalar>
std::vector<std::size_t>
Opm::SatfuncConsistencyChecks<Scalar>::
sortedPointIndices(const ViolationSample& violation,
                   const std::size_t      checkIx) const
{
    auto sortedIdxs = std::vector<std::size_t>
        (this->numPoints(violation, checkIx));

    std::iota(sortedIdxs.begin(), sortedIdxs.end(), std::size_t{0});

    std::sort(sortedIdxs.begin(), sortedIdxs.end(),
              [pointIDs = violation.pointID.data() + (checkIx * this->numSamplePoints_)]
              (const std::size_t i1, const std::size_t i2)
              {
                  return pointIDs[i1] < pointIDs[i2];
              });

    return sortedIdxs;
}

template <typename Scalar>
std::size_t
Opm::SatfuncConsistencyChecks<Scalar>::
numPoints(const ViolationSample& violation,
          const std::size_t      checkIx) const
{
    return std::min(this->numSamplePoints_,
                    violation.count[checkIx]);
}

template <typename Scalar>
std::size_t
Opm::SatfuncConsistencyChecks<Scalar>::
getSampleIndex(const std::size_t sampleSize)
{
    assert (sampleSize > 0);

    this->ensureRandomBitGeneratorIsInitialised();

    return std::uniform_int_distribution<std::size_t>
        { 0, sampleSize - 1 }(*this->urbg_);
}

template <typename Scalar>
void Opm::SatfuncConsistencyChecks<Scalar>::ensureRandomBitGeneratorIsInitialised()
{
    if (this->urbg_ != nullptr) {
        return;
    }

    const auto k = static_cast<std::size_t>
        (std::log2(RandomBitGenerator::modulus) / 32) + 1;

    auto state = std::vector<typename RandomBitGenerator::result_type>(k + 3);

    std::random_device rd{};
    std::generate(state.begin(), state.end(), std::ref(rd));

    std::seed_seq seeds(state.begin(), state.end());
    this->urbg_ = std::make_unique<RandomBitGenerator>(seeds);
}

template <typename Scalar>
typename std::vector<Scalar>::size_type
Opm::SatfuncConsistencyChecks<Scalar>::
violationValueStart(const std::size_t checkIx,
                    const std::size_t sampleIx) const
{
    return this->startCheckValues_[checkIx]
        + (sampleIx * this->battery_[checkIx]->numExportedCheckValues());
}

template <typename Scalar>
bool
Opm::SatfuncConsistencyChecks<Scalar>::
anyFailedChecks(const ViolationLevel level) const
{
    const auto& violation = this->violations_[this->index(level)];

    return std::any_of(violation.count.begin(), violation.count.end(),
                       [](const std::size_t n) { return n > 0; });
}

template <typename Scalar>
template <typename Body>
void Opm::SatfuncConsistencyChecks<Scalar>::checkLoop(Body&& body)
{
    const auto numChecks = this->battery_.size();

    for (auto checkIx = 0*numChecks; checkIx < numChecks; ++checkIx) {
        body(this->battery_[checkIx].get(), checkIx);
    }
}

template <typename Scalar>
template <typename Body>
void Opm::SatfuncConsistencyChecks<Scalar>::checkLoop(Body&& body) const
{
    const auto numChecks = this->battery_.size();

    for (auto checkIx = 0*numChecks; checkIx < numChecks; ++checkIx) {
        body(this->battery_[checkIx].get(), checkIx);
    }
}

// ===========================================================================
// Explicit Specialisations of SatfuncConsistencyChecks Template
//
// No other code below this separator
// ===========================================================================

template class Opm::SatfuncConsistencyChecks<float>;
template class Opm::SatfuncConsistencyChecks<double>;
