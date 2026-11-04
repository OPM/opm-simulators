/*
  Copyright 2026 Equinor ASA

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

#include <opm/simulators/flow/flux/FluxDumper.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace Opm {

FluxDumper::FluxDumper(std::string parentCaseName,
                       const int regionId,
                       const std::array<int, 3>& parentDims,
                       const FluxRegions::Region& region,
                       const Mode mode,
                       const Sampling sampling,
                       const int phaseMask,
                       const std::vector<double>& boundaryTransmissibilities)
{
    if (!boundaryTransmissibilities.empty()
        && boundaryTransmissibilities.size() != region.boundaryFaces.size()) {
        OPM_THROW(std::invalid_argument,
                  "FluxDumper: transmissibility vector must match boundary face count");
    }

    this->data_.header.parentNx = parentDims[0];
    this->data_.header.parentNy = parentDims[1];
    this->data_.header.parentNz = parentDims[2];
    this->data_.header.boxI1 = region.box.i1;
    this->data_.header.boxJ1 = region.box.j1;
    this->data_.header.boxK1 = region.box.k1;
    this->data_.header.boxNx = region.box.nx();
    this->data_.header.boxNy = region.box.ny();
    this->data_.header.boxNz = region.box.nz();
    this->data_.header.numCells = static_cast<int>(region.localToGlobal.size());
    this->data_.header.numBoundaryFaces = static_cast<int>(region.boundaryFaces.size());
    this->data_.header.numReportSteps = 0;
    this->data_.header.mode = mode;
    this->data_.header.sampling = sampling;
    this->data_.header.phaseMask = phaseMask;
    this->data_.header.numPhases = countPhases(phaseMask);

    if (this->data_.header.numPhases <= 0) {
        OPM_THROW(std::invalid_argument, "FluxDumper: phase mask must contain at least one phase");
    }

    this->data_.names = {
        std::move(parentCaseName),
        std::to_string(regionId),
    };

    this->data_.localToGlobal = region.localToGlobal;
    this->regionBoundaryFaces_ = region.boundaryFaces;
    this->data_.boundaryFaces.reserve(region.boundaryFaces.size());

    for (std::size_t i = 0; i < region.boundaryFaces.size(); ++i) {
        const auto& src = region.boundaryFaces[i];

        double transmissibility = 0.0;
        if (!boundaryTransmissibilities.empty()) {
            transmissibility = boundaryTransmissibilities[i];
        }

        this->data_.boundaryFaces.push_back(
            EclIO::FluxFile::BoundaryFace{
                src.interiorLocalCell,
                static_cast<int>(src.direction),
                src.exteriorGlobalCell,
                transmissibility,
            });
    }
}

void FluxDumper::setBoundaryTransmissibilities(const std::vector<double>& boundaryTransmissibilities)
{
    if (boundaryTransmissibilities.size() != this->data_.boundaryFaces.size()) {
        OPM_THROW(std::invalid_argument,
                  "FluxDumper: transmissibility vector must match boundary face count");
    }

    for (std::size_t i = 0; i < this->data_.boundaryFaces.size(); ++i) {
        this->data_.boundaryFaces[i].transmissibility = boundaryTransmissibilities[i];
    }
}

const std::vector<FluxRegions::BoundaryFace>& FluxDumper::regionBoundaryFaces() const
{
    return this->regionBoundaryFaces_;
}

void FluxDumper::setBoundaryExteriorPvtRegions(const std::vector<int>& exteriorPvtRegions)
{
    if (exteriorPvtRegions.size() != this->data_.boundaryFaces.size()) {
        OPM_THROW(std::invalid_argument,
                  "FluxDumper: exterior PVT region vector must match boundary face count");
    }

    for (std::size_t i = 0; i < this->data_.boundaryFaces.size(); ++i) {
        this->data_.boundaryFaces[i].exteriorPvtRegion = exteriorPvtRegions[i];
    }
}

void FluxDumper::setBoundaryExteriorDepths(const std::vector<double>& exteriorDepths)
{
    if (exteriorDepths.size() != this->data_.boundaryFaces.size()) {
        OPM_THROW(std::invalid_argument,
                  "FluxDumper: exterior depth vector must match boundary face count");
    }

    for (std::size_t i = 0; i < this->data_.boundaryFaces.size(); ++i) {
        this->data_.boundaryFaces[i].exteriorDepth = exteriorDepths[i];
    }
}

void FluxDumper::setBoundaryThermalTransmissibilities(const std::vector<double>& thermalTransmissibilities)
{
    if (thermalTransmissibilities.size() != this->data_.boundaryFaces.size()) {
        OPM_THROW(std::invalid_argument,
                  "FluxDumper: thermal transmissibility vector must match boundary face count");
    }

    for (std::size_t i = 0; i < this->data_.boundaryFaces.size(); ++i) {
        this->data_.boundaryFaces[i].thermalTransmissibility = thermalTransmissibilities[i];
    }
}

void FluxDumper::setBoundaryExteriorEquilRegions(const std::vector<int>& exteriorEquilRegions)
{
    if (exteriorEquilRegions.size() != this->data_.boundaryFaces.size()) {
        OPM_THROW(std::invalid_argument,
                  "FluxDumper: exterior equilibration region vector must match "
                  "boundary face count");
    }

    for (std::size_t i = 0; i < this->data_.boundaryFaces.size(); ++i) {
        this->data_.boundaryFaces[i].exteriorEquilRegion = exteriorEquilRegions[i];
    }
}

void FluxDumper::setThresholdPressure(std::vector<double> thresholdPressure)
{
    this->data_.thresholdPressure = std::move(thresholdPressure);
}

void FluxDumper::setParentWells(EclIO::FluxFile::ParentWells parentWells)
{
    this->data_.parentWells = std::move(parentWells);
}

void FluxDumper::setParentActions(EclIO::FluxFile::ParentActions parentActions)
{
    this->data_.parentActions = std::move(parentActions);
}

std::vector<double> FluxDumper::aggregateRates(
    const Sampling sampling,
    const std::vector<std::vector<double>>& rateSnapshots,
    const std::vector<double>& timeWeights)
{
    if (rateSnapshots.empty()) {
        OPM_THROW(std::invalid_argument, "FluxDumper::aggregateRates: rateSnapshots must not be empty");
    }

    const auto expectedSize = rateSnapshots.front().size();
    for (const auto& snapshot : rateSnapshots) {
        if (snapshot.size() != expectedSize) {
            OPM_THROW(std::invalid_argument,
                      "FluxDumper::aggregateRates: all snapshots must have identical size");
        }
    }

    if (sampling == Sampling::Instant) {
        return rateSnapshots.back();
    }

    if (rateSnapshots.size() != timeWeights.size()) {
        OPM_THROW(std::invalid_argument,
                  "FluxDumper::aggregateRates: timeWeights size must match number of snapshots");
    }

    double totalWeight = 0.0;
    for (const auto weight : timeWeights) {
        if (weight < 0.0) {
            OPM_THROW(std::invalid_argument,
                      "FluxDumper::aggregateRates: time weights must be non-negative");
        }
        totalWeight += weight;
    }

    if (totalWeight <= 0.0) {
        OPM_THROW(std::invalid_argument,
                  "FluxDumper::aggregateRates: total time weight must be positive");
    }

    std::vector<double> out(expectedSize, 0.0);
    for (std::size_t t = 0; t < rateSnapshots.size(); ++t) {
        const auto weight = timeWeights[t] / totalWeight;
        for (std::size_t i = 0; i < expectedSize; ++i) {
            out[i] += rateSnapshots[t][i] * weight;
        }
    }

    return out;
}

std::vector<double> FluxDumper::makeFaceMajorRates(const FaceFluxAccessor& getFaceFlux) const
{
    if (!hasFluxMode(this->data_.header.mode)) {
        return {};
    }

    const std::array phases{
        EclIO::FluxFile::Phase::Oil,
        EclIO::FluxFile::Phase::Water,
        EclIO::FluxFile::Phase::Gas,
    };

    std::vector<double> out;
    out.reserve(this->expectedRateSize());

    for (const auto& face : this->regionBoundaryFaces_) {
        for (const auto phase : phases) {
            if (this->data_.header.hasPhase(phase)) {
                out.push_back(getFaceFlux(face, phase));
            }
        }
    }

    return out;
}

FluxDumper::ReportStepData FluxDumper::makeZeroFluxStep(const int reportStep,
                                                        const int simStep,
                                                        const double startTime,
                                                        const double stepLength) const
{
    ReportStepData step;
    step.reportStep = reportStep;
    step.simStep = simStep;
    step.startTime = startTime;
    step.stepLength = stepLength;

    if (hasFluxMode(this->data_.header.mode)) {
        step.massRates.assign(this->expectedRateSize(), 0.0);
    }

    if (hasPressureMode(this->data_.header.mode)) {
        const auto numFaces = this->numBoundaryFaces();
        step.pressures.assign(numFaces, 0.0);
        step.swat.assign(numFaces, 0.0);
        step.sgas.assign(numFaces, 0.0);
        step.rs.assign(numFaces, 0.0);
        step.rv.assign(numFaces, 0.0);
    }

    return step;
}

void FluxDumper::appendReportStep(const ReportStepData& stepData)
{
    this->validateStepData(stepData);

    EclIO::FluxFile::ReportStep step;
    step.reportStep = stepData.reportStep;
    step.simStep = stepData.simStep;
    step.startTime = stepData.startTime;
    step.stepLength = stepData.stepLength;
    step.pressures = stepData.pressures;
    step.swat = stepData.swat;
    step.sgas = stepData.sgas;
    step.rs = stepData.rs;
    step.rv = stepData.rv;
    step.temperature = stepData.temperature;
    step.massRates = stepData.massRates;
    step.relPerm = stepData.relPerm;
    step.capPressure = stepData.capPressure;
    step.externalRegionSums = stepData.externalRegionSums;

    if (!step.temperature.empty()) {
        this->data_.header.hasTemperature = true;
    }

    // Several records sharing a report step means the boundary data was
    // written at sub-report-step resolution.
    if (this->lastReportStep_.has_value()
        && (*this->lastReportStep_ == step.reportStep))
    {
        this->data_.header.boundaryPerTimestep = true;
    }
    this->lastReportStep_ = step.reportStep;

    this->pendingRecords_.push_back(std::move(step));
    this->data_.header.numReportSteps =
        static_cast<int>(this->data_.header.numReportSteps + 1);
}

void FluxDumper::flush(const std::string& filename, const bool formatted)
{
    if (this->writer_ == nullptr) {
        // Deferred to here rather than done in the constructor because the
        // summary keys and the sample intervals are set afterwards, and the
        // static section carries them.
        this->writer_ = std::make_unique<EclIO::FluxFile::Writer>(filename, formatted,
                                                                  this->data_);
    }

    this->writer_->appendRecords(this->pendingRecords_);
    this->writer_->appendSummarySamples(this->pendingSamples_);
    this->writer_->close();

    this->pendingRecords_.clear();
    this->pendingSamples_.clear();
}

void FluxDumper::setSummaryKeys(std::vector<std::string> summaryKeys)
{
    this->data_.summaryKeys = std::move(summaryKeys);
    this->data_.header.numSummaryKeys = static_cast<int>(this->data_.summaryKeys.size());
}

void FluxDumper::appendSummarySample(const double time, std::vector<double> values)
{
    if (values.size() != this->data_.summaryKeys.size()) {
        OPM_THROW(std::invalid_argument,
                  fmt::format("FluxDumper summary sample must contain {} values, got {}",
                              this->data_.summaryKeys.size(), values.size()));
    }

    auto& sample = this->pendingSamples_.emplace_back();
    sample.time = time;
    sample.values = std::move(values);

    this->data_.header.numSummarySamples =
        static_cast<int>(this->data_.header.numSummarySamples + 1);
    this->data_.header.summaryPerTimestep = true;
}

void FluxDumper::setSummaryMinSampleInterval(const double interval)
{
    this->data_.header.summaryMinSampleInterval = interval;
}

void FluxDumper::setBoundaryMinSampleInterval(const double interval)
{
    this->data_.header.boundaryMinSampleInterval = interval;
}

const EclIO::FluxFile::Data& FluxDumper::data() const
{
    return this->data_;
}

int FluxDumper::countPhases(int phaseMask)
{
    int count = 0;
    while (phaseMask != 0) {
        count += (phaseMask & 1);
        phaseMask >>= 1;
    }
    return count;
}

bool FluxDumper::hasFluxMode(const Mode mode)
{
    return (static_cast<int>(mode) & static_cast<int>(Mode::Flux)) != 0;
}

bool FluxDumper::hasPressureMode(const Mode mode)
{
    return (static_cast<int>(mode) & static_cast<int>(Mode::Pressure)) != 0;
}

std::size_t FluxDumper::numBoundaryFaces() const
{
    return this->data_.boundaryFaces.size();
}

std::size_t FluxDumper::expectedRateSize() const
{
    return this->numBoundaryFaces() * static_cast<std::size_t>(this->data_.header.numPhases);
}

void FluxDumper::validateStepData(const ReportStepData& stepData) const
{
    const auto numFaces = this->numBoundaryFaces();

    auto throwSizeError = [](const std::string& name, const std::size_t expected, const std::size_t actual) {
        std::ostringstream os;
        os << "FluxDumper: invalid " << name << " size. expected " << expected << ", got " << actual;
        OPM_THROW(std::invalid_argument, os.str());
    };

    if (hasFluxMode(this->data_.header.mode)) {
        if (stepData.massRates.size() != this->expectedRateSize()) {
            throwSizeError("massRates", this->expectedRateSize(), stepData.massRates.size());
        }
    }
    else if (!stepData.massRates.empty()) {
        throwSizeError("massRates", 0, stepData.massRates.size());
    }

    if (hasPressureMode(this->data_.header.mode)) {
        if (stepData.pressures.size() != numFaces) {
            throwSizeError("pressures", numFaces, stepData.pressures.size());
        }
        if (this->data_.header.hasPhase(EclIO::FluxFile::Phase::Water)) {
            if (stepData.swat.size() != numFaces) {
                throwSizeError("swat", numFaces, stepData.swat.size());
            }
        }
        else if (!stepData.swat.empty()) {
            throwSizeError("swat", 0, stepData.swat.size());
        }

        if (this->data_.header.hasPhase(EclIO::FluxFile::Phase::Gas)) {
            if (stepData.sgas.size() != numFaces) {
                throwSizeError("sgas", numFaces, stepData.sgas.size());
            }
        }
        else if (!stepData.sgas.empty()) {
            throwSizeError("sgas", 0, stepData.sgas.size());
        }
        if (!stepData.rs.empty() && stepData.rs.size() != numFaces) {
            throwSizeError("rs", numFaces, stepData.rs.size());
        }
        if (!stepData.rv.empty() && stepData.rv.size() != numFaces) {
            throwSizeError("rv", numFaces, stepData.rv.size());
        }
        if (!stepData.temperature.empty() && stepData.temperature.size() != numFaces) {
            throwSizeError("temperature", numFaces, stepData.temperature.size());
        }
    }
    else {
        if (!stepData.pressures.empty()) {
            throwSizeError("pressures", 0, stepData.pressures.size());
        }
        if (!stepData.swat.empty()) {
            throwSizeError("swat", 0, stepData.swat.size());
        }
        if (!stepData.sgas.empty()) {
            throwSizeError("sgas", 0, stepData.sgas.size());
        }
        if (!stepData.rs.empty()) {
            throwSizeError("rs", 0, stepData.rs.size());
        }
        if (!stepData.rv.empty()) {
            throwSizeError("rv", 0, stepData.rv.size());
        }
        if (!stepData.temperature.empty()) {
            throwSizeError("temperature", 0, stepData.temperature.size());
        }
    }
}

} // namespace Opm
