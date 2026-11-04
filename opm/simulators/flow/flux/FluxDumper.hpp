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

#ifndef OPM_FLUX_DUMPER_HPP
#define OPM_FLUX_DUMPER_HPP

#include <opm/io/eclipse/FluxFile.hpp>
#include <opm/simulators/flow/flux/FluxRegions.hpp>

#include <array>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace Opm {

class FluxDumper
{
public:
    using Mode = EclIO::FluxFile::Mode;
    using Sampling = EclIO::FluxFile::Sampling;

    struct ReportStepData {
        int reportStep = 0;
        int simStep = 0;
        double startTime = 0.0;
        double stepLength = 0.0;

        std::vector<double> pressures;
        std::vector<double> swat;
        std::vector<double> sgas;
        std::vector<double> rs;
        std::vector<double> rv;
        std::vector<double> temperature;

        /// Component mass rates across each boundary face, face-major over the
        /// active phases. The whole of the Flux-mode payload.
        std::vector<double> massRates;
        std::vector<double> relPerm;
        std::vector<double> capPressure;
        std::vector<double> externalRegionSums;
    };

    FluxDumper(std::string parentCaseName,
               int regionId,
               const std::array<int, 3>& parentDims,
               const FluxRegions::Region& region,
               Mode mode,
               Sampling sampling,
               int phaseMask,
               const std::vector<double>& boundaryTransmissibilities = {});

    static std::vector<double> aggregateRates(
        Sampling sampling,
        const std::vector<std::vector<double>>& rateSnapshots,
        const std::vector<double>& timeWeights);

    using FaceFluxAccessor = std::function<double(const FluxRegions::BoundaryFace&, EclIO::FluxFile::Phase)>;
    std::vector<double> makeFaceMajorRates(const FaceFluxAccessor& getFaceFlux) const;

    ReportStepData makeZeroFluxStep(int reportStep,
                                    int simStep,
                                    double startTime,
                                    double stepLength) const;

    void appendReportStep(const ReportStepData& stepData);

    /// Append everything buffered since the last call to the FLUX file.
    ///
    /// The file grows by one block per call rather than being rewritten, so
    /// this is cheap to do often and nothing accumulates in memory between
    /// calls. The first call also lays down the static section, so the summary
    /// keys and the sample intervals must be set before it.
    void flush(const std::string& filename, bool formatted = false);

    void setSummaryKeys(std::vector<std::string> summaryKeys);

    /// Append one snapshot of the parent summary vectors.
    ///
    /// Samples are independent of the report-step sequence.  \p time is the
    /// end of the interval the sample represents, in seconds, and \p values
    /// must have one entry per registered summary key.  Rate-type entries
    /// are expected to be time-averaged over the interval since the previous
    /// sample; see EclIO::FluxFile::SummarySample.
    void appendSummarySample(double time, std::vector<double> values);

    /// Record the minimum interval, in seconds, enforced between consecutive
    /// summary samples.  Diagnostic only.
    void setSummaryMinSampleInterval(double interval);

    /// Record the minimum interval, in seconds, enforced between consecutive
    /// boundary records.  Diagnostic only.
    void setBoundaryMinSampleInterval(double interval);

    void setBoundaryTransmissibilities(const std::vector<double>& boundaryTransmissibilities);

    //! \brief PVT region of the cell on the far side of each boundary face.
    void setBoundaryExteriorPvtRegions(const std::vector<int>& exteriorPvtRegions);

    //! \brief Record how deep the cell on the far side of each face is.
    //!
    //! \details The pressures written per report step are that cell's, taken
    //!   at its centre, while a reduced run imposes them at the face. Without
    //!   the depth it cannot carry the one to the other.
    void setBoundaryExteriorDepths(const std::vector<double>& exteriorDepths);

    //! \brief Geometric thermal transmissibility of each boundary face, both
    //!        half-cells in series.
    void setBoundaryThermalTransmissibilities(const std::vector<double>& thermalTransmissibilities);

    //! \brief Record which equilibration region the cell on the far side of
    //!        each face belongs to.
    void setBoundaryExteriorEquilRegions(const std::vector<int>& exteriorEquilRegions);

    //! \brief Record the producing run's threshold pressures.
    //!
    //! \details A defaulted THPRES entry is the largest initial potential
    //!   difference along a whole region boundary, so a sector holding part of
    //!   that boundary cannot arrive at the same number and has to be handed
    //!   this one.
    void setThresholdPressure(std::vector<double> thresholdPressure);
    void setParentWells(EclIO::FluxFile::ParentWells parentWells);
    void setParentActions(EclIO::FluxFile::ParentActions parentActions);

    const EclIO::FluxFile::Data& data() const;
    const std::vector<FluxRegions::BoundaryFace>& regionBoundaryFaces() const;

private:
    static int countPhases(int phaseMask);
    static bool hasFluxMode(Mode mode);
    static bool hasPressureMode(Mode mode);

    std::size_t numBoundaryFaces() const;
    std::size_t expectedRateSize() const;
    void validateStepData(const ReportStepData& stepData) const;

    std::vector<FluxRegions::BoundaryFace> regionBoundaryFaces_;

    EclIO::FluxFile::Data data_;

    /// Records and samples waiting to go into the next block. Cleared by
    /// flush(), so the run does not carry its whole history around.
    std::vector<EclIO::FluxFile::ReportStep> pendingRecords_;
    std::vector<EclIO::FluxFile::SummarySample> pendingSamples_;

    /// Created on the first flush, once the summary keys and intervals are
    /// settled, and kept so that later flushes append rather than truncate.
    std::unique_ptr<EclIO::FluxFile::Writer> writer_;

    /// Report step of the record appended last. Two records sharing one is what
    /// sub-report-step boundary output looks like; the accumulated series used
    /// to be inspected for this, but it no longer survives a flush.
    std::optional<int> lastReportStep_;
};

} // namespace Opm

#endif // OPM_FLUX_DUMPER_HPP
