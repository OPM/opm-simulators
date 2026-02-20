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

#include <opm/simulators/utils/satfunc/SatfuncConsistencyCheckManager.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/EndpointScaling.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldProps.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/EclipseState/Grid/SatfuncPropertyInitializers.hpp>
#include <opm/input/eclipse/EclipseState/Phase.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>

#include <opm/material/fluidmatrixinteractions/EclEpsGridProperties.hpp>
#include <opm/material/fluidmatrixinteractions/EclEpsScalingPoints.hpp>

#include <opm/simulators/utils/satfunc/SatfuncCheckPointInterface.hpp>
#include <opm/simulators/utils/satfunc/SatfuncConsistencyChecks.hpp>
#include <opm/simulators/utils/satfunc/ScaledSatfuncCheckPoint.hpp>
#include <opm/simulators/utils/satfunc/UnscaledSatfuncCheckPoint.hpp>

#include <opm/simulators/utils/satfunc/GasPhaseConsistencyChecks.hpp>
#include <opm/simulators/utils/satfunc/OilPhaseConsistencyChecks.hpp>
#include <opm/simulators/utils/satfunc/ThreePointHorizontalConsistencyChecks.hpp>
#include <opm/simulators/utils/satfunc/WaterPhaseConsistencyChecks.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace {
    Opm::satfunc::RawTableEndPoints
    rawTableEndpoints(const Opm::EclipseState& eclipseState)
    {
        const auto& rspec = eclipseState.runspec();

        return Opm::satfunc::getRawTableEndpoints
            (eclipseState.getTableManager(),
             rspec.phases(),
             rspec.saturationFunctionControls()
             .minimumRelpermMobilityThreshold());
    }

    Opm::satfunc::RawFunctionValues
    rawFunctionValues(const Opm::EclipseState&               eclipseState,
                      const Opm::satfunc::RawTableEndPoints& rtep)
    {
        return Opm::satfunc::getRawFunctionValues
            (eclipseState.getTableManager(),
             eclipseState.runspec().phases(), rtep);
    }
} // Anonymous namespace

// ---------------------------------------------------------------------------

template <typename Scalar>
Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
SatfuncConsistencyCheckManager(const std::size_t    numSamplePoints,
                               const EclipseState&  eclipseState,
                               const LocalToGlobal& localToGlobal)
    : eclipseState_  { std::cref(eclipseState) }
    , localToGlobal_ { localToGlobal }
    , rtep_          { rawTableEndpoints(eclipseState) }
    , rfunc_         { rawFunctionValues(eclipseState, rtep_) }
{
    // Note: This setup is limited to
    //   1. Drainage only--no hysteresis
    //   2. Non-directional relative permeability
    //   3. Relative permeability only--no capillary pressure

    this->configureCurveChecks(numSamplePoints);
    this->addChecks();
}

template <typename Scalar>
bool Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
anyFailedStandardChecks() const
{
    return std::any_of(this->curves_.begin(), this->curves_.end(),
                       [](const auto& curve)
                       { return curve.checks.anyFailedStandardChecks(); });
}

template <typename Scalar>
bool Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
anyFailedCriticalChecks() const
{
    return std::any_of(this->curves_.begin(), this->curves_.end(),
                       [](const auto& curve)
                       { return curve.checks.anyFailedCriticalChecks(); });
}

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
reportFailures(const ViolationLevel      level,
               const ReportRecordOutput& emitReportRecord) const
{
    if (! this->isRoot_) {
        return;
    }

    this->curveLoop([level, &emitReportRecord](const auto& curve)
    {
        curve.checks.reportFailures(level, emitReportRecord);
    });
}

// ===========================================================================
// Private member functions for SatfuncConsistencyCheckManager template
// ===========================================================================

template <typename Scalar>
Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
CurveCollection::CurveCollection(std::unique_ptr<SatfuncCheckPointInterface<Scalar>> point_arg,
                                 std::string_view  pointName,
                                 const std::size_t numSamplePoints)
    : point  { std::move(point_arg) }
    , checks { pointName, numSamplePoints }
{}

// ---------------------------------------------------------------------------

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
warnIfDirectionalOrIrreversibleEPS() const
{
    if (! this->isRoot_) { return; }

    if (const auto& eps = this->eclipseState_.get().runspec().endpointScaling(); !eps) {
        // End-point scaling not active in run.  Don't need to check
        // anything else.
        return;
    }
    else if (eps.directional() || eps.irreversible()) {
        OpmLog::warning("Directional and/or irreversible end-point "
                        "scaling is currently not included in the "
                        "saturation function consistency checks");
    }
}

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
runCellChecks(const int cellIdx)
{
    this->curveLoop([cellIdx, endPoints = EclEpsScalingPointsInfo<Scalar>{}]
                    (auto& curve) mutable
    {
        const auto pointID = curve.point->pointID(cellIdx);
        if (! pointID.has_value()) {
            // Check does not apply to this cell for 'curve'.  Might be
            // because it's a region based check and we already ran the
            // checks for this particular underlying region.
            return;
        }

        curve.point->populateCheckPoint(cellIdx, endPoints);
        curve.checks.checkEndpoints(*pointID, endPoints);
    });
}

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
configureCurveChecks(const std::size_t numSamplePoints)
{
    const auto unscaledChecks =
        this->configureUnscaledCurveChecks("SATNUM", numSamplePoints);

    if (unscaledChecks == nullptr) {
        // SATNUM array does not exist (unexpected), or end-point scaling is
        // not active in the current run.  There's no need to configure
        // consistency checks for the scaled curves.
        return;
    }

    const auto useImbibition = false;
    this->configureScaledCurveChecks(*unscaledChecks,
                                     useImbibition,
                                     numSamplePoints);
}

template <typename Scalar>
std::unique_ptr<Opm::Satfunc::PhaseChecks::UnscaledSatfuncCheckPoint<Scalar>>
Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
configureUnscaledCurveChecks(const std::string& regionName,
                             const std::size_t  numSamplePoints)
{
    const auto& fp = this->eclipseState_.get().fieldProps();

    if (! fp.has_int(regionName)) {
        // Region property array (SATNUM, IMBNUM, &c) not available.
        // Nothing to do.
        return {};
    }

    using UEP = typename UnscaledSatfuncCheckPoint<Scalar>::UnscaledEndPoints;

    const auto regIdxOffset = 1; // regIdx contains one-based region indices.
    auto unscaledChecks = std::make_unique<UnscaledSatfuncCheckPoint<Scalar>>
        (&fp.get_int(regionName), regIdxOffset, UEP { &this->rtep_, &this->rfunc_ });

    if (! this->eclipseState_.get().runspec().endpointScaling()) {
        // Include consistency checks for the unscaled/input/tabulated
        // saturation functions only if end-point scaling is NOT enabled.
        this->curves_.emplace_back
            (std::move(unscaledChecks), regionName, numSamplePoints);

        // Return nullptr because there are no scaled curves in this run and
        // we therefore do not need to configure consistency checks for such
        // curves.
        return {};
    }

    // If we get here then the run includes end-point scaling.  Return
    // sampling points on the unscaled curve as the fall-back points for the
    // scaled curve.  Returning a non-null pointer here also lets the caller
    // know that we need to configure the associate consistency checks for
    // the scaled curves.
    return unscaledChecks;
}

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
configureScaledCurveChecks(const UnscaledSatfuncCheckPoint<Scalar>& unscaledChecks,
                           const bool                               useImbibition,
                           const std::size_t                        numSamplePoints)
{
    this->gridProps_.emplace_back(this->eclipseState_, useImbibition);

    const auto& gdims = this->eclipseState_.get().gridDims();

    auto& curve = this->curves_.emplace_back
        (std::make_unique<ScaledSatfuncCheckPoint<Scalar>>
         (unscaledChecks, &this->eclipseState_.get(),
          &this->gridProps_.back(), this->localToGlobal_),
         "Grid Block", numSamplePoints);

    const auto nchar = std::max({
            fmt::formatted_size("{}", gdims.getNX()),
            fmt::formatted_size("{}", gdims.getNY()),
            fmt::formatted_size("{}", gdims.getNZ()),
        });

    curve.checks.setPointIDFormatCallback([nchar, gdims](const std::size_t globalCell)
    {
        const auto ijk = gdims.getIJK(globalCell);

        return fmt::format("({1:>{0}}, {2:>{0}}, {3:>{0}})", nchar,
                           ijk[0] + 1, ijk[1] + 1, ijk[2] + 1);
    });
}

namespace {

    /// Factory for creating individual end-point checks.
    ///
    /// \tparam Scalar Element type.  Typically \c float or \c double.
    template <typename Scalar>
    class CheckCreationFactory
    {
    public:
        /// Constructor
        ///
        /// \param[in] phases Run's active phases.  Needed to determine
        /// which end-point checks to include in the test set.
        ///
        /// \param[in] threePointScaling Whether or not run uses the
        /// alternative, three-point method for horizontal saturation
        /// function end-point scaling ("SCALECRS = YES").
        explicit CheckCreationFactory(const Opm::Phases& phases,
                                      const bool threePointScaling);

        /// Start of sequence of end-point check creation functions.
        auto begin() const { return this->creationFunctions_.begin(); }

        /// End of sequence of end-point check creation functions.
        auto end() const { return this->creationFunctions_.end(); }

    private:
        /// Convenience type alias for individual checks.
        using Check = typename Opm::SatfuncConsistencyChecks<Scalar>::Check;

        /// Type alias for a check creation function.
        using CreationFunction = std::function<std::unique_ptr<Check>()>;

        /// Collection of pertinent test creation functions.
        std::vector<CreationFunction> creationFunctions_{};

        /// Incorporate end-point checks for an active oil phase.
        ///
        /// \param[in] phases Run's active phases.  Needed to determine
        /// which of the two-phase G/O and/or O/W end-point checks to
        /// include in the test set.
        void addOilChecks(const Opm::Phases& phases);

        /// Incorporate end-point checks for the two-phase G/O system.
        void addGasOilChecks();

        /// Incorporate end-point checks for the two-phase O/W system.
        void addOilWaterChecks();

        /// Incorporate end-point checks for an active gas phase.
        void addGasChecks();

        /// Incorporate end-point checks for an active water phase.
        void addWaterChecks();

        /// Incorporate end-point checks for the alternative, three-point
        /// scaling method ("SCALECRS = YES").
        ///
        /// \param[in] phases Run's active phases.  Needed to determine
        /// which of the two-phase G/O and/or O/W end-point checks to
        /// include in the test set.
        void addThreePointChecks(const Opm::Phases& phases);
    };

    template <typename Scalar>
    CheckCreationFactory<Scalar>::CheckCreationFactory(const Opm::Phases& phases,
                                                       const bool threePointScaling)
    {
        if (phases.active(Opm::Phase::OIL)) {
            this->addOilChecks(phases);
        }

        if (phases.active(Opm::Phase::GAS)) {
            this->addGasChecks();
        }

        if (phases.active(Opm::Phase::WATER)) {
            this->addWaterChecks();
        }

        if (threePointScaling && phases.active(Opm::Phase::OIL)) {
            this->addThreePointChecks(phases);
        }
    }

    template <typename Scalar>
    void CheckCreationFactory<Scalar>::addOilChecks(const Opm::Phases& phases)
    {
        if (phases.active(Opm::Phase::GAS)) {
            this->addGasOilChecks();
        }

        if (phases.active(Opm::Phase::WATER)) {
            this->addOilWaterChecks();
        }
    }

    template <typename Scalar>
    void CheckCreationFactory<Scalar>::addGasOilChecks()
    {
        namespace OChecks = Opm::Satfunc::PhaseChecks::Oil;

        this->creationFunctions_.insert(this->creationFunctions_.end(), {
                CreationFunction { []() { return std::make_unique<OChecks::SOcr_GO<Scalar>>(); } },
                CreationFunction { []() { return std::make_unique<OChecks::SOmin_GO<Scalar>>(); } },
                CreationFunction { []() { return std::make_unique<OChecks::MobileOil_GO_SGmin<Scalar>>(); } },
                CreationFunction { []() { return std::make_unique<OChecks::MobileOil_GO_SGcr<Scalar>>(); } },
            });
    }

    template <typename Scalar>
    void CheckCreationFactory<Scalar>::addOilWaterChecks()
    {
        namespace OChecks = Opm::Satfunc::PhaseChecks::Oil;

        this->creationFunctions_.insert(this->creationFunctions_.end(), {
                CreationFunction { []() { return std::make_unique<OChecks::SOcr_OW<Scalar>>(); } },
                CreationFunction { []() { return std::make_unique<OChecks::SOmin_OW<Scalar>>(); } },
                CreationFunction { []() { return std::make_unique<OChecks::MobileOil_OW_SWmin<Scalar>>(); } },
                CreationFunction { []() { return std::make_unique<OChecks::MobileOil_OW_SWcr<Scalar>>(); } },
            });
    }

    template <typename Scalar>
    void CheckCreationFactory<Scalar>::addGasChecks()
    {
        namespace GChecks = Opm::Satfunc::PhaseChecks::Gas;

        this->creationFunctions_.insert(this->creationFunctions_.end(), {
                CreationFunction { []() { return std::make_unique<GChecks::SGmin<Scalar>>(); } },
                CreationFunction { []() { return std::make_unique<GChecks::SGmax<Scalar>>(); } },
                CreationFunction { []() { return std::make_unique<GChecks::SGcr <Scalar>>(); } },
            });
    }

    template <typename Scalar>
    void CheckCreationFactory<Scalar>::addWaterChecks()
    {
        namespace WChecks = Opm::Satfunc::PhaseChecks::Water;

        this->creationFunctions_.insert(this->creationFunctions_.end(), {
                CreationFunction { []() { return std::make_unique<WChecks::SWmin<Scalar>>(); } },
                CreationFunction { []() { return std::make_unique<WChecks::SWmax<Scalar>>(); } },
                CreationFunction { []() { return std::make_unique<WChecks::SWcr <Scalar>>(); } },
            });
    }

    template <typename Scalar>
    void CheckCreationFactory<Scalar>::addThreePointChecks(const Opm::Phases& phases)
    {
        namespace TChecks = Opm::Satfunc::PhaseChecks::ThreePointHorizontal;

        if (phases.active(Opm::Phase::GAS)) {
            this->creationFunctions_.emplace_back
                ([]() { return std::make_unique<TChecks::DisplacingOil_GO<Scalar>>(); });
        }

        if (phases.active(Opm::Phase::WATER)) {
            this->creationFunctions_.emplace_back
                ([]() { return std::make_unique<TChecks::DisplacingOil_OW<Scalar>>(); });
        }
    }

} // Anonymous namespace

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::addChecks()
{
    const auto& rspec = this->eclipseState_.get().runspec();

    const auto checkCreationFactory = CheckCreationFactory<Scalar> {
        rspec.phases(),
        [&eps = rspec.endpointScaling()]() {
            return eps && eps.threepoint();
        }()
    };

    this->curveLoop([&checkCreationFactory](auto& curve)
    {
        curve.checks.resetCheckSet();

        for (const auto& makeCheck : checkCreationFactory) {
            curve.checks.addCheck(makeCheck());
        }

        curve.checks.finaliseCheckSet();
    });
}

template <typename Scalar>
void Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
collectFailures(const Parallel::Communication& comm)
{
    this->curveLoop([root = this->root_, comm](auto& curve)
    {
        curve.checks.collectFailures(root, comm);
    });
}

template <typename Scalar>
template <typename Body>
void Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
curveLoop(Body&& body)
{
    std::ranges::for_each(this->curves_, std::forward<Body>(body));
}

template <typename Scalar>
template <typename Body>
void Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<Scalar>::
curveLoop(Body&& body) const
{
    std::ranges::for_each(this->curves_, std::forward<Body>(body));
}

// ===========================================================================
// Explicit Specialisations
//
// No other code below this separator
// ===========================================================================

template class Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<float>;
template class Opm::Satfunc::PhaseChecks::SatfuncConsistencyCheckManager<double>;
