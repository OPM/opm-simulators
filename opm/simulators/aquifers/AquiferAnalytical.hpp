/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS

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

#ifndef OPM_AQUIFERANALYTICAL_HEADER_INCLUDED
#define OPM_AQUIFERANALYTICAL_HEADER_INCLUDED

#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/EclipseState/Aquifer/Aquancon.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/utils/basicproperties.hh>

#include <opm/output/data/Aquifer.hpp>

#include <opm/simulators/aquifers/AquiferInterface.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <optional>
#include <vector>

namespace Opm
{
template <typename TypeTag>
class AquiferAnalytical : public AquiferInterface<TypeTag>
{
public:
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BlackoilIndices = GetPropType<TypeTag, Properties::Indices>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;

    enum { enableTemperature = getPropValue<TypeTag, Properties::EnableTemperature>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableBrine = getPropValue<TypeTag, Properties::EnableBrine>() };
    enum { enableVapwat = getPropValue<TypeTag, Properties::EnableVapwat>() };
    enum { has_disgas_in_water = getPropValue<TypeTag, Properties::EnableDisgasInWater>() };

    enum { enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>() };

    static constexpr int numEq = BlackoilIndices::numEq;
    using Scalar = double;

    using Eval = DenseAd::Evaluation<double, /*size=*/numEq>;

    using FluidState = BlackOilFluidState<Eval,
                                          FluidSystem,
                                          enableTemperature,
                                          enableEnergy,
                                          BlackoilIndices::gasEnabled,
                                          enableVapwat,
                                          enableBrine,
                                          enableSaltPrecipitation,
                                          has_disgas_in_water,
                                          BlackoilIndices::numPhases>;

    // Constructor
    AquiferAnalytical(const int aqID,
                      const std::vector<Aquancon::AquancCell>& connections,
                      const Simulator& ebosSimulator)
        : AquiferInterface<TypeTag>(aqID, ebosSimulator)
        , connections_(connections)
    {
        this->initializeConnectionMappings();
    }

    // Destructor
    virtual ~AquiferAnalytical()
    {}

    void computeFaceAreaFraction(const std::vector<double>& total_face_area) override
    {
        assert (total_face_area.size() >= static_cast<std::vector<double>::size_type>(this->aquiferID()));

        const auto tfa = total_face_area[this->aquiferID() - 1];
        const auto eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());

        if (tfa < eps_sqrt) {
            this->alphai_.assign(this->size(), Scalar{0});
        }
        else {
            std::transform(this->faceArea_connected_.begin(),
                           this->faceArea_connected_.end(),
                           this->alphai_.begin(),
                           [tfa](const Scalar area)
                           {
                               return area / tfa;
                           });
        }

        this->area_fraction_ = this->totalFaceArea() / tfa;
    }

    double totalFaceArea() const override
    {
        return this->total_face_area_;
    }

    void initFromRestart(const data::Aquifers& aquiferSoln) override
    {
        auto xaqPos = aquiferSoln.find(this->aquiferID());
        if (xaqPos == aquiferSoln.end())
            return;

        this->assignRestartData(xaqPos->second);

        this->W_flux_ = xaqPos->second.volume * this->area_fraction_;
        this->pa0_ = xaqPos->second.initPressure;

        this->solution_set_from_restart_ = true;
    }

    void initialSolutionApplied() override
    {
        initQuantities();
    }

    void beginTimeStep() override
    {
        ElementContext elemCtx(this->ebos_simulator_);
        OPM_BEGIN_PARALLEL_TRY_CATCH();

        for (const auto& elem : elements(this->ebos_simulator_.gridView())) {
            elemCtx.updatePrimaryStencil(elem);

            const int cellIdx = elemCtx.globalSpaceIndex(0, 0);
            const int idx = cellToConnectionIdx_[cellIdx];
            if (idx < 0)
                continue;

            elemCtx.updateIntensiveQuantities(0);
            const auto& iq = elemCtx.intensiveQuantities(0, 0);
            pressure_previous_[idx] = getValue(iq.fluidState().pressure(this->phaseIdx_()));
        }

        OPM_END_PARALLEL_TRY_CATCH("AquiferAnalytical::beginTimeStep() failed: ",
                                   this->ebos_simulator_.vanguard().grid().comm());
    }

    void addToSource(RateVector& rates,
                     const unsigned cellIdx,
                     const unsigned timeIdx) override
    {
        const auto& model = this->ebos_simulator_.model();

        const int idx = this->cellToConnectionIdx_[cellIdx];
        if (idx < 0)
            return;

        const auto& intQuants = model.intensiveQuantities(cellIdx, timeIdx);

        // This is the pressure at td + dt
        this->updateCellPressure(this->pressure_current_, idx, intQuants);
        this->calculateInflowRate(idx, this->ebos_simulator_);

        rates[BlackoilIndices::conti0EqIdx + compIdx_()]
            += this->Qai_[idx] / model.dofTotalVolume(cellIdx);

        if constexpr (enableEnergy) {
            auto fs = intQuants.fluidState();
            if (this->Ta0_.has_value() && this->Qai_[idx] > 0)
            {
                fs.setTemperature(this->Ta0_.value());
                typedef typename std::decay<decltype(fs)>::type::Scalar FsScalar;
                typename FluidSystem::template ParameterCache<FsScalar> paramCache;
                const unsigned pvtRegionIdx = intQuants.pvtRegionIndex();
                paramCache.setRegionIndex(pvtRegionIdx);
                paramCache.setMaxOilSat(this->ebos_simulator_.problem().maxOilSaturation(cellIdx));
                paramCache.updatePhase(fs, this->phaseIdx_());
                const auto& h = FluidSystem::enthalpy(fs, paramCache, this->phaseIdx_());
                fs.setEnthalpy(this->phaseIdx_(), h);
            }
            rates[BlackoilIndices::contiEnergyEqIdx]
            += this->Qai_[idx] *fs.enthalpy(this->phaseIdx_()) * FluidSystem::referenceDensity( this->phaseIdx_(), intQuants.pvtRegionIndex()) / model.dofTotalVolume(cellIdx);

        }
    }

    std::size_t size() const
    {
        return this->connections_.size();
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(pressure_previous_);
        serializer(pressure_current_);
        serializer(Qai_);
        serializer(rhow_);
        serializer(W_flux_);
    }

    bool operator==(const AquiferAnalytical& rhs) const
    {
        return this->pressure_previous_ == rhs.pressure_previous_ &&
               this->pressure_current_ == rhs.pressure_current_ &&
               this->Qai_ == rhs.Qai_ &&
               this->rhow_ == rhs.rhow_ &&
               this->W_flux_ == rhs.W_flux_;
    }

protected:
    virtual void assignRestartData(const data::AquiferData& xaq) = 0;
    virtual void calculateInflowRate(int idx, const Simulator& simulator) = 0;
    virtual void calculateAquiferCondition() = 0;
    virtual void calculateAquiferConstants() = 0;
    virtual Scalar aquiferDepth() const = 0;

    Scalar gravity_() const
    {
        return this->ebos_simulator_.problem().gravity()[2];
    }

    int compIdx_() const
    {
        if (this->co2store_or_h2store_())
            return FluidSystem::oilCompIdx;

        return FluidSystem::waterCompIdx;
    }

    void initQuantities()
    {
        // We reset the cumulative flux at the start of any simulation, so, W_flux = 0
        if (! this->solution_set_from_restart_) {
            W_flux_ = Scalar{0};
        }

        this->initializeConnectionDepths();
        this->calculateAquiferCondition();
        this->calculateAquiferConstants();

        this->pressure_previous_.resize(this->size(), Scalar{0});
        this->pressure_current_.resize(this->size(), Scalar{0});
        this->Qai_.resize(this->size(), Scalar{0});
    }

    void updateCellPressure(std::vector<Eval>& pressure_water,
                            const int idx,
                            const IntensiveQuantities& intQuants)
    {
        const auto& fs = intQuants.fluidState();
        pressure_water.at(idx) = fs.pressure(this->phaseIdx_());
    }

    void updateCellPressure(std::vector<Scalar>& pressure_water,
                            const int idx,
                            const IntensiveQuantities& intQuants)
    {
        const auto& fs = intQuants.fluidState();
        pressure_water.at(idx) = fs.pressure(this->phaseIdx_()).value();
    }

    void initializeConnectionMappings()
    {
        this->alphai_.resize(this->size(), 1.0);
        this->faceArea_connected_.resize(this->size(), Scalar{0});

        // total_face_area_ is the sum of the areas connected to an aquifer
        this->total_face_area_ = Scalar{0};
        this->cellToConnectionIdx_.resize(this->ebos_simulator_.gridView().size(/*codim=*/0), -1);
        const auto& gridView = this->ebos_simulator_.vanguard().gridView();
        for (std::size_t idx = 0; idx < this->size(); ++idx) {
            const auto global_index = this->connections_[idx].global_index;
            const int cell_index = this->ebos_simulator_.vanguard().compressedIndex(global_index);
            if (cell_index < 0) {
                continue;
            }

            auto elemIt = gridView.template begin</*codim=*/ 0>();
            std::advance(elemIt, cell_index);

            // The global_index is not part of this grid
            if (elemIt->partitionType() != Dune::InteriorEntity) {
                continue;
            }

            this->cellToConnectionIdx_[cell_index] = idx;
        }

        // Translate the C face tag into the enum used by opm-parser's TransMult class
        FaceDir::DirEnum faceDirection;

        // Get areas for all connections
        const auto& elemMapper = this->ebos_simulator_.model().dofMapper();
        for (const auto& elem : elements(gridView)) {
            const unsigned cell_index = elemMapper.index(elem);
            const int idx = this->cellToConnectionIdx_[cell_index];

            // Only deal with connections given by the aquifer
            if (idx < 0) {
                continue;
            }

            for (const auto& intersection : intersections(gridView, elem)) {
                // Only deal with grid boundaries
                if (! intersection.boundary()) {
                    continue;
                }

                switch (intersection.indexInInside()) {
                case 0:
                    faceDirection = FaceDir::XMinus;
                    break;
                case 1:
                    faceDirection = FaceDir::XPlus;
                    break;
                case 2:
                    faceDirection = FaceDir::YMinus;
                    break;
                case 3:
                    faceDirection = FaceDir::YPlus;
                    break;
                case 4:
                    faceDirection = FaceDir::ZMinus;
                    break;
                case 5:
                    faceDirection = FaceDir::ZPlus;
                    break;
                default:
                    OPM_THROW(std::logic_error,
                              "Internal error in initialization of aquifer.");
                }

                if (faceDirection == this->connections_[idx].face_dir) {
                    this->faceArea_connected_[idx] = this->connections_[idx].influx_coeff;
                    break;
                }
            }

            this->total_face_area_ += this->faceArea_connected_.at(idx);
        }
    }

    void initializeConnectionDepths()
    {
        this->cell_depth_.resize(this->size(), this->aquiferDepth());

        const auto& gridView = this->ebos_simulator_.vanguard().gridView();
        for (std::size_t idx = 0; idx < this->size(); ++idx) {
            const int cell_index = this->ebos_simulator_.vanguard()
                .compressedIndex(this->connections_[idx].global_index);
            if (cell_index < 0) {
                continue;
            }

            auto elemIt = gridView.template begin</*codim=*/ 0>();
            std::advance(elemIt, cell_index);

            // The global_index is not part of this grid
            if (elemIt->partitionType() != Dune::InteriorEntity) {
                continue;
            }

            this->cell_depth_.at(idx) = this->ebos_simulator_.vanguard().cellCenterDepth(cell_index);
        }
    }

    // This function is for calculating the aquifer properties from equilibrium state with the reservoir
    Scalar calculateReservoirEquilibrium()
    {
        // Since the global_indices are the reservoir index, we just need to extract the fluidstate at those indices
        std::vector<Scalar> pw_aquifer;
        Scalar water_pressure_reservoir;

        ElementContext elemCtx(this->ebos_simulator_);
        const auto& gridView = this->ebos_simulator_.gridView();
        for (const auto& elem : elements(gridView)) {
            elemCtx.updatePrimaryStencil(elem);

            const auto cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto idx = this->cellToConnectionIdx_[cellIdx];
            if (idx < 0)
                continue;

            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            const auto& iq0 = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = iq0.fluidState();

            water_pressure_reservoir = fs.pressure(this->phaseIdx_()).value();
            const auto water_density = fs.density(this->phaseIdx_());

            const auto gdz =
                this->gravity_() * (this->cell_depth_[idx] - this->aquiferDepth());

            pw_aquifer.push_back(this->alphai_[idx] *
                (water_pressure_reservoir - water_density.value()*gdz));
        }

        // We take the average of the calculated equilibrium pressures.
        const auto& comm = this->ebos_simulator_.vanguard().grid().comm();

        Scalar vals[2];
        vals[0] = std::accumulate(this->alphai_.begin(), this->alphai_.end(), Scalar{0});
        vals[1] = std::accumulate(pw_aquifer.begin(), pw_aquifer.end(), Scalar{0});

        comm.sum(vals, 2);

        return vals[1] / vals[0];
    }

    const std::vector<Aquancon::AquancCell> connections_;

    // Grid variables
    std::vector<Scalar> faceArea_connected_;
    std::vector<int> cellToConnectionIdx_;

    // Quantities at each grid id
    std::vector<Scalar> cell_depth_;
    std::vector<Scalar> pressure_previous_;
    std::vector<Eval> pressure_current_;
    std::vector<Eval> Qai_;
    std::vector<Scalar> alphai_;

    Scalar Tc_{}; // Time constant
    Scalar pa0_{}; // initial aquifer pressure
    std::optional<Scalar> Ta0_{}; // initial aquifer temperature
    Scalar rhow_{};

    Scalar total_face_area_{};
    Scalar area_fraction_{Scalar{1}};

    Eval W_flux_;

    bool solution_set_from_restart_ {false};
    bool has_active_connection_on_proc_{false};
};

} // namespace Opm

#endif
