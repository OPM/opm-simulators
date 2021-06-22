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

#ifndef OPM_AQUIFERINTERFACE_HEADER_INCLUDED
#define OPM_AQUIFERINTERFACE_HEADER_INCLUDED

#include <opm/common/utility/numeric/linearInterpolation.hpp>
#include <opm/parser/eclipse/EclipseState/Aquifer/Aquancon.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace Opm
{
template <typename TypeTag>
class AquiferInterface
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

    static const int numEq = BlackoilIndices::numEq;
    typedef double Scalar;

    typedef DenseAd::Evaluation<double, /*size=*/numEq> Eval;

    typedef BlackOilFluidState<Eval,
                               FluidSystem,
                               enableTemperature,
                               enableEnergy,
                               BlackoilIndices::gasEnabled,
                               enableBrine,
                               BlackoilIndices::numPhases>
        FluidState;

    static const auto waterCompIdx = FluidSystem::waterCompIdx;
    static const auto waterPhaseIdx = FluidSystem::waterPhaseIdx;

    // Constructor
    AquiferInterface(int aqID,
                     const std::vector<Aquancon::AquancCell>& connections,
                     const Simulator& ebosSimulator)
        : aquiferID_(aqID)
        , connections_(connections)
        , ebos_simulator_(ebosSimulator)
    {
    }

    // Destructor
    virtual ~AquiferInterface()
    {
    }

    void initFromRestart(const data::Aquifers& aquiferSoln)
    {
        auto xaqPos = aquiferSoln.find(this->aquiferID());
        if (xaqPos == aquiferSoln.end())
            return;

        this->assignRestartData(xaqPos->second);

        this->W_flux_ = xaqPos->second.volume;
        this->pa0_ = xaqPos->second.initPressure;
        this->solution_set_from_restart_ = true;
    }

    void initialSolutionApplied()
    {
        initQuantities();
    }

    void beginTimeStep()
    {
        ElementContext elemCtx(ebos_simulator_);
        auto elemIt = ebos_simulator_.gridView().template begin<0>();
        const auto& elemEndIt = ebos_simulator_.gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);

            const int cellIdx = elemCtx.globalSpaceIndex(0, 0);
            const int idx = cellToConnectionIdx_[cellIdx];
            if (idx < 0)
                continue;

            elemCtx.updateIntensiveQuantities(0);
            const auto& iq = elemCtx.intensiveQuantities(0, 0);
            pressure_previous_[idx] = getValue(iq.fluidState().pressure(waterPhaseIdx));
        }
    }

    template <class Context>
    void addToSource(RateVector& rates,
                     const Context& context,
                     const unsigned spaceIdx,
                     const unsigned timeIdx)
    {
        const unsigned cellIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

        const int idx = this->cellToConnectionIdx_[cellIdx];
        if (idx < 0)
            return;

        // We are dereferencing the value of IntensiveQuantities because
        // cachedIntensiveQuantities return a const pointer to
        // IntensiveQuantities of that particular cell_id
        const auto& intQuants = context.intensiveQuantities(spaceIdx, timeIdx);

        // This is the pressure at td + dt
        this->updateCellPressure(this->pressure_current_, idx, intQuants);
        this->calculateInflowRate(idx, context.simulator());

        rates[BlackoilIndices::conti0EqIdx + FluidSystem::waterCompIdx]
            += this->Qai_[idx] / context.dofVolume(spaceIdx, timeIdx);
    }

    std::size_t size() const {
        return this->connections_.size();
    }

    int aquiferID() const { return this->aquiferID_; }

protected:
    inline Scalar gravity_() const
    {
        return ebos_simulator_.problem().gravity()[2];
    }

    inline void initQuantities()
    {
        // We reset the cumulative flux at the start of any simulation, so, W_flux = 0
        if (!this->solution_set_from_restart_) {
            W_flux_ = Scalar{0};
        }

        // We next get our connections to the aquifer and initialize these quantities using the initialize_connections
        // function
        initializeConnections();
        calculateAquiferCondition();
        calculateAquiferConstants();

        pressure_previous_.resize(this->connections_.size(), Scalar{0});
        pressure_current_.resize(this->connections_.size(), Scalar{0});
        Qai_.resize(this->connections_.size(), Scalar{0});
    }

    inline void
    updateCellPressure(std::vector<Eval>& pressure_water, const int idx, const IntensiveQuantities& intQuants)
    {
        const auto& fs = intQuants.fluidState();
        pressure_water.at(idx) = fs.pressure(waterPhaseIdx);
    }

    inline void
    updateCellPressure(std::vector<Scalar>& pressure_water, const int idx, const IntensiveQuantities& intQuants)
    {
        const auto& fs = intQuants.fluidState();
        pressure_water.at(idx) = fs.pressure(waterPhaseIdx).value();
    }

    virtual void endTimeStep() = 0;

    const int aquiferID_{};
    const std::vector<Aquancon::AquancCell> connections_;
    const Simulator& ebos_simulator_;

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
    Scalar rhow_{};

    Eval W_flux_;

    bool solution_set_from_restart_ {false};

    void initializeConnections()
    {
        this->cell_depth_.resize(this->size(), this->aquiferDepth());
        this->alphai_.resize(this->size(), 1.0);
        this->faceArea_connected_.resize(this->size(), Scalar{0});

        // Translate the C face tag into the enum used by opm-parser's TransMult class
        FaceDir::DirEnum faceDirection;

        bool has_active_connection_on_proc = false;

        // denom_face_areas is the sum of the areas connected to an aquifer
        Scalar denom_face_areas{0};
        this->cellToConnectionIdx_.resize(this->ebos_simulator_.gridView().size(/*codim=*/0), -1);
        const auto& gridView = this->ebos_simulator_.vanguard().gridView();
        for (std::size_t idx = 0; idx < this->size(); ++idx) {
            const auto global_index = this->connections_[idx].global_index;
            const int cell_index = this->ebos_simulator_.vanguard().compressedIndex(global_index);
            auto elemIt = gridView.template begin</*codim=*/ 0>();
            if (cell_index > 0)
                std::advance(elemIt, cell_index);

           //the global_index is not part of this grid
            if ( cell_index < 0 || elemIt->partitionType() != Dune::InteriorEntity)
                continue;

            has_active_connection_on_proc = true;

            this->cellToConnectionIdx_[cell_index] = idx;
            this->cell_depth_.at(idx) = this->ebos_simulator_.vanguard().cellCenterDepth(cell_index);
        }
        // get areas for all connections
        ElementMapper elemMapper(gridView, Dune::mcmgElementLayout());
        auto elemIt = gridView.template begin</*codim=*/ 0>();
        const auto& elemEndIt = gridView.template end</*codim=*/ 0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            unsigned cell_index = elemMapper.index(elem);
            int idx = this->cellToConnectionIdx_[cell_index];

            // only deal with connections given by the aquifer
            if( idx < 0)
                continue;

            auto isIt = gridView.ibegin(elem);
            const auto& isEndIt = gridView.iend(elem);
            for (; isIt != isEndIt; ++ isIt) {
                // store intersection, this might be costly
                const auto& intersection = *isIt;

                // only deal with grid boundaries
                if (!intersection.boundary())
                    continue;

                int insideFaceIdx  = intersection.indexInInside();
                switch (insideFaceIdx) {
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
            denom_face_areas += this->faceArea_connected_.at(idx);
        }

        const auto& comm = this->ebos_simulator_.vanguard().grid().comm();
        comm.sum(&denom_face_areas, 1);
        const double eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());
        for (std::size_t idx = 0; idx < this->size(); ++idx) {
            // Protect against division by zero NaNs.
            this->alphai_.at(idx) = (denom_face_areas < eps_sqrt)
                ? Scalar{0}
                : this->faceArea_connected_.at(idx) / denom_face_areas;
        }

        if (this->solution_set_from_restart_) {
            this->rescaleProducedVolume(has_active_connection_on_proc);
        }
    }

    void rescaleProducedVolume(const bool has_active_connection_on_proc)
    {
        // Needed in parallel restart to approximate influence of aquifer
        // being "owned" by a subset of the parallel processes.  If the
        // aquifer is fully owned by a single process--i.e., if all cells
        // connecting to the aquifer are on a single process--then this_area
        // is tot_area on that process and zero elsewhere.

        const auto this_area = has_active_connection_on_proc
            ? std::accumulate(this->alphai_.begin(),
                              this->alphai_.end(),
                              Scalar{0})
            : Scalar{0};

        const auto tot_area = this->ebos_simulator_.vanguard()
            .grid().comm().sum(this_area);

        this->W_flux_ *= this_area / tot_area;
    }

    virtual void assignRestartData(const data::AquiferData& xaq) = 0;

    virtual void calculateInflowRate(int idx, const Simulator& simulator) = 0;

    virtual void calculateAquiferCondition() = 0;

    virtual void calculateAquiferConstants() = 0;

    virtual Scalar aquiferDepth() const = 0;

    // This function is for calculating the aquifer properties from equilibrium state with the reservoir
    virtual Scalar calculateReservoirEquilibrium()
    {
        // Since the global_indices are the reservoir index, we just need to extract the fluidstate at those indices
        std::vector<Scalar> pw_aquifer;
        Scalar water_pressure_reservoir;

        ElementContext elemCtx(this->ebos_simulator_);
        const auto& gridView = this->ebos_simulator_.gridView();
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            elemCtx.updatePrimaryStencil(elem);

            const auto cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto idx = this->cellToConnectionIdx_[cellIdx];
            if (idx < 0)
                continue;

            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            const auto& iq0 = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = iq0.fluidState();

            water_pressure_reservoir = fs.pressure(waterPhaseIdx).value();
            const auto water_density = fs.density(waterPhaseIdx);

            const auto gdz =
                this->gravity_() * (this->cell_depth_[idx] - this->aquiferDepth());

            pw_aquifer.push_back(this->alphai_[idx] *
                (water_pressure_reservoir - water_density.value()*gdz));
        }

        // We take the average of the calculated equilibrium pressures.
        const auto& comm = ebos_simulator_.vanguard().grid().comm();

        Scalar vals[2];
        vals[0] = std::accumulate(this->alphai_.begin(), this->alphai_.end(), Scalar{0});
        vals[1] = std::accumulate(pw_aquifer.begin(), pw_aquifer.end(), Scalar{0});

        comm.sum(vals, 2);

        return vals[1] / vals[0];
    }

    // This function is used to initialize and calculate the alpha_i for each grid connection to the aquifer
};
} // namespace Opm
#endif
