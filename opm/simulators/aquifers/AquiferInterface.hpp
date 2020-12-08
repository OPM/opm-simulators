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
#include <opm/parser/eclipse/EclipseState/Aquancon.hpp>
#include <opm/parser/eclipse/EclipseState/AquiferCT.hpp>
#include <opm/parser/eclipse/EclipseState/Aquifetp.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Evaluation.hpp>
#include <opm/material/densead/Math.hpp>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>

#include <algorithm>
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

    typedef Opm::BlackOilFluidState<Eval,
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
        : aquiferID(aqID)
        , connections_(connections)
        , ebos_simulator_(ebosSimulator)
    {
    }

    // Deconstructor
    virtual ~AquiferInterface()
    {
    }

    void initFromRestart(const std::vector<data::AquiferData>& aquiferSoln)
    {
        auto xaqPos
            = std::find_if(aquiferSoln.begin(), aquiferSoln.end(), [this](const data::AquiferData& xaq) -> bool {
                   return xaq.aquiferID == this->aquiferID;
              });

        if (xaqPos == aquiferSoln.end())
            return;

        this->assignRestartData(*xaqPos);
        this->W_flux_ = xaqPos->volume;
        this->pa0_ = xaqPos->initPressure;
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

            int cellIdx = elemCtx.globalSpaceIndex(0, 0);
            int idx = cellToConnectionIdx_[cellIdx];
            if (idx < 0)
                continue;

            elemCtx.updateIntensiveQuantities(0);
            const auto& iq = elemCtx.intensiveQuantities(0, 0);
            pressure_previous_[idx] = Opm::getValue(iq.fluidState().pressure(waterPhaseIdx));
        }
    }

    template <class Context>
    void addToSource(RateVector& rates, const Context& context, unsigned spaceIdx, unsigned timeIdx)
    {
        unsigned cellIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

        int idx = cellToConnectionIdx_[cellIdx];
        if (idx < 0)
            return;

        // We are dereferencing the value of IntensiveQuantities because cachedIntensiveQuantities return a const
        // pointer to IntensiveQuantities of that particular cell_id
        const IntensiveQuantities intQuants = context.intensiveQuantities(spaceIdx, timeIdx);
        // This is the pressure at td + dt
        updateCellPressure(pressure_current_, idx, intQuants);
        updateCellDensity(idx, intQuants);
        calculateInflowRate(idx, context.simulator());
        rates[BlackoilIndices::conti0EqIdx + FluidSystem::waterCompIdx]
            += Qai_[idx] / context.dofVolume(spaceIdx, timeIdx);
    }


    std::size_t size() const {
        return this->connections_.size();
    }


protected:
    inline Scalar gravity_() const
    {
        return ebos_simulator_.problem().gravity()[2];
    }

    inline void initQuantities()
    {
        // We reset the cumulative flux at the start of any simulation, so, W_flux = 0
        if (!this->solution_set_from_restart_) {
            W_flux_ = 0.;
        }

        // We next get our connections to the aquifer and initialize these quantities using the initialize_connections
        // function
        initializeConnections();
        calculateAquiferCondition();
        calculateAquiferConstants();

        pressure_previous_.resize(this->connections_.size(), 0.);
        pressure_current_.resize(this->connections_.size(), 0.);
        Qai_.resize(this->connections_.size(), 0.0);
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

    inline void updateCellDensity(const int idx, const IntensiveQuantities& intQuants)
    {
        const auto& fs = intQuants.fluidState();
        rhow_.at(idx) = fs.density(waterPhaseIdx);
    }

    template <class Intersection>
    inline double getFaceArea(const Intersection& intersection,
                              unsigned idx) const
    {
        const auto& geometry = intersection.geometry();
        const auto defaultFaceArea = geometry.volume();
        return (!this->connections_[idx].influx_coeff.first) ? defaultFaceArea : this->connections_[idx].influx_coeff.second;
    }

    virtual void endTimeStep() = 0;

    const int aquiferID;
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
    std::vector<Eval> rhow_;
    std::vector<Scalar> alphai_;

    Scalar Tc_; // Time constant
    Scalar pa0_; // initial aquifer pressure

    Eval W_flux_;

    bool solution_set_from_restart_ {false};

    virtual void initializeConnections() = 0;

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

            size_t cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            int idx = this->cellToConnectionIdx_[cellIdx];
            if (idx < 0)
                continue;

            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            const auto& iq0 = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = iq0.fluidState();

            water_pressure_reservoir = fs.pressure(waterPhaseIdx).value();
            this->rhow_[idx] = fs.density(waterPhaseIdx);
            pw_aquifer.push_back(
                (water_pressure_reservoir
                 - this->rhow_[idx].value() * this->gravity_() * (this->cell_depth_[idx] - this->aquiferDepth()))
                * this->alphai_[idx]);
        }

        // We take the average of the calculated equilibrium pressures.
        const Scalar sum_alpha = std::accumulate(this->alphai_.begin(), this->alphai_.end(), 0.);
        const Scalar aquifer_pres_avg = std::accumulate(pw_aquifer.begin(), pw_aquifer.end(), 0.) / sum_alpha;
        return aquifer_pres_avg;
    }

    // This function is used to initialize and calculate the alpha_i for each grid connection to the aquifer
};
} // namespace Opm
#endif
