/*
  Copyright (C) 2020 Equinor ASA
  Copyright (C) 2020 SINTEF Digital

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

#ifndef OPM_AQUIFERNUMERICAL_HEADER_INCLUDED
#define OPM_AQUIFERNUMERICAL_HEADER_INCLUDED

#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/SingleNumericalAquifer.hpp>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/densead/Evaluation.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <opm/simulators/aquifers/AquiferInterface.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

namespace Opm
{
template <typename TypeTag>
class AquiferNumerical : public AquiferInterface<TypeTag>
{
public:
    using BlackoilIndices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = FluidSystem::numPhases };
    static constexpr int numEq = BlackoilIndices::numEq;

    using Eval =  DenseAd::Evaluation<double, numEq>;
    using Toolbox = MathToolbox<Eval>;

    using typename AquiferInterface<TypeTag>::RateVector;

    // Constructor
    AquiferNumerical(const SingleNumericalAquifer& aquifer,
                     const Simulator& ebos_simulator)
        : AquiferInterface<TypeTag>(aquifer.id(), ebos_simulator)
        , flux_rate_      (0.0)
        , cumulative_flux_(0.0)
        , init_pressure_  (aquifer.numCells(), 0.0)
    {
        this->cell_to_aquifer_cell_idx_.resize(this->ebos_simulator_.gridView().size(/*codim=*/0), -1);

        auto aquifer_on_process = false;
        for (std::size_t idx = 0; idx < aquifer.numCells(); ++idx) {
            const auto* cell = aquifer.getCellPrt(idx);

            // Due to parallelisation, the cell might not exist in the current process
            const int compressed_idx = ebos_simulator.vanguard().compressedIndexForInterior(cell->global_index);
            if (compressed_idx >= 0) {
                this->cell_to_aquifer_cell_idx_[compressed_idx] = idx;
                aquifer_on_process = true;
            }
        }

        if (aquifer_on_process) {
            this->checkConnectsToReservoir();
        }
    }

    static AquiferNumerical serializationTestObject(const Simulator& ebos_simulator)
    {
        AquiferNumerical result({}, ebos_simulator);
        result.flux_rate_ = 1.0;
        result.cumulative_flux_ = 2.0;
        result.init_pressure_ = {3.0, 4.0};
        result.pressure_ = 5.0;

        return result;
    }

    void initFromRestart(const data::Aquifers& aquiferSoln) override
    {
        auto xaqPos = aquiferSoln.find(this->aquiferID());
        if (xaqPos == aquiferSoln.end())
            return;

        if (this->connects_to_reservoir_) {
            this->cumulative_flux_ = xaqPos->second.volume;
        }

        if (const auto* aqData = xaqPos->second.typeData.template get<data::AquiferType::Numerical>();
            aqData != nullptr)
        {
            this->init_pressure_ = aqData->initPressure;
        }

        this->solution_set_from_restart_ = true;
    }

    void beginTimeStep() override {}
    void addToSource(RateVector&, const unsigned, const unsigned) override {}

    void endTimeStep() override
    {
        this->pressure_ = this->calculateAquiferPressure();
        this->flux_rate_ = this->calculateAquiferFluxRate();
        this->cumulative_flux_ += this->flux_rate_ * this->ebos_simulator_.timeStepSize();
    }

    data::AquiferData aquiferData() const override
    {
        data::AquiferData data;
        data.aquiferID = this->aquiferID();
        data.pressure = this->pressure_;
        data.fluxRate = this->flux_rate_;
        data.volume = this->cumulative_flux_;

        auto* aquNum = data.typeData.template create<data::AquiferType::Numerical>();
        aquNum->initPressure = this->init_pressure_;

        return data;
    }

    void initialSolutionApplied() override
    {
        if (this->solution_set_from_restart_) {
            return;
        }

        this->pressure_ = this->calculateAquiferPressure(this->init_pressure_);
        this->flux_rate_ = 0.;
        this->cumulative_flux_ = 0.;
    }

    void computeFaceAreaFraction(const std::vector<double>& /*total_face_area*/) override
    {}

    double totalFaceArea() const override
    {
        return 1.0;
    }

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(flux_rate_);
        serializer(cumulative_flux_);
        serializer(init_pressure_);
        serializer(pressure_);
    }

    bool operator==(const AquiferNumerical& rhs) const
    {
        return this->flux_rate_ == rhs.flux_rate_ &&
               this->cumulative_flux_ == rhs.cumulative_flux_ &&
               this->init_pressure_ == rhs.init_pressure_ &&
               this->pressure_ == rhs.pressure_;
    }

    double cumulativeFlux() const
    {
        return this->cumulative_flux_;
    }

private:
    void checkConnectsToReservoir()
    {
        ElementContext elem_ctx(this->ebos_simulator_);
        auto elemIt = std::find_if(this->ebos_simulator_.gridView().template begin</*codim=*/0>(),
                                   this->ebos_simulator_.gridView().template end</*codim=*/0>(),
            [&elem_ctx, this](const auto& elem) -> bool
        {
            elem_ctx.updateStencil(elem);

            const auto cell_index = elem_ctx
                .globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);

            return this->cell_to_aquifer_cell_idx_[cell_index] == 0;
        });

        assert ((elemIt != this->ebos_simulator_.gridView().template end</*codim=*/0>())
                && "Internal error locating numerical aquifer's connecting cell");

        this->connects_to_reservoir_ =
            elemIt->partitionType() == Dune::InteriorEntity;
    }

    double calculateAquiferPressure() const
    {
        auto capture = std::vector<double>(this->init_pressure_.size(), 0.0);
        return this->calculateAquiferPressure(capture);
    }

    double calculateAquiferPressure(std::vector<double>& cell_pressure) const
    {
        double sum_pressure_watervolume = 0.;
        double sum_watervolume = 0.;

        ElementContext  elem_ctx(this->ebos_simulator_);
        const auto& gridView = this->ebos_simulator_.gridView();
        OPM_BEGIN_PARALLEL_TRY_CATCH();

        for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
            elem_ctx.updatePrimaryStencil(elem);

            const std::size_t cell_index = elem_ctx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const int idx = this->cell_to_aquifer_cell_idx_[cell_index];
            if (idx < 0) {
                continue;
            }

            elem_ctx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            const auto& iq0 = elem_ctx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = iq0.fluidState();

            // TODO: the porosity of the cells are still wrong for numerical aquifer cells
            // Because the dofVolume still based on the grid information.
            // The pore volume is correct. Extra efforts will be done to get sensible porosity value here later.
            const double water_saturation = fs.saturation(this->phaseIdx_()).value();
            const double porosity = iq0.porosity().value();
            const double volume = elem_ctx.dofTotalVolume(0, 0);
            // TODO: not sure we should use water pressure here
            const double water_pressure_reservoir = fs.pressure(this->phaseIdx_()).value();
            const double water_volume = volume * porosity * water_saturation;
            sum_pressure_watervolume += water_volume * water_pressure_reservoir;
            sum_watervolume += water_volume;

            cell_pressure[idx] = water_pressure_reservoir;
        }
        OPM_END_PARALLEL_TRY_CATCH("AquiferNumerical::calculateAquiferPressure() failed: ", this->ebos_simulator_.vanguard().grid().comm());
        const auto& comm = this->ebos_simulator_.vanguard().grid().comm();
        comm.sum(&sum_pressure_watervolume, 1);
        comm.sum(&sum_watervolume, 1);

        // Ensure all processes have same notion of the aquifer cells' pressure values.
        comm.sum(cell_pressure.data(), cell_pressure.size());

        return sum_pressure_watervolume / sum_watervolume;
    }

    template <class ElemCtx>
    double getWaterFlux(const ElemCtx& elem_ctx, unsigned face_idx) const
    {
        const auto& exQuants = elem_ctx.extensiveQuantities(face_idx, /*timeIdx*/ 0);
        const double water_flux = Toolbox::value(exQuants.volumeFlux(this->phaseIdx_()));
        return water_flux;
    }

    double calculateAquiferFluxRate() const
    {
        double aquifer_flux = 0.0;

        if (! this->connects_to_reservoir_) {
            return aquifer_flux;
        }

        ElementContext elem_ctx(this->ebos_simulator_);
        const auto& gridView = this->ebos_simulator_.gridView();
        for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
            // elem_ctx.updatePrimaryStencil(elem);
            elem_ctx.updateStencil(elem);

            const std::size_t cell_index = elem_ctx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const int idx = this->cell_to_aquifer_cell_idx_[cell_index];
            // we only need the first aquifer cell
            if (idx != 0) {
                continue;
            }

            const std::size_t num_interior_faces = elem_ctx.numInteriorFaces(/*timeIdx*/ 0);
            // const auto &problem = elem_ctx.problem();
            const auto& stencil = elem_ctx.stencil(0);
            // const auto& inQuants = elem_ctx.intensiveQuantities(0, /*timeIdx*/ 0);

            for (std::size_t face_idx = 0; face_idx < num_interior_faces; ++face_idx) {
                const auto& face = stencil.interiorFace(face_idx);
                // dof index
                const std::size_t i = face.interiorIndex();
                const std::size_t j = face.exteriorIndex();
                // compressed index
                // const std::size_t I = stencil.globalSpaceIndex(i);
                const std::size_t J = stencil.globalSpaceIndex(j);

                assert(stencil.globalSpaceIndex(i) == cell_index);

                // we do not consider the flux within aquifer cells
                // we only need the flux to the connections
                if (this->cell_to_aquifer_cell_idx_[J] > 0) {
                    continue;
                }
                elem_ctx.updateAllIntensiveQuantities();
                elem_ctx.updateAllExtensiveQuantities();

                const double water_flux = getWaterFlux(elem_ctx,face_idx);
                const std::size_t up_id = water_flux >= 0.0 ? i : j;
                const auto& intQuantsIn = elem_ctx.intensiveQuantities(up_id, 0);
                const double invB = Toolbox::value(intQuantsIn.fluidState().invB(this->phaseIdx_()));
                const double face_area = face.area();
                aquifer_flux += water_flux * invB * face_area;
            }

            // we only need to handle the first aquifer cell, we can exit loop here
            break;
        }

        return aquifer_flux;
    }

    double flux_rate_; // aquifer influx rate
    double cumulative_flux_; // cumulative aquifer influx
    std::vector<double> init_pressure_{};
    double pressure_; // aquifer pressure
    bool solution_set_from_restart_ {false};
    bool connects_to_reservoir_ {false};

    // TODO: maybe unordered_map can also do the work to save memory?
    std::vector<int> cell_to_aquifer_cell_idx_;
};

} // namespace Opm

#endif
