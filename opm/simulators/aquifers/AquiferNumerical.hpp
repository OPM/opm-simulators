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

#include <opm/output/data/Aquifer.hpp>

#include <opm/parser/eclipse/EclipseState/NumericalAquifer.hpp>

namespace Opm
{
template <typename TypeTag>
class AquiferNumerical
{
public:
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BlackoilIndices = GetPropType<TypeTag, Properties::Indices>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;

    enum { dimWorld = GridView::dimensionworld };

    static const auto waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static const int numEq = BlackoilIndices::numEq;

    using Eval =  DenseAd::Evaluation<double, numEq>;
    using Toolbox = Opm::MathToolbox<Eval>;

    // Constructor
    AquiferNumerical(const SingleNumericalAquifer& aquifer,
                     const std::unordered_map<int, int>& cartesian_to_compressed,
                     const Simulator& ebos_simulator)
    : aquifer_(aquifer)
    , ebos_simulator_(ebos_simulator)
    , flux_rate_(0.)
    , cumulative_flux_(0.)
    {

        this->cell_to_aquifer_cell_idx_.resize(this->ebos_simulator_.gridView().size(/*codim=*/0), -1);
        const auto& cells = this->aquifer_.cells();

        // TODO: here, the parallisam is obviously ignored
        for (size_t idx = 0; idx < cells.size(); ++idx) {
            const int global_idx = cells[idx].global_index;
            const int cell_idx = cartesian_to_compressed.at(global_idx);
            this->cell_to_aquifer_cell_idx_[cell_idx] = idx;
        }
    }

    void initFromRestart(const std::vector<data::AquiferData>& aquiferSoln)
    {
        // NOT handling Restart for now
    }

    void endTimeStep()
    {
        this->pressure_ = this->calculateAquiferPressure();
        this->flux_rate_ = this->calculateAquiferFluxRate();
        this->cumulative_flux_ += this->flux_rate_ * this->ebos_simulator_.timeStepSize();
    }

    Opm::data::AquiferData aquiferData() const
    {
        data::AquiferData data;
        data.aquiferID = this->aquifer_.id();
        data.initPressure = this->init_pressure_;
        data.pressure = this->pressure_;
        data.fluxRate = this->flux_rate_;
        data.volume = this->cumulative_flux_;
        data.type = Opm::data::AquiferType::Numerical;
        return data;
    }

    void initialSolutionApplied()
    {
        this->init_pressure_ = this->calculateAquiferPressure();
        this->pressure_ = this->init_pressure_;
        this->flux_rate_ = 0.;
        this->cumulative_flux_ = 0.;
    }

private:
    const Opm::SingleNumericalAquifer& aquifer_;
    const Simulator& ebos_simulator_;
    double flux_rate_; // aquifer influx rate
    double cumulative_flux_; // cumulative aquifer influx
    double init_pressure_;
    double pressure_; // aquifer pressure

    // TODO: maybe unordered_map can also do the work to save memory?
    std::vector<int> cell_to_aquifer_cell_idx_;

    double calculateAquiferPressure() const
    {
        double sum_pv_pressure = 0.;
        double sum_pv = 0.;

        ElementContext  elem_ctx(this->ebos_simulator_);
        const auto& gridView = this->ebos_simulator_.gridView();
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            elem_ctx.updatePrimaryStencil(elem);

            const size_t cell_index = elem_ctx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const int idx = this->cell_to_aquifer_cell_idx_[cell_index];
            if (idx < 0) {
                continue;
            }
            elem_ctx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            const auto& iq0 = elem_ctx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = iq0.fluidState();

            const double water_pressure_reservoir = fs.pressure(waterPhaseIdx).value();
            const auto& pvs = this->ebos_simulator_.vanguard().eclState().fieldProps().porv(true);
            // TODO: should get this PV, how to consider the rock compressiblity
            const double pv = pvs[cell_index];
            sum_pv_pressure += pv * water_pressure_reservoir;
            sum_pv += pv;
        }

        assert(sum_pv > 0.);
        return sum_pv_pressure/ sum_pv;
    }

    double calculateAquiferFluxRate() const
    {
        double aquifer_flux = 0.;

        ElementContext  elem_ctx(this->ebos_simulator_);
        const auto& gridView = this->ebos_simulator_.gridView();
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto &elem = *elemIt;
            // elem_ctx.updatePrimaryStencil(elem);
            elem_ctx.updateStencil(elem);

            const size_t cell_index = elem_ctx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const int idx = this->cell_to_aquifer_cell_idx_[cell_index];
            // we only need the first aquifer cell
            if (idx != 0) {
                continue;
            }
            elem_ctx.updateAllIntensiveQuantities();
            elem_ctx.updateAllExtensiveQuantities();

            const size_t num_interior_faces = elem_ctx.numInteriorFaces(/*timeIdx*/ 0);
            const auto &problem = elem_ctx.problem();
            const auto &stencil = elem_ctx.stencil(0);
            // const auto& inQuants = elem_ctx.intensiveQuantities(0, /*timeIdx*/ 0);

            for (size_t face_idx = 0; face_idx < num_interior_faces; ++face_idx) {
                const auto &face = stencil.interiorFace(face_idx);
                // dof index
                const size_t i = face.interiorIndex();
                const size_t j = face.exteriorIndex();
                // compressed index
                const size_t I = stencil.globalSpaceIndex(i);
                const size_t J = stencil.globalSpaceIndex(j);

                assert(I == cell_index);

                // we do not consider the flux within aquifer cells
                // we only need the flux to the connections
                if (this->cell_to_aquifer_cell_idx_[J] > 0) {
                    continue;
                }
                const auto &exQuants = elem_ctx.extensiveQuantities(face_idx, /*timeIdx*/ 0);
                const double water_flux = Toolbox::value(exQuants.volumeFlux(waterPhaseIdx));

                const size_t up_id = water_flux >= 0. ? i : j;
                const auto &intQuantsIn = elem_ctx.intensiveQuantities(up_id, 0);
                const double invB = Toolbox::value(intQuantsIn.fluidState().invB(waterPhaseIdx));
                const double face_area = face.area();
                aquifer_flux += water_flux * invB * face_area;
            }

            // we only need to handle the first aquifer cell, we can exit loop here
            break;
        }

        return aquifer_flux;
    }
};
} // namespace Opm
#endif
