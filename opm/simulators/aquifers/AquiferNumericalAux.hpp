/*
  Copyright (C) 2026 SINTEF Digital

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

#ifndef OPM_AQUIFERNUMERICALAUX_HEADER_INCLUDED
#define OPM_AQUIFERNUMERICALAUX_HEADER_INCLUDED

#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/nonlinear/newtonmethodproperties.hh>

#include <opm/material/common/MathToolbox.hpp>

#include <opm/output/data/Aquifer.hpp>

#include <opm/simulators/aquifers/AquiferInterface.hpp>
#include <opm/simulators/aquifers/NumericalAquiferAuxCells.hpp>

#include <cstddef>
#include <vector>

namespace Opm {

/*!
 * \brief Reporting for a numerical aquifer represented as auxiliary cells.
 *
 * The counterpart of AquiferNumerical, which finds its cells in the grid and reads their
 * fluxes off the element context.  Neither is possible here -- the aquifer is not in the
 * grid -- so both quantities are formed from the degrees of freedom directly:
 *
 *  - the aquifer pressure is the water-volume-weighted water pressure of its cells, the
 *    same average the grid-cell representation forms over the same cells;
 *  - the influx is the water flux across the connections that reach into the reservoir,
 *    computed with the very local residual that assembled them, so the number reported is
 *    the number the equations used rather than a second opinion about it.
 *
 * The intra-aquifer chain is deliberately excluded: what the aquifer reports is what it
 * gives the reservoir, not what moves inside itself.
 */
template <typename TypeTag>
class AquiferNumericalAux : public AquiferInterface<TypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using LocalResidual = GetPropType<TypeTag, Properties::LocalResidual>;
    using Linearizer = GetPropType<TypeTag, Properties::Linearizer>;
    using AuxCells = NumericalAquiferAuxCells<TypeTag>;

    using typename AquiferInterface<TypeTag>::RateVector;

    static constexpr bool conserveSurfaceVolume =
        getPropValue<TypeTag, Properties::BlackoilConserveSurfaceVolume>();

public:
    AquiferNumericalAux(const int aquiferId,
                        const AuxCells& cells,
                        const Simulator& simulator)
        : AquiferInterface<TypeTag>(aquiferId, simulator)
        , cells_(cells)
        , dofs_(cells.aquiferDofs(aquiferId))
        , connections_(cells.reservoirConnections(aquiferId))
        , init_pressure_(cells.initialPressure(aquiferId))
    {}

    void initFromRestart(const data::Aquifers&) override
    {
        // Restart together with auxiliary cells is refused at setup; there is nothing
        // sensible to restore here and pretending otherwise would hide that.
    }

    void initialSolutionApplied() override
    {
        this->pressure_ = this->aquiferPressure();
        this->flux_rate_ = 0.0;
        this->cumulative_flux_ = 0.0;
    }

    void beginTimeStep() override {}

    void endTimeStep() override
    {
        this->pressure_ = this->aquiferPressure();
        this->flux_rate_ = this->aquiferFluxRate();
        this->cumulative_flux_ += this->flux_rate_ * this->simulator_.timeStepSize();
    }

    void addToSource(RateVector&, const unsigned, const unsigned) override
    {
        // The aquifer is a set of degrees of freedom with their own equations, not a
        // source term on a reservoir cell.
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

    void computeFaceAreaFraction(const std::vector<Scalar>&) override {}

    Scalar totalFaceArea() const override
    { return 1.0; }

    //! Pore volume of the aquifer, for the field totals it would otherwise drop out of.
    Scalar poreVolume() const
    {
        const auto& model = this->simulator_.model();

        Scalar pv = 0.0;
        for (const auto dof : this->dofs_) {
            const auto& iq = model.intensiveQuantities(dof, /*timeIdx=*/0);
            pv += model.dofTotalVolume(dof) * getValue(iq.porosity());
        }

        return pv;
    }

    Scalar cumulativeFlux() const
    { return this->cumulative_flux_; }

private:
    Scalar aquiferPressure() const
    {
        const auto& model = this->simulator_.model();
        const auto waterPos = this->phaseIdx_();

        Scalar sumPressureWaterVolume = 0.0;
        Scalar sumWaterVolume = 0.0;

        for (const auto dof : this->dofs_) {
            const auto& iq = model.intensiveQuantities(dof, /*timeIdx=*/0);
            const auto& fs = iq.fluidState();

            const Scalar waterVolume = model.dofTotalVolume(dof)
                * getValue(iq.porosity())
                * getValue(fs.saturation(waterPos));

            sumPressureWaterVolume += waterVolume * getValue(fs.pressure(waterPos));
            sumWaterVolume += waterVolume;
        }

        return (sumWaterVolume > 0.0)
            ? sumPressureWaterVolume / sumWaterVolume
            : Scalar{0};
    }

    /*!
     * \brief Surface water rate from the aquifer into the reservoir.
     *
     * Formed with the same LocalResidual::computeFlux() the linearizer used, given the
     * same neighbour information it cached, so the reported influx is the flux the
     * equations were assembled with -- upwinding, gravity, threshold pressures and all --
     * rather than a reimplementation that would drift from it.
     */
    Scalar aquiferFluxRate() const
    {
        if constexpr (! Linearizer::assemblesAuxiliaryDofEquations) {
            // Cannot happen: auxiliary DOFs carrying equations are refused on such a
            // linearizer when they are registered.
            return Scalar{0};
        }
        else {
            const auto& model = this->simulator_.model();
            const auto& problem = this->simulator_.problem();
            const auto& neighborInfo = model.linearizer().getNeighborInfo();
            const auto waterPos = this->phaseIdx_();

            Scalar rate = 0.0;
            for (const auto& [aquiferDof, reservoirDof] : this->connections_) {
                for (const auto& nbInfo : neighborInfo[aquiferDof]) {
                    if (nbInfo.neighbor != reservoirDof) {
                        continue;
                    }

                    const auto& iqIn = model.intensiveQuantities(aquiferDof, 0);
                    const auto& iqEx = model.intensiveQuantities(reservoirDof, 0);

                    RateVector flux(0.0);
                    RateVector darcy(0.0);
                    LocalResidual::computeFlux(flux, darcy,
                                               aquiferDof, reservoirDof,
                                               iqIn, iqEx,
                                               nbInfo.res_nbinfo,
                                               problem.moduleParams());

                    const auto& fsys = iqIn.fluidState().fluidSystem();
                    const auto waterEqIdx = Indices::conti0EqIdx
                        + fsys.canonicalToActiveCompIdx(fsys.solventComponentIndex(waterPos));

                    // computeFlux() reports a flux per unit area, which the linearizer
                    // scales by the face area; an authored connection carries the whole
                    // geometry in its transmissibility, so that area is unity.
                    Scalar connRate = getValue(flux[waterEqIdx]) * nbInfo.res_nbinfo.faceArea;

                    // ... and it is written in mass unless the model conserves surface
                    // volume, while the aquifer influx is reported as a surface rate.
                    if constexpr (! conserveSurfaceVolume) {
                        connRate /= fsys.referenceDensity(waterPos,
                                                          problem.pvtRegionIndex(aquiferDof));
                    }

                    rate += connRate;
                    break;
                }
            }

            return rate;
        }
    }

    const AuxCells& cells_;
    std::vector<unsigned> dofs_{};
    std::vector<std::pair<unsigned, unsigned>> connections_{};
    std::vector<Scalar> init_pressure_{};

    Scalar pressure_{0.0};
    Scalar flux_rate_{0.0};
    Scalar cumulative_flux_{0.0};
};

} // namespace Opm

#endif // OPM_AQUIFERNUMERICALAUX_HEADER_INCLUDED
