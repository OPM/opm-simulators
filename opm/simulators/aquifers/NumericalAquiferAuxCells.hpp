// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::NumericalAquiferAuxCells
 */
#ifndef OPM_NUMERICAL_AQUIFER_AUX_CELLS_HPP
#define OPM_NUMERICAL_AQUIFER_AUX_CELLS_HPP

#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/simulators/flow/FlowAuxCellModule.hpp>

#include <fmt/format.h>

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace Opm {

/*!
 * \brief Numerical aquifers (AQUNUM/AQUCON) represented as auxiliary cells.
 *
 * The alternative to letting each aquifer cell take over a grid cell.  The aquifer
 * satisfies the same flow equations either way, and with the same coefficients: the pore
 * volume, depth and regions come from the AQUNUM record exactly as they would when
 * overriding a grid cell's field properties, and the connections are the very NNCs the
 * grid-cell representation would have generated -- reused here rather than reimplemented,
 * so that the two representations are the same discrete system and can be compared
 * directly.
 *
 * What differs is only where the unknown lives.  The grid keeps the shape the deck gives
 * it, which means no cell is resurrected to host an aquifer, no field property is
 * overridden, and no non-neighbour connection is created.
 */
template <class TypeTag>
class NumericalAquiferAuxCells : public FlowAuxCellModule<TypeTag>
{
    using ParentType = FlowAuxCellModule<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
    using GlobalEqVector = GetPropType<TypeTag, Properties::GlobalEqVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

    enum { dimWorld = GridView::dimensionworld };

public:
    using Connection = typename ParentType::Connection;

    explicit NumericalAquiferAuxCells(Simulator& simulator)
        : simulator_(simulator)
    {
        const auto& eclState = simulator_.vanguard().eclState();
        const auto& aquifers = eclState.aquifer().numericalAquifers();

        // Aquifer cells are enumerated in the order of allAquiferCells(), which is keyed
        // on the cartesian index the AQUNUM record names.  That index identifies the
        // aquifer cell in the deck and in the NNCs below; it is *not* a grid cell here,
        // since in this mode the grid never took it over.
        for (const auto& [globalIndex, cell] : aquifers.allAquiferCells()) {
            this->cells_.push_back(cell);
            this->cartesianToLocal_.emplace(globalIndex, this->cells_.size() - 1);
        }

        // The cells these records name were deactivated when the grid was built, which
        // is what makes room for the aquifer to live outside it.
        this->checkAquiferCellsAreNotInterior();
    }

    unsigned numDofs() const override
    { return static_cast<unsigned>(this->cells_.size()); }

    Scalar poreVolume(unsigned localIdx) const override
    { return this->cells_.at(localIdx)->poreVolume(); }

    Scalar bulkVolume(unsigned localIdx) const override
    { return this->cells_.at(localIdx)->cellVolume(); }

    Scalar depth(unsigned localIdx) const override
    { return this->cells_.at(localIdx)->depth; }

    unsigned pvtRegionIndex(unsigned localIdx) const override
    { return static_cast<unsigned>(this->cells_.at(localIdx)->pvttable) - 1; }

    unsigned satRegionIndex(unsigned localIdx) const override
    { return static_cast<unsigned>(this->cells_.at(localIdx)->sattable) - 1; }

    /*!
     * \brief The reservoir cell this aquifer cell hangs off.
     *
     * Used to start the aquifer from a sensible fluid state.  All AQUCON connections of
     * an aquifer attach to its first cell, so the whole chain is initialised from that
     * cell's first reservoir neighbour.
     */
    unsigned initialisationPartner(unsigned localIdx) const override
    { return this->initialisationPartner_.at(localIdx); }

    void connections(std::vector<Connection>& conns) const override
    { conns.insert(conns.end(), this->connections_.begin(), this->connections_.end()); }

    /*!
     * \brief Build the connection list.
     *
     * Deferred out of the constructor because it needs the grid's cartesian-to-compressed
     * mapping, which is only meaningful once the grid has been distributed.
     */
    void buildConnections()
    {
        const auto& vanguard = simulator_.vanguard();
        const auto& eclState = vanguard.eclState();
        const auto& aquifers = eclState.aquifer().numericalAquifers();

        this->connections_.clear();
        this->initialisationPartner_.assign(this->cells_.size(), 0);

        // Aquifer cell to aquifer cell: the chain within one aquifer.  Both endpoints
        // are auxiliary cells.
        for (const auto& nnc : aquifers.aquiferCellNNCs()) {
            const auto dof1 = this->auxDofOf(nnc.cell1);
            const auto dof2 = this->auxDofOf(nnc.cell2);
            this->connections_.push_back({dof1, dof2, static_cast<Scalar>(nnc.trans), 0.0, 0.0});
        }

        // Aquifer cell to reservoir cell.  The second endpoint is a real grid cell, so it
        // has to be translated into the local compressed numbering -- and skipped when
        // this rank does not own it.
        // NNCdata normalises the order of its two cells, so which endpoint is the
        // aquifer cell is not fixed; identify it by membership rather than by position.
        std::size_t skipped = 0;
        for (const auto& nnc : aquifers.aquiferConnectionNNCs(eclState.getInputGrid(),
                                                              eclState.fieldProps()))
        {
            const bool firstIsAquifer = this->cartesianToLocal_.count(nnc.cell1) > 0;
            const auto aquiferCartesian = firstIsAquifer ? nnc.cell1 : nnc.cell2;
            const auto reservoirCartesian = firstIsAquifer ? nnc.cell2 : nnc.cell1;

            const auto dof1 = this->auxDofOf(aquiferCartesian);
            const int reservoirCell = vanguard.compressedIndexForInterior(reservoirCartesian);
            if (reservoirCell < 0) {
                ++skipped;
                continue;
            }

            const auto dof2 = static_cast<unsigned>(reservoirCell);
            this->connections_.push_back({dof1, dof2, static_cast<Scalar>(nnc.trans), 0.0, 0.0});

            const auto localIdx = this->localOf(aquiferCartesian);
            if (this->initialisationPartner_[localIdx] == 0) {
                this->initialisationPartner_[localIdx] = dof2;
            }
        }

        // An aquifer cell in the middle of a chain has no reservoir connection of its
        // own; initialise it from the same place as the cell that does.
        unsigned fallback = 0;
        for (auto& partner : this->initialisationPartner_) {
            if (partner != 0) { fallback = partner; }
            else              { partner = fallback; }
        }

        if (skipped > 0) {
            OpmLog::debug(fmt::format("Numerical aquifer: {} connection(s) target cells "
                                      "not owned by this rank", skipped));
        }
    }

    /*!
     * \brief Set the initial state of the aquifer cells.
     *
     * They cannot be equilibrated the ordinary way -- that needs the cell's geometry --
     * so each one starts from the state of the reservoir cell it is connected to, with
     * the phase pressures carried to its own depth and the cell filled with water.  That
     * is the same thing the grid-cell representation arranges for its aquifer cells, and
     * an explicit initial pressure on the AQUNUM record overrides it either way.
     */
    void applyInitial() override
    {
        auto& solution = simulator_.model().solution(/*timeIdx=*/0);
        const auto& problem = simulator_.problem();
        const auto gravity = problem.gravity()[dimWorld - 1];

        for (unsigned localIdx = 0; localIdx < this->numDofs(); ++localIdx) {
            const auto globalIdx = static_cast<unsigned>(this->localToGlobalDof(localIdx));
            const auto partner = this->initialisationPartner_.at(localIdx);

            auto fs = problem.initialFluidState(partner);

            const auto waterPos = FluidSystem::waterPhaseIdx;
            const auto rho = getValue(fs.density(waterPos));
            const auto dz = this->depth(localIdx) - problem.dofCenterDepth(partner);

            const auto& cell = *this->cells_.at(localIdx);
            for (unsigned phase = 0; phase < FluidSystem::numPhases; ++phase) {
                if (!FluidSystem::phaseIsActive(phase)) {
                    continue;
                }

                fs.setSaturation(phase, (phase == waterPos) ? 1.0 : 0.0);
                fs.setPressure(phase, cell.init_pressure.has_value()
                               ? cell.init_pressure.value()
                               : getValue(fs.pressure(phase)) + rho * gravity * dz);
            }

            solution[globalIdx].assignNaive(fs);
            solution[globalIdx].setPvtRegionIndex(this->pvtRegionIndex(localIdx));
        }
    }

    void linearize(SparseMatrixAdapter&, GlobalEqVector&) override
    {
        // Nothing to add: the aquifer cells are assembled by the model's own local
        // residual, like grid cells, and they are always active so there are no dormant
        // rows to condition.
    }

private:
    unsigned auxDofOf(std::size_t cartesianIndex) const
    { return static_cast<unsigned>(this->localToGlobalDof(this->localOf(cartesianIndex))); }

    unsigned localOf(std::size_t cartesianIndex) const
    {
        auto pos = this->cartesianToLocal_.find(cartesianIndex);
        if (pos == this->cartesianToLocal_.end()) {
            OPM_THROW(std::logic_error,
                      fmt::format("Numerical aquifer connection names cell {} which is "
                                  "not an aquifer cell", cartesianIndex));
        }

        return static_cast<unsigned>(pos->second);
    }

    /*!
     * \brief Refuse an aquifer cell buried inside the model.
     *
     * Deactivating the cell an AQUNUM record names is only sound where that cell is at
     * the edge of the model.  With live rock both above and below it, removing it opens
     * a hole in the middle of a column: the neighbours lose a connection they would
     * otherwise have, and with PINCH active the grid processing may bridge across the
     * gap and connect them to each other instead.  Neither is what the deck describes,
     * and neither matches what the grid-cell representation does, so the two modes would
     * quietly stop being comparable.
     *
     * Placing an aquifer cell there is questionable modelling to begin with -- a
     * numerical aquifer is meant to hang off the model, not to sit inside it -- so this
     * is refused rather than approximated.  The proper answer is for a numerical aquifer
     * to be defined independently of the grid and then connected to it, which is the
     * direction the auxiliary-cell representation is going; the restriction can be
     * lifted once the aquifer no longer needs a cell to name at all.
     */
    void checkAquiferCellsAreNotInterior() const
    {
        const auto& grid = simulator_.vanguard().eclState().getInputGrid();
        const auto nz = grid.getNZ();

        for (const auto* cell : this->cells_) {
            if ((cell->K == 0) || (cell->K + 1 >= nz)) {
                continue; // at the top or the bottom of the model
            }

            const auto above = grid.getGlobalIndex(cell->I, cell->J, cell->K - 1);
            const auto below = grid.getGlobalIndex(cell->I, cell->J, cell->K + 1);

            if (grid.cellActive(above) && grid.cellActive(below)) {
                OPM_THROW(std::runtime_error,
                          fmt::format("AQUNUM record for aquifer {} names cell "
                                      "({},{},{}), which has active cells both above and "
                                      "below it. Representing numerical aquifers outside "
                                      "the grid removes the cell they name, which would "
                                      "open a hole inside the model; place the aquifer "
                                      "cell at the edge of the model, or run with the "
                                      "grid-cell representation.",
                                      cell->aquifer_id,
                                      cell->I + 1, cell->J + 1, cell->K + 1));
            }
        }
    }

    Simulator& simulator_;

    //! Aquifer cells, in auxiliary-DOF order.
    std::vector<const NumericalAquiferCell*> cells_{};

    //! Cartesian index named by an AQUNUM record -> auxiliary-DOF-local index.
    std::unordered_map<std::size_t, std::size_t> cartesianToLocal_{};

    std::vector<Connection> connections_{};
    std::vector<unsigned> initialisationPartner_{};
};

} // namespace Opm

#endif // OPM_NUMERICAL_AQUIFER_AUX_CELLS_HPP
