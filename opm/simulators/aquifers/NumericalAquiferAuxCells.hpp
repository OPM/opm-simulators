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

#include <algorithm>
#include <cstddef>
#include <map>
#include <utility>
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

        // Enumerate aquifer by aquifer and, within an aquifer, in the order its cells
        // were declared.  Not allAquiferCells(), which is an unordered map: the degree of
        // freedom numbering would then depend on the hash order, and the chain structure
        // -- which cell of an aquifer carries the reservoir connections -- would be lost.
        for (const auto& [id, aquifer] : aquifers.aquifers()) {
            const auto first = this->cells_.size();

            for (std::size_t i = 0; i < aquifer.numCells(); ++i) {
                const auto* cell = aquifer.getCellPrt(i);
                this->cells_.push_back(cell);
                this->cartesianToLocal_.emplace(cell->global_index, this->cells_.size() - 1);
            }

            this->aquiferRange_.emplace(id, std::make_pair(first, this->cells_.size()));
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
        this->hasReservoirConnection_.assign(this->cells_.size(), false);

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
            if (!this->hasReservoirConnection_[localIdx]) {
                this->initialisationPartner_[localIdx] = dof2;
                this->hasReservoirConnection_[localIdx] = true;
            }
        }

        // Only the cell that carries the AQUCON connections has a reservoir neighbour of
        // its own; the rest of the chain hangs off it.  Give them all that cell's
        // neighbour to start from -- per aquifer, so that two aquifers cannot borrow each
        // other's.  Note this is a fallback for initialisation only: it says where a cell
        // takes its initial state from, not what it is connected to.
        for (const auto& [id, range] : this->aquiferRange_) {
            unsigned connected = 0;
            bool found = false;
            for (auto i = range.first; i < range.second; ++i) {
                if (this->hasReservoirConnection_.at(i)) {
                    connected = this->initialisationPartner_[i];
                    found = true;
                    break;
                }
            }

            if (!found) {
                OPM_THROW(std::runtime_error,
                          fmt::format("Numerical aquifer {} has no connection to any "
                                      "reservoir cell owned by this process", id));
            }

            for (auto i = range.first; i < range.second; ++i) {
                if (!this->hasReservoirConnection_.at(i)) {
                    this->initialisationPartner_[i] = connected;
                }
            }
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

            // Order matters, as it does on the grid path: assignNaive() decides what the
            // primary variables mean and needs the PVT region to do it.
            solution[globalIdx].setPvtRegionIndex(this->pvtRegionIndex(localIdx));
            solution[globalIdx].assignNaive(fs);
        }
    }

    void linearize(SparseMatrixAdapter&, GlobalEqVector&) override
    {
        // Nothing to add: the aquifer cells are assembled by the model's own local
        // residual, like grid cells, and they are always active so there are no dormant
        // rows to condition.
    }

    //! The aquifer identifiers this module represents, in declaration order.
    std::vector<int> aquiferIds() const
    {
        std::vector<int> ids;
        ids.reserve(this->aquiferRange_.size());
        for (const auto& [id, range] : this->aquiferRange_) {
            static_cast<void>(range);
            ids.push_back(static_cast<int>(id));
        }

        return ids;
    }

    //! The global degrees of freedom carrying one aquifer, in declaration order.
    std::vector<unsigned> aquiferDofs(const int aquiferId) const
    {
        std::vector<unsigned> dofs;

        auto pos = this->aquiferRange_.find(static_cast<std::size_t>(aquiferId));
        if (pos == this->aquiferRange_.end()) {
            return dofs;
        }

        dofs.reserve(pos->second.second - pos->second.first);
        for (auto localIdx = pos->second.first; localIdx < pos->second.second; ++localIdx) {
            dofs.push_back(static_cast<unsigned>(this->localToGlobalDof(localIdx)));
        }

        return dofs;
    }

    /*!
     * \brief The connections of one aquifer that reach into the grid.
     *
     * The intra-aquifer chain is left out: what an aquifer reports as its influx is what
     * crosses into the reservoir, not what moves inside itself.  Each entry is
     * (aquifer degree of freedom, reservoir degree of freedom).
     */
    std::vector<std::pair<unsigned, unsigned>> reservoirConnections(const int aquiferId) const
    {
        std::vector<std::pair<unsigned, unsigned>> conns;

        const auto dofs = this->aquiferDofs(aquiferId);
        const auto isMine = [&dofs](const unsigned dof) {
            return std::find(dofs.begin(), dofs.end(), dof) != dofs.end();
        };
        const auto isGridDof = [this](const unsigned dof) {
            return dof < this->simulator_.model().numGridDof();
        };

        for (const auto& conn : this->connections_) {
            if (isMine(conn.dof1) && isGridDof(conn.dof2)) {
                conns.emplace_back(conn.dof1, conn.dof2);
            }
            else if (isMine(conn.dof2) && isGridDof(conn.dof1)) {
                conns.emplace_back(conn.dof2, conn.dof1);
            }
        }

        return conns;
    }

    //! Initial pressure of one aquifer's cells, for the restart record.
    std::vector<Scalar> initialPressure(const int aquiferId) const
    {
        std::vector<Scalar> pressures;

        auto pos = this->aquiferRange_.find(static_cast<std::size_t>(aquiferId));
        if (pos == this->aquiferRange_.end()) {
            return pressures;
        }

        pressures.reserve(pos->second.second - pos->second.first);
        for (auto localIdx = pos->second.first; localIdx < pos->second.second; ++localIdx) {
            const auto& cell = *this->cells_.at(localIdx);
            pressures.push_back(cell.init_pressure.has_value()
                                ? static_cast<Scalar>(cell.init_pressure.value())
                                : Scalar{0});
        }

        return pressures;
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
    std::vector<bool> hasReservoirConnection_{};

    //! aquifer id -> [first, last) range of its cells in the local numbering
    std::map<std::size_t, std::pair<std::size_t, std::size_t>> aquiferRange_{};
};

} // namespace Opm

#endif // OPM_NUMERICAL_AQUIFER_AUX_CELLS_HPP
