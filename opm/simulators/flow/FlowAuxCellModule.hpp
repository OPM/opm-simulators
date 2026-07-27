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
 * \copydoc Opm::FlowAuxCellModule
 */
#ifndef OPM_FLOW_AUX_CELL_MODULE_HPP
#define OPM_FLOW_AUX_CELL_MODULE_HPP

#include <opm/models/discretization/common/baseauxiliarymodule.hh>

#include <cstddef>
#include <vector>

namespace Opm {

/*!
 * \brief Base class for auxiliary modules whose degrees of freedom are *cells*.
 *
 * An auxiliary cell satisfies the same conservation equations as a grid cell, and is
 * assembled by the same local residual.  What it does not have is geometry: its pore
 * volume, depth and regions are authored rather than derived from a grid entity, and so
 * is its connection list.  A numerical aquifer is the simplest example -- it is pure
 * bookkeeping, a volume and a list of connections -- and fracture flow cells are the
 * dynamic one, where the same quantities are recomputed as the aperture changes.
 *
 * This is deliberately not the shape of the well models, which are also auxiliary
 * modules.  Their unknowns are of a different kind and they assemble and scale their own
 * equations in linearize(); they leave carriesModelEquations() false.  A module derived
 * from this class declares the opposite, and in exchange has to supply, per degree of
 * freedom, everything the model would otherwise read off the grid.
 *
 * Indices passed to the accessors below are *module-local*, in [0, numDofs());
 * localToGlobalDof() converts to the model's numbering.  The connections reported by
 * connections() use global numbering, because they may name grid cells.
 */
template <class TypeTag>
class FlowAuxCellModule : public BaseAuxiliaryModule<TypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ParentType = BaseAuxiliaryModule<TypeTag>;

protected:
    using NeighborSet = typename ParentType::NeighborSet;

public:
    /*!
     * \brief A flux connection authored by this module.
     *
     * Either endpoint may be a grid cell or an auxiliary cell; both are global degree of
     * freedom indices.  The transmissibility is the whole of the geometry -- there is no
     * face area to multiply by.
     */
    struct Connection
    {
        unsigned dof1{};
        unsigned dof2{};
        Scalar trans{};

        //! Thermal half transmissibilities, dof1->dof2 and dof2->dof1.  Only consulted
        //! when the model solves an energy equation.
        Scalar thermalHalfTrans12{};
        Scalar thermalHalfTrans21{};
    };

    //! The degrees of freedom of an auxiliary cell module are cells.
    bool carriesModelEquations() const override
    { return true; }

    /*!
     * \brief The connections of this module, each reported exactly once.
     *
     * The discretization enters every connection in the neighbour list of *both* of its
     * endpoints and assembles it from both sides, so reporting one per pair is correct;
     * reporting it twice would double the flux.
     */
    virtual void connections(std::vector<Connection>& conns) const = 0;

    /*!
     * \brief Pore volume of an auxiliary cell.
     */
    virtual Scalar poreVolume(unsigned localIdx) const = 0;

    /*!
     * \brief Bulk volume of an auxiliary cell.
     *
     * The model divides the pore volume by this to obtain a porosity, so the two have to
     * be authored consistently.  A module with no meaningful bulk volume should report
     * the pore volume and accept a porosity of one.
     */
    virtual Scalar bulkVolume(unsigned localIdx) const = 0;

    /*!
     * \brief Datum depth of an auxiliary cell, used for the gravity head.
     */
    virtual Scalar depth(unsigned localIdx) const = 0;

    //! Zero-based PVT region of an auxiliary cell.
    virtual unsigned pvtRegionIndex(unsigned localIdx) const = 0;

    //! Zero-based saturation-function region of an auxiliary cell.
    virtual unsigned satRegionIndex(unsigned localIdx) const = 0;

    /*!
     * \brief A grid cell whose initial state this auxiliary cell is derived from.
     *
     * Auxiliary cells cannot be equilibrated by the ordinary machinery, which needs the
     * cell's geometry.  Reporting a partner lets the problem start from that cell's
     * initial fluid state; a module which sets its own initial state in applyInitial()
     * may ignore this.
     */
    virtual unsigned initialisationPartner(unsigned localIdx) const = 0;

    /*!
     * \brief Whether this auxiliary cell currently takes part in the flow problem.
     *
     * A module may preallocate degrees of freedom it does not use yet -- fracture cells
     * that have not opened.  Such a cell has no volume and no connections, so its row
     * would be empty; the module is responsible for conditioning it in linearize().
     */
    virtual bool isActive(unsigned /*localIdx*/) const
    { return true; }

    /*!
     * \brief Report the sparsity contribution of this module's connections.
     *
     * Both directions, plus the diagonal, so that a cell with no connections still has a
     * block to be conditioned in.
     */
    void addNeighbors(std::vector<NeighborSet>& neighbors) const override
    {
        std::vector<Connection> conns;
        this->connections(conns);

        for (const auto& conn : conns) {
            neighbors[conn.dof1].insert(conn.dof2);
            neighbors[conn.dof2].insert(conn.dof1);
        }

        for (unsigned localIdx = 0; localIdx < this->numDofs(); ++localIdx) {
            const auto globalIdx = static_cast<unsigned>(this->localToGlobalDof(localIdx));
            neighbors[globalIdx].insert(globalIdx);
        }
    }

    /*!
     * \brief Hand the same connections to the discretization, which builds the
     *        neighbour info from them.
     */
    void addConnections(std::vector<typename ParentType::AuxiliaryConnection>& conns) const override
    {
        std::vector<Connection> own;
        this->connections(own);

        conns.reserve(conns.size() + own.size());
        for (const auto& conn : own) {
            conns.push_back({conn.dof1, conn.dof2});
        }
    }
};

} // namespace Opm

#endif // OPM_FLOW_AUX_CELL_MODULE_HPP
