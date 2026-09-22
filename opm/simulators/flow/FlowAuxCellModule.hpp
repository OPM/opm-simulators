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
 * \brief Base class for auxiliary modules whose degrees of freedom are cells without
 *        geometry, assembled by the model's local residual.
 *
 * Accessor indices are module-local; Connection endpoints are global DOFs.
 */
template <class TypeTag>
class FlowAuxCellModule : public BaseAuxiliaryModule<TypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ParentType = BaseAuxiliaryModule<TypeTag>;

protected:
    using NeighborSet = typename ParentType::NeighborSet;

public:
    //! Endpoints are global DOFs, grid or auxiliary; there is no face area.
    struct Connection
    {
        unsigned dof1{};
        unsigned dof2{};
        Scalar trans{};

        //! dof1->dof2 and dof2->dof1; energy equation only.
        Scalar thermalHalfTrans12{};
        Scalar thermalHalfTrans21{};
    };

    bool carriesModelEquations() const override
    { return true; }

    Scalar dofVolume(unsigned localIdx) const override
    { return this->bulkVolume(localIdx); }

    //! Report each connection once; it is assembled from both sides, so twice doubles the flux.
    virtual void connections(std::vector<Connection>& conns) const = 0;

    virtual Scalar poreVolume(unsigned localIdx) const = 0;

    //! Porosity is poreVolume/bulkVolume; without a bulk volume, report the pore volume.
    virtual Scalar bulkVolume(unsigned localIdx) const = 0;

    virtual Scalar depth(unsigned localIdx) const = 0;

    //! Zero-based.
    virtual unsigned pvtRegionIndex(unsigned localIdx) const = 0;

    //! Zero-based.
    virtual unsigned satRegionIndex(unsigned localIdx) const = 0;

    //! Input-grid cell supplying the reporting regions; -1 falls back to the
    //! initialisation partner's.
    virtual int hostCartesianIndex(unsigned /*localIdx*/) const
    { return -1; }

    //! Grid cell to take the initial state from; equilibration needs geometry.
    virtual unsigned initialisationPartner(unsigned localIdx) const = 0;

    //! Inactive (preallocated) cells have empty rows; the module conditions them in linearize().
    virtual bool isActive(unsigned /*localIdx*/) const
    { return true; }

    //! Diagonal too, so an unconnected cell still has a block to condition.
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
