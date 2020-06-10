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
 * \copydoc Opm::BaseVanguard
 */
#ifndef EWOMS_BASE_VANGUARD_HH
#define EWOMS_BASE_VANGUARD_HH

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/parametersystem.hh>

#include <dune/common/version.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/space/common/dofmanager.hh>
#endif

#include <type_traits>
#include <memory>

namespace Opm {

/*!
 * \brief Provides the base class for most (all?) simulator vanguards.
 */
template <class TypeTag>
class BaseVanguard
{
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Implementation = GetPropType<TypeTag, Properties::Vanguard>;

#if HAVE_DUNE_FEM
    using GridPart = GetPropType<TypeTag, Properties::GridPart>;
#endif

public:
    BaseVanguard(Simulator& simulator)
        : simulator_(simulator)
    {}

    BaseVanguard(const BaseVanguard&) = delete;

    /*!
     * \brief Returns a reference to the grid view to be used.
     */
    const GridView& gridView() const
    { return *gridView_; }

#if HAVE_DUNE_FEM
    /*!
     * \brief Returns a reference to the grid part to be used.
     */
    const GridPart& gridPart() const
    { return *gridPart_; }

    /*!
     * \brief Returns a reference to the grid part to be used.
     */
    GridPart& gridPart()
    { return *gridPart_; }
#endif

    /*!
     * \brief Returns the number of times the grid has been changed since its creation.
     *
     * This basically says how often the grid has been adapted in the current simulation
     * run.
     */
    int gridSequenceNumber () const
    {
#if HAVE_DUNE_FEM
        using FemDofManager = Dune::Fem::DofManager< Grid >;
        return FemDofManager::instance( asImp_().grid() ).sequence();
#else
        return 0; // return the same sequence number >= 0 means the grid never changes
#endif
    }


    /*!
     * \brief Distribute the grid (and attached data) over all
     *        processes.
     */
    void loadBalance()
    {
        asImp_().grid().loadBalance();
        updateGridView_();
    }

protected:
    // this method should be called after the grid has been allocated
    void finalizeInit_()
    {
        updateGridView_();
    }

    void updateGridView_()
    {
#if HAVE_DUNE_FEM
        // first delete old grid part
        // this is due to a bug in dune-fem (dangling reference)
        gridPart_.reset();
        gridPart_.reset(new GridPart(asImp_().grid()));
        gridView_.reset(new GridView(static_cast<GridView>(*gridPart_)));
        assert(gridView_->size(0) == asImp_().grid().leafGridView().size(0));
#else
        gridView_.reset(new GridView(asImp_().grid().leafGridView()));
#endif
    }

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    Simulator& simulator_;
#if HAVE_DUNE_FEM
    std::unique_ptr<GridPart> gridPart_;
#endif
    std::unique_ptr<GridView> gridView_;
};

} // namespace Opm

#endif
