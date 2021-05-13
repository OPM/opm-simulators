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
 * \copydoc Opm::EclCpGridVanguard
 */
#ifndef EWOMS_ECL_CP_GRID_GENERIC_VANGUARD_HH
#define EWOMS_ECL_CP_GRID_GENERIC_VANGUARD_HH

#include <ebos/eclgenericvanguard.hh>
#include <opm/grid/CpGrid.hpp>

#include <functional>

namespace Opm {

template<class ElementMapper, class GridView, class Scalar>
class EclGenericCpGridVanguard {
protected:
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Dune::CpGrid>;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    EclGenericCpGridVanguard();

    /*!
     * \brief Return a reference to the simulation grid.
     */
    Dune::CpGrid& grid()
    { return *grid_; }

    /*!
     * \brief Return a reference to the simulation grid.
     */
    const Dune::CpGrid& grid() const
    { return *grid_; }

    /*!
     * \brief Returns a refefence to the grid which should be used by the EQUIL
     *        initialization code.
     *
     * The EQUIL keyword is used to specify the initial condition of the reservoir in
     * hydrostatic equilibrium. Since the code which does this is not accepting arbitrary
     * DUNE grids (the code is part of the opm-core module), this is not necessarily the
     * same as the grid which is used for the actual simulation.
     */
    const Dune::CpGrid& equilGrid() const;

    /*!
     * \brief Indicates that the initial condition has been computed and the memory used
     *        by the EQUIL grid can be released.
     *
     * Depending on the implementation, subsequent accesses to the EQUIL grid lead to
     * crashes.
     */
    void releaseEquilGrid();

    /// \brief Sets a function that returns external load balancing information when passed the grid
    ///
    /// The information is a vector of integers indication the partition index for each cell id.
    static void setExternalLoadBalancer(const std::function<std::vector<int> (const Dune::CpGrid&)>& loadBalancer)
    {
        externalLoadBalancer_ = loadBalancer;
    }

    /*!
     * \brief Returns the object which maps a global element index of the simulation grid
     *        to the corresponding element index of the logically Cartesian index.
     */
    const CartesianIndexMapper& cartesianIndexMapper() const;

    /*!
     * \brief Returns mapper from compressed to cartesian indices for the EQUIL grid
     */
    const CartesianIndexMapper& equilCartesianIndexMapper() const;

protected:
    /*!
     * \brief Distribute the simulation grid over multiple processes
     *
     * (For parallel simulation runs.)
     */
#if HAVE_MPI
    void doLoadBalance_(Dune::EdgeWeightMethod edgeWeightsMethod,
                        bool ownersFirst, bool serialPartitioning,
                        bool enableDistributedWells, double zoltanImbalanceTol,
                        const GridView& gridv, const Schedule& schedule,
                        std::vector<double>& centroids,
                        EclipseState& eclState,
                        EclGenericVanguard::ParallelWellStruct& parallelWells);

    void distributeFieldProps_(EclipseState& eclState);
#endif

    void allocCartMapper();

    void doCreateGrids_(EclipseState& eclState);

    virtual void allocTrans() = 0;
    virtual double getTransmissibility(unsigned I, unsigned J) = 0;

    // removing some connection located in inactive grid cells
    void doFilterConnections_(Schedule& schedule);

    Scalar computeCellThickness(const Element& element) const;

    std::unique_ptr<Dune::CpGrid> grid_;
    std::unique_ptr<Dune::CpGrid> equilGrid_;
    std::unique_ptr<CartesianIndexMapper> cartesianIndexMapper_;
    std::unique_ptr<CartesianIndexMapper> equilCartesianIndexMapper_;

    /// \brief optional functor returning external load balancing information
    ///
    /// If it is set then this will be used during loadbalance.
    static std::optional<std::function<std::vector<int> (const Dune::CpGrid&)>> externalLoadBalancer_;
    int mpiRank;
};

} // namespace Opm

#endif
