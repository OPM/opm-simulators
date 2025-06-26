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
 * \copydoc Opm::CpGridVanguard
 */
#ifndef OPM_GENERIC_CPGRID_VANGUARD_HPP
#define OPM_GENERIC_CPGRID_VANGUARD_HPP

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

#include <opm/simulators/flow/FlowGenericVanguard.hpp>

#include <functional>
#include <memory>
#include <optional>
#include <vector>

#if HAVE_MPI
#include <filesystem>
#endif

namespace Opm {
    class EclipseState;
    class Schedule;
    class Well;
    class ParallelEclipseState;
    class LgrCollection;
}

namespace Opm {

#if HAVE_MPI
namespace details {
class MPIPartitionFromFile
{
public:
    explicit MPIPartitionFromFile(const std::filesystem::path& partitionFile)
        : partitionFile_(partitionFile)
    {}

    std::vector<int> operator()(const Dune::CpGrid& grid) const;

private:
    std::filesystem::path partitionFile_{};
};

} // namespace Opm::details
#endif // HAVE_MPI


/// \brief optional functor returning external load balancing information
///
/// If it is set then this will be used during loadbalance.
extern std::optional<std::function<std::vector<int> (const Dune::CpGrid&)>> externalLoadBalancer;

template<class ElementMapper, class GridView, class Scalar>
class GenericCpGridVanguard {
protected:
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Dune::CpGrid>;
    using LevelCartesianIndexMapper = Opm::LevelCartesianIndexMapper<Dune::CpGrid>;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    GenericCpGridVanguard();

    virtual ~GenericCpGridVanguard() = default;

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
        externalLoadBalancer = loadBalancer;
    }

    /*!
     * \brief Returns the object which maps a global element index of the simulation grid
     *        to the corresponding element index of the logically Cartesian index.
     */
    const CartesianIndexMapper& cartesianIndexMapper() const;

    /*!
     * \brief Returns the object which maps a global element index of the simulation grid
     *        to the corresponding element index of the level logically Cartesian index.
     */
    const LevelCartesianIndexMapper levelCartesianIndexMapper() const;

    /*!
     * \brief Returns mapper from compressed to cartesian indices for the EQUIL grid
     */
    const CartesianIndexMapper& equilCartesianIndexMapper() const;

    const std::vector<int>& cellPartition() const
    {
        return this->cell_part_;
    }

protected:
    /*!
     * \brief Distribute the simulation grid over multiple processes
     *
     * (For parallel simulation runs.)
     */
#if HAVE_MPI
    void doLoadBalance_(const Dune::EdgeWeightMethod             edgeWeightsMethod,
                        const bool                               ownersFirst,
                        const bool                               addCorners,
                        const int                                numOverlap,
                        const Dune::PartitionMethod              partitionMethod,
                        const bool                               serialPartitioning,
                        const bool                               enableDistributedWells,
                        const bool                               allowSplittingInactiveWells,
                        const double                             imbalanceTol,
                        const GridView&                          gridView,
                        const Schedule&                          schedule,
                        EclipseState&                            eclState,
                        FlowGenericVanguard::ParallelWellStruct& parallelWells,
                        const int                                numJacobiBlocks,
                        const bool                               enableEclOutput);

    void distributeFieldProps_(EclipseState& eclState);

private:
    std::vector<double> extractFaceTrans(const GridView& gridView) const;

    void distributeGrid(const Dune::EdgeWeightMethod                          edgeWeightsMethod,
                        const bool                                            ownersFirst,
                        const bool                                            addCorners,
                        const int                                             numOverlap,
                        const Dune::PartitionMethod                           partitionMethod,
                        const bool                                            serialPartitioning,
                        const bool                                            enableDistributedWells,
                        const double                                          imbalanceTol,
                        const bool                                            loadBalancerSet,
                        const std::vector<double>&                            faceTrans,
                        const std::vector<Well>&                              wells,
                        const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                        EclipseState&                                         eclState,
                        FlowGenericVanguard::ParallelWellStruct&              parallelWells);

    void distributeGrid(const Dune::EdgeWeightMethod                          edgeWeightsMethod,
                        const bool                                            ownersFirst,
                        const bool                                            addCorners,
                        const int                                             numOverlap,
                        const Dune::PartitionMethod                           partitionMethod,
                        const bool                                            serialPartitioning,
                        const bool                                            enableDistributedWells,
                        const double                                          imbalanceTol,
                        const bool                                            loadBalancerSet,
                        const std::vector<double>&                            faceTrans,
                        const std::vector<Well>&                              wells,
                        const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
                        ParallelEclipseState*                                 eclState,
                        FlowGenericVanguard::ParallelWellStruct&              parallelWells);

protected:
    virtual const std::string& zoltanParams() const = 0;
    virtual double zoltanPhgEdgeSizeThreshold() const = 0;
    virtual const std::string& metisParams() const = 0;

#endif  // HAVE_MPI

    void allocCartMapper();

    void doCreateGrids_(EclipseState& eclState);
    void addLgrsUpdateLeafView(const LgrCollection& lgrCollection,
                               const int lgrsSize,
                               Dune::CpGrid& grid);

    virtual void allocTrans() = 0;
    virtual double getTransmissibility(unsigned I, unsigned J) const = 0;

    // removing some connection located in inactive grid cells
    void doFilterConnections_(Schedule& schedule);

    Scalar computeCellThickness(const Element& element) const;

    std::unique_ptr<Dune::CpGrid> grid_;
    std::unique_ptr<Dune::CpGrid> equilGrid_;
    std::unique_ptr<CartesianIndexMapper> cartesianIndexMapper_;
    std::unique_ptr<CartesianIndexMapper> equilCartesianIndexMapper_;
    std::unique_ptr<LevelCartesianIndexMapper> levelCartesianIndexMapper_;

    int mpiRank;
    std::vector<int> cell_part_{};
};

} // namespace Opm

#endif // OPM_GENERIC_CPGRID_VANGUARD_HPP
