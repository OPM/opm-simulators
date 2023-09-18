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
 * \copydoc Opm::EclBaseVanguard
 */
#ifndef EWOMS_ECL_BASE_VANGUARD_HH
#define EWOMS_ECL_BASE_VANGUARD_HH

#include <ebos/eclgenericvanguard.hh>

#include <opm/grid/common/GridEnums.hpp>
#include <opm/grid/common/CartesianIndexMapper.hpp>
#include <opm/grid/LookUpCellCentroid.hh>
#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquiferCell.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/models/discretization/common/fvbaseproperties.hh>
#include <opm/models/io/basevanguard.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/BlackoilModelParametersEbos.hpp>



#include <array>
#include <cstddef>
#include <optional>
#include <unordered_map>
#include <vector>

namespace Opm {
template <class TypeTag>
class EclBaseVanguard;
template<typename Grid, typename GridView> struct LookUpCellCentroid;
}

namespace Opm::Properties {

namespace TTag {
struct EclBaseVanguard {};
}

// declare the properties required by the for the ecl simulator vanguard
template<class TypeTag, class MyTypeTag>
struct EquilGrid {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EnableOpmRstFile {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct ParsingStrictness {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct SchedRestart {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EclOutputInterval {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct IgnoreKeywords {
    using type = UndefinedProperty;
};
template<class TypeTag, class MyTypeTag>
struct EdgeWeightsMethod {
    using type = UndefinedProperty;
};

#if HAVE_OPENCL || HAVE_ROCSPARSE
template<class TypeTag, class MyTypeTag>
struct NumJacobiBlocks {
    using type = UndefinedProperty;
};
#endif  // HAVE_OPENCL || HAVE_ROCSPARSE

template<class TypeTag, class MyTypeTag>
struct OwnerCellsFirst {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct SerialPartitioning {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct ZoltanImbalanceTol {
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct ZoltanParams {
    using type = UndefinedProperty;
};

template <class TypeTag, class MyTypeTag>
struct ExternalPartition
{
    using type = UndefinedProperty;
};

template<class TypeTag, class MyTypeTag>
struct AllowDistributedWells {
    using type = UndefinedProperty;
};

template<class TypeTag>
struct IgnoreKeywords<TypeTag, TTag::EclBaseVanguard> {
    static constexpr auto value = "";
};
template<class TypeTag>
struct EclDeckFileName<TypeTag, TTag::EclBaseVanguard> {
    static constexpr auto value = "";
};
template<class TypeTag>
struct EclOutputInterval<TypeTag, TTag::EclBaseVanguard> {
    static constexpr int value = -1;
};
template<class TypeTag>
struct EnableOpmRstFile<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct ParsingStrictness<TypeTag, TTag::EclBaseVanguard> {
    static constexpr auto value = "normal";
};
template<class TypeTag>
struct SchedRestart<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EdgeWeightsMethod<TypeTag, TTag::EclBaseVanguard> {
    static constexpr int value = 1;
};

#if HAVE_OPENCL || HAVE_ROCSPARSE
template<class TypeTag>
struct NumJacobiBlocks<TypeTag, TTag::EclBaseVanguard> {
    static constexpr int value = 0;
};
#endif // HAVE_OPENCL || HAVE_ROCSPARSE

template<class TypeTag>
struct OwnerCellsFirst<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct SerialPartitioning<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct ZoltanImbalanceTol<TypeTag, TTag::EclBaseVanguard> {
    static constexpr double value = 1.1;
};

template<class TypeTag>
struct ZoltanParams<TypeTag,TTag::EclBaseVanguard> {
    static constexpr auto value = "graph";
};

template <class TypeTag>
struct ExternalPartition<TypeTag, TTag::EclBaseVanguard>
{
    static constexpr auto* value = "";
};

template<class TypeTag>
struct AllowDistributedWells<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = false;
};

template<class T1, class T2>
struct UseMultisegmentWell;

// Same as in BlackoilModelParametersEbos.hpp but for here.
template<class TypeTag>
struct UseMultisegmentWell<TypeTag, TTag::EclBaseVanguard> {
    static constexpr bool value = true;
};
} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 */
template <class TypeTag>
class EclBaseVanguard : public BaseVanguard<TypeTag>,
                        public EclGenericVanguard
{
    using ParentType = BaseVanguard<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Vanguard>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;

    enum { enableExperiments = getPropValue<TypeTag, Properties::EnableExperiments>() };

public:
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

protected:
    static const int dimension = Grid::dimension;
    static const int dimensionworld = Grid::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;


public:
    /*!
     * \brief Register the common run-time parameters for all ECL simulator vanguards.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, std::string, EclDeckFileName,
                             "The name of the file which contains the ECL deck to be simulated");
        EWOMS_REGISTER_PARAM(TypeTag, int, EclOutputInterval,
                             "The number of report steps that ought to be skipped between two writes of ECL results");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableOpmRstFile,
                             "Include OPM-specific keywords in the ECL restart file to enable restart of OPM simulators from these files");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, IgnoreKeywords,
                             "List of Eclipse keywords which should be ignored. As a ':' separated string.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, ParsingStrictness,
                             "Set strictness of parsing process. Available options are "
                             "normal (stop for critical errors), "
                             "high (stop for all errors) and "
                             "low (as normal, except do not stop due to unsupported keywords even if marked critical");
        EWOMS_REGISTER_PARAM(TypeTag, bool, SchedRestart,
                             "When restarting: should we try to initialize wells and groups from historical SCHEDULE section.");
        EWOMS_REGISTER_PARAM(TypeTag, int, EdgeWeightsMethod,
                             "Choose edge-weighing strategy: 0=uniform, 1=trans, 2=log(trans).");

#if HAVE_OPENCL || HAVE_ROCSPARSE
        EWOMS_REGISTER_PARAM(TypeTag, int, NumJacobiBlocks,
                             "Number of blocks to be created for the Block-Jacobi preconditioner.");
#endif

        EWOMS_REGISTER_PARAM(TypeTag, bool, OwnerCellsFirst,
                             "Order cells owned by rank before ghost/overlap cells.");
#if HAVE_MPI
        EWOMS_REGISTER_PARAM(TypeTag, bool, SerialPartitioning,
                             "Perform partitioning for parallel runs on a single process.");
        EWOMS_REGISTER_PARAM(TypeTag, double, ZoltanImbalanceTol,
                             "Tolerable imbalance of the loadbalancing provided by Zoltan (default: 1.1).");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, ZoltanParams,
                             "Configuration of Zoltan partitioner. "
                             "Valid options are: graph, hypergraph or scotch. "
                             "Alternatively, you can request a configuration to be read "
                             "from a JSON file by giving the filename here, ending with '.json.' "
                             "See https://sandialabs.github.io/Zoltan/ug_html/ug.html "
                             "for available Zoltan options.");
        EWOMS_HIDE_PARAM(TypeTag, ZoltanParams);
        EWOMS_REGISTER_PARAM(TypeTag, std::string, ExternalPartition,
                             "Name of file from which to load an externally generated "
                             "partitioning of the model's active cells for MPI "
                             "distribution purposes. If empty, the built-in partitioning "
                             "method will be employed.");
        EWOMS_HIDE_PARAM(TypeTag, ExternalPartition);
#endif
        EWOMS_REGISTER_PARAM(TypeTag, bool, AllowDistributedWells,
                             "Allow the perforations of a well to be distributed to interior of multiple processes");
        // register here for the use in the tests without BlackoildModelParametersEbos
        EWOMS_REGISTER_PARAM(TypeTag, bool, UseMultisegmentWell, "Use the well model for multi-segment wells instead of the one for single-segment wells");

    }

    /*!
     * \brief Create the grid for problem data files which use the ECL file format.
     *
     * This is the file format used by the commercial ECLiPSE simulator. Usually it uses
     * a cornerpoint description of the grid.
     */
    EclBaseVanguard(Simulator& simulator)
        : ParentType(simulator)
    {
        fileName_ = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);
        edgeWeightsMethod_   = Dune::EdgeWeightMethod(EWOMS_GET_PARAM(TypeTag, int, EdgeWeightsMethod));

#if HAVE_OPENCL || HAVE_ROCSPARSE
        numJacobiBlocks_ = EWOMS_GET_PARAM(TypeTag, int, NumJacobiBlocks);
#endif

        ownersFirst_ = EWOMS_GET_PARAM(TypeTag, bool, OwnerCellsFirst);
#if HAVE_MPI
        serialPartitioning_ = EWOMS_GET_PARAM(TypeTag, bool, SerialPartitioning);
        zoltanImbalanceTol_ = EWOMS_GET_PARAM(TypeTag, double, ZoltanImbalanceTol);
        zoltanParams_ = EWOMS_GET_PARAM(TypeTag, std::string, ZoltanParams);
        externalPartitionFile_ = EWOMS_GET_PARAM(TypeTag, std::string, ExternalPartition);
#endif
        enableDistributedWells_ = EWOMS_GET_PARAM(TypeTag, bool, AllowDistributedWells);
        ignoredKeywords_ = EWOMS_GET_PARAM(TypeTag, std::string, IgnoreKeywords);
        int output_param = EWOMS_GET_PARAM(TypeTag, int, EclOutputInterval);
        if (output_param >= 0)
            outputInterval_ = output_param;
        useMultisegmentWell_ = EWOMS_GET_PARAM(TypeTag, bool, UseMultisegmentWell);
        enableExperiments_ = enableExperiments;

        init();
    }


    const CartesianIndexMapper& cartesianMapper() const
    {  return asImp_().cartesianIndexMapper(); }

    /*!
     * \brief Returns the number of logically Cartesian cells in each direction
     */
    const std::array<int, dimension>& cartesianDimensions() const
    { return asImp_().cartesianIndexMapper().cartesianDimensions(); }

    /*!
     * \brief Returns the overall number of cells of the logically Cartesian grid
     */
    int cartesianSize() const
    { return asImp_().cartesianIndexMapper().cartesianSize(); }

    /*!
     * \brief Returns the overall number of cells of the logically EquilCartesian grid
     */
    int equilCartesianSize() const
    { return asImp_().equilCartesianIndexMapper().cartesianSize(); }

    /*!
     * \brief Returns the Cartesian cell id for identifaction with ECL data
     */
    unsigned cartesianIndex(unsigned compressedCellIdx) const
    { return asImp_().cartesianIndexMapper().cartesianIndex(compressedCellIdx); }

    /*!
     * \brief Return the index of the cells in the logical Cartesian grid
     */
    unsigned cartesianIndex(const std::array<int,dimension>& coords) const
    {
        unsigned cartIndex = coords[0];
        int factor = cartesianDimensions()[0];
        for (unsigned i = 1; i < dimension; ++i) {
            cartIndex += coords[i]*factor;
            factor *= cartesianDimensions()[i];
        }

        return cartIndex;
    }

    /*!
     * \brief Return compressed index from cartesian index
     *
     * \return compressed index of cell is in interior, -1 otherwise
     *
     */
    int compressedIndex(int cartesianCellIdx) const
    {
        auto index_pair = cartesianToCompressed_.find(cartesianCellIdx);
        if (index_pair!=cartesianToCompressed_.end())
            return index_pair->second;
        else
            return -1;
    }

    /*!
     * \brief Return compressed index from cartesian index only in interior
     *
     * \return compressed index of cell is in interior, -1 otherwise
     *
     */
    int compressedIndexForInterior(int cartesianCellIdx) const
    {
        auto index_pair = cartesianToCompressed_.find(cartesianCellIdx);
        if (index_pair == cartesianToCompressed_.end() ||
            !is_interior_[index_pair->second])
        {
            return -1;
        }
        else
        {
            return index_pair->second;
        }
    }
    /*!
     * \brief Extract Cartesian index triplet (i,j,k) of an active cell.
     *
     * \param [in] cellIdx Active cell index.
     * \param [out] ijk Cartesian index triplet
     */
    void cartesianCoordinate(unsigned cellIdx, std::array<int,3>& ijk) const
    { return asImp_().cartesianIndexMapper().cartesianCoordinate(cellIdx, ijk); }

    /*!
     * \brief Returns the Cartesian cell id given an element index for the grid used for equilibration
     */
    unsigned equilCartesianIndex(unsigned compressedEquilCellIdx) const
    { return asImp_().equilCartesianIndexMapper().cartesianIndex(compressedEquilCellIdx); }

    /*!
     * \brief Extract Cartesian index triplet (i,j,k) of an active cell of the grid used for EQUIL.
     *
     * \param [in] cellIdx Active cell index.
     * \param [out] ijk Cartesian index triplet
     */
    void equilCartesianCoordinate(unsigned cellIdx, std::array<int,3>& ijk) const
    { return asImp_().equilCartesianIndexMapper().cartesianCoordinate(cellIdx, ijk); }


    /*!
     * \brief Returns the depth of a degree of freedom [m]
     *
     * For ECL problems this is defined as the average of the depth of an element and is
     * thus slightly different from the depth of an element's centroid.
     */
    Scalar cellCenterDepth(unsigned globalSpaceIdx) const
    {
        return cellCenterDepth_[globalSpaceIdx];
    }

    const std::vector<Scalar>& cellCenterDepths() const
    {
        return cellCenterDepth_;
    }

    /*!
     * \brief Returns the thickness of a degree of freedom [m]
     *
     * For ECL problems this is defined as the average of the depths of the top surface
     * corners minus the average of the depths of the bottom surface corners
     * The cell thickness is computed only when needed.
     */
    Scalar cellThickness(unsigned globalSpaceIdx) const
    {
        assert(!cellThickness_.empty());
        return cellThickness_[globalSpaceIdx];
    }

    /*!
     * \brief Get the number of cells in the global leaf grid view.
     * \warn This is a collective operation that needs to be called
     * on all ranks.
     */
    std::size_t globalNumCells() const
    {
        const auto& grid = asImp_().grid();
        if (grid.comm().size() == 1)
        {
            return grid.leafGridView().size(0);
        }
        const auto& gridView = grid.leafGridView();
        constexpr int codim = 0;
        constexpr auto Part = Dune::Interior_Partition;
        auto local_cells = std::distance(gridView.template begin<codim, Part>(),
                                         gridView.template end<codim, Part>());
        return grid.comm().sum(local_cells);
    }

    void setupCartesianToCompressed_() {
          this->updateCartesianToCompressedMapping_();
    }

protected:
    /*!
     * \brief Get function to query cell centroids for a distributed grid.
     *
     * Currently this only non-empty for a loadbalanced CpGrid.
     * It is a function return the centroid for the given element
     * index.
     * \param cartMapper The cartesian index mapper for lookup of
     *        cartesian indices
     */
    template<class CartMapper>
    std::function<std::array<double,dimensionworld>(int)>
    cellCentroids_(const CartMapper& cartMapper, const bool& isCpGrid) const
    {
        return [this, cartMapper, isCpGrid](int elemIdx) {
            std::array<double,dimensionworld> centroid;
            const auto rank = this->gridView().comm().rank();
            const auto maxLevel = this->gridView().grid().maxLevel();
            bool useEclipse = !isCpGrid || (isCpGrid && (rank == 0) && (maxLevel == 0));
            if (useEclipse)
            {
                centroid =  this->eclState().getInputGrid().getCellCenter(cartMapper.cartesianIndex(elemIdx));
            }
            else
            {
                LookUpCellCentroid<Grid,GridView> lookUpCellCentroid(this->gridView(), cartMapper, nullptr);
                centroid = lookUpCellCentroid(elemIdx);
            }
            return centroid;
        };
    }

    void callImplementationInit()
    {
        asImp_().createGrids_();
        asImp_().filterConnections_();
        std::string outputDir = EWOMS_GET_PARAM(TypeTag, std::string, OutputDir);
        bool enableEclCompatFile = !EWOMS_GET_PARAM(TypeTag, bool, EnableOpmRstFile);
        asImp_().updateOutputDir_(outputDir, enableEclCompatFile);
        asImp_().finalizeInit_();
    }

    void updateCartesianToCompressedMapping_()
    {
        std::size_t num_cells = asImp_().grid().leafGridView().size(0);
        is_interior_.resize(num_cells);

        ElementMapper elemMapper(this->gridView(), Dune::mcmgElementLayout());
        for (const auto& element : elements(this->gridView()))
        {
            const auto elemIdx = elemMapper.index(element);
            unsigned cartesianCellIdx = cartesianIndex(elemIdx);
            cartesianToCompressed_[cartesianCellIdx] = elemIdx;
            if (element.partitionType() == Dune::InteriorEntity)
            {
                is_interior_[elemIdx] = 1;
            }
            else
            {
                is_interior_[elemIdx] = 0;
            }
        }
    }

    void updateCellDepths_()
    {
        int numCells = this->gridView().size(/*codim=*/0);
        cellCenterDepth_.resize(numCells);

        ElementMapper elemMapper(this->gridView(), Dune::mcmgElementLayout());

        const auto num_aqu_cells = this->allAquiferCells();

        for(const auto& element : elements(this->gridView())) {
            const unsigned int elemIdx = elemMapper.index(element);
            cellCenterDepth_[elemIdx] = cellCenterDepth(element);

            if (!num_aqu_cells.empty()) {
               const unsigned int global_index = cartesianIndex(elemIdx);
               const auto search = num_aqu_cells.find(global_index);
               if (search != num_aqu_cells.end()) {
                    // updating the cell depth using aquifer cell depth
                    cellCenterDepth_[elemIdx] = search->second->depth;
                }
            }
        }
    }
    void updateCellThickness_()
    {
        if (!this->drsdtconEnabled())
            return;

        ElementMapper elemMapper(this->gridView(), Dune::mcmgElementLayout());

        int numElements = this->gridView().size(/*codim=*/0);
        cellThickness_.resize(numElements);

        for (const auto& elem : elements(this->gridView())) {
            const unsigned int elemIdx = elemMapper.index(elem);
            cellThickness_[elemIdx] = computeCellThickness(elem);
        }
    }

private:
    // computed from averaging cell corner depths
    Scalar cellCenterDepth(const Element& element) const
    {
        typedef typename Element::Geometry Geometry;
        static constexpr int zCoord = Element::dimension - 1;
        Scalar zz = 0.0;

        const Geometry& geometry = element.geometry();
        const int corners = geometry.corners();
        for (int i=0; i < corners; ++i)
            zz += geometry.corner(i)[zCoord];

        return zz/Scalar(corners);
    }

    Scalar computeCellThickness(const typename GridView::template Codim<0>::Entity& element) const
    {
        typedef typename Element::Geometry Geometry;
        static constexpr int zCoord = Element::dimension - 1;
        Scalar zz1 = 0.0;
        Scalar zz2 = 0.0;

        const Geometry& geometry = element.geometry();
        // This code only works with CP-grid where the
        // number of corners are 8 and
        // also assumes that the first
        // 4 corners are the top surface and
        // the 4 next are the bottomn.
        assert(geometry.corners() == 8);
        for (int i=0; i < 4; ++i){
            zz1 += geometry.corner(i)[zCoord];
            zz2 += geometry.corner(i+4)[zCoord];
        }
        zz1 /=4;
        zz2 /=4;
        return zz2-zz1;
     }

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

protected:
    /*! \brief Mapping between cartesian and compressed cells.
     *  It is initialized the first time it is called
     */
    std::unordered_map<int,int> cartesianToCompressed_;

    /*! \brief Cell center depths
     */
    std::vector<Scalar> cellCenterDepth_;

    /*! \brief Cell thickness
     */
    std::vector<Scalar> cellThickness_;

    /*! \brief Whether a cells is in the interior.
     */
    std::vector<int> is_interior_;
};

} // namespace Opm

#endif
