/*
  Copyright 2023 Inria, Bretagne–Atlantique Research Center

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
#include <opm/simulators/utils/GridDataOutput.hpp>

#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/io/file/vtk/common.hh>

#include <opm/common/ErrorMacros.hpp>

#include <cassert>
#include <cstddef>
#include <iterator>
#include <ostream>

namespace Opm::GridDataOutput {

template <class GridView, unsigned int partitions>
SimMeshDataAccessor<GridView,partitions>::
 SimMeshDataAccessor(const GridView& gridView, Dune::PartitionSet<partitions> dunePartition)
    : gridView_(gridView)
    , dunePartition_(dunePartition)
{
    dimw_ = GridView::dimension; // this is an enum
    partition_value_ = dunePartition.value;
    countEntities();
}

template <class GridView, unsigned int partitions>
bool SimMeshDataAccessor<GridView,partitions>::polyhedralCellPresent() const
{
    for (const auto& cit : elements(gridView_, dunePartition_)) {
        auto corner_geom = cit.geometry();
        if (Dune::VTK::geometryType(corner_geom.type()) == Dune::VTK::polyhedron) {
            return true;
        }
    }
    return false;
}

template <class GridView, unsigned int partitions>
void SimMeshDataAccessor<GridView,partitions>::countEntities()
{
    // We include all the vertices for this ranks partition
    const auto& vert_partition = vertices(gridView_, Dune::Partitions::all);
    nvertices_ = std::distance(vert_partition.begin(), vert_partition.end());

    const auto& cell_partition = elements(gridView_, dunePartition_);
    ncells_ = 0;
    ncorners_ = 0;
    for (const auto& cit : cell_partition) {
        auto corner_geom = cit.geometry();
        ncorners_ += corner_geom.corners();
        ++ncells_;
    }
}

template <class GridView, unsigned int partitions>
template <typename T>
long SimMeshDataAccessor<GridView,partitions>::
writeGridPoints(T* x_inout, T* y_inout, T* z_inout, long max_size) const
{
    if (max_size < nvertices_) {
        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeGridPoints( T&  x_inout,  T&  "
                  "y_inout, T& z_inout )  "
                      + " Input objects max_size (" + std::to_string(max_size)
                      + ") is not sufficient to fit the nvertices_ values (" + std::to_string(nvertices_) + ")");
    }

    long i = 0;
    if (dimw_ == 3) {
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0); // vertices only have one corner
            x_inout[i] = static_cast<T>(xyz_local[0]);
            y_inout[i] = static_cast<T>(xyz_local[1]);
            z_inout[i] = static_cast<T>(xyz_local[2]);
            i++;
        }
    } else if (dimw_ == 2) {
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0); // vertices only have one corner
            x_inout[i] = static_cast<T>(xyz_local[0]);
            y_inout[i] = static_cast<T>(xyz_local[1]);
            z_inout[i] = static_cast<T>(0.0);
            i++;
        }
    }
    assert(i == nvertices_); // As we are templated on the
                             // Dune::PartitionSet<partitions>, this cannot change
    return i;
}

template <class GridView, unsigned int partitions>
template <typename VectType>
long SimMeshDataAccessor<GridView,partitions>::
writeGridPoints(VectType& x_inout, VectType& y_inout, VectType& z_inout) const
{
    const std::size_t check_size_x = x_inout.size();
    const std::size_t check_size_y = y_inout.size();
    const std::size_t check_size_z = z_inout.size();

    using VT = decltype(x_inout.data()[0]);

    if ((check_size_x < static_cast<std::size_t>(nvertices_)) ||
        (check_size_y < static_cast<std::size_t>(nvertices_)) ||
        (check_size_z < static_cast<std::size_t>(nvertices_))) {
        // assert(check_size >= nvertices_);
        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeGridPoints( VectType&  x_inout,  VectType&  "
                  "y_inout, VectType& z_inout )  At least one of the inputs"
                      + "  object x size " + std::to_string(check_size_x) + "  object y size "
                      + std::to_string(check_size_y) + "  object z size " + std::to_string(check_size_z)
                      + " is not sufficient to fit the nvertices_ values( " + std::to_string(nvertices_) + " )");
    }

    long i = 0;
    if (dimw_ == 3) {
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0); // vertices only have one corner
            x_inout.data()[i] = static_cast<VT>(xyz_local[0]);
            y_inout.data()[i] = static_cast<VT>(xyz_local[1]);
            z_inout.data()[i] = static_cast<VT>(xyz_local[2]);
            i++;
        }
    } else if (dimw_ == 2) {
        double td = 0.0;
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0); // vertices only have one corner
            x_inout.data()[i] = static_cast<VT>(xyz_local[0]);
            y_inout.data()[i] = static_cast<VT>(xyz_local[1]);
            z_inout.data()[i] = static_cast<VT>(td);
            i++;
        }
    }
    assert(i == nvertices_); // As we are templated on the
                             // Dune::PartitionSet<partitions>, this cannot change
    return i;
}

template <class GridView, unsigned int partitions>
template <typename T>
long SimMeshDataAccessor<GridView,partitions>::
writeGridPoints_AOS(T* xyz_inout, long max_size) const
{
    if (max_size < nvertices_ * 3) {
        assert(max_size >= nvertices_ * 3);
        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeGridPoints_AOS( T*  xyz_inout )  " + " Input objects max_size ("
                      + std::to_string(max_size) + ") is not sufficient to fit the nvertices_ * 3 values ("
                      + std::to_string(nvertices_ * 3) + ")");
    }

    long i = 0;
    if (dimw_ == 3) {
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0);
            xyz_inout[i++] = static_cast<T>(xyz_local[0]);
            xyz_inout[i++] = static_cast<T>(xyz_local[1]);
            xyz_inout[i++] = static_cast<T>(xyz_local[2]);
        }
    } else if (dimw_ == 2) {
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0);
            xyz_inout[i++] = static_cast<T>(xyz_local[0]);
            xyz_inout[i++] = static_cast<T>(xyz_local[1]);
            xyz_inout[i++] = static_cast<T>(0.0);
        }
    }
    return i / 3;
}

template <class GridView, unsigned int partitions>
template <typename VectType>
long SimMeshDataAccessor<GridView,partitions>::
writeGridPoints_AOS(VectType& xyz_inout) const
{
    const std::size_t check_size = xyz_inout.size();

    using VT = decltype(xyz_inout.data()[0]);

    if (check_size < static_cast<std::size_t>(nvertices_ * 3)) {
        assert(check_size >= nvertices_ * 3);
        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeGridPoints_AOS( VectType&  xyz_inout )  "
                      + " Input objects check_size (" + std::to_string(check_size)
                      + ") is not sufficient to fit the nvertices_ * 3 values (" + std::to_string(nvertices_ * 3)
                      + ")");
    }

    long i = 0;
    if (dimw_ == 3) {
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0);
            xyz_inout.data()[i++] = static_cast<VT>(xyz_local[0]);
            xyz_inout[i++] = static_cast<VT>(xyz_local[1]);
            xyz_inout[i++] = static_cast<VT>(xyz_local[2]);
        }
    } else if (dimw_ == 2) {
        double td = 0.0;
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0);
            xyz_inout[i++] = static_cast<VT>(xyz_local[0]);
            xyz_inout[i++] = static_cast<VT>(xyz_local[1]);
            xyz_inout[i++] = static_cast<VT>(td);
        }
    }
    return i / 3;
}

template <class GridView, unsigned int partitions>
template <typename T>
long SimMeshDataAccessor<GridView,partitions>::
writeGridPoints_SOA(T* xyz_inout, long max_size) const
{
    if (max_size < nvertices_ * 3) {
        // assert(max_size >= nvertices_ * 3);
        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeGridPoints_SOA( T&  xyz_inout )  " + " Input objects max_size ("
                      + std::to_string(max_size) + ") is not sufficient to fit the nvertices_ * 3 values ("
                      + std::to_string(nvertices_ * 3) + ")");
    }

    long i = 0;
    // Get offsets into structure
    T* xyz_inout_y = xyz_inout + nvertices_;
    T* xyz_inout_z = xyz_inout + (2 * nvertices_);

    if (dimw_ == 3) {
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0);
            xyz_inout[i] = static_cast<T>(xyz_local[0]);
            xyz_inout_y[i] = static_cast<T>(xyz_local[1]);
            xyz_inout_z[i] = static_cast<T>(xyz_local[2]);
            i++;
        }
    } else if (dimw_ == 2) {
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0);
            xyz_inout[i] = static_cast<T>(xyz_local[0]);
            xyz_inout_y[i] = static_cast<T>(xyz_local[1]);
            xyz_inout_z[i] = static_cast<T>(0.0);
            i++;
        }
    }
    return i;
}

template <class GridView, unsigned int partitions>
template <typename VectType>
long SimMeshDataAccessor<GridView,partitions>::
writeGridPoints_SOA(VectType& xyz_inout) const
{
    const std::size_t check_size = xyz_inout.size();

    if (check_size < static_cast<std::size_t>(nvertices_ * 3)) {
        // assert(check_size >= nvertices_ * 3);
        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeGridPoints_SOA( VectType&  xyz_inout )  "
                      + " Input objects check_size (" + std::to_string(check_size)
                      + ") is not sufficient to fit the nvertices_ * 3 values (" + std::to_string(nvertices_ * 3)
                      + ")");
    }

    using VT = decltype(xyz_inout.data()[0]);

    long i = 0;
    // Get offsets into structure
    VT* xyz_inout_y = xyz_inout.data() + nvertices_;
    VT* xyz_inout_z = xyz_inout.data() + (2 * nvertices_);

    if (dimw_ == 3) {
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0);
            xyz_inout.data()[i] = static_cast<VT>(xyz_local[0]);
            xyz_inout_y[i] = static_cast<VT>(xyz_local[1]);
            xyz_inout_z[i] = static_cast<VT>(xyz_local[2]);
            i++;
        }
    } else if (dimw_ == 2) {
        double td = 0.0;
        for (const auto& vit : vertices(gridView_, Dune::Partitions::all)) {
            auto xyz_local = vit.geometry().corner(0);
            xyz_inout.data()[i] = static_cast<VT>(xyz_local[0]);
            xyz_inout_y[i] = static_cast<VT>(xyz_local[1]);
            xyz_inout_z[i] = static_cast<VT>(td);
            i++;
        }
    }
    return i;
}

template <class GridView, unsigned int partitions>
template <typename Integer>
long SimMeshDataAccessor<GridView,partitions>::
writeConnectivity(Integer* connectivity_inout,
                  ConnectivityVertexOrder whichOrder, long max_size) const
{
    if (max_size < ncorners_) {

        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeConnectivity( T*  connectivity_inout,... )  "
                      + " Input max_size value (" + std::to_string(max_size)
                      + ") is not sufficient to fit the ncorners_ values (" + std::to_string(ncorners_) + ")");
    }

    long i = 0;
    if (whichOrder == DUNE) {
        // DUNE order
        for (const auto& cit : elements(gridView_, dunePartition_)) {
            auto cell_corners = cit.geometry().corners();
            for (auto vx = 0; vx < cell_corners; ++vx) {
                const int vxIdx = gridView_.indexSet().subIndex(cit, vx, 3);
                connectivity_inout[i + vx] = vxIdx;
            }
            i += cell_corners;
        }
    } else {
        // VTK order
        for (const auto& cit : elements(gridView_, dunePartition_)) {
            auto cell_corners = cit.geometry().corners();
            for (auto vx = 0; vx < cell_corners; ++vx) {
                const int vxIdx = gridView_.indexSet().subIndex(cit, vx, 3);
                int vtkOrder;
                vtkOrder = Dune::VTK::renumber(cit.type(), vx);
                connectivity_inout[i + vtkOrder] = vxIdx;
            }
            i += cell_corners;
        }
    }
    return i;
}

template <class GridView, unsigned int partitions>
template <typename VectType>
long SimMeshDataAccessor<GridView,partitions>::
writeConnectivity(VectType& connectivity_inout,
                  ConnectivityVertexOrder whichOrder) const
{
    const std::size_t check_size = connectivity_inout.size();

    if (check_size < static_cast<std::size_t>(ncorners_)) {
        // assert(check_size >= ncorners_);
        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeConnectivity( VectType&  "
                  "connectivity_inout )  "
                      + " Input objects size (" + std::to_string(check_size)
                      + ") is not sufficient to fit the ncorners_ values (" + std::to_string(ncorners_) + ")");
    }

    long i = 0;
    if (whichOrder == DUNE) {
        // DUNE order
        for (const auto& cit : elements(gridView_, dunePartition_)) {
            auto cell_corners = cit.geometry().corners();
            for (auto vx = 0; vx < cell_corners; ++vx) {
                const int vxIdx = gridView_.indexSet().subIndex(cit, vx, 3);
                connectivity_inout.data()[i + vx] = vxIdx;
            }
            i += cell_corners;
        }
    } else {
        // VTK order
        for (const auto& cit : elements(gridView_, dunePartition_)) {
            auto cell_corners = cit.geometry().corners();
            for (auto vx = 0; vx < cell_corners; ++vx) {
                const int vxIdx = gridView_.indexSet().subIndex(cit, vx, 3);
                int vtkOrder;
                vtkOrder = Dune::VTK::renumber(cit.type(), vx);
                connectivity_inout.data()[i + vtkOrder] = vxIdx;
            }
            i += cell_corners;
        }
    }
    return i;
}

template <class GridView, unsigned int partitions>
template <typename Integer>
long SimMeshDataAccessor<GridView,partitions>::
writeOffsetsCells(Integer* offsets_inout, long max_size) const
{
    if (max_size < ncells_) {
        // assert(max_size >= ncells_);
        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeOffsetsCells( T*  offsets_inout )  " + " Input objects max_size ("
                      + std::to_string(max_size) + ") is not sufficient to fit the ncells_ values ("
                      + std::to_string(ncells_) + ")");
    }
    long i = 1;
    offsets_inout[0] = 0;
    for (const auto& cit : elements(gridView_, dunePartition_)) {
        auto cell_corners = cit.geometry().corners();
        offsets_inout[i] = offsets_inout[i - 1] + cell_corners;
        i++;
    }
    return i; // This should be 1 greater than ncells_
}

template <class GridView, unsigned int partitions>
template <typename VectType>
long SimMeshDataAccessor<GridView,partitions>::
writeOffsetsCells(VectType& offsets_inout) const
{
    const std::size_t check_size = offsets_inout.size();
    if (check_size < static_cast<std::size_t>(ncells_)) {
        // assert(check_size >= ncells_);
        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeOffsetsCells( VectType& "
                  "offsets_inout )  "
                      + " Input objects size (" + std::to_string(offsets_inout.size())
                      + ") is not sufficient to fit the ncells_ values (" + std::to_string(ncells_) + ")");
    }

    // using VT = decltype(offsets_inout.data()[0]);

    long i = 1;
    offsets_inout.data()[0] = 0;
    for (const auto& cit : elements(gridView_, dunePartition_)) {
        auto cell_corners = cit.geometry().corners();
        offsets_inout.data()[i] = offsets_inout.data()[i - 1] + cell_corners;
        i++;
    }
    return i; // This should be 1 greater than ncells_
}

template <class GridView, unsigned int partitions>
template <typename Integer>
long SimMeshDataAccessor<GridView,partitions>::
writeCellTypes(Integer* types_inout, long max_size) const
{
    if (max_size < ncells_) {
        // assert(max_size >= ncells_);
        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeCellTypes( T*  types_inout )  " + " Input objects max_size ("
                      + std::to_string(max_size) + ") is not sufficient to fit the ncells_ values ("
                      + std::to_string(ncells_) + ")");
    }
    int i = 0;
    for (const auto& cit : elements(gridView_, dunePartition_)) {
        Integer vtktype = static_cast<Integer>(Dune::VTK::geometryType(cit.type()));
        types_inout[i++] = vtktype;
    }
    return i;
}

template <class GridView, unsigned int partitions>
template <typename VectType>
long SimMeshDataAccessor<GridView,partitions>::
writeCellTypes(VectType& types_inout) const
{
    const std::size_t check_size = types_inout.size();

    if (check_size < static_cast<std::size_t>(ncells_)) {
        OPM_THROW(std::runtime_error,
                  "Opm::GridDataOutput::writeCellTypes( VectType&  types_inout )  " + " Input objects check_size ("
                      + std::to_string(check_size) + ") is not sufficient to fit the ncells_ values ("
                      + std::to_string(ncells_) + ")");
    }

    int i = 0;
    for (const auto& cit : elements(gridView_, dunePartition_)) {
        int vtktype = static_cast<int>(Dune::VTK::geometryType(cit.type()));
        types_inout.data()[i++] = vtktype;
    }
    return i;
}

template <class GridView, unsigned int partitions>
std::string SimMeshDataAccessor<GridView,partitions>::
getPartitionTypeString() const
{
    if (this->dunePartition_ == Dune::Partitions::all) {
        return "Dune::Partitions::all";
    }
    if (this->dunePartition_ == Dune::Partitions::interior) {
        return "Dune::Partitions::interior";
    }
    if (this->dunePartition_ == Dune::Partitions::interiorBorder) {
        return "Dune::Partitions::interiorBorder";
    }
    if (this->dunePartition_ == Dune::Partitions::interiorBorderOverlap) {
        return "Dune::Partitions::interiorBorderOverlap";
    }
    if (this->dunePartition_ == Dune::Partitions::front) {
        return "Dune::Partitions::front";
    }
    if (this->dunePartition_ == Dune::Partitions::interiorBorderOverlapFront) {
        return "Dune::Partitions::InteriorBorderOverlapFront";
    }
    if (this->dunePartition_ == Dune::Partitions::border) {
        return "Dune::Partitions::border";
    }
    if (this->dunePartition_ == Dune::Partitions::ghost) {
        return "Dune::Partitions::ghost";
    }

    return "Unknown Dune::PartitionSet<>";
}

template <class GridView, unsigned int partitions>
void SimMeshDataAccessor<GridView,partitions>::
printGridDetails(std::ostream& outstr) const
{
    outstr << "Dune Partition = " << partition_value_ << ", " << getPartitionTypeString() << std::endl;
    outstr << "ncells_: " << getNCells() << std::endl;
    outstr << "nvertices_: " << getNVertices() << std::endl;
    outstr << "ncorners_: " << getNCorners() << std::endl;
}

} // namespace Opm::GridDataOutput
