// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
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

#ifndef OPM_SIM_MESH_DATA2_HPP
#define OPM_SIM_MESH_DATA2_HPP

#include <sstream>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/io/file/vtk/common.hh>

/** @file
    @brief Allows model geometry data to be passed to external code - via a copy direct to input pointers. 
    
    This data extractor provides the full set of vertices (corresponding to Dune::Partition::all) and then 
    allows a user to specify Dune sub-partitions to get the references into the vertex array and element 
    (aka cell) types for the sub-partition. This allows the full set of verticies to be reused for 
    visualisation of the various sub-partitions, at the expense of copying all the vertices. Typically
    a user is interested in the interiorBoarder elements which make use of the bulk (~80%) of the vertices.
    This saves having to renumber the indexes to the vertices for the sub-partitions.
    The vertex data can be retrieved as seperate x, y and z arrays, or as a single array of structures, 
    or as a single structure of arrays based 


    Example:
        // From the opm-simulators repository
        #include <opm/simulators/utils/GridDataOutput.hpp>
        
        // N.B. does not seem to be able to be allocated with new operator.
        Opm::GridDataOutput::SimMeshDataAccessor geomData(gridView,  Dune::Partition::interior ) ;

        geomData.printGridDetails() ;
        
        // example using seperate x, y and z arrays
        int nvert = geomData.getNVertices() ;
        double * x_vert = new double[nvert] ;
        double * y_vert = new double[nvert] ;
        double * z_vert = new double[nvert] ;
        geomData.writeGridPoints(x_vert,y_vert,z_vert) ;

        ... do something with vertex data x_vert, y_vert and z_vert ....

        free [] x_vert;
        free [] y_vert;
        free [] z_vert;
        
        // example using AOS 
        double * xyz_vert_aos = new double[nvert*3] ;
        geomData.writeGridPoints_SOA(xyz_vert_aos) ;

        ... do something with vertex data xyz_vert_aos....

        free [] xyz_vert_aos;
        

*/

namespace Opm::GridDataOutput
{
   /**
   * Allows selection of order of verticies in writeConnectivity()
   */
  enum ConnectivityVertexOrder { DUNE = 0 , VTK = 1 } ;

  template< class GridView, unsigned int partitions >
  class SimMeshDataAccessor {
  public:
    /**
     * @brief Construct a SimMeshDataAccessor working on a specific GridView and specialize to a Dune::PartitionSet<>.
     *
     * @param gridView The gridView 
     * @param PartitionSet<> the set of cells from which to extract geometric data
     *
     *  The PartitionSet of the data can be specified from one of:
     *   Dune::Partitions::all
     *   Dune::Partitions::interior
     *   Dune::Partitions::border
     *   Dune::Partitions::overlap
     *   Dune::Partitions::front
     *   Dune::Partitions::ghost
     *   Dune::Partitions::interiorBorder
     *   Dune::Partitions::interiorBorderOverlap
     *   Dune::Partitions::interiorBorderOverlapFront
     *   Dune::Partitions::all
     *
     * N.B. To visualise 'field' data on the extracted grid mesh then the field variable
     *       should contain at least as many vlaues as the mesh has cells  (ncells_) or vertices (nvertices_)
     *       depending on if data is cell centred or vertex centred, respectively.
     *
     *  As we are templated on the Dune::PartitionSet<partitions>, values for ncorners_, nvertices_ and ncells_ cannot change
     *
     *  This class does not work with grids containing polyhedral cells (well, it has not been tested
     *  with this kind of grid data). The user should call polyhedralCellPresent() to test if polyhedral 
     *  cells are present and decide what they want to do before copying data using the data accessor methods.
     */
    explicit SimMeshDataAccessor ( const GridView &gridView,  
                            Dune::PartitionSet<partitions>  dunePartition)
        : gridView_( gridView ),
          dunePartition_(dunePartition)
    {
        dimw_ = GridView::dimension ; // this is an enum
        partition_value_ = dunePartition.value ;
        countEntities() ;
    }

    //! destructor
    ~SimMeshDataAccessor ()
    {
    }

     /**
       Checks for cells that have polyhedral type within the current partition of cells
       
       Returns true if a polyhedral sell is found. If this is the case then this partition 
       is not going to be available for visualisation as this class does not yet handle 
       polyhedral cells.
    */
    bool polyhedralCellPresent() 
    {
        for (const auto& cit :  elements(gridView_, dunePartition_))
        {
            auto  corner_geom = cit.geometry() ;
            if( Dune::VTK::geometryType( corner_geom.type() ) == Dune::VTK::polyhedron )
            {
                return true ;
            }
       } 
       return false;
    }

    /** 
        Count the vertices, cells and corners.

        Count all the vertices ( the Dune::Partitions::all partition ) as then we do not need to renumber
        the vertices as all the subsets use references to the full set.
    */ 
    void countEntities(  )
    {
        // We include all the vertices for this ranks partition
        const auto& vert_partition_it = vertices(gridView_, Dune::Partitions::all);
        nvertices_ = std::distance(vert_partition_it.begin(), vert_partition_it.end());

        const auto& cell_partition_it = elements(gridView_, dunePartition_);
        ncells_ = std::distance(cell_partition_it.begin(), cell_partition_it.end());

        ncorners_ = 0 ;
        for (const auto& cit : cell_partition_it)
        {
            auto  corner_geom = cit.geometry() ;
            ncorners_ += corner_geom.corners() ;
        }
    }

     /**
       Write the positions of vertices - directly to the pointers given in paramaters 

       Returns the number of vertices written
    */
    template <typename T> 
    long  writeGridPoints( T*  x_inout,  T*  y_inout, T* z_inout )
    {
        long i = 0 ;
        if (dimw_ == 3) {
            for (const auto& vit :   vertices(gridView_, Dune::Partitions::all) )
            {
                // if (i < nvertices_) // As we are templated on the Dune::PartitionSet<partitions>, this cannot change
                auto xyz_local = vit.geometry().corner(0);  // verticies only have one corner
                x_inout[i] = static_cast<T>(xyz_local[0]) ;
                y_inout[i] = static_cast<T>(xyz_local[1]) ;
                z_inout[i] = static_cast<T>(xyz_local[2]) ; 
                i++ ;
            }
        } else if (dimw_ == 2)  {
            for (const auto& vit :   vertices(gridView_, Dune::Partitions::all) )
            {      
                // if (i < nvertices_) // As we are templated on the Dune::PartitionSet<partitions>, this cannot change
                auto xyz_local = vit.geometry().corner(0);  // verticies only have one corner
                x_inout[i] = static_cast<T>(xyz_local[0]) ;
                y_inout[i] = static_cast<T>(xyz_local[1]) ;
                z_inout[i] = static_cast<T>(0.0);
                i++ ;
            }
        }
        return i ;
    }

    /**
     Write positions of vertices as array of structures : x,y,z,x,y,z,x,y,z,...

     Returns the number of vertices written
    */
    template <typename T> 
    long writeGridPoints_AOS( T*  xyz_inout )
    {
        long i = 0 ;
        if (dimw_ == 3) {
            for (const auto& vit :   vertices(gridView_, Dune::Partitions::all))
            {
                auto xyz_local = vit.geometry().corner(0);
                xyz_inout[i++] = static_cast<T>(xyz_local[0]) ;
                xyz_inout[i++] = static_cast<T>(xyz_local[1]) ;
                xyz_inout[i++] = static_cast<T>(xyz_local[2]) ;
                
            }
        } else if (dimw_ == 2)  {
            for (const auto& vit :   vertices(gridView_, Dune::Partitions::all))
            {
                auto xyz_local = vit.geometry().corner(0);
                xyz_inout[i++] = static_cast<T>(xyz_local[0]) ;
                xyz_inout[i++] = static_cast<T>(xyz_local[1]) ;
                xyz_inout[i++] = static_cast<T>(0.0) ;  
            }
        }
        return ( (i) / 3 );
    }

    /**
     Write positions of vertices as structure of arrays  : x,x,x,...,y,y,y,...,z,z,z,...

     Returns the number of vertices written
    */
    template <typename T> 
    long writeGridPoints_SOA( T*  xyz_inout )
    {
        long i = 0 ;
        // Get offsets into structure
        T * xyz_inout_y = xyz_inout + nvertices_ ;
        T * xyz_inout_z = xyz_inout + (2*nvertices_) ;

        if (dimw_ == 3) {
            for (const auto& vit :   vertices(gridView_, Dune::Partitions::all))
            {     
                auto xyz_local = vit.geometry().corner(0);
                xyz_inout[i] = static_cast<T>(xyz_local[0]) ;
                xyz_inout_y[i]= static_cast<T>(xyz_local[1]) ;
                xyz_inout_z[i] = static_cast<T>(xyz_local[2]) ;
                i++ ;
            }
        } else if (dimw_ == 2)  {
            for (const auto& vit :   vertices(gridView_, Dune::Partitions::all))
            {     
                auto xyz_local = vit.geometry().corner(0);
                xyz_inout[i] = static_cast<T>(xyz_local[0]) ;
                xyz_inout_y[i]= static_cast<T>(xyz_local[1]) ;
                xyz_inout_z[i] = static_cast<T>(0.0);  
                i++ ;
            }
        }
        return (i) ;
    }

    /**
    * Write the connectivity array - directly to the pointer given in paramater 1
      Reorders the indecies as selected either in DUNE order or  VTK order.

      Returns the number of corner indices written.
    */
    template <typename I>
    long writeConnectivity(I * connectivity_inout, ConnectivityVertexOrder whichOrder)
    {
        long i = 0 ;
        if ( whichOrder == DUNE )  {
            // DUNE order
            for (const auto& cit :  elements(gridView_, dunePartition_))
            {  
              auto cell_corners = cit.geometry().corners() ;
              for( auto vx = 0; vx < cell_corners; ++ vx )  
              {
                  const int vxIdx = gridView_.indexSet().subIndex( cit, vx, 3 );
                  connectivity_inout[i + vx] = vxIdx ;
              }
              i += cell_corners ; 
            }
        } else {
            // VTK order
            for (const auto& cit :  elements(gridView_, dunePartition_))
            {
              auto cell_corners = cit.geometry().corners() ;
              for( auto vx = 0; vx < cell_corners; ++ vx )  
              {
                  const int vxIdx = gridView_.indexSet().subIndex( cit, vx, 3 );
                  int vtkOrder ; 
                  vtkOrder = Dune::VTK::renumber(cit.type(), vx) ;
                  connectivity_inout[i + vtkOrder] = vxIdx ;
              }
              i += cell_corners ; 
            }
        }
        return (i) ;
    }

    /**
    * Write the offsets values  - directly to the pointer given in paramater 1
      Returns the number of offset values written, which should be 1 greater than ncells_ 
      or -1 if an error was detected
    */
    template <typename I>
    long writeOffsetsCells( I* offsets_inout  )
    {
        long i = 1 ;
        offsets_inout[0] = 0 ;
        for (const auto& cit :  elements(gridView_, dunePartition_))
        {  
            auto cell_corners = cit.geometry().corners() ;
            offsets_inout[i] = offsets_inout[i-1] +  cell_corners ;
            i++ ;
        }
        return (i) ;  // This should be 1 greater than ncells_
    }

    /**
    * Write the Cell types array - directly to the pointer given in paramater 1
    */
    template <typename I>
    long writeCellTypes( I* types_inout)
    {
        int i = 0 ;
        for (const auto& cit :  elements(gridView_, dunePartition_))
        {
              I vtktype = static_cast<I>(Dune::VTK::geometryType(cit.type()));
              types_inout[i++] = vtktype ;
        }
        return (i) ;
    }

    std::string getPartitionTypeString (  )
    {
        if (this->dunePartition_ ==  Dune::Partitions::all)
            return (std::string("Dune::Partitions::all")) ;
        if (this->dunePartition_ ==  Dune::Partitions::interior)
            return (std::string("Dune::Partitions::interior")) ;
        if (this->dunePartition_ ==  Dune::Partitions::interiorBorder)
            return (std::string("Dune::Partitions::interiorBorder")) ;
        if (this->dunePartition_ ==  Dune::Partitions::interiorBorderOverlap)
            return (std::string("Dune::Partitions::interiorBorderOverlap")) ;
        if (this->dunePartition_ ==  Dune::Partitions::front)
            return (std::string("Dune::Partitions::front")) ;
        if (this->dunePartition_ ==  Dune::Partitions::interiorBorderOverlapFront)
            return (std::string("Dune::Partitions::InteriorBorderOverlapFront")) ;
        if (this->dunePartition_ ==  Dune::Partitions::border)
            return (std::string("Dune::Partitions::border")) ;
        if (this->dunePartition_ ==  Dune::Partitions::ghost)
            return (std::string("Dune::Partitions::ghost")) ;
        
        return (std::string("Unknown Dune::PartitionSet<>")) ;
    }

    Dune::PartitionSet<partitions> getPartition ( void )
    {
        return ( this->dunePartition_ ) ;
    }

    void   printGridDetails()
    {
        std::cout << "Dune Partition = " << partition_value_ << ", " << getPartitionTypeString()  << std::endl ;
        printNCells() ;
        printNVertices() ;
        printNCorners() ;
    }

    void printNCells()
    {
        std::cout << "ncells = " << ncells_ << std::endl ;
    }

    void printNVertices()
    {
        std::cout << "nvertices = " << nvertices_ << std::endl ;
    }

    void printNCorners()
    {
        std::cout << "ncorners = " << ncorners_ << std::endl ;
    }
    
    int getNCells()
    {
        return(ncells_) ;
    }

    int getNVertices()
    {
        return(nvertices_) ;
    }

    int getNCorners()
    {
        return(ncorners_) ;
    }

    std::string getError()
    {
        return error_strm_.str() ;
    }

    void clearError()
    {
         error_strm_.str("") ;
    }

    bool hasError()
    {
        if ( error_strm_.str().length() > 0 ) 
            return true ;
        else 
            return false ;
    }

  protected:

    GridView gridView_;  // the grid
    
    Dune::PartitionSet<partitions>  dunePartition_ ;
    unsigned int partition_value_ ;

    /**
    Current partition grid information
    */
    int ncells_;
    /**
    Current partition grid information
    */
    int nvertices_;
    /**
    Current partition grid information
    */
    int ncorners_; 
    
    int dimw_ ;  // dimensions of the input grid
    
  private:
    std::stringstream error_strm_ ;

  };

}

#endif 
