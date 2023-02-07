// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef OPM_SIM_MESH_DATA_HPP
#define OPM_SIM_MESH_DATA_HPP

#include <cstring>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <type_traits>
#include <vector>
#include <list>
#include <map>

#include <dune/common/visibility.hh>
#include <dune/common/typetraits.hh>
// #include <dune/common/exceptions.hh>
#include <dune/common/indent.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/path.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/dataarraywriter.hh>
//#include <dune/grid/io/file/vtk/function.hh>
//#include <dune/grid/io/file/vtk/pvtuwriter.hh>
//#include <dune/grid/io/file/vtk/streams.hh>
//#include <dune/grid/io/file/vtk/vtuwriter.hh>

/** @file
    @author Joshua Bowden
    @brief Allows model geometry data to be passed to external code - via copy in std::vector<T> or copy direct to input pointers. 
    Based off the Dune VTKWriter.hh code. 
    
    Example:
        // From the opm-simulators repository
        #include <opm/simulators/utils/GridDataOutput.hpp>
        
        // N.B. does not seem to be able to be allocated with new operator.
        Opm::GridDataOutput::SimMeshDataAccessor geomData(gridView, Dune::VTK::conforming) ;
        // geomData = new SimMeshDataAccessor(gridView, Dune::VTK::conforming);  // this does not compile
        geomData.printGridDetails() ;
        
        int nvert = geomData.getNVertices() ;
        double * x_vert = new double[nvert] ;
        double * y_vert = new double[nvert] ;
        double * z_vert = new double[nvert] ;
        geomData.writeGridPoints(x_vert,y_vert,z_vert) ;
        
        ... do something with vertex data ....
        
        free [] x_vert;
        free [] y_vert;
        free [] z_vert;
        
        
    TODO: Extend to support Grid field data output.
*/

namespace Opm::GridDataOutput
{


  /**
   * @brief Writer for the ouput of grid functions in the vtk format.
   * @ingroup VTK
   *
   * Writes arbitrary grid functions (living on cells or vertices of a grid)
   * to a file suitable for easy visualization with
   * <a href="http://public.kitware.com/VTK/">The Visualization Toolkit (VTK)</a>.
   */
  template< class GridView >
  class SimMeshDataAccessor {

    // extract types
    typedef typename GridView::Grid Grid;
    typedef typename GridView::ctype DT;
    enum { n = GridView::dimension };
    enum { w = GridView::dimensionworld };

    typedef typename GridView::template Codim< 0 >::Entity Cell;
    typedef typename GridView::template Codim< n >::Entity Vertex;
    typedef Cell Entity;

    typedef typename GridView::IndexSet IndexSet;

    static const Dune::PartitionIteratorType VTK_Partition = Dune::InteriorBorder_Partition;
    //static const PartitionIteratorType VTK_Partition = All_Partition;

    typedef typename GridView::template Codim< 0 >
    ::template Partition< VTK_Partition >::Iterator
    GridCellIterator;
    typedef typename GridView::template Codim< n >
    ::template Partition< VTK_Partition >::Iterator
    GridVertexIterator;

    typedef typename GridCellIterator::Reference EntityReference;

    typedef typename GridView::template Codim< 0 >
    ::Entity::Geometry::LocalCoordinate Coordinate;

    typedef Dune::MultipleCodimMultipleGeomTypeMapper< GridView > VertexMapper;

    // return true if entity should be skipped in Vertex and Corner iterator
    static bool skipEntity( const Dune::PartitionType entityType )
    {
      switch( VTK_Partition )
      {
        // for All_Partition no entity has to be skipped
        case Dune::All_Partition:             return false;
        case Dune::InteriorBorder_Partition:  return ( entityType != Dune::InteriorEntity );
        default: DUNE_THROW(Dune::NotImplemented,"Add check for this partition type");
      }
      return false ;
    }


  public:
    
    //! Iterator over the grids elements
    /**
     * This class iterates over the gridview's elements.  It is the same as
     * the gridview's Codim<0>::Iterator for the InteriorBorder_Partition,
     * except that it add a position() method.
     */
    class CellIterator : public GridCellIterator
    {
    public:
      //! construct a CellIterator from the gridview's Iterator.
      CellIterator(const GridCellIterator & x) : GridCellIterator(x) {}
      //! get the position of the center of the element, in element-local
      //! coordinates
      const Dune::FieldVector<DT,n> position() const
      {
        return Dune::Geo::ReferenceElements<DT,n>::general((*this)->type()).position(0,0);
      }
    };

    CellIterator cellBegin() const
    {
      return gridView_.template begin< 0, VTK_Partition >();
    }

    CellIterator cellEnd() const
    {
      return gridView_.template end< 0, VTK_Partition >();
    }

    //! Iterate over the grid's vertices
    /**
     * This class iterates over the elements, and within the elements over the
     * corners.  If the data mode dm is nonconforming, each vertex is visited
     * once for each element where it is a corner (similar to CornerIterator).
     * If dm is conforming each vertex is visited only once globally, for the
     * first element where it is a corner.  Contrary to CornerIterator, visit
     * the corners of a given element in Dune-ordering.
     *
     * Dereferencing the iterator yields the current entity, and the index of
     * the current corner within that entity is returned by the iterators
     * localindex() method.  Another useful method on the iterator itself is
     * position() which returns the element-local position of the current
     * corner.
     */
    class VertexIterator :
      public Dune::ForwardIteratorFacade<VertexIterator, const Entity, EntityReference, int>
    {
      GridCellIterator git;
      GridCellIterator gend;
      Dune::VTK::DataMode datamode;
      // Index of the currently visited corner within the current element.
      // NOTE: this is in Dune-numbering, in contrast to CornerIterator.
      int cornerIndexDune;
      const VertexMapper & vertexmapper;
      std::vector<bool> visited;
      // in conforming mode, for each vertex id (as obtained by vertexmapper)
      // hold its number in the iteration order (VertexIterator)
      int offset;

      // hide operator ->
      void operator->();
    protected:
      void basicIncrement ()
      {
        if( git == gend )
          return;
        ++cornerIndexDune;
        const int numCorners = git->subEntities(n);
        if( cornerIndexDune == numCorners )
        {
          offset += numCorners;
          cornerIndexDune = 0;

          ++git;
          while( (git != gend) && skipEntity( git->partitionType() ) )
            ++git;
        }
      }
    public:
      VertexIterator(const GridCellIterator & x,
                     const GridCellIterator & end,
                     const Dune::VTK::DataMode & dm,
                     const VertexMapper & vm) :
        git(x), gend(end), datamode(dm), cornerIndexDune(0),
        vertexmapper(vm), visited(vm.size(), false),
        offset(0)
      {
        if (datamode == Dune::VTK::conforming && git != gend)
          visited[vertexmapper.subIndex(*git,cornerIndexDune,n)] = true;
      }
      void increment ()
      {
        switch (datamode)
        {
        case Dune::VTK::conforming :
          while(visited[vertexmapper.subIndex(*git,cornerIndexDune,n)])
          {
            basicIncrement();
            if (git == gend) return;
          }
          visited[vertexmapper.subIndex(*git,cornerIndexDune,n)] = true;
          break;
        case Dune::VTK::nonconforming :
          basicIncrement();
          break;
        }
      }
      bool equals (const VertexIterator & cit) const
      {
        return git == cit.git
               && cornerIndexDune == cit.cornerIndexDune
               && datamode == cit.datamode;
      }
      EntityReference dereference() const
      {
        return *git;
      }
      //! index of vertex within the entity, in Dune-numbering
      int localindex () const
      {
        return cornerIndexDune;
      }
      //! position of vertex inside the entity
      Dune::FieldVector<DT,n> position () const
      {
        return Dune::referenceElement<DT,n>(git->type())
          .position(cornerIndexDune,n);
      }
    };

    VertexIterator vertexBegin () const
    {
      return VertexIterator( gridView_.template begin< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper );
    }

    VertexIterator vertexEnd () const
    {
      return VertexIterator( gridView_.template end< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper );
    }

    //! Iterate over the elements' corners
    /**
     * This class iterates over the elements, and within the elements over the
     * corners.  Each vertex in the grid can be a corner in multiple elements,
     * and is visited once for each element it is associated with.  This class
     * differs from VertexIterator in that it visits the corners of a given
     * element in VTK-ordering, and that it always visits a given vertex once
     * for each element where that vertex is a corner in, independent of the
     * data mode dm.
     *
     * Dereferencing the iterator yields the current entity.  Another useful
     * method on the iterator itself is id(), which returns the number of the
     * current corners associated vertex, in the numbering given by the
     * iteration order of VertexIterator.
     */
    class CornerIterator :
      public Dune::ForwardIteratorFacade<CornerIterator, const Entity, EntityReference, int>
    {
      GridCellIterator git;
      GridCellIterator gend;
      Dune::VTK::DataMode datamode;
      // Index of the currently visited corner within the current element.
      // NOTE: this is in VTK-numbering, in contrast to VertexIterator.
      int cornerIndexVTK;
      const VertexMapper & vertexmapper;
      // in conforming mode, for each vertex id (as obtained by vertexmapper)
      // hold its number in the iteration order of VertexIterator (*not*
      // CornerIterator)
      const std::vector<int> & number;
      // holds the number of corners of all the elements we have seen so far,
      // excluding the current element
      int offset;

      // hide operator ->
      void operator->();
    public:
      CornerIterator(const GridCellIterator & x,
                     const GridCellIterator & end,
                     const Dune::VTK::DataMode & dm,
                     const VertexMapper & vm,
                     const std::vector<int> & num) :
        git(x), gend(end), datamode(dm), cornerIndexVTK(0),
        vertexmapper(vm),
        number(num), offset(0) {}
      void increment ()
      {
        if( git == gend )
          return;
        ++cornerIndexVTK;
        const int numCorners = git->subEntities(n);
        if( cornerIndexVTK == numCorners )
        {
          offset += numCorners;
          cornerIndexVTK = 0;

          ++git;
          while( (git != gend) && skipEntity( git->partitionType() ) )
            ++git;
        }
      }
      bool equals (const CornerIterator & cit) const
      {
        return git == cit.git
               && cornerIndexVTK == cit.cornerIndexVTK
               && datamode == cit.datamode;
      }
      EntityReference dereference() const
      {
        return *git;
      }
      //! Process-local consecutive zero-starting vertex id
      /**
       * This method returns the number of this corners associated vertex, in
       * the numbering given by the iteration order of VertexIterator.
       */
      int id () const
      {
        switch (datamode)
        {
        case Dune::VTK::conforming :
          return
            number[vertexmapper.subIndex(*git,Dune::VTK::renumber(*git,cornerIndexVTK),
                                    n)];
        case Dune::VTK::nonconforming :
          return offset + Dune::VTK::renumber(*git,cornerIndexVTK);
        default :
          DUNE_THROW(Dune::IOError,"SimMeshDataAccessor: unsupported DataMode" << datamode);
        }
      }
    };

    CornerIterator cornerBegin () const
    {
      return CornerIterator( gridView_.template begin< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper, number );
    }

    CornerIterator cornerEnd () const
    {
      return CornerIterator( gridView_.template end< 0, VTK_Partition >(),
                             gridView_.template end< 0, VTK_Partition >(),
                             datamode, *vertexmapper, number );
    }

  public:
    /**
     * @brief Construct a SimMeshDataAccessor working on a specific GridView.
     *
     *
     * @param gridView The gridView the grid functions live on. (E. g. a LevelGridView.)
     * @param dm The data mode.
     * @param coordPrecision the precision with which to write out the coordinates
     */
    explicit SimMeshDataAccessor ( const GridView &gridView,
                         Dune::VTK::DataMode dm = Dune::VTK::conforming)
      : gridView_( gridView ),
        datamode( dm ),
       /* coordPrec (coordPrecision),*/
        polyhedralCellsPresent_( checkForPolyhedralCells() )
    { 
        
        this->setup_called_bool_ = false ;
        this->setupGeomData() ;
    
    }


    //! clear list of registered functions
    void clear ()
    {
      //celldata.clear();
      //vertexdata.clear();
    }

    //! get the precision with which coordinates are written out
    /*Dune::VTK::Precision coordPrecision() const
    { return coordPrec; }
*/
    //! destructor
    ~SimMeshDataAccessor ()
    {
      this->clear();
    }

    

    //! Get the detais of the grid data 
    void setupGeomData ( void )
    {
      Dune::VTK::FileType fileType =
        (n == 1) ? Dune::VTK::polyData : Dune::VTK::unstructuredGrid;

      // VTK::VTUWriter writer(s, outputtype, fileType);

      // Grid characteristics
      vertexmapper = new VertexMapper( gridView_, Dune::mcmgVertexLayout() );
      if (datamode == Dune::VTK::conforming)
      {
        number.resize(vertexmapper->size());
        for (std::vector<int>::size_type i=0; i<number.size(); i++) 
            number[i] = -1;
      }
      countEntities(nvertices, ncells, ncorners);

      this->setup_called_bool_ = true ;

      delete vertexmapper; number.clear();
    }

    void writeGeometryData() {
      // PointData fields
      //writeVertexData(writer);

      // CellData fields
      //writeCellData(writer);

      // x,y,z vertices coordinates
      // writeGridPoints(writer);

      // Cells -connectivity, offested, types, polygon extras if required
      // writeGridCells(writer);
    }


    std::string getTypeString() const
    {
      if (n==1)
        return "PolyData";
      else
        return "UnstructuredGrid";
    }

    //! count the vertices, cells and corners
    void countEntities(int &nvertices_, int &ncells_, int &ncorners_)
    {
      nvertices_ = 0;
      ncells_ = 0;
      ncorners_ = 0;
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        ncells_++;
        // because of the use of vertexmapper->map(), this iteration must be
        // in the order of Dune's numbering.
        const int subEntities = it->subEntities(n);
        for (int i=0; i<subEntities; ++i)
        {
          ncorners_++;
          if (datamode == Dune::VTK::conforming)
          {
            int alpha = vertexmapper->subIndex(*it,i,n);
            if (number[alpha]<0)
              number[alpha] = nvertices_++;
          }
          else
          {
            nvertices_++;
          }
        }
      }
    }

/*
    template<typename T>
    std::tuple<std::string,std::string> getDataNames(const T& data) const
    {
      std::string scalars = "";
      for (auto it = data.begin(),
             end = data.end();
           it != end;
           ++it)
        if (it->fieldInfo().type() == VTK::FieldInfo::Type::scalar)
          {
            scalars = it->name();
            break;
          }

      std::string vectors = "";
      for (auto it = data.begin(),
             end = data.end();
           it != end;
           ++it)
        if (it->fieldInfo().type() == VTK::FieldInfo::Type::vector)
          {
            vectors = it->name();
            break;
          }
      return std::make_tuple(scalars,vectors);
    }

    template<typename Data, typename Iterator>
    void writeData(VTK::VTUWriter& writer, const Data& data, const Iterator begin, const Iterator end, int nentries)
    {
      for (auto it = data.begin(),
             iend = data.end();
           it != iend;
           ++it)
      {
        const auto& f = *it;
        VTK::FieldInfo fieldInfo = f.fieldInfo();
        std::size_t writecomps = fieldInfo.size();
        switch (fieldInfo.type())
          {
          case VTK::FieldInfo::Type::scalar:
            break;
          case VTK::FieldInfo::Type::vector:
            // vtk file format: a vector data always should have 3 comps (with
            // 3rd comp = 0 in 2D case)
            if (writecomps > 3)
              DUNE_THROW(IOError,"Cannot write VTK vectors with more than 3 components (components was " << writecomps << ")");
            writecomps = 3;
            break;
          case VTK::FieldInfo::Type::tensor:
            DUNE_THROW(NotImplemented,"VTK output for tensors not implemented yet");
          }
        std::shared_ptr<VTK::DataArrayWriter> p
          (writer.makeArrayWriter(f.name(), writecomps, nentries, fieldInfo.precision()));
        if(!p->writeIsNoop())
          for (Iterator eit = begin; eit!=end; ++eit)
          {
            const Entity & e = *eit;
            f.bind(e);
            f.write(eit.position(),*p);
            f.unbind();
            // vtk file format: a vector data always should have 3 comps
            // (with 3rd comp = 0 in 2D case)
            for (std::size_t j=fieldInfo.size(); j < writecomps; ++j)
              p->write(0.0);
          }
      }
    }

    //! write cell data
    virtual void writeCellData(VTK::VTUWriter& writer)
    {
      if(celldata.size() == 0)
        return;

      std::string scalars, vectors;
      std::tie(scalars,vectors) = getDataNames(celldata);

      writer.beginCellData(scalars, vectors);
      writeData(writer,celldata,cellBegin(),cellEnd(),ncells);
      writer.endCellData();
    }

    //! write vertex data
    virtual void writeVertexData(VTK::VTUWriter& writer)
    {
      if(vertexdata.size() == 0)
        return;

      std::string scalars, vectors;
      std::tie(scalars,vectors) = getDataNames(vertexdata);

      writer.beginPointData(scalars, vectors);
      writeData(writer,vertexdata,vertexBegin(),vertexEnd(),nvertices);
      writer.endPointData();
    }*/

    //! write the positions of vertices
    template <typename T>
    void writeGridPoints( std::vector<T> & x_inout,  std::vector<T> & y_inout, std::vector<T> & z_inout )
    {
        int dimw=w;
        VertexIterator vEnd = vertexEnd();
        for (VertexIterator vit=vertexBegin(); vit!=vEnd; ++vit)
        {          
          x_inout.push_back(  (T) (*vit).geometry().corner(vit.localindex())[0] );
          y_inout.push_back(  (T) (*vit).geometry().corner(vit.localindex())[1] );
          if (dimw == 3) 
            z_inout.push_back( (T) (*vit).geometry().corner(vit.localindex())[2] ); 
          else 
            z_inout.push_back((T)  0.0);  
        }
    }
    
    //! write the positions of vertices - directly to the pointers given in paramaters 1
    template <typename T> 
    void writeGridPoints( T*  x_inout,  T*  y_inout, T* z_inout )
    {
        int dimw=w;
        VertexIterator vEnd = vertexEnd();
        int i = 0 ;
        T zero_val = 0.0 ;
        for (VertexIterator vit=vertexBegin(); vit!=vEnd; ++vit)
        {          
          x_inout[i] = static_cast<T>((*vit).geometry().corner(vit.localindex())[0]) ;
          y_inout[i] = static_cast<T>((*vit).geometry().corner(vit.localindex())[1]) ;
          if (dimw == 3) 
            z_inout[i] = static_cast<T>((*vit).geometry().corner(vit.localindex())[2]) ; 
          else 
            z_inout[i] = zero_val;  
        
          i++ ;
        }
    }

    //! write the connectivity array
    template <typename T, typename U, typename V, typename W>
    void writeGridCells(std::vector<T> & connectivity_inout,  std::vector<U> & offsets_inout, std::vector<V> & types_inout, std::vector<W> & polyfaces_inout, std::vector<W> & polyfacesoffsets_inout  )
    {
      // connectivity
       for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
              connectivity_inout.push_back(static_cast<T>(it.id())) ;

      // offsets
      U offset = 0;
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        offset += static_cast<U>(it->subEntities(n));
        offsets_inout.push_back(offset) ;
      }


      // types
      if (n>1)
      {
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          int vtktype = Dune::VTK::geometryType(it->type());
            types_inout.push_back(static_cast<V>(vtktype));
        }
        
        // if polyhedron cells found also cell faces need to be written
        if( polyhedralCellsPresent_ )
        {
          writeCellFaces( polyfaces_inout, polyfacesoffsets_inout );
        }
      }
    }
    
    /**
    * write the connectivity array - make a copy into a std::vector<T>
    */
    template <typename T>
    void writeConnectivity(std::vector<T> & connectivity_inout )
    {
    
      // connectivity
       int i = 0 ;
       for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
       {
           T connect_data = static_cast<T>(it.id()) ;
           connectivity_inout.push_back( connect_data );
       }

    }
    
    /**
    * write the offsets values - make a copy into a std::vector<T>
    */
    template <typename T>
    void writeOffsetsCells(std::vector<T> & offsets_inout  )
    {
      // offsets
      T offset = 0;
      int i = 0 ;
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        offset += static_cast<T>(it->subEntities(n));
        offsets_inout.push_back( offset ) ;
      }
    }
    
    /**
    * write the Cell types array - make a copy into a std::vector<T>
    */
    template <typename T>
    void writeCellTypes( std::vector<T> & types_inout )
    {
      int i = 0 ;
      // types
      if (n>1)
      {
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          T vtktype = static_cast<T>(Dune::VTK::geometryType(it->type()));
          types_inout.push_back( vtktype ) ;
        }
      }
    }
    

    /**
    * write the connectivity array - directly to the pointer given in paramater 1
    */
    template <typename T>
    void writeConnectivity(T * connectivity_inout )
    {
    
      // connectivity
       int i = 0 ;
       for (CornerIterator it=cornerBegin(); it!=cornerEnd(); ++it)
       {
           T connect_data = static_cast<T>(it.id()) ;
           connectivity_inout[i++] = connect_data ;
       }

    }
    
    /**
    * write the offsets values  - directly to the pointer given in paramater 1
    */
    template <typename T>
    void writeOffsetsCells( T* offsets_inout  )
    {
      // offsets
      T offset = 0;
      int i = 0 ;
      for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
      {
        offset += static_cast<T>(it->subEntities(n));
        offsets_inout[i++] = offset ;
      }
    }
    
    /**
    * write the Cell types array - directly to the pointer given in paramater 1
    */
    template <typename T>
    void writeCellTypes( T* types_inout )
    {
      int i = 0 ;
      // types
      if (n>1)
      {
        for (CellIterator it=cellBegin(); it!=cellEnd(); ++it)
        {
          T vtktype = static_cast<T>(Dune::VTK::geometryType(it->type()));
          types_inout[i++] = vtktype ;

        }
      }
    }
    
    
    /**
    * write the connectivity array
    * N.B. W* polyfaces_inout, X* polyfacesoffsets_inout are only needed if hasPolyhedralCells() returns true
    */
    template <typename T>
    void writePolyCells( T* polyfaces_inout, T* polyfacesoffsets_inout  )
    {
    
      int i = 0 ;
      // types
      if (n>1)
      {
        // if polyhedron cells found also cell faces need to be written
        if( polyhedralCellsPresent_ )
        {
          writeCellFaces( polyfaces_inout, polyfacesoffsets_inout );
        }
      }

    }
    

    bool checkForPolyhedralCells() const
    {
      // check if polyhedron cells are present
      for( const auto& geomType : gridView_.indexSet().types( 0 ) )
      {
        if( Dune::VTK::geometryType( geomType ) == Dune::VTK::polyhedron )
        {
          return true;
        }
      }
      return false;
    }

    //! write the polyhedral cell faces array - make a copy into a std::vector<T> method arguments
    template <typename T>
    void writeCellFaces( std::vector<T> & polyfaces_inout, std::vector<T> & polyfacesoffsets_inout )
    {
      std::shared_ptr< std::pair< std::vector<int>, std::vector<int> > > faceVertices_;
      if( ! faceVertices_ )
      {
        faceVertices_.reset( new std::pair< std::vector<T>, std::vector<T> > () );
        // fill face vertex structure
        fillFaceVertices( cornerBegin(), cornerEnd(), gridView_.indexSet(),
                          faceVertices_->first, faceVertices_->second );
      }

      std::vector< T >& faces = faceVertices_->first;
      std::vector< T >& faceOffsets = faceVertices_->second;
      assert( int(faceOffsets.size()) == ncells );

      {
          T face_typew ;
          for( const auto& face : faces )
          {
              face_typew = static_cast<T>(face);
              polyfaces_inout.push_back(face) ;
          }
      }

      {
          T faceoffset_typeX ;
          for( const auto& offset : faceOffsets )
          {
              faceoffset_typeX = static_cast<T>(offset);
              polyfacesoffsets_inout.push_back(faceoffset_typeX) ;
          }
          // clear face vertex structure
          faceVertices_.reset();
      }
    }

    template <class CornerIterator, class IndexSet, class T>
    inline void fillFaceVertices( CornerIterator it,
                           const CornerIterator end,
                           const IndexSet& indexSet,
                           std::vector<T>& faces,
                           std::vector<T>& faceOffsets )
    {
      if( n == 3 && it != end )
      {
        // clear output arrays
        faces.clear();
        faces.reserve( 15 * ncells );
        faceOffsets.clear();
        faceOffsets.reserve( ncells );

        int offset = 0;

        Cell element = *it;
        int elIndex = indexSet.index( element );
        std::vector< T > vertices;
        vertices.reserve( 30 );
        for( ; it != end; ++it )
        {
          const Cell& cell = *it ;
          const int cellIndex = indexSet.index( cell ) ;
          if( elIndex != cellIndex )
          {
            fillFacesForElement( element, indexSet, vertices, offset, faces, faceOffsets );

            vertices.clear();
            element = cell ;
            elIndex = cellIndex ;
          }
          vertices.push_back( it.id() );
        }

        // fill faces for last element
        fillFacesForElement( element, indexSet, vertices, offset, faces, faceOffsets );
      }
    }

    template <class Entity, class IndexSet, class T>
    static void fillFacesForElement( const Entity& element,
                                     const IndexSet& indexSet,
                                     const std::vector<T>& vertices,
                                     T& offset,
                                     std::vector<T>& faces,
                                     std::vector<T>& faceOffsets )
    {
      const int dim = n;

      std::map< T, T > vxMap;

      // get number of local faces
      const int nVertices = element.subEntities( dim );
      for( int vx = 0; vx < nVertices; ++ vx )
      {
        const int vxIdx = indexSet.subIndex( element, vx, dim );
        vxMap[ vxIdx ] = vertices[ vx ];
      }

      // get number of local faces
      const int nFaces = element.subEntities( 1 );
      // store number of faces for current element
      faces.push_back( nFaces );
      ++offset;
      // extract each face as a set of vertex indices
      for( int fce = 0; fce < nFaces; ++ fce )
      {
        // obtain face
        const auto face = element.template subEntity< 1 > ( fce );

        // get all vertex indices from current face
        const int nVxFace = face.subEntities( dim );
        faces.push_back( nVxFace );
        ++offset ;
        for( int i=0; i<nVxFace; ++i )
        {
          const T vxIndex = indexSet.subIndex( face, i, dim );
          assert( vxMap.find( vxIndex ) != vxMap.end() );
          faces.push_back( vxMap[ vxIndex ] );
          ++offset ;
        }
      }

      // store face offset for each element
      faceOffsets.push_back( offset );
    }
    
  void   printGridDetails()
  {
      printNCells() ;
      printNVertices() ;
      printNCorners() ;
      std::cout << "Mesh Type = " << getTypeString()  << std::endl ;
  }
    
  void printNCells()
  {
      std::cout << "ncells = " << ncells << std::endl ;
  }
  
  void printNVertices()
  {
      std::cout << "nvertices = " << nvertices << std::endl ;
  }
  
  void printNCorners()
  {
      std::cout << "ncorners = " << ncorners << std::endl ;
  }
  
  int getNCells()
  {
      return(ncells) ;
  }
  
  int getNVertices()
  {
      return(nvertices) ;
  }
  
  int getNCorners()
  {
      return(ncorners) ;
  }
  
  bool hasPolyhedralCells( void )
  {
      return (polyhedralCellsPresent_) ;
  }
  
  
  protected:
    // the list of registered functions
    // std::list<VTKLocalFunction> celldata;
    // std::list<VTKLocalFunction> vertexdata;

    // the grid
    GridView gridView_;

    // temporary grid information
    int ncells;
    int nvertices;
    int ncorners;
    
  private:
    VertexMapper* vertexmapper;
    // in conforming mode, for each vertex id (as obtained by vertexmapper)
    // hold its number in the iteration order (VertexIterator)
    std::vector<int> number;
    Dune::VTK::DataMode datamode;
    
    // Dune::VTK::Precision coordPrec;

    // true if polyhedral cells are present in the grid
    const bool polyhedralCellsPresent_;
    
    // Used to check that class data has been initialized
    bool setup_called_bool_ ;

    // pointer holding face vertex connectivity if needed
    // std::shared_ptr< std::pair< std::vector<int>, std::vector<int> > > faceVertices_;

  protected:
    Dune::VTK::OutputType outputtype;
  };

}

#endif 
