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
#ifndef EWOMS_PARALLELSERIALOUTPUT_HH
#define EWOMS_PARALLELSERIALOUTPUT_HH

//#if HAVE_OPM_GRID
#include <dune/grid/common/p2pcommunicator.hh>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/common/gridenums.hh>
//#else
//#error "This header needs the opm-grid module."
//#endif

#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <stdexcept>

namespace Ewoms
{
    template < class GridManager >
    class CollectDataToIORank
    {
    public:
        typedef typename GridManager :: Grid  Grid;
        typedef typename Grid :: CollectiveCommunication  CollectiveCommunication;

        // global id
        class GlobalCellIndex
        {
            int globalId_;
            int localIndex_;
            bool isInterior_;
        public:
            GlobalCellIndex() : globalId_(-1), localIndex_(-1), isInterior_(true) {}
            void setGhost() { isInterior_ = false; }

            void setId( const int globalId ) { globalId_ = globalId; }
            void setIndex( const int localIndex ) { localIndex_ = localIndex; }

            int localIndex () const { return localIndex_; }
            int id () const { return globalId_; }
            bool isInterior() const { return isInterior_; }
        };

        typedef typename Dune::PersistentContainer< Grid, GlobalCellIndex > GlobalIndexContainer;

        static const int dimension = Grid :: dimension ;

        typedef typename Grid :: LeafGridView             GridView;
        typedef GridView                                  AllGridView;

        typedef Dune :: Point2PointCommunicator< Dune :: SimpleMessageBuffer > P2PCommunicatorType;
        typedef typename P2PCommunicatorType :: MessageBufferType MessageBufferType;

        typedef std::vector< GlobalCellIndex > LocalIndexMapType;

        typedef std::vector<int> IndexMapType;
        typedef std::vector< IndexMapType > IndexMapStorageType;

        class DistributeIndexMapping : public P2PCommunicatorType::DataHandleInterface
        {
        protected:
            const std::vector<int>& distributedGlobalIndex_;
            IndexMapType& localIndexMap_;
            IndexMapStorageType& indexMaps_;
            std::map< const int, const int > globalPosition_;
#ifndef NDEBUG
            std::set< int > checkPosition_;
#endif

        public:
            DistributeIndexMapping( const std::vector<int>& globalIndex,
                                    const std::vector<int>& distributedGlobalIndex,
                                    IndexMapType& localIndexMap,
                                    IndexMapStorageType& indexMaps )
            : distributedGlobalIndex_( distributedGlobalIndex ),
              localIndexMap_( localIndexMap ),
              indexMaps_( indexMaps ),
              globalPosition_()
            {
                const size_t size = globalIndex.size();
                // create mapping globalIndex --> localIndex
                for ( size_t index = 0; index < size; ++index )
                {
                    globalPosition_.insert( std::make_pair( globalIndex[ index ], index ) );
                }

                // on I/O rank we need to create a mapping from local to global
                if( ! indexMaps_.empty() )
                {
                    // for the ioRank create a localIndex to index in global state map
                    IndexMapType& indexMap = indexMaps_.back();
                    const size_t localSize = localIndexMap_.size();
                    indexMap.resize( localSize );
                    for( size_t i=0; i<localSize; ++i )
                    {
                        const int id = distributedGlobalIndex_[ localIndexMap_[ i ] ];
                        indexMap[ i ] = globalPosition_[ id ] ;
#ifndef NDEBUG
                        assert( checkPosition_.find( id ) == checkPosition_.end() );
                        checkPosition_.insert( id );
#endif
                    }
                }
            }

            void pack( const int link, MessageBufferType& buffer )
            {
                // we should only get one link
                if( link != 0 ) {
                    OPM_THROW(std::logic_error,"link in method pack is not 0 as execpted");
                }

                // pack all interior global cell id's
                const int size = localIndexMap_.size();
                buffer.write( size );

                for( int index = 0; index < size; ++index )
                {
                    const int globalIdx = distributedGlobalIndex_[ localIndexMap_[ index ] ];
                    buffer.write( globalIdx );
                }
            }

            void unpack( const int link, MessageBufferType& buffer )
            {
                // get index map for current link
                IndexMapType& indexMap = indexMaps_[ link ];
                assert( ! globalPosition_.empty() );

                // unpack all interior global cell id's
                int numCells = 0;
                buffer.read( numCells );
                indexMap.resize( numCells );
                for( int index = 0; index < numCells; ++index )
                {
                    int globalId = -1;
                    buffer.read( globalId );
                    assert( globalPosition_.find( globalId ) != globalPosition_.end() );
                    indexMap[ index ] = globalPosition_[ globalId ];
#ifndef NDEBUG
                    assert( checkPosition_.find( globalId ) == checkPosition_.end() );
                    checkPosition_.insert( globalId );
#endif
                }
            }
        };

        enum { ioRank = 0 };

        static const bool needsReordering = ! std::is_same<
            typename GridManager::Grid, typename GridManager::EquilGrid > :: value ;

        CollectDataToIORank( const GridManager& gridManager )
            : toIORankComm_( ),
              isIORank_( gridManager.grid().comm().rank() == ioRank ),
              isParallel_( gridManager.grid().comm().size() > 1 )
        {
            // index maps only have to be build when reordering is needed
            if( ! needsReordering && ! isParallel_ )
            {
                return ;
            }

            const CollectiveCommunication& comm = gridManager.grid().comm();

            {
                std::set< int > send, recv;
                // the I/O rank receives from all other ranks
                if( isIORank() )
                {
                    typedef typename GridManager::EquilGrid::LeafGridView EquilGridView;
                    const EquilGridView equilGridView = gridManager.equilGrid().leafGridView() ;

                    typedef Dune::MultipleCodimMultipleGeomTypeMapper< EquilGridView, Dune::MCMGElementLayout>
                        EquilElementMapper;

                    EquilElementMapper equilElemMapper( equilGridView );

                    // the I/O rank needs a picture of the global grid, here we
                    // use equilGrid which represents a view on the global grid
                    const size_t globalSize = gridManager.equilGrid().leafGridView().size( 0 );
                    // reserve memory
                    globalCartesianIndex_.resize(globalSize, -1);

                    // loop over all elements (global grid) and store Cartesian index
                    auto elemIt = gridManager.equilGrid().leafGridView().template begin<0>();
                    const auto& elemEndIt = gridManager.equilGrid().leafGridView().template end<0>();
                    for (; elemIt != elemEndIt; ++elemIt) {
                        int elemIdx = equilElemMapper.index(*elemIt );
                        int cartElemIdx = gridManager.equilCartesianIndexMapper().cartesianIndex(elemIdx);
                        globalCartesianIndex_[elemIdx] = cartElemIdx;
                    }

                    for(int i=0; i<comm.size(); ++i)
                    {
                        if( i != ioRank )
                        {
                            recv.insert( i );
                        }
                    }
                }
                else // all other simply send to the I/O rank
                {
                    send.insert( ioRank );
                }

                localIndexMap_.clear();
                const size_t gridSize = gridManager.grid().size( 0 );
                localIndexMap_.reserve( gridSize );

                // store the local Cartesian index
                IndexMapType distributedCartesianIndex;
                distributedCartesianIndex.resize(gridSize, -1);

                typedef typename GridManager::GridView LocalGridView;
                const LocalGridView localGridView = gridManager.gridView() ;

                typedef Dune::MultipleCodimMultipleGeomTypeMapper< LocalGridView, Dune::MCMGElementLayout>
                    ElementMapper;

                ElementMapper elemMapper( localGridView );

                for( auto it = localGridView.template begin< 0, Dune::Interior_Partition >(),
                     end = localGridView.template end< 0, Dune::Interior_Partition >(); it != end; ++it )
                {
                    const auto element = *it ;
                    int elemIdx = elemMapper.index( element );
                    distributedCartesianIndex[elemIdx] = gridManager.cartesianIndex( elemIdx );

                    // only store interior element for collection
                    assert( element.partitionType() == Dune :: InteriorEntity );

                    localIndexMap_.push_back( elemIdx );
                }

                // insert send and recv linkage to communicator
                toIORankComm_.insertRequest( send, recv );

                if( isIORank() )
                {
                    // need an index map for each rank
                    indexMaps_.clear();
                    indexMaps_.resize( comm.size() );
                }

                // distribute global id's to io rank for later association of dof's
                DistributeIndexMapping distIndexMapping( globalCartesianIndex_, distributedCartesianIndex, localIndexMap_, indexMaps_ );
                toIORankComm_.exchange( distIndexMapping );
            }
        }

        template <class BufferList>
        class PackUnPackOutputBuffers : public P2PCommunicatorType::DataHandleInterface
        {
            BufferList& bufferList_;

            const IndexMapType& localIndexMap_;
            const IndexMapStorageType& indexMaps_;

        public:
            PackUnPackOutputBuffers( BufferList& bufferList,
                                     const IndexMapType& localIndexMap,
                                     const IndexMapStorageType& indexMaps,
                                     const size_t globalSize,
                                     const bool isIORank )
            : bufferList_( bufferList ),
              localIndexMap_( localIndexMap ),
              indexMaps_( indexMaps )
            {
                if( isIORank )
                {
                    MessageBufferType buffer;
                    pack( 0, buffer );
                    // resize all buffers
                    for (auto it = bufferList_.begin(), end = bufferList_.end(); it != end; ++it )
                    {
                      it->second->resize( globalSize );
                    }

                    // the last index map is the local one
                    doUnpack( indexMaps.back(), buffer );
                }
            }

            // pack all data associated with link
            void pack( const int link, MessageBufferType& buffer )
            {
                // we should only get one link
                if( link != 0 ) {
                    OPM_THROW(std::logic_error,"link in method pack is not 0 as expected");
                }

                size_t buffers = bufferList_.size();
                buffer.write( buffers );
                for (auto it = bufferList_.begin(), end = bufferList_.end(); it != end; ++it )
                {
                    write( buffer, localIndexMap_, *(it->second) );
                }
            }

            void doUnpack( const IndexMapType& indexMap, MessageBufferType& buffer )
            {
                size_t buffers = 0;
                buffer.read( buffers );
                assert( buffers == bufferList_.size() );
                for( auto it = bufferList_.begin(), end = bufferList_.end(); it != end; ++it )
                {
                    read( buffer, indexMap, *(it->second) );
                }
            }

            // unpack all data associated with link
            void unpack( const int link, MessageBufferType& buffer )
            {
                doUnpack( indexMaps_[ link ], buffer );
            }

        protected:
            template <class Vector>
            void write( MessageBufferType& buffer, const IndexMapType& localIndexMap, const Vector& data ) const
            {
                const size_t size = localIndexMap.size();
                assert( size <= data.size() );
                buffer.write( size );
                for( size_t i=0; i<size; ++i )
                {
                    buffer.write( data[ localIndexMap[ i ] ] );
                }
            }

            template <class Vector>
            void read( MessageBufferType& buffer, const IndexMapType& indexMap, Vector& data ) const
            {
                size_t size = indexMap.size();
                assert( size <= data.size() );
                buffer.read( size );
                assert( size == indexMap.size() );
                for( size_t i=0; i<size; ++i )
                {
                    buffer.read( data[ indexMap[ i ] ] );
                }
            }

            void writeString( MessageBufferType& buffer, const std::string& s) const
            {
                const int size = s.size();
                buffer.write( size );
                for( int i=0; i<size; ++i )
                {
                    buffer.write( s[ i ] );
                }
            }

            void readString( MessageBufferType& buffer, std::string& s) const
            {
                int size = -1;
                buffer.read( size );
                s.resize( size );
                for( int i=0; i<size; ++i )
                {
                    buffer.read( s[ i ] );
                }
            }
        };

        // gather solution to rank 0 for EclipseWriter
        template <class BufferList>
        void collect( BufferList& bufferList ) const
        {
            // index maps only have to be build when reordering is needed
            if( ! needsReordering && ! isParallel_ )
            {
                return ;
            }

            // this also packs and unpacks the local buffers one ioRank
            PackUnPackOutputBuffers< BufferList >
                packUnpack( bufferList,
                            localIndexMap_,
                            indexMaps_,
                            numCells(),
                            isIORank() );

            if ( ! isParallel_ )
            {
                // no need to collect anything.
                return;
            }

            //toIORankComm_.exchangeCached( packUnpack );
            toIORankComm_.exchange( packUnpack );

#ifndef NDEBUG
            // mkae sure every process is on the same page
            toIORankComm_.barrier();
#endif
        }

        bool isIORank() const
        {
            return isIORank_;
        }

        bool isParallel() const
        {
            return toIORankComm_.size() > 1;
        }

        size_t numCells () const { return globalCartesianIndex_.size(); }

    protected:
        P2PCommunicatorType             toIORankComm_;
        IndexMapType                    globalCartesianIndex_;
        IndexMapType                    localIndexMap_;
        IndexMapStorageType             indexMaps_;
        // true if we are on I/O rank
        bool                            isIORank_;
        /// \brief True if there is more than one MPI process
        bool                            isParallel_;
    };

} // end namespace Opm
#endif
