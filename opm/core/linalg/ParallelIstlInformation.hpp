/*
  Copyright 2014 Dr. Markus Blatt - HPC-Simulation-Software & Services

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
#ifndef OPM_PARALLELISTLINFORMTION_HEADER_INCLUDED
#define OPM_PARALLELISTLINFORMTION_HEADER_INCLUDED

#if HAVE_MPI && HAVE_DUNE_ISTL

#include "mpi.h"
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/common/parallel/interface.hh>
#include <dune/common/parallel/communicator.hh>
#include <dune/common/enumset.hh>

#include<algorithm>

namespace Opm
{

/// \brief Class that encapsulates the parallelization information needed by the
/// ISTL solvers.
class ParallelISTLInformation
{
public:
    /// \brief The type of the parallel index set used.
    typedef Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet ParallelIndexSet;
    /// \brief The type of the remote indices information used.
    typedef Dune::OwnerOverlapCopyCommunication<int, int>::RemoteIndices RemoteIndices;

    /// \brief Constructs an empty parallel information object using MPI_COMM_WORLD
    ParallelISTLInformation()
        : indexSet_(new ParallelIndexSet),
          remoteIndices_(new RemoteIndices(*indexSet_, *indexSet_, MPI_COMM_WORLD)),
          communicator_(MPI_COMM_WORLD)
    {}
    /// \brief Constructs an empty parallel information object using a communicator.
    /// \param communicator The communicator to use.
    ParallelISTLInformation(MPI_Comm communicator)
        : indexSet_(new ParallelIndexSet),
          remoteIndices_(new RemoteIndices(*indexSet_, *indexSet_, communicator)),
          communicator_(communicator)
    {}
    /// \brief Constructs a parallel information object from the specified information.
    /// \param indexSet The parallel index set to use.
    /// \param remoteIndices The remote indices information to use.
    /// \param communicator The communicator to use.
    ParallelISTLInformation(const std::shared_ptr<ParallelIndexSet>& indexSet,
                            const std::shared_ptr<RemoteIndices>& remoteIndices,
                            MPI_Comm communicator)
        : indexSet_(indexSet), remoteIndices_(remoteIndices), communicator_(communicator)
    {}
    /// \brief Copy constructor.
    ///
    /// The information will be shared by the the two objects.
    ParallelISTLInformation(const ParallelISTLInformation& other)
    : indexSet_(other.indexSet_), remoteIndices_(other.remoteIndices_),
      communicator_(other.communicator_)
    {}
    /// \brief Get a pointer to the underlying index set.
    std::shared_ptr<ParallelIndexSet> indexSet() const
    {
        return indexSet_;
    }
    /// \brief Get a pointer to the remote indices information.
    std::shared_ptr<RemoteIndices> remoteIndices() const
    {
        return remoteIndices_;
    }
    /// \brief Get the MPI communicator that we use.
    MPI_Comm communicator() const
    {
        return communicator_;
    }
    /// \brief Copy the information stored to the specified objects.
    /// \param[out] indexSet The object to store the index set in.
    /// \param[out] remoteIndices The object to store the remote indices information in.
    void copyValuesTo(ParallelIndexSet& indexSet, RemoteIndices& remoteIndices) const
    {
        indexSet.beginResize();
        IndexSetInserter<ParallelIndexSet> inserter(indexSet);
        std::for_each(indexSet_->begin(), indexSet_->end(), inserter);
        indexSet.endResize();
        remoteIndices.rebuild<false>();
    }
    /// \brief Communcate the dofs owned by us to the other process.
    ///
    /// Afterwards all associated dofs will contain the same data.
    template<class T>
    void copyOwnerToAll (const T& source, T& dest) const
    {
        typedef Dune::Combine<Dune::EnumItem<Dune::OwnerOverlapCopyAttributeSet::AttributeSet,Dune::OwnerOverlapCopyAttributeSet::owner>,Dune::EnumItem<Dune::OwnerOverlapCopyAttributeSet::AttributeSet,Dune::OwnerOverlapCopyAttributeSet::overlap>,Dune::OwnerOverlapCopyAttributeSet::AttributeSet> OwnerOverlapSet;
        typedef Dune::Combine<OwnerOverlapSet, Dune::EnumItem<Dune::OwnerOverlapCopyAttributeSet::AttributeSet,Dune::OwnerOverlapCopyAttributeSet::copy>,Dune::OwnerOverlapCopyAttributeSet::AttributeSet> AllSet;
      OwnerOverlapSet sourceFlags;
      AllSet destFlags;
      Dune::Interface interface(communicator_);
      if(!remoteIndices_->isSynced())
          remoteIndices_->rebuild<false>();
      interface.build(*remoteIndices_,sourceFlags,destFlags);
      Dune::BufferedCommunicator communicator;
      communicator.template build<T>(interface);
      communicator.template forward<CopyGatherScatter<T> >(source,dest);
      communicator.free();
    }    
private:
        /** \brief gather/scatter callback for communcation */
    template<typename T>
    struct CopyGatherScatter
    {
        typedef typename Dune::CommPolicy<T>::IndexedType V;

      static V gather(const T& a, std::size_t i)
      {
        return a[i];
      }

      static void scatter(T& a, V v, std::size_t i)
      {
        a[i] = v;
      }
    };
    template<class T>
    class IndexSetInserter
    {
    public:
        typedef T ParallelIndexSet;

        IndexSetInserter(ParallelIndexSet& indexSet)
        : indexSet_(&indexSet)
        {}
        void operator()(const typename ParallelIndexSet::IndexPair& pair)
        {
            indexSet_->add(pair.global(), pair.local());
        }

    private:
        ParallelIndexSet* indexSet_;
    };

    std::shared_ptr<ParallelIndexSet> indexSet_;
    std::shared_ptr<RemoteIndices> remoteIndices_;
    MPI_Comm communicator_;
};
} // end namespace Opm
#endif
#endif
