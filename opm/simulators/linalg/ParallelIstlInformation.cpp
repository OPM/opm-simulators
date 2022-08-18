/*
  Copyright 2014, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA
  Copyright 2015 NTNU

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

#include <config.h>

#if HAVE_MPI && HAVE_DUNE_ISTL

#include <opm/simulators/linalg/ParallelIstlInformation.hpp>

#include <dune/common/enumset.hh>

#include <opm/common/ErrorMacros.hpp>

#include <cstddef>
#include <exception>
#include <mpi.h>
#include <numeric>

namespace
{

template<class T>
class IndexSetInserter
{
public:
    using ParallelIndexSet = T;
    using LocalIndex = typename ParallelIndexSet::LocalIndex;
    using GlobalIndex = typename ParallelIndexSet::GlobalIndex;

    IndexSetInserter(ParallelIndexSet& indexSet, const GlobalIndex& component_size,
                     std::size_t local_component_size, std::size_t num_components)
        : indexSet_(&indexSet), component_size_(component_size),
          local_component_size_(local_component_size),
          num_components_(num_components)
    {}

    void operator()(const typename ParallelIndexSet::IndexPair& pair)
    {
        for(std::size_t i = 0; i < num_components_; i++)
            indexSet_->add(i * component_size_ + pair.global(),
                           LocalIndex(i * local_component_size_  + pair.local(),
                                      pair.local().attribute()));
    }
private:
    ParallelIndexSet* indexSet_;
    /// \brief The global number of unknowns per component/equation.
    GlobalIndex component_size_;
    /// \brief The local number of unknowns per component/equation.
    std::size_t local_component_size_;
    /// \brief The number of components/equations.
    std::size_t num_components_;
};

/** \brief gather/scatter callback for communcation */
template<typename T>
struct CopyGatherScatter
{
    using V = typename Dune::CommPolicy<T>::IndexedType;

    static V gather(const T& a, std::size_t i)
    {
      return a[i];
    }

    static void scatter(T& a, V v, std::size_t i)
    {
      a[i] = v;
    }
};

template<int I=0, typename... BinaryOperators, typename... ReturnValues>
typename std::enable_if<I == sizeof...(BinaryOperators), void>::type
computeGlobalReduction(const std::tuple<ReturnValues...>&,
                       std::tuple<BinaryOperators...>&,
                       std::tuple<ReturnValues...>&)
{}

template<int I=0, typename... BinaryOperators, typename... ReturnValues>
typename std::enable_if<I !=sizeof...(BinaryOperators), void>::type
computeGlobalReduction(const std::tuple<ReturnValues...>& receivedValues,
                       std::tuple<BinaryOperators...>& operators,
                       std::tuple<ReturnValues...>& values)
{
    auto& val = std::get<I>(values);
    val = std::get<I>(operators).localOperator()(val, std::get<I>(receivedValues));
    computeGlobalReduction<I+1>(receivedValues, operators, values);
}

template<int I=0, typename... Containers, typename... BinaryOperators, typename... ReturnValues>
typename std::enable_if<I==sizeof...(Containers), void>::type
computeLocalReduction(const std::tuple<Containers...>&,
                      std::tuple<BinaryOperators...>&,
                      std::tuple<ReturnValues...>&,
                      const std::vector<double>&)
{}

template<int I=0, typename... Containers, typename... BinaryOperators, typename... ReturnValues>
typename std::enable_if<I!=sizeof...(Containers), void>::type
computeLocalReduction(const std::tuple<Containers...>& containers,
                      std::tuple<BinaryOperators...>& operators,
                      std::tuple<ReturnValues...>& values,
                      const std::vector<double>& ownerMask)
{
    const auto& container = std::get<I>(containers);
    if (container.size())
    {
        auto& reduceOperator = std::get<I>(operators);
        // Eigen:Block does not support STL iterators!!!!
        // Therefore we need to rely on the harder random-access
        // property of the containers. But this should be save, too.
        // Just commenting out code in the hope that Eigen might improve
        // in this regard in the future.
        //auto newVal = container.begin();
        auto mask   = ownerMask.begin();
        auto& value = std::get<I>(values);
        value = reduceOperator.getInitialValue();

        for (auto endVal = ownerMask.end(); mask != endVal; /*++newVal,*/ ++mask )
        {
            value = reduceOperator(value, container[mask-ownerMask.begin()], *mask);
        }
    }
    computeLocalReduction<I+1>(containers, operators, values, ownerMask);
}

}

namespace Opm
{
namespace
{

    template<class T>
    struct is_tuple
        : std::integral_constant<bool, false>
    {};
    template<typename... T>
    struct is_tuple<std::tuple<T...> >
        : std::integral_constant<bool, true>
    {};
}

ParallelISTLInformation::ParallelISTLInformation()
    : indexSet_(new ParallelIndexSet),
      remoteIndices_(new RemoteIndices(*indexSet_, *indexSet_, MPI_COMM_WORLD)),
      communicator_(MPI_COMM_WORLD)
{}


ParallelISTLInformation::ParallelISTLInformation(MPI_Comm communicator)
    : indexSet_(new ParallelIndexSet),
      remoteIndices_(new RemoteIndices(*indexSet_, *indexSet_, communicator)),
      communicator_(communicator)
{}


ParallelISTLInformation::
ParallelISTLInformation(const std::shared_ptr<ParallelIndexSet>& indexSet,
                        const std::shared_ptr<RemoteIndices>& remoteIndices,
                        MPI_Comm communicator)
    : indexSet_(indexSet), remoteIndices_(remoteIndices), communicator_(communicator)
{}

ParallelISTLInformation::ParallelISTLInformation(const ParallelISTLInformation& other)
    : indexSet_(other.indexSet_), remoteIndices_(other.remoteIndices_),
      communicator_(other.communicator_)
{}

void ParallelISTLInformation::copyValuesTo(ParallelIndexSet& indexSet,
                                           RemoteIndices& remoteIndices,
                                           std::size_t local_component_size,
                                           std::size_t num_components) const
{
    ParallelIndexSet::GlobalIndex global_component_size  = local_component_size;
    if ( num_components > 1 )
    {
        ParallelIndexSet::GlobalIndex max_gi = 0;
        // component the max global index
        for( auto i = indexSet_->begin(), end = indexSet_->end(); i != end; ++i )
        {
            max_gi = std::max(max_gi, i->global());
        }
        global_component_size = max_gi+1;
        global_component_size = communicator_.max(global_component_size);
    }
    indexSet.beginResize();
    IndexSetInserter<ParallelIndexSet> inserter(indexSet, global_component_size,
                                                local_component_size, num_components);
    std::for_each(indexSet_->begin(), indexSet_->end(), inserter);
    indexSet.endResize();
    remoteIndices.rebuild<false>();
}

template<class T>
void ParallelISTLInformation::copyOwnerToAll(const T& source, T& dest) const
{
    using AS = Dune::OwnerOverlapCopyAttributeSet;
    using CopySet = Dune::EnumItem<AS, AS::copy>;
    using OwnerSet = Dune::EnumItem<AS, AS::owner>;
    using OverlapSet = Dune::EnumItem<AS, AS::overlap>;
    using OwnerOverlapSet = Dune::Combine<OwnerSet, OverlapSet, AS::AttributeSet>;
    using AllSet = Dune::Combine<OwnerOverlapSet, CopySet, AS::AttributeSet>;
    OwnerSet sourceFlags;
    AllSet destFlags;
    Dune::Interface interface(communicator_);
    if( !remoteIndices_->isSynced() )
    {
        remoteIndices_->rebuild<false>();
    }
    interface.build(*remoteIndices_,sourceFlags,destFlags);
    Dune::BufferedCommunicator communicator;
    communicator.template build<T>(interface);
    communicator.template forward<CopyGatherScatter<T>>(source,dest);
    communicator.free();
}

template<class T>
const std::vector<double>&
ParallelISTLInformation::updateOwnerMask(const T& container) const
{
    if (!indexSet_)
    {
        OPM_THROW(std::runtime_error, "Trying to update owner mask without parallel information!");
    }
    if (static_cast<std::size_t>(container.size()) != ownerMask_.size())
    {
        ownerMask_.resize(container.size(), 1.);
        for (const auto& i : *indexSet_)
        {
            if (i.local().attribute() != Dune::OwnerOverlapCopyAttributeSet::owner)
            {
                ownerMask_[i.local().local()] = 0.;
            }
        }
    }
    return ownerMask_;
}

template<typename Container, typename BinaryOperator, typename T>
void ParallelISTLInformation::computeReduction(const Container& container,
                                               BinaryOperator binaryOperator,
                                               T& value) const
{
    if constexpr (is_tuple<Container>())
        computeTupleReduction(container, binaryOperator, value);
    else
    {
        std::tuple<const Container&> containers = std::tuple<const Container&>(container);
        auto values = std::make_tuple(value);
        auto operators = std::make_tuple(binaryOperator);
        computeTupleReduction(containers, operators, values);
        value = std::get<0>(values);
    }
}

template<typename... Containers, typename... BinaryOperators, typename... ReturnValues>
void ParallelISTLInformation::computeTupleReduction(const std::tuple<Containers...>& containers,
                                                    std::tuple<BinaryOperators...>& operators,
                                                    std::tuple<ReturnValues...>& values) const
{
    static_assert(std::tuple_size<std::tuple<Containers...> >::value ==
                  std::tuple_size<std::tuple<BinaryOperators...> >::value,
                  "We need the same number of containers and binary operators");
    static_assert(std::tuple_size<std::tuple<Containers...> >::value ==
                  std::tuple_size<std::tuple<ReturnValues...> >::value,
                  "We need the same number of containers and return values");
    if (std::tuple_size<std::tuple<Containers...> >::value == 0)
    {
        return;
    }

    // Copy the initial values.
    std::tuple<ReturnValues...> init = values;
    updateOwnerMask(std::get<0>(containers));
    computeLocalReduction(containers, operators, values, ownerMask_);
    std::vector<std::tuple<ReturnValues...> > receivedValues(communicator_.size());
    communicator_.allgather(&values, 1, &(receivedValues[0]));
    values = init;
    for (auto& rval : receivedValues)
    {
        computeGlobalReduction(rval, operators, values);
    }
}

template<class T1>
auto
accumulateMaskedValues(const T1& container, const std::vector<double>* maskContainer)
    -> decltype(container[0]*(*maskContainer)[0])
{
    decltype(container[0]*(*maskContainer)[0]) initial = 0;

    if (maskContainer)
    {
        return std::inner_product(container.begin(), container.end(),
                                  maskContainer->begin(), initial);
    }
    else
    {
        return std::accumulate(container.begin(), container.end(), initial);
    }
}

template<class T> using C1 = std::vector<T>;
template<class T> using Ops1 = Reduction::MaskIDOperator<std::plus<T>>;

template<class T>
using C2 = std::tuple<std::vector<T>,
                      std::vector<T>,
                      std::vector<T>,
                      std::vector<T>,
                      std::vector<T>>;
template<class T>
using Ops2 = std::tuple<decltype(Reduction::makeGlobalSumFunctor<T>()),
                        decltype(Reduction::makeGlobalMaxFunctor<T>()),
                        decltype(Reduction::makeGlobalMinFunctor<T>()),
                        decltype(Reduction::makeInnerProductFunctor<T>()),
                        decltype(Reduction::makeLInfinityNormFunctor<T>())>;
template<class T>
using Vals2 = std::tuple<T,T,T,T,T>;

#define INSTANCE1(T) \
    template void ParallelISTLInformation::computeReduction<C1<T>,Ops1<T>,T>(const C1<T>&,Ops1<T>,T&) const;

#define INSTANCE2(T) \
    template void ParallelISTLInformation::computeReduction<C2<T>,Ops2<T>,Vals2<T>>(const C2<T>&,Ops2<T>,Vals2<T>&) const;

#define INSTANCE(T) \
    INSTANCE1(T) \
    INSTANCE2(T)

INSTANCE(int)
INSTANCE(float)
INSTANCE(std::size_t)

} // namespace Opm

#endif // HAVE_MPI && HAVE_DUNE_ISTL
