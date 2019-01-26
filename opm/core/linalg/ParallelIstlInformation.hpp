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
#ifndef OPM_PARALLELISTLINFORMTION_HEADER_INCLUDED
#define OPM_PARALLELISTLINFORMTION_HEADER_INCLUDED

#include <opm/grid/UnstructuredGrid.h>
#include <opm/common/ErrorMacros.hpp>
#include <boost/any.hpp>
#include <exception>

#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>
#include <type_traits>
#include <vector>

#if HAVE_MPI && HAVE_DUNE_ISTL

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <mpi.h>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/common/parallel/interface.hh>
#include <dune/common/parallel/communicator.hh>
#include <dune/common/enumset.hh>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

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
    /// \brief Get the Collective MPI communicator that we use.
    Dune::CollectiveCommunication<MPI_Comm> communicator() const
    {
        return communicator_;
    }
    /// \brief Copy the information stored to the specified objects.
    /// \param[out] indexSet The object to store the index set in.
    /// \param[out] remoteIndices The object to store the remote indices information in.
    void copyValuesTo(ParallelIndexSet& indexSet, RemoteIndices& remoteIndices,
                      std::size_t local_component_size = 0, std::size_t num_components = 1) const
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
    /// \brief Communcate the dofs owned by us to the other process.
    ///
    /// Afterwards all associated dofs will contain the same data.
    template<class T>
    void copyOwnerToAll (const T& source, T& dest) const
    {
        typedef Dune::Combine<Dune::EnumItem<Dune::OwnerOverlapCopyAttributeSet::AttributeSet,Dune::OwnerOverlapCopyAttributeSet::owner>,Dune::EnumItem<Dune::OwnerOverlapCopyAttributeSet::AttributeSet,Dune::OwnerOverlapCopyAttributeSet::overlap>,Dune::OwnerOverlapCopyAttributeSet::AttributeSet> OwnerOverlapSet;
        typedef Dune::EnumItem<Dune::OwnerOverlapCopyAttributeSet::AttributeSet,Dune::OwnerOverlapCopyAttributeSet::owner> OwnerSet;
        typedef Dune::Combine<OwnerOverlapSet, Dune::EnumItem<Dune::OwnerOverlapCopyAttributeSet::AttributeSet,Dune::OwnerOverlapCopyAttributeSet::copy>,Dune::OwnerOverlapCopyAttributeSet::AttributeSet> AllSet;
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
      communicator.template forward<CopyGatherScatter<T> >(source,dest);
      communicator.free();
    }
    template<class T>
    const std::vector<double>& updateOwnerMask(const T& container) const
    {
        if( ! indexSet_ )
        {
            OPM_THROW(std::runtime_error, "Trying to update owner mask without parallel information!");
        }
        if( static_cast<std::size_t>(container.size())!= ownerMask_.size() )
        {
            ownerMask_.resize(container.size(), 1.);
            for( auto i=indexSet_->begin(), end=indexSet_->end(); i!=end; ++i )
            {
                if (i->local().attribute()!=Dune::OwnerOverlapCopyAttributeSet::owner)
                {
                    ownerMask_[i->local().local()] = 0.;
                }
            }
        }
        return ownerMask_;
    }

    /// \brief Get the owner Mask.
    ///
    /// \return A vector with entries 0, and 1. 0 marks an index that we cannot
    ///         compute correct results for. 1 marks an index that this process
    ///         is responsible for and computes correct results in parallel.
    const std::vector<double>& getOwnerMask() const
    {
        return ownerMask_;
    }

    /// \brief Compute one or more global reductions.
    ///
    /// This function can either be used with a container, an operator, and an initial value
    /// to compute a reduction. Or with tuples of them to compute multiple reductions with only
    /// one global communication.
    /// The possible functors needed can be constructed with Opm::Reduction::makeGlobalMaxFunctor(),
    /// Opm::Reduction::makeLInfinityNormFunctor(),
    /// Opm::Reduction::makeGlobalMinFunctor(), and 
    /// Opm::Reduction::makeGlobalSumFunctor().
    /// \tparam type of the container or the tuple of  containers.
    /// \tparam tyoe of the operator or a tuple of operators, examples are e.g. 
    /// Reduction::MaskIDOperator, Reduction::MaskToMinOperator,
    /// and Reduction::MaskToMaxOperator. Has to provide an operator() that takes three
    /// arguments (the last one is the mask value: 1 for a dof that we own, 0 otherwise),
    /// a method maskValue that takes a value and mask value, and localOperator that
    /// returns the underlying binary operator.
    /// \param container A container or tuple of containers.
    /// \param binaryOperator An operator doing the reduction of two values.
    /// \param value The initial value or a tuple of them.
    template<typename Container, typename BinaryOperator, typename T>
    void computeReduction(const Container& container, BinaryOperator binaryOperator,
                          T& value) const
    {
        computeReduction(container, binaryOperator, value, is_tuple<Container>());
    }
private:
    /// \brief compute the reductions for tuples.
    ///
    /// This is a helper function to prepare for calling computeTupleReduction.
    template<typename Container, typename BinaryOperator, typename T>
    void computeReduction(const Container& container, BinaryOperator binaryOperator,
                          T& value, std::integral_constant<bool,true>) const
    {
        computeTupleReduction(container, binaryOperator, value);
    }
    /// \brief compute the reductions for non-tuples.
    ///
    /// This is a helper function to prepare for calling computeTupleReduction.
    template<typename Container, typename BinaryOperator, typename T>
    void computeReduction(const Container& container, BinaryOperator binaryOperator,
                          T& value, std::integral_constant<bool,false>) const
    {
        std::tuple<const Container&> containers=std::tuple<const Container&>(container);
        auto values=std::make_tuple(value);
        auto operators=std::make_tuple(binaryOperator);
        computeTupleReduction(containers, operators, values);
        value=std::get<0>(values);
    }
    /// \brief Compute the reductions for tuples.
    template<typename... Containers, typename... BinaryOperators, typename... ReturnValues>
    void computeTupleReduction(const std::tuple<Containers...>& containers,
                               std::tuple<BinaryOperators...>& operators,
                               std::tuple<ReturnValues...>& values) const
    {
        static_assert(std::tuple_size<std::tuple<Containers...> >::value==
                      std::tuple_size<std::tuple<BinaryOperators...> >::value,
                      "We need the same number of containers and binary operators");
        static_assert(std::tuple_size<std::tuple<Containers...> >::value==
                      std::tuple_size<std::tuple<ReturnValues...> >::value,
                      "We need the same number of containers and return values");
        if( std::tuple_size<std::tuple<Containers...> >::value==0 )
        {
            return;
        }
        // Copy the initial values.
        std::tuple<ReturnValues...> init=values;
        updateOwnerMask(std::get<0>(containers));
        computeLocalReduction(containers, operators, values);
        std::vector<std::tuple<ReturnValues...> > receivedValues(communicator_.size());
        communicator_.allgather(&values, 1, &(receivedValues[0]));
        values=init;
        for( auto rvals=receivedValues.begin(), endvals=receivedValues.end(); rvals!=endvals;
             ++rvals )
        {
            computeGlobalReduction(*rvals, operators, values);
        }
    }
    /// \brief TMP for computing the the global reduction after receiving the local ones.
    ///
    /// End of recursion.
    template<int I=0, typename... BinaryOperators, typename... ReturnValues>
    typename std::enable_if<I == sizeof...(BinaryOperators), void>::type
    computeGlobalReduction(const std::tuple<ReturnValues...>&,
                                std::tuple<BinaryOperators...>&,
                                std::tuple<ReturnValues...>&) const
    {}
    /// \brief TMP for computing the the global reduction after receiving the local ones.
    template<int I=0, typename... BinaryOperators, typename... ReturnValues>
    typename std::enable_if<I !=sizeof...(BinaryOperators), void>::type
    computeGlobalReduction(const std::tuple<ReturnValues...>& receivedValues,
                           std::tuple<BinaryOperators...>& operators,
                           std::tuple<ReturnValues...>& values) const
    {
        auto& val=std::get<I>(values);
        val = std::get<I>(operators).localOperator()(val, std::get<I>(receivedValues));
        computeGlobalReduction<I+1>(receivedValues, operators, values);
    }
    /// \brief TMP for computing the the local reduction on the DOF that the process owns.
    ///
    /// End of recursion.
    template<int I=0, typename... Containers, typename... BinaryOperators, typename... ReturnValues>
    typename std::enable_if<I==sizeof...(Containers), void>::type
    computeLocalReduction(const std::tuple<Containers...>&,
                          std::tuple<BinaryOperators...>&,
                          std::tuple<ReturnValues...>&) const
    {}
    /// \brief TMP for computing the the local reduction on the DOF that the process owns.
    template<int I=0, typename... Containers, typename... BinaryOperators, typename... ReturnValues>
    typename std::enable_if<I!=sizeof...(Containers), void>::type
    computeLocalReduction(const std::tuple<Containers...>& containers,
                          std::tuple<BinaryOperators...>& operators,
                          std::tuple<ReturnValues...>& values) const
    {
        const auto& container = std::get<I>(containers);
        if( container.size() )
        {
            auto& reduceOperator  = std::get<I>(operators);
            // Eigen:Block does not support STL iterators!!!!
            // Therefore we need to rely on the harder random-access
            // property of the containers. But this should be save, too.
            // Just commenting out code in the hope that Eigen might improve
            // in this regard in the future.
            //auto newVal = container.begin();
            auto mask   = ownerMask_.begin();
            auto& value = std::get<I>(values);
            value =  reduceOperator.getInitialValue();

            for( auto endVal=ownerMask_.end(); mask!=endVal;
                 /*++newVal,*/ ++mask )
            {
                value = reduceOperator(value, container[mask-ownerMask_.begin()], *mask);
            }
        }
        computeLocalReduction<I+1>(containers, operators, values);
    }
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
        typedef typename ParallelIndexSet::LocalIndex LocalIndex;
        typedef typename ParallelIndexSet::GlobalIndex GlobalIndex;

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
    std::shared_ptr<ParallelIndexSet> indexSet_;
    std::shared_ptr<RemoteIndices> remoteIndices_;
    Dune::CollectiveCommunication<MPI_Comm> communicator_;
    mutable std::vector<double> ownerMask_;
};

    namespace Reduction
    {
    /// \brief An operator that only uses values where mask is 1.
    ///
    /// Could be used to compute a global sum
    /// \tparam BinaryOperator The wrapped binary operator that specifies
    // the reduction operation.
    template<typename BinaryOperator>
    struct MaskIDOperator
    {
        // This is a real nice one: numeric limits needs a type without const
        // or reference qualifier. Otherwise we get complete nonesense.
        typedef typename std::remove_cv<
            typename std::remove_reference<typename BinaryOperator::result_type>::type
            >::type Result;
        /// \brief Apply the underlying binary operator according to the mask.
        ///
        /// The BinaryOperator will be called with t1, and mask*t2.
        /// \param t1 first value
        /// \param t2 second value (might be modified).
        /// \param mask The mask (0 or 1).
        template<class T, class T1>
        T operator()(const T& t1, const T& t2, const T1& mask)
        {
            return b_(t1, maskValue(t2, mask));
        }
        template<class T, class T1>
        T maskValue(const T& t, const T1& mask)
        {
            return t*mask;
        }
        BinaryOperator& localOperator()
        {
            return b_;
        }
        Result getInitialValue()
        {
            return Result();
        }
    private:
        BinaryOperator b_;
    };

    /// \brief An operator for computing a parallel inner product.
    template<class T>
    struct InnerProductFunctor
    {
        /// \brief Apply the underlying binary operator according to the mask.
        ///
        /// The BinaryOperator will be called with t1, and mask*t2.
        /// \param t1 first value
        /// \param t2 second value (might be modified).
        /// \param mask The mask (0 or 1).
        template<class T1>
        T operator()(const T& t1, const T& t2, const T1& mask)
        {
            T masked =  maskValue(t2, mask);
            return t1 + masked * masked;
        }
        template<class T1>
        T maskValue(const T& t, const T1& mask)
        {
            return t*mask;
        }
        std::plus<T> localOperator()
        {
            return std::plus<T>();
        }
        T getInitialValue()
        {
            return T();
        }
    };

    /// \brief An operator that converts the values where mask is 0 to the minimum value
    ///
    /// Could be used to compute a global maximum.
    /// \tparam BinaryOperator The wrapped binary operator that specifies
    // the reduction operation.
    template<typename BinaryOperator>
    struct MaskToMinOperator
    {
        // This is a real nice one: numeric limits has to a type without const
        // or reference. Otherwise we get complete nonesense.
        typedef typename std::remove_reference<
            typename std::remove_const<typename BinaryOperator::result_type>::type
            >::type Result;

        MaskToMinOperator(BinaryOperator b)
        : b_(b)
        {}
        /// \brief Apply the underlying binary operator according to the mask.
        ///
        /// If mask is 0 then t2 will be substituted by the lowest value,
        /// else t2 will be used.
        /// \param t1 first value
        /// \param t2 second value (might be modified).
        template<class T, class T1>
        T operator()(const T& t1, const T& t2, const T1& mask)
        {
            return b_(t1, maskValue(t2, mask));
        }
        template<class T, class T1>
        T maskValue(const T& t, const T1& mask)
        {
            if( mask )
            {
                return t;
            }
            else
            {
                return getInitialValue();
            }
        }
        Result getInitialValue()
        {
            //g++-4.4 does not support std::numeric_limits<T>::lowest();
            // we rely on IEE 754 for floating point values and use min()
            // for integral types.
            if( std::is_integral<Result>::value )
            {
                return std::numeric_limits<Result>::min();
            }
            else
            {
                return -std::numeric_limits<Result>::max();
            }
        }
        /// \brief Get the underlying binary operator.
        ///
        /// This might be needed to compute the reduction after each processor
        /// has computed its local one.
        BinaryOperator& localOperator()
        {
            return b_;
        }
    private:
        BinaryOperator b_;
    };

    /// \brief An operator that converts the values where mask is 0 to the maximum value
    ///
    /// Could be used to compute a global minimum.
    template<typename BinaryOperator>
    struct MaskToMaxOperator
    {
        // This is a real nice one: numeric limits has to a type without const
        // or reference. Otherwise we get complete nonesense.
        typedef typename std::remove_cv<
            typename std::remove_reference<typename BinaryOperator::result_type>::type
            >::type Result;

        MaskToMaxOperator(BinaryOperator b)
        : b_(b)
        {}
        /// \brief Apply the underlying binary operator according to the mask.
        ///
        /// If mask is 0 then t2 will be substituted by the maximum value,
        /// else t2 will be used.
        /// \param t1 first value
        /// \param t2 second value (might be modified).
        template<class T, class T1>
        T operator()(const T& t1, const T& t2, const T1& mask)
        {
            return b_(t1, maskValue(t2, mask));
        }
        template<class T, class T1>
        T maskValue(const T& t, const T1& mask)
        {
            if( mask )
            {
                return t;
            }
            else
            {
                return std::numeric_limits<T>::max();
            }
        }
        BinaryOperator& localOperator()
        {
            return b_;
        }
        Result getInitialValue()
        {
            return std::numeric_limits<Result>::max();
        }
    private:
        BinaryOperator b_;
    };
    /// \brief Create a functor for computing a global sum.
    ///
    /// To be used with ParallelISTLInformation::computeReduction.
    template<class T>
    MaskIDOperator<std::plus<T> >
    makeGlobalSumFunctor()
    {
        return MaskIDOperator<std::plus<T> >();
    }
    /// \brief Create a functor for computing a global maximum.
    ///
    /// To be used with ParallelISTLInformation::computeReduction.
    template<class T>
    MaskToMinOperator<std::pointer_to_binary_function<const T&,const T&,const T&> >
    makeGlobalMaxFunctor()
    {
        return MaskToMinOperator<std::pointer_to_binary_function<const T&,const T&,const T&> >
            (std::pointer_to_binary_function<const T&,const T&,const T&>
             ((const T&(*)(const T&, const T&))std::max<T>));
    }

    namespace detail
    {
        /// \brief Computes the maximum of the absolute values of two values.
        template<typename T, typename Enable = void>
        struct MaxAbsFunctor
        {
            using result_type = T;
            result_type operator()(const T& t1,
                                   const T& t2)
            {
                return std::max(std::abs(t1), std::abs(t2));
            }
        };

        // Specialization for unsigned integers. They need their own
        // version since abs(x) is ambiguous (as well as somewhat
        // meaningless).
        template<typename T>
        struct MaxAbsFunctor<T, typename std::enable_if<std::is_unsigned<T>::value>::type>
        {
            using result_type = T;
            result_type operator()(const T& t1,
                                   const T& t2)
            {
                return std::max(t1, t2);
            }
        };
    }

    /// \brief Create a functor for computing a global L infinity norm
    ///
    /// To be used with ParallelISTLInformation::computeReduction.
    template<class T>
    MaskIDOperator<detail::MaxAbsFunctor<T> >
    makeLInfinityNormFunctor()
    {
        return MaskIDOperator<detail::MaxAbsFunctor<T> >();
    }
    /// \brief Create a functor for computing a global minimum.
    ///
    /// To be used with ParallelISTLInformation::computeReduction.
    template<class T>
    MaskToMaxOperator<std::pointer_to_binary_function<const T&,const T&,const T&> >
    makeGlobalMinFunctor()
    {
        return MaskToMaxOperator<std::pointer_to_binary_function<const T&,const T&,const T&> >
            (std::pointer_to_binary_function<const T&,const T&,const T&>
             ((const T&(*)(const T&, const T&))std::min<T>));
    }
    template<class T>
    InnerProductFunctor<T>
    makeInnerProductFunctor()
    {
        return InnerProductFunctor<T>();
    }
    } // end namespace Reduction
} // end namespace Opm

#endif

namespace Opm
{
/// \brief Extracts the information about the data decomposition from the grid for dune-istl
/// 
/// In the case that grid is a parallel grid this method will query it to get the information
/// about the data decompoisition and convert it to the format expected by the linear algebra
/// of dune-istl.
/// \warn for UnstructuredGrid this function doesn't do anything.
/// \param anyComm The handle to store the information in. If grid is a parallel grid
/// then this will ecapsulate an instance of ParallelISTLInformation.
/// \param grid The grid to inspect.

inline void extractParallelGridInformationToISTL(boost::any& anyComm, const UnstructuredGrid& grid)
{
    (void)anyComm; (void)grid;
}

/// \brief Accumulates entries masked with 1.
/// \param container The container whose values to accumulate.
/// \param maskContainer null pointer or a pointer to a container
/// with entries 0 and 1. Only values at indices with a 1 stored
/// will be accumulated. If null then all values will be accumulated
/// \return the summ of all entries that should be represented.
template<class T1>
auto
accumulateMaskedValues(const T1& container, const std::vector<double>* maskContainer)
    -> decltype(container[0]*(*maskContainer)[0])
{
    decltype(container[0]*(*maskContainer)[0]) initial = 0;

    if( maskContainer )
    {
        return std::inner_product(container.begin(), container.end(), maskContainer->begin(),
                                  initial);
    }else
    {
        return std::accumulate(container.begin(), container.end(), initial);
    }
}   
} // end namespace Opm

#endif
