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
#ifndef OPM_PARALLELISTLINFORMATION_HEADER_INCLUDED
#define OPM_PARALLELISTLINFORMATION_HEADER_INCLUDED

#include <vector>

#if HAVE_MPI && HAVE_DUNE_ISTL

#include <algorithm>
#include <limits>
#include <memory>
#include <tuple>
#include <type_traits>

#include <dune/istl/owneroverlapcopy.hh>

#include <opm/simulators/utils/ParallelCommunication.hpp>

namespace Opm
{

/// \brief Class that encapsulates the parallelization information needed by the
/// ISTL solvers.
class ParallelISTLInformation
{
public:
    /// \brief The type of the parallel index set used.
    using ParallelIndexSet = Dune::OwnerOverlapCopyCommunication<int, int>::ParallelIndexSet;
    /// \brief The type of the remote indices information used.
    using RemoteIndices = Dune::OwnerOverlapCopyCommunication<int, int>::RemoteIndices;

    /// \brief Constructs an empty parallel information object using MPI_COMM_WORLD
    ParallelISTLInformation();

    /// \brief Constructs an empty parallel information object using a communicator.
    /// \param communicator The communicator to use.
    ParallelISTLInformation(MPI_Comm communicator);

    /// \brief Constructs a parallel information object from the specified information.
    /// \param indexSet The parallel index set to use.
    /// \param remoteIndices The remote indices information to use.
    /// \param communicator The communicator to use.
    ParallelISTLInformation(const std::shared_ptr<ParallelIndexSet>& indexSet,
                            const std::shared_ptr<RemoteIndices>& remoteIndices,
                            MPI_Comm communicator);

    /// \brief Copy constructor.
    ///
    /// The information will be shared by the the two objects.
    ParallelISTLInformation(const ParallelISTLInformation& other);

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
    Parallel::Communication communicator() const
    {
        return communicator_;
    }
    /// \brief Copy the information stored to the specified objects.
    /// \param[out] indexSet The object to store the index set in.
    /// \param[out] remoteIndices The object to store the remote indices information in.
    void copyValuesTo(ParallelIndexSet& indexSet, RemoteIndices& remoteIndices,
                      std::size_t local_component_size = 0,
                      std::size_t num_components = 1) const;

    /// \brief Communcate the dofs owned by us to the other process.
    ///
    /// Afterwards all associated dofs will contain the same data.
    template<class T>
    void copyOwnerToAll (const T& source, T& dest) const;

    template<class T>
    const std::vector<double>& updateOwnerMask(const T& container) const;

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
                          T& value) const;

private:
    template<typename... Containers, typename... BinaryOperators, typename... ReturnValues>
    void computeTupleReduction(const std::tuple<Containers...>& containers,
                               std::tuple<BinaryOperators...>& operators,
                               std::tuple<ReturnValues...>& values) const;

    std::shared_ptr<ParallelIndexSet> indexSet_;
    std::shared_ptr<RemoteIndices> remoteIndices_;
    Parallel::Communication communicator_;
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
    auto makeGlobalMaxFunctor()
    {
        struct MaxOp
        {
            using result_type = T;
            const result_type& operator()(const T& t1, const T& t2)
            {
                return std::max(t1, t2);
            }
        };
        return MaskToMinOperator(MaxOp());
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
    auto
    makeGlobalMinFunctor()
    {
        struct MinOp
        {
            using result_type = T;
            const result_type& operator()(const T& t1, const T& t2)
            {
                return std::min(t1, t2);
            }
        };
        return MaskToMaxOperator(MinOp());
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
/// \brief Accumulates entries masked with 1.
/// \param container The container whose values to accumulate.
/// \param maskContainer null pointer or a pointer to a container
/// with entries 0 and 1. Only values at indices with a 1 stored
/// will be accumulated. If null then all values will be accumulated
/// \return the summ of all entries that should be represented.
template<class T1>
auto
accumulateMaskedValues(const T1& container, const std::vector<double>* maskContainer)
    -> decltype(container[0]*(*maskContainer)[0]);

} // end namespace Opm

#endif
