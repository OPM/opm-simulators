/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_RESERVOIR_COUPLING_MPI_TRAITS_HPP
#define OPM_RESERVOIR_COUPLING_MPI_TRAITS_HPP

#include <dune/common/parallel/mpitraits.hh>

#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>

#include <array>
#include <mutex>  // For std::call_once
#include <tuple>   // for std::tuple_size

#include <array>
#include <mutex>  // For std::call_once
#include <mpi.h>

namespace Dune {
namespace detail {

/// @brief Generic MPI traits implementation for structs
///
/// This class template provides automatic MPI datatype creation for C++ structs
/// by analyzing the struct's member pointers. It generates an MPI_Datatype that
/// correctly represents the memory layout of the struct, allowing it to be sent
/// via MPI communication functions.
///
/// The implementation handles:
/// - Scalar fields of MPI-compatible types
/// - C-style arrays (T[N])
/// - std::array<T, N> containers
/// - Enum fields (using their underlying type)
///
/// **Thread Safety:**
/// The MPI datatype is created exactly once using std::call_once() to ensure
/// thread-safe initialization in case this code is called from multiple threads.
///
/// **Usage Example:**
/// ```cpp
/// struct MyData {
///     double value;
///     std::array<int, 3> counts;
/// };
///
/// // Define the trait by listing member pointers:
/// template<>
/// struct MPITraits<MyData>
///     : detail::StructMPITraits<MyData, &MyData::value, &MyData::counts>
/// { };
///
/// // Now MyData can be sent via MPI:
/// MPI_Send(&data, 1, MPITraits<MyData>::getType(), dest, tag, comm);
/// ```
///
/// @tparam Struct The struct type for which to create an MPI datatype
/// @tparam Members Variadic pack of pointers-to-members listing all fields to include
///
/// **Assumptions:**
/// Each field of the struct must be either:
/// - An already defined MPI type (int, double, float, etc.)
/// - A C array or std::array of an already defined MPI type
/// - An enum whose underlying type is an already defined MPI type
///
/// @note The order of member pointers in the template argument list determines
///       the field order in the MPI datatype
/// @note The datatype is automatically resized to account for struct padding
/// @see MPITraits
template<class Struct, auto... Members>
struct StructMPITraitsImpl
{
    static MPI_Datatype getType()
    {
        // One-time, thread-safe construction: As a precaution, in case this code will be
        //   run by threads in the future, we want to prevent threads entering the code below
        //   simultaneously to build the datatype. This is done by using std::call_once() with
        //   a static flag.
        std::call_once(flag_, [] {

            constexpr std::size_t N = sizeof...(Members);

            // Array of block lengths (number of elements per field)
            //  NOTE: The block length is set to 1 for each field (blk.fill(1)),
            //        but is strictly not necessary currently, as the processMember()
            //        function will always set the block length anyway.
            //        However, it is kept here for clarity and future-proofing.
            std::array<int,      N> blk{};           blk.fill(1);
            // Array of displacements for each field
            std::array<MPI_Aint, N> disp{};
            // Array of MPI data types for each field
            std::array<MPI_Datatype, N> types{};

            // Create a dummy instance of the struct to get the base address
            Struct dummy{};
            // Get the base address of the dummy instance
            MPI_Aint base{};
            MPI_Get_address(&dummy, &base);

            std::size_t i = 0;
            // fold expression over the member pointer pack
            ( processMember<Members>(dummy, base, disp, blk, types, i++), ... );

            MPI_Datatype tmp;
            // Create the MPI datatype
            // The below assertion ensures that the displacements are in declaration order.
            // It is strictly not necessary, since MPI_Type_create_struct does not require
            // the displacements to be in declaration order, but it is a good practice
            // to ensure that the displacements are in declaration order for better readability
            // and maintainability of the code.
            for (std::size_t k = 1; k < N; ++k)
                assert(disp[k-1] < disp[k] && "StructMPITraits member list not in declaration order");
            MPI_Type_create_struct(N, blk.data(), disp.data(), types.data(), &tmp);
            // Resize the datatype to account for possible padding issues
            MPI_Type_create_resized(tmp, 0, sizeof(Struct), &type_);
            MPI_Type_commit(&type_);
            // Free the temporary datatype
            // Note: This is necessary to avoid memory leaks, as the temporary datatype
            //       created by MPI_Type_create_struct is not automatically freed.
            MPI_Type_free(&tmp);
        });
        return type_;
    }

    static constexpr bool is_intrinsic = false;

  private:
    // --- helper to recognise std::array ----------------------------------------
    template<class T> struct is_std_array : std::false_type {};
    template<class T, std::size_t N>
    struct is_std_array<std::array<T, N>> : std::true_type {};
    // Default case: not an enum
    template <typename T, typename Enable = void>
    struct MpiDispatch {
        using Type = MPITraits<T>;
    };
    // Specialization for enums
    template <typename T>
    struct MpiDispatch<T, typename std::enable_if<std::is_enum<T>::value>::type> {
        using Type = MPITraits<typename std::underlying_type<T>::type>;
    };

    template<auto Member, class Dummy>
    static void processMember(Dummy& d, MPI_Aint base,
                              std::array<MPI_Aint,     sizeof...(Members)>& disp,
                              std::array<int,          sizeof...(Members)>& blk,
                              std::array<MPI_Datatype, sizeof...(Members)>& types,
                              std::size_t idx)
    {
        using MemberT = std::remove_reference_t<decltype(d.*Member)>;

        MPI_Get_address(&(d.*Member), &disp[idx]);
        disp[idx] -= base;

        if constexpr (std::is_array_v<MemberT>) {
            // C array  T[N]
            using Elem = std::remove_extent_t<MemberT>;
            blk  [idx] = std::extent_v<MemberT>;
            types[idx] = MPITraits<Elem>::getType();
        }
        else if constexpr (is_std_array<MemberT>::value) {
            // std::array<T,N>
            using Elem = typename MemberT::value_type;
            blk [idx]  = std::tuple_size<MemberT>::value;
            types[idx] = MPITraits<Elem>::getType();
        }
        else {
            // scalar or enum
            blk  [idx] = 1;
            using MPIType = typename MpiDispatch<MemberT>::Type;
            types[idx] = MPIType::getType();
        }
    }

    // Initial value of MPI_DATATYPE_NULL is used to indicate that the type
    //  has not been created yet
    static inline MPI_Datatype type_ = MPI_DATATYPE_NULL;
    static inline std::once_flag flag_;
};

// Convenience alias for StructMPITraitsImpl
template<class Struct, auto... Members>
using StructMPITraits = StructMPITraitsImpl<Struct, Members...>;

} // namespace Dune::detail

// -----------------------------------------------------------------------
// MPI traits for ReservoirCoupling structs
// The fields of the structs should be listed in the order they are declared
// in the struct definition.
// -----------------------------------------------------------------------

// Trait for InjectionGroupTarget
template<class Scalar>
struct MPITraits<
        ::Opm::ReservoirCoupling::InjectionGroupTarget<Scalar>>
  : detail::StructMPITraits<
        ::Opm::ReservoirCoupling::InjectionGroupTarget<Scalar>,
        &::Opm::ReservoirCoupling::InjectionGroupTarget<Scalar>::group_name_idx,
        &::Opm::ReservoirCoupling::InjectionGroupTarget<Scalar>::target,
        &::Opm::ReservoirCoupling::InjectionGroupTarget<Scalar>::cmode,
        &::Opm::ReservoirCoupling::InjectionGroupTarget<Scalar>::phase>  { };

// Trait for ProductionGroupTarget
template<class Scalar>
struct MPITraits<
        ::Opm::ReservoirCoupling::ProductionGroupTarget<Scalar>>
  : detail::StructMPITraits<
        ::Opm::ReservoirCoupling::ProductionGroupTarget<Scalar>,
        &::Opm::ReservoirCoupling::ProductionGroupTarget<Scalar>::group_name_idx,
        &::Opm::ReservoirCoupling::ProductionGroupTarget<Scalar>::target,
        &::Opm::ReservoirCoupling::ProductionGroupTarget<Scalar>::cmode> { };

/// @brief MPI datatype trait for SlaveGroupProductionData structure
///
/// This specialization enables SlaveGroupProductionData to be sent via MPI
/// by defining how to construct an MPI_Datatype that matches the memory layout
/// of the structure. The trait uses StructMPITraits to automatically generate
/// the MPI datatype based on the listed member pointers.
///
/// @tparam Scalar Floating-point type (typically double or float)
///
/// @note This trait is used internally by MPI communication functions when
///       sending/receiving SlaveGroupProductionData between processes
/// @see StructMPITraits
/// @see SlaveGroupProductionData
template<class Scalar>
struct MPITraits<::Opm::ReservoirCoupling::SlaveGroupProductionData<Scalar>>
    : detail::StructMPITraits<
          ::Opm::ReservoirCoupling::SlaveGroupProductionData<Scalar>,
          &::Opm::ReservoirCoupling::SlaveGroupProductionData<Scalar>::potentials,
          &::Opm::ReservoirCoupling::SlaveGroupProductionData<Scalar>::surface_rates,
          &::Opm::ReservoirCoupling::SlaveGroupProductionData<Scalar>::network_surface_rates,
          &::Opm::ReservoirCoupling::SlaveGroupProductionData<Scalar>::reservoir_rates,
          &::Opm::ReservoirCoupling::SlaveGroupProductionData<Scalar>::voidage_rate,
          &::Opm::ReservoirCoupling::SlaveGroupProductionData<Scalar>::gas_reinjection_rate
        >  { };

/// @brief MPI datatype trait for SlaveGroupInjectionData structure
///
/// This specialization enables SlaveGroupInjectionData to be sent via MPI
/// by defining how to construct an MPI_Datatype that matches the memory layout
/// of the structure. The trait uses StructMPITraits to automatically generate
/// the MPI datatype based on the listed member pointers.
///
/// @tparam Scalar Floating-point type (typically double or float)
///
/// @note This trait is used internally by MPI communication functions when
///       sending/receiving SlaveGroupInjectionData between processes
/// @see StructMPITraits
/// @see SlaveGroupInjectionData
template<class Scalar>
struct MPITraits<::Opm::ReservoirCoupling::SlaveGroupInjectionData<Scalar>>
    : detail::StructMPITraits<
          ::Opm::ReservoirCoupling::SlaveGroupInjectionData<Scalar>,
          &::Opm::ReservoirCoupling::SlaveGroupInjectionData<Scalar>::surface_rates,
          &::Opm::ReservoirCoupling::SlaveGroupInjectionData<Scalar>::reservoir_rates
        >  { };
} // namespace Dune

#endif // OPM_RESERVOIR_COUPLING_MPI_TRAITS_HPP
