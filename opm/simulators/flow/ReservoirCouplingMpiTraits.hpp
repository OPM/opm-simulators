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

#include <opm/simulators/flow/ReservoirCoupling.hpp>

#include <array>
#include <mutex>  // For std::call_once
#include <tuple>   // for std::tuple_size

#include <mpi.h>

namespace Dune {

namespace detail {

// -----------------------------------------------------------------------
//  StructMPITraitsImpl  --  generic MPI traits implementation for structs
//
// Assumes that each field of the struct is either
//   - an already defined MPI type, or
//   - is an array or std::array of an already defined MPI type, or
//   - is an enum with underlying type that is an already defined MPI type
// ------------------------------------------------------------------------
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

// Trait for Potentials
template<>
struct MPITraits<::Opm::ReservoirCoupling::Potentials>
    : detail::StructMPITraits<
          ::Opm::ReservoirCoupling::Potentials,
          &::Opm::ReservoirCoupling::Potentials::rate>  { };

} // namespace Dune

#endif // OPM_RESERVOIR_COUPLING_MPI_TRAITS_HPP
