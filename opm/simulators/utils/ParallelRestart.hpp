/*
  Copyright 2019 Equinor AS.

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
#ifndef PARALLEL_RESTART_HPP
#define PARALLEL_RESTART_HPP

#if HAVE_MPI
#include <mpi.h>
#endif

#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/Summary.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <vector>
#include <map>
#include <unordered_map>

namespace Opm
{

class EDITNNC;
class NNC;
struct NNCdata;
class Rock2dTable;
class Rock2dtrTable;
class ThresholdPressure;

namespace Mpi
{
template<class T>
std::size_t packSize(const T*, std::size_t, Dune::MPIHelper::MPICommunicator,
                     std::integral_constant<bool, false>);

template<class T>
std::size_t packSize(const T*, std::size_t l, Dune::MPIHelper::MPICommunicator comm,
                     std::integral_constant<bool, true>);

template<class T>
std::size_t packSize(const T* data, std::size_t l, Dune::MPIHelper::MPICommunicator comm);

template<class T>
std::size_t packSize(const T&, Dune::MPIHelper::MPICommunicator,
                     std::integral_constant<bool, false>);

template<class T>
std::size_t packSize(const T&, Dune::MPIHelper::MPICommunicator comm,
                     std::integral_constant<bool, true>);

template<class T>
std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2>
std::size_t packSize(const std::pair<T1,T2>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T, class A>
std::size_t packSize(const std::vector<T,A>& data, Dune::MPIHelper::MPICommunicator comm);

std::size_t packSize(const char* str, Dune::MPIHelper::MPICommunicator comm);

std::size_t packSize(const std::string& str, Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class C, class A>
std::size_t packSize(const std::map<T1,T2,C,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
std::size_t packSize(const std::unordered_map<T1,T2,H,P,A>& data, Dune::MPIHelper::MPICommunicator comm);

////// pack routines

template<class T>
void pack(const T*, std::size_t, std::vector<char>&, int&,
          Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>);

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T>
void pack(const T&, std::vector<char>&, int&,
          Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>);


template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2>
void pack(const std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T, class A>
void pack(const std::vector<T,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class C, class A>
void pack(const std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
void pack(const std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

void pack(const char* str, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

/// unpack routines

template<class T>
void unpack(T*, const std::size_t&, std::vector<char>&, int&,
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm,
            std::integral_constant<bool, true>);

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T>
void unpack(T&, std::vector<char>&, int&,
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>);

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2>
void unpack(std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T, class A>
void unpack(std::vector<T,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class C, class A>
void unpack(std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
void unpack(std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

void unpack(char* str, std::size_t length, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

/// prototypes for complex types

#define ADD_PACK_PROTOTYPES(T) \
  std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm); \
  void pack(const T& data, std::vector<char>& buffer, int& position, \
          Dune::MPIHelper::MPICommunicator comm); \
  void unpack(T& data, std::vector<char>& buffer, int& position, \
              Dune::MPIHelper::MPICommunicator comm);

ADD_PACK_PROTOTYPES(data::CellData)
ADD_PACK_PROTOTYPES(data::Connection)
ADD_PACK_PROTOTYPES(data::Rates)
ADD_PACK_PROTOTYPES(data::Segment)
ADD_PACK_PROTOTYPES(data::Solution)
ADD_PACK_PROTOTYPES(data::Well)
ADD_PACK_PROTOTYPES(data::WellRates)
ADD_PACK_PROTOTYPES(EDITNNC)
ADD_PACK_PROTOTYPES(NNC)
ADD_PACK_PROTOTYPES(NNCdata)
ADD_PACK_PROTOTYPES(RestartKey)
ADD_PACK_PROTOTYPES(RestartValue)
ADD_PACK_PROTOTYPES(Rock2dTable)
ADD_PACK_PROTOTYPES(Rock2dtrTable)
ADD_PACK_PROTOTYPES(std::string)
ADD_PACK_PROTOTYPES(ThresholdPressure)

} // end namespace Mpi
RestartValue loadParallelRestart(const EclipseIO* eclIO, SummaryState& summaryState,
                                 const std::vector<Opm::RestartKey>& solutionKeys,
                                 const std::vector<Opm::RestartKey>& extraKeys,
                                 Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm);

} // end namespace Opm
#endif // PARALLEL_RESTART_HPP
