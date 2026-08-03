/*
  Copyright Equinor ASA 2026

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
#ifndef OPM_SYSTEMPRECONDITIONERFACTORY_HEADER_INCLUDED
#define OPM_SYSTEMPRECONDITIONERFACTORY_HEADER_INCLUDED

#include <opm/simulators/linalg/system/MultiComm.hpp>
#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/system/SystemTypes.hpp>

#include <dune/istl/operators.hh>
#include <dune/istl/paamg/pinfo.hh>

namespace Opm
{

template <class Operator, class Comm, typename>
struct StandardPreconditioners;

template<typename Scalar>
using SystemSeqOp = Dune::MatrixAdapter<SystemMatrix<Scalar>, SystemVector<Scalar>, SystemVector<Scalar>>;

#if HAVE_MPI
using SystemComm = Dune::MultiCommunicator<const Dune::OwnerOverlapCopyCommunication<int, int>&,
                                           const Dune::JacComm&>;
template<typename Scalar>
using SystemParOp = Dune::OverlappingSchwarzOperator<SystemMatrix<Scalar>, SystemVector<Scalar>,
                                                      SystemVector<Scalar>, SystemComm>;
#endif


template <>
struct StandardPreconditioners<SystemSeqOp<double>, Dune::Amg::SequentialInformation, void>
{
    static void add();
};

#if FLOW_INSTANTIATE_FLOAT
template <>
struct StandardPreconditioners<SystemSeqOp<float>, Dune::Amg::SequentialInformation, void>
{
    static void add();
};
#endif

#if HAVE_MPI
template <>
struct StandardPreconditioners<SystemParOp<double>, Dune::Amg::SequentialInformation, void>
{
    static void add();
};

template <>
struct StandardPreconditioners<SystemParOp<double>, SystemComm, void>
{
    static void add();
};

#if FLOW_INSTANTIATE_FLOAT
template <>
struct StandardPreconditioners<SystemParOp<float>, Dune::Amg::SequentialInformation, void>
{
    static void add();
};

template <>
struct StandardPreconditioners<SystemParOp<float>, SystemComm, void>
{
    static void add();
};

#endif
#endif

} // namespace Opm

#endif // OPM_SYSTEMPRECONDITIONERFACTORY_HEADER_INCLUDED
