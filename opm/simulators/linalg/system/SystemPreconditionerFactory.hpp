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
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>
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

// Full specialisations of StandardPreconditioners for the coupled system
// operators.  Partial specialisations would be ambiguous with the generic
// serial factory (!is_gpu_operator_v guard), so full specialisations are
// used; detail:: helpers factor out the shared logic.

namespace detail {

template<typename Scalar>
void addSystemCprSeq()
{
    using O = SystemSeqOp<Scalar>;
    using F = PreconditionerFactory<O, Dune::Amg::SequentialInformation>;
    using V = SystemVector<Scalar>;
    using P = PropertyTree;

    F::addCreator("system_cpr",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex) {
                      std::function<ResVector<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      return std::make_shared<SystemPreconditioner<Scalar, SeqResOperator<Scalar>>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm);
                  });
}

#if HAVE_MPI
// Register a sequential (non-MPI) version of the system_cpr preconditioner
// for the parallel operator.  This allows each MPI rank to apply a local,
// communication‑free CPR preconditioner inside the overlapping Schwarz
// framework.  It is a lightweight alternative to the fully parallel
// preconditioner registered by addSystemCprPar().
template<typename Scalar>
void addSystemCprParSeq()
{
    using O = SystemParOp<Scalar>;
    using F = PreconditionerFactory<O, Dune::Amg::SequentialInformation>;
    using V = SystemVector<Scalar>;
    using P = PropertyTree;

    F::addCreator("system_cpr",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex) {
                      std::function<ResVector<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      return std::make_shared<SystemPreconditioner<Scalar, SeqResOperator<Scalar>>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm);
                  });
}

template<typename Scalar>
void addSystemCprPar()
{
    using O = SystemParOp<Scalar>;
    using F = PreconditionerFactory<O, SystemComm>;
    using V = SystemVector<Scalar>;
    using P = PropertyTree;

    F::addCreator("system_cpr",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex,
                     const SystemComm& comm) {
                      std::function<ResVector<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      const auto& resComm = comm[Dune::Indices::_0];
                      return std::make_shared<SystemPreconditioner<Scalar, ParResOperator<Scalar>, ParResComm>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm, resComm);
                  });
}
#endif

} // namespace detail

template <>
struct StandardPreconditioners<SystemSeqOp<double>, Dune::Amg::SequentialInformation, void> {
    static void add() { detail::addSystemCprSeq<double>(); }
};

template <>
struct StandardPreconditioners<SystemSeqOp<float>, Dune::Amg::SequentialInformation, void> {
    static void add() { detail::addSystemCprSeq<float>(); }
};

#if HAVE_MPI
template <>
struct StandardPreconditioners<SystemParOp<double>, Dune::Amg::SequentialInformation, void> {
    static void add() { detail::addSystemCprParSeq<double>(); }
};

template <>
struct StandardPreconditioners<SystemParOp<float>, Dune::Amg::SequentialInformation, void> {
    static void add() { detail::addSystemCprParSeq<float>(); }
};

template <>
struct StandardPreconditioners<SystemParOp<double>, SystemComm, void> {
    static void add() { detail::addSystemCprPar<double>(); }
};

template <>
struct StandardPreconditioners<SystemParOp<float>, SystemComm, void> {
    static void add() { detail::addSystemCprPar<float>(); }
};
#endif

} // namespace Opm

#endif // OPM_SYSTEMPRECONDITIONERFACTORY_HEADER_INCLUDED