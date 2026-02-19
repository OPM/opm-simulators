#pragma once

#include "MultiComm.hpp"
#include "SystemPreconditioner.hpp"
#include "SystemTypes.hpp"

#include <opm/simulators/linalg/PreconditionerFactory.hpp>

#include <dune/istl/operators.hh>
#include <dune/istl/paamg/pinfo.hh>

namespace Opm
{

template <class Operator, class Comm, typename>
struct StandardPreconditioners;

using SystemSeqOp = Dune::MatrixAdapter<SystemMatrix, SystemVector, SystemVector>;

#if HAVE_MPI
using SystemComm = Dune::MultiCommunicator<const Dune::OwnerOverlapCopyCommunication<int, int>&,
                                           const Dune::JacComm&>;
using SystemParOp
    = Dune::OverlappingSchwarzOperator<SystemMatrix, SystemVector, SystemVector, SystemComm>;
#endif

template <>
struct StandardPreconditioners<SystemSeqOp, Dune::Amg::SequentialInformation, void> {
    static void add()
    {
        using O = SystemSeqOp;
        using C = Dune::Amg::SequentialInformation;
        using F = PreconditionerFactory<O, C>;
        using V = SystemVector;
        using P = PropertyTree;

        F::addCreator("system_cpr",
                      [](const O& op,
                         const P& prm,
                         const std::function<V()>& sysWeightCalc,
                         std::size_t pressureIndex) {
                          std::function<ResVector()> resWeightCalc;
                          if (sysWeightCalc) {
                              resWeightCalc = [sysWeightCalc]() {
                                  return sysWeightCalc()[Dune::Indices::_0];
                              };
                          }
                          return std::make_shared<SystemPreconditioner<SeqResOperator>>(
                              op.getmat(), resWeightCalc, pressureIndex, prm);
                      });
    }
};

#if HAVE_MPI
// Sequential initOpPrecSp overload instantiates PreconditionerFactory<SystemParOp, SeqInfo>
template <>
struct StandardPreconditioners<SystemParOp, Dune::Amg::SequentialInformation, void> {
    static void add()
    {
        using O = SystemParOp;
        using C = Dune::Amg::SequentialInformation;
        using F = PreconditionerFactory<O, C>;
        using V = SystemVector;
        using P = PropertyTree;

        F::addCreator("system_cpr",
                      [](const O& op,
                         const P& prm,
                         const std::function<V()>& sysWeightCalc,
                         std::size_t pressureIndex) {
                          std::function<ResVector()> resWeightCalc;
                          if (sysWeightCalc) {
                              resWeightCalc = [sysWeightCalc]() {
                                  return sysWeightCalc()[Dune::Indices::_0];
                              };
                          }
                          return std::make_shared<SystemPreconditioner<SeqResOperator>>(
                              op.getmat(), resWeightCalc, pressureIndex, prm);
                      });
    }
};

template <>
struct StandardPreconditioners<SystemParOp, SystemComm, void> {
    static void add()
    {
        using O = SystemParOp;
        using F = PreconditionerFactory<O, SystemComm>;
        using V = SystemVector;
        using P = PropertyTree;

        F::addCreator("system_cpr",
                      [](const O& op,
                         const P& prm,
                         const std::function<V()>& sysWeightCalc,
                         std::size_t pressureIndex,
                         const SystemComm& comm) {
                          std::function<ResVector()> resWeightCalc;
                          if (sysWeightCalc) {
                              resWeightCalc = [sysWeightCalc]() {
                                  return sysWeightCalc()[Dune::Indices::_0];
                              };
                          }
                          const auto& resComm = comm[Dune::Indices::_0];
                          return std::make_shared<SystemPreconditioner<ParResOperator, ParResComm>>(
                              op.getmat(), resWeightCalc, pressureIndex, prm, resComm);
                      });
    }
};
#endif

} // namespace Opm
