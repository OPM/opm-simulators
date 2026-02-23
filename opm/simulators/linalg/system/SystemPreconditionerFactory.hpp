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

template<typename Scalar>
using SystemSeqOpT = Dune::MatrixAdapter<SystemMatrixT<Scalar>, SystemVectorT<Scalar>, SystemVectorT<Scalar>>;

#if HAVE_MPI
using SystemComm = Dune::MultiCommunicator<const Dune::OwnerOverlapCopyCommunication<int, int>&,
                                           const Dune::JacComm&>;
template<typename Scalar>
using SystemParOpT = Dune::OverlappingSchwarzOperator<SystemMatrixT<Scalar>, SystemVectorT<Scalar>,
                                                      SystemVectorT<Scalar>, SystemComm>;
#endif

// Full specialisations of StandardPreconditioners for the coupled system
// operators.  Partial specialisations would be ambiguous with the generic
// serial factory (!is_gpu_operator_v guard), so full specialisations are
// used; detail:: helpers factor out the shared logic.

namespace detail {

template<typename Scalar>
void addSystemCprSeq()
{
    using O = SystemSeqOpT<Scalar>;
    using F = PreconditionerFactory<O, Dune::Amg::SequentialInformation>;
    using V = SystemVectorT<Scalar>;
    using P = PropertyTree;

    F::addCreator("system_cpr",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex) {
                      std::function<ResVectorT<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      return std::make_shared<SystemPreconditioner<Scalar, SeqResOperatorT<Scalar>>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm);
                  });
}

#if HAVE_MPI
template<typename Scalar>
void addSystemCprParSeq()
{
    using O = SystemParOpT<Scalar>;
    using F = PreconditionerFactory<O, Dune::Amg::SequentialInformation>;
    using V = SystemVectorT<Scalar>;
    using P = PropertyTree;

    F::addCreator("system_cpr",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex) {
                      std::function<ResVectorT<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      return std::make_shared<SystemPreconditioner<Scalar, SeqResOperatorT<Scalar>>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm);
                  });
}

template<typename Scalar>
void addSystemCprPar()
{
    using O = SystemParOpT<Scalar>;
    using F = PreconditionerFactory<O, SystemComm>;
    using V = SystemVectorT<Scalar>;
    using P = PropertyTree;

    F::addCreator("system_cpr",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex,
                     const SystemComm& comm) {
                      std::function<ResVectorT<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      const auto& resComm = comm[Dune::Indices::_0];
                      return std::make_shared<SystemPreconditioner<Scalar, ParResOperatorT<Scalar>, ParResComm>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm, resComm);
                  });
}
#endif

} // namespace detail

template <>
struct StandardPreconditioners<SystemSeqOpT<double>, Dune::Amg::SequentialInformation, void> {
    static void add() { detail::addSystemCprSeq<double>(); }
};

template <>
struct StandardPreconditioners<SystemSeqOpT<float>, Dune::Amg::SequentialInformation, void> {
    static void add() { detail::addSystemCprSeq<float>(); }
};

#if HAVE_MPI
template <>
struct StandardPreconditioners<SystemParOpT<double>, Dune::Amg::SequentialInformation, void> {
    static void add() { detail::addSystemCprParSeq<double>(); }
};

template <>
struct StandardPreconditioners<SystemParOpT<float>, Dune::Amg::SequentialInformation, void> {
    static void add() { detail::addSystemCprParSeq<float>(); }
};

template <>
struct StandardPreconditioners<SystemParOpT<double>, SystemComm, void> {
    static void add() { detail::addSystemCprPar<double>(); }
};

template <>
struct StandardPreconditioners<SystemParOpT<float>, SystemComm, void> {
    static void add() { detail::addSystemCprPar<float>(); }
};
#endif

} // namespace Opm
