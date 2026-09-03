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
#include <config.h>
#include <opm/simulators/linalg/system/SystemPreconditionerFactory.hpp>

#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>

#include <functional>
#include <memory>

namespace {

template<typename Scalar>
void addSystemCprSeq()
{
    using O = Opm::SystemSeqOp<Scalar>;
    using F = Opm::PreconditionerFactory<O, Dune::Amg::SequentialInformation>;
    using V = Opm::SystemVector<Scalar>;
    using P = Opm::PropertyTree;

    F::addCreator("system_cpr",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex) {
                      std::function<Opm::ResVector<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      return std::make_shared<Opm::SystemPreconditioner<Scalar, Opm::SeqResOperator<Scalar>>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm);
                  });
}

template<typename Scalar>
void addSystemCprwSeq()
{
    using O = Opm::SystemSeqOp<Scalar>;
    using F = Opm::PreconditionerFactory<O, Dune::Amg::SequentialInformation>;
    using V = Opm::SystemVector<Scalar>;
    using P = Opm::PropertyTree;

    F::addCreator("system_cprw",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex) {
                      std::function<Opm::ResVector<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      return std::make_shared<
                          Opm::SystemPreconditioner<Scalar, Opm::SeqResOperator<Scalar>,
                                                    Dune::Amg::SequentialInformation, true>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm);
                  });
}

#if HAVE_MPI
// Register a sequential (non-MPI) version of the system_cpr preconditioner
// for the parallel operator. This allows each MPI rank to apply a local,
// communication-free CPR preconditioner inside the overlapping Schwarz
// framework. It is a lightweight alternative to the fully parallel
// preconditioner registered by addSystemCprPar().
template<typename Scalar>
void addSystemCprParSeq()
{
    using O = Opm::SystemParOp<Scalar>;
    using F = Opm::PreconditionerFactory<O, Dune::Amg::SequentialInformation>;
    using V = Opm::SystemVector<Scalar>;
    using P = Opm::PropertyTree;

    F::addCreator("system_cpr",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex) {
                      std::function<Opm::ResVector<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      return std::make_shared<Opm::SystemPreconditioner<Scalar, Opm::SeqResOperator<Scalar>>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm);
                  });
}

template<typename Scalar>
void addSystemCprwParSeq()
{
    using O = Opm::SystemParOp<Scalar>;
    using F = Opm::PreconditionerFactory<O, Dune::Amg::SequentialInformation>;
    using V = Opm::SystemVector<Scalar>;
    using P = Opm::PropertyTree;

    F::addCreator("system_cprw",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex) {
                      std::function<Opm::ResVector<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      return std::make_shared<
                          Opm::SystemPreconditioner<Scalar, Opm::SeqResOperator<Scalar>,
                                                    Dune::Amg::SequentialInformation, true>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm);
                  });
}

template<typename Scalar>
void addSystemCprPar()
{
    using O = Opm::SystemParOp<Scalar>;
    using F = Opm::PreconditionerFactory<O, Opm::SystemComm>;
    using V = Opm::SystemVector<Scalar>;
    using P = Opm::PropertyTree;

    F::addCreator("system_cpr",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex,
                     const Opm::SystemComm& comm) {
                      std::function<Opm::ResVector<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      const auto& resComm = comm[Dune::Indices::_0];
                      return std::make_shared<
                          Opm::SystemPreconditioner<Scalar, Opm::ParResOperator<Scalar>, Opm::ParResComm>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm, resComm);
                  });
}

template<typename Scalar>
void addSystemCprwPar()
{
    using O = Opm::SystemParOp<Scalar>;
    using F = Opm::PreconditionerFactory<O, Opm::SystemComm>;
    using V = Opm::SystemVector<Scalar>;
    using P = Opm::PropertyTree;

    F::addCreator("system_cprw",
                  [](const O& op, const P& prm,
                     const std::function<V()>& sysWeightCalc,
                     std::size_t pressureIndex,
                     const Opm::SystemComm& comm) {
                      std::function<Opm::ResVector<Scalar>()> resWeightCalc;
                      if (sysWeightCalc) {
                          resWeightCalc = [sysWeightCalc]() {
                              return sysWeightCalc()[Dune::Indices::_0];
                          };
                      }
                      const auto& resComm = comm[Dune::Indices::_0];
                      return std::make_shared<
                          Opm::SystemPreconditioner<Scalar, Opm::ParResOperator<Scalar>, Opm::ParResComm, true>>(
                          op.getmat(), resWeightCalc, pressureIndex, prm, resComm);
                  });
}
#endif

} // namespace

namespace Opm {

void StandardPreconditioners<SystemSeqOp<double>, Dune::Amg::SequentialInformation, void>::add()
{
    addSystemCprSeq<double>();
    addSystemCprwSeq<double>();
}

#if FLOW_INSTANTIATE_FLOAT
void StandardPreconditioners<SystemSeqOp<float>, Dune::Amg::SequentialInformation, void>::add()
{
    addSystemCprSeq<float>();
    addSystemCprwSeq<float>();
}
#endif

#if HAVE_MPI
void StandardPreconditioners<SystemParOp<double>, Dune::Amg::SequentialInformation, void>::add()
{
    addSystemCprParSeq<double>();
    addSystemCprwParSeq<double>();
}

void StandardPreconditioners<SystemParOp<double>, SystemComm, void>::add()
{
    addSystemCprPar<double>();
    addSystemCprwPar<double>();
}

#if FLOW_INSTANTIATE_FLOAT
void StandardPreconditioners<SystemParOp<float>, Dune::Amg::SequentialInformation, void>::add()
{
    addSystemCprParSeq<float>();
    addSystemCprwParSeq<float>();
}

void StandardPreconditioners<SystemParOp<float>, SystemComm, void>::add()
{
    addSystemCprPar<float>();
    addSystemCprwPar<float>();
}
#endif
#endif

} // namespace Opm
