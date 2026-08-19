// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 NORCE AS

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
#ifndef OPM_TPSA_PRECONDITIONER_FACTORY_HPP
#define OPM_TPSA_PRECONDITIONER_FACTORY_HPP

#include <dune/istl/operators.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/system/MultiComm.hpp>
#include <opm/simulators/linalg/tpsa/TpsaPreconditioner.hpp>
#include <opm/simulators/linalg/tpsa/TpsaTypes.hpp>

#include <cstddef>
#include <functional>
#include <memory>

namespace Opm
{

template <class Operator, class Comm, typename>
struct StandardPreconditioners;

template <typename Scalar>
using TpsaSeqOperator = Dune::MatrixAdapter<Linear::TpsaMatrixView<Scalar>,
                                            Linear::TpsaMultiVector<Scalar>,
                                            Linear::TpsaMultiVector<Scalar> >;

#if HAVE_MPI
//! \brief One communicator per field; in practice the same one five times.
using TpsaComm = Dune::MultiCommunicator<const Dune::OwnerOverlapCopyCommunication<int, int>&,
                                         const Dune::OwnerOverlapCopyCommunication<int, int>&,
                                         const Dune::OwnerOverlapCopyCommunication<int, int>&,
                                         const Dune::OwnerOverlapCopyCommunication<int, int>&,
                                         const Dune::OwnerOverlapCopyCommunication<int, int>&>;

template <typename Scalar>
using TpsaParOperator = Dune::OverlappingSchwarzOperator<Linear::TpsaMatrixView<Scalar>,
                                                         Linear::TpsaMultiVector<Scalar>,
                                                         Linear::TpsaMultiVector<Scalar>,
                                                         TpsaComm>;
#endif

namespace detail
{

    template <typename Scalar>
    void
    addTpsaBlockSeq()
    {
        using O = TpsaSeqOperator<Scalar>;
        using F = PreconditionerFactory<O, Dune::Amg::SequentialInformation>;
        using V = Linear::TpsaMultiVector<Scalar>;
        using P = PropertyTree;

        F::addCreator("tpsa_block",
                      [](const O& op,
                         const P& prm,
                         const std::function<V()>&,
                         std::size_t) {
                          using PreCond = TpsaPreconditioner<Scalar,
                                                             SeqDispDispOperatorT<Scalar>,
                                                             SeqRotRotOperatorT<Scalar>,
                                                             SeqSPresSPresOperatorT<Scalar> >;
                          return std::make_shared<PreCond>(op.getmat(), prm);
                      });
    }

#if HAVE_MPI
    template <typename Scalar>
    void
    addTpsaBlockParSeq()
    {
        // An MPI build running on a single rank still instantiates the parallel
        // operator, but with sequential information.
        using O = TpsaParOperator<Scalar>;
        using F = PreconditionerFactory<O, Dune::Amg::SequentialInformation>;
        using V = Linear::TpsaMultiVector<Scalar>;
        using P = PropertyTree;

        F::addCreator("tpsa_block",
                      [](const O& op,
                         const P& prm,
                         const std::function<V()>&,
                         std::size_t) {
                          using PreCond = TpsaPreconditioner<Scalar,
                                                             SeqDispDispOperatorT<Scalar>,
                                                             SeqRotRotOperatorT<Scalar>,
                                                             SeqSPresSPresOperatorT<Scalar> >;
                          return std::make_shared<PreCond>(op.getmat(), prm);
                      });
    }

    template <typename Scalar>
    void
    addTpsaBlockPar()
    {
        using O = TpsaParOperator<Scalar>;
        using F = PreconditionerFactory<O, TpsaComm>;
        using V = Linear::TpsaMultiVector<Scalar>;
        using P = PropertyTree;

        F::addCreator("tpsa_block",
                      [](const O& op,
                         const P& prm,
                         const std::function<V()>&,
                         std::size_t,
                         const TpsaComm& comm) {
                          // All five fields share one dof partition.
                          const auto& fieldComm = comm[Dune::Indices::_0];
                          using PreCond = TpsaPreconditioner<Scalar,
                                                             ParDispDispOperatorT<Scalar>,
                                                             ParRotRotOperatorT<Scalar>,
                                                             ParSPresSPresOperatorT<Scalar>,
                                                             TpsaParComm>;
                          return std::make_shared<PreCond>(op.getmat(), prm, fieldComm);
                      });
    }
#endif

} // namespace detail

template <>
struct StandardPreconditioners<TpsaSeqOperator<double>, Dune::Amg::SequentialInformation, void>
{
    static void add()
    {
        detail::addTpsaBlockSeq<double>();
    }
};

template <>
struct StandardPreconditioners<TpsaSeqOperator<float>, Dune::Amg::SequentialInformation, void>
{
    static void add()
    {
        detail::addTpsaBlockSeq<float>();
    }
};

#if HAVE_MPI
template <>
struct StandardPreconditioners<TpsaParOperator<double>, Dune::Amg::SequentialInformation, void>
{
    static void add()
    {
        detail::addTpsaBlockParSeq<double>();
    }
};

template <>
struct StandardPreconditioners<TpsaParOperator<float>, Dune::Amg::SequentialInformation, void>
{
    static void add()
    {
        detail::addTpsaBlockParSeq<float>();
    }
};

template <>
struct StandardPreconditioners<TpsaParOperator<double>, TpsaComm, void>
{
    static void add()
    {
        detail::addTpsaBlockPar<double>();
    }
};

template <>
struct StandardPreconditioners<TpsaParOperator<float>, TpsaComm, void>
{
    static void add()
    {
        detail::addTpsaBlockPar<float>();
    }
};
#endif

} // namespace Opm

#endif // OPM_TPSA_PRECONDITIONER_FACTORY_HPP
