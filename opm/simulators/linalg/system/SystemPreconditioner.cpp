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
#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>
#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>
#include <opm/simulators/linalg/system/SystemPreconditionerFactory.hpp>

#define INSTANTIATE_SYSTEM_PF_SEQ(T)                                                                  \
    template class Opm::SystemPreconditioner<T, Opm::SeqResOperator<T>>;                             \
    template class Dune::FlexibleSolver<                                                               \
        Dune::MatrixAdapter<Opm::WWMatrix<T>, Opm::WellVector<T>, Opm::WellVector<T>>>;           \
    template class Dune::FlexibleSolver<Opm::SystemSeqOp<T>>;                                        \
    template class Opm::PreconditionerFactory<Opm::SystemSeqOp<T>, Dune::Amg::SequentialInformation>;

#if HAVE_MPI
#define INSTANTIATE_SYSTEM_PF_PAR(T)                                                                  \
    template class Opm::SystemPreconditioner<T, Opm::ParResOperator<T>, Opm::ParResComm>;            \
    template class Opm::SystemPreconditioner<T, Opm::ParResOperator<T>, Opm::ParResComm, true>;      \
    template class Dune::FlexibleSolver<Opm::SystemParOp<T>>;                                        \
    template Dune::FlexibleSolver<Opm::SystemParOp<T>>::FlexibleSolver(                               \
        Opm::SystemParOp<T>& op,                                                                     \
        const Opm::SystemComm& comm,                                                                  \
        const Opm::PropertyTree& prm,                                                                 \
        const std::function<Opm::SystemVector<T>()>& weightsCalculator,                              \
        std::size_t pressureIndex);                                                                    \
    template class Opm::PreconditionerFactory<Opm::SystemParOp<T>, Opm::SystemComm>;                 \
    template class Opm::PreconditionerFactory<Opm::SystemParOp<T>, Dune::Amg::SequentialInformation>;

#define INSTANTIATE_SYSTEM_PF(T) \
    INSTANTIATE_SYSTEM_PF_PAR(T) \
    INSTANTIATE_SYSTEM_PF_SEQ(T)
#else
#define INSTANTIATE_SYSTEM_PF(T) INSTANTIATE_SYSTEM_PF_SEQ(T)
#endif

INSTANTIATE_SYSTEM_PF(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_SYSTEM_PF(float)
#endif
