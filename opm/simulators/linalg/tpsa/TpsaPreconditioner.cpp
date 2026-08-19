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
#include <config.h>

#include <opm/simulators/linalg/tpsa/TpsaPreconditioner_impl.hpp>
#include <opm/simulators/linalg/tpsa/TpsaPreconditionerFactory.hpp>

#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>
#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>

#define INSTANTIATE_TPSA_PF_SEQ(T)                                          \
    template class Opm::TpsaPreconditioner<T,                               \
                                           Opm::SeqDispDispOperatorT<T>,    \
                                           Opm::SeqRotRotOperatorT<T>,      \
                                           Opm::SeqSPresSPresOperatorT<T>>; \
    template class Dune::FlexibleSolver<Opm::TpsaSeqOperator<T>>;           \
    template class Opm::PreconditionerFactory<Opm::TpsaSeqOperator<T>,      \
                                              Dune::Amg::SequentialInformation>;

#if HAVE_MPI
#define INSTANTIATE_TPSA_PF_PAR(T)                                                    \
    template class Opm::TpsaPreconditioner<T,                                         \
                                           Opm::ParDispDispOperatorT<T>,              \
                                           Opm::ParRotRotOperatorT<T>,                \
                                           Opm::ParSPresSPresOperatorT<T>,            \
                                           Opm::TpsaParComm>;                         \
    template class Dune::FlexibleSolver<Opm::TpsaParOperator<T>>;                     \
    template Dune::FlexibleSolver<Opm::TpsaParOperator<T>>::FlexibleSolver(           \
        Opm::TpsaParOperator<T>& op,                                                  \
        const Opm::TpsaComm& comm,                                                    \
        const Opm::PropertyTree& prm,                                                 \
        const std::function<Opm::Linear::TpsaMultiVector<T>()>& weightsCalculator,    \
        std::size_t pressureIndex);                                                   \
    template class Opm::PreconditionerFactory<Opm::TpsaParOperator<T>, Opm::TpsaComm>;\
    template class Opm::PreconditionerFactory<Opm::TpsaParOperator<T>,                \
                                              Dune::Amg::SequentialInformation>;

#define INSTANTIATE_TPSA_PF(T) \
    INSTANTIATE_TPSA_PF_PAR(T) \
    INSTANTIATE_TPSA_PF_SEQ(T)

#else
#define INSTANTIATE_TPSA_PF(T) INSTANTIATE_TPSA_PF_SEQ(T)
#endif

INSTANTIATE_TPSA_PF(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TPSA_PF(float)
#endif
