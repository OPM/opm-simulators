#include <config.h>
#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>
#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>
#include <opm/simulators/linalg/system/SystemPreconditionerFactory.hpp>

// Explicit instantiations

template class Opm::SystemPreconditioner<Opm::SeqResOperator>;

template class Dune::FlexibleSolver<
    Dune::MatrixAdapter<Opm::WWMatrix, Opm::WellVector, Opm::WellVector>>;

template class Dune::FlexibleSolver<Opm::SystemSeqOp>;
template class Opm::PreconditionerFactory<Opm::SystemSeqOp, Dune::Amg::SequentialInformation>;

#if HAVE_MPI
template class Opm::SystemPreconditioner<Opm::ParResOperator, Opm::ParResComm>;

template class Dune::FlexibleSolver<Opm::SystemParOp>;
template Dune::FlexibleSolver<Opm::SystemParOp>::FlexibleSolver(
    Opm::SystemParOp& op,
    const Opm::SystemComm& comm,
    const Opm::PropertyTree& prm,
    const std::function<Opm::SystemVector()>& weightsCalculator,
    std::size_t pressureIndex);
template class Opm::PreconditionerFactory<Opm::SystemParOp, Opm::SystemComm>;
template class Opm::PreconditionerFactory<Opm::SystemParOp, Dune::Amg::SequentialInformation>;
#endif
