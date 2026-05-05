#include <config.h>
#include <opm/simulators/linalg/FlexibleSolver_impl.hpp>
#include <opm/simulators/linalg/PreconditionerFactory_impl.hpp>
#include <opm/simulators/linalg/system/SystemPreconditioner.hpp>
#include <opm/simulators/linalg/system/SystemPreconditionerFactory.hpp>

#define INSTANTIATE_SYSTEM_PF_SEQ(T)                                                                  \
    template class Opm::SystemPreconditioner<T, Opm::SeqResOperatorT<T>>;                             \
    template class Dune::FlexibleSolver<                                                               \
        Dune::MatrixAdapter<Opm::WWMatrixT<T>, Opm::WellVectorT<T>, Opm::WellVectorT<T>>>;           \
    template class Dune::FlexibleSolver<Opm::SystemSeqOpT<T>>;                                        \
    template class Opm::PreconditionerFactory<Opm::SystemSeqOpT<T>, Dune::Amg::SequentialInformation>;

#if HAVE_MPI
#define INSTANTIATE_SYSTEM_PF_PAR(T)                                                                  \
    template class Opm::SystemPreconditioner<T, Opm::ParResOperatorT<T>, Opm::ParResComm>;            \
    template class Dune::FlexibleSolver<Opm::SystemParOpT<T>>;                                        \
    template Dune::FlexibleSolver<Opm::SystemParOpT<T>>::FlexibleSolver(                               \
        Opm::SystemParOpT<T>& op,                                                                     \
        const Opm::SystemComm& comm,                                                                  \
        const Opm::PropertyTree& prm,                                                                 \
        const std::function<Opm::SystemVectorT<T>()>& weightsCalculator,                              \
        std::size_t pressureIndex);                                                                    \
    template class Opm::PreconditionerFactory<Opm::SystemParOpT<T>, Opm::SystemComm>;                 \
    template class Opm::PreconditionerFactory<Opm::SystemParOpT<T>, Dune::Amg::SequentialInformation>;

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
