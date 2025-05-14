/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/simulators/linalg/PreconditionerFactory.hpp>

#include <opm/simulators/linalg/DILU.hpp>
#include <opm/simulators/linalg/ExtraSmoothers.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/OwningBlockPreconditioner.hpp>
#include <opm/simulators/linalg/OwningTwoLevelPreconditioner.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/PressureBhpTransferPolicy.hpp>
#include <opm/simulators/linalg/PressureTransferPolicy.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/WellOperators.hpp>
#include <opm/simulators/linalg/amgcpr.hh>
#include <opm/simulators/linalg/ilufirstelement.hh>
#include <opm/simulators/linalg/matrixblock.hh>

#include <dune/common/unused.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/paamg/kamg.hh>
#include <dune/istl/preconditioners.hh>

// Include all cuistl/GPU preconditioners inside of this headerfile
#include <opm/simulators/linalg/PreconditionerFactoryGPUIncludeWrapper.hpp>

#if HAVE_HYPRE
#include <opm/simulators/linalg/HyprePreconditioner.hpp>
#endif

#if HAVE_AMGX
#include <opm/simulators/linalg/AmgxPreconditioner.hpp>
#endif

#include <cassert>

#include <opm/simulators/linalg/StandardPreconditioners.hpp>

namespace Opm {



template <class Operator, class Comm>
PreconditionerFactory<Operator, Comm>::PreconditionerFactory()
{
}


template <class Operator, class Comm>
PreconditionerFactory<Operator, Comm>&
PreconditionerFactory<Operator, Comm>::instance()
{
    static PreconditionerFactory singleton;
    return singleton;
}

template <class Operator, class Comm>
typename PreconditionerFactory<Operator, Comm>::PrecPtr
PreconditionerFactory<Operator, Comm>::doCreate(const Operator& op,
                                                const PropertyTree& prm,
                                                const std::function<Vector()> weightsCalculator,
                                                std::size_t pressureIndex)
{
    if (!defAdded_) {
        StandardPreconditioners<Operator, Comm>::add();
        defAdded_ = true;
    }
    const std::string& type = prm.get<std::string>("type", "ParOverILU0");
    auto it = creators_.find(type);
    if (it == creators_.end()) {
        std::ostringstream msg;
        msg << "Preconditioner type " << type << " is not registered in the factory. Available types are: ";
        for (const auto& prec : creators_) {
            msg << prec.first << ' ';
        }
        msg << std::endl;
        OPM_THROW(std::invalid_argument, msg.str());
    }
    return it->second(op, prm, weightsCalculator, pressureIndex);
}

template <class Operator, class Comm>
typename PreconditionerFactory<Operator, Comm>::PrecPtr
PreconditionerFactory<Operator, Comm>::doCreate(const Operator& op,
                                                const PropertyTree& prm,
                                                const std::function<Vector()> weightsCalculator,
                                                std::size_t pressureIndex,
                                                const Comm& comm)
{
    if (!defAdded_) {
        StandardPreconditioners<Operator, Comm>::add();
        defAdded_ = true;
    }
    const std::string& type = prm.get<std::string>("type", "ParOverILU0");
    auto it = parallel_creators_.find(type);
    if (it == parallel_creators_.end()) {
        std::ostringstream msg;
        msg << "Parallel preconditioner type " << type << " is not registered in the factory. Available types are: ";
        for (const auto& prec : parallel_creators_) {
            msg << prec.first << ' ';
        }
        msg << std::endl;
        OPM_THROW(std::invalid_argument, msg.str());
    }
    return it->second(op, prm, weightsCalculator, pressureIndex, comm);
}

template <class Operator, class Comm>
void
PreconditionerFactory<Operator, Comm>::doAddCreator(const std::string& type, Creator c)
{
    creators_[type] = c;
}

template <class Operator, class Comm>
void
PreconditionerFactory<Operator, Comm>::doAddCreator(const std::string& type, ParCreator c)
{
    parallel_creators_[type] = c;
}

template <class Operator, class Comm>
typename PreconditionerFactory<Operator, Comm>::PrecPtr
PreconditionerFactory<Operator, Comm>::create(const Operator& op,
                                              const PropertyTree& prm,
                                              const std::function<Vector()>& weightsCalculator,
                                              std::size_t pressureIndex)
{
    return instance().doCreate(op, prm, weightsCalculator, pressureIndex);
}

template <class Operator, class Comm>
typename PreconditionerFactory<Operator, Comm>::PrecPtr
PreconditionerFactory<Operator, Comm>::create(const Operator& op,
                                              const PropertyTree& prm,
                                              const std::function<Vector()>& weightsCalculator,
                                              const Comm& comm,
                                              std::size_t pressureIndex)
{
    return instance().doCreate(op, prm, weightsCalculator, pressureIndex, comm);
}


template <class Operator, class Comm>
typename PreconditionerFactory<Operator, Comm>::PrecPtr
PreconditionerFactory<Operator, Comm>::create(const Operator& op,
                                              const PropertyTree& prm,
                                              const Comm& comm,
                                              std::size_t pressureIndex)
{
    return instance().doCreate(op, prm, std::function<Vector()>(), pressureIndex, comm);
}

template <class Operator, class Comm>
void
PreconditionerFactory<Operator, Comm>::addCreator(const std::string& type, Creator creator)
{
    instance().doAddCreator(type, creator);
}

template <class Operator, class Comm>
void
PreconditionerFactory<Operator, Comm>::addCreator(const std::string& type, ParCreator creator)
{
    instance().doAddCreator(type, creator);
}

using CommSeq = Dune::Amg::SequentialInformation;

template<class Scalar, int Dim>
using OpFSeq = Dune::MatrixAdapter<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, Dim, Dim>>,
                                   Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>,
                                   Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>>;
template<class Scalar, int Dim>
using OpBSeq = Dune::MatrixAdapter<Dune::BCRSMatrix<MatrixBlock<Scalar, Dim, Dim>>,
                                   Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>,
                                   Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>>;

template<class Scalar, int Dim, bool overlap>
using OpW = WellModelMatrixAdapter<Dune::BCRSMatrix<MatrixBlock<Scalar, Dim, Dim>>,
                                   Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>,
                                   Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>>;

template<class Scalar, int Dim, bool overlap>
using OpWG = WellModelGhostLastMatrixAdapter<Dune::BCRSMatrix<MatrixBlock<Scalar, Dim, Dim>>,
                                             Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>,
                                             Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>,
                                             overlap>;

#if HAVE_MPI
using CommPar = Dune::OwnerOverlapCopyCommunication<int, int>;

template<class Scalar, int Dim>
using OpFPar = Dune::OverlappingSchwarzOperator<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, Dim, Dim>>,
                                                Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>,
                                                Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>,
                                                CommPar>;

template<class Scalar, int Dim>
using OpBPar = Dune::OverlappingSchwarzOperator<Dune::BCRSMatrix<MatrixBlock<Scalar, Dim, Dim>>,
                                                Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>,
                                                Dune::BlockVector<Dune::FieldVector<Scalar, Dim>>,
                                                CommPar>;
template<class Scalar, int Dim>
using OpGLFPar = Opm::GhostLastMatrixAdapter<Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,Dim,Dim>>,
                                             Dune::BlockVector<Dune::FieldVector<Scalar,Dim>>,
                                             Dune::BlockVector<Dune::FieldVector<Scalar,Dim>>,
                                             CommPar>;

template<class Scalar, int Dim>
using OpGLBPar = Opm::GhostLastMatrixAdapter<Dune::BCRSMatrix<MatrixBlock<Scalar,Dim,Dim>>,
                                             Dune::BlockVector<Dune::FieldVector<Scalar,Dim>>,
                                             Dune::BlockVector<Dune::FieldVector<Scalar,Dim>>,
                                             CommPar>;

#define INSTANTIATE_PF_PAR(T,Dim)                                     \
    template class PreconditionerFactory<OpBSeq<T,Dim>, CommPar>;     \
    template class PreconditionerFactory<OpFPar<T,Dim>, CommPar>;     \
    template class PreconditionerFactory<OpBPar<T,Dim>, CommPar>;     \
    template class PreconditionerFactory<OpGLFPar<T,Dim>, CommPar>;   \
    template class PreconditionerFactory<OpGLBPar<T,Dim>, CommPar>;   \
    template class PreconditionerFactory<OpW<T,Dim, false>, CommPar>; \
    template class PreconditionerFactory<OpWG<T,Dim, true>, CommPar>; \
    template class PreconditionerFactory<OpBPar<T,Dim>, CommSeq>;     \
    template class PreconditionerFactory<OpGLBPar<T,Dim>, CommSeq>;
#endif

#define INSTANTIATE_PF_SEQ(T,Dim)                                     \
    template class PreconditionerFactory<OpFSeq<T,Dim>, CommSeq>;     \
    template class PreconditionerFactory<OpBSeq<T,Dim>, CommSeq>;     \
    template class PreconditionerFactory<OpW<T,Dim, false>, CommSeq>; \
    template class PreconditionerFactory<OpWG<T,Dim, true>, CommSeq>;

#if HAVE_MPI
#define INSTANTIATE_PF(T,Dim) \
    INSTANTIATE_PF_PAR(T,Dim) \
    INSTANTIATE_PF_SEQ(T,Dim)
#else
#define INSTANTIATE_PF(T,Dim) INSTANTIATE_PF_SEQ(T,Dim)
#endif

} // namespace Opm
