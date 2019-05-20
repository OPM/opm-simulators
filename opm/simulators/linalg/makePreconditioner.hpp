/*
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

#ifndef OPM_MAKEPRECONDITIONER_HEADER_INCLUDED
#define OPM_MAKEPRECONDITIONER_HEADER_INCLUDED

#include <opm/simulators/linalg/OwningBlockPreconditioner.hpp>
#include <opm/simulators/linalg/OwningTwoLevelPreconditioner.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
#include <opm/simulators/linalg/amgcpr.hh>

#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/fastamg.hh>
#include <dune/istl/preconditioners.hh>

#include <boost/property_tree/ptree.hpp>

namespace Dune
{

template <class OperatorType, class VectorType>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeSeqPreconditioner(const OperatorType& linearoperator, const boost::property_tree::ptree& prm)
{
    using MatrixType = typename OperatorType::matrix_type;
    auto& matrix = linearoperator.getmat();
    double w = prm.get<double>("w");
    int n = prm.get<int>("n");
    std::string precond(prm.get<std::string>("preconditioner"));
    if (precond == "ILU0") {
        return wrapPreconditioner<Dune::SeqILU0<MatrixType, VectorType, VectorType>>(matrix, w);
    } else if (precond == "ParOverILU0") {
        return wrapPreconditioner<Opm::ParallelOverlappingILU0<MatrixType, VectorType, VectorType>>(
            matrix, n, w, Opm::MILU_VARIANT::ILU);
    } else if (precond == "Jac") {
        return wrapPreconditioner<Dune::SeqJac<MatrixType, VectorType, VectorType>>(matrix, n, w);
    } else if (precond == "GS") {
        return wrapPreconditioner<Dune::SeqGS<MatrixType, VectorType, VectorType>>(matrix, n, w);
    } else if (precond == "SOR") {
        return wrapPreconditioner<Dune::SeqSOR<MatrixType, VectorType, VectorType>>(matrix, n, w);
    } else if (precond == "SSOR") {
        return wrapPreconditioner<Dune::SeqSSOR<MatrixType, VectorType, VectorType>>(matrix, n, w);
    } else if (precond == "ILUn") {
        return wrapPreconditioner<Dune::SeqILUn<MatrixType, VectorType, VectorType>>(matrix, n, w);
    } else {
        std::string msg("No such preconditioner : ");
        msg += precond;
        throw std::runtime_error(msg);
    }
}

template <class OperatorType, class VectorType, class Comm>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeParPreconditioner(const OperatorType& linearoperator, const boost::property_tree::ptree& prm, const Comm& comm)
{
    auto& matrix = linearoperator.getmat();
    using MatrixType = typename OperatorType::matrix_type;
    double w = prm.get<double>("w");
    int n = prm.get<int>("n");
    std::string precond(prm.get<std::string>("preconditioner"));
    if (precond == "ILU0") {
        return wrapBlockPreconditioner<DummyUpdatePreconditioner<Dune::SeqILU0<MatrixType, VectorType, VectorType>>>(
            comm, matrix, w);
    } else if (precond == "ParOverILU0") {
        // Already a parallel preconditioner. Need to pass comm, but no need to wrap it in a BlockPreconditioner.
        return wrapPreconditioner<
            DummyUpdatePreconditioner<Opm::ParallelOverlappingILU0<MatrixType, VectorType, VectorType, Comm>>>(
            matrix, comm, n, w, Opm::MILU_VARIANT::ILU);
    } else if (precond == "Jac") {
        return wrapBlockPreconditioner<DummyUpdatePreconditioner<Dune::SeqJac<MatrixType, VectorType, VectorType>>>(
            comm, matrix, n, w);
    } else if (precond == "GS") {
        return wrapBlockPreconditioner<DummyUpdatePreconditioner<Dune::SeqGS<MatrixType, VectorType, VectorType>>>(
            comm, matrix, n, w);
    } else if (precond == "SOR") {
        return wrapBlockPreconditioner<DummyUpdatePreconditioner<Dune::SeqSOR<MatrixType, VectorType, VectorType>>>(
            comm, matrix, n, w);
    } else if (precond == "SSOR") {
        return wrapBlockPreconditioner<DummyUpdatePreconditioner<Dune::SeqSSOR<MatrixType, VectorType, VectorType>>>(
            comm, matrix, n, w);
    } else if (precond == "ILUn") {
        return wrapBlockPreconditioner<DummyUpdatePreconditioner<Dune::SeqILUn<MatrixType, VectorType, VectorType>>>(
            comm, matrix, n, w);
    } else {
        std::string msg("No such preconditioner : ");
        msg += precond;
        throw std::runtime_error(msg);
    }
}

template <class Smoother, class OperatorType, class VectorType>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeAmgPreconditioner(OperatorType& linearoperator, const boost::property_tree::ptree& global_prm)
{
    boost::property_tree::ptree prm = global_prm.get_child("amg");
    using MatrixType = typename OperatorType::matrix_type;
    using CriterionBase
        = Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricMatrixDependency<MatrixType, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;
    int coarsenTarget = prm.get<int>("coarsenTarget");
    int ml = prm.get<int>("maxlevel");
    Criterion criterion(15, coarsenTarget);
    criterion.setDefaultValuesIsotropic(2);
    criterion.setAlpha(prm.get<double>("alpha"));
    criterion.setBeta(prm.get<double>("beta"));
    criterion.setMaxLevel(ml);
    criterion.setSkipIsolated(false);
    criterion.setDebugLevel(prm.get<int>("verbosity"));
    if (global_prm.get<std::string>("preconditioner") == "famg") {
        Dune::Amg::Parameters parms;
        parms.setNoPreSmoothSteps(1);
        parms.setNoPostSmoothSteps(1);
        return wrapPreconditioner<Dune::Amg::FastAMG<OperatorType, VectorType>>(linearoperator, criterion, parms);
    } else {
        typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = prm.get<int>("n");
        // smootherArgs.overlap=SmootherArgs::vertex;
        // smootherArgs.overlap=SmootherArgs::none;
        // smootherArgs.overlap=SmootherArgs::aggregate;
        smootherArgs.relaxationFactor = prm.get<double>("w");
        return std::make_shared<Dune::Amg::AMGCPR<OperatorType, VectorType, Smoother>>(
            linearoperator, criterion, smootherArgs);
    }
    return std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>(nullptr);
}

template <class Smoother, class OperatorType, class VectorType, class Comm>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeParAmgPreconditioner(OperatorType& linearoperator, const boost::property_tree::ptree& global_prm, const Comm& comm)
{
    boost::property_tree::ptree prm = global_prm.get_child("amg");
    using MatrixType = typename OperatorType::matrix_type;
    using CriterionBase
        = Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricMatrixDependency<MatrixType, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;
    int coarsenTarget = prm.get<int>("coarsenTarget");
    int ml = prm.get<int>("maxlevel");
    Criterion criterion(15, coarsenTarget);
    criterion.setDefaultValuesIsotropic(2);
    criterion.setAlpha(prm.get<double>("alpha"));
    criterion.setBeta(prm.get<double>("beta"));
    criterion.setMaxLevel(ml);
    criterion.setSkipIsolated(false);
    criterion.setDebugLevel(prm.get<int>("verbosity"));
    if (global_prm.get<std::string>("preconditioner") == "famg") {
        throw std::runtime_error("The FastAMG preconditioner cannot be used in parallel.");
    } else {
        typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
        SmootherArgs smootherArgs;
        smootherArgs.iterations = prm.get<int>("n");
        // smootherArgs.overlap=SmootherArgs::vertex;
        // smootherArgs.overlap=SmootherArgs::none;
        // smootherArgs.overlap=SmootherArgs::aggregate;
        smootherArgs.relaxationFactor = prm.get<double>("w");
        return std::make_shared<Dune::Amg::AMGCPR<OperatorType, VectorType, Smoother, Comm>>(
            linearoperator, criterion, smootherArgs, comm);
    }
}





template <class OperatorType, class VectorType>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeAmgPreconditioners(OperatorType& linearoperator, const boost::property_tree::ptree& prm)
{
    using MatrixType = typename OperatorType::matrix_type;
    if (prm.get<std::string>("preconditioner") == "famg") {
        // smoother type should not be used
        return makeAmgPreconditioner<Dune::SeqILU0<MatrixType, VectorType, VectorType>, OperatorType, VectorType>(
            linearoperator, prm);
    }

    std::string precond = prm.get<std::string>("amg.smoother");
    if (precond == "ILU0") {
        return makeAmgPreconditioner<Dune::SeqILU0<MatrixType, VectorType, VectorType>, OperatorType, VectorType>(
            linearoperator, prm);
    } else if (precond == "Jac") {
        return makeAmgPreconditioner<Dune::SeqJac<MatrixType, VectorType, VectorType>, OperatorType, VectorType>(
            linearoperator, prm);
        // } else if (precond == "GS") {
        //   return makeAmgPreconditioner<
        // 	Dune::SeqGS<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(linearoperator, prm);
    } else if (precond == "SOR") {
        return makeAmgPreconditioner<Dune::SeqSOR<MatrixType, VectorType, VectorType>, OperatorType, VectorType>(
            linearoperator, prm);
    } else if (precond == "SSOR") {
        return makeAmgPreconditioner<Dune::SeqSSOR<MatrixType, VectorType, VectorType>, OperatorType, VectorType>(
            linearoperator, prm);
    } else if (precond == "ILUn") {
        return makeAmgPreconditioner<Dune::SeqILUn<MatrixType, VectorType, VectorType>, OperatorType, VectorType>(
            linearoperator, prm);
    } else {
        std::string msg("No such sequential preconditioner : ");
        msg += precond;
        throw std::runtime_error(msg);
    }
}


template <class OperatorType, class VectorType, class Comm>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeParAmgPreconditioners(OperatorType& linearoperator, const boost::property_tree::ptree& prm, const Comm& comm)
{
    if (prm.get<std::string>("preconditioner") == "famg") {
        throw std::runtime_error("The FastAMG preconditioner cannot be used in parallel.");
    }

    using MatrixType = typename OperatorType::matrix_type;
    std::string precond = prm.get<std::string>("amg.smoother");
    if (precond == "ILU0") {
        using SmootherType = Opm::ParallelOverlappingILU0<MatrixType, VectorType, VectorType, Comm>;
        return makeParAmgPreconditioner<SmootherType, OperatorType, VectorType, Comm>(linearoperator, prm, comm);
        /*
        return makeParAmgPreconditioner<Dune::SeqILU0<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm, comm);
    } else if (precond == "Jac") {
        return makeParAmgPreconditioner<Dune::SeqJac<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm, comm);
        // } else if (precond == "GS") {
        //   return makeParAmgPreconditioner<
        // 	Dune::SeqGS<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(linearoperator, prm, comm);
    } else if (precond == "SOR") {
        return makeParAmgPreconditioner<Dune::SeqSOR<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm, comm);
    } else if (precond == "SSOR") {
        return makeParAmgPreconditioner<Dune::SeqSSOR<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm, comm);
    } else if (precond == "ILUn") {
        return makeParAmgPreconditioner<Dune::SeqILUn<MatrixType, VectorType, VectorType>, MatrixType, VectorType>(
                linearoperator, prm, comm);
        */
    } else {
        std::string msg("No such parallel preconditioner : ");
        msg += precond;
        throw std::runtime_error(msg);
    }
}



template <class OperatorType, class VectorType>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeTwoLevelPreconditioner(OperatorType& linearoperator, const boost::property_tree::ptree& global_prm)
{
    boost::property_tree::ptree prm = global_prm.get_child("cpr");
    if (global_prm.get<std::string>("preconditioner") == "cpr") {
        return std::make_shared<OwningTwoLevelPreconditioner<OperatorType, VectorType, false>>(linearoperator, prm);
    } else if (global_prm.get<std::string>("preconditioner") == "cprt") {
        return std::make_shared<OwningTwoLevelPreconditioner<OperatorType, VectorType, true>>(linearoperator, prm);
    } else {
        std::string msg("Wrong cpr Should not happen");
        throw std::runtime_error(msg);
    }
}

template <class OperatorType, class VectorType, class Comm>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makeParTwoLevelPreconditioner(OperatorType& linearoperator,
                              const boost::property_tree::ptree& global_prm,
                              const Comm& comm)
{
    boost::property_tree::ptree prm = global_prm.get_child("cpr");
    if (global_prm.get<std::string>("preconditioner") == "cpr") {
        return std::make_shared<OwningTwoLevelPreconditioner<OperatorType, VectorType, false, Comm>>(
            linearoperator, prm, comm);
    } else if (global_prm.get<std::string>("preconditioner") == "cprt") {
        return std::make_shared<OwningTwoLevelPreconditioner<OperatorType, VectorType, true, Comm>>(
            linearoperator, prm, comm);
    } else {
        std::string msg("Wrong cpr Should not happen");
        throw std::runtime_error(msg);
    }
}

template <class OperatorType, class VectorType>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makePreconditioner(OperatorType& linearoperator, const boost::property_tree::ptree& prm)
{
    if ((prm.get<std::string>("preconditioner") == "famg") or (prm.get<std::string>("preconditioner") == "amg")) {
        return makeAmgPreconditioners<OperatorType, VectorType>(linearoperator, prm);
    } else if ((prm.get<std::string>("preconditioner") == "cpr")
               or (prm.get<std::string>("preconditioner") == "cprt")) {
        return makeTwoLevelPreconditioner<OperatorType, VectorType>(linearoperator, prm);
    } else {
        return makeSeqPreconditioner<OperatorType, VectorType>(linearoperator, prm);
    }
}

template <class OperatorType, class VectorType, class Comm>
std::shared_ptr<Dune::PreconditionerWithUpdate<VectorType, VectorType>>
makePreconditioner(OperatorType& linearoperator, const boost::property_tree::ptree& prm, const Comm& comm)
{
    if ((prm.get<std::string>("preconditioner") == "famg") or (prm.get<std::string>("preconditioner") == "amg")) {
        return makeParAmgPreconditioners<OperatorType, VectorType, Comm>(linearoperator, prm, comm);
    } else if ((prm.get<std::string>("preconditioner") == "cpr")
               or (prm.get<std::string>("preconditioner") == "cprt")) {
        return makeParTwoLevelPreconditioner<OperatorType, VectorType, Comm>(linearoperator, prm, comm);
    } else {
        return makeParPreconditioner<OperatorType, VectorType, Comm>(linearoperator, prm, comm);
    }
}


} // namespace Dune


#endif // OPM_MAKEPRECONDITIONER_HEADER_INCLUDED
