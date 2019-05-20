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

#ifndef OPM_OWNINGTWOLEVELPRECONDITIONER_HEADER_INCLUDED
#define OPM_OWNINGTWOLEVELPRECONDITIONER_HEADER_INCLUDED

#include <opm/simulators/linalg/PressureSolverPolicy.hpp>
#include <opm/simulators/linalg/PressureTransferPolicy.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/amg.hh>

#include <boost/property_tree/ptree.hpp>

#include <fstream>

namespace Dune
{

// Circular dependency between makePreconditioner() [which can make an OwningTwoLevelPreconditioner]
// and OwningTwoLevelPreconditioner [which uses makePreconditioner() to choose the fine-level smoother]
// must be broken, accomplished by forward-declaration here.
template <class MatrixType, class VectorType>
std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>>
makePreconditioner(const Dune::MatrixAdapter<MatrixType, VectorType, VectorType>& linearoperator,
                   const boost::property_tree::ptree& prm);

// Must forward-declare FlexibleSolver as we want to use it as solver for the pressure system.
template <class MatrixTypeT, class VectorTypeT>
class FlexibleSolver;

template <class MatrixTypeT, class VectorTypeT, bool transpose = false>
class OwningTwoLevelPreconditioner : public Dune::Preconditioner<VectorTypeT, VectorTypeT>
{
public:
    using pt = boost::property_tree::ptree;
    using MatrixType = MatrixTypeT;
    using VectorType = VectorTypeT;
    using OperatorType = Dune::MatrixAdapter<MatrixType, VectorType, VectorType>;

    OwningTwoLevelPreconditioner(OperatorType& linearoperator, pt& prm)
        : finesmoother_(makePreconditioner<MatrixType, VectorType>(linearoperator, prm.get_child("finesmoother")))
        , comm_()
        , weights_(Opm::Amg::getQuasiImpesWeights<MatrixType, VectorType>(
              linearoperator.getmat(), prm.get<int>("pressure_var_index"), transpose))
        , levelTransferPolicy_(comm_, weights_, prm.get<int>("pressure_var_index"))
        , coarseSolverPolicy_(prm.get_child("coarsesolver"))
        , twolevel_method_(linearoperator,
                           finesmoother_,
                           levelTransferPolicy_,
                           coarseSolverPolicy_,
                           transpose ? 1 : 0,
                           transpose ? 0 : 1)
    {
        if (prm.get<int>("verbosity") > 10) {
            std::ofstream outfile(prm.get<std::string>("weights_filename"));
            if (!outfile) {
                throw std::runtime_error("Could not write weights");
            }
            Dune::writeMatrixMarket(weights_, outfile);
        }
    }

    virtual void pre(VectorType& x, VectorType& b) override
    {
        twolevel_method_.pre(x, b);
    }
    virtual void apply(VectorType& v, const VectorType& d) override
    {
        twolevel_method_.apply(v, d);
    }
    virtual void post(VectorType& x) override
    {
        twolevel_method_.post(x);
    }
    virtual Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    using PressureMatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using PressureVectorType = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using CoarseOperatorType = Dune::MatrixAdapter<PressureMatrixType, PressureVectorType, PressureVectorType>;
    using Communication = Dune::Amg::SequentialInformation;
    using LevelTransferPolicy = Opm::PressureTransferPolicy<OperatorType, CoarseOperatorType, Communication, transpose>;
    using CoarseSolverPolicy
        = Dune::Amg::PressureSolverPolicy<CoarseOperatorType, FlexibleSolver<PressureMatrixType, PressureVectorType>>;
    using TwoLevelMethod
        = Dune::Amg::TwoLevelMethod<OperatorType, CoarseSolverPolicy, Dune::Preconditioner<VectorType, VectorType>>;

    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> finesmoother_;
    Communication comm_;
    VectorType weights_;
    LevelTransferPolicy levelTransferPolicy_;
    CoarseSolverPolicy coarseSolverPolicy_;
    TwoLevelMethod twolevel_method_;
};

} // namespace Dune




#endif // OPM_OWNINGTWOLEVELPRECONDITIONER_HEADER_INCLUDED
