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

#ifndef OPM_PRESSURE_SOLVER_POLICY_HEADER_INCLUDED
#define OPM_PRESSURE_SOLVER_POLICY_HEADER_INCLUDED

#include <opm/simulators/linalg/PressureTransferPolicy.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <dune/istl/solver.hh>
#include <dune/istl/owneroverlapcopy.hh>

namespace Dune
{
namespace Amg
{
    template <class OperatorType, class Solver, class LevelTransferPolicy>
    class PressureSolverPolicy
    {
    public:
        /** @brief The type of the linear operator used. */
        using Operator = OperatorType;
        /**
         * @brief Constructs the coarse solver policy.
         * @param prm Parameter tree specifying the solver details.
         */
        explicit PressureSolverPolicy(const Opm::PropertyTree& prm)
            : prm_(prm)
        {
        }

    private:
        using X = typename Operator::range_type;
        /**
         * @brief A wrapper that makes an inverse operator out of AMG.
         *
         * The operator will use one step of AMG to approximately solve
         * the coarse level system.
         */
        struct PressureInverseOperator : public Dune::InverseOperator<X, X>
        {
#if HAVE_MPI
            template <typename GlobalIndex, typename LocalIndex>
            PressureInverseOperator(Operator& op,
                                    const Opm::PropertyTree& prm,
                                    const Dune::OwnerOverlapCopyCommunication<GlobalIndex, LocalIndex>& comm)
                : linsolver_()
            {
                assert(op.category() == Dune::SolverCategory::overlapping);
                // Assuming that we do not use Cpr as Pressure solver and use hard
                // coded pressure index that might be wrong but should be unused.
                linsolver_ = std::make_unique<Solver>(op, comm, prm, std::function<X()>(),
                                                      /* pressureIndex = */ 1);
            }
#endif // HAVE_MPI

            PressureInverseOperator(Operator& op,
                                    const Opm::PropertyTree& prm,
                                    const SequentialInformation&)
                : linsolver_()
            {
                assert(op.category() != Dune::SolverCategory::overlapping);
                // Assuming that we do not use Cpr as Pressure solver and use hard
                // coded pressure index that might be wrong but should be unused.
                linsolver_ = std::make_unique<Solver>(op, prm, std::function<X()>(),
                                                      /* pressureIndex = */ 1);
            }


            Dune::SolverCategory::Category category() const override
            {
                return linsolver_->category();
            }

            void apply(X& x, X& b, double reduction, Dune::InverseOperatorResult& res) override
            {
                linsolver_->apply(x, b, reduction, res);
            }

            void apply(X& x, X& b, Dune::InverseOperatorResult& res) override
            {
                linsolver_->apply(x, b, res);
            }

            void updatePreconditioner()
            {
                linsolver_->preconditioner().update();
            }

        private:
            std::unique_ptr<Solver> linsolver_;
        };

    public:
        /** @brief The type of solver constructed for the coarse level. */
        using CoarseLevelSolver = PressureInverseOperator;

        /**
         * @brief Constructs a coarse level solver.
         *
         * @param transferPolicy The policy describing the transfer between levels.
         * @return A pointer to the constructed coarse level solver.
         */
        template <class LTP>
        void setCoarseOperator(LTP& transferPolicy)
        {
            coarseOperator_ = transferPolicy.getCoarseLevelOperator();
        }
        template <class LTP>
        CoarseLevelSolver* createCoarseLevelSolver(LTP& transferPolicy)
        {
            coarseOperator_ = transferPolicy.getCoarseLevelOperator();
            auto& tp = dynamic_cast<LevelTransferPolicy&>(transferPolicy); // TODO: make this unnecessary.
            PressureInverseOperator* inv
                = new PressureInverseOperator(*coarseOperator_, prm_, tp.getCoarseLevelCommunication());
            return inv;
        }

    private:
        /** @brief The coarse level operator. */
        std::shared_ptr<Operator> coarseOperator_;
        Opm::PropertyTree prm_;
    };
} // namespace Amg
} // namespace Dune


#endif // OPM_PRESSURE_SOLVER_POLICY_HEADER_INCLUDED
