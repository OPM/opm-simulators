/*
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services

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
#ifndef OPM_AMGCPR_HEADER_INCLUDED
#define OPM_AMGCPR_HEADER_INCLUDED

#include <opm/simulators/linalg/twolevelmethodcpr.hh>
#include <ewoms/linear/matrixblock.hh>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/CPRPreconditioner.hpp>
#include <opm/simulators/linalg/amgcpr.hh>
#include <dune/istl/paamg/twolevelmethod.hh>
#include <dune/istl/paamg/aggregates.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Opm
{

    /**
     * \brief An algebraic twolevel or multigrid approach for solving blackoil (supports CPR with and without AMG)
     *
     * This preconditioner first decouples the component used for coarsening using a simple scaling
     * approach (e.g. Scheichl, Masson 2013,\see scaleMatrixDRS). Then it constructs the
     * coarse level system. The coupling is defined by the weights corresponding to the element located at
     * (COMPONENT_INDEX, VARIABLE_INDEX) in the block matrix. Then the coarse level system is constructed
     * either by extracting these elements, or by applying aggregation to them directly. This coarse level
     * can be solved either by AMG or by ILU. The preconditioner is configured using CPRParameter.
     * \tparam O The type of the operator (encapsulating a BCRSMatrix).
     * \tparam S The type of the smoother.
     * \tparam C The type of coarsening criterion to use.
     * \tparam P The type of the class describing the parallelization.
     * \tparam COMPONENT_INDEX The index of the component to use for coarsening (usually water).
     * \tparam VARIABLE_INDEX The index of the variable to use for coarsening (usually pressure).
     */
    template<typename O, typename S, typename SC, typename C,
             typename P, std::size_t COMPONENT_INDEX, std::size_t VARIABLE_INDEX>
    class BlackoilAmgCpr
        : public Dune::Preconditioner<typename O::domain_type, typename O::range_type>
    {
    public:
        /** \brief The type of the operator (encapsulating a BCRSMatrix). */
        using Operator = O;
        /** \brief The type of coarsening criterion to use. */
        using Criterion = C;
        /** \brief The type of the class describing the parallelization. */
        using Communication = P;
        /** \brief The type of the smoother. */
        using Smoother = S;
        /** \brief The type of the smoother arguments for construction. */
        using SmootherArgs   = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;

    protected:
        using Matrix = typename Operator::matrix_type;
        using CoarseOperator = typename Detail::ScalarType<Operator>::value;
        using CoarseSmoother = typename Detail::ScalarType<SC>::value;
        using FineCriterion  =
            typename Detail::OneComponentCriterionType<Criterion,COMPONENT_INDEX, VARIABLE_INDEX>::value;
        using CoarseCriterion =  typename Detail::ScalarType<Criterion>::value;
        using LevelTransferPolicy =
            OneComponentAggregationLevelTransferPolicy<Operator,
                                                       FineCriterion,
                                                       Communication,
                                                       COMPONENT_INDEX,
                                                       VARIABLE_INDEX>;
        using CoarseSolverPolicy   =
            Detail::OneStepAMGCoarseSolverPolicy<CoarseOperator,
                                                 CoarseSmoother,
                                                 CoarseCriterion,
                                                 LevelTransferPolicy>;
        using TwoLevelMethod =
            Dune::Amg::TwoLevelMethodCpr<Operator,
                                         CoarseSolverPolicy,
                                         Smoother>;
    public:
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
        Dune::SolverCategory::Category category() const override
        {
            return std::is_same<Communication, Dune::Amg::SequentialInformation>::value ?
                Dune::SolverCategory::sequential : Dune::SolverCategory::overlapping;
        }
#else
        // define the category
        enum {
            //! \brief The category the precondtioner is part of.
            category = Operator::category
        };
#endif

        /**
         * \brief Constructor.
         * \param param The parameters used for configuring the solver.
         * \param fineOperator The operator of the fine level.
         * \param criterion The criterion describing the coarsening approach.
         * \param smargs The arguments for constructing the smoother.
         * \param comm The information about the parallelization.
         */
        BlackoilAmgCpr(const CPRParameter& param,
                       const typename TwoLevelMethod::FineDomainType& weights,
                       const Operator& fineOperator, const Criterion& criterion,
                       const SmootherArgs& smargs, const Communication& comm)
            : param_(param),
              weights_(weights),
              scaledMatrix_(Detail::scaleMatrixDRS(fineOperator, COMPONENT_INDEX, weights_, param)),
              scaledMatrixOperator_(Detail::createOperator(fineOperator, *scaledMatrix_, comm)),
              smoother_(Detail::constructSmoother<Smoother>(*scaledMatrixOperator_,
                                                            smargs, comm)),
              levelTransferPolicy_(criterion, comm, param.cpr_pressure_aggregation_),
              coarseSolverPolicy_(&param, smargs, criterion),
              twoLevelMethod_(*scaledMatrixOperator_,
                              smoother_,
                              levelTransferPolicy_,
                              coarseSolverPolicy_, 0, 1)
        {
        }

        void updatePreconditioner(const Operator& fineOperator,
                                  const SmootherArgs& smargs,
                                  const Communication& comm)
        {
            *scaledMatrix_ = *Detail::scaleMatrixDRS(fineOperator, COMPONENT_INDEX, weights_, param_);
            smoother_ = Detail::constructSmoother<Smoother>(*scaledMatrixOperator_, smargs, comm);
            twoLevelMethod_.updatePreconditioner(*scaledMatrixOperator_,
                                                 smoother_,
                                                 coarseSolverPolicy_);
        }

        void pre(typename TwoLevelMethod::FineDomainType& x,
                 typename TwoLevelMethod::FineRangeType& b) override
        {
            twoLevelMethod_.pre(x,b);
        }

        void post(typename TwoLevelMethod::FineDomainType& x) override
        {
            twoLevelMethod_.post(x);
        }

        void apply(typename TwoLevelMethod::FineDomainType& v,
                   const typename TwoLevelMethod::FineRangeType& d) override
        {
            auto scaledD = d;
            Detail::scaleVectorDRS(scaledD, COMPONENT_INDEX, param_, weights_);
            twoLevelMethod_.apply(v, scaledD);
        }

    private:
        const CPRParameter& param_;
        const typename TwoLevelMethod::FineDomainType& weights_;
        std::unique_ptr<Matrix> scaledMatrix_;
        std::unique_ptr<Operator> scaledMatrixOperator_;
        std::shared_ptr<Smoother> smoother_;
        LevelTransferPolicy levelTransferPolicy_;
        CoarseSolverPolicy coarseSolverPolicy_;
        TwoLevelMethod twoLevelMethod_;
    };

} // end namespace Opm
#endif
