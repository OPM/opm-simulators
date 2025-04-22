// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTL_GPUTWOLEVELMETHODCPR_HPP
#define DUNE_ISTL_GPUTWOLEVELMETHODCPR_HPP

// NOTE: This file is a modified version of twolevelmethodcpr.hh for GPU acceleration.


// Include the original twolevelmethodcpr.hh to use LevelTransferPolicyCpr
#include <opm/simulators/linalg/twolevelmethodcpr.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/galerkin.hh>
#include <dune/istl/solver.hh>

#include <dune/common/unused.hh>
#include <dune/common/version.hh>

#include <cstddef>
#include <memory>

/**
 * @addtogroup ISTL_PAAMG
 * @{
 * @file
 * @author Markus Blatt
 * @brief Algebraic twolevel methods for GPU.
 */
namespace Dune
{
namespace Amg
{

/**
 * @brief A policy class for solving the coarse level system using one step of AMG.
 * @tparam O The type of the linear operator used.
 * @tparam S The type of the smoother used in AMG.
 * @tparam C The type of the crition used for the aggregation within AMG.
 */
template<class O, class S, class C>
class GpuOneStepAMGCoarseSolverPolicyCpr
{
public:
  /** @brief The type of the linear operator used. */
  typedef O Operator;
  /** @brief The type of the range and domain of the operator. */
  typedef typename O::range_type X;
  /** @brief The type of the crition used for the aggregation within AMG.*/
  typedef C Criterion;
  /** @brief The type of the smoother used in AMG. */
  typedef S Smoother;
  /** @brief The type of the arguments used for constructing the smoother. */
  typedef typename Dune::Amg::SmootherTraits<S>::Arguments SmootherArgs;
  /** @brief The type of the AMG construct on the coarse level.*/
  typedef AMG<Operator,X,Smoother> AMGType;
  /**
   * @brief Constructs the coarse solver policy.
   * @param args The arguments used for constructing the smoother.
   * @param c The crition used for the aggregation within AMG.
   */
  GpuOneStepAMGCoarseSolverPolicyCpr(const SmootherArgs& args, const Criterion& c)
    : smootherArgs_(args), criterion_(c)
  {}
  /** @brief Copy constructor. */
  GpuOneStepAMGCoarseSolverPolicyCpr(const GpuOneStepAMGCoarseSolverPolicyCpr& other)
  : coarseOperator_(other.coarseOperator_), smootherArgs_(other.smootherArgs_),
    criterion_(other.criterion_)
  {}
private:
  /**
   * @brief A wrapper that makes an inverse operator out of AMG.
   *
   * The operator will use one step of AMG to approximately solve
   * the coarse level system.
   */
  struct AMGInverseOperator : public InverseOperator<X,X>
  {
    AMGInverseOperator(const typename AMGType::Operator& op,
                       const Criterion& crit,
                       const typename AMGType::SmootherArgs& args)
      : amg_(op, crit,args), first_(true)
    {}

    void apply(X& x, X& b, [[maybe_unused]] double reduction, [[maybe_unused]] InverseOperatorResult& res)
    {
      if(first_)
      {
        amg_.pre(x,b);
        first_=false;
        x_=x;
      }
      amg_.apply(x,b);
    }

    void apply(X& x, X& b, InverseOperatorResult& res)
    {
      return apply(x,b,1e-8,res);
    }

    virtual SolverCategory::Category category() const
    {
      return amg_.category();
    }

    ~AMGInverseOperator()
    {
      if(!first_)
        amg_.post(x_);
    }
    AMGInverseOperator(const AMGInverseOperator& other)
    : x_(other.x_), amg_(other.amg_), first_(other.first_)
    {
    }
  private:
    X x_;
    AMGType amg_;
    bool first_;
  };

public:
  /** @brief The type of solver constructed for the coarse level. */
  typedef AMGInverseOperator CoarseLevelSolver;

  /**
    * @brief Constructs a coarse level solver.
    *
    * @param transferPolicy The policy describing the transfer between levels.
    * @return A pointer to the constructed coarse level solver.
    * @tparam P The type of the level transfer policy.
    */
  template<class P>
  CoarseLevelSolver* createCoarseLevelSolver(P& transferPolicy)
  {
    coarseOperator_=transferPolicy.getCoarseLevelOperator();
    AMGInverseOperator* inv = new AMGInverseOperator(*coarseOperator_,
                                                     criterion_,
                                                     smootherArgs_);

    return inv; //std::shared_ptr<InverseOperator<X,X> >(inv);

  }

private:
  /** @brief The coarse level operator. */
  std::shared_ptr<Operator> coarseOperator_;
  /** @brief The arguments used to construct the smoother. */
  SmootherArgs smootherArgs_;
  /** @brief The coarsening criterion. */
  Criterion criterion_;
};

/**
 * @tparam FO The type of the fine level linear operator.
 * @tparam CSP The type of the coarse level solver policy.
 * @tparam S The type of the fine level smoother used.
 */
template<class FO, class CSP, class S>
class GpuTwoLevelMethodCpr :
    public Preconditioner<typename FO::domain_type, typename FO::range_type>
{
public:
  /** @brief The type of the policy for constructing the coarse level solver. */
  typedef CSP CoarseLevelSolverPolicy;
  /** @brief The type of the coarse level solver. */
  typedef typename CoarseLevelSolverPolicy::CoarseLevelSolver CoarseLevelSolver;
  /**
   * @brief The linear operator of the finel level system. Has to be
   * derived from AssembledLinearOperator.
   */
  typedef FO FineOperatorType;
  /**
   * @brief The type of the range of the fine level operator.
   */
  typedef typename FineOperatorType::range_type FineRangeType;
  /**
   * @brief The type of the domain of the fine level operator.
   */
  typedef typename FineOperatorType::domain_type FineDomainType;
  /**
   * @brief The linear operator of the finel level system. Has to be
   * derived from AssembledLinearOperator.
   */
  typedef typename CSP::Operator CoarseOperatorType;
  /**
   * @brief The type of the range of the coarse level operator.
   */
  typedef typename CoarseOperatorType::range_type CoarseRangeType;
  /**
   * @brief The type of the domain of the coarse level operator.
   */
  typedef typename CoarseOperatorType::domain_type CoarseDomainType;
  /**
   * @brief The type of the fine level smoother.
   */
  typedef S SmootherType;
  
  /**
   * @brief Constructs a two level method.
   *
   * @tparam CoarseSolverPolicy The policy for constructing the coarse
   * solver, e.g. OneStepAMGCoarseSolverPolicy
   * @param op The fine level operator.
   * @param smoother The fine level smoother.
   * @param policy The level transfer policy.
   * @param coarsePolicy The policy for constructing the coarse level solver.
   * @param preSteps The number of smoothing steps to apply before the coarse
   * level correction.
   * @param preSteps The number of smoothing steps to apply after the coarse
   * level correction.
   */
  GpuTwoLevelMethodCpr(const FineOperatorType& op,
                    std::shared_ptr<SmootherType> smoother,
                    const Dune::Amg::LevelTransferPolicyCpr<FineOperatorType,
                                                 CoarseOperatorType>& policy,
                    CoarseLevelSolverPolicy& coarsePolicy,
                    std::size_t preSteps=1, std::size_t postSteps=1)
    : operator_(&op), smoother_(smoother),
      preSteps_(preSteps), postSteps_(postSteps)
  {
    policy_ = policy.clone();
    policy_->createCoarseLevelSystem(*operator_);
    coarseSolver_=coarsePolicy.createCoarseLevelSolver(*policy_);
  }

  GpuTwoLevelMethodCpr(const GpuTwoLevelMethodCpr& other)
  : operator_(other.operator_), coarseSolver_(new CoarseLevelSolver(*other.coarseSolver_)),
    smoother_(other.smoother_), policy_(other.policy_->clone()),
    preSteps_(other.preSteps_), postSteps_(other.postSteps_)
  {}

  ~GpuTwoLevelMethodCpr()
  {
    // Each instance has its own policy.
    delete policy_;
    delete coarseSolver_;
  }

  void updatePreconditioner(FineOperatorType& /* op */,
                            std::shared_ptr<SmootherType> smoother,
                            CoarseLevelSolverPolicy& coarsePolicy)
  {
    updatePreconditioner(smoother, coarsePolicy);
  }

  void updatePreconditioner(std::shared_ptr<SmootherType> smoother,
                            CoarseLevelSolverPolicy& coarsePolicy)
  {
    smoother_ = smoother;
    // We assume the matrix has the same sparsity structure and memory locations throughout, only the values are changed.
    smoother_->update();
    if (coarseSolver_) {
      policy_->calculateCoarseEntries(*operator_);
      coarsePolicy.setCoarseOperator(*policy_);
      coarseSolver_->updatePreconditioner();
    } else {
      // we should probably not be heere
      policy_->createCoarseLevelSystem(*operator_);
      coarseSolver_ = coarsePolicy.createCoarseLevelSolver(*policy_);
    }
  }

  void pre(FineDomainType& x, FineRangeType& b)
  {
    smoother_->pre(x,b);
  }

  void post([[maybe_unused]] FineDomainType& x)
  {
  }

  bool hasPerfectUpdate() const
  {
    // The two-level method has perfect update if both the finesmoother and coarse solver do.
    return smoother_->hasPerfectUpdate() && coarseSolver_->hasPerfectUpdate();
  }

  void apply(FineDomainType& v, const FineRangeType& d)
  {
    FineDomainType u(v);
    FineRangeType rhs(d);
    LevelContext context;
    SequentialInformation info;
    context.pinfo=&info;
    context.lhs=&u;
    context.update=&v;
    context.smoother=smoother_;
    context.rhs=&rhs;
    context.matrix=operator_;
    // Presmoothing
    presmooth(context, preSteps_);
    //Coarse grid correction
    policy_->moveToCoarseLevel(*context.rhs);
    InverseOperatorResult res;
    coarseSolver_->apply(policy_->getCoarseLevelLhs(), policy_->getCoarseLevelRhs(), res);
    *context.lhs=0;
    policy_->moveToFineLevel(*context.lhs);
    *context.update += *context.lhs;
    // Postsmoothing
    postsmooth(context, postSteps_);

  }
//   //! Category of the preconditioner (see SolverCategory::Category)
    virtual SolverCategory::Category category() const
    {
      return SolverCategory::sequential;
    }

private:
  /**
   * @brief Struct containing the level information.
   */
  struct LevelContext
  {
    /** @brief The type of the smoother used. */
    typedef S SmootherType;
    /** @brief A pointer to the smoother. */
    std::shared_ptr<SmootherType> smoother;
    /** @brief The left hand side passed to the and returned by the smoother. */
    FineDomainType* lhs;
    /*
     * @brief The right hand side holding the current residual.
     *
     * This is passed to the smoother as the right hand side.
     */
    FineRangeType* rhs;
    /**
     * @brief The total update calculated by the preconditioner.
     *
     * I.e. all update from smoothing and coarse grid correction summed up.
     */
    FineDomainType* update;
    /** @parallel information */
    SequentialInformation* pinfo;
    /**
     * @brief The matrix that we are solving.
     *
     * Needed to update the residual.
     */
    const FineOperatorType* matrix;
  };
  const FineOperatorType* operator_;
  /** @brief The coarse level solver. */
  CoarseLevelSolver* coarseSolver_;
  /** @brief The fine level smoother. */
  std::shared_ptr<S> smoother_;
  /** @brief Policy for prolongation, restriction, and coarse level system creation. */
  Dune::Amg::LevelTransferPolicyCpr<FO,typename CSP::Operator>* policy_;
  /** @brief The number of presmoothing steps to apply. */
  std::size_t preSteps_;
  /** @brief The number of postsmoothing steps to apply. */
  std::size_t postSteps_;
};
}// end namespace Amg
}// end namespace Dune

//#endif // !DUNE_ISTL_GPUTWOLEVELMETHODCPR_HPP

/** @} */
#endif // DUNE_ISTL_GPUTWOLEVELMETHODCPR_HPP
