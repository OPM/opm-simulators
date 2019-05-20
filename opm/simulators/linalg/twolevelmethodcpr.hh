// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTL_TWOLEVELMETHODCPR_HH
#define DUNE_ISTL_TWOLEVELMETHODCPR_HH

// NOTE: This file is a modified version of dune/istl/paamg/twolevelmethod.hh from
// dune-istl release 2.6.0. Modifications have been kept as minimal as possible.

#include <tuple>

#include<dune/istl/operators.hh>
//#include "amg.hh"
//#include"galerkin.hh"
#include<dune/istl/paamg/amg.hh>
#include<dune/istl/paamg/galerkin.hh>
#include<dune/istl/solver.hh>

#include<dune/common/unused.hh>
#include<dune/common/version.hh>

/**
 * @addtogroup ISTL_PAAMG
 * @{
 * @file
 * @author Markus Blatt
 * @brief Algebraic twolevel methods.
 */
namespace Dune
{
namespace Amg
{

/**
 * @brief Abstract base class for transfer between levels and creation
 * of the coarse level system.
 *
 * @tparam FO The type of the linear operator of the finel level system. Has to be
 * derived from AssembledLinearOperator.
 * @tparam CO The type of the linear operator of the coarse level system. Has to be
 * derived from AssembledLinearOperator.
 */
template<class FO, class CO>
class LevelTransferPolicyCpr
{
public:
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
  typedef CO CoarseOperatorType;
  /**
   * @brief The type of the range of the coarse level operator.
   */
  typedef typename CoarseOperatorType::range_type CoarseRangeType;
  /**
   * @brief The type of the domain of the coarse level operator.
   */
  typedef typename CoarseOperatorType::domain_type CoarseDomainType;
  /**
   * @brief Get the coarse level operator.
   * @return A shared pointer to the coarse level system.
   */
  std::shared_ptr<CoarseOperatorType>& getCoarseLevelOperator()
  {
    return operator_;
  }
  /**
   * @brief Get the coarse level right hand side.
   * @return The coarse level right hand side.
   */
  CoarseRangeType& getCoarseLevelRhs()
  {
    return rhs_;
  }

  /**
   * @brief Get the coarse level left hand side.
   * @return The coarse level leftt hand side.
   */
  CoarseDomainType& getCoarseLevelLhs()
  {
    return lhs_;
  }
  /**
   * @brief Transfers the data to the coarse level.
   *
   * Restricts the residual to the right hand side of the
   * coarse level system and initialies the left hand side
   * of the coarse level system. These can afterwards be accessed
   * usinf getCoarseLevelRhs() and getCoarseLevelLhs().
   * @param fineDefect The current residual of the fine level system.
   */
  virtual void moveToCoarseLevel(const FineRangeType& fineRhs)=0;
  /**
   * @brief Updates the fine level linear system after the correction
   * of the coarse levels system.
   *
   * After returning from this function the coarse level correction
   * will have been added to fine level system.
   * @param[inout] fineLhs The left hand side of the fine level to update
   * with the coarse level correction.
   */
  virtual void moveToFineLevel(FineDomainType& fineLhs)=0;
  /**
   * @brief Algebraically creates the coarse level system.
   *
   * After returning from this function the coarse level operator
   * can be accessed using getCoarseLevelOperator().
   * @param fineOperator The operator of the fine level system.
   */
  virtual void createCoarseLevelSystem(const FineOperatorType& fineOperator)=0;

  /**
   * @brief ???.
   */
  virtual void calculateCoarseEntries(const FineOperatorType& fineOperator) = 0;

  /** @brief Clone the current object. */
  virtual LevelTransferPolicyCpr* clone() const =0;

  /** @brief Destructor. */
  virtual ~LevelTransferPolicyCpr(){}

  protected:
  /** @brief The coarse level rhs. */
  CoarseRangeType rhs_;
  /** @brief The coarse level lhs. */
  CoarseDomainType lhs_;
  /** @brief the coarse level linear operator. */
  std::shared_ptr<CoarseOperatorType> operator_;
};

/**
 * @brief A LeveTransferPolicy that used aggregation to construct the coarse level system.
 * @tparam O The type of the fine and coarse level operator.
 * @tparam C The criterion that describes the aggregation procedure.
 */
template<class O, class C>
class AggregationLevelTransferPolicyCpr
  : public LevelTransferPolicyCpr<O,O>
{
  typedef Dune::Amg::AggregatesMap<typename O::matrix_type::size_type> AggregatesMap;
public:
  typedef LevelTransferPolicyCpr<O,O> FatherType;
  typedef C Criterion;
  typedef SequentialInformation ParallelInformation;

  AggregationLevelTransferPolicyCpr(const Criterion& crit)
  : criterion_(crit)
  {}

  void createCoarseLevelSystem(const O& fineOperator)
  {
    prolongDamp_ = criterion_.getProlongationDampingFactor();
    GalerkinProduct<ParallelInformation> productBuilder;
    typedef typename Dune::Amg::MatrixGraph<const typename O::matrix_type> MatrixGraph;
    typedef typename Dune::Amg::PropertiesGraph<MatrixGraph,Dune::Amg::VertexProperties,
      Dune::Amg::EdgeProperties,Dune::IdentityMap,Dune::IdentityMap> PropertiesGraph;
    MatrixGraph mg(fineOperator.getmat());
    PropertiesGraph pg(mg,Dune::IdentityMap(),Dune::IdentityMap());
    typedef NegateSet<typename ParallelInformation::OwnerSet> OverlapFlags;

    aggregatesMap_.reset(new AggregatesMap(pg.maxVertex()+1));

    int noAggregates, isoAggregates, oneAggregates, skippedAggregates;

    std::tie(noAggregates, isoAggregates, oneAggregates, skippedAggregates) =
       aggregatesMap_->buildAggregates(fineOperator.getmat(), pg, criterion_, true);
     std::cout<<"no aggregates="<<noAggregates<<" iso="<<isoAggregates<<" one="<<oneAggregates<<" skipped="<<skippedAggregates<<std::endl;
    // misuse coarsener to renumber aggregates
    Dune::Amg::IndicesCoarsener<Dune::Amg::SequentialInformation,int> renumberer;
    typedef std::vector<bool>::iterator Iterator;
    typedef Dune::IteratorPropertyMap<Iterator, Dune::IdentityMap> VisitedMap;
    std::vector<bool> excluded(fineOperator.getmat().N(), false);
    VisitedMap vm(excluded.begin(), Dune::IdentityMap());
    ParallelInformation pinfo;
    std::size_t aggregates = renumberer.coarsen(pinfo, pg, vm,
                                                *aggregatesMap_, pinfo,
                                                noAggregates);
    std::vector<bool>& visited=excluded;

    typedef std::vector<bool>::iterator Iterator;

    for(Iterator iter= visited.begin(), end=visited.end();
        iter != end; ++iter)
          *iter=false;
    matrix_.reset(productBuilder.build(mg, vm,
                                       SequentialInformation(),
                                       *aggregatesMap_,
                                       aggregates,
                                       OverlapFlags()));
    productBuilder.calculate(fineOperator.getmat(), *aggregatesMap_, *matrix_, pinfo, OverlapFlags());
    this->lhs_.resize(this->matrix_->M());
    this->rhs_.resize(this->matrix_->N());
    this->operator_.reset(new O(*matrix_));
  }

  void moveToCoarseLevel(const typename FatherType::FineRangeType& fineRhs)
  {
    Transfer<std::size_t,typename FatherType::FineRangeType,ParallelInformation>
      ::restrictVector(*aggregatesMap_, this->rhs_, fineRhs, ParallelInformation());
    this->lhs_=0;
  }

  void moveToFineLevel(typename FatherType::FineDomainType& fineLhs)
  {
    Transfer<std::size_t,typename FatherType::FineRangeType,ParallelInformation>
      ::prolongateVector(*aggregatesMap_, this->lhs_, fineLhs,
                         prolongDamp_, ParallelInformation());
  }

  AggregationLevelTransferPolicyCpr* clone() const
  {
    return new AggregationLevelTransferPolicyCpr(*this);
  }

private:
  typename O::matrix_type::field_type prolongDamp_;
  std::shared_ptr<AggregatesMap> aggregatesMap_;
  Criterion criterion_;
  std::shared_ptr<typename O::matrix_type> matrix_;
};

/**
 * @brief A policy class for solving the coarse level system using one step of AMG.
 * @tparam O The type of the linear operator used.
 * @tparam S The type of the smoother used in AMG.
 * @tparam C The type of the crition used for the aggregation within AMG.
 */
template<class O, class S, class C>
class OneStepAMGCoarseSolverPolicyCpr
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
  OneStepAMGCoarseSolverPolicyCpr(const SmootherArgs& args, const Criterion& c)
    : smootherArgs_(args), criterion_(c)
  {}
  /** @brief Copy constructor. */
  OneStepAMGCoarseSolverPolicyCpr(const OneStepAMGCoarseSolverPolicyCpr& other)
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

    void apply(X& x, X& b, double reduction, InverseOperatorResult& res)
    {
      DUNE_UNUSED_PARAMETER(reduction);
      DUNE_UNUSED_PARAMETER(res);
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

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
    virtual SolverCategory::Category category() const
    {
      return amg_.category();
    }
#endif
    
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
class TwoLevelMethodCpr :
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
  
  // define the category
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
  
#else  
  enum {
    //! \brief The category the preconditioner is part of.
    category=SolverCategory::sequential
  };
#endif

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
  TwoLevelMethodCpr(const FineOperatorType& op,
                    std::shared_ptr<SmootherType> smoother,
                    const LevelTransferPolicyCpr<FineOperatorType,
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

  TwoLevelMethodCpr(const TwoLevelMethodCpr& other)
  : operator_(other.operator_), coarseSolver_(new CoarseLevelSolver(*other.coarseSolver_)),
    smoother_(other.smoother_), policy_(other.policy_->clone()),
    preSteps_(other.preSteps_), postSteps_(other.postSteps_)
  {}

  ~TwoLevelMethodCpr()
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
    //assume new matrix is not reallocated the new precondition should anyway be made
    smoother_ = smoother;
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

  void post(FineDomainType& x)
  {
    DUNE_UNUSED_PARAMETER(x);
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
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)  
//   //! Category of the preconditioner (see SolverCategory::Category)
    virtual SolverCategory::Category category() const
    {
      return SolverCategory::sequential;
    }
#endif
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
  LevelTransferPolicyCpr<FO,typename CSP::Operator>* policy_;
  /** @brief The number of presmoothing steps to apply. */
  std::size_t preSteps_;
  /** @brief The number of postsmoothing steps to apply. */
  std::size_t postSteps_;
};
}// end namespace Amg
}// end namespace Dune

/** @} */
#endif
