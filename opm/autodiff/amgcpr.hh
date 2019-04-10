// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_AMG_AMG_CPR_HH
#define DUNE_AMG_AMG_CPR_HH

// NOTE: This file is a modified version of dune/istl/paamg/amg.hh from
// dune-istl release 2.6.0. Modifications have been kept as minimal as possible.

#include <memory>
#include <dune/common/exceptions.hh>
#include <dune/istl/paamg/smoother.hh>
#include <dune/istl/paamg/transfer.hh>
#include <dune/istl/paamg/hierarchy.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/solvertype.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

namespace Dune
{
  namespace Amg
  {
    /**
     * @defgroup ISTL_PAAMG Parallel Algebraic Multigrid
     * @ingroup ISTL_Prec
     * @brief A Parallel Algebraic Multigrid based on Agglomeration.
     */

    /**
     * @addtogroup ISTL_PAAMG
     *
     * @{
     */

    /** @file
     * @author Markus Blatt
     * @brief The AMG preconditioner.
     */

    template<class M, class X, class S, class P, class K, class A>
    class KAMG;

    template<class T>
    class KAmgTwoGrid;

    /**
     * @brief Parallel algebraic multigrid based on agglomeration.
     *
     * \tparam M The matrix type
     * \tparam X The vector type
     * \tparam S The smoother type
     * \tparam A An allocator for X
     *
     * \todo drop the smoother template parameter and replace with dynamic construction
     */
    template<class M, class X, class S, class PI=SequentialInformation,
        class A=std::allocator<X> >
    class AMGCPR : public Preconditioner<X,X>
    {
      template<class M1, class X1, class S1, class P1, class K1, class A1>
      friend class KAMG;

      friend class KAmgTwoGrid<AMGCPR>;

    public:
      /** @brief The matrix operator type. */
      typedef M Operator;
      /**
       * @brief The type of the parallel information.
       * Either OwnerOverlapCommunication or another type
       * describing the parallel data distribution and
       * providing communication methods.
       */
      typedef PI ParallelInformation;
      /** @brief The operator hierarchy type. */
      typedef MatrixHierarchy<M, ParallelInformation, A> OperatorHierarchy;
      /** @brief The parallal data distribution hierarchy type. */
      typedef typename OperatorHierarchy::ParallelInformationHierarchy ParallelInformationHierarchy;

      /** @brief The domain type. */
      typedef X Domain;
      /** @brief The range type. */
      typedef X Range;
      /** @brief the type of the coarse solver. */
      typedef InverseOperator<X,X> CoarseSolver;
      /**
       * @brief The type of the smoother.
       *
       * One of the preconditioners implementing the Preconditioner interface.
       * Note that the smoother has to fit the ParallelInformation.*/
      typedef S Smoother;

      /** @brief The argument type for the construction of the smoother. */
      typedef typename SmootherTraits<Smoother>::Arguments SmootherArgs;

      /**
       * @brief Construct a new amg with a specific coarse solver.
       * @param matrices The already set up matix hierarchy.
       * @param coarseSolver The set up solver to use on the coarse
       * grid, must match the coarse matrix in the matrix hierarchy.
       * @param smootherArgs The  arguments needed for thesmoother to use
       * for pre and post smoothing.
       * @param parms The parameters for the AMG.
       */
      AMGCPR(const OperatorHierarchy& matrices, CoarseSolver& coarseSolver,
          const SmootherArgs& smootherArgs, const Parameters& parms);

      /**
       * @brief Construct an AMG with an inexact coarse solver based on the smoother.
       *
       * As coarse solver a preconditioned CG method with the smoother as preconditioner
       * will be used. The matrix hierarchy is built automatically.
       * @param fineOperator The operator on the fine level.
       * @param criterion The criterion describing the coarsening strategy. E. g. SymmetricCriterion
       * or UnsymmetricCriterion, and providing the parameters.
       * @param smootherArgs The arguments for constructing the smoothers.
       * @param pinfo The information about the parallel distribution of the data.
       */
      template<class C>
      AMGCPR(const Operator& fineOperator, const C& criterion,
          const SmootherArgs& smootherArgs=SmootherArgs(),
          const ParallelInformation& pinfo=ParallelInformation());

      /**
       * @brief Copy constructor.
       */
      AMGCPR(const AMGCPR& amg);

      ~AMGCPR();

      /** \copydoc Preconditioner::pre */
      void pre(Domain& x, Range& b);

      /** \copydoc Preconditioner::apply */
      void apply(Domain& v, const Range& d);

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
      //! Category of the preconditioner (see SolverCategory::Category)
      virtual SolverCategory::Category category() const
      {
        return category_;
      }
#else
      enum {
            //! \brief The category the preconditioner is part of.
            category = std::is_same<PI,Dune::Amg::SequentialInformation>::value?
            Dune::SolverCategory::sequential:Dune::SolverCategory::overlapping
        };
#endif

      /** \copydoc Preconditioner::post */
      void post(Domain& x);

      /**
       * @brief Get the aggregate number of each unknown on the coarsest level.
       * @param cont The random access container to store the numbers in.
       */
      template<class A1>
      void getCoarsestAggregateNumbers(std::vector<std::size_t,A1>& cont);

      std::size_t levels();

      std::size_t maxlevels();

      /**
       * @brief Recalculate the matrix hierarchy.
       *
       * It is assumed that the coarsening for the changed fine level
       * matrix would yield the same aggregates. In this case it suffices
       * to recalculate all the Galerkin products for the matrices of the
       * coarser levels.
       */
      void recalculateHierarchy()
      {
        matrices_->recalculateGalerkin(NegateSet<typename PI::OwnerSet>());
      }

      /**
       * @brief Update the coarse solver and the hierarchies.
       */
      template<class C>  
      void updateSolver(C& criterion, Operator& /* matrix */, const PI& pinfo);

      /**
       * @brief Check whether the coarse solver used is a direct solver.
       * @return True if the coarse level solver is a direct solver.
       */
      bool usesDirectCoarseLevelSolver() const;

    private:
      /**
       * @brief Create matrix and smoother hierarchies.
       * @param criterion The coarsening criterion.
       * @param matrix The fine level matrix operator.
       * @param pinfo The fine level parallel information.
       */
      template<class C>
      void createHierarchies(C& criterion, Operator& matrix,
                             const PI& pinfo);

      void setupCoarseSolver();

      /**
       * @brief A struct that holds the context of the current level.
       *
       * These are the iterators to the smoother, matrix, parallel information,
       * and so on needed for the computations on the current level.
       */
      struct LevelContext
      {
        typedef Smoother SmootherType;
        /**
         * @brief The iterator over the smoothers.
         */
        typename Hierarchy<Smoother,A>::Iterator smoother;
        /**
         * @brief The iterator over the matrices.
         */
        typename OperatorHierarchy::ParallelMatrixHierarchy::ConstIterator matrix;
        /**
         * @brief The iterator over the parallel information.
         */
        typename ParallelInformationHierarchy::Iterator pinfo;
        /**
         * @brief The iterator over the redistribution information.
         */
        typename OperatorHierarchy::RedistributeInfoList::const_iterator redist;
        /**
         * @brief The iterator over the aggregates maps.
         */
        typename OperatorHierarchy::AggregatesMapList::const_iterator aggregates;
        /**
         * @brief The iterator over the left hand side.
         */
        typename Hierarchy<Domain,A>::Iterator lhs;
        /**
         * @brief The iterator over the updates.
         */
        typename Hierarchy<Domain,A>::Iterator update;
        /**
         * @brief The iterator over the right hand sided.
         */
        typename Hierarchy<Range,A>::Iterator rhs;
        /**
         * @brief The level index.
         */
        std::size_t level;
      };


      /**
       * @brief Multigrid cycle on a level.
       * @param levelContext the iterators of the current level.
       */
      void mgc(LevelContext& levelContext);

      void additiveMgc();

      /**
       * @brief Move the iterators to the finer level
       * @param levelContext the iterators of the current level.
       * @param processedFineLevel Whether the process computed on
       *         fine level or not.
       */
      void moveToFineLevel(LevelContext& levelContext,bool processedFineLevel);

      /**
       * @brief Move the iterators to the coarser level.
       * @param levelContext the iterators of the current level
       */
      bool moveToCoarseLevel(LevelContext& levelContext);

      /**
       * @brief Initialize iterators over levels with fine level.
       * @param levelContext the iterators of the current level
       */
      void initIteratorsWithFineLevel(LevelContext& levelContext);

      /**  @brief The matrix we solve. */
      std::shared_ptr<OperatorHierarchy> matrices_;
      /** @brief The arguments to construct the smoother */
      SmootherArgs smootherArgs_;
      /** @brief The hierarchy of the smoothers. */
      std::shared_ptr<Hierarchy<Smoother,A> > smoothers_;
      /** @brief The solver of the coarsest level. */
      std::shared_ptr<CoarseSolver> solver_;
      /** @brief The right hand side of our problem. */
      Hierarchy<Range,A>* rhs_;
      /** @brief The left approximate solution of our problem. */
      Hierarchy<Domain,A>* lhs_;
      /** @brief The total update for the outer solver. */
      Hierarchy<Domain,A>* update_;
      /** @brief The type of the scalar product for the coarse solver. */
      using ScalarProduct = Dune::ScalarProduct<X>;
      /** @brief Scalar product on the coarse level. */
      std::shared_ptr<ScalarProduct> scalarProduct_;
      /** @brief Gamma, 1 for V-cycle and 2 for W-cycle. */
      std::size_t gamma_;
      /** @brief The number of pre and postsmoothing steps. */
      std::size_t preSteps_;
      /** @brief The number of postsmoothing steps. */
      std::size_t postSteps_;
      bool buildHierarchy_;
      bool additive;
      bool coarsesolverconverged;
      std::shared_ptr<Smoother> coarseSmoother_;
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
      /** @brief The solver category. */
      SolverCategory::Category category_;
#endif
      /** @brief The verbosity level. */
      std::size_t verbosity_;
    };

    template<class M, class X, class S, class PI, class A>
    inline AMGCPR<M,X,S,PI,A>::AMGCPR(const AMGCPR& amg)
    : matrices_(amg.matrices_), smootherArgs_(amg.smootherArgs_),
      smoothers_(amg.smoothers_), solver_(amg.solver_),
      rhs_(), lhs_(), update_(),
      scalarProduct_(amg.scalarProduct_), gamma_(amg.gamma_),
      preSteps_(amg.preSteps_), postSteps_(amg.postSteps_),
      buildHierarchy_(amg.buildHierarchy_),
      additive(amg.additive), coarsesolverconverged(amg.coarsesolverconverged),
      coarseSmoother_(amg.coarseSmoother_),
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
      category_(amg.category_),
#endif
      verbosity_(amg.verbosity_)
    {
      if(amg.rhs_)
        rhs_=new Hierarchy<Range,A>(*amg.rhs_);
      if(amg.lhs_)
        lhs_=new Hierarchy<Domain,A>(*amg.lhs_);
      if(amg.update_)
        update_=new Hierarchy<Domain,A>(*amg.update_);
    }

    template<class M, class X, class S, class PI, class A>
    AMGCPR<M,X,S,PI,A>::AMGCPR(const OperatorHierarchy& matrices, CoarseSolver& coarseSolver,
                         const SmootherArgs& smootherArgs,
                         const Parameters& parms)
      : matrices_(stackobject_to_shared_ptr(matrices)), smootherArgs_(smootherArgs),
        smoothers_(new Hierarchy<Smoother,A>), solver_(&coarseSolver),
        rhs_(), lhs_(), update_(), scalarProduct_(0),
        gamma_(parms.getGamma()), preSteps_(parms.getNoPreSmoothSteps()),
        postSteps_(parms.getNoPostSmoothSteps()), buildHierarchy_(false),
        additive(parms.getAdditive()), coarsesolverconverged(true),
        coarseSmoother_(),
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
// #warning should category be retrieved from matrices?
        category_(SolverCategory::category(*smoothers_->coarsest())),
#endif
        verbosity_(parms.debugLevel())
    {
      assert(matrices_->isBuilt());

      // build the necessary smoother hierarchies
      matrices_->coarsenSmoother(*smoothers_, smootherArgs_);
    }

    template<class M, class X, class S, class PI, class A>
    template<class C>
    AMGCPR<M,X,S,PI,A>::AMGCPR(const Operator& matrix,
                         const C& criterion,
                         const SmootherArgs& smootherArgs,
                         const PI& pinfo)
      : smootherArgs_(smootherArgs),
        smoothers_(new Hierarchy<Smoother,A>), solver_(),
        rhs_(), lhs_(), update_(), scalarProduct_(),
        gamma_(criterion.getGamma()), preSteps_(criterion.getNoPreSmoothSteps()),
        postSteps_(criterion.getNoPostSmoothSteps()), buildHierarchy_(true),
        additive(criterion.getAdditive()), coarsesolverconverged(true),
        coarseSmoother_(),
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
        category_(SolverCategory::category(pinfo)),
#endif
        verbosity_(criterion.debugLevel())
    {
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
      if(SolverCategory::category(matrix) != SolverCategory::category(pinfo))
        DUNE_THROW(InvalidSolverCategory, "Matrix and Communication must have the same SolverCategory!");
#endif
      createHierarchies(criterion, const_cast<Operator&>(matrix), pinfo);
    }


    template<class M, class X, class S, class PI, class A>
    AMGCPR<M,X,S,PI,A>::~AMGCPR()
    {
      if(buildHierarchy_) {
        if(solver_)
          solver_.reset();
        if(coarseSmoother_)
          coarseSmoother_.reset();
      }
      if(lhs_)
        delete lhs_;
      lhs_=nullptr;
      if(update_)
        delete update_;
      update_=nullptr;
      if(rhs_)
        delete rhs_;
      rhs_=nullptr;
    }

    template<class M, class X, class S, class PI, class A>
    template<class C>
    void AMGCPR<M,X,S,PI,A>::updateSolver(C& /* criterion */, Operator& /* matrix */, const PI& /* pinfo */)
    {
      Timer watch;
      smoothers_.reset(new Hierarchy<Smoother,A>);
      solver_.reset();
      coarseSmoother_.reset();
      scalarProduct_.reset();
      buildHierarchy_= true;
      coarsesolverconverged = true;
      smoothers_.reset(new Hierarchy<Smoother,A>);
      matrices_->recalculateGalerkin(NegateSet<typename PI::OwnerSet>());
      matrices_->coarsenSmoother(*smoothers_, smootherArgs_);
      setupCoarseSolver();
      if (verbosity_>0 && matrices_->parallelInformation().finest()->communicator().rank()==0) {
        std::cout << "Recalculating galerkin and coarse somothers "<< matrices_->maxlevels() << " levels "
                  << watch.elapsed() << " seconds." << std::endl;
      }
    }

    template<class M, class X, class S, class PI, class A>
    template<class C>
    void AMGCPR<M,X,S,PI,A>::createHierarchies(C& criterion, Operator& matrix,
                                            const PI& pinfo)
    {
      Timer watch;
      matrices_.reset(new OperatorHierarchy(matrix, pinfo));

      matrices_->template build<NegateSet<typename PI::OwnerSet> >(criterion);

      // build the necessary smoother hierarchies
      matrices_->coarsenSmoother(*smoothers_, smootherArgs_);
      setupCoarseSolver();
      if(verbosity_>0 && matrices_->parallelInformation().finest()->communicator().rank()==0)
        std::cout<<"Building hierarchy of "<<matrices_->maxlevels()<<" levels "
                 <<"(inclusive coarse solver) took "<<watch.elapsed()<<" seconds."<<std::endl;
    }

    template<class M, class X, class S, class PI, class A>
    void AMGCPR<M,X,S,PI,A>::setupCoarseSolver()
    {
      // test whether we should solve on the coarse level. That is the case if we
      // have that level and if there was a redistribution on this level then our
      // communicator has to be valid (size()>0) as the smoother might try to communicate
      // in the constructor.
      if(buildHierarchy_ && matrices_->levels()==matrices_->maxlevels()
         && ( ! matrices_->redistributeInformation().back().isSetup() ||
              matrices_->parallelInformation().coarsest().getRedistributed().communicator().size() ) )
      {
        // We have the carsest level. Create the coarse Solver
        SmootherArgs sargs(smootherArgs_);
        sargs.iterations = 1;

        typename ConstructionTraits<Smoother>::Arguments cargs;
        cargs.setArgs(sargs);
        if(matrices_->redistributeInformation().back().isSetup()) {
          // Solve on the redistributed partitioning
          cargs.setMatrix(matrices_->matrices().coarsest().getRedistributed().getmat());
          cargs.setComm(matrices_->parallelInformation().coarsest().getRedistributed());
        }else{
          cargs.setMatrix(matrices_->matrices().coarsest()->getmat());
          cargs.setComm(*matrices_->parallelInformation().coarsest());
        }

        coarseSmoother_.reset(ConstructionTraits<Smoother>::construct(cargs));

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
        scalarProduct_ = createScalarProduct<X>(cargs.getComm(),category());
#else
        typedef Dune::ScalarProductChooser<X,ParallelInformation,category>
          ScalarProductChooser;
        // the scalar product.
        scalarProduct_.reset(ScalarProductChooser::construct(cargs.getComm()));
#endif


        typedef DirectSolverSelector< typename M::matrix_type, X > SolverSelector;

        // Use superlu if we are purely sequential or with only one processor on the coarsest level.
        if( SolverSelector::isDirectSolver &&
            (std::is_same<ParallelInformation,SequentialInformation>::value // sequential mode
           || matrices_->parallelInformation().coarsest()->communicator().size()==1 //parallel mode and only one processor
           || (matrices_->parallelInformation().coarsest().isRedistributed()
               && matrices_->parallelInformation().coarsest().getRedistributed().communicator().size()==1
               && matrices_->parallelInformation().coarsest().getRedistributed().communicator().size()>0) )
          )
        { // redistribute and 1 proc
          if(matrices_->parallelInformation().coarsest().isRedistributed())
          {
            if(matrices_->matrices().coarsest().getRedistributed().getmat().N()>0)
            {
              // We are still participating on this level
              solver_.reset(SolverSelector::create(matrices_->matrices().coarsest().getRedistributed().getmat(), false, false));
            }
            else
              solver_.reset();
          }
          else
          {
            solver_.reset(SolverSelector::create(matrices_->matrices().coarsest()->getmat(), false, false));
          }
          if(verbosity_>0 && matrices_->parallelInformation().coarsest()->communicator().rank()==0)
            std::cout<< "Using a direct coarse solver (" << SolverSelector::name() << ")" << std::endl;
        }
        else
        {
          if(matrices_->parallelInformation().coarsest().isRedistributed())
          {
            if(matrices_->matrices().coarsest().getRedistributed().getmat().N()>0)
              // We are still participating on this level

              // we have to allocate these types using the rebound allocator
              // in order to ensure that we fulfill the alignement requirements
              solver_.reset(new BiCGSTABSolver<X>(const_cast<M&>(matrices_->matrices().coarsest().getRedistributed()),
                                                  // Cast needed for Dune <=2.5
                                                  reinterpret_cast<typename
                                                                   std::conditional<std::is_same<PI,SequentialInformation>::value,
                                                                                    Dune::SeqScalarProduct<X>,
                                                                                    Dune::OverlappingSchwarzScalarProduct<X,PI> >::type&>(*scalarProduct_),
                                                  *coarseSmoother_, 1E-2, 1000, 0));
            else
              solver_.reset();
          }else
          {
              solver_.reset(new BiCGSTABSolver<X>(const_cast<M&>(*matrices_->matrices().coarsest()),
                                                  // Cast needed for Dune <=2.5
                                                  reinterpret_cast<typename
                                                                   std::conditional<std::is_same<PI,SequentialInformation>::value,
                                                                                    Dune::SeqScalarProduct<X>,
                                                                                    Dune::OverlappingSchwarzScalarProduct<X,PI> >::type&>(*scalarProduct_),
                                                  *coarseSmoother_, 1E-2, 1000, 0));
            // // we have to allocate these types using the rebound allocator
            // // in order to ensure that we fulfill the alignement requirements
            // using Alloc = typename A::template rebind<BiCGSTABSolver<X>>::other;
            // Alloc alloc;
            // auto p = alloc.allocate(1);
            // alloc.construct(p,
            //   const_cast<M&>(*matrices_->matrices().coarsest()),
            //   *scalarProduct_,
            //   *coarseSmoother_, 1E-2, 1000, 0);
            // solver_.reset(p,[](BiCGSTABSolver<X>* p){
            //     Alloc alloc;
            //     alloc.destroy(p);
            //     alloc.deallocate(p,1);
            //   });
          }
        }
      }
    }

    template<class M, class X, class S, class PI, class A>
    void AMGCPR<M,X,S,PI,A>::pre(Domain& x, Range& b)
    {
      // Detect Matrix rows where all offdiagonal entries are
      // zero and set x such that  A_dd*x_d=b_d
      // Thus users can be more careless when setting up their linear
      // systems.
      typedef typename M::matrix_type Matrix;
      typedef typename Matrix::ConstRowIterator RowIter;
      typedef typename Matrix::ConstColIterator ColIter;
      typedef typename Matrix::block_type Block;
      Block zero;
      zero=typename Matrix::field_type();

      const Matrix& mat=matrices_->matrices().finest()->getmat();
      for(RowIter row=mat.begin(); row!=mat.end(); ++row) {
        bool isDirichlet = true;
        bool hasDiagonal = false;
        Block diagonal;
        for(ColIter col=row->begin(); col!=row->end(); ++col) {
          if(row.index()==col.index()) {
            diagonal = *col;
            hasDiagonal = false;
          }else{
            if(*col!=zero)
              isDirichlet = false;
          }
        }
        if(isDirichlet && hasDiagonal)
          diagonal.solve(x[row.index()], b[row.index()]);
      }

      if(smoothers_->levels()>0)
        smoothers_->finest()->pre(x,b);
      else
        // No smoother to make x consistent! Do it by hand
        matrices_->parallelInformation().coarsest()->copyOwnerToAll(x,x);
      Range* copy = new Range(b);
      if(rhs_)
        delete rhs_;
      rhs_ = new Hierarchy<Range,A>(copy);
      Domain* dcopy = new Domain(x);
      if(lhs_)
        delete lhs_;
      lhs_ = new Hierarchy<Domain,A>(dcopy);
      dcopy = new Domain(x);
      if(update_)
        delete update_;
      update_ = new Hierarchy<Domain,A>(dcopy);
      matrices_->coarsenVector(*rhs_);
      matrices_->coarsenVector(*lhs_);
      matrices_->coarsenVector(*update_);

      // Preprocess all smoothers
      typedef typename Hierarchy<Smoother,A>::Iterator Iterator;
      typedef typename Hierarchy<Range,A>::Iterator RIterator;
      typedef typename Hierarchy<Domain,A>::Iterator DIterator;
      Iterator coarsest = smoothers_->coarsest();
      Iterator smoother = smoothers_->finest();
      RIterator rhs = rhs_->finest();
      DIterator lhs = lhs_->finest();
      if(smoothers_->levels()>0) {

        assert(lhs_->levels()==rhs_->levels());
        assert(smoothers_->levels()==lhs_->levels() || matrices_->levels()==matrices_->maxlevels());
        assert(smoothers_->levels()+1==lhs_->levels() || matrices_->levels()<matrices_->maxlevels());

        if(smoother!=coarsest)
          for(++smoother, ++lhs, ++rhs; smoother != coarsest; ++smoother, ++lhs, ++rhs)
            smoother->pre(*lhs,*rhs);
        smoother->pre(*lhs,*rhs);
      }


      // The preconditioner might change x and b. So we have to
      // copy the changes to the original vectors.
      x = *lhs_->finest();
      b = *rhs_->finest();

    }
    template<class M, class X, class S, class PI, class A>
    std::size_t AMGCPR<M,X,S,PI,A>::levels()
    {
      return matrices_->levels();
    }
    template<class M, class X, class S, class PI, class A>
    std::size_t AMGCPR<M,X,S,PI,A>::maxlevels()
    {
      return matrices_->maxlevels();
    }

    /** \copydoc Preconditioner::apply */
    template<class M, class X, class S, class PI, class A>
    void AMGCPR<M,X,S,PI,A>::apply(Domain& v, const Range& d)
    {
      LevelContext levelContext;

      if(additive) {
        *(rhs_->finest())=d;
        additiveMgc();
        v=*lhs_->finest();
      }else{
        // Init all iterators for the current level
        initIteratorsWithFineLevel(levelContext);


        *levelContext.lhs = v;
        *levelContext.rhs = d;
        *levelContext.update=0;
        levelContext.level=0;

        mgc(levelContext);

        if(postSteps_==0||matrices_->maxlevels()==1)
          levelContext.pinfo->copyOwnerToAll(*levelContext.update, *levelContext.update);

        v=*levelContext.update;
      }

    }

    template<class M, class X, class S, class PI, class A>
    void AMGCPR<M,X,S,PI,A>::initIteratorsWithFineLevel(LevelContext& levelContext)
    {
      levelContext.smoother = smoothers_->finest();
      levelContext.matrix = matrices_->matrices().finest();
      levelContext.pinfo = matrices_->parallelInformation().finest();
      levelContext.redist =
        matrices_->redistributeInformation().begin();
      levelContext.aggregates = matrices_->aggregatesMaps().begin();
      levelContext.lhs = lhs_->finest();
      levelContext.update = update_->finest();
      levelContext.rhs = rhs_->finest();
    }

    template<class M, class X, class S, class PI, class A>
    bool AMGCPR<M,X,S,PI,A>
    ::moveToCoarseLevel(LevelContext& levelContext)
    {

      bool processNextLevel=true;

      if(levelContext.redist->isSetup()) {
        levelContext.redist->redistribute(static_cast<const Range&>(*levelContext.rhs),
                             levelContext.rhs.getRedistributed());
        processNextLevel = levelContext.rhs.getRedistributed().size()>0;
        if(processNextLevel) {
          //restrict defect to coarse level right hand side.
          typename Hierarchy<Range,A>::Iterator fineRhs = levelContext.rhs++;
          ++levelContext.pinfo;
          Transfer<typename OperatorHierarchy::AggregatesMap::AggregateDescriptor,Range,ParallelInformation>
          ::restrictVector(*(*levelContext.aggregates), *levelContext.rhs,
                           static_cast<const Range&>(fineRhs.getRedistributed()),
                           *levelContext.pinfo);
        }
      }else{
        //restrict defect to coarse level right hand side.
        typename Hierarchy<Range,A>::Iterator fineRhs = levelContext.rhs++;
        ++levelContext.pinfo;
        Transfer<typename OperatorHierarchy::AggregatesMap::AggregateDescriptor,Range,ParallelInformation>
        ::restrictVector(*(*levelContext.aggregates),
                         *levelContext.rhs, static_cast<const Range&>(*fineRhs),
                         *levelContext.pinfo);
      }

      if(processNextLevel) {
        // prepare coarse system
        ++levelContext.lhs;
        ++levelContext.update;
        ++levelContext.matrix;
        ++levelContext.level;
        ++levelContext.redist;

        if(levelContext.matrix != matrices_->matrices().coarsest() || matrices_->levels()<matrices_->maxlevels()) {
          // next level is not the globally coarsest one
          ++levelContext.smoother;
          ++levelContext.aggregates;
        }
        // prepare the update on the next level
        *levelContext.update=0;
      }
      return processNextLevel;
    }

    template<class M, class X, class S, class PI, class A>
    void AMGCPR<M,X,S,PI,A>
    ::moveToFineLevel(LevelContext& levelContext, bool processNextLevel)
    {
      if(processNextLevel) {
        if(levelContext.matrix != matrices_->matrices().coarsest() || matrices_->levels()<matrices_->maxlevels()) {
          // previous level is not the globally coarsest one
          --levelContext.smoother;
          --levelContext.aggregates;
        }
        --levelContext.redist;
        --levelContext.level;
        //prolongate and add the correction (update is in coarse left hand side)
        --levelContext.matrix;

        //typename Hierarchy<Domain,A>::Iterator coarseLhs = lhs--;
        --levelContext.lhs;
        --levelContext.pinfo;
      }
      if(levelContext.redist->isSetup()) {
        // Need to redistribute during prolongateVector
        levelContext.lhs.getRedistributed()=0;
        Transfer<typename OperatorHierarchy::AggregatesMap::AggregateDescriptor,Range,ParallelInformation>
        ::prolongateVector(*(*levelContext.aggregates), *levelContext.update, *levelContext.lhs,
                           levelContext.lhs.getRedistributed(),
                           matrices_->getProlongationDampingFactor(),
                           *levelContext.pinfo, *levelContext.redist);
      }else{
        *levelContext.lhs=0;
        Transfer<typename OperatorHierarchy::AggregatesMap::AggregateDescriptor,Range,ParallelInformation>
        ::prolongateVector(*(*levelContext.aggregates), *levelContext.update, *levelContext.lhs,
                           matrices_->getProlongationDampingFactor(),
                           *levelContext.pinfo);
      }


      if(processNextLevel) {
        --levelContext.update;
        --levelContext.rhs;
      }

      *levelContext.update += *levelContext.lhs;
    }

    template<class M, class X, class S, class PI, class A>
    bool AMGCPR<M,X,S,PI,A>::usesDirectCoarseLevelSolver() const
    {
      return IsDirectSolver< CoarseSolver>::value;
    }

    template<class M, class X, class S, class PI, class A>
    void AMGCPR<M,X,S,PI,A>::mgc(LevelContext& levelContext){
      if(levelContext.matrix == matrices_->matrices().coarsest() && levels()==maxlevels()) {
        // Solve directly
        InverseOperatorResult res;
        res.converged=true; // If we do not compute this flag will not get updated
        if(levelContext.redist->isSetup()) {
          levelContext.redist->redistribute(*levelContext.rhs, levelContext.rhs.getRedistributed());
          if(levelContext.rhs.getRedistributed().size()>0) {
            // We are still participating in the computation
            levelContext.pinfo.getRedistributed().copyOwnerToAll(levelContext.rhs.getRedistributed(),
                                                    levelContext.rhs.getRedistributed());
            solver_->apply(levelContext.update.getRedistributed(),
                           levelContext.rhs.getRedistributed(), res);
          }
          levelContext.redist->redistributeBackward(*levelContext.update, levelContext.update.getRedistributed());
          levelContext.pinfo->copyOwnerToAll(*levelContext.update, *levelContext.update);
        }else{
          levelContext.pinfo->copyOwnerToAll(*levelContext.rhs, *levelContext.rhs);
          solver_->apply(*levelContext.update, *levelContext.rhs, res);
        }

        if (!res.converged)
          coarsesolverconverged = false;
      }else{
        // presmoothing
        presmooth(levelContext, preSteps_);

#ifndef DUNE_AMG_NO_COARSEGRIDCORRECTION
        bool processNextLevel = moveToCoarseLevel(levelContext);

        if(processNextLevel) {
          // next level
          for(std::size_t i=0; i<gamma_; i++)
            mgc(levelContext);
        }

        moveToFineLevel(levelContext, processNextLevel);
#else
        *lhs=0;
#endif

        if(levelContext.matrix == matrices_->matrices().finest()) {
          coarsesolverconverged = matrices_->parallelInformation().finest()->communicator().prod(coarsesolverconverged);
          if(!coarsesolverconverged){
            //DUNE_THROW(MathError, "Coarse solver did not converge");
          }
        }
        // postsmoothing
        postsmooth(levelContext, postSteps_);

      }
    }

    template<class M, class X, class S, class PI, class A>
    void AMGCPR<M,X,S,PI,A>::additiveMgc(){

      // restrict residual to all levels
      typename ParallelInformationHierarchy::Iterator pinfo=matrices_->parallelInformation().finest();
      typename Hierarchy<Range,A>::Iterator rhs=rhs_->finest();
      typename Hierarchy<Domain,A>::Iterator lhs = lhs_->finest();
      typename OperatorHierarchy::AggregatesMapList::const_iterator aggregates=matrices_->aggregatesMaps().begin();

      for(typename Hierarchy<Range,A>::Iterator fineRhs=rhs++; fineRhs != rhs_->coarsest(); fineRhs=rhs++, ++aggregates) {
        ++pinfo;
        Transfer<typename OperatorHierarchy::AggregatesMap::AggregateDescriptor,Range,ParallelInformation>
        ::restrictVector(*(*aggregates), *rhs, static_cast<const Range&>(*fineRhs), *pinfo);
      }

      // pinfo is invalid, set to coarsest level
      //pinfo = matrices_->parallelInformation().coarsest
      // calculate correction for all levels
      lhs = lhs_->finest();
      typename Hierarchy<Smoother,A>::Iterator smoother = smoothers_->finest();

      for(rhs=rhs_->finest(); rhs != rhs_->coarsest(); ++lhs, ++rhs, ++smoother) {
        // presmoothing
        *lhs=0;
        smoother->apply(*lhs, *rhs);
      }

      // Coarse level solve
#ifndef DUNE_AMG_NO_COARSEGRIDCORRECTION
      InverseOperatorResult res;
      pinfo->copyOwnerToAll(*rhs, *rhs);
      solver_->apply(*lhs, *rhs, res);

      if(!res.converged)
        DUNE_THROW(MathError, "Coarse solver did not converge");
#else
      *lhs=0;
#endif
      // Prologate and add up corrections from all levels
      --pinfo;
      --aggregates;

      for(typename Hierarchy<Domain,A>::Iterator coarseLhs = lhs--; coarseLhs != lhs_->finest(); coarseLhs = lhs--, --aggregates, --pinfo) {
        Transfer<typename OperatorHierarchy::AggregatesMap::AggregateDescriptor,Range,ParallelInformation>
        ::prolongateVector(*(*aggregates), *coarseLhs, *lhs, 1.0, *pinfo);
      }
    }


    /** \copydoc Preconditioner::post */
    template<class M, class X, class S, class PI, class A>
    void AMGCPR<M,X,S,PI,A>::post(Domain& x)
    {
      DUNE_UNUSED_PARAMETER(x);
      // Postprocess all smoothers
      typedef typename Hierarchy<Smoother,A>::Iterator Iterator;
      typedef typename Hierarchy<Domain,A>::Iterator DIterator;
      Iterator coarsest = smoothers_->coarsest();
      Iterator smoother = smoothers_->finest();
      DIterator lhs = lhs_->finest();
      if(smoothers_->levels()>0) {
        if(smoother != coarsest  || matrices_->levels()<matrices_->maxlevels())
          smoother->post(*lhs);
        if(smoother!=coarsest)
          for(++smoother, ++lhs; smoother != coarsest; ++smoother, ++lhs)
            smoother->post(*lhs);
        smoother->post(*lhs);
      }
      //delete &(*lhs_->finest());
      delete lhs_;
      lhs_=nullptr;
      //delete &(*update_->finest());
      delete update_;
      update_=nullptr;
      //delete &(*rhs_->finest());
      delete rhs_;
      rhs_=nullptr;
    }

    template<class M, class X, class S, class PI, class A>
    template<class A1>
    void AMGCPR<M,X,S,PI,A>::getCoarsestAggregateNumbers(std::vector<std::size_t,A1>& cont)
    {
      matrices_->getCoarsestAggregatesOnFinest(cont);
    }

  } // end namespace Amg
} // end namespace Dune

#endif
