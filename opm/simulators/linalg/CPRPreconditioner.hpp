/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 IRIS AS.
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU
  Copyright 2015 Statoil AS

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

#ifndef OPM_CPRPRECONDITIONER_HEADER_INCLUDED
#define OPM_CPRPRECONDITIONER_HEADER_INCLUDED

#include <memory>
#include <type_traits>

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/io.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/kamg.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>

#include <opm/simulators/linalg/ParallelRestrictedAdditiveSchwarz.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
namespace Opm
{

template<typename O, typename S, typename C,
         typename P, std::size_t COMPONENT_INDEX, std::size_t VARIABLE_INDEX>
class BlackoilAmg;

namespace Amg
{
    template<int Row, int Column>
    class Element;
}

namespace ISTLUtility
{

template<class T>
void setILUParameters(Opm::ParallelOverlappingILU0Args<T>& args,
                      const CPRParameter& params)
{
    args.setN(params.cpr_ilu_n_);
    args.setMilu(params.cpr_ilu_milu_);
}

template<class T>
void setILUParameters(Opm::ParallelOverlappingILU0Args<T>& args,
                      MILU_VARIANT milu, int n=0)
{
    args.setN(n);
    args.setMilu(milu);
}

template<class S, class P>
void setILUParameters(S&, const P&)
{}

template<class S, class P>
void setILUParameters(S&, bool, int)
{}

///
/// \brief A traits class for selecting the types of the preconditioner.
///
/// \tparam M The type of the matrix.
/// \tparam X The type of the domain of the linear problem.
/// \tparam Y The type of the range of the linear problem.
/// \tparam P The type of the parallel information.
////
template<class M, class X, class Y, class P>
struct CPRSelector
{
    /// \brief The information about the parallelization and communication
    using ParallelInformation = P;
    /// \brief The operator type;
    typedef Dune::OverlappingSchwarzOperator<M,X,X,ParallelInformation> Operator;
    typedef ParallelOverlappingILU0<M,X, X, ParallelInformation>
    EllipticPreconditioner;
    /// \brief The type of the unique pointer to the preconditioner of the elliptic part.
    typedef std::unique_ptr<EllipticPreconditioner>
    EllipticPreconditionerPointer;

    typedef EllipticPreconditioner Smoother;
    typedef Dune::Amg::AMG<Operator, X, Smoother, ParallelInformation> AMG;

    /// \brief creates an Operator from the matrix
    /// \param M The matrix to use.
    /// \param p The parallel information to use.
    static Operator* makeOperator(const M& m, const ParallelInformation& p)
    {
        return new Operator(m, p);
    }
};

template<class M, class X, class Y>
struct CPRSelector<M,X,Y,Dune::Amg::SequentialInformation>
{
    /// \brief The information about the parallelization and communication
    typedef Dune::Amg::SequentialInformation ParallelInformation;
    /// \brief The operator type;
    typedef Dune::MatrixAdapter<M, X, Y> Operator;
    /// \brief The type of the preconditioner used for the elliptic part.
    typedef ParallelOverlappingILU0<M,X, X, ParallelInformation>
    EllipticPreconditioner;
    /// \brief The type of the unique pointer to the preconditioner of the elliptic part.
    typedef std::unique_ptr<EllipticPreconditioner> EllipticPreconditionerPointer;

    /// \brief type of AMG used to precondition the elliptic system.
    typedef EllipticPreconditioner Smoother;
    typedef Dune::Amg::AMG<Operator, X, Smoother, ParallelInformation> AMG;

    /// \brief creates an Operator from the matrix
    /// \param M The matrix to use.
    /// \param p The parallel information to use.
    static Operator* makeOperator(const M& m, const ParallelInformation&)
    {
        return new Operator(m);
    }

};
//! \brief Creates and initializes a unique pointer to an sequential ILU0 preconditioner.
//! \param A     The matrix of the linear system to solve.
//! \param relax The relaxation factor to use.
template<class M, class X, class C>
std::shared_ptr<ParallelOverlappingILU0<M,X,X,C> >
createILU0Ptr(const M& A, const C& comm, double relax, MILU_VARIANT milu)
{
    return std::make_shared<ParallelOverlappingILU0<M,X,X,C> >(A, comm, relax, milu);
}
//! \brief Creates and initializes a shared pointer to an ILUn preconditioner.
//! \param A     The matrix of the linear system to solve.
//! \param ilu_n The n parameter for the extension of the nonzero pattern.
//! \param relax The relaxation factor to use.
template<class M, class X, class C>
std::shared_ptr<ParallelOverlappingILU0<M,X,X,C> >
createILUnPtr(const M& A, const C& comm, int ilu_n, double relax, MILU_VARIANT milu)
{
    return std::make_shared<ParallelOverlappingILU0<M,X,X,C> >( A, comm, ilu_n, relax, milu );
}
/// \brief Creates the elliptic preconditioner (ILU0)
/// \param Ae    The matrix of the elliptic system.
/// \param relax The relaxation parameter for ILU0.
/// \param milu  If true, the modified ilu approach is used. Dropped elements
///              will get added to the diagonal of U to preserve the row sum
///              for constant vectors (Ae = LUe).
/// \param comm  The object describing the parallelization information and communication.
template<class M, class X=typename M::range_type, class P>
typename CPRSelector<M,X,X,P>::EllipticPreconditionerPointer
createEllipticPreconditionerPointer(const M& Ae, double relax,
                                    MILU_VARIANT milu, const P& comm)
{
    typedef typename CPRSelector<M,X,X,P >
        ::EllipticPreconditioner ParallelPreconditioner;

    typedef typename CPRSelector<M,X,X,P>
        ::EllipticPreconditionerPointer EllipticPreconditionerPointer;
    return EllipticPreconditionerPointer(new ParallelPreconditioner(Ae, comm, relax, milu));
}

template < class C, class Op, class P, class S, std::size_t PressureEqnIndex, std::size_t PressureVarIndex, class Vector>
inline void
createAMGPreconditionerPointer(Op& opA, const double relax, const P& comm,
                               std::unique_ptr< BlackoilAmg<Op,S,C,P,PressureEqnIndex,PressureVarIndex> >& amgPtr,
                               const CPRParameter& params,
                               const Vector& weights)
{
    using AMG = BlackoilAmg<Op,S,C,P,PressureEqnIndex,PressureVarIndex>;
    int verbosity = 0;
    if (comm.communicator().rank() == 0) {
        verbosity = params.cpr_solver_verbose_;
    }
    // TODO: revise choice of parameters
    int coarsenTarget=1200;
    using Criterion = C;
    Criterion criterion(15, coarsenTarget);
    criterion.setDebugLevel( verbosity ); // no debug information, 1 for printing hierarchy information
    criterion.setDefaultValuesIsotropic(2);
    criterion.setNoPostSmoothSteps( 1 );
    criterion.setNoPreSmoothSteps( 1 );

    // Since DUNE 2.2 we also need to pass the smoother args instead of steps directly
    typedef typename AMG::Smoother Smoother;
    typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments  SmootherArgs;
    SmootherArgs  smootherArgs;
    smootherArgs.iterations = 1;
    smootherArgs.relaxationFactor = relax;
    setILUParameters(smootherArgs, params);

    amgPtr.reset( new AMG( params, weights, opA, criterion, smootherArgs, comm ) );
}

template < class C, class Op, class P, class AMG >
inline void
createAMGPreconditionerPointer(Op& opA, const double relax, const MILU_VARIANT milu, const P& comm, std::unique_ptr< AMG >& amgPtr)
{
    // TODO: revise choice of parameters
    int coarsenTarget=1200;
    using Criterion = C;
    Criterion criterion(15, coarsenTarget);
    criterion.setDebugLevel( 0 ); // no debug information, 1 for printing hierarchy information
    criterion.setDefaultValuesIsotropic(2);
    criterion.setNoPostSmoothSteps( 1 );
    criterion.setNoPreSmoothSteps( 1 );

    // for DUNE 2.2 we also need to pass the smoother args
    typedef typename AMG::Smoother Smoother;
    typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments  SmootherArgs;
    SmootherArgs  smootherArgs;
    smootherArgs.iterations = 1;
    smootherArgs.relaxationFactor = relax;
    setILUParameters(smootherArgs, milu);

    amgPtr.reset( new AMG(opA, criterion, smootherArgs, comm ) );
}

/// \brief Creates the elliptic preconditioner (ILU0)
/// \param opA     The operator representing the matrix of the system.
/// \param relax   The relaxation parameter for ILU0.
/// \param comm    The object describing the parallelization information and communication.
//  \param amgPtr  The unique_ptr to be filled (return)
template < int PressureEqnIndex, int PressureVarIndex, class Op, class P, class AMG >
inline void
createAMGPreconditionerPointer( Op& opA, const double relax, const MILU_VARIANT milu, const P& comm, std::unique_ptr< AMG >& amgPtr )
{
    // type of matrix
    typedef typename Op::matrix_type  M;

    // The coupling metric used in the AMG
    typedef Opm::Amg::Element<PressureEqnIndex, PressureVarIndex> CouplingMetric;

    // The coupling criterion used in the AMG
    typedef Dune::Amg::SymmetricCriterion<M, CouplingMetric> CritBase;

    // The coarsening criterion used in the AMG
    typedef Dune::Amg::CoarsenCriterion<CritBase> Criterion;

    createAMGPreconditionerPointer<Criterion>(opA, relax, milu, comm, amgPtr);
}

} // end namespace ISTLUtility

    /*!
      \brief CPR preconditioner.

      This is a two-stage preconditioner, combining an elliptic-type
      partial solution with ILU0 for the whole system.

      \tparam M The matrix type to operate on
      \tparam X Type of the update
      \tparam Y Type of the defect
      \tparam P Type of the parallel information. If not provided
                this will be Dune::Amg::SequentialInformation.
                The preconditioner is parallel if this is
                Dune::OwnerOverlapCopyCommunication<int,int>
    */
    template<class M, class X, class Y,
             class P=Dune::Amg::SequentialInformation>
    class CPRPreconditioner : public Dune::Preconditioner<X,Y>
    {
        // prohibit copying for now
        CPRPreconditioner( const CPRPreconditioner& );

    public:
        //! \brief The type describing the parallel information
        typedef P ParallelInformation;
        //! \brief The matrix type the preconditioner is for.
        typedef typename std::remove_const<M>::type matrix_type;
        //! \brief The domain type of the preconditioner.
        typedef X domain_type;
        //! \brief The range type of the preconditioner.
        typedef Y range_type;
        //! \brief The field type of the preconditioner.
        typedef typename X::field_type field_type;

        // define the category
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
        Dune::SolverCategory::Category category() const override
        {
          return std::is_same<P,Dune::Amg::SequentialInformation>::value ?
                 Dune::SolverCategory::sequential : Dune::SolverCategory::overlapping;
        }
#else
        enum {
            //! \brief The category the preconditioner is part of.
            category = std::is_same<P,Dune::Amg::SequentialInformation>::value?
            Dune::SolverCategory::sequential:Dune::SolverCategory::overlapping
        };
#endif

        typedef ISTLUtility::CPRSelector<M,X,X,P>  CPRSelectorType ;

        //! \brief Elliptic Operator
        typedef typename CPRSelectorType::Operator Operator;

        //! \brief preconditioner for the whole system (here either ILU(0) or ILU(n)
        typedef Dune::Preconditioner<X,X> WholeSystemPreconditioner;

        //! \brief type of the unique pointer to the ilu-0 preconditioner
        //! used the for the elliptic system
        typedef typename CPRSelectorType::EllipticPreconditionerPointer
        EllipticPreconditionerPointer;

        //! \brief amg preconditioner for the elliptic system
        typedef typename CPRSelectorType::AMG  AMG;

        /*! \brief Constructor.

          Constructor gets all parameters to operate the prec.
          \param A       The matrix to operate on.
          \param Ae      The top-left elliptic part of A.
          \param relax   The ILU0 relaxation factor.
          \param useAMG  if true, AMG is used as a preconditioner for the elliptic sub-system, otherwise ilu-0 (default)
          \param useBiCG if true, BiCG solver is used (default), otherwise CG solver
          \param paralleInformation The information about the parallelization, if this is a
                                    parallel run
        */
        CPRPreconditioner (const CPRParameter& param, const M& A, const M& Ae,
                           const ParallelInformation& comm=ParallelInformation(),
                           const ParallelInformation& commAe=ParallelInformation())
            : param_( param ),
              A_(A),
              Ae_(Ae),
              de_( Ae_.N() ),
              ve_( Ae_.M() ),
              dmodified_( A_.N() ),
              opAe_(CPRSelectorType::makeOperator(Ae_, commAe)),
              precond_(), // ilu0 preconditioner for elliptic system
              amg_(),     // amg  preconditioner for elliptic system
              pre_(), // copy A will be made be the preconditioner
              vilu_( A_.N() ),
              comm_(comm),
              commAe_(commAe)
        {
            // create appropriate preconditioner for elliptic system
            createEllipticPreconditioner( param_.cpr_use_amg_, commAe_ );

            // create the preconditioner for the whole system.
            if( param_.cpr_ilu_n_ == 0 ) {
                pre_ = ISTLUtility::createILU0Ptr<M,X>( A_, comm, param_.cpr_relax_, param_.cpr_ilu_milu_ );
            }
            else {
                pre_ = ISTLUtility::createILUnPtr<M,X>( A_, comm, param_.cpr_ilu_n_, param_.cpr_relax_, param_.cpr_ilu_milu_);
            }
        }

        /*!
          \brief Prepare the preconditioner.

          \copydoc Preconditioner::pre(X&,Y&)
        */
        virtual void pre (X& /*x*/, Y& /*b*/)
        {
        }

        /*!
          \brief Apply the preconditoner.

          \copydoc Preconditioner::apply(X&,const Y&)
        */
        virtual void apply (X& v, const Y& d)
        {
            // Extract part of d corresponding to elliptic part.
            // Note: Assumes that the elliptic part comes first.
            std::copy_n(d.begin(), de_.size(), de_.begin());

            // Solve elliptic part, extend solution to full.
            // reset result
            ve_ = 0;
            solveElliptic( ve_, de_ );

            //reset return value
            v = 0.0;
            // Again assuming that the elliptic part comes first.
            std::copy(ve_.begin(), ve_.end(), v.begin());

            // Subtract elliptic residual from initial residual.
            // dmodified = d - A * vfull
            dmodified_ = d;
            A_.mmv(v, dmodified_);
            // A is not parallel, do communication manually.
            comm_.copyOwnerToAll(dmodified_, dmodified_);

            // Apply Preconditioner for whole system (relax will be applied already)
            pre_->apply( vilu_, dmodified_);

            // don't apply relaxation if relax_ == 1
            if( std::abs( param_.cpr_relax_ - 1.0 ) < 1e-12 ) {
                v += vilu_;
            }
            else {
                v *= param_.cpr_relax_;
                v += vilu_;
            }
        }

        /*!
          \brief Clean up.

          \copydoc Preconditioner::post(X&)
        */
        virtual void post (X& /*x*/)
        {
        }

     protected:
        void solveElliptic(Y& x, Y& de)
        {
            // Linear solver parameters
            const double tolerance = param_.cpr_solver_tol_;
            const int maxit        = param_.cpr_max_ell_iter_;
            const int verbosity    = ( param_.cpr_solver_verbose_ &&
                                       comm_.communicator().rank()==0 ) ? 1 : 0;

            // operator result containing iterations etc.
            Dune::InverseOperatorResult result;

            // the scalar product chooser
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
            auto sp = Dune::createScalarProduct<X,ParallelInformation>(commAe_, category());
#else
            typedef Dune::ScalarProductChooser<X,ParallelInformation,category>
                ScalarProductChooser;
            // the scalar product.
            std::unique_ptr<typename ScalarProductChooser::ScalarProduct>
                sp(ScalarProductChooser::construct(commAe_));
#endif

            if( amg_ )
            {
                // Solve system with AMG
                if( param_.cpr_use_bicgstab_ ) {
                    Dune::BiCGSTABSolver<X> linsolve(*opAe_, *sp, (*amg_), tolerance, maxit, verbosity);
                    linsolve.apply(x, de, result);
                }
                else {
                    Dune::CGSolver<X> linsolve(*opAe_, *sp, (*amg_), tolerance, maxit, verbosity);
                    linsolve.apply(x, de, result);
                }
            }
            else
            {
                assert( precond_ );
                // Solve system with ILU-0
                if( param_.cpr_use_bicgstab_ ) {
                    Dune::BiCGSTABSolver<X> linsolve(*opAe_, *sp, (*precond_), tolerance, maxit, verbosity);
                    linsolve.apply(x, de, result);
                }
                else {
                    Dune::CGSolver<X> linsolve(*opAe_, *sp, (*precond_), tolerance, maxit, verbosity);
                    linsolve.apply(x, de, result);
                }

            }

            if (!result.converged) {
                OPM_THROW(LinearSolverProblem, "CPRPreconditioner failed to solve elliptic subsystem.");
            }
        }

        //! \brief Parameter collection for CPR
        const CPRParameter& param_;

        //! \brief The matrix for the full linear problem.
        const matrix_type& A_;
        //! \brief The elliptic part of the matrix.
        const matrix_type& Ae_;

        //! \brief temporary variables for elliptic solve
        Y de_, ve_, dmodified_;

        //! \brief elliptic operator
        std::unique_ptr<Operator> opAe_;

        //! \brief ILU0 preconditioner for the elliptic system
        EllipticPreconditionerPointer precond_;
        //! \brief AMG preconditioner with ILU0 smoother
        std::unique_ptr< AMG > amg_;

        //! \brief The preconditioner for the whole system
        //!
        //! We have to use a shared_ptr instead of a unique_ptr
        //! as we need to use a custom allocator based on dynamic
        //! information. But for unique_ptr the type of this deleter
        //! has to be available at coompile time.
        std::shared_ptr< WholeSystemPreconditioner > pre_;

        //! \brief temporary variables for ILU solve
        Y vilu_;

        //! \brief The information about the parallelization of the whole system.
        const P& comm_;
        //! \brief The information about the parallelization of the elliptic part
        //! of the system
        const P& commAe_;
     protected:
        void createEllipticPreconditioner( const bool amg, const P& comm )
        {
            if( amg )
            {
                ISTLUtility::createAMGPreconditionerPointer( *opAe_ , param_.cpr_relax_, param_.cpr_ilu_milu_, comm, amg_ );
            }
            else
            {
                precond_ = ISTLUtility::createEllipticPreconditionerPointer<M,X>( Ae_, param_.cpr_relax_, param_.cpr_ilu_milu_, comm);
            }
       }
    };


} // namespace Opm

#endif // OPM_CPRPRECONDITIONER_HEADER_INCLUDED
