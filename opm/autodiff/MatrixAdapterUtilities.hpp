#ifndef MATRIXADAPTERUTILITIES_HPP
#define MATRIXADAPTERUTILITIES_HPP

#include "transpose.hh"
#include <dune/istl/operators.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/timer.hh>
#include <dune/common/unused.hh>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/any.hpp>

namespace Opm {


template<class M, class X, class Y, class WellModel, class Grid, bool overlapping >
class WellModelMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
{
    typedef Dune::AssembledLinearOperator<M,X,Y> BaseType;

public:
    typedef M matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

#if HAVE_MPI
    typedef Dune::OwnerOverlapCopyCommunication<int,int> communication_type;
#else
    typedef Dune::CollectiveCommunication< Grid > communication_type;
#endif

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
    Dune::SolverCategory::Category category() const override
    {
        return overlapping ?
                    Dune::SolverCategory::overlapping : Dune::SolverCategory::sequential;
    }
#else
    enum {
        //! \brief The solver category.
        category = overlapping ?
        Dune::SolverCategory::overlapping :
        Dune::SolverCategory::sequential
    };
#endif

    //! constructor: just store a reference to a matrix
    WellModelMatrixAdapter (const M& A, const WellModel& wellMod, const boost::any& parallelInformation = boost::any() )
        : A_( A ), wellMod_( wellMod ), comm_()
    {
#if HAVE_MPI
        if( parallelInformation.type() == typeid(ParallelISTLInformation) )
        {
            const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>( parallelInformation);
            comm_.reset( new communication_type( info.communicator() ) );
        }
#endif
    }

    virtual void apply( const X& x, Y& y ) const
    {
        A_.mv( x, y );
        // add well model modification to y
        wellMod_.apply(x, y );

#if HAVE_MPI
        if( comm_ )
            comm_->project( y );
#endif
    }

    // y += \alpha * A * x
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
        A_.usmv(alpha,x,y);
        // add scaled well model modification to y
        wellMod_.applyScaleAdd( alpha, x, y );

#if HAVE_MPI
        if( comm_ )
            comm_->project( y );
#endif
    }

    virtual const matrix_type& getmat() const { return A_; }

    communication_type* comm()
    {
        return comm_.operator->();
    }

protected:
    const matrix_type& A_ ;
    const WellModel& wellMod_;
    std::unique_ptr< communication_type > comm_;
};

template<class M, class X, class Y, class WellModel, class Grid, bool overlapping >
class WellModelTransposeMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
{
    typedef Dune::AssembledLinearOperator<M,X,Y> BaseType;

public:
    typedef M matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

#if HAVE_MPI
    typedef Dune::OwnerOverlapCopyCommunication<int,int> communication_type;
#else
    typedef Dune::CollectiveCommunication< Grid > communication_type;
#endif

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
    Dune::SolverCategory::Category category() const override
    {
        return overlapping ?
                    Dune::SolverCategory::overlapping : Dune::SolverCategory::sequential;
    }
#else
    enum {
        //! \brief The solver category.
        category = overlapping ?
        Dune::SolverCategory::overlapping :
        Dune::SolverCategory::sequential
    };
#endif

    //! constructor: just store a reference to a matrix
    WellModelTransposeMatrixAdapter (const M& A,
                                     const M& A_for_precond,
                                     const WellModel& wellMod,
                                     const boost::any& parallelInformation = boost::any() )
        : A_( A ),  A_for_precond_(A_for_precond), wellMod_( wellMod ), comm_()
    {
       Dune::MatrixVector::transpose<matrix_type>(A_for_precond_, AT_precond_);
  //      AT_precond_ = A_for_precond_.transpose();
#if HAVE_MPI
        if( parallelInformation.type() == typeid(ParallelISTLInformation) )
        {
            const ParallelISTLInformation& info =
                    boost::any_cast<const ParallelISTLInformation&>( parallelInformation);
            comm_.reset( new communication_type( info.communicator() ) );
        }
#endif
    }

    virtual void apply( const X& x, Y& y ) const
    {
        A_.mtv( x, y );
        // add well model modification to y
        wellMod_.applyt(x, y );

#if HAVE_MPI
        if( comm_ )
            comm_->project( y );
#endif
    }

    // y += \alpha * A * x
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
        A_.usmtv(alpha,x,y);
        // add scaled well model modification to y
        wellMod_.applyScaleAdd( alpha, x, y );

#if HAVE_MPI
        if( comm_ )
            comm_->project( y );
#endif
    }

    virtual const matrix_type& getmat() const
    {
        //Dune::MatrixVector::Transposed<matrix_type> AT;
        //Dune::MatrixVector::transpose<matrix_type>(A_, AT);
        return AT_precond_;
        //return Dune::MatrixVector::TransposeHelper<matrix_type>_;
    }

    communication_type* comm()
    {
        return comm_.operator->();
    }

protected:
    const matrix_type& A_ ;
    const matrix_type& A_for_precond_;
    Dune::MatrixVector::Transposed<matrix_type> AT_precond_ ;
    const WellModel& wellMod_;
    std::unique_ptr< communication_type > comm_;
};



}

#endif // MATRIXADAPTERUTILITIES_HPP
