/*
  Copyright 2015 Dr. Blatt - HPC-Simulation-Software & Services
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
#ifndef OPM_PARALLELOVERLAPPINGILU0_HEADER_INCLUDED
#define OPM_PARALLELOVERLAPPINGILU0_HEADER_INCLUDED

#include <dune/istl/preconditioner.hh>
#include <dune/istl/paamg/smoother.hh>

namespace Opm
{

template<class M, class X, class Y, class C>
class ParallelOverlappingILU0;

} // end namespace Opm

namespace Dune
{

namespace Amg
{

/// \brief Tells AMG how to construct the Opm::ParallelOverlappingILU0 smoother
/// \tparam Matrix The type of the Matrix.
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
template<class Matrix, class Domain, class Range, class ParallelInfo>
struct ConstructionTraits<Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo> >
{
    typedef Dune::SeqILU0<Matrix,Domain,Range> T;
    typedef DefaultParallelConstructionArgs<T,ParallelInfo> Arguments;
    typedef ConstructionTraits<T> SeqConstructionTraits;
    static inline Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo>* construct(Arguments& args)
    {
        return new Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo>(args.getMatrix(),
                                                         args.getComm(),
                                                         args.getArgs().relaxationFactor);
    }

    static inline void deconstruct(Opm::ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo>* bp)
    {
        delete bp;
    }

};

} // end namespace Amg

} // end namespace Dune

namespace Opm
{

/// \brief A two-step version of an overlapping Schwarz preconditioner using one step ILU0 as
///
/// This preconditioner differs from a ParallelRestrictedOverlappingSchwarz with
/// Dune:SeqILU0 in the follwing way:
/// During apply we make sure that the current residual is consistent (i.e.
/// each process knows the same value for each index. The we solve
/// Ly= d for y and make y consistent again. Last we solve Ux = y and
/// make sure that x is consistent.
/// In contrast for ParallelRestrictedOverlappingSchwarz we solve (LU)x = d for x
/// without forcing consistency between the two steps.
/// \tparam Matrix The type of the Matrix.
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
template<class Matrix, class Domain, class Range, class ParallelInfo>
class ParallelOverlappingILU0
    : public Dune::Preconditioner<Domain,Range> {
public:
    //! \brief The matrix type the preconditioner is for.
    typedef typename Dune::remove_const<Matrix>::type matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef Domain domain_type;
    //! \brief The range type of the preconditioner.
    typedef Range range_type;
    //! \brief The field type of the preconditioner.
    typedef typename Domain::field_type field_type;

    // define the category
    enum {
        //! \brief The category the preconditioner is part of.
        category=Dune::SolverCategory::overlapping
    };

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \param w The relaxation factor.
    */
    ParallelOverlappingILU0 (const Matrix& A, const ParallelInfo& comm,
                             field_type w)
        : ilu_(A), comm_(comm), w_(w)
    {
        int ilu_setup_successful = 1;
        std::string message;
        try
        {
            bilu0_decomposition(ilu_);
        }
        catch ( Dune::MatrixBlockError error )
        {
            message = error.what();
            std::cerr<<"Exception occured on process " <<
                comm.communicator().rank() << " during " <<
                "setup of ILU0 preconditioner with message: " <<
                message<<std::endl;
            ilu_setup_successful = 0;
        }
        // Check whether there was a problem on some process
        if ( comm.communicator().min(ilu_setup_successful) == 0 )
        {
            throw Dune::MatrixBlockError();
        }
    }

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (Domain& x, Range& b)
    {
        DUNE_UNUSED_PARAMETER(x);
        DUNE_UNUSED_PARAMETER(b);
    }

    /*!
      \brief Apply the preconditoner.

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (Domain& v, const Range& d)
    {
        Range& md = const_cast<Range&>(d);
        comm_.copyOwnerToAll(md,md);
        auto endrow=ilu_.end();
        for ( auto row = ilu_.begin(); row != endrow; ++row )
        {
            auto rhs(d[row.index()]);
            for ( auto col = row->begin(); col.index() < row.index(); ++col )
            {
                col->mmv(v[col.index()],rhs);
            }
            v[row.index()] = rhs;
        }
        comm_.copyOwnerToAll(v, v);
        auto rendrow = ilu_.beforeBegin();
        for( auto row = ilu_.beforeEnd(); row != rendrow; --row)
        {
            auto rhs(v[row.index()]);
            auto col = row->beforeEnd();
            for( ; col.index() > row.index(); --col)
                col->mmv(v[col.index()], rhs);
            v[row.index()] = 0;
            col->umv(rhs, v[row.index()]);
        }
        comm_.copyOwnerToAll(v, v);
        v *= w_;
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (Range& x)
    {
        DUNE_UNUSED_PARAMETER(x);
    }

private:
    //! \brief The ILU0 decomposition of the matrix.
    Matrix ilu_;
    const ParallelInfo& comm_;
    //! \brief The relaxation factor to use.
    field_type w_;

};

} // end namespace Opm
#endif
