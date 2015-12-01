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

template<class M, class X, class Y, class C>
struct ConstructionTraits<Opm::ParallelOverlappingILU0<M,X,Y,C> >
{
    typedef Dune::SeqILU0<M,X,Y> T;
    typedef DefaultParallelConstructionArgs<T,C> Arguments;
    typedef ConstructionTraits<T> SeqConstructionTraits;
    static inline Opm::ParallelOverlappingILU0<M,X,Y,C>* construct(Arguments& args)
    {
        return new Opm::ParallelOverlappingILU0<M,X,Y,C>(args.getMatrix(),
                                                         args.getComm(),
                                                         args.getArgs().relaxationFactor);
    }

    static inline void deconstruct(Opm::ParallelOverlappingILU0<M,X,Y,C>* bp)
    {
        delete bp;
    }

};

} // end namespace Amg

} // end namespace Dune

namespace Opm
{

template<class M, class X, class Y, typename C>
class ParallelOverlappingILU0 : public Dune::Preconditioner<X,Y> {
public:
    //! \brief The matrix type the preconditioner is for.
    typedef typename Dune::remove_const<M>::type matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;

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
    ParallelOverlappingILU0 (const M& A, const C& comm, field_type w)
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
    virtual void pre (X& x, Y& b)
    {
        DUNE_UNUSED_PARAMETER(x);
        DUNE_UNUSED_PARAMETER(b);
    }

    /*!
      \brief Apply the preconditoner.

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (X& v, const Y& d)
    {
        Y& md = const_cast<Y&>(d);
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
    virtual void post (X& x)
    {
        DUNE_UNUSED_PARAMETER(x);
    }

private:
    //! \brief The ILU0 decomposition of the matrix.
    matrix_type ilu_;
    const C& comm_;
    //! \brief The relaxation factor to use.
    field_type w_;

};

} // end namespace Opm
#endif
