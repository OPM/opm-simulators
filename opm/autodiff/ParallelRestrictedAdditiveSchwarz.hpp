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
#ifndef OPM_PARALLELRESTRICTEDADDITIVESCHWARZ_HEADER_INCLUDED
#define OPM_PARALLELRESTRICTEDADDITIVESCHWARZ_HEADER_INCLUDED

#include <dune/istl/preconditioner.hh>
#include <dune/istl/paamg/smoother.hh>
namespace Opm
{
template<class X, class Y, class C, class T>
class ParallelRestrictedOverlappingSchwarz;
}
namespace Dune
{
namespace Amg
{
template<class X, class Y, class C, class T>
struct ConstructionTraits<Opm::ParallelRestrictedOverlappingSchwarz<X,Y,C,T> >
{
    typedef DefaultParallelConstructionArgs<T,C> Arguments;
    typedef ConstructionTraits<T> SeqConstructionTraits;
    static inline Opm::ParallelRestrictedOverlappingSchwarz<X,Y,C,T>* construct(Arguments& args)
    {
        return new Opm::ParallelRestrictedOverlappingSchwarz<X,Y,C,T>(*SeqConstructionTraits::construct(args),
                                                                      args.getComm());
    }

    static inline void deconstruct(Opm::ParallelRestrictedOverlappingSchwarz<X,Y,C,T>* bp)
    {
        SeqConstructionTraits::deconstruct(static_cast<T*>(&bp->preconditioner));
        delete bp;
    }

};
template<class X, class Y, class C, class T>
struct SmootherTraits<Opm::ParallelRestrictedOverlappingSchwarz<X,Y,C,T> >
{
    typedef DefaultSmootherArgs<typename T::matrix_type::field_type> Arguments;

};

} // end namespace Amg
} // end namespace Dune

namespace Opm{

/**
 * @brief Block parallel preconditioner.
 *
 * This is essentially a wrapper that takes a sequential
 * preconditioner. In each step the sequential preconditioner
 * is applied and then all owner data points are updated on
 * all other processes.
 */
template<class X, class Y, class C, class T=Dune::Preconditioner<X,Y> >
class ParallelRestrictedOverlappingSchwarz : public Dune::Preconditioner<X,Y> {
    friend class Dune::Amg::ConstructionTraits<ParallelRestrictedOverlappingSchwarz<X,Y,C,T> >;
public:
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;
    //! \brief The type of the communication object.
    typedef C communication_type;

    // define the category
    enum {
        //! \brief The category the precondtioner is part of.
        category=Dune::SolverCategory::overlapping
    };

    /*! \brief Constructor.

      constructor gets all parameters to operate the prec.
      \param p The sequential preconditioner.
      \param c The communication object for syncing overlap and copy
      data points. (E.~g. OwnerOverlapCommunication )
    */
    ParallelRestrictedOverlappingSchwarz (T& p, const communication_type& c)
        : preconditioner(p), communication(c)
    {   }

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (X& x, Y& b)
    {
        communication.copyOwnerToAll(x,x);     // make dirichlet values consistent
        preconditioner.pre(x,b);
    }

    /*!
      \brief Apply the preconditioner

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (X& v, const Y& d)
    {
        Y& md = const_cast<Y&>(d);
        communication.copyOwnerToAll(md,md);
        preconditioner.apply(v,d);
        communication.copyOwnerToAll(v,v);
        communication.project(md);
    }

    template<bool forward>
    void apply (X& v, const Y& d)
    {
        Y& md = const_cast<Y&>(d);
        communication.copyOwnerToAll(md,md);
        preconditioner.template apply<forward>(v,d);
        communication.copyOwnerToAll(v,v);
        communication.project(md);
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (X& x)
    {
        preconditioner.post(x);
    }

private:
    //! \brief a sequential preconditioner
    T& preconditioner;

    //! \brief the communication object
    const communication_type& communication;
};


} // end namespace OPM
#endif
