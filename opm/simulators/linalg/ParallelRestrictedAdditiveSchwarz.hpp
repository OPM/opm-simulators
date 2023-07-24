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

#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/paamg/smoother.hh>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>

namespace Opm
{

template<class X, class Y, class C, class T>
class ParallelRestrictedOverlappingSchwarz;

} // end namespace Opm

namespace Dune
{

namespace Amg
{

/// \brief Tells AMG how to construct the Opm::ParallelOverlappingILU0 smoother
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
/// \tparam SeqPreconditioner The underlying sequential preconditioner to use.
template<class Range, class Domain, class ParallelInfo, class SeqPreconditioner>
struct ConstructionTraits<Opm::ParallelRestrictedOverlappingSchwarz<Range,
                                                                    Domain,
                                                                    ParallelInfo,
                                                                    SeqPreconditioner> >
{
    typedef DefaultParallelConstructionArgs<SeqPreconditioner,ParallelInfo> Arguments;
    typedef ConstructionTraits<SeqPreconditioner> SeqConstructionTraits;

    /// \brief Construct a parallel restricted overlapping schwarz preconditioner.
    typedef std::shared_ptr< Opm::ParallelRestrictedOverlappingSchwarz<Range,
                                                                       Domain,
                                                                       ParallelInfo,
                                                                       SeqPreconditioner> > ParallelRestrictedOverlappingSchwarzPointer;

    static inline ParallelRestrictedOverlappingSchwarzPointer
    construct(Arguments& args)
    {
        using PROS =
            Opm::ParallelRestrictedOverlappingSchwarz<Range,Domain,
                                                      ParallelInfo,SeqPreconditioner>;
        return std::make_shared<PROS>(*SeqConstructionTraits::construct(args),
                                      args.getComm());
    }

    /// \brief Deconstruct and free a parallel restricted overlapping schwarz preconditioner.
    static inline void deconstruct(Opm::ParallelRestrictedOverlappingSchwarz
                                   <Range,Domain,ParallelInfo,SeqPreconditioner>* bp)
    {
        SeqConstructionTraits
            ::deconstruct(static_cast<SeqPreconditioner*>(&bp->preconditioner));
        delete bp;
    }

};

/// \brief Tells AMG how to use Opm::ParallelOverlappingILU0 smoother
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
/// \tparam SeqPreconditioner The underlying sequential preconditioner to use.
template<class Range, class Domain, class ParallelInfo, class SeqPreconditioner>
struct SmootherTraits<Opm::ParallelRestrictedOverlappingSchwarz<Range,
                                                                Domain,
                                                                ParallelInfo,
                                                                SeqPreconditioner> >
{
    typedef DefaultSmootherArgs<typename SeqPreconditioner::matrix_type::field_type> Arguments;

};

} // end namespace Amg

} // end namespace Dune

namespace Opm{

/// \brief Block parallel preconditioner.
///
/// This is essentially a wrapper that takes a sequential
/// preconditioner. In each step the sequential preconditioner
/// is applied to the whole subdomain and then all owner data
/// points are updated on all other processes from the processor
/// that knows the complete matrix row for this data point (in dune-istl
/// speak that is the one that owns the data).
///
/// Note that this is different from the usual approach in dune-istl where
/// the application of the sequential preconditioner only takes place on
/// the (owner) partition of the process disregarding any overlap/ghost region.
///
/// For more information see https://www.cs.colorado.edu/~cai/papers/rash.pdf
///
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
/// \tparam SeqPreconditioner The underlying sequential preconditioner to use.
template<class Range, class Domain, class ParallelInfo, class SeqPreconditioner=Dune::Preconditioner<Range,Domain> >
class ParallelRestrictedOverlappingSchwarz
    : public Dune::Preconditioner<Range,Domain> {
    friend class Dune::Amg
    ::ConstructionTraits<ParallelRestrictedOverlappingSchwarz<Range,
                                                              Domain,
                                                              ParallelInfo,
                                                              SeqPreconditioner> >;
public:
    //! \brief The domain type of the preconditioner.
    typedef Domain domain_type;
    //! \brief The range type of the preconditioner.
    typedef Range range_type;
    //! \brief The field type of the preconditioner.
    typedef typename Domain::field_type field_type;
    //! \brief The type of the communication object.
    typedef ParallelInfo communication_type;

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
    ParallelRestrictedOverlappingSchwarz (SeqPreconditioner& p, const communication_type& c)
        : preconditioner_(p), communication_(c)
    {   }

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (Domain& x, Range& b)
    {
        OPM_TIMEBLOCK(pre);
        communication_.copyOwnerToAll(x,x);     // make dirichlet values consistent
        preconditioner_.pre(x,b);
    }

    /*!
      \brief Apply the preconditioner

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (Domain& v, const Range& d)
    {
        apply<true>(v, d);
    }

    template<bool forward>
    void apply (Domain& v, const Range& d)
    {
        OPM_TIMEBLOCK(apply);
        // hack us a mutable d to prevent copying.
        Range& md = const_cast<Range&>(d);
        communication_.copyOwnerToAll(md,md);
        preconditioner_.template apply<forward>(v,d);
        communication_.copyOwnerToAll(v,v);
        // Make sure that d is the same as at the beginning of apply.
        communication_.project(md);
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (Range& x)
    {
        OPM_TIMEBLOCK(post);
        preconditioner_.post(x);
    }

private:
    //! \brief a sequential preconditioner
    SeqPreconditioner& preconditioner_;

    //! \brief the communication object
    const communication_type& communication_;
};


} // end namespace OPM
#endif
