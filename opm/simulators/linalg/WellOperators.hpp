/*
  Copyright 2016 IRIS AS
  Copyright 2019, 2020 Equinor ASA
  Copyright 2020 SINTEF

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

#ifndef OPM_WELLOPERATORS_HEADER_INCLUDED
#define OPM_WELLOPERATORS_HEADER_INCLUDED

#include <dune/common/parallel/communication.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <opm/common/TimingMacros.hpp>

#include <opm/simulators/linalg/matrixblock.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/istl/paamg/smoother.hh>

#include <cstddef>

namespace Opm {

//=====================================================================
// Implementation for ISTL-matrix based operators
// Note: the classes WellModelMatrixAdapter and
// WellModelGhostLastMatrixAdapter were moved from ISTLSolver.hpp
// and subsequently modified.
//=====================================================================

/// Linear operator wrapper for well model.
///
/// This class is intended to hide the actual type of the well model
/// (which depends on a TypeTag) by making a simple linear operator
/// wrapper. That way the WellModelMatrixAdapter template does not need
/// the concrete WellModel type, and we can avoid instantiating
/// WellModelMatrixAdapter separately for each TypeTag, as it will only
/// depend on the matrix and vector types involved, which typically are
/// just one for each block size with block sizes 1-4.
template <class X, class Y>
class LinearOperatorExtra : public Dune::LinearOperator<X, Y>
{
public:
    using field_type = typename X::field_type;
    using PressureMatrix = Dune::BCRSMatrix<MatrixBlock<field_type, 1, 1>>;
    virtual void addWellPressureEquations(PressureMatrix& jacobian,
                                          const X& weights,
                                          const bool use_well_weights) const = 0;
    virtual void addWellPressureEquationsStruct(PressureMatrix& jacobian) const = 0;
    virtual int getNumberOfExtraEquations() const = 0;
};

template <class WellModel, class X, class Y>
class WellModelAsLinearOperator : public LinearOperatorExtra<X, Y>
{
public:
    using Base = LinearOperatorExtra<X, Y>;
    using field_type = typename Base::field_type;
    using PressureMatrix = typename Base::PressureMatrix;
    explicit WellModelAsLinearOperator(const WellModel& wm)
        : wellMod_(wm)
    {
    }

    /*! \brief apply operator to x:  \f$ y = A(x) \f$
       The input vector is consistent and the output must also be
       consistent on the interior+border partition.
     */
    void apply(const X& x, Y& y) const override
    {
        OPM_TIMEBLOCK(apply);
        wellMod_.apply(x, y);
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    void applyscaleadd(field_type alpha, const X& x, Y& y) const override
    {
        OPM_TIMEBLOCK(applyscaleadd);
        wellMod_.applyScaleAdd(alpha, x, y);
    }

    /// Category for operator.
    /// This is somewhat tricky, I consider this operator sequential
    /// since (unlike WellModelMatrixAdapter) it does not do any
    /// projections etc. but only forwards the calls to the sequential
    /// well model.
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    void addWellPressureEquations(PressureMatrix& jacobian,
                                  const X& weights,
                                  const bool use_well_weights) const override
    {
        OPM_TIMEBLOCK(addWellPressureEquations);
        wellMod_.addWellPressureEquations(jacobian, weights, use_well_weights);
    }

    void addWellPressureEquationsStruct(PressureMatrix& jacobian) const override
    {
        OPM_TIMEBLOCK(addWellPressureEquationsStruct);
        wellMod_.addWellPressureEquationsStruct(jacobian);
    }

    int getNumberOfExtraEquations() const override
    {
        return wellMod_.numLocalWellsEnd();
    }

protected:
    const WellModel& wellMod_;
};


template <class WellModel, class X, class Y>
class DomainWellModelAsLinearOperator : public WellModelAsLinearOperator<WellModel, X, Y>
{
public:
    using WBase = WellModelAsLinearOperator<WellModel, X, Y>;
    using WBase::WBase;  // inherit all constructors from the base class
    using field_type = typename WBase::field_type;
    using PressureMatrix = typename WBase::PressureMatrix;

    void setDomainIndex(int index) { domainIndex_ = index; }

    void apply(const X& x, Y& y) const override
    {
        OPM_TIMEBLOCK(apply);
        this->wellMod_.applyDomain(x, y, domainIndex_);
    }

    void applyscaleadd(field_type alpha, const X& x, Y& y) const override
    {
        OPM_TIMEBLOCK(applyscaleadd);
        this->wellMod_.applyScaleAddDomain(alpha, x, y, domainIndex_);
    }

    void addWellPressureEquations(PressureMatrix& jacobian,
                                  const X& weights,
                                  const bool use_well_weights) const override
    {
        OPM_TIMEBLOCK(addWellPressureEquations);
        this->wellMod_.addWellPressureEquationsDomain(jacobian, weights, use_well_weights, domainIndex_);
    }

private:
    int domainIndex_ = -1;
};

/*!
   \brief Adapter to combine a matrix and another linear operator into
   a combined linear operator.

   Adapts a matrix A plus another linear operator W (typically from
   wells) to the assembled linear operator interface by returning S
   from getmat() and making apply() and applyscaleadd() apply both A
   and W to the input vector. In addition this is a parallel-aware
   adapter, that does not require the W operator to be parallel, but
   makes it into one by making the proper projections.
 */
template<class M, class X, class Y, bool overlapping >
class WellModelMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
{
public:
    using matrix_type = M;
    using domain_type = X;
    using range_type = Y;
    using field_type = typename X::field_type;
    using PressureMatrix = Dune::BCRSMatrix<MatrixBlock<field_type, 1, 1>>;
#if HAVE_MPI
    using communication_type = Dune::OwnerOverlapCopyCommunication<int,int>;
#else
    using communication_type = Dune::Communication<int>;
#endif

    Dune::SolverCategory::Category category() const override
    {
        return overlapping ?
               Dune::SolverCategory::overlapping : Dune::SolverCategory::sequential;
    }

    //! constructor: just store a reference to a matrix
    WellModelMatrixAdapter (const M& A,
                            const LinearOperatorExtra<X, Y>& wellOper,
                            const std::shared_ptr<communication_type>& comm = {})
        : A_( A ), wellOper_( wellOper ), comm_(comm)
    {}

    void apply( const X& x, Y& y ) const override
    {
      OPM_TIMEBLOCK(apply);
      A_.mv(x, y);

      // add well model modification to y
      wellOper_.apply(x, y);

  #if HAVE_MPI
      if (comm_) {
        comm_->project(y);
      }
  #endif
    }

    // y += \alpha * A * x
    void applyscaleadd (field_type alpha, const X& x, Y& y) const override
    {
      OPM_TIMEBLOCK(applyscaleadd);
      A_.usmv(alpha, x, y);

      // add scaled well model modification to y
      wellOper_.applyscaleadd(alpha, x, y);

  #if HAVE_MPI
      if (comm_) {
          comm_->project( y );
      }
  #endif
    }

    const matrix_type& getmat() const override { return A_; }

    void addWellPressureEquations(PressureMatrix& jacobian,
                                  const X& weights,
                                  const bool use_well_weights) const
    {
        OPM_TIMEBLOCK(addWellPressureEquations);
        wellOper_.addWellPressureEquations(jacobian, weights, use_well_weights);
    }

    void addWellPressureEquationsStruct(PressureMatrix& jacobian) const
    {
        OPM_TIMEBLOCK(addWellPressureEquations);
        wellOper_.addWellPressureEquationsStruct(jacobian);
    }

    int getNumberOfExtraEquations() const
    {
        return wellOper_.getNumberOfExtraEquations();
    }

protected:
    const matrix_type& A_ ;
    const LinearOperatorExtra<X, Y>& wellOper_;
    std::shared_ptr<communication_type> comm_;
};

/*!
   \brief Adapter to combine a matrix and another linear operator into
   a combined linear operator.

   This is similar to WellModelMatrixAdapter, with the difference that
   here we assume a parallel ordering of rows, where ghost rows are
   located after interior rows.
 */
template<class M, class X, class Y, bool overlapping >
class WellModelGhostLastMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
{
public:
    using matrix_type = M;
    using domain_type = X;
    using range_type = Y;
    using field_type = typename X::field_type;
    using PressureMatrix = Dune::BCRSMatrix<MatrixBlock<field_type, 1, 1>>;
#if HAVE_MPI
    using communication_type = Dune::OwnerOverlapCopyCommunication<int,int>;
#else
    using communication_type = Dune::Communication<int>;
#endif

    Dune::SolverCategory::Category category() const override
    {
        return overlapping ?
            Dune::SolverCategory::overlapping : Dune::SolverCategory::sequential;
    }

    //! constructor: just store a reference to a matrix
    WellModelGhostLastMatrixAdapter (const M& A,
                                     const LinearOperatorExtra<X, Y>& wellOper,
                                     const std::size_t interiorSize )
        : A_( A ), wellOper_( wellOper ), interiorSize_(interiorSize)
    {}

    void apply(const X& x, Y& y) const override
    {
        OPM_TIMEBLOCK(apply);
        for (auto row = A_.begin(); row.index() < interiorSize_; ++row)
        {
            y[row.index()]=0;
            auto endc = (*row).end();
            for (auto col = (*row).begin(); col != endc; ++col)
                (*col).umv(x[col.index()], y[row.index()]);
        }

        // add well model modification to y
        wellOper_.apply(x, y);

        ghostLastProject(y);
    }

    // y += \alpha * A * x
    void applyscaleadd (field_type alpha, const X& x, Y& y) const override
    {
        OPM_TIMEBLOCK(applyscaleadd);
        for (auto row = A_.begin(); row.index() < interiorSize_; ++row)
        {
            auto endc = (*row).end();
            for (auto col = (*row).begin(); col != endc; ++col)
                (*col).usmv(alpha, x[col.index()], y[row.index()]);
        }
        // add scaled well model modification to y
        wellOper_.applyscaleadd(alpha, x, y);

        ghostLastProject(y);
    }

    const matrix_type& getmat() const override { return A_; }

    void addWellPressureEquations(PressureMatrix& jacobian,
                                  const X& weights,
                                  const bool use_well_weights) const
    {
        OPM_TIMEBLOCK(addWellPressureEquations);
        wellOper_.addWellPressureEquations(jacobian, weights, use_well_weights);
    }

    void addWellPressureEquationsStruct(PressureMatrix& jacobian) const
    {
        OPM_TIMEBLOCK(addWellPressureEquationsStruct);
        wellOper_.addWellPressureEquationsStruct(jacobian);
    }

    int getNumberOfExtraEquations() const
    {
        return wellOper_.getNumberOfExtraEquations();
    }

protected:
    void ghostLastProject(Y& y) const
    {
        std::size_t end = y.size();
        for (std::size_t i = interiorSize_; i < end; ++i)
            y[i] = 0;
    }

    const matrix_type& A_ ;
    const LinearOperatorExtra<X, Y>& wellOper_;
    std::size_t interiorSize_;
};

/*!
   \brief Dune linear operator that assumes ghost rows are ordered after 
   interior rows. Avoids some computations because of this.
 
   This is similar to WellModelGhostLastMatrixAdapter, with the difference that
   here we do not have a well model, and also do calcilate the interiorSize
   using the parallel index set. Created for use in AMG/CPR smoothers.
 */
template<class M, class X, class Y, class C>
class GhostLastMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
{
public:
    typedef M matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;


    typedef C communication_type;

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::overlapping;
    }

    //! constructor: just store a reference to a matrix
    GhostLastMatrixAdapter (const M& A,
                            const communication_type& comm)
        : A_( Dune::stackobject_to_shared_ptr(A) ), comm_(comm)
    {
        interiorSize_ = setInteriorSize(comm_);
    }

    GhostLastMatrixAdapter (const std::shared_ptr<M> A,
                            const communication_type& comm)
        : A_( A ), comm_(comm)
    {
        interiorSize_ = setInteriorSize(comm_);
    }

    virtual void apply( const X& x, Y& y ) const override
    {
        for (auto row = A_->begin(); row.index() < interiorSize_; ++row)
        {
            y[row.index()]=0;
            auto endc = (*row).end();
            for (auto col = (*row).begin(); col != endc; ++col)
                (*col).umv(x[col.index()], y[row.index()]);
        }

        ghostLastProject( y );
    }

    // y += \alpha * A * x
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const override
    {
        for (auto row = A_->begin(); row.index() < interiorSize_; ++row)
        {
            auto endc = (*row).end();
            for (auto col = (*row).begin(); col != endc; ++col)
                (*col).usmv(alpha, x[col.index()], y[row.index()]);
        }

        ghostLastProject( y );
    }

    virtual const matrix_type& getmat() const override { return *A_; }

    size_t getInteriorSize() const { return interiorSize_;}

private:
    void ghostLastProject(Y& y) const
    {
        size_t end = y.size();
        for (size_t i = interiorSize_; i < end; ++i)
            y[i] = 0; //project to interiorsize, i.e. ignore ghost
    }

    size_t setInteriorSize(const communication_type& comm) const
    {
        auto indexSet = comm.indexSet();

        size_t is = 0;
        // Loop over index set
        for (auto idx = indexSet.begin(); idx!=indexSet.end(); ++idx) {
            //Only take "owner" indices
            if (idx->local().attribute()==1) {
                //get local index
                auto loc = idx->local().local();
                // if loc is higher than "old interior size", update it
                if (loc > is) {
                    is = loc;
                }
            }
        }
        return is + 1; //size is plus 1 since we start at 0
    }
    const std::shared_ptr<const matrix_type> A_ ;
    const communication_type&  comm_;
    size_t interiorSize_;
};

} // namespace Opm

namespace Dune {
namespace Amg {

template<class M, class X, class Y, class C>
class ConstructionTraits<Opm::GhostLastMatrixAdapter<M,X,Y,C> >
{
public:
    typedef ParallelOperatorArgs<M,C> Arguments;

    static inline std::shared_ptr<Opm::GhostLastMatrixAdapter<M,X,Y,C>> construct(const Arguments& args)
    {
        return std::make_shared<Opm::GhostLastMatrixAdapter<M,X,Y,C>>
            (args.matrix_, args.comm_);
    }
};

} // end namespace Amg
} // end namespace Dune

#endif // OPM_WELLOPERATORS_HEADER_INCLUDED

