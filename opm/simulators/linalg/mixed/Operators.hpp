#ifndef OPM_MIXED_OPERATORS_HEADER_INCLUDED
#define OPM_MIXED_OPERATORS_HEADER_INCLUDED


namespace Opm
{

//! @brief Adapter to take advantage of the fact that all matrix rows
//! associated with ghost cells are located at the end of the matrix
//!
//! @note The underlying mixed-matrix already ignores ghost rows.
//!
//! @param M matrix class
//! @param V vector class
//! @param C communicator class
template<class M, class V, class C>
class MixedGhostLastMatrixAdapter : public Dune::AssembledLinearOperator<M,V,V>
{
public:

    //! constructor: just store a reference to matrix and communicator
    MixedGhostLastMatrixAdapter (const M& A, const C& comm) : A_( A ), comm_(comm) {}

    // y = A * x
    virtual void apply( const V& x, V& y ) const override
    {
        A_.mv(x,y);
        comm_.project(y);
    }

    // y += \alpha * A * x
    virtual void applyscaleadd (double alpha, const V& x, V& y) const override
    {
        A_.usmv(alpha,x,y);
        comm_.project(y);
    }

    // accessor to matix object
    virtual const M& getmat() const override { return A_; }

    // solver category
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::overlapping;
    }

private:
    const M& A_;
    const C& comm_;
};


//! @brief Adapter to combine a matrix with another linear operator while
//! taking advantage of the fact that all matrix rows associated with ghost
//! cells are located at the end of the matrix
//!
//! @note The underlying mixed-matrix already ignores ghost rows.
//!
//! @param M matrix class
//! @param V vector class
//! @param C communicator class
template<class M, class V, class C>
class WellModelMixedGhostLastMatrixAdapter : public Dune::AssembledLinearOperator<M,V,V>
{
public:
    using field_type = typename V::field_type;
    using PressureMatrix = Dune::BCRSMatrix<MatrixBlock<field_type, 1, 1>>;

    //! constructor: just store a reference to a matrix
    WellModelMixedGhostLastMatrixAdapter (const M& A,
                                     const LinearOperatorExtra<V, V>& wellOper,
                                     const C& comm
                                    )
        : A_( A ), wellOper_( wellOper ), comm_ ( comm )
    {}

    // y = A * x
    virtual void apply( const V& x, V& y ) const override
    {
        A_.mv(x,y);
        wellOper_.apply(x, y);
        comm_.project(y);
    }

    // y += \alpha * A * x
    virtual void applyscaleadd (double alpha, const V& x, V& y) const override
    {
        A_.usmv(alpha,x,y);
        wellOper_.applyscaleadd(alpha, x, y);
        comm_.project(y);
    }

    // accessor to matix object
    const M& getmat() const override { return A_; }

    void addWellPressureEquations(PressureMatrix& jacobian,
                                  const V& weights,
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

    // solver category
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::overlapping;
    }


protected:

    const M& A_ ;
    const C& comm_ ;
    const LinearOperatorExtra<V, V>& wellOper_;
};

} // namespace Opm
#endif // OPM_MIXED_OPERATORS_HEADER_INCLUDED
