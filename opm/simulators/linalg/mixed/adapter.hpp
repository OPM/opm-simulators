#ifndef OPM_MIXED_ADAPTER_HEADER_INCLUDED
#define OPM_MIXED_ADAPTER_HEADER_INCLUDED

#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/owneroverlapcopy.hh>

#include <dune/istl/matrixindexset.hh>

#include <opm/simulators/linalg/mixed/MatrixWrapper.hpp>



namespace Dune
{
//#include <dune/istl/solvers.hh>


//! @brief Optimized sequential scalar product.
//!
//! @tparam Vector block-vector class with data stored as contiguous double array
template<class Vector>
class SeqOptmizedProduct : public Dune::SeqScalarProduct<Vector>
{
public:

    // extract block size
    static constexpr auto block_size = Vector::block_type::dimension;

    // Compute the dot product <vx, vy>
    virtual double dot(const Vector& vx, const Vector& vy) const override
    {
        // access underlying data
        double const *x = &vx[0][0];
        double const *y = &vy[0][0];

        // total array length
        int NN = block_size*vx.N();

        // unroll loop in multiples of 8
        int n=NN/8;
        int N=8*n;
        double agg[8];
        for(int i=0;i<8;i++) agg[i]=0.0;
        for(int i=0;i<N;i+=8) for(int j=0;j<8;j++) agg[j]+=x[i+j]*y[i+j];
        for(int j=0;j<4;j++) agg[j]+=agg[j+4];
        for(int j=0;j<2;j++) agg[j]+=agg[j+2];
        for(int j=0;j<1;j++) agg[j]+=agg[j+1];

        // trailing end
        for(int j=N;j<NN;j++) agg[0]+=x[j]*y[j];

        return agg[0];
    }

    // Compute the norm ||x||
    virtual double norm(const Vector& x) const override {
        return std::sqrt(this->dot(x, x));
    }
};

//! @brief Generalized mixed precision operator interface
//!
//! @tparam Matrix the block-matrix used by linear operator
//! @tparam Vector the block-vector used by linear operator
//! @tparam Comm the communicator used by linear operator
template <class Matrix, class Vector, class Comm>
struct MixedOperator
{
    using type = Dune::OverlappingSchwarzOperator<Matrix, Vector, Vector, Comm>;
};

//! @brief Generalized mixed precision operator interface
//!
//! @tparam Matrix the block-matrix used by linear operator
//! @tparam Vector the block-vector used by linear operator
template <class Matrix, class Vector>
struct MixedOperator<Matrix, Vector, Dune::Amg::SequentialInformation>
{
    using type = Dune::MatrixAdapter<Matrix, Vector, Vector>;
};


//! @brief Wraps mixed precision
//!
//! @tparam Comm the communicator passed to FlexibleLinearSolver
//! @tparam Operator the linear operator passed to FlexibleLinearSolver
//! @tparam Vector the block vector passed to FlexibleLinearSolver
template <class Comm, class Operator, class Vector>
class MixedAdapter:public InverseOperator<Vector, Vector>
{
    public:

    using AbstractPrecondType = Dune::PreconditionerWithUpdate<Vector,Vector>;
    using AbstractScalarProductType = Dune::ScalarProduct<Vector>;

    static constexpr auto block_size = Vector::block_type::dimension;
    using MixedMatrixType      = MatrixWrapper<Vector, block_size>;
    using MixedOperatorType    = MixedOperator<MixedMatrixType, Vector, Comm>::type;
    using OptimizedProductType = SeqOptmizedProduct<Vector>;

    //! @brief constructor
    //!
    //! @param op the linear operator (assumed double precision)
    //! @param sp the scalar product
    //! @param prec the preconditioner to use
    //! @param reduction the reduction factor passed to the iterative solver
    //! @param maxit maximum number of iterations for the linear solver
    //! @param verbose verbosity level
    //! @param comm the communication object.
    MixedAdapter(Operator *op,
                 std::shared_ptr<AbstractScalarProductType> sp,
                 std::shared_ptr<AbstractPrecondType> prec,
                 const double& tol,
                 const int& maxiter,
                 const int& verbosity,
                 const Comm &comm)
    {

        // Access matrix data from double precision operator
        auto &A = op->getmat();
        int nrows = A.N();
        int nnz = A.nonzeroes();
        double_data_ = &A[0][0][0][0];

        //allocate mixed matrix
        mixed_matrix_ = std::make_shared<MixedMatrixType>(nrows,nnz);
        //auto &B = *mixed_matrix_;

        // copy sparsity pattern from double precision matrix
        int *rows = mixed_matrix_->rowptr();
        int *cols = mixed_matrix_->colidx();

        int irow = 0;
        int icol = 0;
        rows[0]  = 0;
        for(auto row=A.begin(); row!=A.end(); row++)
        {
            for(unsigned int i=0; i<row->getsize(); i++)
            {
                cols[icol++] = row->getindexptr()[i];
            }
            rows[irow+1]     = rows[irow]+row->getsize();
            irow++;
        }

        //initialize mixed operator and optimized scalar product
        double_operator_ = op;
        if constexpr (std::is_same_v<Comm, Dune::Amg::SequentialInformation>)
        {
            mixed_operator_ = std::make_shared<MixedOperatorType>(*mixed_matrix_);
            scalar_product_ = std::make_shared<OptimizedProductType>();
        }
        else
        {
            mixed_operator_ = std::make_shared<MixedOperatorType>(*mixed_matrix_,comm);
            scalar_product_ = sp;
        }

        //initialize bicgstab solver from Dune
        solver_ = std::make_shared<Dune::BiCGSTABSolver<Vector>>(
                                                              *mixed_operator_,
                                                              *scalar_product_,
                                                              *prec,
                                                              tol, // desired residual reduction factor
                                                              maxiter, // maximum number of iterations
                                                              verbosity);
    }

    virtual void apply(Vector &x, Vector &b, InverseOperatorResult &res) override
    {
        //transpose dense blocks and demote to single precision
        mixed_matrix_->update(double_data_);

        //apply bicgstab solver from Dune
        solver_->apply(x,b,res);
    }

    virtual void apply(Vector &x, Vector &b, double reduction, InverseOperatorResult &res) override
    {
        x=0;
        b=0;
        res.reduction = reduction;
        OPM_THROW(std::invalid_argument, "MixedAdapter::apply(...) not implemented yet.");
    }

    virtual Dune::SolverCategory::Category category() const override{return Dune::SolverCategory::overlapping;};

    private:
    using AbstractSolverType = Dune::InverseOperator<Vector,Vector>;


    Operator *double_operator_;
    std::shared_ptr<AbstractSolverType> solver_;
    std::shared_ptr<MixedOperatorType> mixed_operator_;
    std::shared_ptr<MixedMatrixType> mixed_matrix_;
    std::shared_ptr<AbstractScalarProductType> scalar_product_;
    double const *double_data_;

};

}

#endif // OPM_MIXED_ADAPTER_HEADER_INCLUDED
