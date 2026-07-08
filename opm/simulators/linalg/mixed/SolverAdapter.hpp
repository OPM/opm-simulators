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


template<class Vector, class Comm>
class GhostLastScalarProduct : public ScalarProduct<Vector>
{
    public:

    static constexpr auto block_size = Vector::block_type::dimension;

    //! \brief The type of the vector to compute the scalar product on.
    //!
    //! E.g. BlockVector or another type fulfilling the ISTL
    //! vector interface.

    /*!
    * \param com The communication object for syncing overlap and copy
    * data points.
    * \param cat parallel solver category (nonoverlapping or overlapping)
    */
    GhostLastScalarProduct (std::shared_ptr<const Comm> com, SolverCategory::Category cat)
        : _communication(com), _category(cat)
    {
        count_ = getLocalCount();
        int verify = verifyLocalCount();
        if (count_ != verify) OPM_THROW(std::runtime_error, "Inconsistent local node count!!\n");
    }

    /*!
    * \param com The communication object for syncing overlap and copy
    * data points.
    * \param cat parallel solver category (nonoverlapping or overlapping)
    * \note if you use this constructor you have to make sure com stays alive
    */
    GhostLastScalarProduct (const Comm& com, SolverCategory::Category cat)
        : GhostLastScalarProduct(stackobject_to_shared_ptr(com), cat)
    {}

    /*! \brief Dot product of two vectors.
    It is assumed that the vectors are consistent on the interior+border
    partition.
    */
    void setLocalCount(int count){ count_ = count; }

    virtual double dot (const Vector& vx, const Vector& vy) const override
    {

        // access underlying data
        double const *x = &vx[0][0];
        double const *y = &vy[0][0];

        // total array length
        int NN = block_size*count_;

        //double result(0);

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

        // Global summation
        auto cc = _communication->communicator();
        double result = cc.sum(agg[0]);
        return result;
    }

    /*! \brief Norm of a right-hand side vector.
    The vector must be consistent on the interior+border partition
    */
    virtual double norm (const Vector& x) const override
    {
        return sqrt(dot(x,x));
    }

    //! Category of the scalar product (see SolverCategory::Category)
    virtual SolverCategory::Category category() const override
    {
        return _category;
    }

    private:
    std::shared_ptr<const Comm> _communication;
    SolverCategory::Category _category;
    int count_;

    //int getLocalCount(const Comm& comm) const
    int  getLocalCount() const
    {
        int count = 0;
        // Loop over index set
        auto indexSet = _communication->indexSet();
        for (auto idx = indexSet.begin(); idx!=indexSet.end(); ++idx) {
            if (idx->local().attribute()==1) count++; // count non-local indices
        }
        return count;
    }

    int verifyLocalCount() const
    {
        auto indexSet = _communication->indexSet();

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

};



/*
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
    //using type = Dune::AssembledLinearOperator<Matrix,Vector,Vector>;
    //using type = Opm::WellModelMatrixAdapter<Matrix, Vector, Vector>;
};
*/

//! @brief Wraps mixed precision
//!
//! @tparam Comm the communicator passed to FlexibleLinearSolver
//! @tparam Operator the linear operator passed to FlexibleLinearSolver
//! @tparam Vector the block vector passed to FlexibleLinearSolver
template <class Comm, class Operator, class Vector>
class MixedBiCGSTABSolver:public InverseOperator<Vector, Vector>
{
    public:

    using AbstractPrecondType = Dune::PreconditionerWithUpdate<Vector,Vector>;
    using AbstractScalarProductType = Dune::ScalarProduct<Vector>;

    static constexpr auto block_size = Vector::block_type::dimension;
    using MixedMatrixType      = Opm::MixedMatrixWrapper<Vector, block_size>;
    //using MixedOperatorType    = MixedOperator<MixedMatrixType, Vector, Comm>::type;

    //! @brief constructor
    //!
    //! @param op the linear operator (assumed double precision)
    //! @param sp the scalar product
    //! @param prec the preconditioner to use
    //! @param reduction the reduction factor passed to the iterative solver
    //! @param maxit maximum number of iterations for the linear solver
    //! @param verbose verbosity level
    //! @param comm the communication object.
    MixedBiCGSTABSolver(Operator *op,
                 std::shared_ptr<AbstractScalarProductType> sp,
                 std::shared_ptr<AbstractPrecondType> prec,
                 const double& tol,
                 const int& maxiter,
                 const int& verbosity,
                 const Comm &comm)
    {
#if 1
        int halo;
        int nrows;
        int nnz=0;

        auto &A = op->getmat();
        if constexpr (std::is_same_v<Comm, Dune::Amg::SequentialInformation>)
        {
            halo = 0;
            nrows = A.N();
            nnz = A.nonzeroes();
        }
        else
        {
            local_ = new int[A.N()];

            // number of ghost cells
            halo = getHaloCount(comm);

            // number of local cells
            nrows = A.N() - halo;

            // number of nonzeros for local cells
            int irow=0;
            //for(auto row=A.begin(); row!=A.end(); row++)
            for(auto row=A.begin(); row.index() < nrows; row++)
            {
                if(local_[irow++]==1) for(auto col = row->begin(); col != row->end(); col++) nnz++;
                //else break; // This line is used to verify that all local nodes have indices lower than all ghost nodes. This appears to be true!
                //irow++;
            }
        }

        printf("nnz           = %d\n",nnz);
        printf("A.nonzeroes() = %ld\n",A.nonzeroes());

        printf("local = %d\n",nrows);
        printf("halo  = %d\n",halo);
        printf("total = %ld\n",A.N());
        //getchar();
#endif
        // Access matrix data from double precision operator
#if 0
        auto &A = op->getmat();
        int nrows = A.N();
        int nnz = A.nonzeroes();
#endif
        double_data_ = &A[0][0][0][0];

        //allocate mixed matrix
        mixed_matrix_ = std::make_shared<MixedMatrixType>(nrows,nnz);

        // copy sparsity pattern from double precision matrix
        int *rows = mixed_matrix_->rowptr();
        int *cols = mixed_matrix_->colidx();

        int irow = 0;
        int icol = 0;
        rows[0]  = 0;
        for(auto row=A.begin(); row.index() < nrows; row++)
        {
            for(auto col = row->begin(); col != row->end(); ++col)
            {
                cols[icol++] = col.index();
            }
            rows[irow+1] = icol;
            irow++;
        }

        // The following skeleton is in preparation for better support for various MatrixAdapter. For now, it simply throws an error if
        // the operator provided is not supported.
        using MatrixType = std::remove_const_t<std::remove_reference_t<decltype(op->getmat())>>;
        if constexpr (std::is_same_v<std::remove_pointer_t<Operator>, Dune::MatrixAdapter<MatrixType, Vector, Vector>>)
        {
        }
        else if constexpr (std::is_same_v<std::remove_pointer_t<Operator>, Opm::GhostLastMatrixAdapter<MatrixType, Vector, Vector, Comm>>)
        {
        }
        else
        {
            OPM_THROW(std::invalid_argument, "MixedBiCGSTABSolver only supports Dune::MatrixAdapter and Opm::GhostLastMatrixAdapter\n");
        }

        //initialize mixed operator and optimized scalar product
        double_operator_ = op;
        using MatrixType = std::remove_const_t<std::remove_reference_t<decltype(op->getmat())>>;
        if constexpr (std::is_same_v<std::remove_pointer_t<Operator>, Dune::MatrixAdapter<MatrixType, Vector, Vector>>)
        {
            //OPM_THROW(std::invalid_argument, "Dune::MatrixAdapter\n");
            using MixedOperatorType = Dune::MatrixAdapter<MixedMatrixType, Vector, Vector>;
            mixed_operator_ = std::make_shared<MixedOperatorType>(*mixed_matrix_);
            using OptimizedProductType = SeqOptmizedProduct<Vector>;
            scalar_product_ = std::make_shared<OptimizedProductType>();
        }
        else if constexpr (std::is_same_v<std::remove_pointer_t<Operator>, Opm::GhostLastMatrixAdapter<MatrixType, Vector, Vector, Comm>>)
        {
            //OPM_THROW(std::invalid_argument, "Opm::GhostLastMatrixAdapter\n");
            //using MixedOperatorType = Dune::OverlappingSchwarzOperator<MixedMatrixType, Vector, Vector, Comm>;
            using MixedOperatorType = Opm::MixedGhostLastMatrixAdapter<MixedMatrixType, Vector, Comm>;
            mixed_operator_ = std::make_shared<MixedOperatorType>(*mixed_matrix_,comm);

            using OptimizedScalarProductType = GhostLastScalarProduct<Vector,Comm>;
            scalar_product_ = std::make_shared<OptimizedScalarProductType>(comm,Dune::SolverCategory::overlapping);
            //scalar_product_ = sp;
        }
        else if constexpr (std::is_same_v<std::remove_pointer_t<Operator>, Opm::WellModelMatrixAdapter<MatrixType, Vector, Vector>>)
        {
            //OPM_THROW(std::invalid_argument, "Opm::WellModelMatrixAdapter\n");
            using MixedOperatorType = Opm::WellModelMatrixAdapter<MixedMatrixType, Vector, Vector>;
            using WellOperatorType  = Opm::LinearOperatorExtra<Vector,Vector>;
            const WellOperatorType &wellOper = op->getwellOper();
            mixed_operator_ = std::make_shared<MixedOperatorType>(*mixed_matrix_, wellOper);
            using OptimizedProductType = SeqOptmizedProduct<Vector>;
            scalar_product_ = std::make_shared<OptimizedProductType>();
            //scalar_product_ = sp;
        }
        else if constexpr (std::is_same_v<std::remove_pointer_t<Operator>, Opm::WellModelGhostLastMatrixAdapter<MatrixType, Vector, Vector, true>>)
        {
            //OPM_THROW(std::invalid_argument, "Opm::WellModelGhostLastMatrixAdapter\n");
            using MixedOperatorType = Opm::WellModelMixedGhostLastMatrixAdapter<MixedMatrixType, Vector, Comm>;
            using WellOperatorType  = Opm::LinearOperatorExtra<Vector,Vector>;
            const WellOperatorType &wellOper = op->getwellOper();
            mixed_operator_ = std::make_shared<MixedOperatorType>(*mixed_matrix_, wellOper, comm);

            if constexpr (std::is_same_v<Comm, Dune::Amg::SequentialInformation>)
            {
                scalar_product_ = sp;
            }
            else
            {
                using OptimizedScalarProductType = GhostLastScalarProduct<Vector,Comm>;
                scalar_product_ = std::make_shared<OptimizedScalarProductType>(comm,Dune::SolverCategory::overlapping);
            }
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
        OPM_THROW(std::invalid_argument, "MixedBiCGSTABSolver::apply(...) not implemented yet.");
    }

    virtual Dune::SolverCategory::Category category() const override{return Dune::SolverCategory::overlapping;};

    private:

    int getHaloCount(const Comm& comm) const
    {
        int count = 0;
        // Loop over index set
        auto indexSet = comm.indexSet();
        for (auto idx = indexSet.begin(); idx!=indexSet.end(); ++idx)
        {
            if (idx->local().attribute()!=1) count++; // count ghost indices
/*
            // This code snippet is used to check wheter or not all ghost indices occur last in the index set
            // In general, that is NOT the case
            if (idx->local().attribute()==1) count++; // count local indices
            else
            {
                printf("debug: %d: %ld %d\n",count, idx->local().local(), idx->local().attribute());
                break;
            }
*/
            int i=idx->local().local(); // tag local indices
            local_[i] = (idx->local().attribute()==1) ? 1 : 0;
        }

        return count;
    }

    using AbstractSolverType   = Dune::InverseOperator<Vector,Vector>;
    using AbstractOperatorType = Dune::AssembledLinearOperator<MixedMatrixType,Vector,Vector>;

    Operator *double_operator_;
    std::shared_ptr<AbstractSolverType> solver_;
    std::shared_ptr<AbstractOperatorType> mixed_operator_;
    std::shared_ptr<MixedMatrixType> mixed_matrix_;
    std::shared_ptr<AbstractScalarProductType> scalar_product_;
    double const *double_data_;

    int *local_;
};

}

#endif // OPM_MIXED_ADAPTER_HEADER_INCLUDED
