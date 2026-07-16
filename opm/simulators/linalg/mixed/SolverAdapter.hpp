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
#include <opm/simulators/linalg/mixed/Operators.hpp>
#include <opm/simulators/linalg/ScalarProducts.hpp>


namespace Dune
{

//! @brief Adapts BiCGSTAB to mixed precision
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
    using MixedMatrixType      = Opm::MixedMatrixWrapper<Vector>;

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
        int halo;
        int nrows;
        int nnz=0;

        auto &A = op->getmat();
        // trivially determine size of halo==0 for serial linear operators
        if constexpr (std::is_same_v<Comm, Dune::Amg::SequentialInformation>)
        {
            halo = 0;
            nrows = A.N();
            nnz = A.nonzeroes();
        }
        // Determine size of halo for parallel linear operators
        else
        {
            local_ = new int[A.N()];

            // number of ghost cells
            halo = getHaloCount(comm);

            // number of local cells
            nrows = A.N() - halo;

            // number of nonzeros for local cells
            int irow=0;
            for(auto row=A.begin(); row.index() < nrows; row++)
            {
                if(local_[irow++]==1) for(auto col = row->begin(); col != row->end(); col++) nnz++;
            }
        }
/*
        printf("nnz           = %d\n",nnz);
        printf("A.nonzeroes() = %ld\n",A.nonzeroes());

        printf("local = %d\n",nrows);
        printf("halo  = %d\n",halo);
        printf("total = %ld\n",A.N());
        //getchar();
*/
        // Access matrix data from double precision operator
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

        // initialize mixed operator and scalar product depending on the linear operator type provided to the constructor
        double_operator_ = op;
        using MatrixType = std::remove_const_t<std::remove_reference_t<decltype(op->getmat())>>;

        // serial runs with plain block-sparse matrices, i.e. Dune::MatrixAdapter
        if constexpr (std::is_same_v<std::remove_pointer_t<Operator>, Dune::MatrixAdapter<MatrixType, Vector, Vector>>)
        {
            using MixedOperatorType = Dune::MatrixAdapter<MixedMatrixType, Vector, Vector>;
            mixed_operator_ = std::make_shared<MixedOperatorType>(*mixed_matrix_);

            using ScalarProductType = SeqOptmizedProduct<Vector>;
            scalar_product_ = std::make_shared<ScalarProductType>();
        }
        // serial runs with separate linear operator for wells, i.e. Opm::WellModelMatrixAdapter
        else if constexpr (std::is_same_v<std::remove_pointer_t<Operator>, Opm::WellModelMatrixAdapter<MatrixType, Vector, Vector>>)
        {
            using MixedOperatorType = Opm::WellModelMatrixAdapter<MixedMatrixType, Vector, Vector>;
            using WellOperatorType  = Opm::LinearOperatorExtra<Vector,Vector>;
            const WellOperatorType &wellOper = op->getwellOper();
            mixed_operator_ = std::make_shared<MixedOperatorType>(*mixed_matrix_, wellOper);

            using ScalarProductType = SeqOptmizedProduct<Vector>;
            scalar_product_ = std::make_shared<ScalarProductType>();
        }
        // parallel runs with plain block-sparse matrices and all ghost cells sorted after local cells, i.e. Opm::GhostLastMatrixAdapter
        else if constexpr (std::is_same_v<std::remove_pointer_t<Operator>, Opm::GhostLastMatrixAdapter<MatrixType, Vector, Vector, Comm>>)
        {
            using MixedOperatorType = Opm::MixedGhostLastMatrixAdapter<MixedMatrixType, Vector, Comm>;
            mixed_operator_ = std::make_shared<MixedOperatorType>(*mixed_matrix_,comm);

            using ScalarProductType = GhostLastScalarProduct<Vector,Comm>;
            scalar_product_ = std::make_shared<ScalarProductType>(comm,Dune::SolverCategory::overlapping);
        }
        // parallel runs with separate linear operators for wells  and all ghost cells sorted after local cells, i.e. Opm::WellModelGhostLastMatrixAdapter
        else if constexpr (std::is_same_v<std::remove_pointer_t<Operator>, Opm::WellModelGhostLastMatrixAdapter<MatrixType, Vector, Vector, true>>)
        {
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
                using ScalarProductType = GhostLastScalarProduct<Vector,Comm>;
                scalar_product_ = std::make_shared<ScalarProductType>(comm,Dune::SolverCategory::overlapping);
            }
        }
        // throw an exception for all other linear operator types
        else { OPM_THROW(std::invalid_argument, "MixedBiCGSTABSolver: Unsupported linear operator type!!\n");}

        //initialize bicgstab solver from Dune
        solver_ = std::make_shared<Dune::BiCGSTABSolver<Vector>>(
                                                              *mixed_operator_,
                                                              *scalar_product_,
                                                              *prec,
                                                              tol, // desired residual reduction factor
                                                              maxiter, // maximum number of iterations
                                                              verbosity);

    }

    //! @brief destructor
    ~MixedBiCGSTABSolver()
    {
        if constexpr (std::is_same_v<Comm, Dune::Amg::SequentialInformation>) return;
        delete [] local_;
    }

    //! @brief Solver application
    virtual void apply(Vector &x, Vector &b, InverseOperatorResult &res) override
    {
        //transpose dense blocks and demote to single precision
        mixed_matrix_->update(double_data_);

        //apply bicgstab solver from Dune
        solver_->apply(x,b,res);
    }

    //! @brief Unused variant of solver application
    virtual void apply(Vector &x, Vector &b, double reduction, InverseOperatorResult &res) override
    {
        x=0;
        b=0;
        res.reduction = reduction;
        OPM_THROW(std::invalid_argument, "MixedBiCGSTABSolver::apply(...) not implemented yet.");
    }

    //! @brief Solver category
    virtual Dune::SolverCategory::Category category() const override{return Dune::SolverCategory::overlapping;};

    private:

    //! @brief Count number of ghost cells
    //!
    //! @param comm communicator object
    int getHaloCount(const Comm& comm) const
    {
        int count = 0;
        // Loop over index set
        auto indexSet = comm.indexSet();
        for (auto idx = indexSet.begin(); idx!=indexSet.end(); ++idx)
        {
            if (idx->local().attribute()!=1) count++; // count ghost indices

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
