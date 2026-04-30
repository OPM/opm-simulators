#ifndef OPM_MIXED_ADAPTER_HEADER_INCLUDED
#define OPM_MIXED_ADAPTER_HEADER_INCLUDED

//#include <opm/simulators/flow/BlackoilModelParameters.hpp>
//#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/PreconditionerWithUpdate.hpp>
//#include <dune/istl/solver.hh>

namespace Dune
{
template <class Operator, class X>
class MixedAdapter : public InverseOperator<X,X>
{

    public:
    using AbstractPrecondType = Dune::PreconditionerWithUpdate<X, X>;
    using AbstractScalarProductType = Dune::ScalarProduct<X>;


    MixedAdapter(Operator *op,
                 //Dune::ScalarProduct<X> &sp,
                 std::shared_ptr<AbstractScalarProductType> sp,
                 std::shared_ptr<AbstractPrecondType> prec,
                 const double& tol,
                 const int& maxiter,
                 const int& verbosity)
    {
        //operator_ = op;
        //prec_     = prec;
        //product_  = sp;

        //initialize bicgstab solver from Dune
        solver_   = std::make_shared<Dune::BiCGSTABSolver<X>>(*op,
                                                              *sp,
                                                              *prec,
                                                              tol, // desired residual reduction factor
                                                              maxiter, // maximum number of iterations
                                                              verbosity);
    }

    virtual void apply (X& x, X& b, InverseOperatorResult& res) override
    {
        //apply bicgstab solver from Dune
        solver_->apply(x,b,res);
    }

    virtual void apply (X& x, X& b, double reduction, InverseOperatorResult& res) override
    {
        x=0;
        b=0;
        res.reduction = reduction;
        OPM_THROW(std::invalid_argument, "MixedAdapter::apply(...) not implemented yet.");
    }


    private:
    using AbstractSolverType = Dune::InverseOperator<X, X>;

    //Operator* operator_;
    //std::shared_ptr<AbstractPrecondType> prec_;
    //std::shared_ptr<AbstractScalarProductType> product_;
    std::shared_ptr<AbstractSolverType> solver_;

};

}

#endif // OPM_MIXED_ADAPTER_HEADER_INCLUDED
