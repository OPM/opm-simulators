/**/

#ifndef OPM_FULLYIMPLICITTWOPHASESOLVER_HEADER_INCLUDED
#define OPM_FULLYIMPLICITTWOPHASESOLVER_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/polymer/IncompPropsAdInterface.hpp>
#include <opm/core/pressure/tpfa/trans_tpfa.h>


struct UnstructuredGrid;
namespace Opm {
//    struct HelperOps;
    class LinearSolverInterface;
    class TwophaseState;

    
    class FullyImplicitTwoPhaseSolver
    {
    public:
        FullyImplicitTwoPhaseSolver(const UnstructuredGrid&        grid,
                                    const IncompPropsAdInterface&  fluid,
                                    const LinearSolverInterface&    linsolver);

        void step(const double   dt,
                  TwophaseState& state,
                  const std::vector<double>& src);
    private:
        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;
        typedef Eigen::Array<double,
                             Eigen::Dynamic,
                             Eigen::Dynamic,
                             Eigen::RowMajor> DataBlock;
        struct SolutionState {
            SolutionState(const int np);
            ADB             pressure;
            std::vector<ADB> saturation;
        };
        const UnstructuredGrid&         grid_;
        const IncompPropsAdInterface&   fluid_;
        const LinearSolverInterface&    linsolver_;
        const std::vector<int>          cells_;
        HelperOps                       ops_;
        std::vector<ADB>                residual_;
       

        SolutionState
        constantState(const TwophaseState& x);
        SolutionState
        variableState(const TwophaseState& x);
        void
        assemble(const V&               pvdt,
                 const SolutionState&   old_state,
                 const TwophaseState&  x,
                 const std::vector<double>& src);
        V solveJacobianSystem() const;
        void updateState(const V&             dx,
                         TwophaseState& x) const;
        std::vector<ADB>
        computeRelPerm(const SolutionState& state) const;
        V
        transmissibility() const;
        ADB
        computeFracFlow(int    phase,
                        const std::vector<ADB>& kr);
        ADB
        computeMassFlux(const int               phase,
                        const V&                trans,
                        const std::vector<ADB>& kr,
                        const SolutionState&    state);
        double
        residualNorm() const;

        ADB
        rockPorosity(const ADB& p) const;
        ADB
        rockPermeability(const ADB& p) const;
        const double
        fluidDensity(const int phase) const;
        ADB
        transMult(const ADB& p) const;
    };
} // namespace Opm
#endif// OPM_FULLYIMPLICITTWOPHASESOLVER_HEADER_INCLUDED
