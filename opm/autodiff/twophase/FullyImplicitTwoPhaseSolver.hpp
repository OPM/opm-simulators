/**/

#ifndef OPM_FULLYIMPLICITTWOPHASESOLVER_HEADER_INCLUDED
#define OPM_FULLYIMPLICITTWOPHASESOLVER_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/twophase/IncompPropsAdInterface.hpp>
#include <opm/core/pressure/tpfa/trans_tpfa.h>


struct UnstructuredGrid;
struct Wells;
namespace Opm {
    class LinearSolverInterface;
    class TwophaseState;
    class WellState;
    
    class FullyImplicitTwoPhaseSolver
    {
    public:
        FullyImplicitTwoPhaseSolver(const UnstructuredGrid&        grid,
                                    const IncompPropsAdInterface&  fluid,
                                    const LinearSolverInterface&   linsolver,
                                    const Wells&                   wells,
                                    const double*                  gravity);

        void step(const double   dt,
                  TwophaseState& state,
                  WellState&     well_state);
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
            ADB              pressure;
            std::vector<ADB> saturation;
            ADB              qs;
            ADB              bhp;
        };
        struct WellOps {
            WellOps(const Wells& wells);
            M w2p;              // well -> perf (scatter)
            M p2w;              // perf -> well (gather)
        };
        const UnstructuredGrid&         grid_;
        const IncompPropsAdInterface&   fluid_;
        const LinearSolverInterface&    linsolver_;
        const Wells&                    wells_;
        const double*                   gravity_;
        const std::vector<int>          cells_;
        HelperOps                       ops_;
        const WellOps                   wops_;
        std::vector<ADB>                mob_;
       
        struct {
            std::vector<ADB>    mass_balance;
            ADB                 well_eq;
            ADB                 well_flux_eq;
        } residual_;
       

        SolutionState
        constantState(const TwophaseState& x,
                      const WellState&     xw);
        SolutionState
        variableState(const TwophaseState& x,
                      const WellState&     xw);
        void
        assemble(const V&               pvdt,
                 const SolutionState&   old_state,
                 const TwophaseState&   x,
                 const WellState&       xw);
                // const std::vector<double>& src);
        V solveJacobianSystem() const;
        void updateState(const V&             dx,
                         TwophaseState&       x,
                         WellState&           xw) const;
        std::vector<ADB>
        computeRelPerm(const SolutionState& state) const;
        V
        transmissibility() const;
        ADB
        computeFracFlow(const int phase);
        ADB 
        accumSource(const int phase,
                    const std::vector<ADB>& kr,
                    const std::vector<double>& src) const;
        ADB
        computeMassFlux(const int               phase,
                        const V&                trans,
                        const std::vector<ADB>& kr,
                        const SolutionState&    state);
        double
        residualNorm() const;

        ADB
        fluidDensity(const int phase,
                     const ADB p) const;
        ADB
        rockPorosity(const ADB& p) const;
        ADB
        rockPermeability(const ADB& p) const;
        ADB
        transMult(const ADB& p) const;
    };
} // namespace Opm
#endif// OPM_FULLYIMPLICITTWOPHASESOLVER_HEADER_INCLUDED
