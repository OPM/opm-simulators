#ifndef OPM_FULLYIMPLICITTWOPHASEPOLYMERSOLVER_HEADER_INCLUDED
#define OPM_FULLYIMPLICITTWOPHASEPOLYMERSOLVER_HEADER_INCLUDED

#include <opm/polymer/fullyimplicit/AutoDiffBlock.hpp>
#include <opm/polymer//fullyimplicit/AutoDiffHelpers.hpp>
#include <opm/polymer/fullyimplicit/IncompPropsAdInterface.hpp>
#include <opm/polymer/PolymerProperties.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>
#include <opm/core/pressure/tpfa/trans_tpfa.h>


struct UnstructuredGrid;
struct Wells;
namespace Opm {
    class LinearSolverInterface;
    class PolymerState;
    class WellState;
    
    class FullyImplicitTwophasePolymerSolver
    {
    public:
        FullyImplicitTwophasePolymerSolver(const UnstructuredGrid&        grid,
                                           const IncompPropsAdInterface&  fluid,
                                           const PolymerPropsAd&          polymer_props_ad,
                                           const LinearSolverInterface&    linsolver,
                                           const Wells&                     wells,
                                           const double*                    gravity);

        void step(const double   dt,
                  PolymerState& state,
                  WellState&    well_state,
                  const std::vector<double>& polymer_inflow);
    private:
        typedef AutoDiffBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;
        typedef Eigen::Array<double,
                             Eigen::Dynamic,
                             Eigen::Dynamic,
                             Eigen::RowMajor> DataBlock;

        struct ReservoirResidualQuant {
            ReservoirResidualQuant();
            std::vector<ADB> accum; // Accumulations
            ADB              mflux; // Mass flux (surface conditions)
            ADB              b;     // Reciprocal FVF
            ADB              head;  // Pressure drop across int. interfaces
            ADB              mob;   // Phase mobility (per cell)
        };

        struct SolutionState {
            SolutionState(const int np);
            ADB             pressure;
            std::vector<ADB> saturation;
            ADB             concentration;
            ADB             qs;
            ADB             bhp;
        };

        struct WellOps {
            WellOps(const Wells& wells);
            M w2p;              // well -> perf (scatter)
            M p2w;              // perf -> well (gather)
        };

        const UnstructuredGrid&         grid_;
        const IncompPropsAdInterface&   fluid_;
        const PolymerPropsAd&           polymer_props_ad_;
        const LinearSolverInterface&    linsolver_;
        const Wells&                    wells_;
        const double*                   gravity_;
        const std::vector<int>          cells_;
        HelperOps                       ops_;
        const WellOps                   wops_;
      	V								cmax_;
        std::vector<ReservoirResidualQuant> rq_;
        struct {
            std::vector<ADB>     mass_balance;
            ADB                  well_eq;
            ADB                  well_flux_eq;
        } residual_;

        SolutionState
        constantState(const PolymerState& x,
                      const WellState&    xw);
        SolutionState
        variableState(const PolymerState& x,
                      const WellState&    xw);
        void
        assemble(const V&               pvdt,
                 const SolutionState&   old_state,
                 const PolymerState&    x,
                 const WellState&       xw,
                 const std::vector<double>& polymer_inflow);
        V solveJacobianSystem() const;
        void updateState(const V&             dx,
                         PolymerState& x,
                         WellState&    xw) const;
        std::vector<ADB>
        computeRelPerm(const SolutionState& state) const;
        V
        transmissibility() const;
        
		void
        computeMassFlux(const V&                trans,
                        const ADB&              mc,
                        const ADB&              kro,
                        const ADB&              krw_eff,
                        const SolutionState&    state );
    
        std::vector<ADB>
        accumSource(const ADB&                 kro,
                    const ADB&                 krw_eff,
                    const ADB&                 c,
                    const std::vector<double>& src,
                    const std::vector<double>& polymer_inflow_c) const;

        
        std::vector<ADB>
        computeFracFlow() const;
        double
        residualNorm() const;
        ADB
        polymerSource(const std::vector<ADB>& kr,
                      const std::vector<double>& src,
                      const std::vector<double>& polymer_inflow_c,
                      const SolutionState& state) const;

        void
        computeCmax(PolymerState& state,
					const ADB& c);
    	void
		computeAccum(const SolutionState& state,
                 	const int            aix  );
        ADB 
        computeMc(const SolutionState& state) const;
        ADB
        rockPorosity(const ADB& p) const;
        ADB
        rockPermeability(const ADB& p) const;
        const double
        fluidDensity(const int phase) const;
        ADB
        fluidDensity(const int phase,
                     const ADB p) const;
        ADB
        transMult(const ADB& p) const;
    };


} // namespace Opm
#endif// OPM_FULLYIMPLICITTWOPHASESOLVER_HEADER_INCLUDED
