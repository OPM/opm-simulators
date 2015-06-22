/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 STATOIL ASA.

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

#ifndef OPM_FULLYIMPLICITCOMPRESSIBLEPOLYMERSOLVER_HEADER_INCLUDED
#define OPM_FULLYIMPLICITCOMPRESSIBLEPOLYMERSOLVER_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/NewtonIterationBlackoilInterface.hpp>
#include <opm/autodiff/LinearisedBlackoilResidual.hpp>
#include <opm/polymer/PolymerProperties.hpp>
#include <opm/polymer/fullyimplicit/WellStateFullyImplicitBlackoilPolymer.hpp>
#include <opm/polymer/fullyimplicit/PolymerPropsAd.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

struct UnstructuredGrid;
struct Wells;

namespace Opm {

    class DerivedGeology;
    class RockCompressibility;
    class NewtonIterationBlackoilInterface;
    class PolymerBlackoilState;
    class WellStateFullyImplicitBlackoil;

    /// A fully implicit solver for the oil-water with polymer problem.
    ///
    /// The simulator is capable of handling oil-water-polymer problems
    /// It uses an industry-standard TPFA discretization with per-phase
    /// upwind weighting of mobilities.
    ///
    /// It uses automatic differentiation via the class AutoDiffBlock
    /// to simplify assembly of the jacobian matrix.
    class FullyImplicitCompressiblePolymerSolver
    {
    public:
        /// Construct a solver. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] grid             grid data structure
        /// \param[in] fluid            fluid properties
        /// \param[in] geo              rock properties
        /// \param[in] rock_comp_props  if non-null, rock compressibility properties
        /// \param[in] polymer_props_ad polymer properties
        /// \param[in] wells            well structure
        /// \param[in] linsolver        linear solver
        FullyImplicitCompressiblePolymerSolver(const UnstructuredGrid&         grid ,
        		                               const BlackoilPropsAdInterface& fluid,
                   			                   const DerivedGeology&           geo  ,
                              			       const RockCompressibility*      rock_comp_props,
                                    		   const PolymerPropsAd&           polymer_props_ad,
                                    		   const Wells&                    wells,
                                    		   const NewtonIterationBlackoilInterface&    linsolver);

        /// Take a single forward step, modifiying
        ///   state.pressure()
        ///   state.faceflux()
        ///   state.saturation()
        ///   state.concentration()
        ///   wstate.bhp()
        /// \param[in] dt        time step size
        /// \param[in] state     reservoir state
        /// \param[in] wstate    well state
        /// \param[in] polymer_inflow	polymer influx
        int
        step(const double   			dt,
             PolymerBlackoilState& 		state ,
             WellStateFullyImplicitBlackoilPolymer& wstate);

        int newtonIterations() const;
        int linearIterations() const;

        /// Not used by this class except to satisfy interface requirements.
        typedef parameter::ParameterGroup SolverParameters;

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
			std::vector<ADB> ads;   // Adsorption term.
        };

        struct SolutionState {
            SolutionState(const int np);
            ADB              pressure;
            ADB              temperature;
            std::vector<ADB> saturation;
            ADB              concentration;
            ADB              qs;
            ADB              bhp;
        };

        struct WellOps {
            WellOps(const Wells& wells);
            M w2p;              // well -> perf (scatter)
            M p2w;              // perf -> well (gather)
        };

        enum { Water = BlackoilPropsAdInterface::Water,
               Oil   = BlackoilPropsAdInterface::Oil  };

        // Member data
        const UnstructuredGrid&         grid_;
        const BlackoilPropsAdInterface& fluid_;
        const DerivedGeology&           geo_;
        const RockCompressibility*      rock_comp_props_;
        const PolymerPropsAd&           polymer_props_ad_;
        const Wells&                    wells_;
        const NewtonIterationBlackoilInterface&    linsolver_;
        const std::vector<int>          cells_;  // All grid cells
        HelperOps                       ops_;
        const WellOps                   wops_;
        const M                         grav_;
		V    			 				cmax_;
        std::vector<PhasePresence> phaseCondition_;
        std::vector<ReservoirResidualQuant> rq_;
        // The mass_balance vector has one element for each active phase,
        // each of which has size equal to the number of cells.
        // The well_eq has size equal to the number of wells.
        LinearisedBlackoilResidual  residual_;

        unsigned int newtonIterations_;
        unsigned int linearIterations_;

        // Private methods.
        SolutionState
        constantState(const PolymerBlackoilState& x,
                      const WellStateFullyImplicitBlackoil&     xw);

        SolutionState
        variableState(const PolymerBlackoilState& x,
                      const WellStateFullyImplicitBlackoil&     xw);

        void
        computeAccum(const SolutionState& state,
                     const int            aix  );

        void
        assemble(const double             	 dt,
                 const PolymerBlackoilState& x,
                 const WellStateFullyImplicitBlackoil& xw,  
                 const std::vector<double>& polymer_inflow);

        V solveJacobianSystem() const;

        void updateState(const V& dx,
                         PolymerBlackoilState& state,
                         WellStateFullyImplicitBlackoil& well_state) const;

        std::vector<ADB>
        computeRelPerm(const SolutionState& state) const;

        std::vector<ADB>
        computePressures(const SolutionState& state) const;

        void
        computeMassFlux(const int               actph ,
                        const V&                transi,
                        const std::vector<ADB>& kr    ,
                        const SolutionState&    state );

        void
        computeMassFlux(const V&                trans,
                        const ADB&              mc,
                        const ADB&              kro,
                        const ADB&              krw_eff,
                        const SolutionState&    state);

        std::vector<ADB>
        computeFracFlow(const ADB&              kro,
                        const ADB&              krw_eff,
                        const ADB&              c) const;

        void
        computeCmax(PolymerBlackoilState& state);

        ADB
        computeMc(const SolutionState&  state) const;

        ADB
        rockPorosity(const ADB& p) const;

        ADB
        rockPermeability(const ADB& p) const;

        double
        residualNorm() const;

        ADB
        fluidViscosity(const int                         phase,
                       const ADB&                        p    ,
                       const ADB&                        T    ,
                       const std::vector<PhasePresence>& cond,
                       const std::vector<int>&           cells) const;

        ADB
        fluidReciprocFVF(const int                         phase,
                         const ADB&                        p    ,
                         const ADB&                        T    ,
                         const std::vector<PhasePresence>& cond,
                         const std::vector<int>&           cells) const;

        ADB
        fluidDensity(const int                         phase,
                     const ADB&                        p    ,
                     const ADB&                        T    ,
                     const std::vector<PhasePresence>& cond,
                     const std::vector<int>&           cells) const;

        ADB
        poroMult(const ADB& p) const;

        ADB
        transMult(const ADB& p) const;
        
        const std::vector<PhasePresence>
        phaseCondition() const { return phaseCondition_; }
        
        void
        classifyCondition(const PolymerBlackoilState& state);
    };
} // namespace Opm


#endif // OPM_FULLYIMPLICITCOMPRESSIBLEPOLYMERSOLVER_HEADER_INCLUDED
