/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_TRANSPORTSOLVERTWOPHASECOMPRESSIBLEPOLYMER_HEADER_INCLUDED
#define OPM_TRANSPORTSOLVERTWOPHASECOMPRESSIBLEPOLYMER_HEADER_INCLUDED

#include <opm/polymer/PolymerProperties.hpp>
#include <opm/core/transport/reorder/ReorderSolverInterface.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <vector>
#include <list>

struct UnstructuredGrid;

namespace {
    class ResSOnCurve;
    class ResCOnCurve;
}

namespace Opm
{

    class BlackoilPropertiesInterface;

    /// Implements a reordering transport solver for incompressible two-phase flow
    /// with polymer in the water phase.
    /// \TODO Include permeability reduction effect.
    class TransportSolverTwophaseCompressiblePolymer : public ReorderSolverInterface
    {
    public:

	enum SingleCellMethod { Bracketing, Newton, NewtonC, Gradient};
        enum GradientMethod { Analytic, FinDif }; // Analytic is chosen (hard-coded)

	/// Construct solver.
	/// \param[in] grid       A 2d or 3d grid.
	/// \param[in] props      Rock and fluid properties.
	/// \param[in] polyprops  Polymer properties.
        /// \param[in] rock_comp  Rock compressibility properties
	/// \param[in] method     Bracketing: solve for c in outer loop, s in inner loop,
        ///                                   each solve being bracketed for robustness.
	///                       Newton: solve simultaneously for c and s with Newton's method.
        ///                               (using gradient variant and bracketing as fallbacks).
	/// \param[in] tol        Tolerance used in the solver.
	/// \param[in] maxit      Maximum number of non-linear iterations used.
	TransportSolverTwophaseCompressiblePolymer(const UnstructuredGrid& grid,
                                          const BlackoilPropertiesInterface& props,
                                          const PolymerProperties& polyprops,
                                          const SingleCellMethod method,
                                          const double tol,
                                          const int maxit);

	/// Set the preferred method, Bracketing or Newton.
        void setPreferredMethod(SingleCellMethod method);

	/// Solve for saturation, concentration and cmax at next timestep.
	/// Using implicit Euler scheme, reordered.
	/// \param[in] darcyflux           Array of signed face fluxes.
	/// \param[in] initial_pressure    Array with pressure at start of timestep.
	/// \param[in] pressure            Array with pressure.
	/// \param[in] porevolume0         Array with pore volume at start of timestep.
	/// \param[in] porevolume          Array with pore volume.
	/// \param[in] source              Transport source term, to be interpreted by sign:
        ///                                 (+) Inflow, value is first phase flow (water)
        ///                                     per second, in *surface* volumes (unlike the
        ///                                     incompressible version).
        ///                                 (-) Outflow, value is total flow of all phases
        ///                                     per second, in reservoir volumes.
	/// \param[in] polymer_inflow_c    Array of inflow polymer concentrations per cell.
	/// \param[in] dt                  Time step.
	/// \param[in, out] saturation     Phase saturations.
	/// \param[in, out] surfacevol     Surface volumes.
	/// \param[in, out] concentration  Polymer concentration.
	/// \param[in, out] cmax           Highest concentration that has occured in a given cell.
	void solve(const double* darcyflux,
                   const std::vector<double>& initial_pressure,
                   const std::vector<double>& pressure,
                   const double* porevolume0,
                   const double* porevolume,
		   const double* source,
		   const double* polymer_inflow_c,
		   const double dt,
		   std::vector<double>& saturation,
		   std::vector<double>& surfacevol,
                   std::vector<double>& concentration,
                   std::vector<double>& cmax);

        /// Initialise quantities needed by gravity solver.
        /// \param[in] grav    Gravity vector
        void initGravity(const double* grav);

        /// Solve for gravity segregation.
        /// This uses a column-wise nonlinear Gauss-Seidel approach.
        /// It assumes that the input columns contain cells in a single
        /// vertical stack, that do not interact with other columns (for
        /// gravity segregation.
	/// \param[in] columns             Vector of cell-columns.
	/// \param[in] dt                  Time step.
	/// \param[in, out] saturation     Phase saturations.
	/// \param[in, out] surfacevol     Surface volumes.
	/// \param[in, out] concentration  Polymer concentration.
	/// \param[in, out] cmax           Highest concentration that has occured in a given cell.
        void solveGravity(const std::vector<std::vector<int> >& columns,
                          const double dt,
                          std::vector<double>& saturation,
                          std::vector<double>& surfacevol,
                          std::vector<double>& concentration,
                          std::vector<double>& cmax);

        



    private: 

	const UnstructuredGrid& grid_;
	const BlackoilPropertiesInterface& props_;
	const PolymerProperties& polyprops_;
	const double* darcyflux_;   // one flux per grid face
        const double* porevolume0_; // one volume per cell
        const double* porevolume_;  // one volume per cell
	const double* source_;      // one source per cell
	const double* polymer_inflow_c_;
	double dt_;
	double tol_;
	double maxit_;
	SingleCellMethod method_;
	double adhoc_safety_;

        std::vector<double> saturation_; // one per cell, only water saturation!
        std::vector<int> allcells_;
	double* concentration_;
	double* cmax_;
	std::vector<double> fractionalflow_;  // one per cell
	std::vector<double> mc_;  // one per cell
        std::vector<double> visc_; // viscosity (without polymer, for given pressure)
        std::vector<double> A_;
        std::vector<double> A0_;
	std::vector<double> smin_;
	std::vector<double> smax_;
	
        // For gravity segregation.
        const double* gravity_;
        std::vector<double> trans_;
        std::vector<double> density_;
        std::vector<double> gravflux_;
        std::vector<double> mob_;
        std::vector<double> cmax0_;

        // For gravity segregation, column variables
        std::vector<double> s0_;
        std::vector<double> c0_;

        // Storing the upwind and downwind graphs for experiments.
        std::vector<int> ia_upw_;
        std::vector<int> ja_upw_;
        std::vector<int> ia_downw_;
        std::vector<int> ja_downw_;
        
	struct ResidualC;
	struct ResidualS;

	class ResidualCGrav;
	class ResidualSGrav;

        class ResidualEquation;
        class ResSOnCurve;
        class ResCOnCurve;

	friend class TransportSolverTwophaseCompressiblePolymer::ResidualEquation;
        friend class TransportSolverTwophaseCompressiblePolymer::ResSOnCurve;
        friend class TransportSolverTwophaseCompressiblePolymer::ResCOnCurve;


	virtual void solveSingleCell(const int cell);
	virtual void solveMultiCell(const int num_cells, const int* cells);
	void solveSingleCellBracketing(int cell);
	void solveSingleCellNewton(int cell, bool use_sc, bool use_explicit_step = false);
	void solveSingleCellGradient(int cell);
        void solveSingleCellGravity(const std::vector<int>& cells,
                                    const int pos,
                                    const double* gravflux);
        int solveGravityColumn(const std::vector<int>& cells);

        void initGravityDynamic();

	void fracFlow(double s, double c, double cmax, int cell, double& ff) const;
	void fracFlowWithDer(double s, double c, double cmax, int cell, double& ff,
                               double* dff_dsdc) const;
	void fracFlowBoth(double s, double c, double cmax, int cell, double& ff,
                          double* dff_dsdc, bool if_with_der) const;
	void computeMc(double c, double& mc) const;
	void computeMcWithDer(double c, double& mc, double& dmc_dc) const;
        void mobility(double s, double c, int cell, double* mob) const;
        void scToc(const double* x, double* x_c) const;
        #ifdef PROFILING
        class Newton_Iter {
        public:
            bool res_s;
            int cell;
            double s;
            double c;

            Newton_Iter(bool res_s_val, int cell_val, double s_val, double c_val) {
                res_s = res_s_val;
                cell = cell_val;
                s = s_val;
                c = c_val;
            }
        };

        std::list<Newton_Iter> res_counts;
        #endif
    };

} // namespace Opm

#endif // OPM_TRANSPORTSOLVERTWOPHASECOMPRESSIBLEPOLYMER_HEADER_INCLUDED
