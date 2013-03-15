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

#ifndef OPM_TRANSPORTSOLVERTWOPHASEPOLYMER_HEADER_INCLUDED
#define OPM_TRANSPORTSOLVERTWOPHASEPOLYMER_HEADER_INCLUDED

#include <opm/polymer/PolymerProperties.hpp>
#include <opm/core/transport/reorder/ReorderSolverInterface.hpp>
#include <opm/core/utility/linearInterpolation.hpp>
#include <vector>
#include <list>

class UnstructuredGrid;

namespace Opm
{

    class IncompPropertiesInterface;

    /// Implements a reordering transport solver for incompressible two-phase flow
    /// with polymer in the water phase.
    /// \TODO Include permeability reduction effect.
    class TransportSolverTwophasePolymer : public ReorderSolverInterface
    {
    public:

	enum SingleCellMethod { Bracketing, Newton, Gradient, NewtonSimpleSC, NewtonSimpleC};
        enum GradientMethod { Analytic, FinDif }; // Analytic is chosen (hard-coded)

	/// Construct solver.
	/// \param[in] grid       A 2d or 3d grid.
	/// \param[in] props      Rock and fluid properties.
	/// \param[in] polyprops  Polymer properties.
	/// \param[in] method     Bracketing: solve for c in outer loop, s in inner loop,
        ///                                   each solve being bracketed for robustness.
	///                       Newton: solve simultaneously for c and s with Newton's method.
        ///                               (using gradient variant and bracketing as fallbacks).
	/// \param[in] tol        Tolerance used in the solver.
	/// \param[in] maxit      Maximum number of non-linear iterations used.
	TransportSolverTwophasePolymer(const UnstructuredGrid& grid,
                                       const IncompPropertiesInterface& props,
                                       const PolymerProperties& polyprops,
                                       const SingleCellMethod method,
                                       const double tol,
                                       const int maxit);

	/// Set the preferred method, Bracketing or Newton.
        void setPreferredMethod(SingleCellMethod method);

	/// Solve for saturation, concentration and cmax at next timestep.
	/// Using implicit Euler scheme, reordered.
	/// \param[in] darcyflux           Array of signed face fluxes.
	/// \param[in] porevolume          Array of pore volumes.
	/// \param[in] source              Transport source term, to be interpreted by sign:
        ///                                 (+) Inflow, value is first phase flow (water)
        ///                                     per second, in reservoir volumes.
        ///                                 (-) Outflow, value is total flow of all phases
        ///                                     per second, in reservoir volumes.
	/// \param[in] polymer_inflow_c    Array of inflow polymer concentrations per cell.
	/// \param[in] dt                  Time step.
	/// \param[in, out] saturation     Phase saturations.
	/// \param[in, out] concentration  Polymer concentration.
	/// \param[in, out] cmax           Highest concentration that has occured in a given cell.
	void solve(const double* darcyflux,
                   const double* porevolume,
		   const double* source,
                   const double* polymer_inflow_c,
		   const double dt,
		   std::vector<double>& saturation,
                   std::vector<double>& concentration,
                   std::vector<double>& cmax);

        /// Solve for gravity segregation.
        /// This uses a column-wise nonlinear Gauss-Seidel approach.
        /// It assumes that the input columns contain cells in a single
        /// vertical stack, that do not interact with other columns (for
        /// gravity segregation.
	/// \param[in] columns             Vector of cell-columns.
	/// \param[in] porevolume          Array of pore volumes.
	/// \param[in] dt                  Time step.
	/// \param[in, out] saturation     Phase saturations.
	/// \param[in, out] concentration  Polymer concentration.
	/// \param[in, out] cmax           Highest concentration that has occured in a given cell.
        void solveGravity(const std::vector<std::vector<int> >& columns,
                          const double* porevolume,
                          const double dt,
                          std::vector<double>& saturation,
                          std::vector<double>& concentration,
                          std::vector<double>& cmax);

    public: // But should be made private...
	virtual void solveSingleCell(const int cell);
	virtual void solveMultiCell(const int num_cells, const int* cells);
	void solveSingleCellBracketing(int cell);
	void solveSingleCellNewton(int cell);
	void solveSingleCellGradient(int cell);
	void solveSingleCellNewtonSimple(int cell,bool use_sc);
	class ResidualEquation;

        void initGravity(const double* grav);
        void solveSingleCellGravity(const std::vector<int>& cells,
                                    const int pos,
                                    const double* gravflux);
        int solveGravityColumn(const std::vector<int>& cells);
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


    private:
	const UnstructuredGrid& grid_;
	const double* porosity_;
	const double* porevolume_;  // one volume per cell
	const IncompPropertiesInterface& props_;
	const PolymerProperties& polyprops_;
	std::vector<double> smin_;
	std::vector<double> smax_;
	double tol_;
	double maxit_;

	const double* darcyflux_;   // one flux per grid face
	const double* source_;      // one source per cell
	const double* polymer_inflow_c_;
	double dt_;
        std::vector<double> saturation_; // one per cell, only water saturation!
	double* concentration_;
	double* cmax_;
	std::vector<double> fractionalflow_;  // one per cell
	std::vector<double> mc_;  // one per cell
	const double* visc_;
	SingleCellMethod method_;
	double adhoc_safety_;
	
        // For gravity segregation.
        std::vector<double> gravflux_;
        std::vector<double> mob_;
        std::vector<double> cmax0_;
        // For gravity segregation, column variables
        std::vector<double> s0_;
        std::vector<double> c0_;

	struct ResidualC;
	struct ResidualS;

	class ResidualCGrav;
	class ResidualSGrav;


	void fracFlow(double s, double c, double cmax, int cell, double& ff) const;
	void fracFlowWithDer(double s, double c, double cmax, int cell, double& ff,
                               double* dff_dsdc) const;
	void fracFlowBoth(double s, double c, double cmax, int cell, double& ff,
                          double* dff_dsdc, bool if_with_der) const;
	void computeMc(double c, double& mc) const;
	void computeMcWithDer(double c, double& mc, double& dmc_dc) const;
        void mobility(double s, double c, int cell, double* mob) const;
    };

} // namespace Opm

#endif // OPM_TRANSPORTSOLVERTWOPHASEPOLYMER_HEADER_INCLUDED
