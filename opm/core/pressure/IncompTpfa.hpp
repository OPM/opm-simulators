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

#ifndef OPM_INCOMPTPFA_HEADER_INCLUDED
#define OPM_INCOMPTPFA_HEADER_INCLUDED


#include <vector>

struct UnstructuredGrid;
struct ifs_tpfa_data;
struct FlowBoundaryConditions;

namespace Opm
{

    class LinearSolverInterface;

    /// Encapsulating a tpfa pressure solver for the incompressible case.
    /// Supports gravity and simple sources as driving forces.
    /// Below we use the shortcuts D for the number of dimensions, N
    /// for the number of cells and F for the number of faces.
    /// Note: we intend to add wells in the future.
    class IncompTpfa
    {
    public:
	/// Construct solver.
	/// \param[in] g             A 2d or 3d grid.
	/// \param[in] permeability  Array of permeability tensors, the array
	///                          should have size N*D^2, if D == g.dimensions
	///                          and N == g.number_of_cells.
	/// \param[in] gravity       Gravity vector. If nonzero, the array should
	///                          have D elements.
	IncompTpfa(const UnstructuredGrid& g,
		   const double* permeability,
		   const double* gravity,
                   const LinearSolverInterface& linsolver);

	/// Destructor.
	~IncompTpfa();

	/// Assemble and solve pressure system.
	/// \param[in]  totmob     Must contain N total mobility values (one per cell).
	///                        totmob = \sum_{p} kr_p/mu_p.
	/// \param[in]  omega      Must be empty if constructor gravity argument was null.
	///                        Otherwise must contain N mobility-weighted density values (one per cell).
	///                        omega = \frac{\sum_{p} mob_p rho_p}{\sum_p rho_p}.
	/// \param[in]  src        Must contain N source rates (one per cell).
	///                        Positive values represent total inflow rates,
	///                        negative values represent total outflow rates.
	/// \param[in]  bcs        If non-null, specifies boundary conditions.
	///                        If null, noflow conditions are assumed.
	/// \param[out] pressure   Will contain N cell-pressure values.
	/// \param[out] faceflux   Will contain F signed face flux values.
	void solve(const std::vector<double>& totmob,
		   const std::vector<double>& omega,
		   const std::vector<double>& src,
		   const FlowBoundaryConditions* bcs,
		   std::vector<double>& pressure,
		   std::vector<double>& faceflux);

        /// Expose read-only reference to internal half-transmissibility.
        const ::std::vector<double>& getHalfTrans() const { return htrans_; }

    private:
	const UnstructuredGrid& grid_;
        const LinearSolverInterface& linsolver_;
	::std::vector<double> htrans_;
	::std::vector<double> trans_ ;
	::std::vector<double> gpress_;
	::std::vector<double> gpress_omegaweighted_;

	struct ifs_tpfa_data* h_;
    };

} // namespace Opm

#endif // OPM_INCOMPTPFA_HEADER_INCLUDED
