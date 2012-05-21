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

#ifndef OPM_COMPRESSIBLETPFA_HEADER_INCLUDED
#define OPM_COMPRESSIBLETPFA_HEADER_INCLUDED


#include <vector>

struct UnstructuredGrid;
struct cfs_tpfa_res_data;
struct Wells;
struct FlowBoundaryConditions;

namespace Opm
{

    class BlackoilState;
    class BlackoilPropertiesInterface;
    class LinearSolverInterface;
    class WellState;

    /// Encapsulating a tpfa pressure solver for the compressible-fluid case.
    /// Supports gravity, wells and simple sources as driving forces.
    /// Below we use the shortcuts D for the number of dimensions, N
    /// for the number of cells and F for the number of faces.
    class CompressibleTpfa
    {
    public:
        /// Construct solver.
        /// \param[in] grid          A 2d or 3d grid.
        /// \param[in] props         Rock and fluid properties.
        /// \param[in] linsolver     Linear solver to use.
        /// \param[in] gravity       Gravity vector. If nonzero, the array should
        ///                          have D elements.
        /// \param[in] wells         The wells argument. Will be used in solution, 
        ///                          is ignored if NULL.
        ///                          Note: this class observes the well object, and
        ///                                makes the assumption that the well topology
        ///                                and completions does not change during the
        ///                                run. However, controls (only) are allowed
        ///                                to change.
	CompressibleTpfa(const UnstructuredGrid& grid,
                         const BlackoilPropertiesInterface& props,
                         const LinearSolverInterface& linsolver,
                         const double* gravity,
                         const Wells* wells);

	/// Destructor.
	~CompressibleTpfa();

        void solve(const double dt,
                   BlackoilState& state,
                   WellState& well_state);

    private:
        void computePerSolveDynamicData(const double dt,
                                        const BlackoilState& state,
                                        const WellState& well_state);
        void computeWellPotentials(const BlackoilState& state);
        void computePerIterationDynamicData(const double dt,
                                            const BlackoilState& state,
                                            const WellState& well_state);
        void computeCellDynamicData(const double dt,
                                    const BlackoilState& state,
                                    const WellState& well_state);
        void computeFaceDynamicData(const double dt,
                                    const BlackoilState& state,
                                    const WellState& well_state);
        void computeWellDynamicData(const double dt,
                                    const BlackoilState& state,
                                    const WellState& well_state);
        void assemble(const double dt,
                      const BlackoilState& state,
                      const WellState& well_state);
        void solveIncrement();

	void computeResults(std::vector<double>& pressure,
                            std::vector<double>& faceflux,
                            std::vector<double>& well_bhp,
                            std::vector<double>& well_rate);

        // ------ Data that will remain unmodified after construction. ------
	const UnstructuredGrid& grid_;
        const BlackoilPropertiesInterface& props_;
        const LinearSolverInterface& linsolver_;
        const double* gravity_;
        const Wells* wells_;   // Outside may modify controls (only) between calls to solve().
	std::vector<double> htrans_;
	std::vector<double> trans_ ;
        std::vector<double> porevol_;
        std::vector<int> allcells_;

        // ------ Internal data for the cfs_tpfa_res solver. ------
	struct cfs_tpfa_res_data* h_;

        // ------ Data that will be modified for every solve. ------
        std::vector<double> wellperf_gpot_;

        // ------ Data that will be modified for every solver iteration. ------
        // Gravity and capillary contributions (per face).
        std::vector<double> cell_A_;
        std::vector<double> cell_dA_;
        std::vector<double> cell_viscosity_;
        std::vector<double> cell_phasemob_;
        std::vector<double> cell_voldisc_;
        std::vector<double> face_A_;
        std::vector<double> face_phasemob_;
        std::vector<double> face_gravcap_;
        std::vector<double> wellperf_A_;
        std::vector<double> wellperf_phasemob_;
        // The update to be applied to the pressures (cell and bhp).
        std::vector<double> pressure_increment_;





    };

} // namespace Opm


#endif // OPM_COMPRESSIBLETPFA_HEADER_INCLUDED
