/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_INCOMPTPFASINGLEPHASE_HEADER_INCLUDED
#define OPM_INCOMPTPFASINGLEPHASE_HEADER_INCLUDED


#include <opm/core/pressure/tpfa/ifs_tpfa.h>
#include <vector>

struct UnstructuredGrid;
struct Wells;

namespace Opm
{

    class IncompPropertiesSinglePhase;
    class LinearSolverInterface;

    /// Encapsulating a tpfa pressure solver for the incompressible-fluid case.
    /// Supports gravity, wells controlled by bhp or reservoir rates,
    /// boundary conditions and simple sources as driving forces.
    /// Rock compressibility can be included, and necessary nonlinear
    /// iterations are handled.
    /// Below we use the shortcuts D for the number of dimensions, N
    /// for the number of cells and F for the number of faces.
    class IncompTpfaSinglePhase
    {
    public:
        /// Construct solver for incompressible case.
        /// \param[in] grid             A 2d or 3d grid.
        /// \param[in] props            Rock and fluid properties.
        /// \param[in] linsolver        Linear solver to use.
        /// \param[in] wells            The wells used as driving forces.
        IncompTpfaSinglePhase(const UnstructuredGrid& grid,
                              const IncompPropertiesSinglePhase& props,
                              const LinearSolverInterface& linsolver,
                              const Wells& wells);

        /// Destructor.
        ~IncompTpfaSinglePhase();

        /// Solve the pressure equation.
        void solve(std::vector<double>& press,
                   std::vector<double>& flux,
                   std::vector<double>& bhp,
                   std::vector<double>& wellrates);

    private:
        // Helper functions.
        void computeStaticData();
        void computePerSolveDynamicData();

    protected:
        // ------ Data that will remain unmodified after construction. ------
        const UnstructuredGrid& grid_;
        const IncompPropertiesSinglePhase& props_;
        const LinearSolverInterface& linsolver_;
        const Wells& wells_;
        std::vector<double> htrans_;
        std::vector<double> trans_ ;
        std::vector<double> zeros_;
        std::vector<double> totmob_;
        struct ifs_tpfa_forces forces_;

        // ------ Internal data for the ifs_tpfa solver. ------
        struct ifs_tpfa_data* h_;
    };

} // namespace Opm

#endif // OPM_INCOMPTPFASINGLEPHASE_HEADER_INCLUDED
