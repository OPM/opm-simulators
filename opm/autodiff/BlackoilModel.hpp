/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2014, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services
  Copyright 2015 NTNU

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

#ifndef OPM_BLACKOILMODEL_HEADER_INCLUDED
#define OPM_BLACKOILMODEL_HEADER_INCLUDED

#include <opm/autodiff/BlackoilModelBase.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>

namespace Opm {

    /// A model implementation for three-phase black oil.
    ///
    /// The simulator is capable of handling three-phase problems
    /// where gas can be dissolved in oil and vice versa. It
    /// uses an industry-standard TPFA discretization with per-phase
    /// upwind weighting of mobilities.
    ///
    /// It uses automatic differentiation via the class AutoDiffBlock
    /// to simplify assembly of the jacobian matrix.
    template<class Grid>
    class BlackoilModel : public BlackoilModelBase<Grid, BlackoilModel<Grid> >
    {
    public:
        typedef BlackoilModelBase<Grid, BlackoilModel<Grid> > Base;

        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] fluid            fluid properties
        /// \param[in] geo              rock properties
        /// \param[in] rock_comp_props  if non-null, rock compressibility properties
        /// \param[in] wells            well structure
        /// \param[in] linsolver        linear solver
        /// \param[in] eclState         eclipse state
        /// \param[in] has_disgas       turn on dissolved gas
        /// \param[in] has_vapoil       turn on vaporized oil feature
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilModel(const typename Base::ModelParameters&   param,
                      const Grid&                             grid,
                      const BlackoilPropsAdInterface&         fluid,
                      const DerivedGeology&                   geo,
                      const RockCompressibility*              rock_comp_props,
                      const Wells*                            wells,
                      const NewtonIterationBlackoilInterface& linsolver,
                      Opm::EclipseStateConstPtr               eclState,
                      const bool                              has_disgas,
                      const bool                              has_vapoil,
                      const bool                              terminal_output)
            : Base(param, grid, fluid, geo, rock_comp_props, wells, linsolver,
                   eclState, has_disgas, has_vapoil, terminal_output)
        {
        }
    };


    /// Providing types by template specialisation of ModelTraits for BlackoilModel.
    template <class Grid>
    struct ModelTraits< BlackoilModel<Grid> >
    {
        typedef BlackoilState ReservoirState;
        typedef WellStateFullyImplicitBlackoil WellState;
        typedef BlackoilModelParameters ModelParameters;
        typedef DefaultBlackoilSolutionState SolutionState;
    };

} // namespace Opm


#endif // OPM_BLACKOILMODEL_HEADER_INCLUDED
