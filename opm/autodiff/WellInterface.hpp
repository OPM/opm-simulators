/*
  Copyright 2017 SINTEF ICT, Applied Mathematics.
  Copyright 2017 Statoil ASA.

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


#ifndef OPM_WELLINTERFACE_HEADER_INCLUDED
#define OPM_WELLINTERFACE_HEADER_INCLUDED

#include <opm/common/OpmLog/OpmLog.hpp>


#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/core/wells.h>
#include <opm/core/well_controls.h>
#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/autodiff/VFPProperties.hpp>
#include <opm/autodiff/VFPInjProperties.hpp>
#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/autodiff/WellHelpers.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoilDense.hpp>
#include <opm/autodiff/BlackoilModelParameters.hpp>

#include <string>
#include <memory>
#include <vector>
#include <cassert>

namespace Opm
{


    template<typename TypeTag>
    class WellInterface
    {
    public:

        using WellState = WellStateFullyImplicitBlackoilDense;

        typedef BlackoilModelParameters ModelParameters;
        typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typedef typename GET_PROP_TYPE(TypeTag, Indices) BlackoilIndices;
        typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

        static const int solventCompIdx = 3; //TODO get this from ebos
        static const bool has_solvent = GET_PROP_VALUE(TypeTag, EnableSolvent);

        /// Constructor
        WellInterface(const Well* well, const int time_step, const Wells* wells);

        /// Well name.
        const std::string& name() const;

        /// The index of the well in Wells struct
        // It is used to locate the inforation in Wells and also WellState for now.
        int indexOfWell() const;

        /// Well type, INJECTOR or PRODUCER.
        WellType wellType() const;

        /// number of phases
        int numberOfPhases() const;

        /// Component fractions for each phase for the well
        const std::vector<double>& compFrac() const;

        /// Well controls
        // TODO: later to see whether we need to return const.
        WellControls* wellControls() const;

        /// Number of the perforations
        int numberOfPerforations() const;

        /// Well productivity index for each perforation.
        const std::vector<double>& wellIndex() const;

        /// Depth of perforations
        const std::vector<double>& perfDepth() const;

        /// Indices of the grid cells/blocks that perforations are completed within.
        const std::vector<int>& wellCells() const;

        // TODO: the following function should be able to be removed by refactoring the well class
        // It is probably only needed for StandardWell
        /// the densities of the fluid  in each perforation
        virtual const std::vector<double>& perfDensities() const = 0;
        virtual std::vector<double>& perfDensities() = 0;

        /// the pressure difference between different perforations
        virtual const std::vector<double>& perfPressureDiffs() const = 0;
        virtual std::vector<double>& perfPressureDiffs() = 0;

        // TODO: the parameters need to be optimized/adjusted
        void init(const PhaseUsage* phase_usage_arg,
                  const std::vector<bool>* active_arg,
                  const VFPProperties* vfp_properties_arg,
                  const double gravity_arg,
                  const int num_cells);

        // TODO: temporary
        virtual void setWellVariables(const WellState& well_state) = 0;

        const std::vector<bool>& active() const;

        const PhaseUsage& phaseUsage() const;

        int flowPhaseToEbosCompIdx( const int phaseIdx ) const;

        int flowToEbosPvIdx( const int flowPv ) const;

        int flowPhaseToEbosPhaseIdx( const int phaseIdx ) const;

        int numPhases() const;

        int numComponents() const;

        bool allowCrossFlow() const;

        // TODO: for this kind of function, maybe can make a function with parameter perf
        const std::vector<int>& saturationTableNumber() const;

        const double wsolvent() const;

        virtual bool getWellConvergence(Simulator& ebosSimulator,
                                        std::vector<double>& B_avg,
                                        const ModelParameters& param) const = 0;

    protected:
        // TODO: some variables shared by all the wells should be made static
        // well name
        std::string name_;

        // the index of well in Wells struct
        int index_of_well_;

        // well type
        // INJECTOR or PRODUCER
        enum WellType well_type_;

        // whether the well allows crossflow
        bool allow_cf_;

        // number of phases
        int number_of_phases_;

        // component fractions for each well
        // typically, it should apply to injection wells
        std::vector<double> comp_frac_;

        // controls for this well
        // TODO: later will check whehter to let it stay with pointer
        struct WellControls* well_controls_;

        // number of the perforations for this well
        int number_of_perforations_;

        // record the index of the first perforation
        // TODO: it might not be needed if we refactor WellState to be a vector
        // of states of individual well.
        int first_perf_;

        // well index for each perforation
        std::vector<double> well_index_;

        // depth for each perforation
        std::vector<double> perf_depth_;

        // reference depth for the BHP
        int ref_depth_;

        double well_efficiency_factor_;

        // cell index for each well perforation
        std::vector<int> well_cell_;

        // saturation table nubmer for each well perforation
        std::vector<int> saturation_table_number_;

        const PhaseUsage* phase_usage_;

        const std::vector<bool>* active_;

        const VFPProperties* vfp_properties_;

        double gravity_;
    };

}

#include "WellInterface_impl.hpp"

#endif // OPM_WELLINTERFACE_HEADER_INCLUDED
