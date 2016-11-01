/*
  Copyright (c) 2014 SINTEF ICT, Applied Mathematics.
  Copyright (c) 2015-2016 IRIS AS

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
#include "config.h"

#include "SimulatorFullyImplicitBlackoilOutputEbos.hpp"

#include <opm/common/data/SimulationDataContainer.hpp>

#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/utility/DataMap.hpp>
#include <opm/autodiff/Compat.hpp>
#include <opm/output/vtk/writeVtkData.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/utility/miscUtilities.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>

#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/BackupRestore.hpp>

#include <sstream>
#include <iomanip>
#include <fstream>

#include <boost/filesystem.hpp>

//For OutputWriterHelper
#include <map>
#include <opm/parser/eclipse/Units/UnitSystem.hpp>


#ifdef HAVE_OPM_GRID
#include <opm/common/utility/platform_dependent/disable_warnings.h>
#include <dune/common/version.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <opm/common/utility/platform_dependent/reenable_warnings.h>
#endif
namespace Opm
{
    namespace detail {

        struct WriterCallEbos : public ThreadHandle :: ObjectInterface
        {
            BlackoilOutputWriterEbos& writer_;
            std::unique_ptr< SimulatorTimerInterface > timer_;
            const SimulationDataContainer state_;
            const WellState wellState_;
            data::Solution simProps_;
            const bool substep_;

            explicit WriterCallEbos( BlackoilOutputWriterEbos& writer,
                                 const SimulatorTimerInterface& timer,
                                 const SimulationDataContainer& state,
                                 const WellState& wellState,
                                 const data::Solution& simProps,
                                 bool substep )
                : writer_( writer ),
                  timer_( timer.clone() ),
                  state_( state ),
                  wellState_( wellState ),
                  simProps_( simProps ),
                  substep_( substep )
            {
            }

            // callback to writer's serial writeTimeStep method
            void run ()
            {
                // write data
                writer_.writeTimeStepSerial( *timer_, state_, wellState_, simProps_, substep_ );
            }
        };
    }

    void
    BlackoilOutputWriterEbos::
    writeTimeStepWithoutCellProperties(
                  const SimulatorTimerInterface& timer,
                  const SimulationDataContainer& localState,
                  const WellState& localWellState,
                  bool substep)
    {
        data::Solution noCellProperties;
        writeTimeStepWithCellProperties(timer, localState, localWellState, noCellProperties, substep);
    }





    void
    BlackoilOutputWriterEbos::
    writeTimeStepWithCellProperties(
                  const SimulatorTimerInterface& timer,
                  const SimulationDataContainer& localState,
                  const WellState& localWellState,
                  const data::Solution& sol,
                  bool substep)
    {
        bool isIORank = output_ ;
        if( parallelOutput_ && parallelOutput_->isParallel() )
        {
            // If this is not the initial write and no substep, then the well
            // state used in the computation is actually the one of the last
            // step. We need that well state for the gathering. Otherwise
            // It an exception with a message like "global state does not
            // contain well ..." might be thrown.
            int wellStateStepNumber = ( ! substep && timer.reportStepNum() > 0) ?
                (timer.reportStepNum() - 1) : timer.reportStepNum();
            // collect all solutions to I/O rank
            isIORank = parallelOutput_->collectToIORank( localState, localWellState, wellStateStepNumber );
        }

        const SimulationDataContainer& state = (parallelOutput_ && parallelOutput_->isParallel() ) ? parallelOutput_->globalReservoirState() : localState;
        const WellState& wellState  = (parallelOutput_ && parallelOutput_->isParallel() ) ? parallelOutput_->globalWellState() : localWellState;

        // serial output is only done on I/O rank
        if( isIORank )
        {
            if( asyncOutput_ ) {
                // dispatch the write call to the extra thread
                asyncOutput_->dispatch( detail::WriterCallEbos( *this, timer, state, wellState, sol, substep ) );
            }
            else {
                // just write the data to disk
                writeTimeStepSerial( timer, state, wellState, sol, substep );
            }
        }
    }



    void
    BlackoilOutputWriterEbos::
    writeTimeStepSerial(const SimulatorTimerInterface& timer,
                        const SimulationDataContainer& state,
                        const WellState& wellState,
                        const data::Solution& sol,
                        bool substep)
    {
        // ECL output
        if ( eclWriter_ )
        {
            const auto& initConfig = eclipseState_.getInitConfig();
            if (initConfig.restartRequested() && ((initConfig.getRestartStep()) == (timer.currentStepNum()))) {
                std::cout << "Skipping restart write in start of step " << timer.currentStepNum() << std::endl;
            } else {
                eclWriter_->writeTimeStep(timer.reportStepNum(),
                                          substep,
                                          timer.simulationTimeElapsed(),
                                          simToSolution( state, phaseUsage_ ),
                                          wellState.report(phaseUsage_));
            }
        }

        // write backup file
        if( backupfile_.is_open() )
        {
            int reportStep      = timer.reportStepNum();
            int currentTimeStep = timer.currentStepNum();
            if( (reportStep == currentTimeStep || // true for SimulatorTimer
                 currentTimeStep == 0 || // true for AdaptiveSimulatorTimer at reportStep
                 timer.done() ) // true for AdaptiveSimulatorTimer at reportStep
               && lastBackupReportStep_ != reportStep ) // only backup report step once
            {
                // store report step
                lastBackupReportStep_ = reportStep;
                // write resport step number
                backupfile_.write( (const char *) &reportStep, sizeof(int) );

                try {
                    backupfile_ << state;

                    const WellStateFullyImplicitBlackoil& boWellState = static_cast< const WellStateFullyImplicitBlackoil& > (wellState);
                    backupfile_ << boWellState;
                }
                catch ( const std::bad_cast& e )
                {
                }

                backupfile_ << std::flush;
            }
        } // end backup
    }

    void
    BlackoilOutputWriterEbos::
    restore(SimulatorTimerInterface& timer,
            BlackoilState& state,
            WellStateFullyImplicitBlackoilDense& wellState,
            const std::string& filename,
            const int desiredResportStep )
    {
        std::ifstream restorefile( filename.c_str() );
        if( restorefile )
        {
            std::cout << "============================================================================"<<std::endl;
            std::cout << "Restoring from ";
            if( desiredResportStep < 0 ) {
                std::cout << "last";
            }
            else {
                std::cout << desiredResportStep;
            }
            std::cout << " report step! filename = " << filename << std::endl << std::endl;

            int reportStep;
            restorefile.read( (char *) &reportStep, sizeof(int) );

            const int readReportStep = (desiredResportStep < 0) ?
                std::numeric_limits<int>::max() : desiredResportStep;

            while( reportStep <= readReportStep && ! timer.done() && restorefile )
            {
                restorefile >> state;
                restorefile >> wellState;

                // No per cell data is written for restore steps, but will be
                // for subsequent steps, when we have started simulating
                writeTimeStepWithoutCellProperties( timer, state, wellState );

                // some output
                std::cout << "Restored step " << timer.reportStepNum() << " at day "
                          <<  unit::convert::to(timer.simulationTimeElapsed(),unit::day) << std::endl;

                if( readReportStep == reportStep ) {
                    break;
                }

                // if the stream is not valid anymore we just use the last state read
                if( ! restorefile ) {
                    std::cerr << "Reached EOF, using last state read!" << std::endl;
                    break;
                }

                // try to read next report step
                restorefile.read( (char *) &reportStep, sizeof(int) );

                // if read failed, exit loop
                if( ! restorefile ) {
                    break;
                }

                // next step
                timer.advance();

                if( timer.reportStepNum() != reportStep ) {
                    break;
                }
            }
        }
        else
        {
            std::cerr << "Warning: Couldn't open restore file '" << filename << "'" << std::endl;
        }
    }


    bool BlackoilOutputWriterEbos::isRestart() const {
        const auto& initconfig = eclipseState_.getInitConfig();
        return initconfig.restartRequested();
    }
}
