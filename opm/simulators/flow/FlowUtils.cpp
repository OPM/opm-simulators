/*
  Copyright 2013, 2014, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2015 IRIS AS
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

#include "opm/simulators/flow/BlackoilModelParameters.hpp"
#include "opm/simulators/flow/FlowProblemParameters.hpp"
#include "opm/simulators/flow/VtkTracerModule.hpp"
#include <config.h>
#include <opm/simulators/flow/FlowUtils.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/String.hpp>

#include <opm/models/common/multiphasebaseparameters.hh>
#include <opm/models/discretization/common/fvbaseparameters.hh>
#include <opm/models/io/vtkblackoilmodule.hpp>
#include <opm/models/io/vtkcompositionmodule.hpp>
#include <opm/models/io/vtkdiffusionmodule.hpp>
#include <opm/models/io/vtkmultiphasemodule.hpp>
#include <opm/models/io/vtkprimaryvarsmodule.hpp>
#include <opm/models/io/vtktemperaturemodule.hpp>
#include <opm/models/nonlinear/newtonmethodparams.hpp>
#include <opm/models/utils/basicparameters.hh>
#include <opm/models/utils/parametersystem.hpp>

#include <opm/simulators/flow/ConvergenceOutputConfiguration.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/utils/ParallelFileMerger.hpp>

#if HAVE_MPI
#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#endif

#include <fmt/format.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace Opm::detail {

void mergeParallelLogFiles(std::string_view output_dir,
                           std::string_view deckFilename,
                           bool enableLoggingFalloutWarning)
{
    namespace fs = ::std::filesystem;
    fs::path output_path(output_dir);
    fs::path deck_filename(deckFilename);
    std::string basename;
    // Strip extension "." and ".DATA"
    std::string extension = uppercase(deck_filename.extension().string());
    if (extension == ".DATA" || extension == ".") {
        basename = uppercase(deck_filename.stem().string());
    } else {
        basename = uppercase(deck_filename.filename().string());
    }
    std::for_each(fs::directory_iterator(output_path),
                  fs::directory_iterator(),
                  detail::ParallelFileMerger(output_path, basename,
                                             enableLoggingFalloutWarning));
}

void handleExtraConvergenceOutput(SimulatorReport& report,
                                  std::string_view option,
                                  std::string_view optionName,
                                  std::string_view output_dir,
                                  std::string_view base_name)
{
    const auto extraConvOutput = ConvergenceOutputConfiguration {
        option, optionName
    };

    if (extraConvOutput.want(ConvergenceOutputConfiguration::Option::Steps)) {
      namespace fs = ::std::filesystem;

      const auto infostep = fs::path{output_dir} / fs::path{base_name}.concat(".INFOSTEP");

      std::ofstream os(infostep);
      report.fullReports(os);
    }
}
void checkAllMPIProcesses()
{
#if HAVE_MPI
    const auto& comm = FlowGenericVanguard::comm();
    if (comm.size() > 1)
    {
        // we try to prevent the abort here.
        // For that we need a signal that each process is here.
        // Each process sends  a message to rank 0.
        const int tag = 357912;
        if (comm.rank() == 0)
        {
            // wait for a message from all processes.
            std::vector<MPI_Request> requests(comm.size() - 1, MPI_REQUEST_NULL);
            std::vector<int> data(comm.size()-1);

            for(decltype(comm.size()) i = 1; i < comm.size(); ++i)
            {
                if (auto error = MPI_Irecv(data.data() + (i - 1), 1, MPI_INT, i, tag, comm, requests.data() + (i - 1));
                    error != MPI_SUCCESS) {
                    OpmLog::error(fmt::format("Error: Could not set up MPI receive (error code : {})", error));
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }
            }
            std::size_t msgs = comm.size() - 1;
            for(std::size_t tries = 0; msgs >0 && tries < 3; ++tries)
            {
                sleep(3);
                int flag, idx;
                for(auto left_msgs = msgs; left_msgs > 0; --left_msgs)
                {
                    if( auto error = MPI_Testany(comm.size()-1, requests.data(), &idx, &flag, MPI_STATUS_IGNORE);
                        error != MPI_SUCCESS) {
                        OpmLog::error(fmt::format("Error: Could not test for MPI message (error code : {})", error));
                        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                    }
                    if (flag)
                    {
                        --msgs;
                    }
                }
            }
            if (msgs) {
                // seems like some processes are stuck. Abort just to be save
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }
        else
        {
            int data= 3;
            MPI_Request request = MPI_REQUEST_NULL;
            if (auto error = MPI_Isend(&data, 1, MPI_INT, 0, tag, comm, &request);
                error != MPI_SUCCESS) {
                OpmLog::error(fmt::format("Error: Could send MPI message (error code : {})", error));
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            bool completed = false;
            for(std::size_t tries = 0; !completed && tries < 3; tries++)
            {
                sleep(3);
                int flag;
                if( auto error = MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
                    error != MPI_SUCCESS) {
                    OpmLog::error(fmt::format("Error: Could not test for MPI message (error code : {})", error));
                    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                }
                if (flag)
                {
                    completed = true;
                }
            }
            if (!completed) {
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
        }
    }
#endif
}

template<class Scalar>
void hideUnusedParameters()
{
    // hide the parameters unused by flow. TODO: this is a pain to maintain
    Parameters::Hide<Parameters::EnableGravity>();
    Parameters::Hide<Parameters::EnableGridAdaptation>();

    // this parameter is actually used in eWoms, but the flow well model
    // hard-codes the assumption that the intensive quantities cache is enabled,
    // so flow crashes. Let's hide the parameter for that reason.
    Parameters::Hide<Parameters::EnableIntensiveQuantityCache>();

    // thermodynamic hints are not implemented/required by the eWoms blackoil
    // model
    Parameters::Hide<Parameters::EnableThermodynamicHints>();

    // in flow only the deck file determines the end time of the simulation
    Parameters::Hide<Parameters::EndTime<Scalar>>();

    // time stepping is not done by the eWoms code in flow
    Parameters::Hide<Parameters::InitialTimeStepSize<Scalar>>();
    Parameters::Hide<Parameters::MaxTimeStepDivisions>();
    Parameters::Hide<Parameters::MaxTimeStepSize<Scalar>>();
    Parameters::Hide<Parameters::MinTimeStepSize<Scalar>>();
    Parameters::Hide<Parameters::PredeterminedTimeStepsFile>();

    // flow also does not use the eWoms Newton method
    Parameters::Hide<Parameters::NewtonMaxError<Scalar>>();
    Parameters::Hide<Parameters::NewtonTolerance<Scalar>>();
    Parameters::Hide<Parameters::NewtonTargetIterations>();
    Parameters::Hide<Parameters::NewtonVerbose>();
    Parameters::Hide<Parameters::NewtonWriteConvergence>();

    // the default eWoms checkpoint/restart mechanism does not work with flow
    Parameters::Hide<Parameters::RestartTime<Scalar>>();
    Parameters::Hide<Parameters::RestartWritingInterval>();
    // hide all vtk related it is not currently possible to do this dependet on if the vtk writing is used
    //if(not(Parameters::Get<Parameters::EnableVtkOutput>())){
        Parameters::Hide<Parameters::VtkWriteOilFormationVolumeFactor>();
        Parameters::Hide<Parameters::VtkWriteOilSaturationPressure>();
        Parameters::Hide<Parameters::VtkWriteOilVaporizationFactor>();
        Parameters::Hide<Parameters::VtkWritePorosity>();
        Parameters::Hide<Parameters::VtkWritePotentialGradients>();
        Parameters::Hide<Parameters::VtkWritePressures>();
        Parameters::Hide<Parameters::VtkWritePrimaryVars>();
        Parameters::Hide<Parameters::VtkWritePrimaryVarsMeaning>();
        Parameters::Hide<Parameters::VtkWriteProcessRank>();
        Parameters::Hide<Parameters::VtkWriteRelativePermeabilities>();
        Parameters::Hide<Parameters::VtkWriteSaturatedGasOilVaporizationFactor>();
        Parameters::Hide<Parameters::VtkWriteSaturatedOilGasDissolutionFactor>();
        Parameters::Hide<Parameters::VtkWriteSaturationRatios>();
        Parameters::Hide<Parameters::VtkWriteSaturations>();
        Parameters::Hide<Parameters::VtkWriteTemperature>();
        Parameters::Hide<Parameters::VtkWriteViscosities>();
        Parameters::Hide<Parameters::VtkWriteWaterFormationVolumeFactor>();
        Parameters::Hide<Parameters::VtkWriteGasDissolutionFactor>();
        Parameters::Hide<Parameters::VtkWriteGasFormationVolumeFactor>();
        Parameters::Hide<Parameters::VtkWriteGasSaturationPressure>();
        Parameters::Hide<Parameters::VtkWriteIntrinsicPermeabilities>();
        Parameters::Hide<Parameters::VtkWriteTracerConcentration>();
        Parameters::Hide<Parameters::VtkWriteExtrusionFactor>();
        Parameters::Hide<Parameters::VtkWriteFilterVelocities>();
        Parameters::Hide<Parameters::VtkWriteDensities>();
        Parameters::Hide<Parameters::VtkWriteDofIndex>();
        Parameters::Hide<Parameters::VtkWriteMobilities>();
        //}
    Parameters::Hide<Parameters::VtkWriteAverageMolarMasses>();
    Parameters::Hide<Parameters::VtkWriteFugacities>();
    Parameters::Hide<Parameters::VtkWriteFugacityCoeffs>();
    Parameters::Hide<Parameters::VtkWriteMassFractions>();
    Parameters::Hide<Parameters::VtkWriteMolarities>();
    Parameters::Hide<Parameters::VtkWriteMoleFractions>();
    Parameters::Hide<Parameters::VtkWriteTotalMassFractions>();
    Parameters::Hide<Parameters::VtkWriteTotalMoleFractions>();

    Parameters::Hide<Parameters::VtkWriteTortuosities>();
    Parameters::Hide<Parameters::VtkWriteDiffusionCoefficients>();
    Parameters::Hide<Parameters::VtkWriteEffectiveDiffusionCoefficients>();

    // hide average density option
    Parameters::Hide<Parameters::UseAverageDensityMsWells>();
}

int eclPositionalParameter(std::function<void(const std::string&, const std::string&)> addKey,
                           std::set<std::string>& seenParams,
                           std::string& errorMsg,
                           const char** argv,
                           int paramIdx)
{
    std::string param  = argv[paramIdx];
    std::size_t i = param.find('=');
    if (i != std::string::npos) {
        std::string oldParamName = param.substr(0, i);
        std::string oldParamValue = param.substr(i+1);
        std::string newParamName = "--" + oldParamName;
        std::replace(newParamName.begin(),
                     newParamName.end(), '_' , '-');
        errorMsg =
          "The old syntax to specify parameters on the command line is no longer supported: "
          "Try replacing '" + oldParamName + "=" + oldParamValue + "' with "+
          "'" + newParamName + "=" + oldParamValue + "'!";
        return 0;
    }

    if (seenParams.count("EclDeckFileName") > 0) {
        errorMsg =
            "Parameter 'EclDeckFileName' specified multiple times"
            " as a command line parameter";
        return 0;
    }

    addKey("EclDeckFileName", argv[paramIdx]);
    seenParams.insert("EclDeckFileName");
    return 1;
}

template void hideUnusedParameters<double>();

#if FLOW_INSTANTIATE_FLOAT
template void hideUnusedParameters<float>();
#endif

} // namespace Opm::detail
