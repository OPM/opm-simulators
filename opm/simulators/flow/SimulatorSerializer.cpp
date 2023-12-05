/*
  Copyright 2023 Equinor ASA.

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

#include <config.h>
#include <opm/simulators/flow/SimulatorSerializer.hpp>

#include <dune/common/hash.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/String.hpp>

#include <opm/input/eclipse/EclipseState/IOConfig/IOConfig.hpp>

#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#if HAVE_HDF5
#include <ebos/hdf5serializer.hh>
#endif

#include <algorithm>
#include <filesystem>
#include <stdexcept>

namespace Opm {

SimulatorSerializer::SimulatorSerializer(SerializableSim& simulator,
                                         Parallel::Communication& comm,
                                         const IOConfig& ioconfig,
                                         const std::string& saveSpec,
                                         int loadStep,
                                         const std::string& saveFile,
                                         const std::string& loadFile)
    : simulator_(simulator)
    , comm_(comm)
    , loadStep_(loadStep)
    , saveFile_(saveFile)
    , loadFile_(loadFile)
{
    if (saveSpec == "all") {
        saveStride_ = 1;
    } else if (saveSpec == "last") {
        saveStride_ = -1;
    } else if (!saveSpec.empty() && saveSpec[0] == ':') {
        saveStride_ = std::atoi(saveSpec.c_str()+1);
    } else if (!saveSpec.empty()) {
        saveStep_ = std::atoi(saveSpec.c_str());
    }

#if !HAVE_HDF5
    if (loadStep_ > -1)  {
        OPM_THROW(std::runtime_error, "Loading of serialized state requested, "
                                      "but no HDF5 support available.");
    }
    if (saveStep_ > -1)  {
        OPM_THROW(std::runtime_error, "Saving of serialized state requested, "
                                      "but no HDF5 support available.");
    }
#endif

    if (loadFile_.empty() || saveFile_.empty()) {
        if (saveFile_.empty()) saveFile_ = ioconfig.fullBasePath() + ".OPMRST";
        if (loadFile_.empty()) loadFile_ = saveFile_;
        if (loadStep_ != -1 && !std::filesystem::exists(loadFile_)) {
            std::filesystem::path path(ioconfig.getInputDir() + "/");
            path.replace_filename(ioconfig.getBaseName() + ".OPMRST");
            loadFile_ = path;
            if (!std::filesystem::exists(loadFile_)) {
                OPM_THROW(std::runtime_error, "Error locating serialized restart file " + loadFile_);
            }
        }
    }
}

void SimulatorSerializer::save(SimulatorTimer& timer)
{
    if (saveStride_ == 0 && saveStep_ == -1) {
        return;
    }

    OPM_BEGIN_PARALLEL_TRY_CATCH();

    int nextStep = timer.currentStepNum();
    if ((saveStep_ != -1 && nextStep == saveStep_)  ||
        (saveStride_ != 0 && (nextStep % saveStride_) == 0)) {
#if HAVE_HDF5
        const std::string groupName = "/report_step/" + std::to_string(nextStep);
        if (saveStride_ < 0 || nextStep == saveStride_ || nextStep == saveStep_) {
            std::filesystem::remove(saveFile_);
        }
        HDF5Serializer writer(saveFile_, HDF5File::OpenMode::APPEND, comm_);
        if (saveStride_ < 0 || nextStep == saveStride_ || nextStep == saveStep_) {
            const auto data = simulator_.getHeader();
            writer.writeHeader(data[0], data[1], data[2], data[3], data[4], comm_.size());

            if (comm_.size() > 1) {
                const auto& cellMapping = simulator_.getCellMapping();
                std::size_t hash = Dune::hash_range(cellMapping.begin(), cellMapping.end());
                writer.write(hash, "/", "grid_checksum");
            }
        }
        simulator_.saveState(writer, groupName);
        writer.write(timer, groupName, "simulator_timer",
                     HDF5File::DataSetMode::ROOT_ONLY);
        OpmLog::info("Serialized state written for report step " + std::to_string(nextStep));
#endif
    }

    OPM_END_PARALLEL_TRY_CATCH("Error saving serialized state: ", comm_);
}

    //! \brief Load timer info from serialized state.
void SimulatorSerializer::loadTimerInfo([[maybe_unused]] SimulatorTimer& timer)
{
#if HAVE_HDF5
    OPM_BEGIN_PARALLEL_TRY_CATCH();

    HDF5Serializer reader(loadFile_, HDF5File::OpenMode::READ, comm_);

    if (loadStep_ == 0) {
        loadStep_ = reader.lastReportStep();
    }

    OpmLog::info("Loading serialized state for report step " + std::to_string(loadStep_));
    const std::string groupName = "/report_step/" + std::to_string(loadStep_);
    reader.read(timer, groupName, "simulator_timer", HDF5File::DataSetMode::ROOT_ONLY);

    std::tuple<std::array<std::string,5>,int> header;
    reader.read(header, "/", "simulator_info", HDF5File::DataSetMode::ROOT_ONLY);
    const auto& [strings, procs] = header;

    if (comm_.size() != procs) {
        throw std::runtime_error("Number of processes (procs=" +
                                 std::to_string(comm_.size()) +
                                 ") does not match .OPMRST file (procs=" +
                                 std::to_string(procs) + ")");
    }

    if (comm_.size() > 1) {
        std::size_t stored_hash;
        reader.read(stored_hash, "/", "grid_checksum");
        const auto& cellMapping = simulator_.getCellMapping();
        std::size_t hash = Dune::hash_range(cellMapping.begin(), cellMapping.end());
        if (hash != stored_hash) {
            throw std::runtime_error("Grid hash mismatch, .OPMRST file cannot be used");
        }
    }

    if (comm_.rank() == 0) {
        const auto& curr_header = simulator_.getHeader();
        checkSerializedCmdLine(curr_header[4], strings[4]);
    }

    OPM_END_PARALLEL_TRY_CATCH("Error loading serialized state: ", comm_);
#endif
}

    //! \brief Load simulator state from serialized state.
void SimulatorSerializer::loadState()
{
#if HAVE_HDF5
    OPM_BEGIN_PARALLEL_TRY_CATCH();

    HDF5Serializer reader(loadFile_, HDF5File::OpenMode::READ, comm_);
    const std::string groupName = "/report_step/" + std::to_string(loadStep_);
    simulator_.loadState(reader, groupName);

    OPM_END_PARALLEL_TRY_CATCH("Error loading serialized state: ", comm_);
#endif
    loadStep_ = -1;
}

void SimulatorSerializer::checkSerializedCmdLine(const std::string& current,
                                                 const std::string& stored)
{
    auto filter_strings = [](const std::vector<std::string>& input)
    {
        std::vector<std::string> output;
        output.reserve(input.size());
        std::copy_if(input.begin(), input.end(), std::back_inserter(output),
                     [](const std::string& line)
                     {
                        return line.compare(0, 11, "EclDeckFile") != 0 &&
                               line.compare(0, 8, "LoadStep") != 0 &&
                               line.compare(0, 9, "OutputDir") != 0 &&
                               line.compare(0, 8, "SaveFile") != 0 &&
                               line.compare(0, 8, "SaveStep") != 0;
                     });
        return output;
    };

    auto curr_strings = split_string(current, '\n');
    auto stored_strings = split_string(stored, '\n');
    std::sort(curr_strings.begin(), curr_strings.end());
    std::sort(stored_strings.begin(), stored_strings.end());
    curr_strings = filter_strings(curr_strings);
    stored_strings = filter_strings(stored_strings);

    std::vector<std::string> difference;
    std::set_symmetric_difference(stored_strings.begin(), stored_strings.end(),
                                  curr_strings.begin(), curr_strings.end(),
                                  std::back_inserter(difference));

    std::vector<std::string> only_stored, only_curr;
    if (!difference.empty()) {
        for (std::size_t i = 0; i < difference.size(); ) {
            auto stored_it = std::find(stored_strings.begin(),
                                       stored_strings.end(), difference[i]);
            auto pos = difference[i].find_first_of('=');
            if (i < difference.size() - 1 &&
                difference[i].compare(0, pos, difference[i+1], 0, pos) == 0) {
                if (stored_it == stored_strings.end()) {
                    std::swap(difference[i], difference[i+1]);
                }
                i += 2;
            } else {
                if (stored_it == stored_strings.end()) {
                    only_curr.push_back(difference[i]);
                } else {
                    only_stored.push_back(difference[i]);
                }
                difference.erase(difference.begin() + i);
            }
        }
        std::stringstream str;
        str << "Differences:\n";
        for (std::size_t i = 0; i < difference.size(); ++i) {
            str << '\t' << (i % 2 == 0 ? '-' : '+') << difference[i] << '\n';
        }
        if (!only_stored.empty()) {
            str << "Only in serialized parameters:\n";
            for (const std::string& line : only_stored)
                str << '\t' << line << '\n';
        }
        if (!only_curr.empty()) {
            str << "Only in current parameters:\n";
            for (const std::string& line : only_curr)
                str << '\t' << line << '\n';
        }
        OpmLog::warning("Command line parameters mismatch:\n" + str.str());
    }
}

} // namespace Opm
