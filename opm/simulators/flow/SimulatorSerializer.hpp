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

#ifndef OPM_SIMULATOR_SERIALIZER_HEADER_INCLUDED
#define OPM_SIMULATOR_SERIALIZER_HEADER_INCLUDED

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <array>
#include <string>
#include <vector>

namespace Opm {

class HDF5Serializer;
class IOConfig;
class SimulatorTimer;

//! \brief Abstract interface for simulator serialization ops.
struct SerializableSim {
    //! \brief Load simulator state from file.
    virtual void loadState(HDF5Serializer& serializer,
                           const std::string& groupName) = 0;

    //! \brief Save simulator state to file.
    virtual void saveState(HDF5Serializer& serializer,
                           const std::string& groupName) const = 0;

    //! \brief Get header info to save to file.
    virtual std::array<std::string,5> getHeader() const = 0;

    //! \brief Obtain local-to-global cell mapping.
    virtual const std::vector<int>& getCellMapping() const = 0;
};

//! \brief Class handling simulator serialization.
class SimulatorSerializer {
public:
    //! \brief Constructor inits parameters.
    //! \param saveSpec Specification of steps to save
    //! \param loadStep Step to load
    //! \paramn saveFile File to save to
    //! \param loadFile File to load from
    SimulatorSerializer(SerializableSim& simulator,
                        Parallel::Communication& comm,
                        const IOConfig& ioconfig,
                        const std::string& saveSpec,
                        int loadStep,
                        const std::string& saveFile,
                        const std::string& loadFile);

    //! \brief Returns whether or not a state should be loaded.
    bool shouldLoad() const { return loadStep_ > -1; }

    //! \brief Returns step to load.
    int loadStep() const { return loadStep_; }

    //! \brief Save data to file if appropriate.
    void save(SimulatorTimer& timer);

    //! \brief Loads time step info from file.
    void loadTimerInfo(SimulatorTimer& timer);

    //! \brief Load state from file.
    void loadState();

private:
    //! \brief Checks for differences between command line parameters.
    void checkSerializedCmdLine(const std::string& current,
                                const std::string& stored);

    SerializableSim& simulator_; //!< Reference to simulator to be use
    Parallel::Communication& comm_; //!< Communication to use
    int saveStride_ = 0; //!< Stride to save serialized state at, negative to only keep last
    int saveStep_ = -1; //!< Specific step to save serialized state at
    int loadStep_ = -1; //!< Step to load serialized state from
    std::string saveFile_; //!< File to save serialized state to
    std::string loadFile_; //!< File to load serialized state from
};

} // namespace Opm

#endif // OPM_SIMULATOR_SERIALIZER_HEADER_INCLUDED
