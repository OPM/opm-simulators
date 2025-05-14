// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::Restart
 */
#ifndef OPM_RESTART_HPP
#define OPM_RESTART_HPP

#include <dune/geometry/dimension.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Opm {

/*!
 * \brief Load or save a state of a problem to/from the harddisk.
 */
class Restart
{
    /*!
     * \brief Create a magic cookie for restart files, so that it is
     *        unlikely to load a restart file for an incorrectly.
     */
    template <class GridView>
    static std::string magicRestartCookie_(const GridView& gridView)
    {
        static const std::string gridName = "blubb"; // gridView.grid().name();
        static const int dim = GridView::dimension;

        const int numVertices = gridView.size(dim);
        const int numElements = gridView.size(0);
        const int numEdges = gridView.size(dim - 1);
        const int numCPUs = gridView.comm().size();
        const int rank = gridView.comm().rank();

        std::ostringstream oss;
        oss << "eWoms restart file: "
            << "gridName='" << gridName << "' "
            << "numCPUs=" << numCPUs << " "
            << "myRank=" << rank << " "
            << "numElements=" << numElements << " "
            << "numEdges=" << numEdges << " "
            << "numVertices=" << numVertices;
        return oss.str();
    }

    /*!
     * \brief Return the restart file name.
     */
    static std::string restartFileName_(int rank,
                                        const std::string& outputDir,
                                        const std::string& simName,
                                        double t);

public:
    /*!
     * \brief Returns the name of the file which is (de-)serialized.
     */
    const std::string& fileName() const
    { return fileName_; }

    /*!
     * \brief Write the current state of the model to disk.
     */
    template <class Simulator>
    void serializeBegin(Simulator& simulator)
    {
        fileName_ = restartFileName_(simulator.gridView().comm().rank(),
                                     simulator.problem().outputDir(),
                                     simulator.problem().name(),
                                     simulator.time());

        // open output file and write magic cookie
        openOutputStream(magicRestartCookie_(simulator.gridView()));
    }

    /*!
     * \brief The output stream to write the serialized data.
     */
    std::ostream& serializeStream()
    { return outStream_; }

    /*!
     * \brief Start a new section in the serialized output.
     */
    void serializeSectionBegin(const std::string& cookie);

    /*!
     * \brief End of a section in the serialized output.
     */
    void serializeSectionEnd();

    /*!
     * \brief Serialize all leaf entities of a codim in a gridView.
     *
     * The actual work is done by Serializer::serialize(Entity)
     */
    template <int codim, class Serializer, class GridView>
    void serializeEntities(Serializer& serializer, const GridView& gridView)
    {
        serializeSectionBegin("Entities: Codim " + std::to_string(codim));

        for (const auto& entity : entities(gridView, Dune::Codim<codim>())) {
            serializer.serializeEntity(outStream_, entity);
            outStream_ << "\n";
        }

        serializeSectionEnd();
    }

    /*!
     * \brief Finish the restart file.
     */
    void serializeEnd();

    /*!
     * \brief Start reading a restart file at a certain simulated
     *        time.
     */
    template <class Simulator, class Scalar>
    void deserializeBegin(Simulator& simulator, Scalar t)
    {
        fileName_ = restartFileName_(simulator.gridView().comm().rank(),
                                     simulator.problem().outputDir(),
                                     simulator.problem().name(), t);
        openInputStream(magicRestartCookie_(simulator.gridView()));
    }

    /*!
     * \brief The input stream to read the data which ought to be
     *        deserialized.
     */
    std::istream& deserializeStream()
    { return inStream_; }

    /*!
     * \brief Start reading a new section of the restart file.
     */
    void deserializeSectionBegin(const std::string& cookie);

    /*!
     * \brief End of a section in the serialized output.
     */
    void deserializeSectionEnd();

    /*!
     * \brief Deserialize all leaf entities of a codim in a grid.
     *
     * The actual work is done by Deserializer::deserialize(Entity)
     */
    template <int codim, class Deserializer, class GridView>
    void deserializeEntities(Deserializer& deserializer, const GridView& gridView)
    {
        deserializeSectionBegin("Entities: Codim " + std::to_string(codim));

        std::string curLine;

        // read entity data
        for (const auto& entity : entities(gridView, Dune::Codim<codim>())) {
            if (!inStream_.good()) {
                throw std::runtime_error("Restart file is corrupted");
            }

            std::getline(inStream_, curLine);
            std::istringstream curLineStream(curLine);
            deserializer.deserializeEntity(curLineStream, entity);
        }

        deserializeSectionEnd();
    }

    /*!
     * \brief Stop reading the restart file.
     */
    void deserializeEnd();

private:
    /*!
     * \brief Open input stream and check for magic cookie.
     * \param cookie Magic cookie to check for
     */
    void openInputStream(const std::string& cookie);

    /*!
     * \brief Open output stream and write magic cookie.
     * \param cookie Magic cookie to write
     */
    void openOutputStream(const std::string& cookie);

    std::string fileName_;
    std::ifstream inStream_;
    std::ofstream outStream_;
};

} // namespace Opm

#endif // OPM_RESTART_HPP
