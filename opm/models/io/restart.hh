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
#ifndef EWOMS_RESTART_HH
#define EWOMS_RESTART_HH

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

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
    static const std::string magicRestartCookie_(const GridView& gridView)
    {
        static const std::string gridName = "blubb"; // gridView.grid().name();
        static const int dim = GridView::dimension;

        int numVertices = gridView.size(dim);
        int numElements = gridView.size(0);
        int numEdges = gridView.size(dim - 1);
        int numCPUs = gridView.comm().size();
        int rank = gridView.comm().rank();

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
    template <class GridView, class Scalar>
    static const std::string restartFileName_(const GridView& gridView,
                                              const std::string& outputDir,
                                              const std::string& simName,
                                              Scalar t)
    {
        std::string dir = outputDir;
        if (dir == ".")
            dir = "";
        else if (!dir.empty() && dir.back() != '/')
            dir += "/";

        int rank = gridView.comm().rank();
        std::ostringstream oss;
        oss << dir << simName << "_time=" << t << "_rank=" << rank << ".ers";
        return oss.str();
    }

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
        const std::string magicCookie = magicRestartCookie_(simulator.gridView());
        fileName_ = restartFileName_(simulator.gridView(),
                                     simulator.problem().outputDir(),
                                     simulator.problem().name(),
                                     simulator.time());

        // open output file and write magic cookie
        outStream_.open(fileName_.c_str());
        outStream_.precision(20);

        serializeSectionBegin(magicCookie);
        serializeSectionEnd();
    }

    /*!
     * \brief The output stream to write the serialized data.
     */
    std::ostream& serializeStream()
    { return outStream_; }

    /*!
     * \brief Start a new section in the serialized output.
     */
    void serializeSectionBegin(const std::string& cookie)
    { outStream_ << cookie << "\n"; }

    /*!
     * \brief End of a section in the serialized output.
     */
    void serializeSectionEnd()
    { outStream_ << "\n"; }

    /*!
     * \brief Serialize all leaf entities of a codim in a gridView.
     *
     * The actual work is done by Serializer::serialize(Entity)
     */
    template <int codim, class Serializer, class GridView>
    void serializeEntities(Serializer& serializer, const GridView& gridView)
    {
        std::ostringstream oss;
        oss << "Entities: Codim " << codim;
        std::string cookie = oss.str();
        serializeSectionBegin(cookie);

        // write element data
        typedef typename GridView::template Codim<codim>::Iterator Iterator;

        Iterator it = gridView.template begin<codim>();
        const Iterator& endIt = gridView.template end<codim>();
        for (; it != endIt; ++it) {
            serializer.serializeEntity(outStream_, *it);
            outStream_ << "\n";
        }

        serializeSectionEnd();
    }

    /*!
     * \brief Finish the restart file.
     */
    void serializeEnd()
    { outStream_.close(); }

    /*!
     * \brief Start reading a restart file at a certain simulated
     *        time.
     */
    template <class Simulator, class Scalar>
    void deserializeBegin(Simulator& simulator, Scalar t)
    {
        fileName_ = restartFileName_(simulator.gridView(), simulator.problem().outputDir(), simulator.problem().name(), t);

        // open input file and read magic cookie
        inStream_.open(fileName_.c_str());
        if (!inStream_.good()) {
            throw std::runtime_error("Restart file '"+fileName_+"' could not be opened properly");
        }

        // make sure that we don't open an empty file
        inStream_.seekg(0, std::ios::end);
        auto pos = inStream_.tellg();
        if (pos == 0) {
            throw std::runtime_error("Restart file '"+fileName_+"' is empty");
        }
        inStream_.seekg(0, std::ios::beg);

        const std::string magicCookie = magicRestartCookie_(simulator.gridView());

        deserializeSectionBegin(magicCookie);
        deserializeSectionEnd();
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
    void deserializeSectionBegin(const std::string& cookie)
    {
        if (!inStream_.good())
            throw std::runtime_error("Encountered unexpected EOF in restart file.");
        std::string buf;
        std::getline(inStream_, buf);
        if (buf != cookie)
            throw std::runtime_error("Could not start section '"+cookie+"'");
    }

    /*!
     * \brief End of a section in the serialized output.
     */
    void deserializeSectionEnd()
    {
        std::string dummy;
        std::getline(inStream_, dummy);
        for (unsigned i = 0; i < dummy.length(); ++i) {
            if (!std::isspace(dummy[i])) {
                throw std::logic_error("Encountered unread values while deserializing");
            }
        }
    }

    /*!
     * \brief Deserialize all leaf entities of a codim in a grid.
     *
     * The actual work is done by Deserializer::deserialize(Entity)
     */
    template <int codim, class Deserializer, class GridView>
    void deserializeEntities(Deserializer& deserializer, const GridView& gridView)
    {
        std::ostringstream oss;
        oss << "Entities: Codim " << codim;
        std::string cookie = oss.str();
        deserializeSectionBegin(cookie);

        std::string curLine;

        // read entity data
        typedef typename GridView::template Codim<codim>::Iterator Iterator;
        Iterator it = gridView.template begin<codim>();
        const Iterator& endIt = gridView.template end<codim>();
        for (; it != endIt; ++it) {
            if (!inStream_.good()) {
                throw std::runtime_error("Restart file is corrupted");
            }

            std::getline(inStream_, curLine);
            std::istringstream curLineStream(curLine);
            deserializer.deserializeEntity(curLineStream, *it);
        }

        deserializeSectionEnd();
    }

    /*!
     * \brief Stop reading the restart file.
     */
    void deserializeEnd()
    { inStream_.close(); }

private:
    std::string fileName_;
    std::ifstream inStream_;
    std::ofstream outStream_;
};
} // namespace Opm

#endif
