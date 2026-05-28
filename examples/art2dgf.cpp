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

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <opm/common/utility/String.hpp>

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace Ewoms {
/*!
 * \brief Reads in mesh files in the ART format.
 *
 * This file format is used to specify grids with fractures.
 */

 struct Art2DGF
 {
    using Scalar = double;
    using GlobalPosition = Dune::FieldVector<Scalar, 2>;

    enum class ParseMode { Vertex, Edge, Element, Finished };

    //! \brief Static entry point for conversion.
    //! \param artFileName File name for art file
    //! \param dgfFile Stream to output dgf file to
    //! \param precision Precision for floating point numbers
    static void convert(const std::string& artFileName,
                        std::ostream& dgfFile,
                        const unsigned precision = 16)
    {
        Art2DGF c;
        c.doConvert(artFileName, dgfFile, precision);
    }

private:
    //! \brief Converts a given .art file to a .dgf file
    //! \param artFileName File name for art file
    //! \param dgfFile Stream to output dgf file to
    //! \param precision Precision for floating point numbers
    void doConvert(const std::string& artFileName,
                   std::ostream& dgfFile,
                   const unsigned precision)
    {
        std::ifstream inStream(artFileName);
        if (!inStream.is_open()) {
            throw std::runtime_error("File '"+artFileName
                                     +"' does not exist or is not readable");
        }
        std::string curLine;
        ParseMode curParseMode = ParseMode::Vertex;
        while (inStream) {
            std::getline(inStream, curLine);

            // remove comments
            auto commentPos = curLine.find("%");
            if (commentPos != curLine.npos) {
                curLine.erase(commentPos);
            }

            curLine = Opm::trim_copy(curLine);

            // a section of the file is finished, go to the next one
            if (curLine == "$") {
                if (curParseMode == ParseMode::Vertex)
                    curParseMode = ParseMode::Edge;
                else if (curParseMode == ParseMode::Edge)
                    curParseMode = ParseMode::Element;
                else if (curParseMode == ParseMode::Element)
                    curParseMode = ParseMode::Finished;
                continue;
            }

            // skip empty lines
            if (curLine.empty())
                continue;

            if (curParseMode == ParseMode::Vertex) {
                vertexPos.emplace_back(parseVertex(curLine), 0);
            }
            else if (curParseMode == ParseMode::Edge) {
                // read an edge and update the fracture mapper

                // read the data attached to the edge
                std::istringstream iss(curLine);
                int dataVal;
                std::string tmp;
                iss >> dataVal;
                iss >> tmp;
                assert(tmp == ":");

                // read the vertex indices of an edge
                std::vector<unsigned int> vertIndices;
                while (iss) {
                    unsigned int tmp2;
                    iss >> tmp2;
                    if (!iss)
                        break;
                    vertIndices.push_back(tmp2);
                    assert(tmp2 < vertexPos.size());
                }

                // an edge always has two indices!
                assert(vertIndices.size() == 2);

                std::pair<unsigned, unsigned> edge(vertIndices[0], vertIndices[1]);
                edges.push_back(edge);

                // add the edge to the fracture mapper if it is a fracture
                if (dataVal < 0) {
                    fractureEdges.push_back(edge);
                    vertexPos[ edge.first  ].second = 1;
                    vertexPos[ edge.second ].second = 1;
                }
            }
            else if (curParseMode == ParseMode::Element) {
                // skip the data attached to an element
                std::istringstream iss(curLine);
                int dataVal;
                std::string tmp;
                iss >> dataVal;
                iss >> tmp;
                assert(tmp == ":");

                // read the edge indices of an element
                std::vector<unsigned> edgeIndices;
                while (iss) {
                    unsigned tmp2;
                    iss >> tmp2;
                    if (!iss)
                        break;
                    edgeIndices.push_back(tmp2);
                    assert(tmp2 < edges.size());
                }

                // so far, we only support triangles
                assert(edgeIndices.size() == 3);

                // extract the vertex indices of the element
                std::vector<unsigned> vertIndices;
                for (unsigned i = 0; i < 3; ++i) {
                    bool haveFirstVertex = false;
                    for (unsigned j = 0; j < vertIndices.size(); ++j) {
                        assert(edgeIndices[i] < edges.size());
                        if (vertIndices[j] == edges[edgeIndices[i]].first) {
                            haveFirstVertex = true;
                            break;
                        }
                    }
                    if (!haveFirstVertex)
                        vertIndices.push_back(edges[edgeIndices[i]].first);

                    bool haveSecondVertex = false;
                    for (unsigned j = 0; j < vertIndices.size(); ++j) {
                        assert(edgeIndices[i] < edges.size());
                        if (vertIndices[j] == edges[edgeIndices[i]].second) {
                            haveSecondVertex = true;
                            break;
                        }
                    }
                    if (!haveSecondVertex)
                        vertIndices.push_back(edges[edgeIndices[i]].second);
                }

                // check whether the element's vertices are given in
                // mathematically positive direction. if not, swap the
                // first two.
                Dune::FieldMatrix<Scalar, 2, 2> mat;
                mat[0] = vertexPos[vertIndices[1]].first;
                mat[0] -= vertexPos[vertIndices[0]].first;
                mat[1] = vertexPos[vertIndices[2]].first;
                mat[1] -= vertexPos[vertIndices[0]].first;
                assert(std::abs(mat.determinant()) > 1e-50);
                if (mat.determinant() < 0)
                    std::swap(vertIndices[2], vertIndices[1]);

                elements.push_back( vertIndices );
            }
            else if (curParseMode == ParseMode::Finished) {
                assert(curLine.size() == 0);
            }
        }

        write(dgfFile, precision);
    }

    //! \brief Parses a vertex from a string.
    //! \param curLine String to parse
    static GlobalPosition parseVertex(const std::string& curLine)
    {
        GlobalPosition coord;
        std::istringstream iss(curLine);
        // parse only the first two numbers as the vertex
        // coordinate. the last number is the Z coordinate
        // which we ignore (so far)
        iss >> coord[0] >> coord[1];

        return coord;
    }

    //! \brief Writes out the data to the output stream.
    //! \param dgfFile Stream to write to
    //! \param precision Precision to write at
    void write(std::ostream& dgfFile,
                const unsigned precision)
    {
        dgfFile << "DGF" << std::endl << std::endl;

        dgfFile << "GridParameter" << std::endl
                << "overlap 1" << std::endl
                << "closure green" << std::endl
                << "#" << std::endl << std::endl;

        dgfFile << "Vertex" << std::endl;
        const bool hasFractures = fractureEdges.size() > 0;
        if (hasFractures) {
            dgfFile << "parameters 1" << std::endl;
        }
        dgfFile << std::scientific;
        dgfFile.precision( precision );
        for (const auto& vtx : vertexPos) {
            dgfFile << vtx.first;
            if (hasFractures) {
                dgfFile << " " << vtx.second;
            }
            dgfFile << std::endl;
        }

        dgfFile << "#" << std::endl << std::endl;

        dgfFile << "Simplex" << std::endl;
        for (const auto& element : elements) {
            for (const auto idx : element) {
                dgfFile << idx << " ";
            }
            dgfFile << std::endl;
        }

        dgfFile << "#" << std::endl << std::endl;
        dgfFile << "BoundaryDomain" << std::endl;
        dgfFile << "default 1" << std::endl;
        dgfFile << "#" << std::endl << std::endl;
        dgfFile << "#" << std::endl;
    }

    std::vector< std::pair<GlobalPosition, unsigned> > vertexPos;
    std::vector<std::pair<unsigned, unsigned> > edges;
    std::vector<std::pair<unsigned, unsigned> > fractureEdges;
    std::vector<std::vector<unsigned> > elements;
 };

} // namespace Ewoms

int main( int argc, char** argv )
try {
    if (argc != 2) {
        std::cout << "Converts a grid file from the ART file format to DGF (Dune grid format)\n"
                  << "\n"
                  << "Usage: " << argv[0] << " ART_FILENAME\n"
                  << "\n"
                  << "The result will be written to the file $ART_FILENAME.dgf\n";
        return 1;
    }

    std::string filename(argv[1]);
    std::string dgfname(filename);
    dgfname += ".dgf";

    std::cout << "Converting ART file \"" << filename << "\" to DGF file \"" << dgfname << "\"\n";
    std::ofstream dgfFile(dgfname);
    Ewoms::Art2DGF::convert(filename, dgfFile);

    return 0;
}
catch (const std::exception& e) {
    std::cerr << "Exception thrown " << e.what() << std::endl;
    return 2;
}
catch (...) {
    std::cerr << "Unknown exception thrown" << std::endl;
    return 2;
}
