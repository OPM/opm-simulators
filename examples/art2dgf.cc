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
#include <opm/material/common/Exceptions.hpp>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <string>
#include <memory>

namespace Ewoms {
/*!
 * \brief Reads in mesh files in the ART format.
 *
 * This file format is used to specify grids with fractures.
 */

 struct Art2DGF
 {
    /*!
     * \brief Create the Grid
     */
    static void convert( const std::string& artFileName,
                         std::ostream& dgfFile,
                         const unsigned precision = 16 )
    {
        using Scalar = double;
        using GlobalPosition = Dune::FieldVector< Scalar, 2 >;
        enum ParseMode { Vertex, Edge, Element, Finished };
        std::vector< std::pair<GlobalPosition, unsigned> > vertexPos;
        std::vector<std::pair<unsigned, unsigned> > edges;
        std::vector<std::pair<unsigned, unsigned> > fractureEdges;
        std::vector<std::vector<unsigned> > elements;
        std::ifstream inStream(artFileName);
        if (!inStream.is_open()) {
            throw std::runtime_error("File '"+artFileName
                                     +"' does not exist or is not readable");
        }
        std::string curLine;
        ParseMode curParseMode = Vertex;
        while (inStream) {
            std::getline(inStream, curLine);

            // remove comments
            auto commentPos = curLine.find("%");
            if (commentPos != curLine.npos) {
                curLine = curLine.substr(0, commentPos);
            }

            // remove leading whitespace
            unsigned numLeadingSpaces = 0;
            while (curLine.size() > numLeadingSpaces
                   && std::isspace(curLine[numLeadingSpaces]))
                ++numLeadingSpaces;
            curLine = curLine.substr(numLeadingSpaces,
                                     curLine.size() - numLeadingSpaces);

            // remove trailing whitespace
            unsigned numTrailingSpaces = 0;
            while (curLine.size() > numTrailingSpaces
                   && std::isspace(curLine[curLine.size() - numTrailingSpaces]))
                ++numTrailingSpaces;
            curLine = curLine.substr(0, curLine.size() - numTrailingSpaces);

            // a section of the file is finished, go to the next one
            if (curLine == "$") {
                if (curParseMode == Vertex)
                    curParseMode = Edge;
                else if (curParseMode == Edge)
                    curParseMode = Element;
                else if (curParseMode == Element)
                    curParseMode = Finished;
                continue;
            }

            // skip empty lines
            if (curLine.empty())
                continue;

            if (curParseMode == Vertex) {
                GlobalPosition coord;
                std::istringstream iss(curLine);
                // parse only the first two numbers as the vertex
                // coordinate. the last number is the Z coordinate
                // which we ignore (so far)
                iss >> coord[0] >> coord[1];
                vertexPos.push_back( std::make_pair( coord, 0 ) );
            }
            else if (curParseMode == Edge) {
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
            else if (curParseMode == Element) {
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
            else if (curParseMode == Finished) {
                assert(curLine.size() == 0);
            }
        }

        dgfFile << "DGF" << std::endl << std::endl;

        dgfFile << "GridParameter" << std::endl
                << "overlap 1" << std::endl
                << "closure green" << std::endl
                << "#" << std::endl << std::endl;

        dgfFile << "Vertex" << std::endl;
        const bool hasFractures = fractureEdges.size() > 0;
        if( hasFractures )
        {
            dgfFile << "parameters 1" << std::endl;
        }
        dgfFile << std::scientific;
        dgfFile.precision( precision );
        const size_t vxSize = vertexPos.size();
        for( size_t i=0; i<vxSize; ++i)
        {
            dgfFile << vertexPos[ i ].first;
            if( hasFractures )
            {
                dgfFile << " " << vertexPos[ i ].second;
            }
            dgfFile << std::endl;
        }

        dgfFile << "#" << std::endl << std::endl;

        dgfFile << "Simplex" << std::endl;
        const size_t elSize = elements.size();
        for( size_t i=0; i<elSize; ++i )
        {
            const size_t elVx = elements[ i ].size();
            for( size_t j=0; j<elVx; ++j )
                dgfFile << elements[ i ][ j ] << " ";
            dgfFile << std::endl;
        }

        dgfFile << "#" << std::endl << std::endl;
        dgfFile << "BoundaryDomain" << std::endl;
        dgfFile << "default 1" << std::endl;
        dgfFile << "#" << std::endl << std::endl;
        dgfFile << "#" << std::endl;
    }
 };

} // namespace Ewoms

int main( int argc, char** argv )
{
    if (argc != 2) {
        std::cout << "Converts a grid file from the ART file format to DGF (Dune grid format)\n"
                  << "\n"
                  << "Usage: " << argv[0] << " ART_FILENAME\n"
                  << "\n"
                  << "The result will be written to the file $ART_FILENAME.dgf\n";
        return 1;
    }

    std::string filename( argv[ 1 ] );
    std::string dgfname( filename );
    dgfname += ".dgf";

    std::cout << "Converting ART file \"" << filename << "\" to DGF file \"" << dgfname << "\"\n";
    std::ofstream dgfFile( dgfname );
    Ewoms::Art2DGF::convert( filename, dgfFile );

    return 0;
}
