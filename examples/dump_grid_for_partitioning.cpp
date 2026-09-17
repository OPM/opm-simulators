/*
  Copyright 2025 Equinor ASA

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
#include <opm/grid/CpGrid.hpp>

#include <cstddef>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>

#include <opm/grid/common/ZoltanGraphFunctions.hpp>
#include <dune/common/parallel/mpihelper.hh>


size_t getNumberOfEdges(const Dune::CpGrid& grid) {
    size_t totalEdges = 0;
    for(int i = 0; i < grid.numCells();  i++ ) {
        totalEdges += Dune::cpgrid::getNumberOfEdgesForSpecificCell(grid, i);
    }

    return totalEdges;
}

void createMetisGraph(const Dune::CpGrid& grid, std::ostream& os) {
     
    os << grid.numCells() << " " << getNumberOfEdges(grid)/2 << "\n";

    for (int cell = 0; cell < grid.numCells(); ++cell) {
        for (int localFace = 0 ; localFace < grid.numCellFaces(cell); ++localFace ) {
            const int face  = grid.cellFace(cell, localFace);
            int otherCell   = grid.faceCell(face, 0);
            if (otherCell == cell || otherCell == -1) {
                otherCell = grid.faceCell(face, 1);
                if (otherCell == cell || otherCell == -1) {
                    continue;
                }
            }
            // In Metis, the vertices are 1-based, so we add 1
            os << otherCell + 1 << " ";
        }
        os << std::endl;
    }
}
    


int
main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc,argv);

    
    if (argc < 3) {
        std::cerr << "A simple program to dump a grid in a format suitable for partitioning with METIS." << std::endl;
        std::cerr << "The output file can be paritioned with the command:" << std::endl << std::endl;
        std::cerr << "  $ gpmetis <outputfile> <num_partitions>" << std::endl << std::endl;;
        std::cerr << "The resulting partitioning can be used with flow by specifying" << std::endl << std::endl;
        std::cerr << "   --external-partition=<outputfile>.part.<num_partitions>" << std::endl << std::endl;
        std::cerr << "Usage: " << argv[0] << " <input deck file> <outputfile>" << std::endl;
        return 1;
    }

    std::string inputFilename = argv[1];
    std::string outputFilename = argv[2];

    {
        std::ofstream outputFile(outputFilename);
        if (!outputFile) {
            std::cerr << "Error opening output file: " << outputFilename << std::endl;
            return 1;
        }
    }
    Opm::Parser parser;
    const auto deck = parser.parseFile(inputFilename);

    Dune::CpGrid grid;
    Opm::EclipseGrid eclGrid(deck);

    // TODO: Handle transmissibilities and wells if needed
    grid.processEclipseFormat(&eclGrid, nullptr, false, false, false);

    std::ofstream outputFile(outputFilename);
    createMetisGraph(grid, outputFile);
}