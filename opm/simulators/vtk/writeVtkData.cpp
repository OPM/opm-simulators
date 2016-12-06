/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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
#include <opm/simulators/vtk/writeVtkData.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/core/grid.h>
#include <set>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>



namespace Opm
{

    void writeVtkData(const std::array<int, 3>& dims,
                      const std::array<double, 3>& cell_size,
                      const std::map< std::string, const std::vector< double >* >& data,
                      std::ostream& os)
    {
        // Dimension is hardcoded in the prototype and the next two lines,
        // but the rest is flexible (allows dimension == 2 or 3).
        int dimension = 3;
        int num_cells = dims[0]*dims[1]*dims[2];

        assert(dimension == 2 || dimension == 3);
        assert(num_cells == dims[0]*dims[1]* (dimension == 2 ? 1 : dims[2]));

        os << "# vtk DataFile Version 2.0\n";
        os << "Structured Grid\n \n";
        os << "ASCII \n";
        os << "DATASET STRUCTURED_POINTS\n";

        os << "DIMENSIONS "
                 << dims[0] + 1 << " "
                 << dims[1] + 1 << " ";
        if (dimension == 3) {
            os << dims[2] + 1;
        } else {
            os << 1;
        }
        os << "\n";

        os << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << "\n";

        os << "SPACING " << cell_size[0] << " " << cell_size[1];
        if (dimension == 3) {
            os << " " << cell_size[2];
        } else {
            os << " " << 0.0;
        }
        os << "\n";

        os << "\nCELL_DATA " << num_cells << '\n';
        for (auto dit = data.begin(); dit != data.end(); ++dit) {
            std::string name = dit->first;
            os << "SCALARS " << name << " float" << '\n';
            os << "LOOKUP_TABLE " << name << "_table " << '\n';
            const std::vector<double>& field = *(dit->second);
            // We always print only the first data item for every
            // cell, using 'stride'.
            // This is a hack to get water saturation nicely.
            // \TODO: Extend to properly printing vector data.
            const int stride = field.size()/num_cells;
            const int num_per_line = 5;
            for (int c = 0; c < num_cells; ++c) {
                os << field[stride*c] << ' ';
                if (c % num_per_line == num_per_line - 1
                    || c == num_cells - 1) {
                    os << '\n';
                }
            }
        }
    }

    typedef std::map<std::string, std::string> PMap;


    struct Tag
    {
        Tag(const std::string& tag, const PMap& props, std::ostream& os)
            : name_(tag), os_(os)
        {
            indent(os);
            os << "<" << tag;
            for (PMap::const_iterator it = props.begin(); it != props.end(); ++it) {
                os << " " << it->first << "=\"" << it->second << "\"";
            }
            os << ">\n";
            ++indent_;
        }
        Tag(const std::string& tag, std::ostream& os)
            : name_(tag), os_(os)
        {
            indent(os);
            os << "<" << tag << ">\n";
            ++indent_;
        }
        ~Tag()
        {
            --indent_;
            indent(os_);
            os_ << "</" << name_ << ">\n";
        }
        static void indent(std::ostream& os)
        {
            for (int i = 0; i < indent_; ++i) {
                os << "  ";
            }
        }
    private:
        static int indent_;
        std::string name_;
        std::ostream& os_;
    };

    int Tag::indent_ = 0;


    void writeVtkData(const UnstructuredGrid& grid,
                      const std::map< std::string, const std::vector< double >* >& data,
                      std::ostream& os)
    {
       if (grid.dimensions != 3) {
           OPM_THROW(std::runtime_error, "Vtk output for 3d grids only");
       }
       os.precision(12);
       os << "<?xml version=\"1.0\"?>\n";
       PMap pm;
       pm["type"] = "UnstructuredGrid";
       Tag vtkfiletag("VTKFile", pm, os);
       Tag ugtag("UnstructuredGrid", os);
       int num_pts = grid.number_of_nodes;
       int num_cells = grid.number_of_cells;
       pm.clear();
       pm["NumberOfPoints"] = std::to_string(num_pts);
       pm["NumberOfCells"] = std::to_string(num_cells);
       Tag piecetag("Piece", pm, os);
       {
           Tag pointstag("Points", os);
           pm.clear();
           pm["type"] = "Float64";
           pm["Name"] = "Coordinates";
           pm["NumberOfComponents"] = "3";
           pm["format"] = "ascii";
           Tag datag("DataArray", pm, os);
           for (int i = 0; i < num_pts; ++i) {
               Tag::indent(os);
               os << grid.node_coordinates[3*i + 0] << ' '
                  << grid.node_coordinates[3*i + 1] << ' '
                  << grid.node_coordinates[3*i + 2] << '\n';
           }
       }
       {
           Tag cellstag("Cells", os);
           pm.clear();
           pm["type"] = "Int32";
           pm["NumberOfComponents"] = "1";
           pm["format"] = "ascii";
           std::vector<int> cell_numpts;
           cell_numpts.reserve(num_cells);
           {
               pm["Name"] = "connectivity";
               Tag t("DataArray", pm, os);
               int hf = 0;
               for (int c = 0; c < num_cells; ++c) {
                   std::set<int> cell_pts;
                   for (; hf < grid.cell_facepos[c+1]; ++hf) {
                       int f = grid.cell_faces[hf];
                       const int* fnbeg = grid.face_nodes + grid.face_nodepos[f];
                       const int* fnend = grid.face_nodes + grid.face_nodepos[f+1];
                       cell_pts.insert(fnbeg, fnend);
                   }
                   cell_numpts.push_back(cell_pts.size());
                   Tag::indent(os);
                   std::copy(cell_pts.begin(), cell_pts.end(),
                             std::ostream_iterator<int>(os, " "));
                   os << '\n';
               }
           }
           {
               pm["Name"] = "offsets";
               Tag t("DataArray", pm, os);
               int offset = 0;
               const int num_per_line = 10;
               for (int c = 0; c < num_cells; ++c) {
                   if (c % num_per_line == 0) {
                       Tag::indent(os);
                   }
                   offset += cell_numpts[c];
                   os << offset << ' ';
                   if (c % num_per_line == num_per_line - 1
                       || c == num_cells - 1) {
                       os << '\n';
                   }
               }
           }
           std::vector<int> cell_foffsets;
           cell_foffsets.reserve(num_cells);
           {
               pm["Name"] = "faces";
               Tag t("DataArray", pm, os);
               const int* fp = grid.cell_facepos;
               int offset = 0;
               for (int c = 0; c < num_cells; ++c) {
                   Tag::indent(os);
                   os << fp[c+1] - fp[c] << '\n';
                   ++offset;
                   for (int hf = fp[c]; hf < fp[c+1]; ++hf) {
                       int f = grid.cell_faces[hf];
                       const int* np = grid.face_nodepos;
                       int f_num_pts = np[f+1] - np[f];
                       Tag::indent(os);
                       os << f_num_pts << ' ';
                       ++offset;
                       std::copy(grid.face_nodes + np[f],
                                 grid.face_nodes + np[f+1],
                                 std::ostream_iterator<int>(os, " "));
                       os << '\n';
                       offset += f_num_pts;
                   }
                   cell_foffsets.push_back(offset);
               }
           }
           {
               pm["Name"] = "faceoffsets";
               Tag t("DataArray", pm, os);
               const int num_per_line = 10;
               for (int c = 0; c < num_cells; ++c) {
                   if (c % num_per_line == 0) {
                       Tag::indent(os);
                   }
                   os << cell_foffsets[c] << ' ';
                   if (c % num_per_line == num_per_line - 1
                       || c == num_cells - 1) {
                       os << '\n';
                   }
               }
           }
           {
               pm["type"] = "UInt8";
               pm["Name"] = "types";
               Tag t("DataArray", pm, os);
               const int num_per_line = 10;
               for (int c = 0; c < num_cells; ++c) {
                   if (c % num_per_line == 0) {
                       Tag::indent(os);
                   }
                   os << "42 ";
                   if (c % num_per_line == num_per_line - 1
                       || c == num_cells - 1) {
                       os << '\n';
                   }
               }
           }
       }
       {
           pm.clear();
           if (data.find("saturation") != data.end()) {
               pm["Scalars"] = "saturation";
           } else if (data.find("pressure") != data.end()) {
               pm["Scalars"] = "pressure";
           }
           Tag celldatatag("CellData", pm, os);
           pm.clear();
           pm["NumberOfComponents"] = "1";
           pm["format"] = "ascii";
           pm["type"] = "Float64";
           for (auto dit = data.begin(); dit != data.end(); ++dit) {
               pm["Name"] = dit->first;
               const std::vector<double>& field = *(dit->second);
               const int num_comps = field.size()/grid.number_of_cells;
               pm["NumberOfComponents"] = std::to_string(num_comps);
               Tag ptag("DataArray", pm, os);
               const int num_per_line = num_comps == 1 ? 5 : num_comps;
               for (int item = 0; item < num_cells*num_comps; ++item) {
                   if (item % num_per_line == 0) {
                       Tag::indent(os);
                   }
                   double value = field[item];
                   if (std::fabs(value) < std::numeric_limits<double>::min()) {
                       // Avoiding denormal numbers to work around
                       // bug in Paraview.
                       value = 0.0;
                   }
                   os << value << ' ';
                   if (item % num_per_line == num_per_line - 1
                       || item == num_cells - 1) {
                       os << '\n';
                   }
               }
           }
       }
    }

} // namespace Opm

