/*
  Copyright 2022 2023 Inria, Bretagne–Atlantique Research Center
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.

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

#include <string>

namespace Opm::DamarisOutput
{


/*
    Below is the XML file for Damaris that is supported by Damaris.

    The entries in the map below will be filled by corresponding Damaris
    Keywords.
    
    N.B. Ensure all text items that are to be replaced are quoted with double quotes
    e.g. unit="_REPLACE_UNIT_"   
    *not* unit=_REPLACE_UNIT_
*/
std::string initDamarisXmlFile()
{
    std::string init_damaris = R"V0G0N(<?xml version="1.0"?>
<simulation name="opm-flow-sim" language="c" xmlns="http://damaris.gforge.inria.fr/damaris/model">
<architecture>
    <domains count="1"/>
    <dedicated cores="_DC_REGEX_" nodes="_DN_REGEX_"/>
    <buffer name="buffer" size="_SHMEM_BUFFER_BYTES_REGEX_" />
    <placement />
    <queue  name="queue" size="300" />
</architecture>

<data>
    <parameter name="n_elements_total"     type="int" value="1" />
    <parameter name="n_elements_local"     type="int" value="1" />
    <parameter name="n"     type="int" value="1" />

    <layout   name="zonal_layout_usmesh_integer"             type="int" dimensions="n_elements_local"   global="n_elements_total"   comment="For the field data e.g. Pressure"  />
    <variable name="GLOBAL_CELL_INDEX"    layout="zonal_layout_usmesh_integer"     type="scalar"  visualizable="false"  time-varying="false"  centering="zonal" />
    <layout   name="zonal_layout_usmesh"             type="double" dimensions="n_elements_local"   global="n_elements_total"   comment="For the field data e.g. Pressure"  />
    <variable name="PRESSURE"    layout="zonal_layout_usmesh"     type="scalar"  visualizable="false"     unit="bar"   centering="zonal"  store="MyStore"  script="PythonScript" />
    _MORE_VARIABLES_REGEX_
    
    
    <parameter name="n_coords_local"     type="int" value="1" />
    <layout    name="n_coords_layout"    type="double" dimensions="n_coords_local"   comment="For the individual x, y and z coordinates of the mesh vertices, these values are referenced in the topologies/topo/subelements/connectivity_pg data"  />
    <group name="coordset/coords/values"> 
        <variable name="x"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonScript" time-varying="false" />
        <variable name="y"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonScript" time-varying="false" />
        <variable name="z"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="PythonScript" time-varying="false" />
    </group>
    
    
    <parameter name="n_connectivity_ph"        type="int"  value="1" />
    <layout    name="n_connections_layout_ph"  type="int"  dimensions="n_connectivity_ph"   comment="Layout for connectivities "  />
    <parameter name="n_offsets_types_ph"       type="int"  value="1" />
    <layout    name="n_offsets_layout_ph"      type="int"  dimensions="n_offsets_types_ph"  comment="Layout for the offsets_ph"  />
    <layout    name="n_types_layout_ph"        type="char" dimensions="n_offsets_types_ph"  comment="Layout for the types_ph "  />
    <group name="topologies/topo/elements">
        <variable name="connectivity" layout="n_connections_layout_ph"  type="scalar"  visualizable="false"    script="PythonScript" time-varying="false" />
        <variable name="offsets"      layout="n_offsets_layout_ph"    type="scalar"  visualizable="false"     script="PythonScript" time-varying="false" />
        <variable name="types"        layout="n_types_layout_ph"    type="scalar"  visualizable="false"    script="PythonScript" time-varying="false" />
    </group>
    
    
    <mesh name="unstructured_mesh" type="unstructured" topology="3" time-varying="false" 
               comment="This Mesh definition is for connection with Paraview.
                        This definition only references the actual variables that define an unstructured mesh available above.
                        Current issues: 
                          1. x,y,z vertex coordinates are separate - need to test if 3 individual arrays works;
                          2. GID is not defined for ASCENT data above; 
                          3. Paraview is expecting sizes and not offsets" >
          <coord                name="coordset/coords/values/x"  unit="m"    />
          <coord                name="coordset/coords/values/y"  unit="m"    />
          <coord                name="coordset/coords/values/z"  unit="m"    />
          <vertex_global_id     name=""                          offset="-1" />
          <section_types        name="topologies/topo/elements/types"        />
          <section_sizes        name="topologies/topo/elements/offsets"      />
          <section_connectivity name="topologies/topo/elements/connectivity" />
    </mesh>
    
</data>

<storage>
    <store name="MyStore" type="HDF5">
        <option key="FileMode">_File_Mode</option>
        <option key="XDMFMode">NoIteration</option>
        <option key="FilesPath">_PATH_REGEX_/</option>
    </store>
</storage>

<scripts>
    <pyscript name="PythonScript" file="polygonal_mesh_conduit.py" language="python" frequency="1" scheduler-file="" nthreads="0" keep-workers="no" />
</scripts>

<paraview update-frequency="1" >
        <script></script>
</paraview>

<actions>
</actions>

<log FileName="_PATH_REGEX_/damaris_log/opm-flow" RotationSize="5" LogFormat="[%TimeStamp%]: %Message%"  Flush="True"  LogLevel="info" />

</simulation>)V0G0N";

    return init_damaris;
}

} // namespace Opm::DamarisOutput
