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
<simulation name="_SIM_NAME_" language="c" xmlns="http://damaris.gforge.inria.fr/damaris/model">
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
    <variable name="PRESSURE"    layout="zonal_layout_usmesh"     type="scalar"  visualizable="true"  mesh="us_mesh"   unit="_PRESSURE_UNIT_"   centering="zonal" select-file="GLOBAL_CELL_INDEX" store="_MYSTORE_OR_EMPTY_REGEX_"  script="_MAKE_AVAILABLE_IN_PYTHON_" />

    _MORE_VARIABLES_REGEX_
    <variable name="MPI_RANK"  layout="zonal_layout_usmesh_integer"   type="scalar"  visualizable="true" mesh="us_mesh" unit="rank"  centering="zonal"  store="_MYSTORE_OR_EMPTY_REGEX_" time-varying="false"  select-file="GLOBAL_CELL_INDEX"  script="_MAKE_AVAILABLE_IN_PYTHON_" comment="The MPI rank of each cell"/>
    
    <variable name="KRNSW_GO"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"     time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"   script="#" />
    <variable name="KRNSW_OW"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"     time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"   script="#" />
    <variable name="PCSWM_GO"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"     time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"   script="#" />
    <variable name="PCSWM_OW"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"     time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"   script="#" />
    <variable name="PPCW"      layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit="Bar"  centering="zonal"  time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"   script="#" />    
    <variable name="RS"        layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit="Bar"  centering="zonal"  time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"   script="#" comment="Dissolved Gas units Gas Oil Ratio" />
    <variable name="RV"        layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit="Bar"  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX" store="#"  script="#" />
    <variable name="SOMAX"     layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""     centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX" store="#"  script="#" />
    <variable name="1OVERBG"   layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />
    <variable name="1OVERBO"   layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />
    <variable name="1OVERBW"   layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />
    <variable name="GASKR"     layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />
    <variable name="GAS_DEN"   layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />
    <variable name="GAS_VISC"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />
    <variable name="OILKR"     layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />
    <variable name="OIL_DEN"   layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />
    <variable name="OIL_VISC"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />
    <variable name="WATKR"     layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />
    <variable name="WAT_DEN"   layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />
    <variable name="WAT_VISC"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />

    <variable name="2FBF"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"   time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="4FBF"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="DFBF"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="GCDI"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="GCDM"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="GIP"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"   time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />      
    <variable name="GIPG"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="GIPL"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"   time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="GIPR"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="HTOF"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="OIP"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"   time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />      
    <variable name="OIPG"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="OIPL"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="OIPR"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="RPV"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"   time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />      
    <variable name="S36F"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="SEAF"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="SGAS"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="SIP"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />      
    <variable name="SWAT"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="TFBF"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="WCD"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />      
    <variable name="WIP"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"   time-varying="true" select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />      
    <variable name="WIPG"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="WIPL"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     
    <variable name="WIPR"  layout="zonal_layout_usmesh"  type="scalar"  visualizable="true" mesh="#"  unit=""  centering="zonal"  time-varying="true"  select-file="GLOBAL_CELL_INDEX"  store="#"  script="#" />     

    <parameter name="n_coords_local"     type="int" value="1" />
    <layout    name="n_coords_layout"    type="double" dimensions="n_coords_local"   comment="For the individual x, y and z coordinates of the mesh vertices, these values are referenced in the topologies/topo/subelements/connectivity_pg data"  />
    <group name="coordset/coords/values"> 
        <variable name="x"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="_MAKE_AVAILABLE_IN_PYTHON_" time-varying="false" />
        <variable name="y"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="_MAKE_AVAILABLE_IN_PYTHON_" time-varying="false" />
        <variable name="z"    layout="n_coords_layout"  type="scalar"  visualizable="false"  unit="m"   script="_MAKE_AVAILABLE_IN_PYTHON_" time-varying="false" />
    </group>

    <parameter name="n_connectivity_ph"        type="int"  value="1" />
    <layout    name="n_connections_layout_ph"  type="int"  dimensions="n_connectivity_ph"   comment="Layout for connectivities "  />
    <parameter name="n_offsets_types_ph"       type="int"  value="1" />
    <layout    name="n_offsets_layout_ph"      type="int"  dimensions="n_offsets_types_ph + 1"  comment="Layout for the offsets_ph"  />
    <layout    name="n_types_layout_ph"        type="char" dimensions="n_offsets_types_ph"  comment="Layout for the types_ph "  />
    <group name="topologies/topo/elements">
        <variable name="connectivity" layout="n_connections_layout_ph"  type="scalar"  visualizable="false"    script="_MAKE_AVAILABLE_IN_PYTHON_" time-varying="false" />
        <variable name="offsets"      layout="n_offsets_layout_ph"    type="scalar"  visualizable="false"     script="_MAKE_AVAILABLE_IN_PYTHON_" time-varying="false" />
        <variable name="types"        layout="n_types_layout_ph"    type="scalar"  visualizable="false"    script="_MAKE_AVAILABLE_IN_PYTHON_" time-varying="false" />
    </group>

    <mesh name="us_mesh" type="unstructured" topology="3" time-varying="false" 
             comment="This Mesh definition is for connection with Paraview.
                      This definition references the variables that define an unstructured mesh specified above." >
        <coord                name="coordset/coords/values/x"  unit="m"    />
        <coord                name="coordset/coords/values/y"  unit="m"    />
        <coord                name="coordset/coords/values/z"  unit="m"    />
        <vertex_global_id     name="#"                         offset="0"  />
        <element_offsets      name="topologies/topo/elements/offsets"      />
        <section_types        name="topologies/topo/elements/types"        />
        <section_sizes        name="#"                                     />
        <section_connectivity name="topologies/topo/elements/connectivity" />
        <polyhedral_cell_faces_connectivity  name="#"                      /> 
        <polyhedral_cell_faces_offsets       name="#"                      />
        <polyhedral_n_faces_per_cell         name="#"                      />
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
    <_DISABLEPYTHONSTART_pyscript name="PythonScript" file="_PYTHON_SCRIPT_" language="python" frequency="1" scheduler-file="" nthreads="0" keep-workers="no" /_DISABLEPYTHONFIN_>
</scripts>

<_DISABLEPARAVIEWSTART_paraview update-frequency="1" write-vtk="0" write-vtk-binary="false" >
        <script>_PARAVIEW_PYTHON_SCRIPT_</script>
</paraview _DISABLEPARAVIEWFIN_>

<actions>
</actions>

<log FileName="_PATH_REGEX_/damaris_log/_SIM_NAME_" RotationSize="5" LogFormat="[%TimeStamp%]: %Message%"  Flush="_LOG_FLUSH_"  LogLevel="_LOG_LEVEL_" />

</simulation>)V0G0N";

    return init_damaris;
}

} // namespace Opm::DamarisOutput
