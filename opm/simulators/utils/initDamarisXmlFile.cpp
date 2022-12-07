/*
  Copyright 2022 KerData Research Team, Inria Rennes, Bretagne–Atlantique Research Center
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
*/
std::string initDamarisXmlFile()
{
    std::string init_damaris = R"V0G0N(<?xml version="1.0"?>
<simulation name="opm-flow" language="c" xmlns="http://damaris.gforge.inria.fr/damaris/model">
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
    <variable name="PRESSURE"    layout="zonal_layout_usmesh"     type="scalar"  visualizable="false"     unit="Pa"   centering="zonal"  store="_MYSTORE_OR_EMPTY_REGEX_" />
    _MORE_VARIABLES_REGEX_
</data>

<storage>
    <store name="MyStore" type="HDF5">
        <option key="FileMode">_File_Mode</option>
        <option key="XDMFMode">NoIteration</option>
        <option key="FilesPath">_PATH_REGEX_/</option>
    </store>
</storage>

<actions>
</actions>

<log FileName="_PATH_REGEX_/damaris_log/exa_dbg" RotationSize="5" LogFormat="[%TimeStamp%]: %Message%"  Flush="True"  LogLevel="debug" />

</simulation>)V0G0N";

    return init_damaris;
}

} // namespace Opm::DamarisOutput
