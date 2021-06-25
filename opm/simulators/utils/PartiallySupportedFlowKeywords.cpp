/*
  Copyright 2021 Equinor.

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

#include <opm/simulators/utils/PartiallySupportedFlowKeywords.hpp>

using namespace Opm::KeywordValidation;
namespace Opm::FlowKeywordValidation
{

template <>
const PartiallySupportedKeywords<std::string>&
partiallySupported()
{
   static const PartiallySupportedKeywords<std::string> partially_supported_keywords_strings = {
         {
            "BRINE",
            {
               {1,{false, allow_values<std::string> {}, "SALTS Multi-Component Brine model not supported"}}, // NO_NAME
            },
         },
         {
            "COMPORD",
            {
               {2,{false, allow_values<std::string> {"DEPTH", "INPUT", "TRACK"}, "COMPORD should be DEPTH/INPUT/TRACK"}}, // ORDER_TYPE
            },
         },
         {
            "EDITNNC",
            {
               {12,{false, allow_values<std::string> {}, "FACE1 not supported use 1* - will continue"}}, // FACE_FLOW12
               {13,{false, allow_values<std::string> {}, "FACE2 not supported use 1* - will continue"}}, // FACE_FLOW21
            },
         },
         {
            "EHYSTR",
            {
               {5,{false, allow_values<std::string> {"BOTH"}, "HYSTOPT only default value of BOTH is supported – will STOP"}}, // limiting_hyst_flag
               {6,{false, allow_values<std::string> {"RETR"}, "HYSTSCAN is not supported and is ignored"}}, // shape_cap_press_flag
               {7,{false, allow_values<std::string> {"DRAIN"}, "HYSTMOB is not supported and is ignored"}}, // init_fluid_mob_flag
               {8,{false, allow_values<std::string> {"OIL"}, "HYSTWET  is not supported and is ignored"}}, // wetting_phase_flag
               {9,{false, allow_values<std::string> {"NO"}, "NOTUSED and is ignored"}}, // baker_flag_oil
               {10,{false, allow_values<std::string> {"NO"}, "NOTUSED and is ignored"}}, // baker_flag_gas
               {11,{false, allow_values<std::string> {"NO"}, "NOTUSED and is ignored"}}, // baker_flag_water
            },
         },
         {
            "ENDSCALE",
            {
               {1,{false, allow_values<std::string> {"NODIR"}, "DIRECT only default value of NODIR supported"}}, // DIRECT
               {2,{false, allow_values<std::string> {"REVERS"}, "IRREVERS only default value of REVERS supported"}}, // IRREVERS
            },
         },
         {
            "EQLOPTS",
            {
               {1,{false, allow_values<std::string> {"THPRES"}, "MOBILE fluid critical saturation end point correction option not supported – value ignored"}}, // OPTION1
               {2,{false, allow_values<std::string> {"THPRES"}, "QUIESC quiescence option not supported – value ignored"}}, // OPTION2
               {3,{false, allow_values<std::string> {"THPRES"}, "THPRES invalid value for THPRES"}}, // OPTION3
               {4,{false, allow_values<std::string> {"THPRES"}, "IRREVER irreversible inter-region equilibration flow option not supported  - value ignored"}}, // OPTION4
            },
         },
         {
            "GRIDOPTS",
            {
               {1,{false, allow_values<std::string> {"NO"}, "TRANMULT only the NO option is supported – value ignored"}}, // TRANMULT
            },
         },
         {
            "MISCIBLE",
            {
               {3,{false, allow_values<std::string> {"NONE"}, "MISOPT only option NONE is supported – value ignored"}}, // TWOPOINT
            },
         },
         {
            "MULTIREG",
            {
               {4,{false, allow_values<std::string> {"F", "M", "O"}, "REGION_NAME must equal to F/M or O"}}, // NOT_DEFINED
            },
         },
         {
            "MULTREGP",
            {
               {3,{false, allow_values<std::string> {"F", "M", "O"}, "REGION_NAME must equal to F/M or O"}}, // REGION_TYPE
            },
         },
         {
            "MULTREGT",
            {
               {6,{false, allow_values<std::string> {"F", "M", "O"}, "REGION_NAME must equal to F/M or O"}}, // REGION_DEF
            },
         },
         {
            "NNC",
            {
               {12,{false, allow_values<std::string> {}, "FACE1 not supported use 1* - will continue"}}, // VE_FACE1
               {13,{false, allow_values<std::string> {}, "FACE2 not supported use 1* - will continue"}}, // VE_FACE2
            },
         },
         {
            "ROCKCOMP",
            {
               {1,{false, allow_values<std::string> {"REVERS"}, "ROCKOPT only the REVERS option is supported"}}, // HYSTERESIS
               {3,{false, allow_values<std::string> {"YES"}, "WATINOPT only equal to YES is supported"}}, // WATER_COMPACTION
               {4,{false, allow_values<std::string> {}, "PORTXROP transmissibility dependent on porosity model is not supported"}}, // PORTXROP
            },
         },
         {
            "SATOPTS",
            {
               {1,{false, allow_values<std::string> {"HYSTER"}, "DIRECT directional relative permeability assignment option not supported - value ignored"}}, // options
               {2,{false, allow_values<std::string> {"HYSTER"}, "IRREVERS reversible directional relative permeability assignment option not supported – value ignored"}}, // IRREVERS
               {3,{false, allow_values<std::string> {"HYSTER"}, "HYSTER hysteresis directional relative permeability assignment option not supported - value ignored"}}, // HYSTER
               {4,{false, allow_values<std::string> {"HYSTER"}, "SURFTENS capillary pressure surface tension pressure dependency option not supported – value ignored"}}, // SURFTENS
            },
         },
         {
            "SPECGRID",
            {
               {5,{true, allow_values<std::string> {"F"}, "TYPE only option F (Cartesian grids supported) supported – will STOP"}}, // COORD_TYPE
            },
         },
         {
            "TABDIMS",
            {
               {20,{false, allow_values<std::string> {}, "NOTUSED should be defaulted (1*) - ignored as not used"}}, // ITEM20_NOT_USED
               {25,{false, allow_values<std::string> {}, "RESVED should be defaulted (1*) - ignored as not used"}}, // RESERVED
            },
         },
         {
            "UDQDIMS",
            {
               {11,{false, allow_values<std::string> {"N"}, "RSEED option is not supported – value ignored"}}, // RESTART_NEW_SEED
            },
         },
   };

   return partially_supported_keywords_strings;
}

template <>
const KeywordValidation::PartiallySupportedKeywords<int>&
partiallySupported()
{
   static const KeywordValidation::PartiallySupportedKeywords<int>partially_supported_keywords_int = {
         {
            "EDITNNC",
            {
               {8,{false, allow_values<int> {0}, "ISATNUM1 only default value of 0 supported – will continue"}}, // SAT_TABLE12
               {9,{false, allow_values<int> {0}, "ISATNUM2 only default value of 0 supported – will continue"}}, // SAT_TABLE21
               {10,{false, allow_values<int> {0}, "IPRSNUM1 only default value of 0 supported – will continue"}}, // PRESS_TABLE12
               {11,{false, allow_values<int> {0}, "IPRSNUM2 only default value of 0 supported – will continue"}}, // PRESS_TABLE21
            },
         },
         {
            "EHYSTR",
            {
               {2,{false, allow_values<int> {0}, "HYSTMOD only default value of 0 (Carlson Hysteresis Model) supported"}}, // relative_perm_hyst
               {13,{false, allow_values<int> {0}, "NOTUSED and is ignored"}}, // FLAG_SOMETHING
            },
         },
         {
            "ENDSCALE",
            {
               {3,{false, allow_values<int> {1}, "NTENDP depth end-point scaling not supported – value ignored"}}, // NTENDP
               {4,{false, allow_values<int> {20}, "NNODES depth end-point scaling not supported – value ignored"}}, // NSENDP
               {5,{false, allow_values<int> {0}, "MODE depth temperature end-point scaling not supported – value ignored"}}, // COMP_MODE
            },
         },
         {
            "EQLDIMS",
            {
               {4,{false, allow_values<int> {1}, "NTTRVD tracer end-point depth scaling not supported – value ignored"}}, // NTTRVD
               {5,{false, allow_values<int> {20}, "NSTRVD tracer end-point depth scaling not supported – value ignored"}}, // NSTRVD
            },
         },
         {
            "GRIDFILE",
            {
               {1,{false, allow_values<int> {0}, "NGRID only default value of 0 supported – will continue"}}, // GRID
               {2,{false, allow_values<int> {1}, "NEGRID only default value of 1 supported – will continue"}}, // EGRID
            },
         },
         {
            "NETWORK",
            {
               {3,{false, allow_values<int> {20}, "NBCMAX option is not used and should be defaulted– value ignored"}}, // NBCMAX
            },
         },
         {
            "NNC",
            {
               {8,{false, allow_values<int> {0}, "ISATNUM1 only default value of 0 supported – will continue"}}, // SIM_DEPENDENT1
               {9,{false, allow_values<int> {0}, "ISATNUM2 only default value of 0 supported – will continue"}}, // SIM_DEPENDENT2
               {10,{false, allow_values<int> {0}, "IPRSNUM1 only default value of 0 supported – will continue"}}, // PRESSURE_TABLE1
               {11,{false, allow_values<int> {0}, "IPRSNUM2 only default value of 0 supported – will continue"}}, // PRESSURE_TABLE2
            },
         },
         {
            "NUMRES",
            {
               {1,{true, allow_values<int> {1}, "NUMRES only a value of one is supported – will STOP"}}, // NUM
            },
         },
         {
            "REGDIMS",
            {
               {5,{false, allow_values<int> {0}, "NUSREG compositional TRACK regions not supported - value ignored"}}, // "MAX_ETRACK
               {6,{false, allow_values<int> {1}, "NTCREG COAL regions not supported - value ignored"}}, // NTCREG
               {8,{false, allow_values<int> {0}, "NWKDREG should be equal to 0 - value ignored"}}, // MAX_OPERATE_DWORK
               {9,{false, allow_values<int> {0}, "NWKIREG should be equal to 0 - value ignored"}}, // MAX_OPERATE_IWORK
            },
         },
         {
            "SPECGRID",
            {
               {4,{false, allow_values<int> {1}, "NUMRES must be greater than or equal to 1"}}, // NUMRES
            },
         },
         {
            "TABDIMS",
            {
               {7,{false, allow_values<int> {20}, "NRVPVT should be defaulted (20) – ignored as not used"}}, // MAX_RV_NODES"
               {9,{false, allow_values<int> {1}, "NMEOSR must be greater than or equal to 1"}}, // NUM_EOS_RES
               {10,{false, allow_values<int> {1}, "NMEOSS should be equal to 1 - ignored as not used"}}, // NUM_EOS_SURFACE
               {12,{false, allow_values<int> {1}, "MXNTHR should be equal to 1 - ignored as not used"}}, // MAX_THERMAL_REGIONS
               {14,{false, allow_values<int> {0}, "MXNPMR should be equal to 0 - ignored as not used"}}, // MAX_PRESSURE_MAINTAINANCE_REGIONS
               {15,{false, allow_values<int> {0}, "NTABKT should be defaulted (0) – ignored as not used"}}, // MAX_KVALUE_TABLES
               {16,{false, allow_values<int> {0}, "NTALPHA should be defaulted (0) - ignored as not used"}}, // NTALPHA
               {17,{false, allow_values<int> {10}, "NASPKA should be defaulted (10) - ignored as not used"}}, // ASPHALTENE_ASPKDAM_MAX_ROWS
               {18,{false, allow_values<int> {10}, "MXRAWG should be defaulted (10) - ignored as not used"}}, // ASPHALTENE_ASPREWG_MAX_ROWS
               {19,{false, allow_values<int> {10}, "MXRASO should be defaulted (10) - ignored as not used"}}, // ASPHALTENE_ASPVISO_MAX_ROWS
               {21,{false, allow_values<int> {5}, "MCASPP should be defaulted (5) - ignored as not used"}}, // ASPHALTENE_ASPPW2D_MAX_COLUMNS
               {22,{false, allow_values<int> {5}, "MRASPP should be defaulted (5) - ignored as not used"}}, // ASPHALTENE_ASPPW2D_MAX_ROWS
               {23,{false, allow_values<int> {5}, "MXRATF should be defaulted (5) - ignored as not used"}}, // ASPHALTENE_ASPWETF_MAX_ROWS
               {24,{false, allow_values<int> {0}, "MXNKVT should be defaulted (0) - ignored as not used"}}, // NUM_KVALUE_TABLES
            },
         },
         {
            "UDADIMS",
            {
               {2,{false, allow_values<int> {0}, "IGNORED should be defaulted (0) – ignored as not used"}}, // IGNORED
            },
         },
         {
            "UDQPARAM",
            {
               {1,{false, allow_values<int> {1}, "RSEED option not supported – value ignored"}}, // RANDOM_SEED
            },
         },
         {
            "WELLDIMS",
            {
               {5,{false, allow_values<int> {5}, "MXSTAGE option not supported – value ignored"}}, // MAX_STAGES
               {6,{false, allow_values<int> {10}, "MXSTRMS option not supported – value ignored"}}, // MAX_STREAMS
               {7,{false, allow_values<int> {5}, "MXMIS option not supported – value ignored"}}, // MAX_MIXTURES
               {8,{false, allow_values<int> {4}, "MXSEPS option not supported – value ignored"}}, // MAX_SEPARATORS
               {9,{false, allow_values<int> {3}, "MXCOMPS option not supported – value ignored"}}, // MAX_MIXTURE_ITEMS
               {10,{false, allow_values<int> {0}, "MXDOCOMP option not supported – value ignored"}}, // MAX_COMPLETION_X
               {11,{false, allow_values<int> {1}, "MXWSLIST option not supported – value ignored"}}, // MAX_WELLIST_PR_WELL
               {12,{false, allow_values<int> {1}, "MXWLISTS option not supported – value ignored"}}, // MAX_DYNAMIC_WELLIST
               {13,{false, allow_values<int> {10}, "MXWSECD option not supported – value ignored"}}, // MAX_SECONDARY_WELLS
               {14,{false, allow_values<int> {201}, "MXNGPP option not supported – value ignored"}}, // NO_JASON_ENTRY
            },
         },
   };

   return partially_supported_keywords_int;
}

template <>
const KeywordValidation::PartiallySupportedKeywords<double>&
partiallySupported()
{
   static const KeywordValidation::PartiallySupportedKeywords<double> partially_supported_keywords_double = {
         {
            "AQUCON",
            {
               {12,{false, allow_values<double> {1}, "VEOPT1 Vertical Equilibrium option one– not used"}}, // VEFRAC
               {13,{false, allow_values<double> {1}, "VEOPI2 Vertical Equilibrium option two– not used"}}, // VEFRACP
            },
         },
         {
            "AQUCT",
            {
               {12,{false, allow_values<double> {0}, "SALTCON option is not used and should be defaulted – value ignored"}}, // INI_SALT
               {13,{false, allow_values<double> {}, "TEMP option is not used and should be defaulted – value ignored"}}, // TEMP_AQUIFER
            },
         },
         {
            "EDITNNC",
            {
               {14,{false, allow_values<double> {0}, "DIFFNNC not supported – will continue"}}, // DIFFM
            },
         },
         {
            "EHYSTR",
            {
               {1,{false, allow_values<double> {0.1}, "HYSTRCP option not supported – value ignored"}}, // curvature_caplillary_pressure_hyst
               {3,{false, allow_values<double> {1.00}, "HYSTREL is not supported and is ignored"}}, // curvature_param_killough_wetting
               {4,{false, allow_values<double> {0.1}, "HYSTSGR is not supported and is ignored"}}, // mod_param_trapped
               {12,{false, allow_values<double> {0}, "NOTUSED and is ignored"}}, // threshold_saturation
            },
         },
         {
            "MULTFLT",
            {
               {3,{false, allow_values<double> {}, "FLT-DIF the diffusivity multiplier option is not supported – will continue"}}, // NOT_DEFINED
            },
         },
         {
            "NNC",
            {
               {14,{false, allow_values<double> {0}, "DIFFNNC not supported – will continue"}}, // DIFFUSIVITY
               {15,{false, allow_values<double> {0}, "DISPNNC not supported – will continue"}}, // SIM_DEPENDENT3
               {16,{false, allow_values<double> {0}, "AREANNC not supported – will continue"}}, // VDFLOW_AREA
               {17,{false, allow_values<double> {0}, "PERMNNC only default value of 0 supported – will continue"}}, // VDFLOW_PERM
            },
         },
         {
            "ROCKCOMP",
            {
               {4,{false, allow_values<double> {0}, "CARKZEXP transmissibility dependent on porosity model is not supported"}}, // CARKZEXP
            },
         },
   };

   return partially_supported_keywords_double;
}

} // namespace Opm::FlowKeywordValidation
