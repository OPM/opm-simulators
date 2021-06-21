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
 
  Generated : keyw-code.py 
  Date      : 2021-06-21 14:08:30
*/ 
 
#include <opm/simulators/utils/PartiallySupportedFlowKeywords.hpp> 
 
using namespace Opm::KeywordValidation;
namespace Opm::FlowKeywordValidation
{ 
 
template<>  
const PartiallySupportedKeywords<std::string>&  
partiallySupported()   
{ 
   static const PartiallySupportedKeywords<std::string> partially_supported_keywords_strings = {
         { 
            "AQUANCON", 
            { 
               {8,{false, allow_values<std::string> {"X+", "Y+", "Z+", "X-", "Y-", "Z-", "I+", "J+", "K+", "I-", "J-", "K-"}, "AQUFACE has an invalid aquifer face"}}, // FACE
               {11,{false, allow_values<std::string> {"YES", "NO"}, "AQUOPT must be either YES or NO"}}, // CONNECT_ADJOINING_ACTIVE_CELL
            },
         },
         { 
            "AQUCON", 
            { 
               {8,{false, allow_values<std::string> {"X+", "Y+", "Z+", "X-", "Y-", "Z-", "I+", "J+", "K+", "I-", "J-", "K-"}, "AQUFACE has an invalid aquifer face"}}, // CONNECT_FACE
               {11,{false, allow_values<std::string> {"YES", "NO"}, "AQUOPT2 must be either YES or NO"}}, // ALLOW_INTERNAL_CELLS
            },
         },
         { 
            "BRINE", 
            { 
               {1,{false, allow_values<std::string> {}, "SALTS Multi-Component Brine model not supported"}}, // NO_NAME
            },
         },
         { 
            "COMPORD", 
            { 
               {1,{false, allow_values<std::string> {}, "WELNAME must be a character string"}}, // WELL
               {2,{false, allow_values<std::string> {"DEPTH", "INPUT", "TRACK"}, "COMPORD should be DEPTH/INPUT/TRACK"}}, // ORDER_TYPE
            },
         },
         { 
            "EHYSTR", 
            { 
               {5,{false, allow_values<std::string> {"BOTH"}, "HYSTOPT only default value of BOTH is supported"}}, // limiting_hyst_flag
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
               {1,{false, allow_values<std::string> {}, "MOBILE fluid critical saturation end point correction option not supported – value ignored"}}, // OPTION1
               {2,{false, allow_values<std::string> {}, "QUIESC quiescence option not supported – value ignored"}}, // OPTION2
               {3,{false, allow_values<std::string> {"THPRESS"}, "THPRES invalid value for THPRES"}}, // OPTION3
               {4,{false, allow_values<std::string> {}, "IRREVER irreversible inter-region equilibration flow option not supported  - value ignored"}}, // OPTION4
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
            "PINCH", 
            { 
               {2,{false, allow_values<std::string> {"NOGAP"}, "PINCHOPT only NOGAP option supported"}}, // CONTROL_OPTION
               {4,{false, allow_values<std::string> {"TOPBOT", "ALL"}, "PINCHCAL should equal TOPBOT or ALL"}}, // PINCHOUT_OPTIONL
               {5,{false, allow_values<std::string> {"TOP", "ALL"}, "PINCHMUL should equal TOP or ALL"}}, // MULTZ_OPTION
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
               {1,{false, allow_values<std::string> {}, "DIRECT directional relative permeability assignment option not supported - value ignored"}}, // options
               {2,{false, allow_values<std::string> {}, "IRREVERS reversible directional relative permeability assignment option not supported – value ignored"}}, // IRREVERS
               {3,{false, allow_values<std::string> {}, "HYSTER hysteresis directional relative permeability assignment option not supported - value ignored"}}, // HYSTER
               {4,{false, allow_values<std::string> {}, "SURFTENS capillary pressure surface tension pressure dependency option not supported – value ignored"}}, // SURFTENS
            },
         },
         { 
            "UDQDIMS", 
            { 
               {11,{false, allow_values<std::string> {"Y", "N"}, "RSEED should be equal to Y or N"}}, // RESTART_NEW_SEED
            },
         },
   }; 
 
   return partially_supported_keywords_strings; 
} 
 
template<>  
const KeywordValidation::PartiallySupportedKeywords<int>&  
partiallySupported()   
{ 
   static const KeywordValidation::PartiallySupportedKeywords<int> partially_supported_keywords_int = { 
         { 
            "ACTDIMS", 
            { 
               {1,{false, [](int x) { return x >= 0; }, "MXACTNS must be greater than or equal to 0"}}, // MAX_ACTION
               {2,{false, [](int x) { return x >= 0; }, "MXLINES must be greater than or equal to 0"}}, // MAX_ACTION_LINES
               {3,{false, [](int x) { return x >= 0; }, "MXCHARS must be greater than or equal to 0"}}, // MAX_ACTION_LINE_CHARACTERS
               {4,{false, [](int x) { return x >= 0; }, "MXSTATMS must be greater than or equal to 0"}}, // MAX_ACTION_COND
            },
         },
         { 
            "AQUANCON", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "AQUID must be greater than zero"}}, // AQUIFER_ID
               {2,{false, [](int x) { return x >= 1; }, "I1 must be greater than zero"}}, // I1
               {3,{false, [](int x) { return x >= 1; }, "I2 must be greater than zero"}}, // I2
               {4,{false, [](int x) { return x >= 1; }, "J1 must be greater than zero"}}, // J1
               {5,{false, [](int x) { return x >= 1; }, "J2 must be greater than zero"}}, // J1
               {6,{false, [](int x) { return x >= 1; }, "K1 must be greater than zero"}}, // K1
               {7,{false, [](int x) { return x >= 1; }, "K2 must be greater than zero"}}, // K2
            },
         },
         { 
            "AQUCON", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "AQUID must be greater than zero"}}, // ID
               {2,{false, [](int x) { return x >= 1; }, "I1 must be greater than zero"}}, // I1
               {3,{false, [](int x) { return x >= 1; }, "I2 must be greater than zero"}}, // I2
               {4,{false, [](int x) { return x >= 1; }, "J1 must be greater than zero"}}, // J1
               {5,{false, [](int x) { return x >= 1; }, "J2 must be greater than zero"}}, // J1
               {6,{false, [](int x) { return x >= 1; }, "K1 must be greater than zero"}}, // K1
               {7,{false, [](int x) { return x >= 1; }, "K2 must be greater than zero"}}, // K2
               {10,{false, [](int x) { return x >= 0 && x <= 1; }, "AQUOPT1 must be equal to one or zero"}}, // TRANS_OPTION
            },
         },
         { 
            "AQUCT", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "AQUID must be greater than zero"}}, // AQUIFER_ID
               {10,{false, [](int x) { return x >= 1; }, "PVTNUM must be greater than one"}}, // TABLE_NUM_WATER_PRESS
               {11,{false, [](int x) { return x > 0; }, "AQUTAB must be greater than zero and less than AQUDIMS(NIFTBL)"}}, // TABLE_NUM_INFLUENCE_FN
            },
         },
         { 
            "AQUDIMS", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "MXNAQN must be greater than or equal to 1"}}, // MXNAQN
               {2,{false, [](int x) { return x >= 1; }, "MXNAQC must be greater than or equal to 1"}}, // MXNAQC
               {3,{false, [](int x) { return x >= 1; }, "NIFTB must be greater than or equal to 1"}}, // NIFTBL
               {4,{false, [](int x) { return x >= 36; }, "NRIFTB must be greater than or equal to 36"}}, // NRIFTB
               {5,{false, [](int x) { return x >= 1; }, "NANAQ must be greater than or equal to 1"}}, // NANAQU
               {6,{false, [](int x) { return x >= 1; }, "NCAMAX must be greater than or equal to 1"}}, // NCAMAX
               {7,{false, [](int x) { return x >= 0; }, "MXNALI must be greater than or equal to 0"}}, // MXNALI
               {8,{false, [](int x) { return x >= 0; }, "MXAAQL must be greater than or equal to 0"}}, // MXAAQL
            },
         },
         { 
            "DIMENS", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "NX must be greater than or equal to 1"}}, // NX
               {2,{false, [](int x) { return x >= 1; }, "NY must be greater than or equal to 1"}}, // NY
               {3,{false, [](int x) { return x >= 1; }, "NZ must be greater than or equal to 1"}}, // NZ
            },
         },
         { 
            "EHYSTR", 
            { 
               {2,{false, [](int x) { return x ==0; }, "HYSTMOD only default value of 0 (Carlson Hysteresis Model) supported"}}, // relative_perm_hyst
               {13,{false, [](int x) { return x ==0; }, "NOTUSED and is ignored"}}, // FLAG_SOMETHING
            },
         },
         { 
            "ENDSCALE", 
            { 
               {3,{false, [](int x) { return x == 1; }, "NTENDP depth end-point scaling not supported – value ignored"}}, // NTENDP
               {4,{false, [](int x) { return x == 20; }, "NNODES depth end-point scaling not supported – value ignored"}}, // NSENDP
               {5,{false, [](int x) { return x == 0; }, "MODE depth temperature end-point scaling not supported"}}, // COMP_MODE
            },
         },
         { 
            "EQLDIMS", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "NTEQUL must be greater then or equal to one"}}, // NTEQUL
               {2,{false, [](int x) { return x >= 2; }, "NPRSVD must be greater then or equal to two"}}, // "DEPTH_NODES_P
               {3,{false, [](int x) { return x >= 2; }, "NDRXVD must be greater then or equal to two"}}, // DEPTH_NODES_TAB
               {4,{false, [](int x) { return x == 1; }, "NTTRVD tracer end-point depth scaling not supported – value ignored"}}, // NTTRVD
               {5,{false, [](int x) { return x == 20; }, "NSTRVD tracer end-point depth scaling not supported – value ignored"}}, // NSTRVD
            },
         },
         { 
            "GRIDOPTS", 
            { 
               {2,{false, [](int x) { return x >= 0; }, "NRMULT must be greater than or equal to 0"}}, // NRMULT
               {3,{false, [](int x) { return x >= 0; }, "NRPINC must be greater than or equal to 0"}}, // NRPINC
            },
         },
         { 
            "MISCIBLE", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "NTMISC must be greater than or equal to 1"}}, // NTMISC
               {2,{false, [](int x) { return x >= 1; }, "NSMISC must be greater than or equal to 1"}}, // NSMISC
            },
         },
         { 
            "NETWORK", 
            { 
               {1,{false, [](int x) { return x >= 0; }, "NODMAX must be greater than zero"}}, // NODMAX
               {2,{false, [](int x) { return x >= 0; }, "NBRMAX must be greater than zero"}}, // NBRMAX
               {3,{false, allow_values<int>{{}}, "NBCMAX option is not used – value ignored"}}, // NBCMAX
            },
         },
         { 
            "NUMRES", 
            { 
               {1,{true, [](int x) { return x == 1; }, "NUMRES only a value of one is supported – will STOP"}}, // NUM
            },
         },
         { 
            "NUPCOL", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "NUPCOL must be greater than or equal to 1"}}, // NUM_ITER
            },
         },
         { 
            "PIMTDIMS", 
            { 
               {1,{false, [](int x) { return x >= 0; }, "NTPIMT must be greater than or equal to 0"}}, // NTPIMT
               {2,{false, [](int x) { return x >= 0; }, "NRPIMT must be greater than or equal to 0"}}, // NRPIMT
            },
         },
         { 
            "REGDIMS", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "NTFIP must be greater than or equal to 1"}}, // NTFIP
               {2,{false, [](int x) { return x >= 1; }, "NMFIPR must be greater than or equal to 1"}}, // NMFIPR
               {3,{false, [](int x) { return x >= 0; }, "NRFREG should be greater than 0"}}, // NRFREG
               {4,{false, [](int x) { return x >= 0; }, "MXNFLX must be greater than 0"}}, // MXNFLX
               {5,{false, [](int x) { return x == 0; }, "NUSREG compositional TRACK regions not supported - value ignored"}}, // "MAX_ETRACK
               {6,{false, [](int x) { return x == 1; }, "NTCREG COAL regions not supported - value ignored"}}, // NTCREG
               {7,{false, [](int x) { return x >= 0; }, "NOPREG must be greater than or equal to 0"}}, // MAX_OPERNUM
               {8,{false, [](int x) { return x == 0; }, "NWKDREG should be equal to 0 - value ignored"}}, // MAX_OPERATE_DWORK
               {9,{false, [](int x) { return x == 0; }, "NWKIREG should be equal to 0 - value ignored"}}, // MAX_OPERATE_IWORK
               {10,{false, [](int x) { return x >= 1; }, "NPLMIX must be greater than or equal to 1"}}, // NPLMIX
            },
         },
         { 
            "ROCKCOMP", 
            { 
               {2,{false, [](int x) { return x >= 1; }, "NTROCC must be greater than or equal to 1"}}, // NTROCC
            },
         },
         { 
            "SATOPTS", 
            { 
               {5,{false, [](int x) { return x > 0; }, "MODE depth temperature end-point scaling not supported"}}, // MODE
            },
         },
         { 
            "SMRYDIMS", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "Must be greater than or equal to 1"}}, // DIMS
            },
         },
         { 
            "TABDIMS", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "NTSFUN must be greater than or equal to 1"}}, // NTSFUN
               {2,{false, [](int x) { return x >= 1; }, "NTPVT must be greater than or equal to 1"}}, // NTPVT
               {3,{false, [](int x) { return x >= 1; }, "NSSFUN must be greater than or equal to 1"}}, // NSSFUN
               {4,{false, [](int x) { return x >= 1; }, "NPPVT must be greater than or equal to 1"}}, // NPPVT
               {5,{false, [](int x) { return x >= 1; }, "NTFIP must be greater than or equal to 1"}}, // NTFIP
               {6,{false, [](int x) { return x >= 1; }, "NRPVT must be greater than or equal to 1"}}, // NRPVT
               {7,{false, [](int x) { return x == 1; }, "NRVPVT should be equal to 1 – ignored as not used"}}, // MAX_RV_NODES"
               {8,{false, [](int x) { return x >= 1; }, "NTENDP must be greater than or equal to 1"}}, // NTENDP
               {9,{false, [](int x) { return x == 1; }, "NMEOSR must be greater than or equal to 1"}}, // NUM_EOS_RES
               {10,{false, [](int x) { return x == 1; }, "NMEOSS should be equal to 1 - ignored as not used"}}, // NUM_EOS_SURFACE
               {11,{false, [](int x) { return x >= 1; }, "MXNFLN must be greater than or equal to 1"}}, // MAX_FLUX_REGIONS
               {12,{false, [](int x) { return x == 1; }, "MXNTHR should be equal to 1 - ignored as not used"}}, // MAX_THERMAL_REGIONS
               {13,{false, [](int x) { return x >= 1; }, "NTROCC must be greater than or equal to 1"}}, // NTROCC
               {14,{false, [](int x) { return x == 0; }, "MXNPMR should be equal to 0 - ignored as not used"}}, // MAX_PRESSURE_MAINTAINANCE_REGIONS
               {15,{false, [](int x) { return x == 0; }, "NTABKT should be equal to 0 - ignored as not used"}}, // MAX_KVALUE_TABLES
               {16,{false, [](int x) { return x == 0; }, "NTALPHA should be equal to 0 - ignored as not used"}}, // NTALPHA
               {17,{false, [](int x) { return x == 0; }, "NASPKA should be equal to 0 - ignored as not used"}}, // ASPHALTENE_ASPKDAM_MAX_ROWS
               {18,{false, [](int x) { return x == 0; }, "MXRAWG should be equal to 0 - ignored as not used"}}, // ASPHALTENE_ASPREWG_MAX_ROWS
               {19,{false, [](int x) { return x == 0; }, "MXRASO should be equal to 0 - ignored as not used"}}, // ASPHALTENE_ASPVISO_MAX_ROWS
               {20,{false, [](int x) { return x == 0; }, "NOTUSED should be equal to 0 - ignored as not used"}}, // ITEM20_NOT_USED
               {21,{false, [](int x) { return x == 0; }, "MCASPP should be equal to 0 - ignored as not used"}}, // ASPHALTENE_ASPPW2D_MAX_COLUMNS
               {22,{false, [](int x) { return x == 0; }, "MRASPP should be equal to 0 - ignored as not used"}}, // ASPHALTENE_ASPPW2D_MAX_ROWS
               {23,{false, [](int x) { return x == 0; }, "MXRATF should be equal to 0 - ignored as not used"}}, // ASPHALTENE_ASPWETF_MAX_ROWS
               {24,{false, [](int x) { return x == 0; }, "MXNKVT should be equal to 0 - ignored as not used"}}, // NUM_KVALUE_TABLES
               {25,{false, [](int x) { return x == 0; }, "RESVED should be equal to 0 - ignored as not used"}}, // RESERVED
            },
         },
         { 
            "UDADIMS", 
            { 
               {1,{false, [](int x) { return x >= 0; }, "NMUDA must be greater than or equal to 0"}}, // NUM_UDQ_REPLACE
               {2,{false, allow_values<int>{{}}, "IGNORED should be equal to 0 – ignored as not used"}}, // IGNORED
               {3,{false, [](int x) { return x >= 0; }, "NMUDA must be greater than or equal to 0"}}, // TOTAL_UDQ_UNIQUE
            },
         },
         { 
            "UDQDIMS", 
            { 
               {1,{false, [](int x) { return x >= 1; }, "MXFUN must be greater than or equal to 1"}}, // MAX_FUNCTIONS
               {2,{false, [](int x) { return x >= 1; }, "MXITEMS must be greater than or equal to 1"}}, // MAX_ITEMS
               {3,{false, [](int x) { return x >= 0; }, "MXUDC must be greater than or equal to 0"}}, // MAX_CONNECTIONS
               {4,{false, [](int x) { return x >= 0; }, "MXUDF must be greater than or equal to 0"}}, // MAX_FIELDS
               {5,{false, [](int x) { return x >= 0; }, "MXUDG must be greater than or equal to 0"}}, // MAX_GROUP
               {6,{false, [](int x) { return x >= 0; }, "MXUDR must be greater than or equal to 0"}}, // MAX_REGION
               {7,{false, [](int x) { return x >= 0; }, "MXUDS must be greater than or equal to 0"}}, // MAX_SEGMENT
               {8,{false, [](int x) { return x >= 0; }, "MXUDW must be greater than or equal to 0"}}, // MAX_WELL
               {9,{false, [](int x) { return x >= 0; }, "MXUDA must be greater than or equal to 0"}}, // MAX_AQUIFER
               {10,{false, [](int x) { return x >= 0; }, "MXUDB must be greater than or equal to 0"}}, // MAX_BLOCK
            },
         },
         { 
            "UDQPARAM", 
            { 
               {1,{false, [](int x) { return x > 0; }, "RSEED must be greater than 0"}}, // RANDOM_SEED
            },
         },
   }; 
 
   return partially_supported_keywords_int; 
} 
 
template<>  
const KeywordValidation::PartiallySupportedKeywords<double>&  
partiallySupported()   
{ 
   static const KeywordValidation::PartiallySupportedKeywords<double> partially_supported_keywords_double = { 
         { 
            "AQUANCON", 
            { 
               {9,{false, [](double x) { return x >= 0.0; }, "AQUFLUX must be greater than or equal to zero"}}, // INFLUX_COEFF
               {10,{false, [](double x) { return x >= 0.0; }, "AQUCOEF must be greater than or equal to zero"}}, // CONNECT_ADJOINING_ACTIVE_CELL
            },
         },
         { 
            "AQUCON", 
            { 
               {9,{false, [](double x) { return x >= 0.0; }, "AQUMULT must be greater than or equal to zero"}}, // TRANS_MULT
            },
         },
         { 
            "AQUCT", 
            { 
               {2,{false, [](double x) { return x >= 0.0; }, "DATUM must be greater than zero"}}, // DAT_DEPTH
               {3,{false, [](double x) { return x > 0.0; }, "PRESS must be greater than or equal to zero"}}, // "P_INI
               {4,{false, [](double x) { return x > 0.0; }, "PERM must be greater than zero"}}, // PERM_AQ
               {5,{false, [](double x) { return x > 0.0; }, "PORO must be greater then zero"}}, // PORO_AQ
               {6,{false, [](double x) { return x > 0.0 && x <= 1.0; }, "RCOMP must greater than zero and less than one"}}, // C_T
               {7,{false, [](double x) { return x > 0.0; }, "RE must be greater than zero"}}, // RAD
               {8,{false, [](double x) { return x > 0.0; }, "DZ must be greater than zero"}}, // THICKNESS_AQ
               {9,{false, [](double x) { return x > 0.0 && x <= 360.0; }, "ANGLE must be greater than zero"}}, // INFLUENCE_ANGLE
               {12,{false, allow_values<double>{0.0}, "SALTCON option is not used – value ignored"}}, // INI_SALT
               {13,{false, allow_values<double>{{}}, "TEMP option is not used – value ignored"}}, // TEMP_AQUIFER
            },
         },
         { 
            "EHYSTR", 
            { 
               {1,{false, [](double x) { return x ==0.1; }, "HYSTRCP option not supported – value ignored"}}, // curvature_caplillary_pressure_hyst
               {3,{false, [](double x) { return x ==1.00; }, "HYSTREL is not supported and is ignored"}}, // curvature_param_killough_wetting
               {4,{false, [](double x) { return x ==0.1; }, "HYSTSGR is not supported and is ignored"}}, // mod_param_trapped
               {12,{false, [](double x) { return x ==0.0; }, "NOTUSED and is ignored"}}, // threshold_saturation
            },
         },
         { 
            "PINCH", 
            { 
               {1,{false, [](double x) { return x >= 0.0; }, "PINCHTHK must be greater than or equal to zero"}}, // THRESHOLD_THICKNESS
               {3,{false, [](double x) { return x <= 1E+20; }, "PINCHGAP value should be less than equal to 1E+20"}}, // MAX_EMPTY_GAP
            },
         },
         { 
            "UDQPARAM", 
            { 
               {2,{false, [](double x) { return x >= 1.0 && x <= 1.0E20; }, "RANGE must be greater than or equal to 1 and less than or equal to 1.0E20"}}, // RANGE
               {3,{false, [](double x) { return x >= -1.0E20 && x <= 1.0E20; }, "DEFAULT must be in the range -1.0E20 to 1.0E20 - normally set to 0.0"}}, // UNDEFINED_VALUE
               {4,{false, [](double x) { return x >= 0.0  && x <= 1.0; }, "TOLUDQ must be in the range 0.0 to 1.0"}}, // CMP_EPSILON
            },
         },
   }; 
 
   return partially_supported_keywords_double; 
} 
 
} // namespace Opm::FlowKeywordValidation
