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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

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
            "BRANPROP",
            {
               {5,{true, allow_values<std::string> {"NONE"}, "BRANPROP(ALQ-DEN): option is not supported. Default is supported."}}, // ALQ_SURFACE_DENSITY
            },
         },
         {
            "EDITNNC",
            {
               {12,{true, allow_values<std::string> {"NONE"}, "EDITNNC(FACE1): not supported use default."}}, // FACE_FLOW12
               {13,{true, allow_values<std::string> {"NONE"}, "EDITNNC(FACE2): not supported use default."}}, // FACE_FLOW21
            },
         },
         {
            "EHYSTR",
            {
               {5,{true, allow_values<std::string> {"KR", "PC", "BOTH"}, "EHYSTR(HYSTOPT): relative permeability hysteresis option should equal BOTH, KR, or PC"}}, // limiting_hyst_flag
               {6,{true, allow_values<std::string> {"RETR"}, "EHYSTR(HYSTSCAN): only RETR supported"}}, // shape_cap_press_flag
               {7,{false, allow_values<std::string> {"DRAIN"}, "EHYSTR(HYSTMOB): mobility option not supported and is ignored"}}, // init_fluid_mob_flag
               {8,{false, allow_values<std::string> {"OIL"}, "EHYSTR(HYSTWET): only OIL supported value is ignored"}}, // wetting_phase_flag
               {9,{false, allow_values<std::string> {"NO"}, "EHYSTR(HYBAKOIL): not used and is ignored"}}, // baker_flag_oil
               {10,{false, allow_values<std::string> {"NO"}, "EHYSTR(HYBAKGAS): not used and is ignored"}}, // baker_flag_gas
               {11,{false, allow_values<std::string> {"NO"}, "EHYSTR(HYBAKWAT): not used and is ignored"}}, // baker_flag_water
            },
         },
         {
            "ENDSCALE",
            {
               {1,{true, allow_values<std::string> {"NODIR"}, "ENDSCALE(DIRECT): only default value of NODIR supported"}}, // DIRECT
               {2,{true, allow_values<std::string> {"REVERS"}, "ENDSCALE(IRREVERS): only default value of REVERS supported"}}, // IRREVERS
            },
         },
         {
            "EQLOPTS",
            {
               {1,{true, allow_values<std::string> {"THPRES"}, "EQLOPTS(MOBILE/QUIESC/IRREVER): options not supported"}}, // OPTION1
               {2,{true, allow_values<std::string> {"THPRES"}, "EQLOPTS(MOBILE/QUIESC/IRREVER): options not supported"}}, // OPTION2
               {3,{true, allow_values<std::string> {"THPRES"}, "EQLOPTS(MOBILE/QUIESC/IRREVER): options not supported"}}, // OPTION3
               {4,{true, allow_values<std::string> {"THPRES"}, "EQLOPTS(MOBILE/QUIESC/IRREVER): options not supported"}}, // OPTION4
            },
         },
         {
            "FOAMOPTS",
            {
               {2,{true, allow_values<std::string> {"TAB"}, "FOAMOPTS(FOAMOPT2): only the default option of TAB is supported"}}, // MODEL
            },
         },
         {
             "FAULTS",
             {
                 {1, {false, [](const std::string& val){ return val.size()<=8;},
                      "FAULTS(FLTNAME): Only names of faults up to 8 characters are supported. Will ignore excess characters."
                     }
                 },},
         },
         {
            "GCONINJE",
            {
               {10,{true, allow_values<std::string> {"RATE", "NETV", "RESV", "VOID"}, "GCONINJE(GUIPHASE): only RATE/NETV/RESV/VOID are supported - will STOP"}}, // GUIDE_RATE_DEF
            },
         },
         {
            "GCONPROD",
            {
               {2,{true, allow_values<std::string> {"NONE", "FLD", "ORAT", "WRAT", "GRAT", "LRAT", "RESV"}, "GCONPROD(TARGET): valid option should be NONE/FLD/ORAT/WRAT/GRAT/LRAT or RESV"}}, // CONTROL_MODE
               {7,{true, allow_values<std::string> {"NONE", "RATE"}, "GCONPROD(ACTION): Only NONE and RATE are supported"}}, 
               {11,{true, allow_values<std::string> {"NONE", "RATE"}, "GCONPROD(ACTWAT): Only NONE and RATE are supported"}}, // WATER_EXCEED_PROCEDURE
               {12,{true, allow_values<std::string> {"NONE", "RATE"}, "GCONPROD(ACTGAS): Only NONE and RATE are supported"}}, // GAS_EXCEED_PROCEDURE
               {13,{true, allow_values<std::string> {"NONE", "RATE"}, "GCONPROD(ACTLIQ): Only NONE and RATE are supported"}}, // LIQUID_EXCEED_PROCEDURE
               {21,{true, allow_values<std::string> {"NONE"}, "GCONPROD(COMBPROC): linearly combined procedure is not used and should be defaulted (1*)"}}, // LIN_TARGET_EXCEED_PROCEDURE
            },
         },
         {
            "GECON",
            {
               {7,{true, allow_values<std::string> {"NONE"}, "GECON(WORKOVER): Workover procedures not implemented"}},
               {8,{true, allow_values<std::string> {"NO"}, "GECON(ENDRUN): End run not implemented"}},
            },
         },
         {
            "GEFAC",
            {
               {3,{true, allow_values<std::string> {"YES"}, "GEFAC(GRPNETWK): Extended Network Model efficiency NO option not implemented"}}, // TRANSFER_EXT_NET
            },
         },
         {
            "GRIDOPTS",
            {
               {1,{true, allow_values<std::string> {"NO", "YES"}, "GRIDOPTS(TRANMULT): should be set to either NO or YES"}}, // TRANMULT
            },
         },
         {
            "GRUPNET",
            {
               {5,{true, allow_values<std::string> {"NO"}, "GRUPNET(SUBSEAMANIFOLD): only option NO is supported"}}, // SUB_SEA_MANIFOLD
               {6,{true, allow_values<std::string> {"NO", "FLO"}, "GRUPNET(LIFTGAS): only option NO and FLO are supported"}}, // ADD_GAS_LIFT_GAS
               {7,{true, allow_values<std::string> {"NONE"}, "GRUPNET(ALQ-DEN): only option NONE is supported"}}, // ALQ_SURFACE_DENSITY
            }
         },
         {
            "GUIDERAT",
            {
               {2,{true, allow_values<std::string> {"OIL", "LIQ", "GAS", "RES", "NONE"}, "GUIDERAT(PHASE): unsupported option must be OIL LIQ GAS RES or NONE"}}, // NOMINATED_PHASE
               {9,{true, allow_values<std::string> {"YES"}, "GUIDERAT(GROPT01): only the default option of YES is supported"}}, // ALLOW_INCREASE
            },
         },
         {
            "MISCIBLE",
            {
               {3,{true, allow_values<std::string> {"NONE"}, "MISCIBLE(MISOPT): only option NONE is supported"}}, // TWOPOINT
            },
         },
         {
             "MULTFLT",
             {
                 {1, {false, [](const std::string& val){ return val.size()<=8;},
                      "MLTFLT(FLTNAME): Only names of faults up to 8 characters are supported. Will ignore excess characters."
                     }
                 },
             },
         },
         {
            "MULTIREG",
            {
               {4,{true, allow_values<std::string> {"F", "M", "O"}, "MULTIREG(REGION_NAME): must equal to F/M or O"}}, // NOT_DEFINED
            },
         },
         {
            "MULTREGP",
            {
               {3,{true, allow_values<std::string> {"F", "M", "O"}, "MULTREGP(REGION_NAME): must equal to F/M or O"}}, // REGION_TYPE
            },
         },
         {
            "MULTREGT",
            {
               {6,{true, allow_values<std::string> {"F", "M", "O"}, "MULTREGT(REGION_NAME): must equal to F/M or O"}}, // REGION_DEF
            },
         },
         {
            "NNC",
            {
               {12,{true, allow_values<std::string> {}, "NNC(FACE1): not supported use 1*."}}, // VE_FACE1
               {13,{true, allow_values<std::string> {}, "NNC(FACE2): not supported use 1*."}}, // VE_FACE2
            },
         },
         {
            "NODEPROP",
            {
               {6,{true, allow_values<std::string> {}, "NODEPROP(GRPNAME): not supported use default"}}, // SOURCE_SINK_GROUP
               {7,{true, allow_values<std::string> {}, "NODEPROP(NETTYPE): not supported use default"}}, // NETWORK_VALUE_TYPE
            },
         },
         {
            "RESTART",
            {
               {3,{true, allow_values<std::string> {}, "RESTART(RSTYPE): restart from SAVE file not supported"}}, // SAVEFILE
               {4,{true, allow_values<std::string> {"UNFORMATTED"}, "RESTART(RSFORMAT): restart from SAVE file not supported"}}, // SAVEFILE_FORMAT
            },
         },
         {
            "ROCKCOMP",
            {
               {1,{true, allow_values<std::string> {"REVERS", "IRREVERS"}, "ROCKCOMP(ROCKOPT): only the REVERS and IRREVERS options are supported – will STOP"}}, // HYSTERESIS
               {3,{true, allow_values<std::string> {"YES", "NO"}, "ROCKCOMP(WATINOPT): only YES and NO are supported"}}, // WATER_COMPACTION
               {4,{false, allow_values<std::string> {}, "ROCKCOMP(PORTXROP): transmissibility dependent on porosity model is not supported"}}, // PORTXROP
            },
         },
         {
            "ROCKOPTS",
            {
               {1,{true, allow_values<std::string> {"PRESSURE"}, "ROCKOPTS: only the PRESSURE options are supported"}}, 
               {2,{true, allow_values<std::string> {"NOSTORE"}, "ROCKOPTS: only the NOSTORE options are supported"}}, 
               {3,{true, allow_values<std::string> {"PVTNUM", "ROCKNUM"}, "ROCKOPTS: only PVTNUM and ROCKNUM are supported"}}, 
            },
         },
         {
            "RPTRST",
            {
               {1,{false, allow_values<std::string> {"ALLPROPS", "BASIC=1", "BASIC=2", "BASIC=3", "BASIC=4", "BASIC=5", "BASIC=6", "DEN", "KRG", "KRO", "KRW", "RSSAT", "RVSAT", "VISC"}, "RPTRST(RPTRST): invalid option or unsupported integer control format"}}, // MNEMONIC_LIST
            },
         },
         {
            "SATOPTS",
            {
               {1,{true, allow_values<std::string> {"HYSTER", "DIRECT"}, "SATOPTS(IRREVERS/SURFTENS): options not supported"}}, // IRREVERS
               {2,{true, allow_values<std::string> {"HYSTER", "DIRECT"}, "SATOPTS(IRREVERS/SURFTENS): options not supported"}}, // IRREVERS
               {3,{true, allow_values<std::string> {"HYSTER", "DIRECT"}, "SATOPTS(IRREVERS/SURFTENS): options not supported"}}, // IRREVERS
               {4,{true, allow_values<std::string> {"HYSTER", "DIRECT"}, "SATOPTS(IRREVERS/SURFTENS): options not supported"}}, // IRREVERS
            },
         },
         {
            "SPECGRID",
            {
               {5,{true, allow_values<std::string> {"F"}, "SPECGRID(TYPE): only option F (Cartesian grids supported) supported"}}, // COORD_TYPE
            },
         },
         {
            "TABDIMS",
            {
               {20,{false, allow_values<std::string> {}, "TABDIMS(NOTUSED): should be defaulted (1*) - ignored as not used"}}, // ITEM20_NOT_USED
               {25,{false, allow_values<std::string> {}, "TABDIMS(RESVED): should be defaulted (1*) - ignored as not used"}}, // RESERVED
            },
         },
         {
             "THPRESFT",
             {
                 {1, {false, [](const std::string& val){ return val.size()<=8;},
                      "THPRESFT(FLTNAME): Only names of faults up to 8 characters are supported. Will ignore excess characters."
                     }
                 },
             },
         },
         {
            "TRACER",
            {
               {1,{true, [](const std::string& val){ return val.size()<=3;}, "TRACER(NAME): Only names of tracers up to 3 characters are supported."}},
               {4,{true, allow_values<std::string> {}, "TRACER(SOLPHASE): partitioned tracer model not supported use default"}}, // SOLUTION_PHASE
               {6,{true, allow_values<std::string> {}, "TRACER(ADSPHASE): partitioned tracer model not supported use default"}}, // ADSORB_PHASE
            },
         },
         {
            "TRACERS",
            {
               {5,{true, allow_values<std::string> {"NODIFF"}, "TRACERS(DIFFOPT): numerical diffusion control for tracers not implemented."}}, // NUMERIC_DIFF
               {8,{true, allow_values<std::string> {"NO"}, "TRACERS(NONLIN): only linear option NO supported"}}, // PASSIVE_NONLINEAR
            },
         },
         {
            "UDQDIMS",
            {
               {11,{true, allow_values<std::string> {"N"}, "UDQDIMS(RSEED): option is not supported use default"}}, // RESTART_NEW_SEED
            },
         },
         {
            "WAGHYSTR",
            {
               {3,{true, allow_values<std::string> {"YES"}, "WAGHYSTR(GAS_MODEL): only the YES option is supported – will STOP"}}, // GAS_MODEL
               {4,{true, allow_values<std::string> {"NO"}, "WAGHYSTR(RES_OIL): only the NO option is supported – will STOP"}}, // RES_OIL
               {5,{true, allow_values<std::string> {"NO"}, "WAGHYSTR(WATER_MODEL): only the NO option is supported – will STOP"}}, // WATER_MODEL
            },
         },

         {
            "WCONHIST",
            {
               {3,{true, allow_values<std::string> {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP"}, "WCONHIST(TARGET): should be set to ORAT/WRAT/GRAT/LRAT/RESV or BHP"}}, // CMODE
            },
         },
         {
            "WEFAC",
            {
               {3,{true, allow_values<std::string> {"YES"}, "WEFAC(WELNETWK): only the YES option is supported"}}, // EXTENDED_NETWORK_OPT
            },
         },
         {
            "WELSPECS",
            {
               {8,{true, allow_values<std::string> {"STD", "NO"}, "WELSPECS(WELNETWK): only the STD and NO options are supported"}}, // INFLOW_EQ
               {12,{true, allow_values<std::string> {"SEG"}, "WELSPECS(DENOPT): only the SEG option is supported"}}, // DENSITY_CALC
               {14,{true, allow_values<std::string> {}, "WELSPECS(STRMLIN1): not used and should be defaulted"}}, // FRONTSIM1
               {15,{true, allow_values<std::string> {}, "WELSPECS(STRMLIN2): not used and should be defaulted"}}, // FRONTSIM2
               {16,{true, allow_values<std::string> {"STD"}, "WELSPECS(TYPECOMP): not used and should be defaulted"}}, // well_model
            },
         },
         {
            "WELTARG",
            {
               {2,{true, allow_values<std::string> {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP", "THP", "VFP", "LIFT", "GUID"}, "WELTARG(TARGET): invalid option"}}, // CMODE
            },
         },
         {
            "WGRUPCON",
            {
               {4,{true, allow_values<std::string> {"OIL", "WAT", "GAS", "LIQ", "RES", "RAT"}, "WGRUPCON(TARGET): only OIL WAT GAS LIQ RES RAT options are supported"}}, // PHASE
            },
         },
         {
            "WHISTCTL",
            {
               {2,{true, allow_values<std::string> {"NO"}, "WHISTCTL(END): only the NO option is supported"}}, // BPH_TERMINATE
            },
         },
         {
            "WLIFTOPT",
            {
               {7,{true, allow_values<std::string> {"NO"}, "WLIFTOPT(OPTLIMIT): only the default NO option is supported"}}, // ALLOCATE_EXTRA_LIFT_GAS
            },
         },
         {
            "WPAVE",
            {
                {4,{false, allow_values<std::string> {"OPEN"}, "WPAVE(WPAVE4) Connection flag not really supported. Should be OPEN or defaulted."}}, // CONNECTION
            },
         },
         {
            "WTEST",
            {
               {3,{true, allow_values<std::string> {"E", "P", "EP", "PE", ""}, "WTEST(TEST): only the E (economic) and P (physical) reason is currently supported"}}, // REASON
            },
         },
         {
            "WVFPEXP",
            {
               {5,{false, allow_values<std::string> {"WG"}, "WVFPEXP(EXTRAP): only linear extrapolation is support "}}, // EXTRAPOLATION_CONTROL
            },
         },
         {
            "WWPAVE",
            {
                {5,{false, allow_values<std::string> {"OPEN"}, "WWPAVE(WPAVE4) Connection flag not really supported. Should be OPEN or defaulted."}}, // CONNECTION
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
               {8,{true, allow_values<int> {0}, "EDITNNC(ISATNUM1): only default value of 0 supported"}}, // SAT_TABLE12
               {9,{true, allow_values<int> {0}, "EDITNNC(ISATNUM2): only default value of 0 supported"}}, // SAT_TABLE21
               {10,{true, allow_values<int> {0}, "EDITNNC(IPRSNUM1): only default value of 0 supported"}}, // PRESS_TABLE12
               {11,{true, allow_values<int> {0}, "EDITNNC(IPRSNUM2): only default value of 0 supported"}}, // PRESS_TABLE21
            },
         },
         {
            "EHYSTR",
            {
               {2,{true, allow_values<int> {0, 1, 2, 3}, "EHYSTR(HYSTMOD): only Carlson or Killough Hysteresis Models supported (0,1 or 2,3)"}}, // relative_perm_hyst
               {13,{true, allow_values<int> {0}, "EHYSTR(HYSWETRP): Killough’s option not supported and should be defaulted"}}, // FLAG_SOMETHING
            },
         },
         {
            "ENDSCALE",
            {
               {3,{false, allow_values<int> {1}, "ENDSCALE(NTENDP): depth end-point scaling not supported – value ignored"}}, // NTENDP
               {4,{false, allow_values<int> {20}, "ENDSCALE(NNODES): depth end-point scaling not supported – value ignored"}}, // NSENDP
               {5,{true, allow_values<int> {0}, "ENDSCALE(MODE): depth temperature end-point scaling not supported"}}, // COMP_MODE
            },
         },
         {
            "EQLDIMS",
            {
               {4,{true, allow_values<int> {1}, "EQLDIMS(NTTRVD): tracer regions (TNUM) not supported. This item should be defaulted to 1."}}, // NTTRVD
            },
         },
         {
            "EQUIL",
            {
               {9,{true, [](int x) { return x >= -20 && x <= 0; }, "EQUIL(EQLOPT3): only values less than or equal to zero are supported (default is -5)"}}, // OIP_INIT
               {10,{false, allow_values<int> {}, "EQUIL(EQLOPT4): compositional option not used, should be defaulted"}}, // EQLOPT4
               {11,{false, allow_values<int> {}, "EQUIL(EQLOPT5): compositional option not used, should be defaulted"}}, // EQLOPT5
            },
         },
         {
            "FOAMROCK",
            {
               {1,{true, allow_values<int> {1}, "FOAMROCK(ADINDX): only the default(1) value is supported"}}, // ADSORPTION_INDEX
            },
         },
         {
            "GRIDFILE",
            {
               {1,{false, allow_values<int> {0}, "GRIDFILE(NGRID): output of GRID file is not supported."}}, // GRID
               {2,{false, allow_values<int> {0,1}, "GRIDFILE(NEGRID): only generate (1) or not (0) an EGRID file is supported."}}, // EGRID
            },
         },
         {
            "MESSAGES",
            {
               {7,{false, allow_values<int> {1000000}, "MESSAGES(STOPMESG): option is not supported"}}, // MESSAGE_STOP_LIMIT
               {8,{false, allow_values<int> {1000000}, "MESSAGES(STOPCOMT): option is not supported"}}, // COMMENT_STOP_LIMIT
               {9,{false, allow_values<int> {10000}, "MESSAGES(STOPWARN): option is not supported"}}, // WARNING_STOP_LIMIT
               {10,{false, allow_values<int> {100}, "MESSAGES(STOPPROB): option is not supported"}}, // PROBLEM_STOP_LIMIT
               {11,{false, allow_values<int> {10}, "MESSAGES(STOPERRS): option is not supported"}}, // ERROR_STOP_LIMIT
               {12,{false, allow_values<int> {1}, "MESSAGES(STOPBUGS): option is not supported"}}, // BUG_STOP_LIMIT
               {13,{false, allow_values<int> {10}, "MESSAGES(PRTGRPMS): option is not supported"}}, // GROUP_PRINT_LIMIT
            },
         },
         {
            "NETBALAN",
            {
               {5,{false, allow_values<int> {10}, "NETBALAN(THPMXITE): option is not supported"}}, // MAX_ITER_THP
            },
         },
         {
            "NETWORK",
            {
               {3,{false, allow_values<int> {20}, "NETWORK(NBCMAX): option is not used and should be defaulted– value ignored"}}, // NBCMAX
            },
         },
         {
            "NNC",
            {
               {8,{true, allow_values<int> {0}, "NNC(ISATNUM1): only default value of 0 supported"}}, // SIM_DEPENDENT1
               {9,{true, allow_values<int> {0}, "NNC(ISATNUM2): only default value of 0 supported"}}, // SIM_DEPENDENT2
               {10,{true, allow_values<int> {0}, "NNC(IPRSNUM1): only default value of 0 supported"}}, // PRESSURE_TABLE1
               {11,{true, allow_values<int> {0}, "NNC(IPRSNUM2): only default value of 0 supported"}}, // PRESSURE_TABLE2
            },
         },
         {
            "NUMRES",
            {
               {1,{true, allow_values<int> {1}, "NUMRES(NUMRES): only a value of one is supported – will STOP"}}, // NUM
            },
         },
         {
            "REGDIMS",
            {
               {5,{false, allow_values<int> {0}, "REGDIMS(NUSREG): compositional TRACK regions not supported - value ignored"}}, // "MAX_ETRACK
               {6,{false, allow_values<int> {1}, "REGDIMS(NTCREG): COAL regions not supported - value ignored"}}, // NTCREG
               {8,{false, allow_values<int> {0}, "REGDIMS(NWKDREG): should be equal to 0 - value ignored"}}, // MAX_OPERATE_DWORK
               {9,{false, allow_values<int> {0}, "REGDIMS(NWKIREG): should be equal to 0 - value ignored"}}, // MAX_OPERATE_IWORK
            },
         },
         {
            "SPECGRID",
            {
               {4,{true, allow_values<int> {1}, "SPECGRID(NUMRES): must be equal to 1"}}, // NUMRES
            },
         },
         {
            "TABDIMS",
            {
               {7,{false, allow_values<int> {20}, "TABDIMS(NRVPVT): should be defaulted (20) – ignored as not used"}}, // MAX_RV_NODES"
               {9,{false, allow_values<int> {1}, "TABDIMS(NMEOSR): must be greater than or equal to 1"}}, // NUM_EOS_RES
               {10,{false, allow_values<int> {1}, "TABDIMS(NMEOSS): should be equal to 1 - ignored as not used"}}, // NUM_EOS_SURFACE
               {12,{false, allow_values<int> {1}, "TABDIMS(MXNTHR): should be equal to 1 - ignored as not used"}}, // MAX_THERMAL_REGIONS
               {14,{false, allow_values<int> {0}, "TABDIMS(MXNPMR): should be equal to 0 - ignored as not used"}}, // MAX_PRESSURE_MAINTAINANCE_REGIONS
               {15,{false, allow_values<int> {0}, "TABDIMS(NTABKT): should be defaulted (0) – ignored as not used"}}, // MAX_KVALUE_TABLES
               {16,{false, allow_values<int> {0}, "TABDIMS(NTALPHA): should be defaulted (0) - ignored as not used"}}, // NTALPHA
               {17,{false, allow_values<int> {10}, "TABDIMS(NASPKA): should be defaulted (10) - ignored as not used"}}, // ASPHALTENE_ASPKDAM_MAX_ROWS
               {18,{false, allow_values<int> {10}, "TABDIMS(MXRAWG): should be defaulted (10) - ignored as not used"}}, // ASPHALTENE_ASPREWG_MAX_ROWS
               {19,{false, allow_values<int> {10}, "TABDIMS(MXRASO): should be defaulted (10) - ignored as not used"}}, // ASPHALTENE_ASPVISO_MAX_ROWS
               {21,{false, allow_values<int> {5}, "TABDIMS(MCASPP): should be defaulted (5) - ignored as not used"}}, // ASPHALTENE_ASPPW2D_MAX_COLUMNS
               {22,{false, allow_values<int> {5}, "TABDIMS(MRASPP): should be defaulted (5) - ignored as not used"}}, // ASPHALTENE_ASPPW2D_MAX_ROWS
               {23,{false, allow_values<int> {5}, "TABDIMS(MXRATF): should be defaulted (5) - ignored as not used"}}, // ASPHALTENE_ASPWETF_MAX_ROWS
               {24,{false, allow_values<int> {0}, "TABDIMS(MXNKVT): should be defaulted (0) - ignored as not used"}}, // NUM_KVALUE_TABLES
            },
         },
         {
            "TRACER",
            {
               {5,{true, allow_values<int> {0}, "TRACER(KPNUM): partitioned tracer model not supported use default"}}, // NUM_PART_TABLE
            },
         },
         {
            "TRACERS",
            {
               {4,{false, allow_values<int> {0}, "TRACERS(MXENVTR):  passive environmental tracers not supported - ignored as not used"}}, // MAX_ENV_TRACERS
               {6,{false, allow_values<int> {12}, "TRACERS(MXITRTR): not supported - ignored as not used"}}, // MAX_ITER
               {7,{false, allow_values<int> {1}, "TRACERS(MNITRTR): not supported - ignored as not used"}}, // MIN_ITER
               {9,{false, allow_values<int> {}, "TRACERS(LNCONFAC): not supported - ignored as not used"}}, // ONEOFF_LIN_TIGHT
               {10,{false, allow_values<int> {}, "TRACERS(NLCONFAC): not supported - ignored as not used"}}, // ONEOFF_NLIN_TIGHT
               {12,{false, allow_values<int> {0}, "TRACERS(NUMCONF): not supported - ignored as not used"}}, // NTIGHTFACTORS
            },
         },
         {
            "UDADIMS",
            {
               {2,{false, allow_values<int> {0}, "UDADIMS(IGNORED): should be defaulted (0) – ignored as not used"}}, // IGNORED
            },
         },
         {
            "UDQPARAM",
            {
               {1,{false, allow_values<int> {1}, "UDQPARAM(RSEED): option not supported – value ignored"}}, // RANDOM_SEED
            },
         },
         {
            "WELLDIMS",
            {
               {5,{false, allow_values<int> {5}, "WELLDIMS(MXSTAGE): option not supported – value ignored"}}, // MAX_STAGES
               {6,{false, allow_values<int> {10}, "WELLDIMS(MXSTRMS): option not supported – value ignored"}}, // MAX_STREAMS
               {7,{false, allow_values<int> {5}, "WELLDIMS(MXMIS): option not supported – value ignored"}}, // MAX_MIXTURES
               {8,{false, allow_values<int> {4}, "WELLDIMS(MXSEPS): option not supported – value ignored"}}, // MAX_SEPARATORS
               {9,{false, allow_values<int> {3}, "WELLDIMS(MXCOMPS): option not supported – value ignored"}}, // MAX_MIXTURE_ITEMS
               {10,{false, allow_values<int> {0}, "WELLDIMS(MXDOCOMP): option not supported – value ignored"}}, // MAX_COMPLETION_X
               {11,{false, allow_values<int> {1}, "WELLDIMS(MXWSLIST): option not supported – value ignored"}}, // MAX_WELLIST_PR_WELL
               {12,{false, allow_values<int> {1}, "WELLDIMS(MXWLISTS): option not supported – value ignored"}}, // MAX_DYNAMIC_WELLIST
               {13,{false, allow_values<int> {10}, "WELLDIMS(MXWSECD): option not supported – value ignored"}}, // MAX_SECONDARY_WELLS
               {14,{false, allow_values<int> {201}, "WELLDIMS(MXNGPP): option not supported – value ignored"}}, // NO_JASON_ENTRY
            },
         },
         {
            "WELSPECS",
            {
               {17,{true, allow_values<int> {0}, "WELSPECS(POLYTAB): only the default value of zero is supported"}}, // POLYMER_TABLE
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
               {12,{false, allow_values<double> {1}, "AQUCON(VEOPT1): Vertical Equilibrium option one– not used"}}, // VEFRAC
               {13,{false, allow_values<double> {1}, "AQUCON(VEOPI2): Vertical Equilibrium option two– not used"}}, // VEFRACP
            },
         },
         {
            "AQUCT",
            {
               {12,{false, allow_values<double> {0}, "AQUCT(SALTCON): option is not used and should be defaulted – value ignored"}}, // INI_SALT
            },
         },
         {
            "AQUFETP",
            {
               {8,{false, allow_values<double> {0}, "AQUFETP(SALTCON): option is not used and should be defaulted – value ignored"}}, // SALINITY
               {9,{false, allow_values<double> {}, "AQUFETP(TEMP): option is not used and should be defaulted – value ignored"}}, // TEMP
            },
         },
         {
            "DIFFC",
            {
               {7,{false, allow_values<double> {0}, "DIFFC(ISATNUM1): only default value of 0 supported – will continue"}}, // GAS_OIL_CROSS_DIFF_COEFF
               {8,{false, allow_values<double> {0}, "DIFFC(ISATNUM2): only default value of 0 supported – will continue"}}, // OIL_OIL_CROSS_DIFF_COEFF
            },
         },
         {
            "EDITNNC",
            {
               {14,{true, allow_values<double> {0}, "EDITNNC(DIFFNNC): not supported and should be defaulted"}}, // DIFFM
            },
         },
         {
            "EHYSTR",   
            {
               {3,{false, allow_values<double> {1.0}, "EHYSTR(HYSTREL): Killough’s option not supported and should be defaulted"}}, // curvature_param_killough_wetting
               {12,{true, allow_values<double> {0}, "EHYSTR(HYTHRESH): Killough’s option not supported and should be defaulted"}}, // threshold_saturation
            },
         },
         {
            "FOAMFSC",
            {
               {4,{false, allow_values<double> {1E-6}, "FOAMFSC(MINSWAT): option is not supported – value ignored"}}, // MIN_WAT_SAT
            },
         },
         {
            "GCONINJE",
            {
               {13,{false, allow_values<double> {}, "GCONINJE(WGASRATE): wet gas rate is not used and should be defaulted (1*)"}}, // WETGAS_TARGET
            },
         },
         {
            "GCONPROD",
            {
               {15,{true, allow_values<double> {}, "GCONPROD(RESVFRAC): reservoir volume fraction is not supported and should be defaulted (1*)"}}, // RESERVOIR_VOLUME_BALANCE
               {16,{false, allow_values<double> {}, "GCONPROD(WGASRATE): wet gas rate is not used and should be defaulted (1*)"}}, // WETGAS_TARGET
               {17,{false, allow_values<double> {}, "GCONPROD(CALRATE): calorific rate is not used and should be defaulted (1*)"}}, // CALORIFIC_TARGET
               {18,{false, allow_values<double> {}, "GCONPROD(GASFRAC): gas production fraction is not used and should be defaulted (1*)"}}, // SURFACE_GAS_FRACTION
               {19,{false, allow_values<double> {}, "GCONPROD(WATFRAC): water production fraction is not used and should be defaulted (1*)"}}, // SURFACE_WAT_FRACTION
               {20,{true, allow_values<double> {}, "GCONPROD(COMBRATE): linearly combined rate is not used and should be defaulted (1*)"}}, // LINEAR_COMBINED_TARGET
            },
         },
         {
            "GUIDERAT",
            {
               {10,{true, [](double x) { return x >= 0; }, "GUIDERAT(GROPT01): only only positive values allowed – will STOP"}}, // ALLOW_INCREASE
            },
         },
         {
            "MULTFLT",
            {
               {3,{true, allow_values<double> {}, "MULTFLT(FLT-DIF): the diffusivity multiplier option is not supported"}}, // NOT_DEFINED
            },
         },
         {
            "NETBALAN",
            {
               {1,{true, [](const double value) { return value <= 0.0; }, "NETBALAN(NSTEP): only negative values or 0 supported"}}, // TIME_INTERVAL
               {4,{false, allow_values<double> {0.01}, "NETBALAN(GRPCNV): not supported"}}, // THP_CONVERGENCE_LIMIT
               {6,{false, allow_values<double> {1e20}, "NETBALAN(NTRGERR): not supported"}}, // TARGET_BALANCE_ERROR
               {7,{false, allow_values<double> {1e20}, "NETBALAN(NMAXERR): not supported"}}, // MAX_BALANCE_ERROR
               {8,{false, allow_values<double> {}, "NETBALAN(NTSMIN): not supported"}}, // MIN_TIME_STEP
            },
         },
         {
            "NNC",
            {
               {14,{true, allow_values<double> {0}, "NNC(DIFFNNC): not supported"}}, // DIFFUSIVITY
               {15,{true, allow_values<double> {0}, "NNC(DISPNNC): not supported"}}, // SIM_DEPENDENT3
               {16,{true, allow_values<double> {0}, "NNC(AREANNC): not supported"}}, // VDFLOW_AREA
               {17,{true, allow_values<double> {0}, "NNC(PERMNNC): not supported"}}, // VDFLOW_PERM
            },
         },
         {
            "PLYMAX",  
            {
               {2,{false, allow_values<double> {}, "PLYMAX(SALTCON): option is ignored since BRINE and POLYMER combination is not implemented in OPM Flow"}}, // MAX_SALT_CONCENTRATION
            },
         },
         {
            "ROCKCOMP",
            {
               {5,{false, allow_values<double> {0}, "ROCKCOMP(CARKZEXP): transmissibility dependent on porosity model is not supported"}}, // CARKZEXP
            },
         },
         {
            "TRACERS",
            {
               {11,{false, allow_values<double> {1.0}, "TRACERS(CONFAC): not supported - ignored as not used"}}, // TIGHTENING_FACTORS
            },
         },
         {
            "VISCREF",
            {
               {3,{true, allow_values<double> {}, "VISCREF(API): API tracking option is not supported"}}, // API_GRAVITY
            },
         },
         {
            "WAGHYSTR",
            {
               {8,{false, allow_values<double> {}, "WAGHYSTR(RES_OIL_MOD_FRACTION): Residual oil modification for STONE1 not supported - value ignored"}}, // RES_OIL_MOD_FRACTION
            },
         },
         {
            "WCONHIST",
            {
               {11,{false, [](double x) { return x == 0; }, "WCONHIST(WGRA): wet gas rate is not supported use default"}}, // WGASRAT
               {12,{false, [](double x) { return x == 0; }, "WCONHIST(NGL): natural gas rate is not supported use default"}}, // NGLRAT
            },
         },
         {
            "WCONINJE",
            {
               {11,{true, [](double x) { return x == 0; }, "WCONINJE(RSSTEAM): is not supported and is ignored"}}, // GAS_STEAM_RATIO
               {12,{true, [](double x) { return x == 0; }, "WCONINJE(OILFRAC): is not supported and is ignored"}}, // SURFACE_OIL_FRACTION
               {13,{true, [](double x) { return x == 0; }, "WCONINJE(WATFRAC): is not supported and is ignored"}}, // SURFACE_WATER_FRACTION
               {14,{true, [](double x) { return x == 0; }, "WCONINJE(GASFRAC): is not supported and is ignored"}}, // SURFACE_GAS_FRACTION
               {15,{true, [](double x) { return x == 0; }, "WCONINJE(OILSTEAM): is not supported and is ignored"}}, // OIL_STEAM_RATIO
            },
         },
         {
            "WCONINJH",
            {
               {8,{true, [](double x) { return x == 0; }, "WCONINJH(RSRVINJ): is not supported and is ignored"}}, // VAPOIL_C
               {9,{true, [](double x) { return x == 0; }, "WCONINJH(OILFRAC): is not supported and is ignored"}}, // SURFACE_OIL_FRACTION
               {10,{true, [](double x) { return x == 0; }, "WCONINJH(WATFRAC): is not supported and is ignored"}}, // SURFACE_WATER_FRACTION
               {11,{true, [](double x) { return x == 0; }, "WCONINJH(GASFRAC): is not supported and is ignored"}}, // SURFACE_GAS_FRACTION
            },
         },
         {
            "WCONPROD",
            {
               {13,{false, allow_values<double> {}, "WCONPROD(WGASRATE): wet gas rate is not used and should be defaulted (1*)"}}, // E300_ITEM13
               {14,{false, allow_values<double> {}, "WCONPROD(MOLARATE): molar rate not used and should be defaulted (1*)"}}, // E300_ITEM14
               {15,{false, allow_values<double> {}, "WCONPROD(STEAMRAT): steam rate is not used and should be defaulted (1*)"}}, // E300_ITEM15
               {16,{false, allow_values<double> {}, "WCONPROD(DELTAP): pressure offset not used and should be defaulted (1*)"}}, // E300_ITEM16
               {17,{false, allow_values<double> {}, "WCONPROD(DELTAT): temperature offset not used and should be defaulted (1*)"}}, // E300_ITEM17
               {18,{false, allow_values<double> {}, "WCONPROD(CALRATE): calorific rate not used and should be defaulted (1*)"}}, // E300_ITEM18
               {19,{true, allow_values<double> {}, "WCONPROD(COMBPROC): linearly combined rate not used and should be defaulted (1*)"}}, // E300_ITEM19
               {20,{false, allow_values<double> {}, "WCONPROD(NGL): natural gas liquid rate  is not used and should be defaulted (1*)"}}, // E300_ITEM20
            },
         },
         {
            "WINJTEMP",
            {
               {2,{true, allow_values<double> {1}, "WINJTEMP(STEAMQAL): steam injection is not supported, this item should be defaulted"}}, // STEAM_QUALITY
               {5,{true, allow_values<double> {0}, "WINJTEMP(ENTHALPY): enthalpy of injected fluid is not supported"}}, // ENTHALPY
            },
         },
         {
            "WLIFTOPT",
            {
               {6,{true, allow_values<double> {0}, "WLIFTOPT(OPTGAS): incremental gas weighting not supported"}}, // DELTA_GAS_RATE_WEIGHT_FACTOR
            },
         },
   };

   return partially_supported_keywords_double;
}

} // namespace Opm::FlowKeywordValidation
