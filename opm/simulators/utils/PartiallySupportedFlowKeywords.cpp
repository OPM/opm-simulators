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
            "BRANPROP",
            {
               {5,{false, allow_values<std::string> {"NONE"}, "BRANPROP(ALQ-DEN): option is not supported and should be defaulted (1*)"}}, // ALQ_SURFACE_DENSITY
            },
         },
         {
            "COMPORD",
            {
               {2,{true, allow_values<std::string> {"DEPTH", "INPUT", "TRACK"}, "COMPORD(COMPORD): should equal DEPTH INPUT or TRACK – will STOP"}}, // ORDER_TYPE
            },
         },
         {
            "EDITNNC",
            {
               {12,{false, allow_values<std::string> {}, "EDITNNC(FACE1): not supported use 1* - will continue"}}, // FACE_FLOW12
               {13,{false, allow_values<std::string> {}, "EDITNNC(FACE2): not supported use 1* - will continue"}}, // FACE_FLOW21
            },
         },
         {
            "EHYSTR",
            {
               {5,{true, allow_values<std::string> {"KR"}, "EHYSTR(HYSTOPT): only the KR relative permeability hysteresis option is supported – will STOP"}}, // limiting_hyst_flag
               {6,{false, allow_values<std::string> {}, "EHYSTR(HYSTSCAN): Killough’s option not supported and is ignored"}}, // shape_cap_press_flag
               {7,{false, allow_values<std::string> {}, "EHYSTR(HYSTMOB): mobility option not supported and is ignored"}}, // init_fluid_mob_flag
               {8,{false, allow_values<std::string> {"OIL"}, "EHYSTR(HYSTWET): only OIL supported value is ignored"}}, // wetting_phase_flag
               {9,{false, allow_values<std::string> {"NO"}, "EHYSTR(HYBAKOIL): not used and is ignored"}}, // baker_flag_oil
               {10,{false, allow_values<std::string> {"NO"}, "EHYSTR(HYBAKGAS): not used and is ignored"}}, // baker_flag_gas
               {11,{false, allow_values<std::string> {"NO"}, "EHYSTR(HYBAKWAT): not used and is ignored"}}, // baker_flag_water
            },
         },
         {
            "ENDSCALE",
            {
               {1,{false, allow_values<std::string> {"NODIR"}, "ENDSCALE(DIRECT): only default value of NODIR supported"}}, // DIRECT
               {2,{false, allow_values<std::string> {"REVERS"}, "ENDSCALE(IRREVERS): only default value of REVERS supported"}}, // IRREVERS
            },
         },
         {
            "EQLOPTS",
            {
               {1,{false, allow_values<std::string> {"THPRES"}, "EQLOPTS(MOBILE/QUIESC/IRREVER): options not supported – value ignored"}}, // OPTION1
               {2,{false, allow_values<std::string> {"THPRES"}, "EQLOPTS(MOBILE/QUIESC/IRREVER): options not supported – value ignored"}}, // OPTION2
               {3,{false, allow_values<std::string> {"THPRES"}, "EQLOPTS(MOBILE/QUIESC/IRREVER): options not supported – value ignored"}}, // OPTION3
               {4,{false, allow_values<std::string> {"THPRES"}, "EQLOPTS(MOBILE/QUIESC/IRREVER): options not supported – value ignored"}}, // OPTION4
            },
         },
         {
            "FOAMOPTS",
            {
               {1,{false, allow_values<std::string> {"GAS"}, "FOAMOPTS(FOAMOPT1): only the default option of GAS is supported – value ignored"}}, // TRANSPORT_PHASE
               {2,{false, allow_values<std::string> {"TAB"}, "FOAMOPTS(FOAMOPT2): only the default option of TAB is supported – value ignored"}}, // MODEL
            },
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
               {2,{true, allow_values<std::string> {"NONE", "FLD", "ORAT", "WRAT", "GRAT", "LRAT", "RESV"}, "GCONPROD(TARGET): valid option should be NONE/FLD/ORAT/WRAT/GRAT/LRAT or RESV – will STOP"}}, // CONTROL_MODE
               {11,{true, allow_values<std::string> {}, "GCONPROD(ACTWAT): water violation procedure not implemented – will STOP"}}, // WATER_EXCEED_PROCEDURE
               {12,{true, allow_values<std::string> {}, "GCONPROD(ACTGAS): gas violation procedure not implemented – will STOP"}}, // GAS_EXCEED_PROCEDURE
               {13,{true, allow_values<std::string> {}, "GCONPROD(ACTLIQ): liquid violation procedure not implemented – will STOP"}}, // LIQUID_EXCEED_PROCEDURE
               {21,{false, allow_values<std::string> {}, "GCONPROD(COMBPROC): linearly combined procedure is not used and should be defaulted (1*) – will continue"}}, // LIN_TARGET_EXCEED_PROCEDURE
            },
         },
         {
            "GEFAC",
            {
               {3,{true, allow_values<std::string> {"YES"}, "GEFAC(GRPNETWK): Extended Network Model efficiency NO option not implemented – will STOP"}}, // TRANSFER_EXT_NET
            },
         },
         {
            "GRIDOPTS",
            {
               {1,{true, allow_values<std::string> {"NO", "YES"}, "GRIDOPTS(TRANMULT): should be set to either NO or YES – will STOP"}}, // TRANMULT
            },
         },
         {
            "GUIDERAT",
            {
               {2,{true, allow_values<std::string> {"OIL", "LIQ", "GAS", "RES", "NONE"}, "GUIDERAT(PHASE): unsupported option must be OIL LIQ GAS RES or NONE – will STOP"}}, // NOMINATED_PHASE
               {9,{true, allow_values<std::string> {"YES"}, "GUIDERAT(GROPT01): only the default option of YES supported – will STOP"}}, // ALLOW_INCREASE
            },
         },
         {
            "MISCIBLE",
            {
               {3,{false, allow_values<std::string> {"NONE"}, "MISCIBLE(MISOPT): only option NONE is supported – value ignored"}}, // TWOPOINT
            },
         },
         {
            "MULTIREG",
            {
               {4,{false, allow_values<std::string> {"F", "M", "O"}, "MULTIREG(REGION_NAME): must equal to F/M or O"}}, // NOT_DEFINED
            },
         },
         {
            "MULTREGP",
            {
               {3,{false, allow_values<std::string> {"F", "M", "O"}, "MULTREGP(REGION_NAME): must equal to F/M or O"}}, // REGION_TYPE
            },
         },
         {
            "MULTREGT",
            {
               {6,{false, allow_values<std::string> {"F", "M", "O"}, "MULTREGT(REGION_NAME): must equal to F/M or O"}}, // REGION_DEF
            },
         },
         {
            "NNC",
            {
               {12,{false, allow_values<std::string> {}, "NNC(FACE1): not supported use 1* - will continue"}}, // VE_FACE1
               {13,{false, allow_values<std::string> {}, "NNC(FACE2): not supported use 1* - will continue"}}, // VE_FACE2
            },
         },
         {
            "NODEPROP",
            {
               {6,{false, allow_values<std::string> {}, "NODEPROP(GRPNAME): not used – will continue"}}, // SOURCE_SINK_GROUP
               {7,{false, allow_values<std::string> {}, "NODEPROP(NETTYPE): not used – will continue"}}, // NETWORK_VALUE_TYPE
            },
         },
         {
            "RESTART",
            {
               {3,{true, allow_values<std::string> {}, "RESTART(RSTYPE): restart from SAVE file not supported – will STOP"}}, // SAVEFILE
               {4,{true, allow_values<std::string> {"UNFORMATTED"}, "RESTART(RSFORMAT): restart from SAVE file not supported – will STOP"}}, // SAVEFILE_FORMAT
            },
         },
         {
            "ROCKCOMP",
            {
               {1,{false, allow_values<std::string> {"REVERS"}, "ROCKCOMP(ROCKOPT): only the REVERS option is supported"}}, // HYSTERESIS
               {3,{false, allow_values<std::string> {"YES"}, "ROCKCOMP(WATINOPT): only equal to YES is supported"}}, // WATER_COMPACTION
               {4,{false, allow_values<std::string> {}, "ROCKCOMP(PORTXROP): transmissibility dependent on porosity model is not supported"}}, // PORTXROP
            },
         },
         {
            "RPTRST",
            {
               {1,{false, allow_values<std::string> {"ALLPROPS", "BASIC=1", "BASIC=2", "BASIC=3", "BASIC=4", "BASIC=5", "BASIC=6", "DEN", "KRG", "KRO", "KRW", "RSSAT", "RVSAT", "VISC"}, "RPTRST(RPTRST): invalid option or unsupported integer control format – will continue"}}, // MNEMONIC_LIST 
            },
         },
         {
            "SATOPTS",
            {
               {1,{false, allow_values<std::string> {"HYSTER"}, "SATOPTS(DIRECT): directional relative permeability assignment option not supported - value ignored"}}, // options
               {2,{false, allow_values<std::string> {"HYSTER"}, "SATOPTS(IRREVERS): reversible directional relative permeability option not supported – value ignored"}}, // IRREVERS
               {3,{false, allow_values<std::string> {"HYSTER"}, "SATOPTS(HYSTER): hysteresis directional relative permeability option not supported - value ignored"}}, // HYSTER
               {4,{false, allow_values<std::string> {"HYSTER"}, "SATOPTS(SURFTENS): capillary pressure surface tension pressure dependency option not supported – value ignored"}}, // SURFTENS
            },
         },
         {
            "SPECGRID",
            {
               {5,{true, allow_values<std::string> {"F"}, "SPECGRID(TYPE): only option F (Cartesian grids supported) supported – will STOP"}}, // COORD_TYPE
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
            "TRACER",
            {
               {4,{false, allow_values<std::string> {}, "TRACER(SOLPHASE): partitioned tracer model not supported - ignored as not used"}}, // SOLUTION_PHASE
               {6,{false, allow_values<std::string> {}, "TRACER(ADSPHASE): partitioned tracer model not supported - ignored as not used"}}, // ADSORB_PHASE
            },
         },
         {
            "TRACERS",
            {
               {5,{false, allow_values<std::string> {"NODIFF"}, "TRACERS(DIFFOPT): numerical diffusion control not implemented. - ignored as not used"}}, // NUMERIC_DIFF
               {8,{false, allow_values<std::string> {"NO"}, "TRACERS(NONLIN): only linear option NO supported – will continue"}}, // PASSIVE_NONLINEAR
            },
         },
         {
            "UDQDIMS",
            {
               {11,{false, allow_values<std::string> {"N"}, "UDQDIMS(RSEED): option is not supported – value ignored"}}, // RESTART_NEW_SEED
            },
         },
         {
            "WCONHIST",
            {
               {3,{true, allow_values<std::string> {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP"}, "WCONHIST(TARGET): should be set to ORAT/WRAT/GRAT/LRAT/RESV or BHP – will STOP"}}, // CMODE
            },
         },
         {
            "WEFAC",
            {
               {3,{true, allow_values<std::string> {"YES"}, "WEFAC(WELNETWK): only the YES option is supported – will STOP"}}, // EXTENDED_NETWORK_OPT
            },
         },
         {
            "WELSPECS",
            {
               {8,{true, allow_values<std::string> {"STD", "NO"}, "WELSPECS(WELNETWK): only the STD and NO options are supported – will STOP"}}, // INFLOW_EQ
               {12,{false, allow_values<std::string> {"SEG"}, "WELSPECS(DENOPT): only the SEG option is supported – will continue"}}, // DENSITY_CALC
               {14,{false, allow_values<std::string> {}, "WELSPECS(STRMLIN1): not used  – will continue"}}, // FRONTSIM1
               {15,{false, allow_values<std::string> {}, "WELSPECS(STRMLIN2): not used – will continue"}}, // FRONTSIM2
               {16,{false, allow_values<std::string> {"STD"}, "WELSPECS(TYPECOMP): not used – will continue"}}, // well_model
            },
         },
         {
            "WELTARG",
            {
               {2,{true, allow_values<std::string> {"ORAT", "WRAT", "GRAT", "LRAT", "RESV", "BHP", "THP", "VFP", "LIFT", "GUID"}, "WELTARG(TARGET): invalid option – will STOP"}}, // CMODE
            },
         },
         {
            "WGRUPCON",
            {
               {4,{true, allow_values<std::string> {"OIL", "WAT", "GAS", "LIQ", "RES", "RAT"}, "WGRUPCON(TARGET): only OIL WAT GAS LIQ RES RAT options are supported – will STOP"}}, // PHASE
            },
         },
         {
            "WHISTCTL",
            {
               {2,{false, allow_values<std::string> {"NO"}, "WHISTCTL(END): only the NO option is supported – will continue"}}, // BPH_TERMINATE
            },
         },
         {
            "WLIFTOPT",
            {
               {7,{false, allow_values<std::string> {"NO"}, "WLIFTOPT(OPTLIMIT): only the default NO option is supported – will continue"}}, // ALLOCATE_EXTRA_LIFT_GAS
            },
         },
         {
            "WRFTPLT",
            {
               {4,{false, allow_values<std::string> {"NO"}, "WRFTPLT(MULTISEG): output is currently not supported – will continue"}}, // OUTPUT_SEGMENT
            },
         },
         {
            "WTEST",
            {
               {3,{false, allow_values<std::string> {"E", "P", "EP", "PE", ""}, "WTEST(TEST): only the E (economic) option is currently supported – will continue"}}, // REASON
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
               {8,{false, allow_values<int> {0}, "EDITNNC(ISATNUM1): only default value of 0 supported – will continue"}}, // SAT_TABLE12
               {9,{false, allow_values<int> {0}, "EDITNNC(ISATNUM2): only default value of 0 supported – will continue"}}, // SAT_TABLE21
               {10,{false, allow_values<int> {0}, "EDITNNC(IPRSNUM1): only default value of 0 supported – will continue"}}, // PRESS_TABLE12
               {11,{false, allow_values<int> {0}, "EDITNNC(IPRSNUM2): only default value of 0 supported – will continue"}}, // PRESS_TABLE21
            },
         },
         {
            "EHYSTR",
            {
               {2,{false, allow_values<double> {0, 1}, "EHYSTR(HYSTMOD): only Carlson Hysteresis Models supported (0 or 1)"}}, // relative_perm_hyst
               {13,{false, allow_values<int> {0}, "EHYSTR(HYSWETRP): Killough’s option not supported and is ignored"}}, // FLAG_SOMETHING
            },
         },
         {
            "ENDSCALE",
            {
               {3,{false, allow_values<int> {1}, "ENDSCALE(NTENDP): depth end-point scaling not supported – value ignored"}}, // NTENDP
               {4,{false, allow_values<int> {20}, "ENDSCALE(NNODES): depth end-point scaling not supported – value ignored"}}, // NSENDP
               {5,{false, allow_values<int> {0}, "ENDSCALE(MODE): depth temperature end-point scaling not supported – value ignored"}}, // COMP_MODE
            },
         },
         {
            "EQLDIMS",
            {
               {4,{false, allow_values<int> {1}, "EQLDIMS(NTTRVD): tracer end-point depth scaling not supported – value ignored"}}, // NTTRVD
               {5,{false, allow_values<int> {20}, "EQLDIMS(NSTRVD): tracer end-point depth scaling not supported – value ignored"}}, // NSTRVD
            },
         },
         {
            "EQUIL",
            {
               {9,{true, [](int x) { return x >= -20 && x <= 0; }, "EQUIL(EQLOPT3): only values less than or equal to zero are supported (default is -5) - will STOP"}}, // OIP_INIT
               {10,{false, allow_values<int> {}, "EQUIL(EQLOPT4): compositional option not used – will continue"}}, // EQLOPT4
               {11,{false, allow_values<int> {}, "EQUIL(EQLOPT5): compositional option not used – will continue"}}, // EQLOPT5
            },
         },
         {
            "FOAMROCK",
            {
               {1,{false, allow_values<int> {1}, "FOAMROCK(ADINDX): only the default(1) value is supported – value ignored"}}, // ADSORPTION_INDEX
            },
         },
         {
            "GRIDFILE",
            {
               {1,{false, allow_values<int> {0}, "GRIDFILE(NGRID): only default value of 0 supported – will continue"}}, // GRID
               {2,{false, allow_values<int> {1}, "GRIDFILE(NEGRID): only default value of 1 supported – will continue"}}, // EGRID
            },
         },
         {
            "MESSAGES",
            {
               {7,{false, allow_values<int> {1000000}, "MESSAGES(STOPMESG): option is not supported – will continue"}}, // MESSAGE_STOP_LIMIT
               {8,{false, allow_values<int> {1000000}, "MESSAGES(STOPCOMT): option is not supported – will continue"}}, // COMMENT_STOP_LIMIT
               {9,{false, allow_values<int> {10000}, "MESSAGES(STOPWARN): option is not supported – will continue"}}, // WARNING_STOP_LIMIT
               {10,{false, allow_values<int> {100}, "MESSAGES(STOPPROB): option is not supported – will continue"}}, // PROBLEM_STOP_LIMIT
               {11,{false, allow_values<int> {10}, "MESSAGES(STOPERRS): option is not supported – will continue"}}, // ERROR_STOP_LIMIT
               {12,{false, allow_values<int> {1}, "MESSAGES(STOPBUGS): option is not supported – will continue"}}, // BUG_STOP_LIMIT
               {13,{false, allow_values<int> {10}, "MESSAGES(PRTGRPMS): option is not supported – will continue"}}, // GROUP_PRINT_LIMIT
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
               {8,{false, allow_values<int> {0}, "NNC(ISATNUM1): only default value of 0 supported – will continue"}}, // SIM_DEPENDENT1
               {9,{false, allow_values<int> {0}, "NNC(ISATNUM2): only default value of 0 supported – will continue"}}, // SIM_DEPENDENT2
               {10,{false, allow_values<int> {0}, "NNC(IPRSNUM1): only default value of 0 supported – will continue"}}, // PRESSURE_TABLE1
               {11,{false, allow_values<int> {0}, "NNC(IPRSNUM2): only default value of 0 supported – will continue"}}, // PRESSURE_TABLE2
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
               {4,{false, allow_values<int> {1}, "SPECGRID(NUMRES): must be equal to 1"}}, // NUMRES
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
               {5,{false, allow_values<int> {0}, "TRACER(KPNUM): partitioned tracer model not supported - ignored as not used"}}, // NUM_PART_TABLE
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
               {17,{true, allow_values<int> {0}, "WELSPECS(POLYTAB): only the default value of zero is supported – will STOP"}}, // POLYMER_TABLE
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
            "COMPDAT",
            {
               {12,{false, allow_values<double> {}, "COMPDAT(DFACT): non-Darcy D factor not supported and should be defaulted (1*) – value ignored"}}, // D_FACTOR
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
               {14,{false, allow_values<double> {0}, "EDITNNC(DIFFNNC): not supported – will continue"}}, // DIFFM
            },
         },
         {
            "EHYSTR",
            {
               {1,{false, allow_values<double> {0.1}, "EHYSTR(HYSTRCP): option is not supported – value ignored"}}, // curvature_capillary_pressure_hyst
               {3,{false, allow_values<double> {1.0}, "EHYSTR(HYSTREL): Killough’s option not supported and is ignored"}}, // curvature_param_killough_wetting
               {4,{false, allow_values<double> {0.1}, "EHYSTR(HYSTSGR): Killough’s option not supported and is ignored"}}, // mod_param_trapped
               {12,{false, allow_values<double> {0}, "EHYSTR(HYTHRESH): Killough’s option not supported and is ignored"}}, // threshold_saturation
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
               {13,{false, allow_values<double> {}, "GCONINJE(WGASRATE): wet gas rate is not used and should be defaulted (1*) – will continue"}}, // WETGAS_TARGET
            },
         },
         {
            "GCONPROD",
            {
               {15,{false, allow_values<double> {}, "GCONPROD(RESVFRAC): reservoir volume fraction is not supported and should be defaulted (1*) – will STOP"}}, // RESERVOIR_VOLUME_BALANCE
               {16,{false, allow_values<double> {}, "GCONPROD(WGASRATE): wet gas rate is not used and should be defaulted (1*) – will continue"}}, // WETGAS_TARGET
               {17,{false, allow_values<double> {}, "GCONPROD(CALRATE): calorific rate is not used and should be defaulted (1*) – will continue"}}, // CALORIFIC_TARGET
               {18,{false, allow_values<double> {}, "GCONPROD(GASFRAC): gas production fraction is not used and should be defaulted (1*) – will continue"}}, // SURFACE_GAS_FRACTION
               {19,{false, allow_values<double> {}, "GCONPROD(WATFRAC): water production fraction is not used and should be defaulted (1*) – will continue"}}, // SURFACE_WAT_FRACTION
               {20,{false, allow_values<double> {}, "GCONPROD(COMBRATE): linearly combined rate is not used and should be defaulted (1*) – will continue"}}, // LINEAR_COMBINED_TARGET
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
               {3,{false, allow_values<double> {}, "MULTFLT(FLT-DIF): the diffusivity multiplier option is not supported – will continue"}}, // NOT_DEFINED
            },
         },
         {
            "NNC",
            {
               {14,{false, allow_values<double> {0}, "NNC(DIFFNNC): not supported – will continue"}}, // DIFFUSIVITY
               {15,{false, allow_values<double> {0}, "NNC(DISPNNC): not supported – will continue"}}, // SIM_DEPENDENT3
               {16,{false, allow_values<double> {0}, "NNC(AREANNC): not supported – will continue"}}, // VDFLOW_AREA
               {17,{false, allow_values<double> {0}, "NNC(PERMNNC): only default value of 0 supported – will continue"}}, // VDFLOW_PERM
            },
         },
         {
            "PLYMAX",
            {
               {2,{false, allow_values<double> {}, "PLYMAX(SALTCON): option is ignored if the BRINE keyword is absent – will continue"}}, // MAX_SALT_CONCENTRATION
            },
         },
         {
            "ROCKCOMP",
            {
               {4,{false, allow_values<double> {0}, "ROCKCOMP(CARKZEXP): transmissibility dependent on porosity model is not supported"}}, // CARKZEXP
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
               {3,{false, allow_values<double> {}, "VISCREF(API): API tracking option is not supported - value ignored"}}, // API_GRAVITY
            },
         },
         {
            "WCONHIST",
            {
               {11,{false, [](double x) { return x == 0; }, "WCONHIST(WGRA): wet gas rate is not supported and is ignored – will continue"}}, // WGASRAT
               {12,{false, [](double x) { return x == 0; }, "WCONHIST(NGL): natural gas rate is not supported and is ignored – will continue"}}, // NGLRAT
            },
         },
         {
            "WCONINJE",
            {
               {10,{false, [](double x) { return x == 0; }, "WCONINJE(RSRVINJ): is not supported and is ignored – will continue"}}, // VAPOIL_C
               {11,{false, [](double x) { return x == 0; }, "WCONINJE(RSSTEAM): is not supported and is ignored – will continue"}}, // GAS_STEAM_RATIO
               {12,{false, [](double x) { return x == 0; }, "WCONINJE(OILFRAC): is not supported and is ignored – will continue"}}, // SURFACE_OIL_FRACTION
               {13,{false, [](double x) { return x == 0; }, "WCONINJE(WATFRAC): is not supported and is ignored – will continue"}}, // SURFACE_WATER_FRACTION
               {14,{false, [](double x) { return x == 0; }, "WCONINJE(GASFRAC): is not supported and is ignored – will continue"}}, // SURFACE_GAS_FRACTION
               {15,{false, [](double x) { return x == 0; }, "WCONINJE(OILSTEAM): is not supported and is ignored – will continue"}}, // OIL_STEAM_RATIO
            },
         },
         {
            "WCONINJH",
            {
               {8,{false, [](double x) { return x == 0; }, "WCONINJH(RSRVINJ): is not supported and is ignored – will continue"}}, // VAPOIL_C
               {9,{false, [](double x) { return x == 0; }, "WCONINJH(OILFRAC): is not supported and is ignored – will continue"}}, // SURFACE_OIL_FRACTION
               {10,{false, [](double x) { return x == 0; }, "WCONINJH(WATFRAC): is not supported and is ignored – will continue"}}, // SURFACE_WATER_FRACTION
               {11,{false, [](double x) { return x == 0; }, "WCONINJH(GASFRAC): is not supported and is ignored – will continue"}}, // SURFACE_GAS_FRACTION
            },
         },
         {
            "WCONPROD",
            {
               {13,{false, allow_values<double> {}, "WCONPROD(WGASRATE): wet gas rate is not used and should be defaulted (1*) – will continue"}}, // E300_ITEM13
               {14,{false, allow_values<double> {}, "WCONPROD(MOLARATE): molar rate not used and should be defaulted (1*) – will continue"}}, // E300_ITEM14
               {15,{false, allow_values<double> {}, "WCONPROD(STEAMRAT): steam rate is not used and should be defaulted (1*) – will continue"}}, // E300_ITEM15
               {16,{false, allow_values<double> {}, "WCONPROD(DELTAP): pressure offset not used and should be defaulted (1*) – will continue"}}, // E300_ITEM16
               {17,{false, allow_values<double> {}, "WCONPROD(DELTAT): temperature offset not used and should be defaulted (1*) – will continue"}}, // E300_ITEM17
               {18,{false, allow_values<double> {}, "WCONPROD(CALRATE): calorific rate not used and should be defaulted (1*) – will continue"}}, // E300_ITEM18
               {19,{false, allow_values<double> {}, "WCONPROD(COMBPROC): linearly combined rate not used and should be defaulted (1*) – will continue"}}, // E300_ITEM19
               {20,{false, allow_values<double> {}, "WCONPROD(NGL): natural gas liquid rate  is not used and should be defaulted (1*) – will continue"}}, // E300_ITEM20
            },
         },
         {
            "WINJTEMP",
            {
               {2,{false, allow_values<double> {1}, "WINJTEMP(STEAMQAL): steam injection is not supported – will continue"}}, // STEAM_QUALITY
               {5,{false, allow_values<double> {0}, "WINJTEMP(ENTHALPY): enthalpy of injected fluid is not supported – will continue"}}, // ENTHALPY
            },
         },
         {
            "WLIFTOPT",
            {
               {6,{false, allow_values<double> {0}, "WLIFTOPT(OPTGAS): incremental gas weighting not supported – will continue"}}, // DELTA_GAS_RATE_WEIGHT_FACTOR
            },
         },
   };

   return partially_supported_keywords_double;
}

} // namespace Opm::FlowKeywordValidation
