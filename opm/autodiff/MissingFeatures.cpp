/*
  Copyright 2016 Statoil ASA.

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
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/parser/eclipse/Deck/DeckKeyword.hpp>
#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/autodiff/MissingFeatures.hpp>
#include <unordered_set>
#include <string>


namespace Opm {

namespace MissingFeatures {

    void checkKeywords(DeckConstPtr deck, ParserConstPtr parser)
    {
        // These keywords are supported by opm-parser, but are not supported
        // by flow. For some of them, only part of the options are supported.
        // The list is used to output messages only.
        std::unordered_set<std::string> unsupported_keywords = {
            "ACTDIMS", "ADSALNOD", "API", "AQUCON", "AQUDIMS", "AQUNUM"
            "COMPLUMP", "COMPSEGS", "CONNECTION", "CPR", 
            "DATE", "ECHO", "EDITNNC", "ENDNUM",
            "ENDSKIP", "ENKSRVD", "ENPTVD", "EQLNUM", "EQUALREG",
            "EXCEL", "EXTRAPMS", "FILLEPS", "FIPNUM", "FMTIN",
            "FMTOUT", "FULLIMP", "GDORIENT", "GECON", "GEFAC", "GRIDUNIT", 
            "GRUPNET", "IMKRVD", "IMPES", "IMPTVD", "MAPUNITS",
            "MAXVALUE", "MESSAGES", "MINVALUE", "MONITOR", "MSGFILE",
            "MULT_XYZ", "NETBALAN", "NEXTSTEP", "NOCASC", "NOECHO",
            "NOGGF", "NOINSPEC", "NOMONITO", "NONNC", "NORSSPEC", "NOSIM",
            "NSTACK", "NUMRES", "NUPCOL", "OILVISCT", "OLDTRAN", "OPTIONS", 
            "PARALLEL", "PBVD", "PCG", "PERMXY", "PERMYZ", 
            "PERMZX", "PIMTDIMS", "PIMULTAB", "PLMIXPAR", "PLYADSS", "PLYDHFLF",
            "RADFIN4", "REGDIMS", "RKTRMDIR", "ROCKCOMP", "ROCKOPTS",
            "ROCKTAB", "RPTGRID", "RPTONLY", "RPTONLYO", "RPTPROS", "PRTRST", "RPTRUNSP",
            "RPTSCHED", "RPTSOL", "RTEMPVD", "RUNSUM", "SATOPTS", "SAVE", "SEPARATE",
            "SKIP", "SKIP100", "SKIP300", "SKIPREST", "SMRYDIMS", "SPECGRID", "SSOL",
            "SUMTHIN", "TEMP", "THCONR", "TRACER", "TRACERS", "UDADIMS", "UDQDIMS",
            "VAPPARS", "VISCREF", "WATVISCT",
            "WPAVE", "WPIMULT", "WPITAB", "WSEGDIMS", "WTEMP",
            "WTEST", "WTRACER", "ZIPPY2" };
        
        // check deck and keyword for flow and parser.
        for (size_t idx = 0; idx < deck->size(); ++idx) {
            const auto& keyword = deck->getKeyword(idx);
            std::unordered_set<std::string>::const_iterator it;
            it = unsupported_keywords.find(keyword.name());
            if (it != unsupported_keywords.end()) {
                std::string msg = "Keyword '" + keyword.name() + "' is not supported by flow.\n"
                    + "In file " + keyword.getFileName() + ", line " + std::to_string(keyword.getLineNumber()) + "\n"; 
                    OpmLog::error(msg);
            }
        }
    }
} // namespace MissingFeatures
} // namespace Opm
