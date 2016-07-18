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
#include <map>


namespace Opm {

namespace MissingFeatures {

    template<typename T>
    struct PartiallySupported {
        int item;
        T item_value;
    };

    std::multimap<std::string, PartiallySupported<std::string> >
    string_options = { {"COMPORD", PartiallySupported<std::string>{1, "TRACK"}},
                       {"ENDSCALE",PartiallySupported<std::string>{0, "NODIR"}},
                       {"ENDSCALE",PartiallySupported<std::string>{1, "REVER"}},
                       {"PINCH",   PartiallySupported<std::string>{1, "GAP"}},
                       {"PINCH",   PartiallySupported<std::string>{3, "TOPBOT"}}
    };

    std::multimap<std::string, PartiallySupported<int> >
    int_options = { {"EHYSTR", PartiallySupported<int>{1, 0}}};

    void checkKeywords(const Deck& deck)
    {
        // These keywords are supported by opm-parser, but are not supported
        // by flow. For some of them, only part of the options are supported.
        // The list is used to output messages only.
        std::unordered_set<std::string> unsupported_keywords = {
            "ADSALNOD", "API", "AQUCON", "AQUNUM"
            "COMPLUMP", "CONNECTION", "CPR",
            "DATE", "ECHO", "EDITNNC", "ENDNUM",
            "ENDSKIP", "ENKSRVD", "ENPTVD", "EQLNUM", "EQUALREG",
            "EXCEL", "EXTRAPMS", "FILLEPS", "FIPNUM",
            "FULLIMP", "GDORIENT", "GECON", "GEFAC", "GRIDUNIT",
            "GRUPNET", "IMKRVD", "IMPES", "IMPTVD", "MAPUNITS",
            "MAXVALUE", "MESSAGES", "MINVALUE", "MONITOR", "MSGFILE",
            "MULT_XYZ", "NETBALAN", "NEXTSTEP", "NOCASC", "NOECHO",
            "NOGGF", "NOINSPEC", "NOMONITO", "NONNC", "NORSSPEC", "NOSIM",
            "NSTACK", "NUMRES", "NUPCOL", "OILVISCT", "OLDTRAN", "OPTIONS",
            "PARALLEL", "PBVD", "PCG", "PERMXY", "PERMYZ",
            "PERMZX", "PIMULTAB", "PLYADSS", "PLYDHFLF",
            "RADFIN4", "RKTRMDIR", "ROCKCOMP", "ROCKOPTS",
            "ROCKTAB", "RPTGRID", "RPTONLY", "RPTONLYO", "RPTPROS", "PRTRST", "RPTRUNSP",
            "RPTSCHED", "RPTSOL", "RTEMPVD", "RUNSUM", "SATOPTS", "SAVE", "SEPARATE",
            "SKIP", "SKIP100", "SKIP300", "SKIPREST", "SPECGRID", "SSOL",
            "SUMTHIN", "TEMP", "THCONR", "TRACER", "TRACERS",
            "VAPPARS", "VISCREF", "WATVISCT",
            "WPAVE", "WPIMULT", "WPITAB", "WTEMP",
            "WTEST", "WTRACER", "ZIPPY2" };

        // check deck and keyword for flow and parser.
        for (size_t idx = 0; idx < deck.size(); ++idx) {
            const auto& keyword = deck.getKeyword(idx);
            std::unordered_set<std::string>::const_iterator it;
            it = unsupported_keywords.find(keyword.name());
            if (it != unsupported_keywords.end()) {
                std::string msg = "Keyword '" + keyword.name() + "' is not supported by flow.\n"
                    + "In file " + keyword.getFileName() + ", line " + std::to_string(keyword.getLineNumber()) + "\n";
                OpmLog::error(msg);
            }

            // check for partially supported keywords.
            std::multimap<std::string, PartiallySupported<std::string>>::iterator string_it, string_low, string_up;
            string_low = string_options.lower_bound(keyword.name());
            string_up  = string_options.upper_bound(keyword.name());
            for (string_it = string_low; string_it != string_up; ++string_it) {
                const auto& record = keyword.getRecord(0);
                if (record.getItem(string_it->second.item).get<std::string>(0) != string_it->second.item_value) {
                    std::string msg = "For keyword '" + string_it->first + "' only value " + string_it->second.item_value + " in item " +
                        std::to_string(string_it->second.item + 1) + " is supported by flow.\n"
                        + "In file " + keyword.getFileName() + ", line " + std::to_string(keyword.getLineNumber()) + "\n";
                    OpmLog::error(msg);
                }
            }

            std::multimap<std::string, PartiallySupported<int>>::iterator int_it, int_low, int_up;
            int_low = int_options.lower_bound(keyword.name());
            int_up  = int_options.upper_bound(keyword.name());
            for (int_it = int_low; int_it != int_up; ++int_it) {
                const auto& record = keyword.getRecord(0);
                if (record.getItem(int_it->second.item).get<int>(0) != int_it->second.item_value) {
                    std::string msg = "For keyword '" + int_it->first + "' only value " + std::to_string(int_it->second.item_value) + " in item " +
                        std::to_string(int_it->second.item + 1) + " is supported by flow.\n"
                        + "In file " + keyword.getFileName() + ", line " + std::to_string(keyword.getLineNumber()) + "\n";
                    OpmLog::error(msg);
                }
            }
        }
    }
} // namespace MissingFeatures
} // namespace Opm
