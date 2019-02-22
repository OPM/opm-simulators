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
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/Parser/ErrorGuard.hpp>
#include <opm/parser/eclipse/Parser/ParserKeywords/C.hpp>
#include <opm/parser/eclipse/Parser/ParserKeywords/E.hpp>
#include <opm/parser/eclipse/Parser/ParserKeywords/P.hpp>

#include <opm/autodiff/MissingFeatures.hpp>

#include <unordered_set>
#include <string>
#include <map>
#include <boost/lexical_cast.hpp>

namespace Opm {

namespace MissingFeatures {


    template <typename Keyword, typename Item, typename T>
    void addSupported(std::multimap<std::string, PartiallySupported<T> >& map, T itemValue)
    {
        std::pair<std::string,PartiallySupported<T> > pair({Keyword::keywordName, PartiallySupported<T>{Item::itemName , itemValue}});
        map.insert(pair);
    }


    template <typename T>
    void checkOptions(const DeckKeyword& keyword, std::multimap<std::string , PartiallySupported<T> >& map, const ParseContext& parseContext, ErrorGuard& errorGuard)
    {
        // check for partially supported keywords.
        typename std::multimap<std::string, PartiallySupported<T> >::iterator it, itlow, itup;
        itlow = map.lower_bound(keyword.name());
        itup  = map.upper_bound(keyword.name());
        for (it = itlow; it != itup; ++it) {
            const auto& record = keyword.getRecord(0);
            if (record.getItem(it->second.item).template get<T>(0) != it->second.item_value) {
                std::string msg = "For keyword '" + it->first + "' only value " + boost::lexical_cast<std::string>(it->second.item_value)
                    + " in item " + it->second.item + " is supported by flow.\n"
                    + "In file " + keyword.getFileName() + ", line " + std::to_string(keyword.getLineNumber()) + "\n";
                parseContext.handleError(ParseContext::SIMULATOR_KEYWORD_ITEM_NOT_SUPPORTED, msg, errorGuard);
            }
        }
    }

    template <typename T>
    void checkKeywords(const Deck& deck, const ParseContext& parseContext, T&& errorGuard) {
        checkKeywords(deck, parseContext, errorGuard);
    }

    void checkKeywords(const Deck& deck) {
        checkKeywords(deck, ParseContext(), ErrorGuard());
    }

    void checkKeywords(const Deck& deck, const ParseContext& parseContext, ErrorGuard& errorGuard)
    {
        // These keywords are supported by opm-parser, but are not supported
        // by flow. For some of them, only part of the options are supported.
        // The list is used to output messages only.
        std::unordered_set<std::string> unsupported_keywords = {
            "ACTION",
            "ACTIONX",
            "ADSALNOD",
            "API",
            "APIGROUP",
            "AQUCON",
            "AQUCT",
            "AQUTAB",
            "AQUNUM"
            "CARFIN"
            "COMPDATL",
            "COMPLUMP",
            "CONNECTION",
            "CPR",
            "DATE",
            "ECHO",
            "EDITNNC",
            "ENDACTIO",
            "ENDFIN"
            "ENDNUM",
            "ENDSKIP",
            "ENKSRVD",
            "ENPTVD",
            "EQLNUM",
            "EQUALREG",
            "EXCEL",
            "EXTRAPMS",
            "FILLEPS",
            "FIPNUM",
            "FULLIMP",
            "GDORIENT",
            "GECON",
            "GLIFTOPT",
            "GNETINJE",
            "GRIDUNIT",
            "GRUPNET",
            "GSATPROD",
            "IMKRVD",
            "IMPES",
            "IMPTVD",
            "LGR",
            "LIFTOPT",
            "MAPUNITS",
            "MAXVALUE",
            "MESSAGES",
            "MINVALUE",
            "MONITOR",
            "MSGFILE",
            "MULT_XYZ",
            "NETBALAN",
            "NEXTSTEP",
            "NOCASC",
            "NOECHO",
            "NOGGF",
            "NOINSPEC",
            "NOMONITO",
            "NONNC",
            "NORSSPEC",
            "NOWARN",
            "NSTACK",
            "NUMRES",
            "NUPCOL",
            "OILVISCT",
            "OLDTRAN",
            "OPERATER",
            "OPTIONS",
            "PARALLEL",
            "PBVD",
            "PCG",
            "PERMR",
            "PERMTHT",
            "PERMXY",
            "PERMYZ",
            "PERMZX",
            "PIMULTAB",
            "PLYADSS",
            "PLYDHFLF",
            "PPCWMAX",
            "REFINE",
            "RADFIN4",
            "RHO",
            "RKTRMDIR",
            "ROCKCOMP",
            "ROCKOPTS",
            "ROCKTAB",
            "RPTGRID",
            "RPTONLY",
            "RPTONLYO",
            "RPTPROS",
            "PRTRST",
            "RPTRUNSP",
            "RPTSCHED",
            "RPTSMRY",
            "RPTSOL",
            "RSCONST",
            "RSCONSTT",
            "RTEMP",
            "RTEMPA",
            "RTEMPVD",
            "RUNSUM",
            "SATOPTS",
            "SAVE",
            "SEPARATE",
            "SKIP",
            "SKIP100",
            "SKIP300",
            "SKIPREST",
            "SPECGRID",
            "SUMTHIN",
            "TEMP",
            "THCONR",
            "TRACER",
            "TRACERS",
            "VAPPARS",
            "VISCREF",
            "WARN",
            "WATVISCT",
            "WELPI",
            "WELSPECL",
            "WGASPROD",
            "WINJMULT",
            "WLIMTOL",
            "WORKTHP",
            "WPAVE",
            "WPITAB",
            "WTEMP",
            "WTEST",
            "WTRACER",
            "ZIPPY2" };
        std::multimap<std::string, PartiallySupported<std::string> > string_options;
        std::multimap<std::string, PartiallySupported<int> > int_options;
        addSupported<ParserKeywords::COMPORD, ParserKeywords::COMPORD::ORDER_TYPE, std::string>(string_options , "DEPTH");
        addSupported<ParserKeywords::ENDSCALE, ParserKeywords::ENDSCALE::DIRECT, std::string>(string_options, "NODIR");
        addSupported<ParserKeywords::ENDSCALE, ParserKeywords::ENDSCALE::IRREVERS, std::string>(string_options, "REVER");
        addSupported<ParserKeywords::PINCH, ParserKeywords::PINCH::CONTROL_OPTION, std::string>(string_options, "GAP");
        addSupported<ParserKeywords::PINCH, ParserKeywords::PINCH::PINCHOUT_OPTION, std::string>(string_options, "TOPBOT");
        addSupported<ParserKeywords::EHYSTR, ParserKeywords::EHYSTR::relative_perm_hyst, int>(int_options , 0);

        // check deck and keyword for flow and parser.
        for (size_t idx = 0; idx < deck.size(); ++idx) {
            const auto& keyword = deck.getKeyword(idx);
            std::unordered_set<std::string>::const_iterator it;
            it = unsupported_keywords.find(keyword.name());
            if (it != unsupported_keywords.end()) {
                std::string msg = "Keyword '" + keyword.name() + "' is not supported by flow.\n"
                    + "In file " + keyword.getFileName() + ", line " + std::to_string(keyword.getLineNumber()) + "\n";
                parseContext.handleError(ParseContext::SIMULATOR_KEYWORD_NOT_SUPPORTED, msg, errorGuard);
            }
            checkOptions<std::string>(keyword, string_options, parseContext, errorGuard);
            checkOptions<int>(keyword, int_options, parseContext, errorGuard);
        }
    }
} // namespace MissingFeatures
} // namespace Opm
