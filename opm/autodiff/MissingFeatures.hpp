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

#ifndef OPM_MISSINGFEATURES_HEADER_INCLUDED
#define OPM_MISSINGFEATURES_HEADER_INCLUDED


namespace Opm {

class ErrorGuard;
class ParseContext;

namespace MissingFeatures {

    template <typename T>
    struct PartiallySupported {
        std::string item;
        T item_value;
    };

    template <typename Keyword, typename Item, typename T>
    void addSupported(std::multimap<std::string, PartiallySupported<T> >& map, T itemValue);

    template <typename T>
    void checkOptions(const DeckKeyword& keyword, std::multimap<std::string , PartiallySupported<T> >& map, const ParseContext& parseContext, ErrorGuard& errorGuard);

    void checkKeywords(const Deck& deck, const ParseContext& parseContext, ErrorGuard& errorGuard);

    template<typename T>
    void checkKeywords(const Deck& deck, const ParseContext& parseContext, T&& errorGuard);

    void checkKeywords(const Deck& deck);
}

}


#endif // OPM_MISSINGFEATURES_HEADER_INCLUDED
