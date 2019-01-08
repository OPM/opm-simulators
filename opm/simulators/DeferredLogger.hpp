/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.

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


#ifndef OPM_DEFERREDLOGGER_HEADER_INCLUDED
#define OPM_DEFERREDLOGGER_HEADER_INCLUDED

#include <opm/common/OpmLog/OpmLog.hpp>

#include <string>
#include <vector>

namespace Opm
{

    class DeferredLogger
    {
    public:
        void info(const std::string& tag, const std::string& message);
        void warning(const std::string& tag, const std::string& message);
        void error(const std::string& tag, const std::string& message);
        void problem(const std::string& tag, const std::string& message);
        void bug(const std::string& tag, const std::string& message);
        void debug(const std::string& tag, const std::string& message);
        void note(const std::string& tag, const std::string& message);

        void info(const std::string& message);
        void warning(const std::string& message);
        void error(const std::string& message);
        void problem(const std::string& message);
        void bug(const std::string& message);
        void debug(const std::string& message);
        void note(const std::string& message);

        void logMessages();

    private:
        struct Message
        {
            int64_t flag;
            std::string tag;
            std::string text;
        };
        std::vector<Message> messages_;
    };

} // namespace Opm

#endif // OPM_DEFERREDLOGGER_HEADER_INCLUDED
