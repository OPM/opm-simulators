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

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <string>
#include <vector>

namespace Opm
{
    /** This class implements a deferred logger:
     * 1) messages can be pushed back to a vector
     * 2) a call to logMessages adds the messages to OpmLog backends
     * */

namespace ExceptionType
{


/*
  If an exception has been raised on more than processor simultaneously the
  highest number type will be thrown in the subsequent global rethrow.
*/

enum ExcEnum {
    NONE = 0,
    RUNTIME_ERROR = 1,
    INVALID_ARGUMENT = 2,
    LOGIC_ERROR = 3,
    DEFAULT = 4,   // will throw std::logic_error()
    NUMERICAL_ISSUE = 5
};
}


    class DeferredLogger
    {
    public:

        struct Message
        {
            int64_t flag;
            std::string tag;
            std::string text;
        };

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

        /// Log all messages to the OpmLog backends,
        /// and clear the message container.
        void logMessages();

        /// Clear the message container without logging them.
        void clearMessages();

    private:
        std::vector<Message> messages_;
        friend DeferredLogger gatherDeferredLogger(const DeferredLogger& local_deferredlogger,
                                                   Parallel::Communication mpi_communicator);
    };

} // namespace Opm

#endif // OPM_DEFERREDLOGGER_HEADER_INCLUDED
