// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>
#include <opm/models/utils/terminal.hpp>

#if HAVE_MPI
#include <mpi.h>
#endif

#include <csignal>
#include <iostream>
#include <string.h>             // strsignal()
#include <sys/ioctl.h>
#include <unistd.h>

namespace Opm {

std::string breakLines(const std::string& msg,
                       int indentWidth,
                       int maxWidth)
{
    std::string result;
    int startInPos = 0;
    int inPos = 0;
    int lastBreakPos = 0;
    int ttyPos = 0;
    for (; inPos < int(msg.size()); ++ inPos, ++ ttyPos) {
        if (msg[inPos] == '\n') {
            result += msg.substr(startInPos, inPos - startInPos + 1);
            startInPos = inPos + 1;
            lastBreakPos = startInPos + 1;

            // we need to use -1 here because ttyPos is incremented after the loop body
            ttyPos = -1;
            continue;
        }

        if (std::isspace(msg[inPos])) {
            lastBreakPos = inPos;
        }

        if (ttyPos >= maxWidth) {
            if (lastBreakPos > startInPos) {
                result += msg.substr(startInPos, lastBreakPos - startInPos);
                startInPos = lastBreakPos + 1;
                lastBreakPos = startInPos;
                inPos = startInPos;
            }
            else {
                result += msg.substr(startInPos, inPos - startInPos);
                startInPos = inPos;
                lastBreakPos = startInPos;
                inPos = startInPos;
            }

            result += "\n";
            for (int i = 0; i < indentWidth; ++i) {
                result += " ";
            }
            ttyPos = indentWidth;
        }
    }

    result += msg.substr(startInPos);

    return result;
}

int getTtyWidth()
{
    int ttyWidth = 10*1000; // effectively do not break lines at all.
    if (isatty(STDOUT_FILENO)) {
#if defined TIOCGWINSZ
        // This is a bit too linux specific, IMO. let's do it anyway
        struct winsize ttySize;
        ioctl(STDOUT_FILENO, TIOCGWINSZ, &ttySize);
        ttyWidth = std::max<int>(80, ttySize.ws_col);
#else
        // default for systems that do not implement the TIOCGWINSZ ioctl
        ttyWidth = 100;
#endif
    }

    return ttyWidth;
}

void assignResetTerminalSignalHandlers()
{
    // set the signal handlers to reset the TTY to a well defined state on unexpected
    // program aborts
    if (isatty(STDIN_FILENO)) {
        signal(SIGINT, resetTerminal);
        signal(SIGHUP, resetTerminal);
        signal(SIGABRT, resetTerminal);
        signal(SIGFPE, resetTerminal);
        signal(SIGSEGV, resetTerminal);
        signal(SIGPIPE, resetTerminal);
        signal(SIGTERM, resetTerminal);
    }
}

void resetTerminal()
{
    // make sure stderr and stderr do not contain any unwritten data and make sure that
    // the TTY does not see any unfinished ANSI escape sequence.
    std::cerr << "    \r\n";
    std::cerr.flush();
    std::cout << "    \r\n";
    std::cout.flush();

    // it seems like some terminals sometimes takes their time to react, so let's
    // accommodate them.
    usleep(/*usec=*/500*1000);

    // this requires the 'stty' command to be available in the command search path. on
    // most linux systems, is the case. (but even if the system() function fails, the
    // worst thing which can happen is that the TTY stays potentially choked up...)
    if (system("stty sane") != 0) {
        std::cout << "Executing the 'stty' command failed."
                  << " Terminal might be left in an undefined state!\n";
    }
}

void resetTerminal(int signum)
{
    // first thing to do when a nuke hits: restore the default signal handler
    signal(signum, SIG_DFL);

#if HAVE_MPI
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) {
        // re-raise the signal
        raise(signum);

        return;
    }
#endif

    if (isatty(fileno(stdout)) && isatty(fileno(stdin))) {
        std::cout << "\n\nReceived signal " << signum
                  << " (\"" << strsignal(signum) << "\")."
                  << " Trying to reset the terminal.\n";

        resetTerminal();
    }

    // after we did our best to clean the pedestrian way, re-raise the signal
    raise(signum);
}

} // namespace Opm
