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
#ifndef OPM_TERMINAL_HPP
#define OPM_TERMINAL_HPP

#include <string>

namespace Opm {

/*!
 * \brief Break up a string in lines suitable for terminal output.
 * \param msg String to print
 * \param indentWidth Size of indent
 * \param maxWidth Maximum with of terminal
 * \return
 */
std::string breakLines(const std::string& msg,
                       int indentWidth,
                       int maxWidth);

/*!
 * \brief Get the width of the tty we are attached to.
 * \return Width of tty
 */
int getTtyWidth();

/*!
 * \brief Resets the current TTY to a usable state if the program was aborted.
 *
 * This is intended to be called as part of a generic exception handler
 */
void resetTerminal();

/*!
 * \brief Resets the current TTY to a usable state if the program was interrupted by
 *        SIGABRT or SIGINT.
 */
void resetTerminal(int signum);

} // namespace Opm

#endif // OPM_TERMINAL_HPP
