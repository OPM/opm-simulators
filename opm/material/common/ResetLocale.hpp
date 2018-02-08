// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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
/*!
 * \file
 * \brief Provides a free function to reset the localization settings
 *
 * Under some circumstances, some environments seem to set a locale which they do not
 * install. In turn this leads to std::runtime_errror being thrown by some parts of Boost
 * (for some versions) which causes unsolicited program aborts.
 *
 * This issue asside, it looks pretty weird if the e.g. the number format is inconsistent
 * with the language used by rest of the simulation.
 */
#ifndef OPM_RESET_LOCALE_HH
#define OPM_RESET_LOCALE_HH

#include <stdlib.h>

namespace Opm {

inline void resetLocale()
{
#ifndef WIN32
    // this probably only works for POSIX compatible operating systems. for all others,
    // unsetting a few environment variables should not hurt, though.
    unsetenv("LC_ALL");
    unsetenv("LANG");
    unsetenv("LANGUAGE");
    unsetenv("LC_ADDRESS");
    unsetenv("LC_COLLATE");
    unsetenv("LC_CTYPE");
    unsetenv("LC_IDENTIFICATION");
    unsetenv("LC_MEASUREMENT");
    unsetenv("LC_MESSAGES");
    unsetenv("LC_MONETARY");
    unsetenv("LC_NAME");
    unsetenv("LC_NUMERIC");
    unsetenv("LC_PAPER");
    unsetenv("LC_TELEPHONE");
    unsetenv("LC_TIME");
#endif // !WIN32
}

} // namespace Opm

#endif
