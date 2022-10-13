/*
  Copyright 2020 Equinor ASA

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

#if HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <sys/time.h>
#include <fstream>

namespace Opm
{
namespace Accelerator
{

double second(void)
{
    struct timeval tv;
    gettimeofday(&tv, nullptr);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

bool even(int n)
{
    if (n % 2 == 0) {
        return true;
    } else {
        return false;
    }
}

int roundUpTo(int i, int n)
{
    if (i % n == 0) {
        return i;
    } else {
        return (i / n + 1) * n;
    }
}

bool fileExists(const char *filename)
{
    FILE *fin;
    fin = fopen(filename, "r");
    if (fin == nullptr) {
        return false;
    } else {
        fclose(fin);
        return true;
    }
}

} // namespace Accelerator
} // namespace Opm
