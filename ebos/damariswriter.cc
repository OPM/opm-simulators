// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2022 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2023 Inria, Bretagne–Atlantique Research Center

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
#include <ebos/damariswriter.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <Damaris.h>
#include <fmt/format.h>

namespace Opm::DamarisOutput {

int setPosition(const char* field, int rank, int64_t pos)
{
    int dam_err = damaris_set_position(field, &pos);
    if (dam_err != DAMARIS_OK) {
        OpmLog::error(fmt::format("damariswriter::writeOutput()       : ( rank:{})"
                                  "damaris_set_position({}, ...), Damaris Error: {}  ",
                                  rank, field, damaris_error_string(dam_err)));
    }

    return dam_err;
}

int setParameter(const char* field, int rank, int value)
{
    int dam_err = damaris_parameter_set(field, &value, sizeof(int));
    if (dam_err != DAMARIS_OK) {
        OpmLog::error(fmt::format("(rank:{}) Damaris library produced an error result for "
                                  "damaris_parameter_set(\"{}\",...)", rank, field));
    }

    return dam_err;
}

int write(const char* field, int rank, const void* data)
{
    int dam_err = damaris_write(field, data);
    if (dam_err != DAMARIS_OK) {
        OpmLog::error(fmt::format("damariswriter::writeOutput()       : ( rank:{}) "
                                  "damaris_write({}, ...), Damaris Error: {}  ",
                                  rank, field, damaris_error_string(dam_err)));
    }

    return dam_err;
}

int endIteration(int rank)
{
    int dam_err =  damaris_end_iteration();
    if (dam_err != DAMARIS_OK) {
        OpmLog::error(fmt::format("damariswriter::writeOutput()       : ( rank:{}) "
                                  "damaris_end_iteration(), Damaris Error: {}  ",
                                  rank, damaris_error_string(dam_err)));
    }

    return dam_err;
}

}
