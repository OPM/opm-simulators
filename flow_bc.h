/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#ifndef FLOW_BC_H_INCLUDED
#define FLOW_BC_H_INCLUDED

#include <stddef.h>


#ifdef __cplusplus
extern "C" {
#endif

enum flowbc_type { UNSET, PRESSURE, FLUX };

typedef struct {
    enum flowbc_type *type;
    double           *bcval;
} flowbc_t;


flowbc_t *
allocate_flowbc(size_t nf);

void
deallocate_flowbc(flowbc_t *fbc);


#ifdef __cplusplus
}
#endif

#endif  /* FLOW_BC_H_INCLUDED */
