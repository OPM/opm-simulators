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

#ifndef OPM_FLOW_BC_HEADER_INCLUDED
#define OPM_FLOW_BC_HEADER_INCLUDED

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

enum FlowBCType { BC_NOFLOW      ,
                  BC_PRESSURE    ,
                  BC_FLUX_TOTVOL };

/* Boundary condition structure.
 *
 * Condition i (in [0 .. nbc-1]) affects (outer) interface face[i], is
 * of type type[i], and specifies a target value of value[i].
 *
 * The field 'cpty' is for internal use by the implementation. */
struct FlowBoundaryConditions {
    size_t           nbc;       /* Current number of bdry. conditions */

    size_t           face_cpty; /* Internal management.  Do not touch */
    size_t           cond_cpty; /* Internal management.  Do not touch */

    size_t          *cond_pos;  /* Indirection pointer into '.face' */

    enum FlowBCType *type;      /* Condition type */
    double          *value;     /* Condition value (target) */
    int             *face;      /* Outer faces affected by ind. target */
};


/* Allocate a 'FlowBoundaryConditions' structure, initially capable of
 * managing 'nbc' individual boundary conditions.  */
struct FlowBoundaryConditions *
flow_conditions_construct(size_t nbc);


/* Release memory resources managed by 'fbc', including the containing
 * 'struct' pointer, 'fbc'. */
void
flow_conditions_destroy(struct FlowBoundaryConditions *fbc);


/* Append a new boundary condition to existing set.
 *
 * Return one (1) if successful, and zero (0) otherwise. */
int
flow_conditions_append(enum FlowBCType                type ,
                       int                            face ,
                       double                         value,
                       struct FlowBoundaryConditions *fbc  );

/* Append a new boundary condition that affects multiple interfaces.
 *
 * Return one (1) if successful, and zero (0) otherwise. */
int
flow_conditions_append_multi(enum FlowBCType                type  ,
                             size_t                         nfaces,
                             const int                     *faces ,
                             double                         value ,
                             struct FlowBoundaryConditions *fbc   );


/* Clear existing set of boundary conditions */
void
flow_conditions_clear(struct FlowBoundaryConditions *fbc);

#ifdef __cplusplus
}
#endif

#endif  /* OPM_FLOW_BC_HEADER_INCLUDED */
