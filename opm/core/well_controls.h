/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_WELL_CONTROLS_H_INCLUDED
#define OPM_WELL_CONTROLS_H_INCLUDED

#include <stdbool.h>
#include <opm/core/well_control_type.h>

#ifdef __cplusplus
extern "C" {
#endif




/**
 * Controls for a single well.
 * Each control specifies a well rate or bottom-hole pressure. Only
 * one control can be active at a time, indicated by current. The
 * meaning of each control's target value depends on the control type:
 *
 *  - BHP            -> target pressure in Pascal.
 *  - RESERVOIR_RATE -> target reservoir volume rate in cubic(meter)/second
 *  - SURFACE_RATE   -> target surface volume rate in cubic(meter)/second
 *
 * The sign convention for RATE targets is as follows:
 *
 *  - (+) Fluid flowing into reservoir, i.e. injecting.
 *  - (-) Fluid flowing out of reservoir, i.e. producing.
 *
 * For *_RATE controls, the distribution of phases used for the control
 * is also needed. For example, a total rate control should have 1.0
 * for each phase, whereas a control on oil rate should have 1.0 for
 * the oil phase and 0.0 for the rest. For BHP controls, this is unused.
 * The active control acts as an equality constraint, whereas the
 * non-active controls should be interpreted as inequality
 * constraints (upper or lower bounds).  For instance, a PRODUCER's
 * BHP constraint defines a minimum acceptable bottom-hole pressure
 * value for the well.
 */

//#ifdef HAVE_WELLCONTROLS
struct WellControls
{
    /**
     * Number of controls.
     */
    int num;

    int number_of_phases;

    /**
     * Array of control types.
     */
    enum WellControlType *type;

    /**
     * Array of control targets.
     */
    double *target;

    /**
     * Array of rate control distributions,
     * <CODE>number_of_phases</CODE> numbers for each control
     */
    double *distr;

    /**
     * Index of current active control.
     */
    int current;

    /* 
       The capacity allocated.
    */
    int cpty;
};
//#else
//struct WellControls;
//#endif

bool 
well_controls_equal(const struct WellControls *ctrls1, const struct WellControls *ctrls2);

int 
well_controls_reserve(int nctrl, int nphases, struct WellControls *ctrl);

struct WellControls * 
well_controls_create(void);

void
well_controls_destroy(struct WellControls *ctrl);





#ifdef __cplusplus
}
#endif

#endif /* OPM_WELL_CONTROLS_H_INCLUDED */
