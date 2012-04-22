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

#ifndef OPM_NEWWELLS_H_INCLUDED
#define OPM_NEWWELLS_H_INCLUDED


/**
 * \file
 *
 * Main OPM-Core well data structure along with functions
 * to create, populate and destroy it.
 */


#ifdef __cplusplus
extern "C" {
#endif

/** Well type indicates desired/expected well behaviour. */
enum WellType         { INJECTOR, PRODUCER };
/** Type of well control equation or inequality constraint.
 *  BHP  -> Well constrained by bottom-hole pressure target.
 *  RATE -> Well constrained by total reservoir volume flow rate.
 */
enum WellControlType  { BHP     , RATE     };
/** Canonical component names and ordering. */
enum SurfaceComponent { WATER = 0, OIL = 1, GAS = 2 };


/** Controls for a single well.
 *  Each control specifies a well rate or bottom-hole pressure. Only
 *  one control can be active at a time, indicated by current. The
 *  meaning of each control's target value depends on the control
 *  type, for BHP controls it is a pressure in Pascal, for RATE
 *  controls it is a volumetric rate in cubic(meter)/second. The
 *  active control as an equality constraint, whereas the
 *  non-active controls should be interpreted as inequality
 *  constraints (upper or lower bounds).  For instance, a PRODUCER's BHP
 *  constraint defines a minimum acceptable bottom-hole pressure value
 *  for the well.
 */
struct WellControls
{
    int                     num;     /** Number of controls. */
    enum WellControlType   *type;    /** Array of control types. */
    double                 *target;  /** Array of control targets. */
    int                     current; /** Index of current active control. */

    void                   *data;    /** Internal management structure. */
};



/** Data structure aggregating static information about all wells in a scenario. */
struct Wells
{
    int                  number_of_wells; /** Number of wells. */

    enum WellType       *type;            /** Array of well types. */
    double              *depth_ref;       /** Array of well bhp reference depths. */
    double              *zfrac;           /** Component volume fractions for each well, size is (3*number_of_wells).
                                           *  This is intended to be used for injection wells. For production wells
                                           *  the component fractions will vary and cannot be specified a priori.
                                           */
    int                 *well_connpos;    /** Array of indices into well_cells (and WI).
                                          *  For a well w, well_connpos[w] and well_connpos[w+1] yield
                                          *  start and one-beyond-end indices into the well_cells array
                                          *  for accessing w's perforation cell indices.
                                          */
    int                 *well_cells;      /** Array of perforation cell indices.
                                           *  Size is number of perforations (== well_connpos[number_of_wells]).
                                           */
    double              *WI;              /** Well productivity index, same size and structure as well_cells. */
    struct WellControls **ctrls;          /** Well controls, one set of controls for each well. */


    void               *data;             /** Internal management structure. */
};


/** Data structure aggregating dynamic information about all wells in a scenario.
 *  All arrays in this structure contain data for each perforation,
 *  ordered the same as Wells::well_cells and Wells:WI.  The array
 *  sizes are, respectively,
 *
 *     gpot       n*NP
 *     A          n²*NP (matrix in column-major (i.e., Fortran) order).
 *     phasemob   n*NP
 *
 *  in which "n" denotes the number of active fluid phases (and
 *  constituent components) and "NP" is the total number of
 *  perforations, <CODE>well_connpos[ number_of_wells ]</CODE>.
 */
struct CompletionData
{
    double *gpot;     /** Gravity potentials. */
    double *A;        /** Volumes to surface-components matrix, A = RB^{-1}. */
    double *phasemob; /** Phase mobilities. */
};

/**
 * Construct a Wells object initially capable of managing a given
 * number of wells and total number of well connections
 * (perforations).
 *
 * Function add_well() is used to populate the Wells object.  No
 * reallocation occurrs in function add_well() as long as the
 * initially indicated capacites are sufficient.  Call function
 * destroy_wells() to dispose of the Wells object and its allocated
 * memory resources.
 *
 * \param[in] nwells Expected number of wells in simulation scenario.
 *                   Pass zero if the total number of wells is unknown.
 *
 * \param[in] nperf  Expected total number of well connections
 *                   (perforations) for all wells in simulation
 *                   scenario.  Pass zero if the total number of well
 *                   connections is unknown.
 *
 * \return A valid Wells object with no wells if successful, and NULL
 * otherwise.
 */
struct Wells *
create_wells(int nwells, int nperf);


/**
 * Append a new well to an existing Wells object.
 *
 * Increments W->number_of_wells by one if successful.  The new well
 * does not include operational constraints.  Such information is
 * specified using function append_well_controls().  The current
 * control index is set to -1 (invalid).
 *
 * \param[in] type       Type of well.
 * \param[in] depth_ref  Reference depth for bhp.
 * \param[in] nperf      Number of perforations.
 * \param[in] zfrac      Injection fraction (three components) or NULL.
 * \param[in] cells      Perforation cell indices.
 * \param[in] WI         Well production index per perforation, or NULL.
 * \param[inout] W       The Wells object to be modified.
 * \return Non-zero (true) if successful and zero otherwise.
 */
int
add_well(enum WellType  type     ,
         double         depth_ref,
         int            nperf    ,
         const double  *zfrac    ,
         const int     *cells    ,
         const double  *WI       ,
         struct Wells  *W        );


/**
 * Append operational constraint to an existing well.
 *
 * Increments ctrl->num by one if successful.  Introducing a new
 * operational constraint does not affect the well's notion of the
 * currently active constraint represented by ctrl->current.
 *
 * \param[in] type     Control type.
 * \param[in] target   Target value for the control.
 * \param[inout] ctrl  Existing set of well controls.
 * \return Non-zero (true) if successful and zero (false) otherwise.
 */
int
append_well_controls(enum WellControlType type  ,
                     double               target,
                     struct WellControls *ctrl  );

/** Clear all controls from a well. */
void
clear_well_controls(struct WellControls *ctrl);


/**
 * Wells object destructor.
 *
 * Disposes of all resources managed by the Wells object.
 *
 * The Wells object must be built using function create_wells() and
 * subsequently populated using function add_well().
 */
void
destroy_wells(struct Wells *W);


#ifdef __cplusplus
}
#endif

#endif /* OPM_NEWWELLS_H_INCLUDED */
