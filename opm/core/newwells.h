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

#ifdef __cplusplus
extern "C" {
#endif

/** Well type indicates desired/expected well behaviour. */
enum WellType         { INJECTOR, PRODUCER };
/** Type of well control equation or inequality. 
 *  BHP  -> bottom hole pressure is specified.
 *  RATE -> flow rate is specified.
 */
enum WellControlType  { BHP     , RATE     };
/** Canonical component names and ordering. */
enum SurfaceComponent { WATER = 0, OIL = 1, GAS = 2 };


/** Controls for a single well.
 *  Each control specifies a well rate or bottom hole pressure. Only
 *  one control can be active at a time, indicated by current. The
 *  meaning of each control's target value depends on the control
 *  type, for BHP controls it is a pressure in Pascal, for RATE
 *  controls it is a volumetric rate in cubic(meter)/second. The
 *  active control should be interpreted as an equation, whereas the
 *  non-active controls should be interpreted as inequalities
 *  specifying constraints on the solution, where BHP controls yield
 *  minimum pressures, and RATE controls yield maximum rates.
 */
struct WellControls
{
    int                     num;     /** Number of controls. */
    int                     cpty;    /** Allocated capacity, for internal use only. */
    enum WellControlType   *type;    /** Array of control types. */
    double                 *target;  /** Array of control targets. */
    int                     current; /** Index of current active control. */
};



/** Struct encapsulating static information about all wells in a scenario. */
struct Wells
{
    int                  number_of_wells; /** Number of wells. */
    int                  well_cpty;       /** Allocated well capacity, for internal use only. */
    int                  perf_cpty;       /** Allocated perforation capacity, for internal use only. */
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
    struct WellControls **ctrls;          /** Well controls, one struct for each well. */
};


/** Struct encapsulating dynamic information about all wells in a scenario.
 *  All arrays in this struct contain data for each perforation, ordered
 *  the same as Wells::well_cells and Wells:WI. Letting NP be the number
 *  of perforations, the array sizes are:
 *     gpot       3*NP
 *     A          9*NP (matrix in Fortran order).
 *     phasemob   3*NP
 * \TODO: Verify that the sizes are correct, check if we should refactor to handle two phases better.
 */
struct CompletionData
{
    double *gpot;     /** Gravity potentials. */
    double *A;        /** Volumes to surface-components matrix, A = RB^{-1}. */
    double *phasemob; /** Phase mobilities. */
};

/** Contruction function initializing a Wells object.
 *  The arguments may be used to indicate expected capacity needed,
 *  they will be used internally for pre-allocation.
 *  \return NULL upon failure, otherwise a valid Wells object with 0 wells.
 *  Call add_well() to populate the Wells object.
 *  Call destroy_wells() to deallocate and clean up the Wells object.
 */
struct Wells *
create_wells(int nwells_reserve_cap, int nperf_reserve_cap);


/** Append a well to a Wells struct.
 *  If successful, W->number_of_wells is incremented by 1.
 *  The newly added well will have no controls associated with it, add
 *  controls using append_well_controls(). The current control index is set
 *  to -1 (invalid).
 *  \param[in] type       Type of well.
 *  \param[in] depth_ref  Reference depth for bhp.
 *  \param[in] nperf      Number of perforations.
 *  \param[in] zfrac      Injection fraction (three components) or NULL.
 *  \param[in] cells      Perforation cell indices.
 *  \param[in] WI         Well production index per perforation, or NULL.
 *  \param[inout] W       The Wells object to be modified.
 *  \return 1 if successful, 0 if failed.
 */
int
add_well(enum WellType  type     ,
         double         depth_ref,
         int            nperf    ,
         const double  *zfrac    ,
         const int     *cells    ,
         const double  *WI       ,
         struct Wells  *W        );


/** Append a control to a well.
 *  If successful, ctrl->num is incremented by 1.
 *  Note that this function does not change ctrl->current.
 *  To append a control to a well with index w, pass its
 *  controls to this function via wellsptr->ctrls[w].
 *  \param[in] type     Control type.
 *  \param[in] target   Target value for the control.
 *  \param[inout] ctrl  The WellControls object to be modified.
 *  \return 1 if successful, 0 if failed.
 */
int
append_well_controls(enum WellControlType type  ,
                     double               target,
                     struct WellControls *ctrl  );

/** Clear all controls from a well. */
void
clear_well_controls(struct WellControls *ctrl);


/** Destruction function for Wells objects.
 *  Assumes that create_wells() and add_wells() have been used to
 *  build the object.
 */
void
destroy_wells(struct Wells *W);


#ifdef __cplusplus
}
#endif

#endif /* OPM_NEWWELLS_H_INCLUDED */
