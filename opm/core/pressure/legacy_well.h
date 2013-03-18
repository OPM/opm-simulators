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

#ifndef OPM_LEGACY_WELL_HEADER_INCLUDED
#define OPM_LEGACY_WELL_HEADER_INCLUDED

/**
 * \file
 * Deprecated (and obsolescent) well definition.  Still in use by
 * the hybridized pressure solvers.
 */

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Well taxonomy.
 */
enum well_type    { INJECTOR, PRODUCER };

/**
 * Control types recognised in system.
 */
enum well_control { BHP     , RATE     };

/**
 * Compositions recognised in injection wells.
 */
enum surface_component { WATER = 0, OIL = 1, GAS = 2 };

/**
 * Basic representation of well topology.
 */
struct WellCompletions {
    int  number_of_wells; /**< Number of wells. */
    int *well_connpos;    /**< Well topology start pointers. */
    int *well_cells;      /**< Well connections */
};

/**
 * Basic representation of well controls.
 */
struct WellControls {
    enum well_type    *type;    /**< Individual well taxonomy */
    enum well_control *ctrl;    /**< Individual well controls */
    double            *target;  /**< Control target */
    double            *zfrac;   /**< Surface injection composition */
};

/**
 * Dynamic discretisation data relating well to flow in reservoir.
 */
struct completion_data {
    double *WI;       /**< Well indices */
    double *gpot;     /**< Gravity potential */
    double *A;        /**< \f$RB^{-1}\f$ for compressible flows. */
    double *phasemob; /**< Phase mobility, per connection. */
};

/**
 * Convenience type alias to preserve backwards compatibility in
 * well topology definitions used by hybridised pressure solver.
 */
typedef struct WellCompletions well_t;

/**
 * Convenience type alias to preserve backwards compatiblity in
 * well control definitions used by hybridised pressure solver.
 */
typedef struct WellControls    well_control_t;

/**
 * Allocate cell-to-well mapping (as a sparse array).
 *
 * @param[in]  nc      Total number of cells.
 * @param[in]  W       Well topology (well-to-cell mapping).
 * @param[out] cwpos   Indirection array.  Points to array of size
 *                     <CODE>nc + 1</CODE> if successful.
 * @param[out] cwells  Cell-to-well mapping.  Points to array
 *                     of size <CODE>W->well_connpos[
 *                     W->number_of_wells]</CODE> if successful.
 * @return Positive number (size of <CODE>*cwells</CODE>)
 * if successful. Zero in case of allocation failure.
 */
int
allocate_cell_wells(int nc, well_t *W, int **cwpos, int **cwells);

/**
 * Dispose of memory resources allocated using function
 * allocate_cell_wells().
 *
 * Following a call to deallocate_cell_wells(), the input pointers
 * are no longer valid.
 *
 * @param[in,out] cvpos  Cell-to-well start pointers.
 * @param[in,out] cwells Cell-to-well mapping.
 */
void
deallocate_cell_wells(int *cvpos, int *cwells);

/**
 * Construct cell-to-well mapping (i.e., transpose the
 * well-to-cell mapping represented by <CODE>W->well_cells</CODE>).
 *
 * @param[in]  nc      Total number of cells.
 * @param[in]  W       Well topology (well-to-cell mapping).
 * @param[out] cwpos   Cell-to-well start pointers.
 * @param[out] cwells  Cell-to-well mapping.
 */
void
derive_cell_wells(int nc, well_t *W, int *cwpos, int *cwells);


#ifdef __cplusplus
}
#endif

#endif /* OPM_LEGACY_WELL_HEADER_INCLUDED */
