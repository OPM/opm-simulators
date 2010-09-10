#ifndef WELL_H_INCLUDED
#define WELL_H_INCLUDED

enum well_type         { INJECTOR, PRODUCER };
enum well_control_type { BHP     , RATE     };

#define WELL_DESCRIPTOR                         \
    int                     number_of_wells;    \
                                                \
    enum well_type         *type;               \
                                                \
    enum well_control_type *control_type;       \
    double                 *target;

#define WELL_TOPOLOGY                           \
    int *well_connpos;                          \
    int *well_cells;

#define WELL_CONNECTION_DATA                    \
    double *WI;                                 \
    double *depth;

/* ---------------------------------------------------------------------- */

typedef struct {
    WELL_DESCRIPTOR
    WELL_TOPOLOGY
    WELL_CONNECTION_DATA
} well_t;

#endif
