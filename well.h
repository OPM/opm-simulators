#ifndef WELL_H_INCLUDED
#define WELL_H_INCLUDED

enum well_type    { INJECTOR, PRODUCER };
enum well_control { BHP     , RATE     };

typedef struct {
    int  number_of_wells;
    int *well_connpos;
    int *well_cells;
} well_t;

typedef struct {
    enum well_type    *type;
    enum well_control *ctrl;
    double            *target;
} well_control_t;

int
allocate_cell_wells(int nc, well_t *W, int **cwpos, int **cwells);

void
deallocate_cell_wells(int *cvpos, int *cwells);

void
derive_cell_wells(int nc, well_t *W, int *cwpos, int *cwells);

#endif
