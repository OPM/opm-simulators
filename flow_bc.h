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
