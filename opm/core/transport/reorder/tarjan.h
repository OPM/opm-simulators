/* Copyright 2011 (c) Jostein R. Natvig <Jostein.R.Natvig at sintef.no> */

#ifndef TARJAN_H_INCLUDED
#define TARJAN_H_INCLUDED
#ifdef __cplusplus
extern "C" {
#endif

    void tarjan (int size, int *ia, int *ja, int *rowP, int *P, 
                 int *ncomp, int *work);

#ifdef __cplusplus
}
#endif
#endif /* TARJAN_H_INCLUDED */

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
