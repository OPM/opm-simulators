#pragma once

#ifdef __cplusplus
extern "C" {
#endif


void vec_copy(double *y, double const * x, int n)
{
    for(int i=0;i<n;i++) y[i]=x[i];
}

void vec_fill(double *y, double x, int n)
{
    for(int i=0;i<n;i++) y[i]=x;
}


#ifdef __cplusplus
}
#endif
