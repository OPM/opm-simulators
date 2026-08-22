#include "dot.h"

#include <stdio.h>

double vec_dot (double const *x, double const *y, int NN)
{
    // unroll loop in multiples of 8
    int n=NN/8;
    int N=8*n;
    double agg[8];
    for(int i=0;i<8;i++) agg[i]=0.0;
    for(int i=0;i<N;i+=8) for(int j=0;j<8;j++) agg[j]+=x[i+j]*y[i+j];
    for(int j=0;j<4;j++) agg[j]+=agg[j+4];
    for(int j=0;j<2;j++) agg[j]+=agg[j+2];
    for(int j=0;j<1;j++) agg[j]+=agg[j+1];

    // loop-peeling of trailing end
    for(int j=N;j<NN;j++) agg[0]+=x[j]*y[j];

    return agg[0];
}

double vec_bdot (double const *x, double const *y, int NN)
{
    // unroll loop in multiples of 8
    // assumes vectors are padded by zeros
    // to nearest multiple of 8 doubles
    int N=8*((NN+7)/8);
    double agg[8];
    for(int i=0;i<8;i++) agg[i]=0.0;
    for(int i=0;i<N;i+=8) for(int j=0;j<8;j++) agg[j]+=x[i+j]*y[i+j];

    for(int j=0;j<4;j++) agg[j]+=agg[j+4];
    for(int j=0;j<2;j++) agg[j]+=agg[j+2];
    for(int j=0;j<1;j++) agg[j]+=agg[j+1];

    return agg[0];
}
