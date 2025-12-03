#pragma once

#ifdef __cplusplus
extern "C" {
#endif


/**
 * @brief Copy vector.
 *
 * @param y Pointer to output vector.
 * @param x Pointer to inout vector.
 * @param n Length of vectors.
 */
void vec_copy(double *y, double const * x, int n)
{
    for(int i=0;i<n;i++) y[i]=x[i];
}

/**
 * @brief Fill vector.
 *
 * @param y Pointer to output vector.
 * @param x Fill value.
 * @param n Length of vector.
 */
void vec_fill(double *y, double x, int n)
{
    for(int i=0;i<n;i++) y[i]=x;
}


#ifdef __cplusplus
}
#endif
