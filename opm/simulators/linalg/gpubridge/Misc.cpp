#include <config.h>

#include <opm/simulators/linalg/gpubridge/Misc.hpp>

#include <cmath>
#include <algorithm>

namespace Opm::Accelerator {

// divide A by B, and round up: return (int)ceil(A/B)
unsigned int ceilDivision(const unsigned int A,
                          const unsigned int B)
{
    return A / B + (A % B > 0);
}

// return the absolute value of the N elements for which the absolute value is highest
template<class Scalar>
Scalar get_absmax(const Scalar* data,
                  const int N)
{
    return std::abs(*std::max_element(data, data + N,
                                      [](Scalar a, Scalar b)
                                      { return std::fabs(a) < std::fabs(b); }));
}

// solve A^T * x = b
template<class Scalar>
void solve_transposed_3x3(const Scalar* A,
                          const Scalar* b,
                          Scalar* x)
{
    const int B = 3;
    // from dune-common/densematrix.hh, but transposed, so replace [r*B+c] with [r+c*B]
    Scalar t4  = A[0+0*B] * A[1+1*B];
    Scalar t6  = A[0+0*B] * A[1+2*B];
    Scalar t8  = A[0+1*B] * A[1+0*B];
    Scalar t10 = A[0+2*B] * A[1+0*B];
    Scalar t12 = A[0+1*B] * A[2+0*B];
    Scalar t14 = A[0+2*B] * A[2+0*B];

    Scalar d = (t4*A[2+2*B]-t6*A[2+1*B]-t8*A[2+2*B]+
               t10*A[2+1*B]+t12*A[1+2*B]-t14*A[1+1*B]); // determinant

    x[0] = (b[0]*A[1+1*B]*A[2+2*B] - b[0]*A[2+1*B]*A[1+2*B]
          - b[1] *A[0+1*B]*A[2+2*B] + b[1]*A[2+1*B]*A[0+2*B]
          + b[2] *A[0+1*B]*A[1+2*B] - b[2]*A[1+1*B]*A[0+2*B]) / d;

    x[1] = (A[0+0*B]*b[1]*A[2+2*B] - A[0+0*B]*b[2]*A[1+2*B]
          - A[1+0*B] *b[0]*A[2+2*B] + A[1+0*B]*b[2]*A[0+2*B]
          + A[2+0*B] *b[0]*A[1+2*B] - A[2+0*B]*b[1]*A[0+2*B]) / d;

    x[2] = (A[0+0*B]*A[1+1*B]*b[2] - A[0+0*B]*A[2+1*B]*b[1]
          - A[1+0*B] *A[0+1*B]*b[2] + A[1+0*B]*A[2+1*B]*b[0]
          + A[2+0*B] *A[0+1*B]*b[1] - A[2+0*B]*A[1+1*B]*b[0]) / d;
}

#define INSTANTIATE_TYPE(T) \
    template void solve_transposed_3x3<T>(const T* A, const T* b, T* x); \
    template T get_absmax<T>(const T* data, const int N);

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

}
