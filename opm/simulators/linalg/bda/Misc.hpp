#ifndef OPM_MISC_HPP
#define OPM_MISC_HPP

namespace Opm::Accelerator {

unsigned int ceilDivision(const unsigned int A, const unsigned int B);

template<class Scalar>
Scalar get_absmax(const Scalar *data, const int N);

template<class Scalar>
void solve_transposed_3x3(const Scalar *A, const Scalar *b, Scalar *x);

}

#endif
