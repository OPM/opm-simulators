#ifndef OPM_MISC_HPP
#define OPM_MISC_HPP

#ifdef HAVE_ROCSPARSE
#include <hip/hip_runtime_api.h>
#include <hip/hip_version.h>
#include <sstream>

#define HIP_CHECK(STAT)                                  \
    do {                                                 \
        const hipError_t stat = (STAT);                  \
        if(stat != hipSuccess)                           \
        {                                                \
            std::ostringstream oss;                      \
            oss << "rocsparseSolverBackend::hip ";       \
            oss << "error: " << hipGetErrorString(stat); \
            OPM_THROW(std::logic_error, oss.str());      \
        }                                                \
    } while(0)

#define ROCSPARSE_CHECK(STAT)                            \
    do {                                                 \
        const rocsparse_status stat = (STAT);            \
        if(stat != rocsparse_status_success)             \
        {                                                \
            std::ostringstream oss;                      \
            oss << "rocsparseSolverBackend::rocsparse "; \
            oss << "error: " << stat;                    \
            OPM_THROW(std::logic_error, oss.str());      \
        }                                                \
    } while(0)

#define ROCBLAS_CHECK(STAT)                              \
    do {                                                 \
        const rocblas_status stat = (STAT);              \
        if(stat != rocblas_status_success)               \
        {                                                \
            std::ostringstream oss;                      \
            oss << "rocsparseSolverBackend::rocblas ";   \
            oss << "error: " << stat;                    \
            OPM_THROW(std::logic_error, oss.str());      \
        }                                                \
    } while(0)
#endif

namespace Opm::Accelerator {

unsigned int ceilDivision(const unsigned int A,
                          const unsigned int B);

template<class Scalar>
Scalar get_absmax(const Scalar *data,
                  const int N);

template<class Scalar>
void solve_transposed_3x3(const Scalar *A,
                          const Scalar *b,
                          Scalar *x);

}

#endif
