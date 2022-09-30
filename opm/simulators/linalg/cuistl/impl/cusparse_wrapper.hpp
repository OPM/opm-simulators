#include <cusparse.h>
#include <type_traits>
#ifndef OPM_CUSPARSE_WRAPPER_HEADER_INCLUDED
#define OPM_CUSPARSE_WRAPPER_HEADER_INCLUDED
namespace Opm::cuistl::impl
{

inline cusparseStatus_t
cusparseBsrilu02_analysis(cusparseHandle_t handle,
                          cusparseDirection_t dirA,
                          int mb,
                          int nnzb,
                          const cusparseMatDescr_t descrA,
                          double* bsrSortedVal,
                          const int* bsrSortedRowPtr,
                          const int* bsrSortedColInd,
                          int blockDim,
                          bsrilu02Info_t info,
                          cusparseSolvePolicy_t policy,
                          void* pBuffer)
{
    return cusparseDbsrilu02_analysis(handle,
                                      dirA,
                                      mb,
                                      nnzb,
                                      descrA,
                                      bsrSortedVal,
                                      bsrSortedRowPtr,
                                      bsrSortedColInd,
                                      blockDim,
                                      info,
                                      policy,
                                      pBuffer);
}

inline cusparseStatus_t
cusparseBsrsv2_analysis(cusparseHandle_t handle,
                        cusparseDirection_t dirA,
                        cusparseOperation_t transA,
                        int mb,
                        int nnzb,
                        const cusparseMatDescr_t descrA,
                        const double* bsrSortedValA,
                        const int* bsrSortedRowPtrA,
                        const int* bsrSortedColIndA,
                        int blockDim,
                        bsrsv2Info_t info,
                        cusparseSolvePolicy_t policy,
                        void* pBuffer)
{
    return cusparseDbsrsv2_analysis(handle,
                                    dirA,
                                    transA,
                                    mb,
                                    nnzb,
                                    descrA,
                                    bsrSortedValA,
                                    bsrSortedRowPtrA,
                                    bsrSortedColIndA,
                                    blockDim,
                                    info,
                                    policy,
                                    pBuffer);
}

inline cusparseStatus_t
cusparseBsrsv2_analysis(cusparseHandle_t handle,
                        cusparseDirection_t dirA,
                        cusparseOperation_t transA,
                        int mb,
                        int nnzb,
                        const cusparseMatDescr_t descrA,
                        const float* bsrSortedValA,
                        const int* bsrSortedRowPtrA,
                        const int* bsrSortedColIndA,
                        int blockDim,
                        bsrsv2Info_t info,
                        cusparseSolvePolicy_t policy,
                        void* pBuffer)
{
    return cusparseSbsrsv2_analysis(handle,
                                    dirA,
                                    transA,
                                    mb,
                                    nnzb,
                                    descrA,
                                    bsrSortedValA,
                                    bsrSortedRowPtrA,
                                    bsrSortedColIndA,
                                    blockDim,
                                    info,
                                    policy,
                                    pBuffer);
}

inline cusparseStatus_t
cusparseBsrilu02_analysis(cusparseHandle_t handle,
                          cusparseDirection_t dirA,
                          int mb,
                          int nnzb,
                          const cusparseMatDescr_t descrA,
                          float* bsrSortedVal,
                          const int* bsrSortedRowPtr,
                          const int* bsrSortedColInd,
                          int blockDim,
                          bsrilu02Info_t info,
                          cusparseSolvePolicy_t policy,
                          void* pBuffer)
{
    return cusparseSbsrilu02_analysis(handle,
                                      dirA,
                                      mb,
                                      nnzb,
                                      descrA,
                                      bsrSortedVal,
                                      bsrSortedRowPtr,
                                      bsrSortedColInd,
                                      blockDim,
                                      info,
                                      policy,
                                      pBuffer);
}

inline cusparseStatus_t
cusparseBsrsv2_solve(cusparseHandle_t handle,
                     cusparseDirection_t dirA,
                     cusparseOperation_t transA,
                     int mb,
                     int nnzb,
                     const double* alpha,
                     const cusparseMatDescr_t descrA,
                     const double* bsrSortedValA,
                     const int* bsrSortedRowPtrA,
                     const int* bsrSortedColIndA,
                     int blockDim,
                     bsrsv2Info_t info,
                     const double* f,
                     double* x,
                     cusparseSolvePolicy_t policy,
                     void* pBuffer)
{
    return cusparseDbsrsv2_solve(handle,
                                 dirA,
                                 transA,
                                 mb,
                                 nnzb,
                                 alpha,
                                 descrA,
                                 bsrSortedValA,
                                 bsrSortedRowPtrA,
                                 bsrSortedColIndA,
                                 blockDim,
                                 info,
                                 f,
                                 x,
                                 policy,
                                 pBuffer);
}


inline cusparseStatus_t
cusparseBsrsv2_solve(cusparseHandle_t handle,
                     cusparseDirection_t dirA,
                     cusparseOperation_t transA,
                     int mb,
                     int nnzb,
                     const float* alpha,
                     const cusparseMatDescr_t descrA,
                     const float* bsrSortedValA,
                     const int* bsrSortedRowPtrA,
                     const int* bsrSortedColIndA,
                     int blockDim,
                     bsrsv2Info_t info,
                     const float* f,
                     float* x,
                     cusparseSolvePolicy_t policy,
                     void* pBuffer)
{
    return cusparseSbsrsv2_solve(handle,
                                 dirA,
                                 transA,
                                 mb,
                                 nnzb,
                                 alpha,
                                 descrA,
                                 bsrSortedValA,
                                 bsrSortedRowPtrA,
                                 bsrSortedColIndA,
                                 blockDim,
                                 info,
                                 f,
                                 x,
                                 policy,
                                 pBuffer);
}


inline cusparseStatus_t
cusparseBsrilu02_bufferSize(cusparseHandle_t handle,
                            cusparseDirection_t dirA,
                            int mb,
                            int nnzb,
                            const cusparseMatDescr_t descrA,
                            double* bsrSortedVal,
                            const int* bsrSortedRowPtr,
                            const int* bsrSortedColInd,
                            int blockDim,
                            bsrilu02Info_t info,
                            int* pBufferSizeInBytes)
{
    return cusparseDbsrilu02_bufferSize(handle,
                                        dirA,
                                        mb,
                                        nnzb,
                                        descrA,
                                        bsrSortedVal,
                                        bsrSortedRowPtr,
                                        bsrSortedColInd,
                                        blockDim,
                                        info,
                                        pBufferSizeInBytes);
}


inline cusparseStatus_t
cusparseBsrilu02_bufferSize(cusparseHandle_t handle,
                            cusparseDirection_t dirA,
                            int mb,
                            int nnzb,
                            const cusparseMatDescr_t descrA,
                            float* bsrSortedVal,
                            const int* bsrSortedRowPtr,
                            const int* bsrSortedColInd,
                            int blockDim,
                            bsrilu02Info_t info,
                            int* pBufferSizeInBytes)
{
    return cusparseSbsrilu02_bufferSize(handle,
                                        dirA,
                                        mb,
                                        nnzb,
                                        descrA,
                                        bsrSortedVal,
                                        bsrSortedRowPtr,
                                        bsrSortedColInd,
                                        blockDim,
                                        info,
                                        pBufferSizeInBytes);
}

inline cusparseStatus_t
cusparseBsrsv2_bufferSize(cusparseHandle_t handle,
                          cusparseDirection_t dirA,
                          cusparseOperation_t transA,
                          int mb,
                          int nnzb,
                          const cusparseMatDescr_t descrA,
                          double* bsrSortedValA,
                          const int* bsrSortedRowPtrA,
                          const int* bsrSortedColIndA,
                          int blockDim,
                          bsrsv2Info_t info,
                          int* pBufferSizeInBytes)
{
    return cusparseDbsrsv2_bufferSize(handle,
                                      dirA,
                                      transA,
                                      mb,
                                      nnzb,
                                      descrA,
                                      bsrSortedValA,
                                      bsrSortedRowPtrA,
                                      bsrSortedColIndA,
                                      blockDim,
                                      info,
                                      pBufferSizeInBytes);
}
inline cusparseStatus_t
cusparseBsrsv2_bufferSize(cusparseHandle_t handle,
                          cusparseDirection_t dirA,
                          cusparseOperation_t transA,
                          int mb,
                          int nnzb,
                          const cusparseMatDescr_t descrA,
                          float* bsrSortedValA,
                          const int* bsrSortedRowPtrA,
                          const int* bsrSortedColIndA,
                          int blockDim,
                          bsrsv2Info_t info,
                          int* pBufferSizeInBytes)
{
    return cusparseSbsrsv2_bufferSize(handle,
                                      dirA,
                                      transA,
                                      mb,
                                      nnzb,
                                      descrA,
                                      bsrSortedValA,
                                      bsrSortedRowPtrA,
                                      bsrSortedColIndA,
                                      blockDim,
                                      info,
                                      pBufferSizeInBytes);
}

inline cusparseStatus_t
cusparseBsrilu02(cusparseHandle_t handle,
                 cusparseDirection_t dirA,
                 int mb,
                 int nnzb,
                 const cusparseMatDescr_t descrA,
                 double* bsrSortedVal,
                 const int* bsrSortedRowPtr,
                 const int* bsrSortedColInd,
                 int blockDim,
                 bsrilu02Info_t info,
                 cusparseSolvePolicy_t policy,
                 void* pBuffer)
{
    return cusparseDbsrilu02(handle,
                             dirA,
                             mb,
                             nnzb,
                             descrA,
                             bsrSortedVal,
                             bsrSortedRowPtr,
                             bsrSortedColInd,
                             blockDim,
                             info,
                             policy,
                             pBuffer);
}
inline cusparseStatus_t
cusparseBsrilu02(cusparseHandle_t handle,
                 cusparseDirection_t dirA,
                 int mb,
                 int nnzb,
                 const cusparseMatDescr_t descrA,
                 float* bsrSortedVal,
                 const int* bsrSortedRowPtr,
                 const int* bsrSortedColInd,
                 int blockDim,
                 bsrilu02Info_t info,
                 cusparseSolvePolicy_t policy,
                 void* pBuffer)
{
    return cusparseSbsrilu02(handle,
                             dirA,
                             mb,
                             nnzb,
                             descrA,
                             bsrSortedVal,
                             bsrSortedRowPtr,
                             bsrSortedColInd,
                             blockDim,
                             info,
                             policy,
                             pBuffer);
}

inline cusparseStatus_t
cusparseBsrmv(cusparseHandle_t handle,
              cusparseDirection_t dirA,
              cusparseOperation_t transA,
              int mb,
              int nb,
              int nnzb,
              const double* alpha,
              const cusparseMatDescr_t descrA,
              const double* bsrSortedValA,
              const int* bsrSortedRowPtrA,
              const int* bsrSortedColIndA,
              int blockDim,
              const double* x,
              const double* beta,
              double* y)
{
    return cusparseDbsrmv(handle,
                          dirA,
                          transA,
                          mb,
                          nb,
                          nnzb,
                          alpha,
                          descrA,
                          bsrSortedValA,
                          bsrSortedRowPtrA,
                          bsrSortedColIndA,
                          blockDim,
                          x,
                          beta,
                          y);
}

inline cusparseStatus_t
cusparseBsrmv(cusparseHandle_t handle,
              cusparseDirection_t dirA,
              cusparseOperation_t transA,
              int mb,
              int nb,
              int nnzb,
              const float* alpha,
              const cusparseMatDescr_t descrA,
              const float* bsrSortedValA,
              const int* bsrSortedRowPtrA,
              const int* bsrSortedColIndA,
              int blockDim,
              const float* x,
              const float* beta,
              float* y)
{
    return cusparseSbsrmv(handle,
                          dirA,
                          transA,
                          mb,
                          nb,
                          nnzb,
                          alpha,
                          descrA,
                          bsrSortedValA,
                          bsrSortedRowPtrA,
                          bsrSortedColIndA,
                          blockDim,
                          x,
                          beta,
                          y);
}
} // namespace Opm::cuistl::impl
#endif