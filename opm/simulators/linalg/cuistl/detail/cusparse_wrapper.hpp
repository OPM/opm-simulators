/*
  Copyright 2022-2023 SINTEF AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * Contains wrappers to make the CuSPARSE library behave as a modern C++ library with function overlading.
 *
 * In simple terms, this allows one to call say cusparseBsrilu02_analysis on both double and single precisision,
 * instead of calling cusparseDbsrilu02_analysis and cusparseDbsrilu02_analysis respectively.
 */
#include <cusparse.h>
#include <type_traits>
#ifndef OPM_CUSPARSE_WRAPPER_HPP
#define OPM_CUSPARSE_WRAPPER_HPP
namespace Opm::cuistl::detail
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
} // namespace Opm::cuistl::detail
#endif
