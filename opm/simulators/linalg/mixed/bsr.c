#include "bsr.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <immintrin.h>

#pragma GCC push_options
#pragma GCC target("avx2")

void bsr_hello()
{
    printf("HaugenLabs\n");
    getchar();
}

bsr_matrix* bsr_new()
{
    bsr_matrix *A=malloc(sizeof(bsr_matrix));
    A->nrows = 0;
    A->ncols = 0;
    A->nnz   = 0;
    A->b     = 0;

    A->rowptr = NULL;
    A->colidx = NULL;
    A->dbl    = NULL;
    A->flt    = NULL;

    return A;
}

void bsr_init(bsr_matrix *A, int nrows, int nnz, int b)
{
    A->nrows=nrows;
    A->ncols=nrows;
    A->nnz=nnz;
    A->b=b;

    A->rowptr = malloc((nrows+1)*sizeof(int));
    A->colidx = malloc(nnz*sizeof(int));
    A->dbl    = malloc(b*b*nnz*sizeof(double));
    A->flt    = malloc(b*b*nnz*sizeof(float));
}

void bsr_info(bsr_matrix *A)
{
    printf("nrows = %d\n",A->nrows);
    printf("ncols = %d\n",A->ncols);
    printf("nnz   = %d\n",A->nnz);
    printf("b     = %d\n",A->b);

    printf("rowptr= 0x%08lX\n",(uint64_t)A->rowptr);
    printf("colidx= 0x%08lX\n",(uint64_t)A->colidx);
    printf("dbl   = 0x%08lX\n",(uint64_t)A->dbl);

    printf("\n");
}

void bsr_vmspmv3(bsr_matrix *A, const double *x, double *y)
{
    int nrows = A->nrows;
    int *rowptr=A->rowptr;
    int *colidx=A->colidx;
    const float *data=A->flt;

    const int b=3;

    __m256d mm_zeros =_mm256_setzero_pd();
    for(int i=0;i<nrows;i++)
    {
        __m256d vA[3];
        for(int k=0;k<3;k++) vA[k] = mm_zeros;
        for(int k=rowptr[i];k<rowptr[i+1];k++)
        {
            const float *AA=data+9*k;

            int j = colidx[k];
            __m256d vx = _mm256_loadu_pd(x+b*j);

            vA[0] += _mm256_cvtps_pd(_mm_loadu_ps(AA+0))*_mm256_permute4x64_pd(vx,0b00000000); //0b01010101
            vA[1] += _mm256_cvtps_pd(_mm_loadu_ps(AA+3))*_mm256_permute4x64_pd(vx,0b01010101); //0b01010101
            vA[2] += _mm256_cvtps_pd(_mm_loadu_ps(AA+6))*_mm256_permute4x64_pd(vx,0b10101010); //0b01010101
        }

        // sum over columns
        __m256d vy, vz;
        vz = vA[0] + vA[1] + vA[2];

        double *y_i = y+b*i;
        vy = _mm256_loadu_pd(y_i);       // optional blend to keep
        vz =_mm256_blend_pd(vy,vz,0x7);  // 4th element unchanged
        _mm256_storeu_pd(y_i,vz);
    }
}

void bsr_downcast(bsr_matrix *M)
{
    int nnz = M->nnz;
    int b = M->b;

    if(M->flt==NULL) posix_memalign((void**)&(M->flt),64,b*b*nnz*sizeof(float));
    for(int i=0;i<b*b*nnz;i++) M->flt[i]=M->dbl[i];
}



void bsr_sparsity(const bsr_matrix *A, const char *name)
{
    printf("%s =\n[\n",name);
    int count=1;
    int offset=0;
    for(int i=0; i<A->nrows; i++)
    {
        printf("%4d: ",offset);
        for(int j=A->rowptr[i];j<A->rowptr[i+1];j++)
        {
            printf(" %4d",A->colidx[j]);
            offset++;
        }
        printf("\n");
        count++;
        if(count>16) break;
    }
    printf("]\n");
}

void bsr_nonzeros(bsr_matrix *A, const char *name)
{
    printf("%s =\n[\n",name);
    int count=1;
    int b=A->b;
    int bb=b*b;
    for(int i=0; i<A->nrows; i++)
    {
        for(int j=A->rowptr[i];j<A->rowptr[i+1];j++)
        {
            printf("|");
            for(int m=0;m<b;m++)
            {
                for(int n=0;n<b;n++)
                {
                    printf(" %+.4e",A->dbl[j*bb + m*b + n]);
                }
                printf(" |");
            }
            printf("\n");
        }
        count++;
        if(count>6) break;
        printf("\n");
    }
    printf("]\n");
}

#pragma GCC pop_options
