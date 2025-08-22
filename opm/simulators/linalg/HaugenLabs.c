#include "HaugenLabs.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>


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

void vec_show(const double *x, int n, const char *name)
{
    printf("%s =\n[\n",name);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<3;j++) printf(" %+.4e",x[3*i+j]);
        printf("\n");
    }
    printf("\n]\n");
}
