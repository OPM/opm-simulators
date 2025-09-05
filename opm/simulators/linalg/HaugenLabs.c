#include "HaugenLabs.h"

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <immintrin.h>

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
    A->dbl    = malloc(b*b*nnz*sizeof(double));
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

bildu_prec *bildu_new()
{
    bildu_prec *P = malloc(sizeof(bildu_prec));
    P->L=bsr_new();
    P->D=bsr_new();
    P->U=bsr_new();
    return P;
}

void bildu_init(bildu_prec *P, bsr_matrix *A)
{
    int b     = A->b;
    int nrows = A->nrows;
    int nnz   = A->nnz;

    bsr_matrix *L=P->L;
    bsr_matrix *D=P->D;
    bsr_matrix *U=P->U;

    bsr_init(P->L, nrows, (nnz-nrows)/2, b);
    bsr_init(P->D, nrows, nrows, b);
    bsr_init(P->U, nrows, (nnz-nrows)/2, b);

    // Sparsity of L,D, and U components of matrix A.
    // Sparsity of L is expressed for its transpose
    // which is identical to the sparsity of U due to
    // the assumption of structural symmetry.
    nnz=0;
    for (int i=0;i<nrows;i++)
    {
        L->rowptr[i]=nnz;
        D->rowptr[i]=i;
        U->rowptr[i]=nnz;
        for (int k=A->rowptr[i];k<A->rowptr[i+1];k++)
        {
            int j=A->colidx[k];
            if(j>i)
            {
                L->colidx[nnz]=j;
                U->colidx[nnz]=j;
                nnz++;
            }
        }
        D->colidx[i]=i;
    }
    L->rowptr[nrows]=nnz;
    D->rowptr[nrows]=nrows;
    U->rowptr[nrows]=nnz;

    L->nnz=nnz;
    D->nnz=nrows;
    U->nnz=nnz;

    // allocate values arrays
    //int bb=b*b;
    //L->dbl = malloc(bb*nnz*sizeof(double));
    //D->dbl = malloc(bb*nrows*sizeof(double));
    //U->dbl = malloc(bb*nnz*sizeof(double));
}

void vec_copy(double *y, double const * x, int n)
{
    for(int i=0;i<n;i++) y[i]=x[i];
}

void mat3_view(double const *M, char const *name)
{
    printf("%s = \n[\n",name);
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            printf("%+.4e ",M[3*j+i]);
        }
        printf("\n");
    }
    printf("\n]\n");
}

void mat3_T(double *B, double const *A)
{
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) B[3*i+j] = A[i+3*j];
}

void mat3_matmul(double *C, const double *A, const double *B)
{
    // assume 3x3 column-major matrices
    // account for possiblity of C and A referring to same memory location
    double M[9];
    for(int k=0;k<9;k++) M[k]=0.0;
    for(int j=0;j<3;j++)
    {
        double *m_j = M+3*j;        // j-th column of M
        for(int k=0;k<3;k++)
        {
            double b_kj = B[k+3*j]; // kj-th element of B
            double const *a_k = A+3*k;    // k-th column of A
            for(int i=0;i<3;i++) m_j[i] += a_k[i]*b_kj; // |m_j> += |a_k> * b_kj;
        }
    }
    for(int k=0;k<9;k++) C[k]=M[k];
}

void mat3_matfms(double *C, const double *A, const double *B)
{
    // assume 3x3 column-major matrices
    // account for possiblity of C and A referring to same memory location
    double M[9];
    for(int k=0;k<9;k++) M[k]=0.0;
    for(int j=0;j<3;j++)
    {
        double *m_j = M+3*j;        // j-th column of M
        for(int k=0;k<3;k++)
        {
            double b_kj = B[k+3*j]; // kj-th element of B
            double const *a_k = A+3*k;    // k-th column of A
            for(int i=0;i<3;i++) m_j[i] += a_k[i]*b_kj; // |m_j> += |a_k> * b_kj;
        }
    }
    for(int k=0;k<9;k++) C[k]-=M[k];
}


void mat3_inv(double *invA, const double *A)
{
    double M[9];
    M[0] = +(A[4]*A[8]-A[5]*A[7]);
    M[3] = -(A[3]*A[8]-A[5]*A[6]);
    M[6] = +(A[3]*A[7]-A[4]*A[6]);
    M[1] = -(A[1]*A[8]-A[2]*A[7]);
    M[4] = +(A[0]*A[8]-A[2]*A[6]);
    M[7] = -(A[0]*A[7]-A[1]*A[6]);
    M[2] = +(A[1]*A[5]-A[2]*A[4]);
    M[5] = -(A[0]*A[5]-A[2]*A[3]);
    M[8] = +(A[0]*A[4]-A[1]*A[3]);

    double detA = A[0]*M[0]+A[1]*M[3]+A[2]*M[6];
    for(int k=0;k<9;k++) invA[k]=M[k]/detA;
}

void inline vec_copy9(double *y, double const *x)
{
    for(int i=0;i<9;i++) y[i]=x[i];
}

void mat3_rmul(double *A, double const *B)
{
    // load left hand matrix
    __m256d vA[3];
    vA[0] = _mm256_loadu_pd(A+0);
    vA[1] = _mm256_loadu_pd(A+3);
    vA[2] = _mm256_loadu_pd(A+6);

    for(int j=0;j<3;j++)
    {
        // load column j of B matrix
        __m256d vbj   = _mm256_loadu_pd(B+3*j);

        // multiply matrix A with column j of matrix B
        __m256d vAB[3];
        vAB[0] = vA[0]*_mm256_permute4x64_pd(vbj,0b00000000); //0b01010101
        vAB[1] = vA[1]*_mm256_permute4x64_pd(vbj,0b01010101); //0b01010101
        vAB[2] = vA[2]*_mm256_permute4x64_pd(vbj,0b10101010); //0b01010101

        __m256d vz = vAB[0] + vAB[1] + vAB[2];

        // Store result in  column j of matrix A
        double z[4];
        _mm256_store_pd(z,vz);
        for(int k=0;k<3;k++) A[3*j+k]=z[k];
    }
}

void mat3_lmul(double const *A, double *B)
{
    // load left hand matrix
    __m256d vA[3];
    vA[0] = _mm256_loadu_pd(A+0);
    vA[1] = _mm256_loadu_pd(A+3);
    vA[2] = _mm256_loadu_pd(A+6);

    for(int j=0;j<3;j++)
    {
        // load column j of B matrix
        __m256d vbj   = _mm256_loadu_pd(B+3*j);

        // multiply matrix A with column j of matrix B
        __m256d vAB[3];
        vAB[0] = vA[0]*_mm256_permute4x64_pd(vbj,0b00000000); //0b01010101
        vAB[1] = vA[1]*_mm256_permute4x64_pd(vbj,0b01010101); //0b01010101
        vAB[2] = vA[2]*_mm256_permute4x64_pd(vbj,0b10101010); //0b01010101

        __m256d vz = vAB[0] + vAB[1] + vAB[2];

        // Store result in  column j of matrix B
        double z[4];
        _mm256_store_pd(z,vz);
        for(int k=0;k<3;k++) B[3*j+k]=z[k];
    }
}

void mat3_vfms(double *C, double const *A, double const *B)
{
    // load left hand matrix
    __m256d vA[3];
    vA[0] = _mm256_loadu_pd(A+0);
    vA[1] = _mm256_loadu_pd(A+3);
    vA[2] = _mm256_loadu_pd(A+6);

    for(int j=0;j<3;j++)
    {
        // load column j of B matrix
        __m256d vbj   = _mm256_loadu_pd(B+3*j);

        // multiply matrix A with column j of matrix B
        __m256d vAB[3];
        vAB[0] = vA[0]*_mm256_permute4x64_pd(vbj,0b00000000); //0b01010101
        vAB[1] = vA[1]*_mm256_permute4x64_pd(vbj,0b01010101); //0b01010101
        vAB[2] = vA[2]*_mm256_permute4x64_pd(vbj,0b10101010); //0b01010101

        __m256d vz = vAB[0] + vAB[1] + vAB[2];

        // Store result in  column j of matrix A
        double z[4];
        _mm256_store_pd(z,vz);
        for(int k=0;k<3;k++) C[3*j+k]-=z[k];
    }
}


void bildu_factorize(bildu_prec *P, bsr_matrix *A)
{
    int nrows = A->nrows;
    int b     = A->b;
    int bb    = b*b;

    bsr_matrix *L=P->L;
    bsr_matrix *D=P->D;
    bsr_matrix *U=P->U;

    // Splitting values of A into L, D, and U, respectively
    int kU=0;
    for(int i=0;i<nrows;i++)
    {
        for (int k=A->rowptr[i];k<A->rowptr[i+1];k++)
        {
            int j=A->colidx[k];
            if(j<i)       // struct-transpose of L
            {
                int kL = L->rowptr[j];
                vec_copy9(L->dbl + bb*kL, A->dbl + bb*k);
                L->rowptr[j]++;
            }
            else if(j==i) // struct-copy of D
            {
                vec_copy9(D->dbl + bb*i, A->dbl + bb*k);
            }
            else if(j>i) // struct-copy of U
            {
                vec_copy9(U->dbl + bb*kU, A->dbl + bb*k);
                kU++;
            }
        }
    }
    // reset rowptr of L
    for(int i=nrows;i>0;i--) L->rowptr[i]=L->rowptr[i-1];
    L->rowptr[0]=0;

    // Factorizing
    double scale[9]; //hard-coded to 3x3 blocks for now
    for(int i=0;i<A->nrows;i++)
    {
        mat3_inv(scale,D->dbl+i*bb);
        vec_copy9(D->dbl+bb*i, scale); //store inverse instead to simplify application
        for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
        {
            //scale column i of L
            mat3_rmul(L->dbl+k*bb,scale);

            //update diagonal of U
            int j=L->colidx[k];
            mat3_vfms(D->dbl+j*bb,L->dbl+k*bb,U->dbl+k*bb);

            //scale row i of U
            mat3_lmul(scale,U->dbl+k*bb);

            //NOT IMPLEMENTED!
            for(int m=L->rowptr[j];m<L->rowptr[j+1];m++)
            {
                if(L->colidx[m]==j)
                {
                    printf("ILU OFF_DIAGONALS NOT IMPLEMENTED!\n");
                    printf("(%d,%d)",m,j);
                    getchar();
                }
            }
        }
    }
}


inline void mat3_vecmul(const double *A, double *x)
{
    const int b=3;
    double z[3];
    for(int k=0;k<3;k++) z[k]=0;
    for(int c=0;c<b;c++)
    {
        for(int r=0;r<b;r++)
        {
            z[r]+=A[c*b+r]*x[c];
        }
    }
    for(int k=0;k<3;k++) x[k]=z[k];
}


inline void mat3_vecfms(double *y, const double *A, const double *x)
{
    const int b=3;
    double z[3];
    for(int k=0;k<3;k++) z[k]=0;
    for(int c=0;c<b;c++)
    {
        for(int r=0;r<b;r++)
        {
            z[r]+=A[b*c+r]*x[c];
        }
    }
    for(int k=0;k<3;k++) y[k]-=z[k];
}

void bildu_apply3(bildu_prec *restrict P, double *x)
{
    //bsr_matrix *LT = &(P->LT);
    bsr_matrix const *L  = P->L;
    bsr_matrix const *D  = P->D;
    bsr_matrix const *U  = P->U;

    int b=L->b;
    int bb=b*b;

    // Lower triangular solve assuming ones on diagonal
    for(int i=0;i<L->ncols;i++)
    {
        for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
        {
            int j=L->colidx[k];
            mat3_vecfms(x+b*j,L->dbl+k*bb,x+b*i);
        }
    }

    // Muliply by (inverse) diagonal matrix
    for (int i=0;i<D->ncols;i++) mat3_vecmul(D->dbl+bb*i,x+b*i);

    // Upper triangular solve assuming ones on diagonal
    for(int i=U->ncols;i>0;i--)
    {
        for(int k=U->rowptr[i]-1;k>U->rowptr[i-1]-1;k--)
        {
            int j=U->colidx[k];
            mat3_vecfms(x+b*(i-1),U->dbl+k*bb,x+b*j);
        }
    }
}

void bildu_info(bildu_prec *P)
{
    bsr_info(P->L);
    bsr_info(P->D);
    bsr_info(P->U);
}

void vec_show(const double *x, int n, const char *name)
{
    printf("%s =\n[\n",name);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<3;j++) printf(" %+.4e",x[3*i+j]);
        printf("\n");
    }
    printf("]\n\n");
}

void headtail(double *x, int n, char const *name)
{

    int const  depth=9;
    printf("%s =\n[\n",name);
    for(int i=0;i<depth;i++)
    {
        printf(" %+.4e\n",x[i]);
    }
    printf("...\n");
    for(int i=n-depth;i<n;i++)
    {
        printf(" %+.4e\n",x[i]);
    }
    printf("]\n");

}
