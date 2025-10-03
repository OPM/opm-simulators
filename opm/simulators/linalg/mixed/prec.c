#include "prec.h"

#include <stdio.h>
#include <immintrin.h>

#pragma GCC push_options
#pragma GCC target("avx2")


prec_t *prec_new()
{
    prec_t *P = malloc(sizeof(prec_t));
    P->L=bsr_new();
    P->D=bsr_new();
    P->U=bsr_new();

    P->noffsets=-1;
    P->offsets=NULL;

    return P;
}

int prec_analyze(bsr_matrix *M, int (*offsets)[3])
{
    int count=0;
    for(int i=0;i<M->nrows;i++)
    {
        for(int z=M->rowptr[i];z<M->rowptr[i+1];z++)
        {
            int j = M->colidx[z];
            int match=0;
            for(int m=M->rowptr[j];m<M->rowptr[j+1];m++)
            {
                int k = M->colidx[m];
                for(int n=M->rowptr[i];n<M->rowptr[i+1];n++)
                {
                    int jjj=M->colidx[n];
                    match += (k==jjj);
                    if(k==jjj)
                    {
                        if(offsets)
                        {
                            offsets[count][0]=z;//ij
                            offsets[count][1]=n;//ik
                            offsets[count][2]=m;//jk
                        }
                        count++;
                    }
                }
            }
        }
    }
    return count;
}


void prec_init(prec_t *P, bsr_matrix const *A)
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

    // replace with asserts?
    L->nnz=nnz;
    D->nnz=nrows;
    U->nnz=nnz;


    // offsets for off-diagonal updates
    int count;
    count = prec_analyze(U,NULL);
    P->offsets = malloc(3*(count+1)*sizeof(int));
    count = prec_analyze(U,P->offsets);
    P->offsets[count][0]=U->nnz;
    P->noffsets=count;
}


/*
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
*/
/*
void mat3_T(double *B, double const *A)
{
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) B[3*i+j] = A[i+3*j];
}
*/
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

inline void vec_copy9(double *y, double const *x)
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


void prec_factorize(prec_t *P, bsr_matrix *A)
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

void prec_factorize2(prec_t *P, bsr_matrix *A)
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
    int idx=0;
    int next = P->offsets[idx][0];
    double scale[9]; //hard-coded to 3x3 blocks for now
    for(int i=0;i<A->nrows;i++)
    {
        mat3_inv(scale,D->dbl+i*bb);
        vec_copy9(D->dbl+bb*i, scale); //store inverse instead to simplify application
        for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
        {
            //scale column i of L
            mat3_rmul(L->dbl+k*bb,scale);

            //update diagonal D
            int j=L->colidx[k];
            mat3_vfms(D->dbl+j*bb,L->dbl+k*bb,U->dbl+k*bb);
        }

        while(next<U->rowptr[i+1])
        {
            int ij = P->offsets[idx][0];
            int ik = P->offsets[idx][1];
            int jk = P->offsets[idx][2];

            //update off-diagonals L and U
            mat3_vfms(U->dbl+jk*bb,L->dbl+ij*bb,U->dbl+ik*bb);
            mat3_vfms(L->dbl+jk*bb,L->dbl+ik*bb,U->dbl+ij*bb);

            //update marker
            next=P->offsets[++idx][0];
        }


        for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
        {
            //scale row i of U
            mat3_lmul(scale,U->dbl+k*bb);
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

void prec_apply3(prec_t *restrict P, double *x)
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



void prec_mapply3c(prec_t *restrict P, double *x)
{
    bsr_matrix *L  = P->L;
    bsr_matrix *D  = P->D;
    bsr_matrix *U  = P->U;

    int b=L->b;
    int bb=b*b;

    __m256d mm256_zero_pd =_mm256_setzero_pd();

    // Lower triangular solve assuming ones on diagonal
    for(int i=0;i<L->ncols;i++)
    {
        __m256d vA[3], vx[3];

        double *xi = x+b*i;
        __m256d vxi = _mm256_loadu_pd(xi);

        vx[0] = _mm256_permute4x64_pd(vxi,0b00000000); //0b01010101
        vx[1] = _mm256_permute4x64_pd(vxi,0b01010101); //0b01010101
        vx[2] = _mm256_permute4x64_pd(vxi,0b10101010); //0b01010101
        for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
        {
            const float *A = L->flt+k*bb;
            int j=U->colidx[k]; // should be L, but does not matter die to structural symmetry?
            vA[0] = _mm256_cvtps_pd(_mm_loadu_ps(A+0))*vx[0];
            vA[1] = _mm256_cvtps_pd(_mm_loadu_ps(A+3))*vx[1];
            vA[2] = _mm256_cvtps_pd(_mm_loadu_ps(A+6))*vx[2];

            double *xj = x+b*j;
            __m256d vxj = _mm256_loadu_pd(xj);
            __m256d vz = (vxj - vA[0]) - (vA[1] + vA[2]);

            //vz =_mm256_blend_pd(vxj,vz,0x7);  // 4th element unchanged
            //_mm256_storeu_pd(xj,vz);
            double z[4];
            _mm256_store_pd(z,vz);
            for(int n=0;n<3;n++) xj[n]=z[n];
        }

        // Muliply by (inverse) diagonal block
        const float *A = D->flt+i*bb;
        vA[0] = _mm256_cvtps_pd(_mm_loadu_ps(A+0))*vx[0]; //0b01010101
        vA[1] = _mm256_cvtps_pd(_mm_loadu_ps(A+3))*vx[1]; //0b01010101
        vA[2] = _mm256_cvtps_pd(_mm_loadu_ps(A+6))*vx[2]; //0b01010101
        __m256d vz = vA[0] + vA[1] + vA[2];

        //vz =_mm256_blend_pd(vxi,vz,0x7);  // 4th element unchanged
        //_mm256_storeu_pd(xi,vz);
        double z[4];
        _mm256_store_pd(z,vz);
        for(int k=0;k<3;k++) xi[k]=z[k];
    }

    // Upper triangular solve assuming nonzeros stored in original order
    for(int i=U->ncols;i>0;i--)
    {
        __m256d vA[3];
        for(int k=0;k<3;k++) vA[k]=mm256_zero_pd;
        for(int k=U->rowptr[i]-1;k>U->rowptr[i-1]-1;k--)
        {
            const float *A = U->flt+k*bb;
            int j=U->colidx[k];
            __m256d vxj = _mm256_loadu_pd(x+b*j);
            vA[0] += _mm256_cvtps_pd(_mm_loadu_ps(A+0))*_mm256_permute4x64_pd(vxj,0b00000000); //0b01010101
            vA[1] += _mm256_cvtps_pd(_mm_loadu_ps(A+3))*_mm256_permute4x64_pd(vxj,0b01010101); //0b01010101
            vA[2] += _mm256_cvtps_pd(_mm_loadu_ps(A+6))*_mm256_permute4x64_pd(vxj,0b10101010); //0b01010101
        }
        double *xi = x+b*(i-1);
        __m256d vxi = _mm256_loadu_pd(xi);
        __m256d vz = (vxi - vA[0]) - (vA[1] + vA[2]);
        vz =_mm256_blend_pd(vxi,vz,0x7);  // 4th element unchanged
        _mm256_storeu_pd(xi,vz);

        //double z[4];
        //_mm256_store_pd(z,vz);
        //for(int k=0;k<3;k++) xi[k]=z[k];
    }
}

void prec_dapply3c(prec_t *restrict P, double *x)
{
    bsr_matrix *L  = P->L;
    bsr_matrix *D  = P->D;
    bsr_matrix *U  = P->U;

    int b=L->b;
    int bb=b*b;

    __m256d mm256_zero_pd =_mm256_setzero_pd();

    // Lower triangular solve assuming ones on diagonal
    for(int i=0;i<L->ncols;i++)
    {
        __m256d vA[3], vx[3];

        double *xi = x+b*i;
        __m256d vxi = _mm256_loadu_pd(xi);

        vx[0] = _mm256_permute4x64_pd(vxi,0b00000000); //0b01010101
        vx[1] = _mm256_permute4x64_pd(vxi,0b01010101); //0b01010101
        vx[2] = _mm256_permute4x64_pd(vxi,0b10101010); //0b01010101
        for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
        {
            const double *A = L->dbl+k*bb;
            int j=U->colidx[k];
            vA[0] = _mm256_loadu_pd(A+0)*vx[0];
            vA[1] = _mm256_loadu_pd(A+3)*vx[1];
            vA[2] = _mm256_loadu_pd(A+6)*vx[2];

            double *xj = x+b*j;
            __m256d vxj = _mm256_loadu_pd(xj);
            __m256d vz = (vxj - vA[0]) - (vA[1] + vA[2]);

            //vz =_mm256_blend_pd(vxj,vz,0x7);  // 4th element unchanged
            //_mm256_storeu_pd(xj,vz);
            double z[4];
            _mm256_store_pd(z,vz);
            for(int kk=0;kk<3;kk++) xj[kk]=z[kk];
        }

        // Muliply by (inverse) diagonal block
        const double *A = D->dbl+i*bb;
        vA[0] = _mm256_loadu_pd(A+0)*vx[0]; //0b01010101
        vA[1] = _mm256_loadu_pd(A+3)*vx[1]; //0b01010101
        vA[2] = _mm256_loadu_pd(A+6)*vx[2]; //0b01010101
        __m256d vz = vA[0] + vA[1] + vA[2];

        //vz =_mm256_blend_pd(vxi,vz,0x7);  // 4th element unchanged
        //_mm256_storeu_pd(xi,vz);
        double z[4];
        _mm256_store_pd(z,vz);
        for(int k=0;k<3;k++) xi[k]=z[k];
    }

    // Upper triangular solve assuming nonzeros stored in original order
    for(int i=U->ncols;i>0;i--)
    {
        __m256d vA[3];
        for(int k=0;k<3;k++) vA[k]=mm256_zero_pd;
        for(int k=U->rowptr[i]-1;k>U->rowptr[i-1]-1;k--)
        {
            const double *A = U->dbl+k*bb;
            int j=U->colidx[k];
            __m256d vxj = _mm256_loadu_pd(x+b*j);
            vA[0] += _mm256_loadu_pd(A+0)*_mm256_permute4x64_pd(vxj,0b00000000); //0b01010101
            vA[1] += _mm256_loadu_pd(A+3)*_mm256_permute4x64_pd(vxj,0b01010101); //0b01010101
            vA[2] += _mm256_loadu_pd(A+6)*_mm256_permute4x64_pd(vxj,0b10101010); //0b01010101
        }
        double *xi = x+b*(i-1);
        __m256d vxi = _mm256_loadu_pd(xi);
        __m256d vz = (vxi - vA[0]) - (vA[1] + vA[2]);
        vz =_mm256_blend_pd(vxi,vz,0x7);  // 4th element unchanged
        _mm256_storeu_pd(xi,vz);

        //double z[4];
        //_mm256_store_pd(z,vz);
        //for(int k=0;k<3;k++) xi[k]=z[k];
    }
}








void prec_downcast(prec_t *P)
{
    bsr_downcast(P->L);
    bsr_downcast(P->D);
    bsr_downcast(P->U);
}



void prec_info(prec_t *P)
{
    bsr_info(P->L);
    bsr_info(P->D);
    bsr_info(P->U);
}
#pragma GCC pop_options
