#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void mat_show(double const *A, int n, char const *name)
{
    printf("%s = [\n",name);
    for(int i = 0;i<n;i++)
    {
        for(int j = 0;j<n;j++) printf(" %+.4e",A[i+n*j]);
        printf("\n");
    }
    printf("]\n\n");
}
#if 0
void mat2_rmul(double *A, double const *B)
{
    double M[4];
    M[0] = A[0]*B[0] + A[2]*B[1];
    M[1] = A[1]*B[0] + A[3]*B[1];
    M[2] = A[0]*B[2] + A[2]*B[3];
    M[3] = A[1]*B[2] + A[3]*B[3];
    for(int k=0;k<4;k++) A[k]=M[k];
}

void mat2_lmul(double const *A, double *B)
{
    double M[4];
    M[0] = A[0]*B[0] + A[2]*B[1];
    M[1] = A[1]*B[0] + A[3]*B[1];
    M[2] = A[0]*B[2] + A[2]*B[3];
    M[3] = A[1]*B[2] + A[3]*B[3];
    for(int k=0;k<4;k++) B[k]=M[k];
}

void mat2_inv(double *A, double const *B)
{
    double M[4];
    M[0] =  B[3];
    M[1] = -B[1];
    M[2] = -B[2];
    M[3] =  B[0];

    double inv_det = 1.0/(M[0]*M[3]-M[1]*M[2]);
    for(int k=0;k<4;k++) A[k]=inv_det*M[k];
}
#endif

/**
 * @brief Matrix inverse for 4x4 matrix.
 *
 * @param invA Pointer to inverse matrix.
 * @param    A Pointer to input matrix.
 */
void mat4_inv(double *invA, const double *A)
{
    double M[16];
    for(int k=0;k<16;k++) M[k] = A[k];

    for(int k=0;k<4;k++)
    {
        double scale=-1.0/M[5*k];        
        for(int i=0;i<4;i++) M[i+4*k] *= i==k?0:scale; // scale column k
        for(int j=0;j<4;j++)
        {
            if (j==k) continue;
            for(int i=0;i<4;i++) M[i+4*j] += i==k?0:M[i+4*k]*M[k+4*j]; //sweep
        }
        scale=-scale;
        for(int j=0;j<4;j++) M[k+4*j] *= scale; // scale row k
        M[5*k] = scale;
    }

    for(int k=0;k<16;k++) invA[k] = M[k];
}

/**
 * @brief In-place right matrix-matrix multiplication for 4x4 matrices.
 *
 * @note Function is optimized with avx2 intrinsics
 *
 * @param A Pointer to left input and output matrix.
 * @param B Pointer to right input matrix.
 */
void mat4_rmul(double *A, double const *B)
{
    // load left hand matrix
    __m256d vA[4];
    vA[0] = _mm256_loadu_pd(A+ 0);
    vA[1] = _mm256_loadu_pd(A+ 4);
    vA[2] = _mm256_loadu_pd(A+ 8);
    vA[3] = _mm256_loadu_pd(A+12);

    for(int j=0;j<4;j++)
    {
        // load column j of B matrix
        __m256d vbj   = _mm256_loadu_pd(B+4*j);

        // multiply matrix A with column j of matrix B
        __m256d vAB[4];
        vAB[0] = vA[0]*_mm256_permute4x64_pd(vbj,0x00);
        vAB[1] = vA[1]*_mm256_permute4x64_pd(vbj,0x55); // 0b01010101
        vAB[2] = vA[2]*_mm256_permute4x64_pd(vbj,0xAA); // 0b10101010
        vAB[3] = vA[3]*_mm256_permute4x64_pd(vbj,0xFF); // 0b10101010

        __m256d vz = vAB[0] + vAB[1] + vAB[2] + vAB[3];

        // Store result in  column j of matrix A
        _mm256_storeu_pd(A+4*j,vz);
    }
}

/**
 * @brief In-place left matrix-matrix multiplication for 4x4 matrices.
 *
 * @note Function is optimized with avx2 intrinsics
 *
 * @param A Pointer to left input  matrix.
 * @param B Pointer to right input and output matrix.
 */
void mat4_lmul(double const *A, double *B)
{
    // load left hand matrix
    __m256d vA[4];
    vA[0] = _mm256_loadu_pd(A+ 0);
    vA[1] = _mm256_loadu_pd(A+ 4);
    vA[2] = _mm256_loadu_pd(A+ 8);
    vA[3] = _mm256_loadu_pd(A+12);

    for(int j=0;j<4;j++)
    {
        // load column j of B matrix
        __m256d vbj   = _mm256_loadu_pd(B+4*j);

        // multiply matrix A with column j of matrix B
        __m256d vAB[4];
        vAB[0] = vA[0]*_mm256_permute4x64_pd(vbj,0x00);
        vAB[1] = vA[1]*_mm256_permute4x64_pd(vbj,0x55); // 0b01010101
        vAB[2] = vA[2]*_mm256_permute4x64_pd(vbj,0xAA); // 0b10101010
        vAB[3] = vA[3]*_mm256_permute4x64_pd(vbj,0xFF); // 0b10101010

        __m256d vz = vAB[0] + vAB[1] + vAB[2] + vAB[3];

        // Store result in  column j of matrix B
        _mm256_storeu_pd(B+4*j,vz);
    }
}

/**
 * @brief Fused multiply-subtract for 4x4 matrices.
 *
 * @note Function is optimized with avx2 intrinsics
 *
 * @param C Pointer to output matrix.
 * @param A Pointer to left input matrix.
 * @param B Pointer to right input matrix.
 */
void mat4_vfms(double *C, double const *A, double const *B)
{
    // load left hand matrix
    __m256d vA[4];
    vA[0] = _mm256_loadu_pd(A+ 0);
    vA[1] = _mm256_loadu_pd(A+ 4);
    vA[2] = _mm256_loadu_pd(A+ 8);
    vA[3] = _mm256_loadu_pd(A+12);

    for(int j=0;j<4;j++)
    {
        // load column j of B matrix
        __m256d vbj   = _mm256_loadu_pd(B+4*j);

        // multiply matrix A with column j of matrix B
        __m256d vAB[4];
        vAB[0] = vA[0]*_mm256_permute4x64_pd(vbj,0x00);
        vAB[1] = vA[1]*_mm256_permute4x64_pd(vbj,0x55); // 0b01010101
        vAB[2] = vA[2]*_mm256_permute4x64_pd(vbj,0xAA); // 0b10101010
        vAB[3] = vA[3]*_mm256_permute4x64_pd(vbj,0xFF); // 0b10101010

        __m256d vz = vAB[0] + vAB[1] + vAB[2] + vAB[3];

        // Subtract from column j of matrix C
        vz = _mm256_loadu_pd(C+4*j) - vz;

        // Store result in column j of matrix C
        _mm256_store_pd(C+4*j,vz);
    }
}


#ifdef __cplusplus
}
#endif
