#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void mat_fshow(float const *A, int n, char const *name)
{
    printf("%s = [\n",name);
    for(int i = 0;i<n;i++)
    {
        for(int j = 0;j<n;j++) printf(" %+.4e",A[i+n*j]);
        printf("\n");
    }
    printf("]\n\n");
}

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

/**
 * @brief In-place right matrix-matrix multiplication for 2x2 matrices.
 *
 * @param A Pointer to left input and output matrix.
 * @param B Pointer to right input matrix.
 */
void mat2_rmul(double *A, double const *B)
{
    double M[4];
    M[0] = A[0]*B[0] + A[2]*B[1];
    M[1] = A[1]*B[0] + A[3]*B[1];
    M[2] = A[0]*B[2] + A[2]*B[3];
    M[3] = A[1]*B[2] + A[3]*B[3];
    for(int k=0;k<4;k++) A[k]=M[k];
}

/**
 * @brief In-place left matrix-matrix multiplication for 2x2 matrices.
 *
 * @param A Pointer to left input  matrix.
 * @param B Pointer to right input and output matrix.
 */
void mat2_lmul(double const *A, double *B)
{
    double M[4];
    M[0] = A[0]*B[0] + A[2]*B[1];
    M[1] = A[1]*B[0] + A[3]*B[1];
    M[2] = A[0]*B[2] + A[2]*B[3];
    M[3] = A[1]*B[2] + A[3]*B[3];
    for(int k=0;k<4;k++) B[k]=M[k];
}

/**
 * @brief Fused multiply-subtract for 2x2 matrices.
 *
 * @param C Pointer to output matrix.
 * @param A Pointer to left input matrix.
 * @param B Pointer to right input matrix.
 */
void mat2_vfms(double *C, double const *A, double const *B)
{
    double M[4];
    M[0] = A[0]*B[0] + A[2]*B[1];
    M[1] = A[1]*B[0] + A[3]*B[1];
    M[2] = A[0]*B[2] + A[2]*B[3];
    M[3] = A[1]*B[2] + A[3]*B[3];
    for(int k=0;k<4;k++) C[k]-=M[k];
}

/**
 * @brief Matrix inverse for 2x2 matrix.
 *
 * @param invA Pointer to inverse matrix.
 * @param    A Pointer to input matrix.
 */
void mat2_inv(double *invA, double const *A)
{
    double M[4];
    M[0] =  A[3];
    M[1] = -A[1];
    M[2] = -A[2];
    M[3] =  A[0];

    double inv_det = 1.0/(M[0]*M[3]-M[1]*M[2]);
    for(int k=0;k<4;k++) invA[k]=inv_det*M[k];
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
        _mm256_storeu_pd(C+4*j,vz);
    }
}


/**
 * @brief Matrix inverse for 4x4 matrix.
 *
 * @note This is a generic implementation without pivoting
 *       hardcoded for 4x4 matrices. Primary purpose is to
 *       confirm the avx2 implementations below
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
 * @brief Matrix inverse for 4x4 matrix.
 *
 * @note This is an avx2 double-precision translation of Intel's
 *       SSE single-precision 4x4 matrix inverse using Cramer's
 *       rule.
 *
 * @param invA Pointer to inverse matrix.
 * @param    A Pointer to input matrix.
 */
void mat4_vinv(double *invA, double const *A)
{
    __m256d col0, col1, col2, col3;
    __m256d row0, row1, row2, row3;
    __m256d tmp1, det;

    // load columns
    col0 = _mm256_load_pd(A+ 0);
    col1 = _mm256_load_pd(A+ 4);
    col2 = _mm256_load_pd(A+ 8);
    col3 = _mm256_load_pd(A+12);

    tmp1 = _mm256_unpacklo_pd(col0,col1);
    row2 = _mm256_unpacklo_pd(col2,col3);

    row0 = _mm256_permute2f128_pd(tmp1, row2, 0x20);
    row2 = _mm256_permute2f128_pd(row2, tmp1, 0x13);

    tmp1 = _mm256_unpackhi_pd(col0,col1);
    row1 = _mm256_unpackhi_pd(col2,col3);

    row3 = _mm256_permute2f128_pd(tmp1, row1, 0x13);
    row1 = _mm256_permute2f128_pd(row1, tmp1, 0x20);
    // ---------------------------------------------
    tmp1 = _mm256_mul_pd(row2, row3);
    tmp1 = _mm256_permute_pd(tmp1, 0x05);

    col0 = _mm256_mul_pd(row1, tmp1);
    col1 = _mm256_mul_pd(row0, tmp1);

    tmp1 = _mm256_permute4x64_pd(tmp1,0x4E);

    col0 = _mm256_fmsub_pd(row1,tmp1,col0);
    col1 = _mm256_fmsub_pd(row0,tmp1,col1);
    col1 = _mm256_permute4x64_pd(col1,0x4E);
    // -----------------------------------------------
    tmp1 = _mm256_mul_pd(row1, row2);
    tmp1 = _mm256_permute_pd(tmp1, 0x05);

    col0 = _mm256_fmadd_pd(row3,tmp1,col0);
    col3 = _mm256_mul_pd(row0, tmp1);

    tmp1 = _mm256_permute4x64_pd(tmp1,0x4E);

    col0 = _mm256_fnmadd_pd(row3,tmp1,col0);
    col3 = _mm256_fmsub_pd(row0,tmp1,col3);
    col3 = _mm256_permute4x64_pd(col3,0x4E);
    // -----------------------------------------------
    tmp1 = _mm256_mul_pd(_mm256_permute4x64_pd(row1, 0x4E), row3);
    tmp1 = _mm256_permute_pd(tmp1, 0x05);
    row2 = _mm256_permute4x64_pd(row2, 0x4E);

    col0 = _mm256_fmadd_pd(row2,tmp1,col0);
    col2 = _mm256_mul_pd(row0, tmp1);

    tmp1 = _mm256_permute4x64_pd(tmp1,0x4E);

    col0 = _mm256_fnmadd_pd(row2,tmp1,col0);
    col2 = _mm256_fmsub_pd(row0,tmp1,col2);
    col2 = _mm256_permute4x64_pd(col2,0x4E);
    // -----------------------------------------------
    tmp1 = _mm256_mul_pd(row0, row1);
    tmp1 = _mm256_permute_pd(tmp1, 0x05);

    col2 = _mm256_fmadd_pd(row3,tmp1,col2);
    col3 = _mm256_fmsub_pd(row2,tmp1,col3);

    tmp1 = _mm256_permute4x64_pd(tmp1,0x4E);

    col2 = _mm256_fmsub_pd(row3,tmp1,col2);
    col3 = _mm256_fnmadd_pd(row2,tmp1,col3);
    // -----------------------------------------------
    tmp1 = _mm256_mul_pd(row0, row3);
    tmp1 = _mm256_permute_pd(tmp1, 0x05);

    col1 = _mm256_fnmadd_pd(row2,tmp1,col1);
    col2 = _mm256_fmadd_pd(row1,tmp1,col2);

    tmp1 = _mm256_permute4x64_pd(tmp1,0x4E);

    col1 = _mm256_fmadd_pd(row2,tmp1,col1);
    col2 = _mm256_fnmadd_pd(row1,tmp1,col2);
    // -----------------------------------------------
    tmp1 = _mm256_mul_pd(row0, row2);
    tmp1 = _mm256_permute_pd(tmp1, 0x05);

    col1 = _mm256_fmadd_pd(row3,tmp1,col1);
    col3 = _mm256_fnmadd_pd(row1,tmp1,col3);

    tmp1 = _mm256_permute4x64_pd(tmp1,0x4E);

    col1 = _mm256_fnmadd_pd(row3,tmp1,col1);
    col3 = _mm256_fmadd_pd(row1,tmp1,col3);
// -----------------------------------------------
    det = _mm256_mul_pd(row0, col0);
    det = _mm256_add_pd(_mm256_permute4x64_pd(det, 0x4E), det);
    det = _mm256_add_pd(_mm256_permute_pd(det, 0x05), det);
    det = _mm256_permute4x64_pd(det,0x00);
    det = _mm256_div_pd(_mm256_set1_pd(1.0),det);

    _mm256_store_pd(invA+ 0,_mm256_mul_pd(det, col0));
    _mm256_store_pd(invA+ 4,_mm256_mul_pd(det, col1));
    _mm256_store_pd(invA+ 8,_mm256_mul_pd(det, col2));
    _mm256_store_pd(invA+12,_mm256_mul_pd(det, col3));

}

/**
 * @brief Matrix inverse for 4x4 matrix.
 *
 * @note This is a slightly modifed version of Intel's SSE
 *       single-precision 4x4 matrix inverse using Cramer's
 *       rule.
 *
 * @param invA Pointer to inverse matrix.
 * @param    A Pointer to input matrix.
 */
void mat4_inv2(float *A)
{
    __m128 col0, col1, col2, col3;
    __m128 row0, row1, row2, row3;
    __m128 tmp1, det;

    // load columns
    col0 = _mm_load_ps(A+ 0);
    col1 = _mm_load_ps(A+ 4);
    col2 = _mm_load_ps(A+ 8);
    col3 = _mm_load_ps(A+12);

    tmp1 = _mm_unpacklo_ps(col0,col1);
    row1 = _mm_unpacklo_ps(col2,col3);

    row0 = _mm_shuffle_ps(tmp1, row1, 0x44);
    row1 = _mm_shuffle_ps(row1, tmp1, 0xEE);

    tmp1 = _mm_unpackhi_ps(col0,col1);
    row3 = _mm_unpackhi_ps(col2,col3);

    row2 = _mm_shuffle_ps(tmp1, row3, 0x44);
    row3 = _mm_shuffle_ps(row3, tmp1, 0xEE);
    // ---------------------------------------------
    tmp1 = _mm_mul_ps(row2, row3);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);

    col0 = _mm_mul_ps(row1, tmp1);
    col1 = _mm_mul_ps(row0, tmp1);

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    col0 = _mm_sub_ps(_mm_mul_ps(row1, tmp1), col0);
    col1 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), col1);
    col1 = _mm_shuffle_ps(col1, col1, 0x4E);
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row1, row2);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);

    col0 = _mm_add_ps(_mm_mul_ps(row3, tmp1), col0);
    col3 = _mm_mul_ps(row0, tmp1);

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    col0 = _mm_sub_ps(col0, _mm_mul_ps(row3, tmp1));
    col3 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), col3);
    col3 = _mm_shuffle_ps(col3, col3, 0x4E);
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(_mm_shuffle_ps(row1, row1, 0x4E), row3);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);
    row2 = _mm_shuffle_ps(row2, row2, 0x4E);

    col0 = _mm_add_ps(_mm_mul_ps(row2, tmp1), col0);
    col2 = _mm_mul_ps(row0, tmp1);

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    col0 = _mm_sub_ps(col0, _mm_mul_ps(row2, tmp1));
    col2 = _mm_sub_ps(_mm_mul_ps(row0, tmp1), col2);
    col2 = _mm_shuffle_ps(col2, col2, 0x4E);
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row0, row1);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);

    col2 = _mm_add_ps(_mm_mul_ps(row3, tmp1), col2);
    col3 = _mm_sub_ps(_mm_mul_ps(row2, tmp1), col3);

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    col2 = _mm_sub_ps(_mm_mul_ps(row3, tmp1), col2);
    col3 = _mm_sub_ps(col3, _mm_mul_ps(row2, tmp1));
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row0, row3);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);

    col1 = _mm_sub_ps(col1, _mm_mul_ps(row2, tmp1));
    col2 = _mm_add_ps(_mm_mul_ps(row1, tmp1), col2);

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    col1 = _mm_add_ps(_mm_mul_ps(row2, tmp1), col1);
    col2 = _mm_sub_ps(col2, _mm_mul_ps(row1, tmp1));
    // -----------------------------------------------
    tmp1 = _mm_mul_ps(row0, row2);
    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0xB1);

    col1 = _mm_add_ps(_mm_mul_ps(row3, tmp1), col1);
    col3 = _mm_sub_ps(col3, _mm_mul_ps(row1, tmp1));

    tmp1 = _mm_shuffle_ps(tmp1, tmp1, 0x4E);

    col1 = _mm_sub_ps(col1, _mm_mul_ps(row3, tmp1));
    col3 = _mm_add_ps(_mm_mul_ps(row1, tmp1), col3);
// -----------------------------------------------
    det = _mm_mul_ps(row0, col0);
    det = _mm_add_ps(_mm_shuffle_ps(det, det, 0x4E), det);
    det = _mm_add_ss(_mm_shuffle_ps(det, det, 0xB1), det);
    tmp1 = _mm_rcp_ss(det);

    det = _mm_sub_ss(_mm_add_ss(tmp1, tmp1), _mm_mul_ss(det, _mm_mul_ss(tmp1, tmp1)));
    det = _mm_shuffle_ps(det, det, 0x00);

    _mm_store_ps(A+ 0,_mm_mul_ps(det, col0));
    _mm_store_ps(A+ 4,_mm_mul_ps(det, col1));
    _mm_store_ps(A+ 8,_mm_mul_ps(det, col2));
    _mm_store_ps(A+12,_mm_mul_ps(det, col3));
}

#ifdef __cplusplus
}
#endif
