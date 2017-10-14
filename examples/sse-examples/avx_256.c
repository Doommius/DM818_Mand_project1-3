//
// Created by jervelund on 9/28/17.
//

#include <immintrin.h>
#include <stdio.h>
//#include "random"
#include <time.h>


int main() {
//    sub_avx();
//    mul_avx();
//    sub_avx_PS();
    sub_avx_PD();

}

int sub_avx_PS(){

    /* Initialize the two argument vectors */
    __m256 evens = _mm256_set_ps(2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0);
    __m256 odds = _mm256_set_ps(1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0);

    /* Compute the difference between the two vectors */
    __m256 result = _mm256_sub_ps(evens, odds);

    /* Display the elements of the result vector */
    float* f = (float*)&result;
    printf("%f %f %f %f %f %f %f %f\n",
           f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);

    return 0;
}

int sub_avx_PD(){

    /* Initialize the two argument vectors */
                                    //4     3    2     1
    __m256d evens = _mm256_setr_pd(1.0, 2.0, 3.0, 4.0);
    __m256d odds = _mm256_setr_pd(0.5, 0.6, 0.7, 0.8);

    /* Compute the difference between the two vectors */
    __m256d result = _mm256_add_pd(evens, odds);

    /* Display the elements of the result vector */
    double* f = (double*)&result;
    printf("%f %f %f %f\n",
           f[0], f[1], f[2], f[3]);

    return 0;
}


int add_avx(){

    /* Initialize the two argument vectors */
    __m512d evens = _mm512_set_ps(2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0);
    __m512d odds = _mm512_set_ps(1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0);

    /* Compute the difference between the two vectors */
    __m256 result = _mm256_add_ps(evens, odds);

    /* Display the elements of the result vector */
    float* f = (float*)&result;
    printf("%f %f %f %f %f %f %f %f\n",
           f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7]);

    return 0;
}

int mul_avx(){

//    int matrix = makematrix();

    /* Initialize the two argument vectors */
    __m256d vec1 = _mm256_set_pd(1.0, 2.0, 3.0, 4.0);
    __m256d vec2 = _mm256_set_pd(2.0, 2.0, 2.0, 2.0);

    /* Compute the difference between the two vectors */
//    __m256d result = _mm256_mullo_epi64(vec1, vec2);
    __m256d  result =  _mm256_add_pd(vec1,vec2);

    /* Display the elements of the result vector */
    double* f = (double *)&result;
    printf("%f %f %f %f\n", f[0], f[1], f[2], f[3]);

    return 0;





}


__m256d multiply_vec(__m256d vec1, __m256d vec2){
    return _mm256_mul_pd(vec1,vec2);

}




int makematrix()
{
    int random[100][100];
    int i, o;

    srand(time(NULL));
    for(o = 0; o<100; o++)
        for(i = 0; i<100; i++)
            random[o][i] = rand();
    return 0;
}


/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block (int lda, int M, int N, int K, double* A, double* B, double* C)
{
    /* For each row i of A */
    for (int i = 0; i < M; ++i)
        /* For each column j of B */
        for (int j = 0; j < N; ++j)
        {
            /* Compute C(i,j) */
            double cij = C[i+j*lda];
            for (int k = 0; k < K; ++k)
                cij += A[i+k*lda] * B[k+j*lda];
            C[i+j*lda] = cij;
        }
}


/* This auxiliary subroutine performs a smaller dgemm operation
 *  C := C + A * B
 * where C is M-by-N, A is M-by-K, and B is K-by-N. */
static void do_block_avx256d (int lda, int M, int N, int K, double* A, double* B, double* C)
{
    /* For each row i of A */
    for (int i = 0; i < M; ++i)
        /* For each column j of B */
        for (int j = 0; j < N; ++j)
        {
            /* Compute C(i,j) */
            double cij = C[i+j*lda];
            for (int k = 0; k < K; ++k){
                /* Initialize the two argument vectors */
                __m256d vec1 = _mm256_set_pd(1.0, 2.0, 3.0, 4.0);
                __m256d vec2 = _mm256_set_pd(2.0, 2.0, 2.0, 2.0);

                /* Compute the difference between the two vectors */
                //    __m256d result = _mm256_mullo_epi64(vec1, vec2);
                __m256d  result =  _mm256_mul_pd(vec1,vec2);

                /* Display the elements of the result vector */
                double* f = (double *)&result;

                printf("%d %d %d %d\n", f[0], f[1], f[2], f[3]);

                cij += f[0] + f[1] + f[2] +f[3];






                cij += A[i+k*lda] * B[k+j*lda];
            }


            C[i+j*lda] = cij;
        }
}


