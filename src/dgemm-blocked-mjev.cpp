#include <mm_malloc.h>
#include <malloc.h>
#include <algorithm>

/*
  In case you're wondering, dgemm stands for Double-precision, GEneral
  Matrix-Matrix multiplication.
*/

const char *dgemm_desc = "mjerv15 blocked dgemm. based on lecture by jacob and goto paper";

//max blocks
constexpr unsigned int NR = 4;
constexpr unsigned int MR = 4;

//blocksize
constexpr unsigned int KC = 128;
constexpr unsigned int MC = 128;

unsigned int lda; //size lda*lda of matrix

double *Ablock;
double *Bblock;
double *Cslice;

//              From Lecture by Jacob
//		 				                      N
//                                +-------------------------+
//                                |                         |
//                                |                         |
//                                |                         |
//                                |                         |
//                                |                         |
//                              K |            B            |
//                                |                         |
//                                |                         |
//                                |                         |
//                                |                         |
//                                |                         |
//                                |                         |
// 				                  +-------------------------+
//		        K
//  +-------------------------+   +-------------------------+
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  |           A             |   |                         |
// M|                         |   |            C            |
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  |                         |   |                         |
//  +-------------------------+   +-------------------------+


/*
  A is M-by-K
  B is K-by-N
  C is M-by-N

  lda is the leading dimension of the matrix (the M of square_dgemm).
*/

void basic_dgemm(int lda, int M, int N, int K,
                 double *A, double *B, double *C) {
    /*
      To optimize this, think about loop unrolling and software
      pipelining.  Hint:  For the majority of the matmuls, you
      know exactly how many iterations there are (the block size)...
    */
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++) {
            double cij = C[i + j * lda];
            for (int k = 0; k < K; k++)
                cij += A[i + k * lda] * B[k + j * lda];
            C[i + j * lda] = cij;
        }
}

void do_block(int lda, double *A, double *B, double *C,
              int i, int j, int k) {
    /*
      Remember that you need to deal with the fringes in each
      dimension.

      If the matrix is 7x7 and the blocks are 3x3, you'll have 1x3,
      3x1, and 1x1 fringe blocks.

            xxxoooX
            xxxoooX
            xxxoooX
            oooxxxO
            oooxxxO
            oooxxxO
            XXXOOOX

      You won't get this to go fast until you figure out a `better'
      way to handle the fringe blocks.  The better way will be more
      machine-efficient, but very programmer-inefficient.
    */
    int M = std::min(128, lda - i);
    int N = std::min(128, lda - j);
    int K = std::min(4, lda - k);

    basic_dgemm(lda, M, N, K, A + i + k * lda, B + k + j * lda, C + i + j * lda);
}

void square_dgemm(int lda, double *A, double *B, double *C) {

    Ablock = (double *) _mm_malloc(MC * KC * sizeof(double), 16); //128*128
    Bblock = (double *) malloc(lda * KC * sizeof(double)); // M * 128
    Cslice = (double *) _mm_malloc(MC * NR * sizeof(double), 16); //128*4

    for (int k = 0; k < lda; k += KC) {
        //pack B
        for (int i = 0; i < lda; i += 128) {
            //PackA
            do_block(C, std::min(lda-i,MC), std::min(lda-k,MC));
        }




    }

}
