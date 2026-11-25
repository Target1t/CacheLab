/*
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 *
 * The cache has 32 sets and 32-byte lines, i.e. 8 ints per line. A row of
 * a 32 x 32 matrix is 128 bytes, so row i and row i + 8 land on the same
 * four sets, and a row of A and the column of B that receives it always
 * collide: reading a row of A into registers before writing B is what
 * keeps the diagonal blocks from thrashing.
 */
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/*
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded.
 *
 * Only nine int variables plus the loop counters are used: v0 .. v7 hold
 * one row of eight elements, i, j, k and t are loop counters.
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, k, t;
    int v0, v1, v2, v3, v4, v5, v6, v7;

    if (M == 32 && N == 32) {
        /*
         * 8 x 8 blocking. One row of the block is read into registers
         * and then written out, so a diagonal block does not evict the
         * line of A it is about to read again.
         */
        for (i = 0; i < 32; i += 8) {
            for (j = 0; j < 32; j += 8) {
                for (k = i; k < i + 8; k++) {
                    v0 = A[k][j + 0];
                    v1 = A[k][j + 1];
                    v2 = A[k][j + 2];
                    v3 = A[k][j + 3];
                    v4 = A[k][j + 4];
                    v5 = A[k][j + 5];
                    v6 = A[k][j + 6];
                    v7 = A[k][j + 7];
                    B[j + 0][k] = v0;
                    B[j + 1][k] = v1;
                    B[j + 2][k] = v2;
                    B[j + 3][k] = v3;
                    B[j + 4][k] = v4;
                    B[j + 5][k] = v5;
                    B[j + 6][k] = v6;
                    B[j + 7][k] = v7;
                }
            }
        }
    } else if (M == 64 && N == 64) {
        /*
         * A row of a 64 x 64 matrix is 256 bytes, so an 8 x 8 block of A
         * and the 8 x 8 block of B that receives it overlap in only a few
         * sets, but 8 x 8 blocking alone still evicts a lot. The block is
         * handled as four 4 x 4 quadrants: the top-right quadrant of A is
         * parked in the right half of B, then moved into place while the
         * bottom-left quadrant of A is written out; the bottom-right
         * quadrant is copied last.
         */
        for (i = 0; i < 64; i += 8) {
            for (j = 0; j < 64; j += 8) {
                for (k = 0; k < 4; k++) {
                    v0 = A[i + k][j + 0];
                    v1 = A[i + k][j + 1];
                    v2 = A[i + k][j + 2];
                    v3 = A[i + k][j + 3];
                    v4 = A[i + k][j + 4];
                    v5 = A[i + k][j + 5];
                    v6 = A[i + k][j + 6];
                    v7 = A[i + k][j + 7];
                    B[j + 0][i + k] = v0;
                    B[j + 1][i + k] = v1;
                    B[j + 2][i + k] = v2;
                    B[j + 3][i + k] = v3;
                    B[j + 0][i + k + 4] = v4;
                    B[j + 1][i + k + 4] = v5;
                    B[j + 2][i + k + 4] = v6;
                    B[j + 3][i + k + 4] = v7;
                }
                for (k = 0; k < 4; k++) {
                    v0 = A[i + 4][j + k];
                    v1 = A[i + 5][j + k];
                    v2 = A[i + 6][j + k];
                    v3 = A[i + 7][j + k];
                    v4 = B[j + k][i + 4];
                    v5 = B[j + k][i + 5];
                    v6 = B[j + k][i + 6];
                    v7 = B[j + k][i + 7];
                    B[j + k][i + 4] = v0;
                    B[j + k][i + 5] = v1;
                    B[j + k][i + 6] = v2;
                    B[j + k][i + 7] = v3;
                    B[j + k + 4][i + 0] = v4;
                    B[j + k + 4][i + 1] = v5;
                    B[j + k + 4][i + 2] = v6;
                    B[j + k + 4][i + 3] = v7;
                }
                for (k = 4; k < 8; k++) {
                    v0 = A[i + k][j + 4];
                    v1 = A[i + k][j + 5];
                    v2 = A[i + k][j + 6];
                    v3 = A[i + k][j + 7];
                    B[j + 4][i + k] = v0;
                    B[j + 5][i + k] = v1;
                    B[j + 6][i + k] = v2;
                    B[j + 7][i + k] = v3;
                }
            }
        }
    } else {
        /*
         * 61 x 67 and anything else: 16 rows x 8 columns. Eight columns
         * are exactly one cache line, and 16 rows of B are 16 lines in
         * sixteen different sets, so a block of B stays put while its
         * rows are written out. The eight values of a row go through
         * registers; the tail columns use a plain loop.
         */
        for (i = 0; i < N; i += 16) {
            for (j = 0; j < M; j += 8) {
                for (k = i; k < i + 16 && k < N; k++) {
                    if (j + 8 <= M) {
                        v0 = A[k][j + 0];
                        v1 = A[k][j + 1];
                        v2 = A[k][j + 2];
                        v3 = A[k][j + 3];
                        v4 = A[k][j + 4];
                        v5 = A[k][j + 5];
                        v6 = A[k][j + 6];
                        v7 = A[k][j + 7];
                        B[j + 0][k] = v0;
                        B[j + 1][k] = v1;
                        B[j + 2][k] = v2;
                        B[j + 3][k] = v3;
                        B[j + 4][k] = v4;
                        B[j + 5][k] = v5;
                        B[j + 6][k] = v6;
                        B[j + 7][k] = v7;
                    } else {
                        for (t = j; t < M; t++)
                            B[t][k] = A[k][t];
                    }
                }
            }
        }
    }
}

/*
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started.
 */

/*
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc);

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc);

}

/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}
