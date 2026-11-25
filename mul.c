/*
 * mul.c - Matrix multiply C = A * B
 *
 * Each multiply function must have a prototype of the form:
 * void mul(int M, int N, int A[N][M], int B[M][N], int C[N][N]);
 *
 * A multiply function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */
#include <stdio.h>
#include "cachelab.h"

int is_mul(int M, int N, int A[N][M], int B[M][N], int C[N][N]);

/*
 * mul_submit - This is the solution multiply function that you will be
 * graded on for Part C of the assignment. Do not change the
 * description string "Multiply submission", as the driver searches
 * for that string to identify the multiply function to be graded.
 * C starts out holding garbage, so every element of C must be written.
 *
 * Cache facts used here (32 x 32 int matrices, 1 KB direct mapped):
 *   - a row is 128 bytes, i.e. four cache lines, so A[i][k], B[k][j] and
 *     C[i][j] all live in the set (4 * index + column_block) mod 32;
 *   - a line of A and a line of B therefore collide only when their row
 *     indices agree mod 8 and their column blocks agree;
 *   - one line of B is eight columns wide, so a fetched line should pay
 *     for as many products as possible.
 *
 * The loops walk a 3 x 4 tile of C. Three rows share every line of B the
 * tile fetches, and the three A values of a k step are loaded once into
 * locals, which stops A and B from evicting each other inside the tile.
 * Twelve accumulators plus three locals keep within the 16-int allowance.
 */
char mul_submit_desc[] = "Multiply submission";
void mul_submit(int M, int N, int A[N][M], int B[M][N], int C[N][N])
{
    int i, j, k;
    int t0, t1, t2;   /* the A values of the tile's rows, hoisted */
    int a0, a1, a2, a3, b0, b1, b2, b3, c0, c1, c2, c3;

    for (i = 0; i + 2 < N; i += 3) {
        for (j = 0; j + 3 < N; j += 4) {
            a0 = a1 = a2 = a3 = 0;
            b0 = b1 = b2 = b3 = 0;
            c0 = c1 = c2 = c3 = 0;
            for (k = 0; k < M; k++) {
                t0 = A[i][k];
                t1 = A[i + 1][k];
                t2 = A[i + 2][k];
                a0 += t0 * B[k][j + 0];
                a1 += t0 * B[k][j + 1];
                a2 += t0 * B[k][j + 2];
                a3 += t0 * B[k][j + 3];
                b0 += t1 * B[k][j + 0];
                b1 += t1 * B[k][j + 1];
                b2 += t1 * B[k][j + 2];
                b3 += t1 * B[k][j + 3];
                c0 += t2 * B[k][j + 0];
                c1 += t2 * B[k][j + 1];
                c2 += t2 * B[k][j + 2];
                c3 += t2 * B[k][j + 3];
            }
            C[i][j + 0] = a0;
            C[i][j + 1] = a1;
            C[i][j + 2] = a2;
            C[i][j + 3] = a3;
            C[i + 1][j + 0] = b0;
            C[i + 1][j + 1] = b1;
            C[i + 1][j + 2] = b2;
            C[i + 1][j + 3] = b3;
            C[i + 2][j + 0] = c0;
            C[i + 2][j + 1] = c1;
            C[i + 2][j + 2] = c2;
            C[i + 2][j + 3] = c3;
        }
        /* Columns left over when N is not a multiple of four. */
        for (j = N - N % 4; j < N; j++) {
            C[i][j] = C[i + 1][j] = C[i + 2][j] = 0;
            for (k = 0; k < M; k++) {
                C[i][j] += A[i][k] * B[k][j];
                C[i + 1][j] += A[i + 1][k] * B[k][j];
                C[i + 2][j] += A[i + 2][k] * B[k][j];
            }
        }
    }

    /* One or two rows left when N is not a multiple of three. */
    for (i = N - N % 3; i < N; i += 2) {
        if (i + 1 < N) {
            for (j = 0; j + 3 < N; j += 4) {
                a0 = a1 = a2 = a3 = 0;
                b0 = b1 = b2 = b3 = 0;
                for (k = 0; k < M; k++) {
                    t0 = A[i][k];
                    t1 = A[i + 1][k];
                    a0 += t0 * B[k][j + 0];
                    a1 += t0 * B[k][j + 1];
                    a2 += t0 * B[k][j + 2];
                    a3 += t0 * B[k][j + 3];
                    b0 += t1 * B[k][j + 0];
                    b1 += t1 * B[k][j + 1];
                    b2 += t1 * B[k][j + 2];
                    b3 += t1 * B[k][j + 3];
                }
                C[i][j + 0] = a0;
                C[i][j + 1] = a1;
                C[i][j + 2] = a2;
                C[i][j + 3] = a3;
                C[i + 1][j + 0] = b0;
                C[i + 1][j + 1] = b1;
                C[i + 1][j + 2] = b2;
                C[i + 1][j + 3] = b3;
            }
            for (j = N - N % 4; j < N; j++) {
                C[i][j] = C[i + 1][j] = 0;
                for (k = 0; k < M; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                    C[i + 1][j] += A[i + 1][k] * B[k][j];
                }
            }
        } else {
            for (j = 0; j + 3 < N; j += 4) {
                a0 = a1 = a2 = a3 = 0;
                for (k = 0; k < M; k++) {
                    t0 = A[i][k];
                    a0 += t0 * B[k][j + 0];
                    a1 += t0 * B[k][j + 1];
                    a2 += t0 * B[k][j + 2];
                    a3 += t0 * B[k][j + 3];
                }
                C[i][j + 0] = a0;
                C[i][j + 1] = a1;
                C[i][j + 2] = a2;
                C[i][j + 3] = a3;
            }
            for (j = N - N % 4; j < N; j++) {
                C[i][j] = 0;
                for (k = 0; k < M; k++)
                    C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

/*
 * mul - A simple baseline multiply function, not optimized for the cache.
 */
char mul_desc[] = "Simple i-j-k multiply";
void mul(int M, int N, int A[N][M], int B[M][N], int C[N][N])
{
    int i, j, k, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            tmp = 0;
            for (k = 0; k < M; k++)
                tmp += A[i][k] * B[k][j];
            C[i][j] = tmp;
        }
    }
}

/*
 * registerFunctions - This function registers your multiply functions
 *     with the driver. At runtime, the driver will evaluate each of
 *     the registered functions and summarize their performance. This
 *     is a handy way to experiment with different multiply strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerMulFunction(mul_submit, mul_submit_desc);

    /* Register any additional multiply functions */
    registerMulFunction(mul, mul_desc);
}

/*
 * is_mul - This helper function checks if C is the product of A and B.
 *     You can check the correctness of your multiply by calling it
 *     before returning from the multiply function.
 */
int is_mul(int M, int N, int A[N][M], int B[M][N], int C[N][N])
{
    int i, j, k, num;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            num = 0;
            for (k = 0; k < M; k++)
                num += A[i][k] * B[k][j];
            if (num != C[i][j])
                return 0;
        }
    }
    return 1;
}
