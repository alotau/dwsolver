/*
 * test_blas.c — Unit tests for dw_blas.c
 *
 * Covers all five functions: dw_daxpy, dw_ddot, dw_dcopy, dw_ddoti,
 * dw_dcoogemv.  Each function has a normal case, a boundary/edge case
 * (len=0 or nz=0), and a case with a negative or zero coefficient.
 *
 * Build:  make check   (via tests/Makefile.am wired into configure.ac)
 * Run:    ./tests/test_blas
 * Exit:   0 on all-pass, 1 on first failure
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dw_blas.h"

/* Absolute-tolerance floating-point comparison. */
#define EXPECT_APPROX(name, expected, got, tol) \
    do { \
        double _e = (expected), _g = (got), _d = _e - _g; \
        if (_d < -(tol) || _d > (tol)) { \
            fprintf(stderr, "FAIL [%s]: expected %.15g, got %.15g\n", \
                    (name), _e, _g); \
            exit(1); \
        } \
    } while(0)

/* Exact integer comparison. */
#define EXPECT_INT(name, expected, got) \
    do { \
        if ((expected) != (got)) { \
            fprintf(stderr, "FAIL [%s]: expected %d, got %d\n", \
                    (name), (int)(expected), (int)(got)); \
            exit(1); \
        } \
    } while(0)

#define TOL 1e-12

/* -------------------------------------------------------------------------
 * dw_daxpy:  y <- alpha*x + y
 * ---------------------------------------------------------------------- */
static void test_daxpy(void) {
    /* Normal: len=3, alpha=2.0 */
    double x3[3] = {1.0, 2.0, 3.0};
    double y3[3] = {10.0, 20.0, 30.0};
    dw_daxpy(3, 2.0, x3, y3);
    EXPECT_APPROX("dw_daxpy: len=3, alpha=2.0, y[0]", 12.0, y3[0], TOL);
    EXPECT_APPROX("dw_daxpy: len=3, alpha=2.0, y[1]", 24.0, y3[1], TOL);
    EXPECT_APPROX("dw_daxpy: len=3, alpha=2.0, y[2]", 36.0, y3[2], TOL);

    /* Boundary: len=0 — y must be unchanged */
    double xz[1] = {99.0};
    double yz[1] = {7.0};
    dw_daxpy(0, 3.0, xz, yz);
    EXPECT_APPROX("dw_daxpy: len=0, y[0] unchanged", 7.0, yz[0], TOL);

    /* alpha=0.0 — y must be unchanged */
    double xa[3] = {1.0, 2.0, 3.0};
    double ya[3] = {5.0, 6.0, 7.0};
    dw_daxpy(3, 0.0, xa, ya);
    EXPECT_APPROX("dw_daxpy: alpha=0.0, y[0] unchanged", 5.0, ya[0], TOL);
    EXPECT_APPROX("dw_daxpy: alpha=0.0, y[2] unchanged", 7.0, ya[2], TOL);

    /* Negative alpha: alpha=-1.0 */
    double xn[2] = {3.0, 4.0};
    double yn[2] = {10.0, 10.0};
    dw_daxpy(2, -1.0, xn, yn);
    EXPECT_APPROX("dw_daxpy: alpha=-1.0, y[0]", 7.0, yn[0], TOL);
    EXPECT_APPROX("dw_daxpy: alpha=-1.0, y[1]", 6.0, yn[1], TOL);
}

/* -------------------------------------------------------------------------
 * dw_ddot:  return x . y
 * ---------------------------------------------------------------------- */
static void test_ddot(void) {
    /* Normal: len=3, known dot product */
    double x3[3] = {1.0, 2.0, 3.0};
    double y3[3] = {4.0, 5.0, 6.0};
    EXPECT_APPROX("dw_ddot: len=3, normal", 32.0, dw_ddot(3, x3, y3), TOL);

    /* Boundary: len=0 — must return 0.0 */
    double xz[1] = {1.0};
    double yz[1] = {1.0};
    EXPECT_APPROX("dw_ddot: len=0 returns 0.0", 0.0, dw_ddot(0, xz, yz), TOL);

    /* Orthogonal vectors — dot product = 0 */
    double xo[2] = {1.0, 0.0};
    double yo[2] = {0.0, 1.0};
    EXPECT_APPROX("dw_ddot: orthogonal vectors = 0", 0.0, dw_ddot(2, xo, yo), TOL);

    /* Negative component */
    double xneg[3] = {-1.0, 2.0, -3.0};
    double yneg[3] = { 1.0, 1.0,  1.0};
    EXPECT_APPROX("dw_ddot: negative components", -2.0, dw_ddot(3, xneg, yneg), TOL);
}

/* -------------------------------------------------------------------------
 * dw_dcopy:  y <- x
 * ---------------------------------------------------------------------- */
static void test_dcopy(void) {
    /* Normal: len=3 */
    double x3[3] = {1.5, -2.5, 3.5};
    double y3[3] = {0.0,  0.0, 0.0};
    dw_dcopy(3, x3, y3);
    EXPECT_APPROX("dw_dcopy: len=3, y[0]", 1.5,  y3[0], TOL);
    EXPECT_APPROX("dw_dcopy: len=3, y[1]", -2.5, y3[1], TOL);
    EXPECT_APPROX("dw_dcopy: len=3, y[2]", 3.5,  y3[2], TOL);

    /* Boundary: len=0 — target must not be modified */
    double xz[1] = {42.0};
    double yz[1] = {99.0};
    dw_dcopy(0, xz, yz);
    EXPECT_APPROX("dw_dcopy: len=0, target unchanged", 99.0, yz[0], TOL);
}

/* -------------------------------------------------------------------------
 * dw_ddoti:  sparse dot product:  sum_i x[i] * y[ind_x[i]]
 * ---------------------------------------------------------------------- */
static void test_ddoti(void) {
    /* Normal: nz=2 */
    double xs[2] = {2.0, 3.0};
    int    is[2] = {1, 3};
    double y[5]  = {10.0, 20.0, 30.0, 40.0, 50.0};
    /* 2.0 * y[1] + 3.0 * y[3] = 2*20 + 3*40 = 40 + 120 = 160 */
    EXPECT_APPROX("dw_ddoti: nz=2, normal", 160.0, dw_ddoti(2, xs, is, y), TOL);

    /* Boundary: nz=0 — must return 0.0 */
    EXPECT_APPROX("dw_ddoti: nz=0 returns 0.0", 0.0, dw_ddoti(0, xs, is, y), TOL);

    /* Single element: nz=1 */
    double xs1[1] = {-5.0};
    int    is1[1] = {2};
    /* -5.0 * y[2] = -5 * 30 = -150 */
    EXPECT_APPROX("dw_ddoti: nz=1, single element", -150.0, dw_ddoti(1, xs1, is1, y), TOL);
}

/* -------------------------------------------------------------------------
 * dw_dcoogemv:  sparse matrix-vector product A*b = c  (COO format)
 *
 *   A is m x k.  val[j], indx[j], jndx[j] define non-zero (row, col, value).
 *   c is zero-initialised by the function.
 * ---------------------------------------------------------------------- */
static void test_dcoogemv(void) {
    /* Single non-zero entry in a 2x2 matrix:
     *   A = [0  2]        b = [1]     c = A*b = [2]
     *       [0  0]            [1]               [0]
     */
    double val1[1] = {2.0};
    int    ix1[1]  = {0};   /* row 0 */
    int    jx1[1]  = {1};   /* col 1 */
    double b1[2]   = {1.0, 1.0};
    double c1[2]   = {0.0, 0.0};
    dw_dcoogemv(2, 2, val1, ix1, jx1, 1, b1, c1);
    EXPECT_APPROX("dw_dcoogemv: single nnz, c[0]", 2.0, c1[0], TOL);
    EXPECT_APPROX("dw_dcoogemv: single nnz, c[1]", 0.0, c1[1], TOL);

    /* Full 2x2 dense matrix in COO:
     *   A = [1  2]        b = [3]     c = [1*3+2*4, 3*3+4*4] = [11, 25]
     *       [3  4]            [4]
     */
    double val2[4] = {1.0, 2.0, 3.0, 4.0};
    int    ix2[4]  = {0, 0, 1, 1};
    int    jx2[4]  = {0, 1, 0, 1};
    double b2[2]   = {3.0, 4.0};
    double c2[2]   = {99.0, 99.0};  /* pre-dirty to prove zeroing works */
    dw_dcoogemv(2, 2, val2, ix2, jx2, 4, b2, c2);
    EXPECT_APPROX("dw_dcoogemv: 2x2 full, c[0]", 11.0, c2[0], TOL);
    EXPECT_APPROX("dw_dcoogemv: 2x2 full, c[1]", 25.0, c2[1], TOL);

    /* Zero nnz — c must be zeroed */
    double val3[1] = {1.0};
    int    ix3[1]  = {0};
    int    jx3[1]  = {0};
    double b3[2]   = {7.0, 8.0};
    double c3[2]   = {5.0, 6.0};
    dw_dcoogemv(2, 2, val3, ix3, jx3, 0, b3, c3);
    EXPECT_APPROX("dw_dcoogemv: nnz=0, c[0] zeroed", 0.0, c3[0], TOL);
    EXPECT_APPROX("dw_dcoogemv: nnz=0, c[1] zeroed", 0.0, c3[1], TOL);
}

/* -------------------------------------------------------------------------
 * main
 * ---------------------------------------------------------------------- */
int main(void) {
    test_daxpy();
    test_ddot();
    test_dcopy();
    test_ddoti();
    test_dcoogemv();
    printf("ALL TESTS PASSED\n");
    return 0;
}
