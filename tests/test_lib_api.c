/* *****************************************************************************
 *
 * tests/test_lib_api.c — Standalone C99 acceptance tests for libdwsolver API.
 *
 * Tests 1–5: US1 — basic API correctness (T008)
 * Tests 6–8: US3 — option-control scenarios (T014)
 *
 * Build:  make check   (via tests/Makefile.am)
 * Run:    ./tests/test_lib_api
 * Exit:   0 on all-pass, 1 on first failure
 *
 * EXAMPLE_DIR is injected at compile time by tests/Makefile.am as:
 *   -DEXAMPLE_DIR='"$(top_srcdir)/examples"'
 *
 **************************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>

#include "dw_solver.h"

/* Path to the examples directory (injected at compile time). */
#ifndef EXAMPLE_DIR
#  define EXAMPLE_DIR "../examples"
#endif

#define BERTSIMAS_DIR  EXAMPLE_DIR "/book_bertsimas"

/* -------------------------------------------------------------------------
 * Assertion helpers
 * ------------------------------------------------------------------------- */
#define ASSERT(name, cond) \
    do { \
        if (!(cond)) { \
            fprintf(stderr, "FAIL [%s]: assertion failed: %s\n", (name), #cond); \
            exit(1); \
        } \
        fprintf(stdout, "PASS [%s]\n", (name)); \
    } while(0)

#define ASSERT_EQ_INT(name, expected, got) \
    do { \
        int _e = (expected), _g = (got); \
        if (_e != _g) { \
            fprintf(stderr, "FAIL [%s]: expected %d, got %d\n", (name), _e, _g); \
            exit(1); \
        } \
        fprintf(stdout, "PASS [%s]: %d\n", (name), _g); \
    } while(0)

/* -------------------------------------------------------------------------
 * Test 1: dw_options_init sets all 12 fields to documented defaults
 * ------------------------------------------------------------------------- */
static void test_1_defaults(void) {
    dw_options_t opts;
    dw_options_init(&opts);

    ASSERT("T1 verbosity=5",               opts.verbosity               == 5);
    ASSERT("T1 mip_gap=0.01",              fabs(opts.mip_gap - 0.01)    < 1e-10);
    ASSERT("T1 max_phase1_iterations=100", opts.max_phase1_iterations   == 100);
    ASSERT("T1 max_phase2_iterations=3000",opts.max_phase2_iterations   == 3000);
    ASSERT("T1 rounding_flag=0",           opts.rounding_flag           == 0);
    ASSERT("T1 integerize_flag=0",         opts.integerize_flag         == 0);
    ASSERT("T1 enforce_sub_integrality=0", opts.enforce_sub_integrality == 0);
    ASSERT("T1 print_timing_data=0",       opts.print_timing_data       == 0);
    ASSERT("T1 print_final_master=0",      opts.print_final_master      == 0);
    ASSERT("T1 print_relaxed_sol=0",       opts.print_relaxed_sol       == 0);
    ASSERT("T1 perturb=0",                 opts.perturb                 == 0);
    ASSERT("T1 shift=0.0",                 fabs(opts.shift)             < 1e-10);
}

/* -------------------------------------------------------------------------
 * Test 2: dw_solve on book_bertsimas returns DW_STATUS_OK with valid result
 * Test 3: dw_result_free zeroes all fields and x becomes NULL
 * ------------------------------------------------------------------------- */
static void test_2_3_solve_and_free(void) {
    dw_options_t opts;
    dw_result_t  result = {0, 0.0, 0, NULL};
    int rc;

    const char *master_file = BERTSIMAS_DIR "/master.cplex";
    const char *subs[1]     = { BERTSIMAS_DIR "/single_sub.cplex" };

    dw_options_init(&opts);
    opts.verbosity        = DW_OUTPUT_SILENT;
    opts.print_final_master = 0;
    opts.print_relaxed_sol  = 0;

    rc = dw_solve(master_file, subs, 1, &opts, &result);

    ASSERT_EQ_INT("T2 return DW_STATUS_OK", DW_STATUS_OK, rc);
    ASSERT("T2 result.status == DW_STATUS_OK", result.status == DW_STATUS_OK);
    ASSERT("T2 result.num_vars > 0",           result.num_vars > 0);
    ASSERT("T2 result.x != NULL",              result.x != NULL);
    ASSERT("T2 result.objective_value finite",
           isfinite(result.objective_value));

    /* Test 3: dw_result_free */
    dw_result_free(&result);
    ASSERT("T3 result.x == NULL after free",    result.x == NULL);
    ASSERT("T3 result.num_vars == 0 after free", result.num_vars == 0);
}

/* -------------------------------------------------------------------------
 * Test 4: dw_solve with NULL master_file returns DW_STATUS_ERR_BAD_ARGS
 * ------------------------------------------------------------------------- */
static void test_4_null_master(void) {
    dw_options_t opts;
    dw_result_t  result = {0, 0.0, 0, NULL};
    const char *subs[1] = { BERTSIMAS_DIR "/single_sub.cplex" };
    int rc;

    dw_options_init(&opts);
    rc = dw_solve(NULL, subs, 1, &opts, &result);

    ASSERT_EQ_INT("T4 NULL master → ERR_BAD_ARGS",
                  DW_STATUS_ERR_BAD_ARGS, rc);
    ASSERT("T4 result.x == NULL on bad args", result.x == NULL);
}

/* -------------------------------------------------------------------------
 * Test 5: dw_version returns a string starting with "dwsolver"
 * ------------------------------------------------------------------------- */
static void test_5_version(void) {
    const char *v = dw_version();
    ASSERT("T5 dw_version != NULL", v != NULL);
    ASSERT("T5 dw_version starts with 'dwsolver'",
           strncmp(v, "dwsolver", 8) == 0);
    fprintf(stdout, "     version string: \"%s\"\n", v);
}

/* -------------------------------------------------------------------------
 * Test 6 (US3): verbosity=DW_OUTPUT_SILENT suppresses solver output
 *
 * Strategy: redirect stdout/stderr to /dev/null, run solve, restore.
 * ------------------------------------------------------------------------- */
static void test_6_silent_verbosity(void) {
    dw_options_t opts;
    dw_result_t  result = {0, 0.0, 0, NULL};
    int rc;
    int saved_stdout, saved_stderr;
    int devnull;

    const char *master_file = BERTSIMAS_DIR "/master.cplex";
    const char *subs[1]     = { BERTSIMAS_DIR "/single_sub.cplex" };

    dw_options_init(&opts);
    opts.verbosity          = DW_OUTPUT_SILENT;
    opts.print_final_master = 0;
    opts.print_relaxed_sol  = 0;

    /* Redirect stdout and stderr to /dev/null. */
    saved_stdout = dup(STDOUT_FILENO);
    saved_stderr = dup(STDERR_FILENO);
    devnull = open("/dev/null", O_WRONLY);
    dup2(devnull, STDOUT_FILENO);
    dup2(devnull, STDERR_FILENO);
    close(devnull);

    rc = dw_solve(master_file, subs, 1, &opts, &result);

    /* Restore stdout and stderr. */
    dup2(saved_stdout, STDOUT_FILENO);
    dup2(saved_stderr, STDERR_FILENO);
    close(saved_stdout);
    close(saved_stderr);

    dw_result_free(&result);
    ASSERT_EQ_INT("T6 silent verbosity → DW_STATUS_OK", DW_STATUS_OK, rc);
}

/* -------------------------------------------------------------------------
 * Test 7 (US3): custom mip_gap=0.0001 — solve still completes
 * ------------------------------------------------------------------------- */
static void test_7_custom_mip_gap(void) {
    dw_options_t opts;
    dw_result_t  result = {0, 0.0, 0, NULL};
    int rc;

    const char *master_file = BERTSIMAS_DIR "/master.cplex";
    const char *subs[1]     = { BERTSIMAS_DIR "/single_sub.cplex" };

    dw_options_init(&opts);
    opts.verbosity          = DW_OUTPUT_SILENT;
    opts.mip_gap            = 0.0001;
    opts.print_final_master = 0;
    opts.print_relaxed_sol  = 0;

    rc = dw_solve(master_file, subs, 1, &opts, &result);

    dw_result_free(&result);
    ASSERT_EQ_INT("T7 custom mip_gap → DW_STATUS_OK", DW_STATUS_OK, rc);
}

/* -------------------------------------------------------------------------
 * Test 8 (US3): opts=NULL uses internal defaults — solve succeeds
 * ------------------------------------------------------------------------- */
static void test_8_null_opts(void) {
    dw_result_t result = {0, 0.0, 0, NULL};
    int rc;
    int saved_stdout, saved_stderr;
    int devnull;

    const char *master_file = BERTSIMAS_DIR "/master.cplex";
    const char *subs[1]     = { BERTSIMAS_DIR "/single_sub.cplex" };

    /* When opts=NULL the library uses defaults (verbosity=5).
     * Suppress output so the test log stays clean. */
    saved_stdout = dup(STDOUT_FILENO);
    saved_stderr = dup(STDERR_FILENO);
    devnull = open("/dev/null", O_WRONLY);
    dup2(devnull, STDOUT_FILENO);
    dup2(devnull, STDERR_FILENO);
    close(devnull);

    rc = dw_solve(master_file, subs, 1, NULL, &result);

    dup2(saved_stdout, STDOUT_FILENO);
    dup2(saved_stderr, STDERR_FILENO);
    close(saved_stdout);
    close(saved_stderr);

    dw_result_free(&result);
    ASSERT_EQ_INT("T8 opts=NULL → DW_STATUS_OK", DW_STATUS_OK, rc);
}

/* -------------------------------------------------------------------------
 * main
 * ------------------------------------------------------------------------- */
int main(void) {
    fprintf(stdout, "=== libdwsolver API tests ===\n\n");

    fprintf(stdout, "--- Test 1: dw_options_init defaults ---\n");
    test_1_defaults();

    fprintf(stdout, "\n--- Test 2+3: dw_solve + dw_result_free (book_bertsimas) ---\n");
    test_2_3_solve_and_free();

    fprintf(stdout, "\n--- Test 4: NULL master_file argument ---\n");
    test_4_null_master();

    fprintf(stdout, "\n--- Test 5: dw_version ---\n");
    test_5_version();

    fprintf(stdout, "\n--- Test 6 (US3): silent verbosity ---\n");
    test_6_silent_verbosity();

    fprintf(stdout, "\n--- Test 7 (US3): custom mip_gap ---\n");
    test_7_custom_mip_gap();

    fprintf(stdout, "\n--- Test 8 (US3): opts=NULL uses defaults ---\n");
    test_8_null_opts();

    fprintf(stdout, "\n=== All tests passed. ===\n");
    return 0;
}
