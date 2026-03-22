/* *****************************************************************************
 *
 *   DWSOLVER - a general, parallel implementation of the Dantzig-Wolfe
 *               Decomposition algorithm.
 *
 *   Copyright 2010 United States Government National Aeronautics and Space
 *   Administration (NASA).  No copyright is claimed in the United States under
 *   Title 17, U.S. Code. All Other Rights Reserved.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   In accordance with the GNU General Public License Version 3, 29 June 2007
 *   (GPL V3) Section 7. Additional Terms, Additional Permissions are added as
 *   exceptions to the terms of the GPL V3 for this program.  These additional
 *   terms should have been received with this program in a file entitled
 *   "ADDITIONAL_LICENSE_TERMS".  If a copy was not provided, you may request
 *   one from the contact author listed below.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *   Contact: Joseph Rios <Joseph.L.Rios@nasa.gov>
 *
 **************************************************************************** */

/*
 * dw_solver.h — Public C API for libdwsolver.
 *
 * Consumers include only this header.  No other project headers are needed.
 *
 * THREAD SAFETY WARNING:
 *   libdwsolver internally uses POSIX pthreads and shares global state
 *   (GLPK problem objects, mutex/semaphore handles) across a single solve.
 *   Calling dw_solve() concurrently from multiple threads is NOT supported.
 *   Each call to dw_solve() must complete before the next call begins.
 */

#ifndef DWSOLVER_H_
#define DWSOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

/* -------------------------------------------------------------------------
 * Visibility macro (self-contained copy for consumer builds that do not
 * include the internal dw.h header).
 * ------------------------------------------------------------------------- */
#ifndef DWSOLVER_API
#if defined(_WIN32) || defined(__CYGWIN__)
#  ifdef DWSOLVER_BUILDING_LIB
#    define DWSOLVER_API __declspec(dllexport)
#  elif defined(DWSOLVER_STATIC)
#    define DWSOLVER_API
#  else
#    define DWSOLVER_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define DWSOLVER_API __attribute__((visibility("default")))
#else
#  define DWSOLVER_API
#endif
#endif /* DWSOLVER_API */

/* -------------------------------------------------------------------------
 * Verbosity level constants (mirror the internal OUTPUT_* values).
 * Pass as dw_options_t.verbosity.
 * ------------------------------------------------------------------------- */
#define DW_OUTPUT_SILENT   (-1)
#define DW_OUTPUT_QUIET      0
#define DW_OUTPUT_NORMAL     5
#define DW_OUTPUT_VERBOSE   10
#define DW_OUTPUT_ALL       15

/* -------------------------------------------------------------------------
 * Status codes
 * ------------------------------------------------------------------------- */
typedef enum {
    DW_STATUS_OK               = 0,
    DW_STATUS_ERR_BAD_ARGS     = 1,
    DW_STATUS_ERR_FILE         = 2,
    DW_STATUS_ERR_INFEASIBLE   = 3,
    DW_STATUS_ERR_PHASE1_LIMIT = 4,
    DW_STATUS_ERR_PHASE2_LIMIT = 5,
    DW_STATUS_ERR_UNBOUNDED    = 6,
    DW_STATUS_ERR_INTERNAL     = 99
} dw_status_t;

/* -------------------------------------------------------------------------
 * Solver options
 *
 * Always initialise with dw_options_init() before setting individual fields
 * to guarantee forward compatibility when new fields are added.
 * ------------------------------------------------------------------------- */
typedef struct {
    int    verbosity;               /* Output verbosity (DW_OUTPUT_*). Default: 5 */
    double mip_gap;                 /* MIP optimality gap. Default: 0.01 */
    int    max_phase1_iterations;   /* Phase I iteration cap. Default: 100 */
    int    max_phase2_iterations;   /* Phase II iteration cap. Default: 3000 */
    int    rounding_flag;           /* Round convexity vars to nearest int. Default: 0 */
    int    integerize_flag;         /* Enforce binary lambdas after LP solve. Default: 0 */
    int    enforce_sub_integrality; /* Enforce integer subproblem solutions each iter. Default: 0 */
    int    print_timing_data;       /* Print runtime timing info. Default: 0 */
    int    print_final_master;      /* Write final master LP to done.cpxlp. Default: 0 */
    int    print_relaxed_sol;       /* Write relaxed solution to file. Default: 1 */
    int    perturb;                 /* Perturb RHS (experimental). Default: 0 */
    double shift;                   /* Objective constant shift. Default: 0.0 */
} dw_options_t;

/* -------------------------------------------------------------------------
 * Solve result
 *
 * Populated by dw_solve() on success.  Call dw_result_free() to release
 * the heap-allocated x array when done.
 * ------------------------------------------------------------------------- */
typedef struct {
    int     status;           /* dw_status_t code */
    double  objective_value;  /* Optimal LP relaxation objective */
    int     num_vars;         /* Length of the x array */
    double *x;                /* Heap-allocated primal solution (caller frees via dw_result_free) */
} dw_result_t;

/* -------------------------------------------------------------------------
 * Public API
 * ------------------------------------------------------------------------- */

/**
 * dw_options_init — Initialise *opts with library defaults.
 *
 * Always call this before setting individual option fields so that any
 * fields added in future library versions are properly defaulted.
 *
 * opts must not be NULL (passing NULL is undefined behaviour).
 */
DWSOLVER_API void dw_options_init(dw_options_t *opts);

/**
 * dw_solve — Run the Dantzig-Wolfe decomposition algorithm.
 *
 * Parameters:
 *   master_file       Path to the master LP problem file (CPLEX LP format).
 *   subproblem_files  Array of num_subproblems paths to subproblem LP files.
 *   num_subproblems   Number of entries in subproblem_files (must be >= 1).
 *   opts              Solver options. Pass NULL to use library defaults.
 *   result            Output struct. Must not be NULL. Caller allocates.
 *
 * Returns DW_STATUS_OK on success; a dw_status_t error code otherwise.
 * The return value is also stored in result->status.
 *
 * On success, result->x is heap-allocated. The caller must release it with
 * dw_result_free().
 */
DWSOLVER_API int dw_solve(const char        *master_file,
                           const char *const *subproblem_files,
                           int                num_subproblems,
                           const dw_options_t *opts,
                           dw_result_t        *result);

/**
 * dw_result_free — Release heap memory owned by *result.
 *
 * Frees result->x and zeroes all fields. Safe to call on a result that
 * was never successfully populated. No-op if result is NULL.
 */
DWSOLVER_API void dw_result_free(dw_result_t *result);

/**
 * dw_version — Return a static version string.
 *
 * Returns a pointer to a null-terminated string such as "dwsolver 1.2.1".
 * The returned pointer is valid for the lifetime of the process.
 */
DWSOLVER_API const char *dw_version(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* DWSOLVER_H_ */
