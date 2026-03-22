# API Contract: libdwsolver C Public Interface

**Version**: 0.0.0 (initial)  
**Header**: `dwsolver.h` (installed to `$(includedir)/dwsolver.h`)  
**Library**: `libdwsolver` (`.a` + `.so`/`.dylib`)  
**Feature**: 010-callable-library  
**Date**: 2026-03-22

This document is the authoritative specification for the callable C API.
All symbols described here carry the `DWSOLVER_API` visibility attribute and are
exported from the shared library. All other symbols are hidden.

---

## Header Include Guard

```c
#ifndef DWSOLVER_H_
#define DWSOLVER_H_
/* ... */
#endif /* DWSOLVER_H_ */
```

Consumers include only `<dwsolver.h>`. No other project headers need be included.

---

## Visibility Macro

```c
#if defined(_WIN32) || defined(__CYGWIN__)
#  ifdef DWSOLVER_BUILDING_LIB
#    define DWSOLVER_API __declspec(dllexport)
#  else
#    define DWSOLVER_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define DWSOLVER_API __attribute__((visibility("default")))
#else
#  define DWSOLVER_API
#endif
```

---

## Types

See [data-model.md](../data-model.md) for field-level documentation.

```c
typedef enum {
    DW_STATUS_OK            = 0,
    DW_STATUS_ERR_BAD_ARGS  = 1,
    DW_STATUS_ERR_FILE      = 2,
    DW_STATUS_ERR_INFEASIBLE = 3,
    DW_STATUS_ERR_PHASE1_LIMIT = 4,
    DW_STATUS_ERR_PHASE2_LIMIT = 5,
    DW_STATUS_ERR_UNBOUNDED  = 6,
    DW_STATUS_ERR_INTERNAL  = 99
} dw_status_t;

typedef struct {
    int    verbosity;
    double mip_gap;
    int    max_phase1_iterations;
    int    max_phase2_iterations;
    int    rounding_flag;
    int    integerize_flag;
    int    enforce_sub_integrality;
    int    print_timing_data;
    int    print_final_master;
    int    print_relaxed_sol;
    int    perturb;
    double shift;
} dw_options_t;

typedef struct {
    int     status;
    double  objective_value;
    int     num_vars;
    double* x;
} dw_result_t;
```

---

## Functions

### `dw_options_init`

```c
DWSOLVER_API void dw_options_init(dw_options_t *opts);
```

**Purpose**: Initialize `*opts` with documented default values. Always call this before
setting individual fields to guarantee forward compatibility when new fields are added.

**Parameters**:
- `opts` — Non-NULL pointer to a `dw_options_t` struct to initialize.

**Returns**: void.

**Preconditions**: `opts` must not be NULL. Passing NULL is undefined behavior.

**Postconditions**: All fields of `*opts` are set to their defaults as documented in
[data-model.md](../data-model.md).

**Errors**: None. This function always succeeds.

---

### `dw_solve`

```c
DWSOLVER_API int dw_solve(const char        *master_file,
                           const char *const *subproblem_files,
                           int                num_subproblems,
                           const dw_options_t *opts,
                           dw_result_t        *result);
```

**Purpose**: Run the Dantzig-Wolfe decomposition algorithm on the given LP problem
files and populate `*result` with the solution.

**Parameters**:

| Parameter | Type | Description |
|-----------|------|-------------|
| `master_file` | `const char*` | Path to the master LP problem file (CPLEX LP or MPS format as accepted by GLPK). Must not be NULL. |
| `subproblem_files` | `const char *const *` | Array of `num_subproblems` paths to subproblem LP files. Must not be NULL. Each element must not be NULL. |
| `num_subproblems` | `int` | Number of entries in `subproblem_files`. Must be >= 1. |
| `opts` | `const dw_options_t*` | Solver options. If NULL, library-internal defaults equivalent to `dw_options_init()` are used. |
| `result` | `dw_result_t*` | Output struct. Must not be NULL. The caller allocates; the library populates. |

**Returns**: `DW_STATUS_OK` (0) on success; a `dw_status_t` error code on failure.
The return value is also stored in `result->status`.

**Preconditions**:
- `master_file` is not NULL and points to an accessible, readable file.
- `subproblem_files` is not NULL; each element is not NULL and points to an accessible file.
- `num_subproblems` >= 1.
- `result` is not NULL.

**Postconditions on success** (`DW_STATUS_OK`):
- `result->status == DW_STATUS_OK`
- `result->objective_value` contains the optimal LP relaxation objective.
- `result->num_vars` > 0.
- `result->x` points to a heap-allocated array of `result->num_vars` doubles containing
  the primal solution. The caller must free this via `dw_result_free()`.

**Postconditions on failure**:
- `result->status` is set to the appropriate error code.
- `result->x` is NULL.
- `result->num_vars` is 0.
- No resources are leaked by the library.

**Thread safety**: NOT reentrant. Only one call may be active at any time within the
process. See Thread Safety Warning below.

**Side effects**: May write files to the current working directory when
`opts->print_final_master`, `opts->print_relaxed_sol`, or `opts->write_bases` are set.
Diagnostic output is written to stdout/stderr according to `opts->verbosity`.

---

### `dw_result_free`

```c
DWSOLVER_API void dw_result_free(dw_result_t *result);
```

**Purpose**: Release any heap memory allocated inside `*result` by a prior `dw_solve()`
call. Specifically, frees `result->x` and zeroes the struct fields.

**Parameters**:
- `result` — Pointer to the `dw_result_t` to clean up. May be NULL (no-op).

**Returns**: void.

**Preconditions**: If non-NULL, `result` must point to a `dw_result_t` that was
previously passed to `dw_solve()`. Calling `dw_result_free()` on an uninitialized or
already-freed struct is undefined behavior.

**Postconditions**:
- `result->x == NULL`
- `result->num_vars == 0`
- `result->objective_value == 0.0`
- `result->status == 0`

---

### `dw_version`

```c
DWSOLVER_API const char* dw_version(void);
```

**Purpose**: Return a null-terminated string describing the library version.

**Returns**: Pointer to a static string of the form `"dwsolver 1.2.1"`. The caller
must not free or modify this string.

**Thread safety**: Safe to call from any thread at any time.

---

## Thread Safety Warning

```
THREAD SAFETY: This library is NOT reentrant. Only one call to dw_solve()
may be in progress at any time within a process. The library uses process-wide
POSIX synchronization objects (mutexes and semaphores initialized in
dw_globals.c) that are not safe for concurrent re-entry. Callers that need
to run multiple solves from different threads must serialize calls with an
external mutex.
```

---

## Error Handling Contract

- `dw_solve()` guarantees it will not `exit()` or `abort()` on recoverable errors
  (bad arguments, missing files, infeasibility). It returns an error code instead.
- On unrecoverable errors (OOM, internal inconsistency), the library may `exit()` via
  the existing `dw_oom_abort()` mechanism. This behavior matches the current CLI.
- The caller should always check the return value of `dw_solve()`.
- **GLPK State**: This library initializes GLPK internally as part of each `dw_solve()` call. Callers that have independently initialized a GLPK environment prior to calling the library may experience undefined behavior. Callers MUST NOT hold an active GLPK environment when calling `dw_solve()`. This limitation is consistent with the concurrency restriction above.

---

## ABI Stability Promise

- The `dw_options_t` struct may grow (new fields appended at the end) in future releases.
  Always call `dw_options_init()` before use.
- The `dw_result_t` struct may grow similarly; `dw_result_free()` will handle cleanup
  regardless of which version allocated the memory.
- The function signatures above are stable and will not change in a backward-incompatible
  way within a major version.
- Internal types (`faux_globals`, `subprob_struct`, etc.) are not ABI-stable.

---

## pkg-config Usage

```sh
cc $(pkg-config --cflags dwsolver) -o my_app my_app.c $(pkg-config --libs dwsolver)
```

Produces flags equivalent to:
```
-I/usr/local/include -L/usr/local/lib -ldwsolver
```
