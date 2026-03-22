# Data Model: Callable Library Interface

**Phase**: 1 — Design  
**Feature**: 010-callable-library  
**Date**: 2026-03-22

This document describes all public types, structs, and enumerations that form the
stable ABI surface of `libdwsolver`. Internal types (`faux_globals`, `subprob_struct`,
`D_matrix`, etc.) are not listed here — they remain private to the implementation.

---

## Public Types Overview

```
dw_options_t        Solver configuration (caller populates before calling dw_solve)
dw_result_t         Solve output (populated by dw_solve; freed via dw_result_free)
dw_status_t         Integer error code enum for dw_solve return values
```

---

## `dw_options_t` — Solver Options Struct

**Purpose**: Groups all tunable parameters accepted by the solver. The caller
initializes this struct via `dw_options_init()` (which fills defaults) and then
overrides specific fields before calling `dw_solve()`.

**Stability**: Fields may be added in future versions; always initialize with
`dw_options_init()` to guarantee forward compatibility.

### Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `verbosity` | `int` | 5 (`OUTPUT_NORMAL`) | Controls diagnostic output level. Use `DW_OUTPUT_SILENT` (-1) to suppress all output. See verbosity constants in `dw_solver.h`. |
| `mip_gap` | `double` | 0.01 | MIP optimality gap tolerance. The solver stops when the relative gap between the best integer solution and the LP relaxation bound is at most this value. |
| `max_phase1_iterations` | `int` | 100 | Maximum iterations for the Phase I (feasibility) procedure. |
| `max_phase2_iterations` | `int` | 3000 | Maximum iterations for the Phase II (optimality) procedure. |
| `rounding_flag` | `int` | 0 | If 1, compute a rounded integer solution after the LP relaxation completes. |
| `integerize_flag` | `int` | 0 | If 1, enforce binary constraints on lambda variables to obtain a MIP solution. |
| `enforce_sub_integrality` | `int` | 0 | If 1, require subproblems to produce integer solutions. |
| `print_timing_data` | `int` | 0 | If 1, print per-phase wall-clock and CPU timing to stdout. |
| `print_final_master` | `int` | 0 | If 1, write the final master LP to `done.cpxlp` at solve completion. |
| `print_relaxed_sol` | `int` | 0 | If 1, write the LP relaxation solution to a results file. |
| `perturb` | `int` | 0 | If 1, apply perturbation to break degeneracy. |
| `shift` | `double` | 0.0 | Objective shift constant applied during solve. |

### Verbosity Constants (defined in `dw_solver.h`)

| Constant | Value | Meaning |
|----------|-------|---------|
| `DW_OUTPUT_SILENT` | -1 | No output at all |
| `DW_OUTPUT_QUIET` | 0 | Errors only |
| `DW_OUTPUT_NORMAL` | 5 | Standard progress messages (default) |
| `DW_OUTPUT_VERBOSE` | 10 | Detailed per-iteration output |
| `DW_OUTPUT_ALL` | 15 | Full diagnostic output |

---

## `dw_result_t` — Solve Result Struct

**Purpose**: Receives the output of a `dw_solve()` call. The caller allocates this
struct (typically on the stack) and passes a pointer to `dw_solve()`. The `x` field
is heap-allocated by the library; caller must call `dw_result_free()` when done.

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `status` | `int` | Solve status code. 0 = success. Non-zero values are `dw_status_t` constants. Populated even on failure. |
| `objective_value` | `double` | Optimal LP relaxation objective value. Valid only when `status == DW_STATUS_OK`. |
| `num_vars` | `int` | Number of entries in the `x` array (equals number of columns in the final master LP). |
| `x` | `double*` | Primal solution vector of length `num_vars`. Heap-allocated by the library. NULL if `status != DW_STATUS_OK`. Must be freed via `dw_result_free()`, not `free()` directly. |

### Ownership and Lifetime

- The `dw_result_t` struct itself is **caller-owned**; allocate on the stack or heap.
- The `x` pointer inside the struct is **library-allocated** and must be freed by calling
  `dw_result_free(&result)` before the `dw_result_t` goes out of scope.
- After `dw_result_free()`, `result.x` is set to NULL and `result.num_vars` to 0.

---

## `dw_status_t` — Status Codes

Returned by `dw_solve()` as the function return value and also stored in `result.status`.

| Constant | Value | Meaning |
|----------|-------|---------|
| `DW_STATUS_OK` | 0 | Solve completed successfully; solution is valid. |
| `DW_STATUS_ERR_BAD_ARGS` | 1 | NULL or invalid argument passed to `dw_solve()`. |
| `DW_STATUS_ERR_FILE` | 2 | Failed to open master or subproblem LP file. |
| `DW_STATUS_ERR_INFEASIBLE` | 3 | Problem is infeasible; no solution found. |
| `DW_STATUS_ERR_PHASE1_LIMIT` | 4 | Phase I iteration limit reached without feasibility. |
| `DW_STATUS_ERR_PHASE2_LIMIT` | 5 | Phase II iteration limit reached without optimality. |
| `DW_STATUS_ERR_UNBOUNDED` | 6 | Problem is unbounded; LP relaxation has no finite optimum. |
| `DW_STATUS_ERR_INTERNAL` | 99 | Unexpected internal error (e.g., memory allocation failure). |

---

## Entity Relationships

```
caller
  │
  │  stack-allocates
  ▼
dw_result_t ──── (library allocates x pointer) ──► double[] heap buffer
  ▲
  │  passed by pointer to
dw_solve()
  │
  │  reads from
  ▼
dw_options_t ──── (caller populates, often from dw_options_init() + overrides)
```

---

## Internal Types (not public, listed for cross-reference)

The following types remain internal and are **not** included in `dw_solver.h`:

| Internal Type | Location | Notes |
|---------------|----------|-------|
| `faux_globals` | `dw.h` | Wide configuration + state struct; populated internally from `dw_options_t` |
| `subprob_struct` | `dw_subprob.h` | Per-subproblem thread data |
| `D_matrix` | `dw.h` | Master constraint matrix representation |
| `master_data` | `dw.h` | Dual values and objective coefficients for master |
| `signal_data` | `dw.h` | Cross-thread iteration synchronization flags |
| `hook_struct` | `dw.h` | GLPK terminal output redirect data |

These types are not exposed to consumers and may change without affecting the public ABI.
