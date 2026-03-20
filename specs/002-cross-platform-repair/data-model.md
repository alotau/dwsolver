# Data Model: Global State Inventory

**Feature**: `002-cross-platform-repair`
**Date**: 2026-03-19
**Purpose**: Complete inventory of shared global state in DWSOLVER — used to scope the
`dw_globals.c` extraction and to document thread-access patterns for safety verification.

---

## Overview

DWSOLVER maintains all cross-thread shared state as global variables currently
co-located in `src/dw.h`. This document enumerates each symbol, its type, its
initializer (if any), which code owns initialization and destruction, and which
threads read/write it.

After the KD-001 fix, `dw.h` will contain only `extern` declarations; the
definitions below will live in `src/dw_globals.c`.

---

## Synchronization Primitives

| Symbol | Type | Init (function, file:line) | Destroy | Reader threads | Writer threads |
|--------|------|--------------------------|---------|----------------|----------------|
| `attr` | `pthread_attr_t` | `init_pthread_data()` `dw_support.c:181` | `pthread_attr_destroy(&attr)` `dw_main.c:729` | main | main only |
| `master_lp_ready_mutex` | `pthread_mutex_t` | `init_pthread_data()` `:183` | `pthread_mutex_destroy` `:733` | master + all sub | master + all sub |
| `master_lp_ready_cv` | `pthread_cond_t` | `init_pthread_data()` `:188` | `pthread_cond_destroy` `:738` | master + all sub | master + all sub |
| `service_queue_mutex` | `pthread_mutex_t` | `init_pthread_data()` `:184` | `pthread_mutex_destroy` `:734` | master + all sub | master + all sub |
| `next_iteration_mutex` | `pthread_mutex_t` | `init_pthread_data()` `:185` | `pthread_mutex_destroy` `:735` | master + all sub | master + all sub |
| `next_iteration_cv` | `pthread_cond_t` | `init_pthread_data()` `:189` | `pthread_cond_destroy` `:739` | master + all sub | master + all sub |
| `master_mutex` | `pthread_mutex_t` | `init_pthread_data()` `:186` | `pthread_mutex_destroy` `:736` | master + all sub | master + all sub |
| `reduced_cost_mutex` | `pthread_mutex_t` | `init_pthread_data()` `:186` | `pthread_mutex_destroy` `:736` | all sub | all sub |
| `glpk_mutex` | `pthread_mutex_t` | `init_pthread_data()` `:187` | `pthread_mutex_destroy` `:737` | master + all sub | master + all sub |
| `fputs_mutex` | `pthread_mutex_t` | `init_pthread_data()` `:187` | `pthread_mutex_destroy` `:737` | master + all sub | master + all sub |
| `sub_data_mutex` | `pthread_mutex_t*` (array) | `init_pthread_data()` `:191` malloc + init loop | destroy loop + free `dw_main.c:740–742` | all sub | all sub |
| `customers` | `sem_t*` or `sem_t` | `init_pthread_data()` — named: `sem_open`; unnamed: `sem_init` | `sem_close`/`sem_unlink` or `sem_destroy` | master + all sub | master + all sub |

### Sub-data per-thread mutex note

`sub_data_mutex` is a pointer to a heap-allocated array of `num_clients`
mutexes. It is initialized in `init_pthread_data()` after the malloc. The
pointer itself is a global; the array storage is heap-allocated.

---

## LP Problem State

| Symbol | Type | Init | Destroy | Access pattern |
|--------|------|------|---------|----------------|
| `original_master_lp` | `glp_prob*` | `dw_support.c:prepare_md()` via `lpx_read_cpxlp()` | `glp_delete_prob` `dw_main.c` | read-only during DW phases; written only once at startup |
| `master_lp` | `glp_prob*` | `glp_copy_prob(master_lp, original_master_lp, ...)` `dw_support.c:prepare_md()` | `glp_delete_prob` | written by master thread; read by subproblem threads under `glpk_mutex` |
| `parm` | `glp_iocp*` | `dw_support.c:162–168` | `free(parm)` | read-only after init; global GLPK integer control params |
| `simplex_control_params` | `glp_smcp*` | `dw_support.c:162–165` | `free(simplex_control_params)` | read-only after init; global GLPK simplex params |

---

## Structural Data

| Symbol | Type | Init | Contains |
|--------|------|------|---------|
| `D` | `D_matrix*` | `dw_support.c:prepare_D()` | Sparse matrix (CRS format) of linking constraints; rows × cols of the coupling block |
| `signals` | `signal_data*` | `dw_support.c:init_globals()` | `master_lp_ready`, `next_iteration_ready`, `current_iteration` flags shared between master and subproblem threads |

### `D_matrix` fields

```c
typedef struct D_matrix {
    double* values;    // non-zero values (CRS)
    int*    columns;   // column indices (CRS)
    int*    pointerB;  // row start pointers
    int*    pointerE;  // row end pointers
    int     rows;      // number of linking constraint rows
    int     cols;      // number of subproblem variables
    int     rows_plus; // rows + number of subproblems (D rows + convexity rows)
} D_matrix;
```

### `signal_data` fields

```c
typedef struct signal_data {
    int master_lp_ready;       // master has posted a new LP for subs to solve
    int next_iteration_ready;  // master has advanced to the next DW iteration
    int current_iteration;     // which iteration we are on (Phase I or II step)
} signal_data;
```

---

## State Initialization Order

The following initialization sequence must be preserved exactly after the KD-001 fix.
`dw_globals.c` introduces *definitions only* — no initializers — so this order is
unchanged.

```
main()
  → process_cmdline()    [dw_support.c] — CLI → faux_globals struct
  → prepare_D()          [dw_support.c] — builds D from master LP file
  → prepare_md()         [dw_support.c] — loads master LP, sets parm, simplex_control_params
  → init_globals()       [dw_support.c] — sets signals, service_queue fields
  → init_pthread_data()  [dw_support.c] — initializes attr, all mutexes, CVars, semaphore
  → thread create loop   [dw_main.c]    — spawns num_clients subproblem_thread()
  [DW iteration loop]
  → free_globals()       [dw_support.c] — frees parm, simplex_control_params, D arrays
  → thread join loop     [dw_main.c]    — joins all sub threads
  → mutex/cv/sem destroy [dw_main.c]    — destroys attr, mutexes, condvars, semaphore
  → glp_delete_prob()    [dw_main.c]    — frees original_master_lp, master_lp
```

---

## Impact of KD-001 Fix on This Inventory

**Zero runtime impact.** The fix converts the declarations in `dw.h` from
definitions to `extern` declarations, and places matching definitions in
`dw_globals.c`. This changes only the *linker input* (one TU provides storage
instead of five TUs each attempting to provide it). No values, types, init
calls, destroy calls, or access patterns change.
