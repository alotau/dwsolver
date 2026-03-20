# Data Model: Struct Layouts for Code Quality Changes

**Phase**: 1 — Design & Contracts
**Date**: 2026-03-20

This file documents the structural changes to C data types needed to support
the code quality improvements, ordered by priority.

---

## P1 — No struct changes (local variable fixes only)

The P1 bug fixes are purely local:

### Fix 1: `val = NULL` in `phase_2_iteration` (dw_phases.c)

No struct change. The fix is in the local scope of the diagnostic block:

```c
// Before (undefined behaviour):
int*    ind = malloc(sizeof(int)   * glp_get_num_rows(master_lp));
double* val = NULL;
int len = glp_get_mat_col(master_lp, col, ind, val);   // UB: val == NULL
// ... inner loops referencing val ...
free(ind);  // val never freed when allocated

// After (safe):
int     nrows = glp_get_num_rows(master_lp);
int*    ind2  = malloc(sizeof(int)    * (nrows + 1));
double* val2  = malloc(sizeof(double) * (nrows + 1));
int len2 = glp_get_mat_col(master_lp, col, ind2, val2);
for (j = 1; j <= len2; j++) {
    // ... diagnostic output ...
}
free(ind2);
free(val2);
```

*(Renamed to `ind2`/`val2` to avoid shadowing the outer `ind`/`val` arrays.)*

### Fix 2: `const char*` leak in `rounding_thread` (dw_rounding.c)

No struct change. The fix is a declaration change:

```c
// Before (wrong type, leaked allocation):
const char* local_col_name = malloc(sizeof(char)*BUFF_SIZE);
// ... later overwritten without free:
local_col_name = glp_get_col_name(my_data->sub_data.lp, j);

// After (stack buffer, no leak, no type confusion):
char local_col_name[BUFF_SIZE];
// ... later:
strncpy(local_col_name, glp_get_col_name(my_data->sub_data.lp, j), BUFF_SIZE - 1);
local_col_name[BUFF_SIZE - 1] = '\0';
```

---

## P2 — Unified `dw_iteration` function signature

### New function signature replacing `phase_1_iteration` and `phase_2_iteration`

```c
// dw_phases.h — replaces two separate declarations

/* Perform one master iteration of the Dantzig-Wolfe algorithm.
 *
 *   phase      : 1 = Phase I (drive out artificial variables),
 *                2 = Phase II (optimise full objective).
 *   first_run  : non-zero on first call to phase 1 only; triggers
 *                y-accumulator setup and sign correction.
 *   obj_names  : (phase 1 only) array to record column names for
 *                objective coefficient restoration after phase 1 ends.
 *                Pass NULL for phase 2.
 *   obj_coefs  : (phase 1 only) parallel array to obj_names. NULL for phase 2.
 *   obj_count  : (phase 1 only) pointer to running count. NULL for phase 2.
 *
 * Returns: number of columns added (> 0 = continue), 0 = done, -1 = keep
 *          going despite no columns added (simplex still making progress).
 */
int dw_iteration(subprob_struct* sub_data,
                 faux_globals*   fg,
                 master_data*    md,
                 int             phase,
                 int             first_run,
                 char**          obj_names,
                 double*         obj_coefs,
                 int*            obj_count);
```

**Callers** in `dw_main.c` change from:

```c
rc = phase_1_iteration(sub_data, globals, j ? 0 : 1,
                        obj_names, obj_coefs, obj_count, md);
// ...
rc = phase_2_iteration(sub_data, globals, md);
```

to:

```c
rc = dw_iteration(sub_data, globals, md, 1, j == 0,
                  obj_names, obj_coefs, obj_count);
// ...
rc = dw_iteration(sub_data, globals, md, 2, 0, NULL, NULL, NULL);
```

The old names `phase_1_iteration` / `phase_2_iteration` become aliases or are
removed. The `static int iteration_count` and `static int simplex_iterations`
variables that currently live inside both functions should be collapsed to a
single pair inside `dw_iteration`, reset when `phase == 1 && first_run`.

---

## P3 — `sync_state` struct (prerequisite sketching for P4)

While P4 (global state migration) is deferred, the P3 lock-efficiency changes
require identifying which mutex operations are affected.

### `pthread_rwlock_t` for `original_master_lp` access

The `master_mutex` currently protects reads of `original_master_lp` in
`subproblem_thread`. This is a write-once / read-many scenario. The replacement
is:

```c
// dw.h — add alongside existing extern declarations (P3 only, not yet in struct)
extern pthread_rwlock_t  master_data_rwlock;  // guards original_master_lp reads
```

`dw_globals.c` definition:
```c
pthread_rwlock_t master_data_rwlock = PTHREAD_RWLOCK_INITIALIZER;
```

Callers:
- `subproblem_thread` (column mapping loop): `pthread_rwlock_rdlock` /
  `pthread_rwlock_rdunlock`  
- Any future writer (none at present after master setup): would use
  `pthread_rwlock_wrlock`

This is a safe incremental step — fewer contention points, same semantic
guarantee. `init_pthread_data` gains `pthread_rwlock_init(&master_data_rwlock, NULL)` and `free_globals` gains `pthread_rwlock_destroy(&master_data_rwlock)`.

---

## P4 — `faux_globals` extension and `sync_state` struct (deferred)

This section documents the target state for the P4 global-state migration to
guide later implementation.

### `sync_state` struct (new)

```c
/* Groups all synchronisation primitives that are currently true globals.
 * Embedded in faux_globals so the master and all subproblems share one
 * instance without relying on file-scope storage.
 */
typedef struct {
    pthread_attr_t    attr;
    pthread_mutex_t   master_lp_ready_mutex;
    pthread_cond_t    master_lp_ready_cv;
    pthread_mutex_t   service_queue_mutex;
    pthread_mutex_t   next_iteration_mutex;
    pthread_cond_t    next_iteration_cv;
    pthread_mutex_t   master_mutex;
    pthread_mutex_t   reduced_cost_mutex;
    pthread_mutex_t   glpk_mutex;
    pthread_mutex_t   fputs_mutex;
    pthread_mutex_t*  sub_data_mutex;   /* malloc'd array, length = num_clients */
    pthread_rwlock_t  master_data_rwlock;
#ifdef USE_NAMED_SEMAPHORES
    sem_t*            customers;
#else
    sem_t             customers;
#endif
} sync_state;
```

### `faux_globals` additions

```c
typedef struct {
    /* ... existing fields ... */

    /* Currently-global solver objects (migrated from file scope) */
    glp_prob*     original_master_lp;
    glp_prob*     master_lp;
    glp_iocp*     parm;
    glp_smcp*     simplex_control_params;
    D_matrix*     D;
    signal_data*  signals;
    sync_state    sync;   /* embedded, not a pointer */
} faux_globals;
```

### Migration impact

All functions that currently reference the global `master_lp`, `D`, etc.
must accept `faux_globals*` (most already do) and dereference via `fg->master_lp`.
Functions that do not currently accept `faux_globals*` (e.g., `prepare_D`,
`check_degeneracy`, `purge_nonbasics`) must either be updated to accept it or
have it passed as a new parameter.

`dw_globals.c` becomes empty of content (or is removed) once all state is in
`faux_globals`.

### Transition strategy

Recommended incremental approach (not all in one PR):
1. Add `sync_state sync` to `faux_globals` (compile-time check, no behaviour change) — then migrate one mutex at a time.
2. Add `D_matrix* D` to `faux_globals` — update all callers — remove `extern D_matrix* D` from `dw.h`.
3. Add `signal_data* signals` — same pattern.
4. Add `master_lp` / `original_master_lp` — highest impact, most callers.
5. Final: remove `dw_globals.c` stubs and the extern declarations from `dw.h`.

Each step must leave the test suite green before proceeding.
