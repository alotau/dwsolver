# Research: DWSOLVER Code Quality Analysis

**Phase**: 0 — Research & Unknowns Resolution
**Date**: 2026-03-20
**Status**: Complete — no external research required; findings come from direct
codebase inspection of all `src/dw*.{c,h}` files (~4,000 lines).

---

## Finding 1: Undefined Behaviour — `val = NULL` in `phase_2_iteration`

**Decision**: Fix by allocating `val` before calling `glp_get_mat_col`.

**Rationale**: Inside the diagnostic block added when a new column's primal
value is non-zero, the code does:
```c
int* ind = malloc(sizeof(int)*glp_get_num_rows(master_lp));
double* val = NULL;
int len = glp_get_mat_col(master_lp, col, ind, val);
```
GLPK's `glp_get_mat_col` documents that when column indices and values are
requested, both `ind` and `val` must point to arrays of length ≥ `1 +
glp_get_num_rows(lp)`. Passing `NULL` for `val` is undefined behaviour.
On most platforms this path is silent but it can corrupt memory or produce
wrong `len` values. The block also allocates `ind` but never frees it,
and frees `ind` and `val` in an inner loop without matching outer frees.

**Alternatives considered**:
- Guarding the call with `if (val != NULL)`: rejected — the data is needed
  for the diagnostic to be meaningful.
- Removing the diagnostic block entirely: viable but loses useful debug output
  for degenerate columns; prefer a correct fix over deletion.

---

## Finding 2: Memory Leak + Wrong Type — `local_col_name` in `rounding_thread`

**Decision**: Change `const char* local_col_name = malloc(...)` to
`char* local_col_name`, print its initial `malloc` allocation, and ensure it
is freed before `pthread_exit`.

**Rationale**: The variable is declared `const char*` but initialised with
`malloc`. This generates a compiler warning and is semantically wrong — the
allocation is owned by the thread and must be freed. The variable is then
overwritten by `glp_get_col_name(...)` (which returns a pointer into GLPK's
internal storage), permanently discarding the `malloc` pointer and leaking
it. With `N` subproblems and the `--round` flag, this leaks one allocation
per subproblem.

**Alternatives considered**:
- Using a local fixed-size buffer instead of `malloc`: valid and preferred for
  a bounded filename. `char local_col_name[BUFF_SIZE]` avoids the allocation
  entirely and is the correct fix.

---

## Finding 3: `phase_1_iteration` / `phase_2_iteration` Duplication

**Decision**: Unify into a single `dw_iteration(subprob_struct*, faux_globals*,
master_data*, int phase, int first_run, ...)` function.

**Rationale**: The two functions (~250 lines each) share:
- Service-queue consumption loop (identical)
- Column addition to master LP (near-identical with one extra `obj_names` step
  in Phase I)
- `COMMAND_GO` signaling (identical)
- `glp_simplex` call + return-code logging (identical)
- Simplex-iteration monotonicity check (identical)
- Dual-variable extraction into `md->row_duals` (identical)
- `pthread_cond_broadcast` to advance iteration (identical)

The existing divergence is already a bug: Phase I uses `> TOLERANCE` for the
column-acceptance test; Phase II uses `> 0.0 + TOLERANCE`. These are
numerically equivalent but make it appear they may be different tests.

Key differences (must be parametrised or branched):
- Phase I receives `obj_names`, `obj_coefs`, `obj_count` for later objective
  restoration; Phase II does not.
- Phase I has a `first_run` branch that initialises `y_accumulators` and
  corrects the sign of `y_`-columns.
- Phase II has convergence detection (comparing `simplex_iterations` to
  `LPX_K_ITCNT`) with the special `rc = -1` case.

**Alternatives considered**:
- Leaving the functions separate and extracting only the shared loops into a
  helper (`service_queue_consume`, `broadcast_next_iteration`): viable and
  lower-risk, but leaves the column-acceptance logic duplicated.
- Full unification: higher integration effort but cleaner result. Recommended
  as the medium-term target.

---

## Finding 4: `glpk_mutex` Held Over Diagnostic Logging Loop

**Decision**: Guard the logging block in `prepare_column` with the verbosity
check before acquiring the mutex.

**Rationale**: `prepare_column` unconditionally acquires `glpk_mutex`, then
calls `dw_printf(IMPORTANCE_DIAG, ...)` inside the locked region. At runtimes
below `OUTPUT_ALL`, `dw_printf` is a no-op. The mutex is still acquired,
serialising all subproblem threads for a no-op. The fix is:
```c
if (data->globals->verbosity >= OUTPUT_ALL) {
    pthread_mutex_lock(&glpk_mutex);
    // ... logging ...
    pthread_mutex_unlock(&glpk_mutex);
}
```

**Alternatives considered**: Removing the diagnostic logging entirely. Rejected
— the logging is valuable at `OUTPUT_ALL`. The fix is minimal and safe.

---

## Finding 5: `master_mutex` Serialising Read-Only Column Mapping

**Decision**: Replace `pthread_mutex_lock(&master_mutex)` /
`pthread_mutex_unlock(&master_mutex)` in the subproblem column-mapping loop
with `pthread_rwlock_rdlock(&master_rwlock)` / `pthread_rwlock_rdunlock`.

**Rationale**: `original_master_lp` is read-only after the master setup signal
is broadcast. All subproblems are reading the same data concurrently. A
read-write lock allows all subproblem threads to proceed in parallel through
this section, which serialises them today.

**Alternatives considered**:
- Copying data into per-thread structs before launching threads (done for
  some data already via `col_translate[]` / `c[]`). The column-mapping loop
  already does this; the mutex protects reading from `original_master_lp`.
  After the loop, the data lives in per-thread arrays so no ongoing sharing
  exists. A `pthread_rwlock_t` is the minimal correct fix.
- Accepting the serialisation because the mapping loop is one-time at startup:
  valid for small `N`; at large `N` this is a measurable bottleneck.

---

## Finding 6: Dead Commented-Out Code

**Decision**: Remove multi-line blocks of commented-out code in `dw_phases.c`
and `dw_rounding.c`. Retain single-line explanatory comments.

**Rationale**: Git history serves as the authoritative record. Dead code
commented out during development adds noise, makes the files harder to read,
and misleads contributors into thinking the code is intentionally preserved
for future use. Blocks identified:
- `dw_phases.c`: `glp_intopt` per-iteration call, `purge_nonbasics` call,
  `check_aux_vars` call, and the `//if(...)` guards around Phase I branches.
- `dw_rounding.c`: Large `print_zeros` helper that references the NASA ATM
  application's specific variable naming scheme (flight, sector, delay);
  extensive commented `fprintf` blocks in `rounding_thread`.

**`print_zeros` assessment**: This function is application-specific (parses
`(flight,sector,time)` variable names), not called anywhere in the current
codebase, and appears to be dead since the original commit. It should be
removed. `print_zeros_simple` is currently used by `rounding_thread`.

**Alternatives considered**: Moving `print_zeros` to a separate file or an
`examples/` directory. Rejected — it is application-specific, not general, and
its functional correctness cannot be verified without the original ATM problem
inputs.

---

## Finding 7: Mixed `printf` / `dw_printf` Usage

**Decision**: Audit all `printf` calls in `dw_phases.c`, `dw_support.c`, and
`dw_rounding.c`. Route diagnostic/informational calls through `dw_printf`
at the appropriate `IMPORTANCE_*` level. Error/warning calls that are
currently `printf` should be routed to `fprintf(stderr, ...)`.

**Rationale**: Direct `printf` calls bypass the verbosity gate, meaning
messages appear even in `--silent` mode. They also go to stdout, mixing
diagnostic output into the solver's structured output stream.

**Exceptions**: The `print_usage` and `print_timing` functions currently use
raw `printf`; these are intended to always print and can remain as-is or be
wrapped in a thin always-on variant.

---

## Finding 8: Global State (`master_lp`, `D`, `signals`, mutexes)

**Decision**: Defer to P4. Migrate `master_lp`, `original_master_lp`, `D`,
and `signals` into `faux_globals`. Group all `pthread_mutex_t` /
`pthread_cond_t` / `sem_t` objects into a new `sync_state` struct embedded
in `faux_globals`.

**Rationale**: This is the highest-impact architectural change. It makes the
solver library-embeddable and eliminates the multiple-instances-in-one-process
limitation. However, `faux_globals` is already threaded through almost all
call sites; the mechanical change is large but the design is clear.

**Alternatives considered**:
- Thread-local storage for `master_lp`: rejected — the master LP is shared
  between the master thread and all subproblem threads by design.
- Leaving globals as-is: acceptable short-term (P4 is explicitly deferred).

**Prerequisite**: P1–P3 must be complete and all tests green before starting P4.

---

## Summary of Decisions

| # | Finding | Action | Risk |
|---|---------|--------|------|
| 1 | `val=NULL` UB in phase 2 | Allocate `val`, fix `ind` leak | Very low |
| 2 | `const char*` leak in rounding | Use stack buffer `char[BUFF_SIZE]` | Very low |
| 3 | phase_1/phase_2 duplication | Unify with `phase` + `first_run` params | Medium |
| 4 | `glpk_mutex` over logging | Guard with verbosity check | Very low |
| 5 | `master_mutex` read serialisation | Use `pthread_rwlock_t` | Low |
| 6 | Dead commented-out code | Delete; use `git log` as history | Very low |
| 7 | Mixed `printf`/`dw_printf` | Route through `dw_printf`/`fprintf(stderr)` | Low |
| 8 | Global solver state | Migrate to `faux_globals` + `sync_state` | High (deferred) |
