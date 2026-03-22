# GLPK Lock-State Audit — POS52-C Compliance

**Feature**: 008-sei-cert-c-compliance  
**Date**: 2026-03-21  
**CERT Rule**: [POS52-C](https://wiki.sei.cmu.edu/confluence/display/c/POS52-C.+Do+not+perform+operations+that+can+block+while+holding+a+POSIX+lock)

---

## Rule Summary

POS52-C prohibits calling functions that may block or take unbounded time while holding a POSIX mutex. `glp_simplex` and `glp_intopt` are LP/MIP solvers that can take seconds to minutes on large problems — they must never be called while any mutex is held.

---

## Complete Call-Site Audit

| # | File | Line | Function | Mutexes Held at Call Time | Verdict | Notes |
|---|------|------|----------|--------------------------|---------|-------|
| 1 | dw_subprob.c | 112 | `glp_simplex(lp, ...)` | None — `glpk_mutex` released at line 103 | ✅ PASS | Initial subproblem solve before DW loop |
| 2 | dw_subprob.c | 116 | `glp_intopt(lp, ...)` | None | ✅ PASS | Initial subproblem integer solve |
| 3 | dw_subprob.c | 303 | `glp_simplex(lp, ...)` | **`sub_data_mutex[id]`** (locked line 260, unlocked line 343) | ❌ **FAIL** | Main iteration re-solve — FIX REQUIRED |
| 4 | dw_subprob.c | 306 | `glp_intopt(lp, ...)` | **`sub_data_mutex[id]`** | ❌ **FAIL** | Main iteration integer solve — FIX REQUIRED |
| 5 | dw_phases.c | 211 | `glp_simplex(master_lp, ...)` | None — all locks released before call | ✅ PASS | Phase 1 master solve |
| 6 | dw_phases.c | 435 | `glp_simplex(master_lp, ...)` | None — `sub_data_mutex[index]` released at 423; no other lock | ✅ PASS | Phase 2 master solve |
| 7 | dw_main.c | 573 | `glp_simplex(master_lp, ...)` | None | ✅ PASS | Phase 1 cleanup re-solve |
| 8 | dw_main.c | 666 | `glp_simplex(master_lp, ...)` | None | ✅ PASS | Post-phase-2 perturbation solve |
| 9 | dw_main.c | 767 | `glp_intopt(master_lp, ...)` | None | ✅ PASS | Integerization solve |
| 10 | dw_main.c | 778 | `glp_simplex(master_lp, ...)` | None | ✅ PASS | Post-integerization re-solve |

**Total violations**: 2 (sites #3 and #4 in `dw_subprob.c`)

---

## Fix Specification for Sites #3 and #4

### Context

In `dw_subprob.c`, the main iteration loop (lines 255–400) works as follows:
1. Acquire `sub_data_mutex[id]` (line 260)
2. Read master dual values from `my_data->md->row_duals` and `my_data->c`
3. Compute reduced-cost vector `d_vector` using BLAS routines
4. Set new objective coefficients on the private `lp` object via `glp_set_obj_coef`
5. **VIOLATION**: Call `glp_simplex(lp, ...)` (line 303) — holding `sub_data_mutex[id]`
6. **VIOLATION**: Call `glp_intopt(lp, ...)` (line 306) — holding `sub_data_mutex[id]`
7. Read results from GLPK, store into `my_data->*` fields
8. Copy solution to `my_data->current_solution[]` (line ~320–340)
9. Release `sub_data_mutex[id]` (line 343)

### Why the fix is safe

The private LP object `lp` is:
- Allocated by the subproblem thread in `init_subproblem_data()` — one per thread
- Never accessed by the master thread or any other subproblem thread
- Therefore, no synchronization is needed to protect `lp` itself

The mutex `sub_data_mutex[id]` protects `sub_data[id].*` struct fields (shared between this thread and the master). Steps 2–4 read from `sub_data[id]` into local variables and into the private `lp`. After step 4, the `lp` contains a complete self-consistent snapshot of the current objective. The mutex may be safely released before step 5.

After the solve completes (steps 5–6), the results are written back to `my_data->*` fields under the re-acquired mutex.

### Required change

Release `sub_data_mutex[id]` after all `glp_set_obj_coef` calls (step 4) and before `glp_simplex` (step 5). Re-acquire `sub_data_mutex[id]` immediately after the solve completes, before writing back to `my_data->*`.

```c
/* Release lock before long-running GLPK solve (POS52-C). */
pthread_mutex_unlock(&sub_data_mutex[id]);  /* always succeeds */

ret = glp_simplex(lp, simplex_control_params);
double result_obj = glp_get_obj_val(lp);
int result_unbounded = glp_get_status(lp);
if (my_data->globals->enforce_sub_integrality) {
    glp_intopt(lp, int_parm);
    result_obj = glp_mip_obj_val(lp);
}
glp_create_index(lp);
/* ... extract my_solution[] from lp ... */

/* Re-acquire lock to write results back to shared struct. */
DW_PTHREAD_CHECK(pthread_mutex_lock(&sub_data_mutex[id]), "relock sub_data_mutex");
my_data->obj       = result_obj;
my_data->unbounded = result_unbounded;
/* ... copy my_solution[] to my_data->condensed_x / .current_solution ... */
```

### Risk assessment

- **Low risk**: `lp` is entirely thread-private; no other thread touches it.
- **Pattern consistency**: `dw_subprob.c:101–116` (initial solve) already follows this pattern correctly — mutex released before `glp_simplex`, results read without lock.
- **Ordering invariant**: The master only reads `sub_data[id].obj` after the subproblem signals via semaphore (`sem_post` at lines 374/376), which occurs later in the function after the re-acquired mutex is released. The happens-before relationship is preserved.

---

## Verification

After implementation:
1. Re-run this audit to confirm column "Mutexes Held at Call Time" shows "None" for all 10 sites.
2. Run ThreadSanitizer build and confirm no data-race warnings on `lp` access.
3. Run full test suite and confirm all examples reproduce expected optima.
