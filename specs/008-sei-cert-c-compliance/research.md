# Research: SEI CERT C Compliance — Code Audit Findings

**Feature**: 008-sei-cert-c-compliance  
**Date**: 2026-03-21  
**Method**: Direct source inspection of `src/dw_*.c` via grep and read-file audit.

---

## Decision Log

### D1 — Error-check macro pattern for pthread/sem calls

- **Decision**: Introduce `DW_PTHREAD_CHECK(rc, name)` macro (or equivalent inline helper `dw_pthread_abort`) — analogous to the existing `dw_oom_abort` OOM-abort pattern already in the codebase.
- **Rationale**: ~50 call sites need identical checks. A macro centralizes the error message format and ensures no site is accidentally skipped. Inline helper vs. macro: inline function preferred in C99 to avoid double-evaluation; macro acceptable if `rc` is always a simple `int` variable.
- **Alternatives considered**:
  - Per-site `if (rc != 0) { fprintf(stderr, ...); exit(1); }` — too verbose, too easy to miss one.
  - Wrapper functions (e.g., `dw_mutex_lock`) — heavier refactor; deferred to future cleanup work.
- **Reference**: `dw_support.c` already has the `dw_oom_abort` inline helper as a pattern to follow.

### D2 — POS52-C fix strategy for `glp_simplex` under `sub_data_mutex[id]`

- **Decision**: In `dw_subprob.c`, release `sub_data_mutex[id]` after writing objective coefficients to the GLPK LP object (before line 303), call `glp_simplex` / `glp_intopt` without the lock, then re-acquire `sub_data_mutex[id]` to write back results (`my_data->obj`, `my_data->unbounded`, etc.).
- **Rationale**: The mutex protects `sub_data[id]` struct fields. The GLPK LP object `lp` is thread-private (allocated per subproblem thread, never shared), so its manipulation does not require the mutex. The solve itself is therefore safe to perform without the lock.
- **Correctness verification**: All fields read from `my_data` to set up the GLPK objective (row_duals, D_vals, c) are copied/computed into GLPK before releasing. Results that the master reads (`sub_data[id].obj`, `sub_data[id].unbounded`, `sub_data[id].current_solution`) are written back only after re-acquiring.
- **Alternatives considered**:
  - Use a separate mutex for the GLPK LP object — unnecessary since `lp` is thread-private.
  - Use `pthread_mutex_trylock` / polling — adds busy-wait complexity with no benefit.

### D3 — Spurious wakeup bug at `dw_subprob.c:134` (second `if` → `while`)

- **Decision**: Change `if (!signals->master_lp_ready)` to `while (!signals->master_lp_ready)` at line 131.
- **Rationale**: `pthread_cond_wait` may return spuriously. The flag `signals->master_lp_ready` is set to `1` only by the master thread under `master_lp_ready_mutex` before broadcasting. A spurious wakeup that exits the `if` guard would proceed with an unready master LP — undefined reads from master data structures follow immediately. This is the same class of bug as was fixed in `fix/spurious-wakeup-cond-wait` (line 382).
- **Note**: This fix was not included in the original `fix/spurious-wakeup-cond-wait` branch. It should be included in this feature.

### D4 — POS51-C lock acquisition order

- **Decision**: Establish the following total order (lower number acquired first):
  1. `glpk_mutex` (innermost GLPK file-I/O guard — always held for the shortest possible time)
  2. `master_lp_ready_mutex` (one-shot startup barrier — acquired before sub_data_mutex in init)
  3. `service_queue_mutex` (queue head/tail pointer update — brief critical section)
  4. `sub_data_mutex[i]` (per-subproblem data — longest-held non-glpk lock)
  5. `next_iteration_mutex` (iteration counter + condvar broadcast — nested inside sub_data_mutex[i] in one code path; requires audit)
  6. `master_mutex` (master dual vector write — acquired in subproblem thread under sub_data_mutex[id])

- **Inversion risk found**: In `dw_phases.c:97→117` and `343→360`, `sub_data_mutex[index]` is acquired first, then `next_iteration_mutex` is acquired inside. This means the order sub_data_mutex → next_iteration_mutex is established. In `dw_main.c:640→647`, the same ordering appears (sub_data_mutex[i] unlocked at 642, then next_iteration_mutex locked at 644 — NOT nested, sequential). No inversion of this pair found.
- **master_mutex nesting**: In `dw_subprob.c:181→236`, `master_mutex` is acquired without holding `sub_data_mutex[id]`. In `dw_phases.c:511→515`, same. No nesting between master_mutex and sub_data_mutex found — they are always used sequentially, not nested.
- **Conclusion**: No actual lock-order inversions found in existing code. Formalization required by FR-003/FR-004.

### D5 — MEM32-C: malloc null-check coverage

- **Decision**: Focus null-check additions on `init_pthread_data()` and initialization paths where a NULL return would cause an immediate NULL-dereference; defer "always-succeeds in practice" small-allocation sites less rigorously per CERT recommendation priority.
- **Rationale**: Approximately 100+ malloc/calloc sites exist. Many allocate fixed-size structs (always <1 MB); OOM on these is fatal regardless and the existing `dw_oom_abort` pattern handles them. The spec (FR-005) requires ALL sites — implementation will mechanically add null checks using `dw_oom_abort` or equivalent. Highest-risk sites first: `my_mutex_attr` at line 187 (dereferenced 3 lines later for `pthread_mutexattr_init`), `sub_data_mutex` array at line 193 (iterated immediately after).
- **Note**: `dw_support.c:384`, `513`, `530`, `549`, `dw_main.c:170`, `dw_rounding.c:155` — all `fopen` calls already have null checks ✅. `dw_main.c:781`, `dw_support.c:612` — also have null checks ✅.

### D6 — INT30/31-C: integer type hazard scope

- **Decision**: Compile with `-Wall -Wextra -Wsign-compare` to enumerate all warnings before fixing; fix only warnings in `src/dw_*.c`, not in vendored GLPK files (which generate their own warnings).
- **Rationale**: No INT30/31 bugs have been observed causing incorrect results. The scope is warning-elimination only. At GLPK API call boundaries, `(int)` explicit casts with size range assertions are the appropriate resolution where `size_t` → `int` conversions occur.
- **Alternatives considered**: Using `size_t` for all loop counters — too invasive; GLPK API takes `int` parameters throughout, so propagating `size_t` would require casts at every GLPK call anyway.

---

## Call-Site Inventories

### POS54-C — Unchecked pthread/sem return values

**Already correctly checked** (have `rc =` assignment + error print):
- `dw_support.c:217` — `rc = pthread_attr_setstacksize(...)` — has `rc` variable but **NOT tested** ❌ (oversight)
- `dw_main.c:204` — `rc = pthread_create(...)` — tested with fprintf + exit ✅
- `dw_main.c:706` — `rc = pthread_join(...)` — tested with fprintf + exit ✅
- `dw_rounding.c:175,181` — `rc = pthread_create(...)` — tested ✅
- `dw_rounding.c:199` — `rc = pthread_join(...)` — tested ✅

**Unchecked — must add DW_PTHREAD_CHECK**:

`dw_support.c` (init_pthread_data):
- Line 190: `pthread_mutexattr_init(my_mutex_attr)`
- Line 195: `pthread_mutex_init(&sub_data_mutex[i], NULL)` (in loop)
- Line 196: `pthread_mutex_init(&master_lp_ready_mutex, NULL)`
- Line 197: `pthread_cond_init(&master_lp_ready_cv, NULL)`
- Line 198: `pthread_mutex_init(&service_queue_mutex, NULL)`
- Line 199: `pthread_mutex_init(&next_iteration_mutex, my_mutex_attr)`
- Line 200: `pthread_mutex_init(&master_mutex, NULL)`
- Line 201: `pthread_mutex_init(&reduced_cost_mutex, NULL)`
- Line 202: `pthread_mutex_init(&glpk_mutex, NULL)`
- Line 203: `pthread_mutex_init(&fputs_mutex, NULL)`
- Line 204: `pthread_cond_init(&next_iteration_cv, NULL)`
- Line 205: `pthread_attr_init(&attr)`
- Line 209: `sem_init(&customers, 0, 0)`
- Line 212: `pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE)`
- Line 215: `pthread_attr_getstacksize(&attr, &stacksize)`
- Line 217: `rc = pthread_attr_setstacksize(...)` — assigned but NOT tested
- Line 218: `pthread_attr_getstacksize(&attr, &stacksize)`

`dw_main.c`:
- Lines 221, 223: `glpk_mutex` lock/unlock
- Lines 258, 260, 261: `master_lp_ready_mutex` lock, `cond_broadcast`, unlock
- Lines 458, 460: `glpk_mutex` lock/unlock
- Lines 479, 481, 541, 543, 640, 642: `sub_data_mutex[j/i]` lock/unlock (multiple call sites)
- Lines 644, 646, 647: `next_iteration_mutex` lock, `cond_broadcast`, unlock
- Lines 826, 828: `fputs_mutex` lock/unlock

`dw_phases.c`:
- Lines 87, 91: `service_queue_mutex` lock/unlock
- Line 97: `sub_data_mutex[index]` lock
- Lines 117, 124: `next_iteration_mutex` lock/unlock
- Line 178: `sub_data_mutex[index]` unlock
- Lines 271, 273: `sub_data_mutex[i]` lock/unlock (in loop)
- Lines 280, 282, 283: `next_iteration_mutex` lock, broadcast, unlock
- Lines 334, 338: `service_queue_mutex` lock/unlock
- Line 343: `sub_data_mutex[index]` lock
- Lines 360, 367: `next_iteration_mutex` lock/unlock
- Line 423: `sub_data_mutex[index]` unlock
- Lines 511, 515: `master_mutex` lock/unlock
- Lines 521, 523: `sub_data_mutex[i]` lock/unlock
- Lines 530, 532, 533: `next_iteration_mutex` lock, broadcast, unlock

`dw_subprob.c`:
- Lines 101, 103: `glpk_mutex` lock/unlock
- Lines 127, 141: `master_lp_ready_mutex` lock/unlock
- Line 134: `pthread_cond_wait(&master_lp_ready_cv, ...)` (also `if`→`while` fix)
- Lines 181, 236: `master_mutex` lock/unlock
- Lines 260, 343: `sub_data_mutex[id]` lock/unlock
- Lines 368, 378: `service_queue_mutex` lock/unlock
- Lines 374/376: `sem_post(&customers)` / `sem_post(customers)` (platform variant)
- Lines 381, 385: `next_iteration_mutex` lock/unlock
- Line 383: `pthread_cond_wait(&next_iteration_cv, ...)` (already `while` after fix branch)
- Lines 389, 396: `sub_data_mutex[my_data->my_id]` lock/unlock
- Lines 435, 443: `glpk_mutex` lock/unlock

`dw_rounding.c`:
- Lines 518, 524: `glpk_mutex` lock/unlock
- Lines 549, 553: `master_lp_ready_mutex` lock/unlock
- Lines 591, 600: `glpk_mutex` lock/unlock (in loop)

**Annotation exceptions** (always-succeeds, per POSIX specification — annotate with comment, no check needed):
- `pthread_mutex_unlock` when called on a mutex the calling thread provably owns — POSIX says this MUST succeed for owned, initialized, non-robustly-locked mutexes. Add `/* always succeeds: unlocking owned mutex */` comment.
- `pthread_cond_broadcast` — POSIX says always succeeds. Add `/* always succeeds */` comment.
- `pthread_attr_getstacksize` — always succeeds after `pthread_attr_init`.

**Checked count**: 5 sites already checked ✅  
**Annotation-only count**: ~20 unlock + broadcast sites  
**New check required**: ~15 init/lock sites

### POS52-C — GLPK calls while holding mutex

| File | Line | Function | Lock Held | Verdict |
|------|------|----------|-----------|---------|
| dw_subprob.c | 112 | glp_simplex(lp) | None (glpk_mutex released line 103) | ✅ PASS |
| dw_subprob.c | 116 | glp_intopt(lp) | None | ✅ PASS |
| dw_subprob.c | 303 | glp_simplex(lp) | **sub_data_mutex[id]** (locked line 260, unlocked line 343) | ❌ FAIL |
| dw_subprob.c | 306 | glp_intopt(lp) | **sub_data_mutex[id]** | ❌ FAIL |
| dw_phases.c | 211 | glp_simplex(master_lp) | None | ✅ PASS |
| dw_phases.c | 435 | glp_simplex(master_lp) | None (all sub_data_mutex released before call) | ✅ PASS |
| dw_main.c | 573 | glp_simplex(master_lp) | None | ✅ PASS |
| dw_main.c | 666 | glp_simplex(master_lp) | None | ✅ PASS |
| dw_main.c | 767 | glp_intopt(master_lp) | None | ✅ PASS |
| dw_main.c | 778 | glp_simplex(master_lp) | None | ✅ PASS |
| dw_rounding.c | — | (none in rounding solve loop) | — | ✅ PASS |

**Violations requiring fix**: 2 (dw_subprob.c lines 303 and 306)

### POS51-C — Multi-lock acquisition sequences

| Location | Lock 1 | Lock 2 | Order |
|----------|--------|--------|-------|
| dw_phases.c:97→117→124→178 | sub_data_mutex[index] | next_iteration_mutex | sub → next |
| dw_phases.c:343→360→367→423 | sub_data_mutex[index] | next_iteration_mutex | sub → next |
| dw_subprob.c:368→374→378 | service_queue_mutex | (sem_post, not a mutex) | service_queue only |
| dw_main.c:640→642→644→647 | sub_data_mutex[i] (unlocked) then next_iteration_mutex | Sequential (not nested) | sub released BEFORE next acquired |
| dw_subprob.c:181→236 | master_mutex | (no other lock nested) | master_mutex alone |
| dw_phases.c:511→515 | master_mutex | (no other lock nested) | master_mutex alone |

**Inversions found**: NONE. The `sub_data_mutex → next_iteration_mutex` ordering is consistent. In `dw_main.c:640–647`, the sequence is unlock-sub, lock-next (sequential, never nested simultaneously).

**Conclusion**: No code changes required for POS51-C correctness. Documentation of the total order (FR-003) is the deliverable.
