# pthread/sem Call-Site Audit — POS54-C Compliance

**Feature**: 008-sei-cert-c-compliance  
**Date**: 2026-03-21  
**CERT Rule**: [POS54-C](https://wiki.sei.cmu.edu/confluence/display/c/POS54-C.+Detect+and+handle+POSIX+library+errors)

---

## Rule Summary

POS54-C requires that every call to a POSIX library function that can fail must have its return value checked. Failures at mutex lock, thread creation, or semaphore operations should be treated as fatal in dwsolver — silent failures corrupt serialization invariants.

**Important — two distinct error conventions**:
- `pthread_*` functions return **0 on success, the error code directly on failure** (`strerror(rc)` gives the message). Use `DW_PTHREAD_CHECK(rc, name)`.
- `sem_*` functions return **0 on success, -1 on failure** with the error code in `errno` (`strerror(errno)` gives the message). Use `DW_SEM_CHECK(ret, name)` — do NOT use `DW_PTHREAD_CHECK` for `sem_*` calls, as `strerror(-1)` produces a useless "Unknown error" message.

### Annotation exceptions (documented "always-succeeds" calls)

Per POSIX specification, the following are unconditionally safe to leave unchecked **if annotated with a comment**:
- `pthread_mutex_unlock` on a mutex the calling thread owns (POSIX: "shall not fail" on owned non-robustly-locked mutex)
- `pthread_cond_broadcast` (POSIX: "shall not fail")
- `pthread_attr_getstacksize` after `pthread_attr_init`

All other calls require either `DW_PTHREAD_CHECK(return_value, "name")` (for `pthread_*`) or `DW_SEM_CHECK(return_value, "name")` (for `sem_*`), or an explanatory comment.

---

## Call-Site Table

### dw_support.c — `init_pthread_data()`

| Line | Call | Current Status | Required Action |
|------|------|---------------|----------------|
| 190 | `pthread_mutexattr_init(my_mutex_attr)` | Unchecked | `DW_PTHREAD_CHECK` |
| 195 | `pthread_mutex_init(&sub_data_mutex[i], NULL)` (loop) | Unchecked | `DW_PTHREAD_CHECK` |
| 196 | `pthread_mutex_init(&master_lp_ready_mutex, NULL)` | Unchecked | `DW_PTHREAD_CHECK` |
| 197 | `pthread_cond_init(&master_lp_ready_cv, NULL)` | Unchecked | `DW_PTHREAD_CHECK` |
| 198 | `pthread_mutex_init(&service_queue_mutex, NULL)` | Unchecked | `DW_PTHREAD_CHECK` |
| 199 | `pthread_mutex_init(&next_iteration_mutex, my_mutex_attr)` | Unchecked | `DW_PTHREAD_CHECK` |
| 200 | `pthread_mutex_init(&master_mutex, NULL)` | Unchecked | `DW_PTHREAD_CHECK` |
| 201 | `pthread_mutex_init(&reduced_cost_mutex, NULL)` | Unchecked | `DW_PTHREAD_CHECK` |
| 202 | `pthread_mutex_init(&glpk_mutex, NULL)` | Unchecked | `DW_PTHREAD_CHECK` |
| 203 | `pthread_mutex_init(&fputs_mutex, NULL)` | Unchecked | `DW_PTHREAD_CHECK` |
| 204 | `pthread_cond_init(&next_iteration_cv, NULL)` | Unchecked | `DW_PTHREAD_CHECK` |
| 205 | `pthread_attr_init(&attr)` | Unchecked | `DW_PTHREAD_CHECK` |
| 209 | `sem_init(&customers, 0, 0)` | Unchecked | `DW_SEM_CHECK` (sem_* uses errno, not rc) |
| 212 | `pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE)` | Unchecked | `DW_PTHREAD_CHECK` |
| 215 | `pthread_attr_getstacksize(&attr, &stacksize)` | Unchecked | Annotate: always succeeds |
| 217 | `rc = pthread_attr_setstacksize(&attr, stacksize)` | Assigned, NOT tested | Add `DW_PTHREAD_CHECK(rc, ...)` |
| 218 | `pthread_attr_getstacksize(&attr, &stacksize)` | Unchecked | Annotate: always succeeds |

### dw_main.c

| Line | Call | Current Status | Required Action |
|------|------|---------------|----------------|
| 204 | `rc = pthread_create(&threads[i], ...)` | **Checked ✅** | No change |
| 221 | `pthread_mutex_lock(&glpk_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 223 | `pthread_mutex_unlock(&glpk_mutex)` | Unchecked | Annotate: always succeeds |
| 258 | `pthread_mutex_lock(&master_lp_ready_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 260 | `pthread_cond_broadcast(&master_lp_ready_cv)` | Unchecked | Annotate: always succeeds |
| 261 | `pthread_mutex_unlock(&master_lp_ready_mutex)` | Unchecked | Annotate: always succeeds |
| 458 | `pthread_mutex_lock(&glpk_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 460 | `pthread_mutex_unlock(&glpk_mutex)` | Unchecked | Annotate: always succeeds |
| 479 | `pthread_mutex_lock(&sub_data_mutex[j])` | Unchecked | `DW_PTHREAD_CHECK` |
| 481 | `pthread_mutex_unlock(&sub_data_mutex[j])` | Unchecked | Annotate: always succeeds |
| 541 | `pthread_mutex_lock(&sub_data_mutex[j])` | Unchecked | `DW_PTHREAD_CHECK` |
| 543 | `pthread_mutex_unlock(&sub_data_mutex[j])` | Unchecked | Annotate: always succeeds |
| 640 | `pthread_mutex_lock(&sub_data_mutex[i])` | Unchecked | `DW_PTHREAD_CHECK` |
| 642 | `pthread_mutex_unlock(&sub_data_mutex[i])` | Unchecked | Annotate: always succeeds |
| 644 | `pthread_mutex_lock(&next_iteration_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 646 | `pthread_cond_broadcast(&next_iteration_cv)` | Unchecked | Annotate: always succeeds |
| 647 | `pthread_mutex_unlock(&next_iteration_mutex)` | Unchecked | Annotate: always succeeds |
| 706 | `rc = pthread_join(threads[i], &status)` | **Checked ✅** | No change |
| 826 | `pthread_mutex_lock(&fputs_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 828 | `pthread_mutex_unlock(&fputs_mutex)` | Unchecked | Annotate: always succeeds |

### dw_phases.c

| Line | Call | Current Status | Required Action |
|------|------|---------------|----------------|
| 83 | `sem_wait(customers)` (named semaphore path) | **Missing from audit** | `DW_SEM_CHECK` (sem_* uses errno, not rc) |
| 85 | `sem_wait(&customers)` (anonymous semaphore path) | **Missing from audit** | `DW_SEM_CHECK` (sem_* uses errno, not rc) |
| 87 | `pthread_mutex_lock(&service_queue_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 91 | `pthread_mutex_unlock(&service_queue_mutex)` | Unchecked | Annotate: always succeeds |
| 97 | `pthread_mutex_lock(&sub_data_mutex[index])` | Unchecked | `DW_PTHREAD_CHECK` |
| 117 | `pthread_mutex_lock(&next_iteration_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 124 | `pthread_mutex_unlock(&next_iteration_mutex)` | Unchecked | Annotate: always succeeds |
| 178 | `pthread_mutex_unlock(&sub_data_mutex[index])` | Unchecked | Annotate: always succeeds |
| 271 | `pthread_mutex_lock(&sub_data_mutex[i])` (loop) | Unchecked | `DW_PTHREAD_CHECK` |
| 273 | `pthread_mutex_unlock(&sub_data_mutex[i])` (loop) | Unchecked | Annotate: always succeeds |
| 280 | `pthread_mutex_lock(&next_iteration_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 282 | `pthread_cond_broadcast(&next_iteration_cv)` | Unchecked | Annotate: always succeeds |
| 283 | `pthread_mutex_unlock(&next_iteration_mutex)` | Unchecked | Annotate: always succeeds |
| 330 | `sem_wait(customers)` (named semaphore path) | **Missing from audit** | `DW_SEM_CHECK` (sem_* uses errno, not rc) |
| 332 | `sem_wait(&customers)` (anonymous semaphore path) | **Missing from audit** | `DW_SEM_CHECK` (sem_* uses errno, not rc) |
| 334 | `pthread_mutex_lock(&service_queue_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 338 | `pthread_mutex_unlock(&service_queue_mutex)` | Unchecked | Annotate: always succeeds |
| 343 | `pthread_mutex_lock(&sub_data_mutex[index])` | Unchecked | `DW_PTHREAD_CHECK` |
| 360 | `pthread_mutex_lock(&next_iteration_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 367 | `pthread_mutex_unlock(&next_iteration_mutex)` | Unchecked | Annotate: always succeeds |
| 423 | `pthread_mutex_unlock(&sub_data_mutex[index])` | Unchecked | Annotate: always succeeds |
| 511 | `pthread_mutex_lock(&master_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 515 | `pthread_mutex_unlock(&master_mutex)` | Unchecked | Annotate: always succeeds |
| 521 | `pthread_mutex_lock(&sub_data_mutex[i])` (loop) | Unchecked | `DW_PTHREAD_CHECK` |
| 523 | `pthread_mutex_unlock(&sub_data_mutex[i])` (loop) | Unchecked | Annotate: always succeeds |
| 530 | `pthread_mutex_lock(&next_iteration_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 532 | `pthread_cond_broadcast(&next_iteration_cv)` | Unchecked | Annotate: always succeeds |
| 533 | `pthread_mutex_unlock(&next_iteration_mutex)` | Unchecked | Annotate: always succeeds |

### dw_subprob.c

| Line | Call | Current Status | Required Action |
|------|------|---------------|----------------|
| 101 | `pthread_mutex_lock(&glpk_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 103 | `pthread_mutex_unlock(&glpk_mutex)` | Unchecked | Annotate: always succeeds |
| 127 | `pthread_mutex_lock(&master_lp_ready_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 131 | `if` guard on `pthread_cond_wait` | Bug: `if` not `while` | Change to `while` (POS53-C, D3) |
| 134 | `pthread_cond_wait(&master_lp_ready_cv, ...)` | Return unchecked | Annotate: spurious returns handled by `while` guard |
| 141 | `pthread_mutex_unlock(&master_lp_ready_mutex)` | Unchecked | Annotate: always succeeds |
| 181 | `pthread_mutex_lock(&master_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 236 | `pthread_mutex_unlock(&master_mutex)` | Unchecked | Annotate: always succeeds |
| 260 | `pthread_mutex_lock(&sub_data_mutex[id])` | Unchecked | `DW_PTHREAD_CHECK` |
| 343 | `pthread_mutex_unlock(&sub_data_mutex[id])` | Unchecked | Annotate: always succeeds — **NOTE: this unlock moves to before line 303 after POS52-C fix** |
| 368 | `pthread_mutex_lock(&service_queue_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 374 | `sem_post(customers)` (named semaphore variant) | Unchecked | `DW_SEM_CHECK` (sem_* uses errno, not rc) |
| 376 | `sem_post(&customers)` (anonymous semaphore variant) | Unchecked | `DW_SEM_CHECK` (sem_* uses errno, not rc) |
| 378 | `pthread_mutex_unlock(&service_queue_mutex)` | Unchecked | Annotate: always succeeds |
| 381 | `pthread_mutex_lock(&next_iteration_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 383 | `pthread_cond_wait(&next_iteration_cv, ...)` | `while` guard ✅ (fix branch) | Return value: annotate |
| 385 | `pthread_mutex_unlock(&next_iteration_mutex)` | Unchecked | Annotate: always succeeds |
| 389 | `pthread_mutex_lock(&sub_data_mutex[my_data->my_id])` | Unchecked | `DW_PTHREAD_CHECK` |
| 396 | `pthread_mutex_unlock(&sub_data_mutex[my_data->my_id])` | Unchecked | Annotate: always succeeds |
| 435 | `pthread_mutex_lock(&glpk_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 443 | `pthread_mutex_unlock(&glpk_mutex)` | Unchecked | Annotate: always succeeds |

### dw_rounding.c

| Line | Call | Current Status | Required Action |
|------|------|---------------|----------------|
| 175 | `rc = pthread_create(&threads[i], ...)` | **Checked ✅** | No change |
| 181 | `rc = pthread_create(&threads[i], ...)` | **Checked ✅** | No change |
| 199 | `rc = pthread_join(threads[i], &status)` | **Checked ✅** | No change |
| 518 | `pthread_mutex_lock(&glpk_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 524 | `pthread_mutex_unlock(&glpk_mutex)` | Unchecked | Annotate: always succeeds |
| 549 | `pthread_mutex_lock(&master_lp_ready_mutex)` | Unchecked | `DW_PTHREAD_CHECK` |
| 553 | `pthread_mutex_unlock(&master_lp_ready_mutex)` | Unchecked | Annotate: always succeeds |
| 591 | `pthread_mutex_lock(&glpk_mutex)` (loop) | Unchecked | `DW_PTHREAD_CHECK` |
| 600 | `pthread_mutex_unlock(&glpk_mutex)` (loop) | Unchecked | Annotate: always succeeds |

---

## Summary Statistics

| Category | Count |
|----------|-------|
| Already correctly checked | 5 |
| Require `DW_PTHREAD_CHECK` | ~24 |
| Annotate as "always succeeds" | ~28 |
| Bug fix (`if`→`while`) | 1 (dw_subprob.c:131) |
| **Total call sites** | **~57** |
