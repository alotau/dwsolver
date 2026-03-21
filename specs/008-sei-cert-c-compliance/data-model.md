# Data Model: SEI CERT C Compliance

**Feature**: 008-sei-cert-c-compliance  
**Date**: 2026-03-21

This document captures the synchronization primitive model (the "entities" of this feature) and the error-handling patterns introduced to achieve CERT compliance. Since this is a C code hardening effort rather than a data-driven application, the "data model" represents the synchronization architecture and the new code-level constructs.

---

## Synchronization Primitives (Entities)

### Five Mutex Objects

| Variable | Type | Declared In | Protects |
|----------|------|-------------|---------|
| `glpk_mutex` | `pthread_mutex_t` | `dw_globals.c` | GLPK file I/O calls (`lpx_read_cpxlp`, `lpx_write_cpxlp`) — not needed around `glp_simplex` |
| `master_lp_ready_mutex` | `pthread_mutex_t` | `dw_globals.c` | `signals->master_lp_ready` flag + `master_lp_ready_cv` |
| `service_queue_mutex` | `pthread_mutex_t` | `dw_globals.c` | `fg->service_queue[]`, `fg->head_service_queue` |
| `sub_data_mutex[i]` | `pthread_mutex_t*` (array, length `num_clients`) | `dw_globals.c` | `sub_data[i].*` fields (obj, r, unbounded, current_solution, phase_one, command, len, ind, val) |
| `next_iteration_mutex` | `pthread_mutex_t` (ERRORCHECK attr) | `dw_globals.c` | `signals->current_iteration` + `next_iteration_cv` |
| `master_mutex` | `pthread_mutex_t` | `dw_globals.c` | `md->row_duals[]` (master dual vector) — written by master, read by subproblem threads |
| `reduced_cost_mutex` | `pthread_mutex_t` | `dw_globals.c` | Reduced cost vector (currently unused in main loop) |
| `fputs_mutex` | `pthread_mutex_t` | `dw_globals.c` | Thread-safe `fputs` / `dw_printf` output |

### Two Condition Variables

| Variable | Type | Associated Mutex | Signal Condition |
|----------|------|-----------------|-----------------|
| `master_lp_ready_cv` | `pthread_cond_t` | `master_lp_ready_mutex` | `signals->master_lp_ready == 1` — set once at master init |
| `next_iteration_cv` | `pthread_cond_t` | `next_iteration_mutex` | `signals->current_iteration > local_iteration` — incremented every DW iteration |

### One Semaphore

| Variable | Type | Semantics |
|----------|------|-----------|
| `customers` | `sem_t` (anonymous) or named semaphore | Counts subproblems ready to report; posted by each subproblem after solve; waited by master in phase_1_iteration / phase_2_iteration |

### Thread Array

| Variable | Type | Lifecycle |
|----------|------|-----------|
| `threads[i]` | `pthread_t*` | Created in `dw_main.c` (or `dw_rounding.c`); joined after all iterations complete |

---

## Lock Acquisition Total Order

Defined by FR-003. Primitives must be acquired in the order below (lower level acquired before higher):

| Level | Primitive | Held Duration | Notes |
|-------|-----------|--------------|-------|
| 1 | `glpk_mutex` | ~microseconds | Single GLPK file-I/O call only |
| 2 | `master_lp_ready_mutex` | ~microseconds | Startup only; re-used for wait on `master_lp_ready_cv` |
| 3 | `service_queue_mutex` | ~microseconds | Head-pointer read/increment only |
| 4 | `sub_data_mutex[i]` | 1–100 ms (iteration body) | **After POS52-C fix**: no GLPK solves held here |
| 5 | `next_iteration_mutex` | ~microseconds | Broadcast + counter increment only; acquired inside sub_data_mutex in phases |
| 6 | `master_mutex` | ~microseconds | Dual vector copy only; always acquired WITHOUT sub_data_mutex held |

**Rules**:
- Levels 1–3 are never nested with each other.  
- Level 4 (`sub_data_mutex[i]`) may be held when level 5 (`next_iteration_mutex`) is acquired (established pattern; no inversion found).  
- Level 6 (`master_mutex`) is NEVER acquired while level 4 (`sub_data_mutex`) is held.  
- All levels are released BEFORE calling `glp_simplex` or `glp_intopt`.

---

## New Code-Level Constructs (Error-Handling Pattern)

### `DW_PTHREAD_CHECK` — error check helper

**Purpose**: Centralizes the POS54-C error check for fallible pthread/sem calls.  
**Placement**: `src/dw_support.h` (inline helper or macro alongside `dw_oom_abort`).

```c
/* POS54-C: Check return value of fallible pthread/sem primitives.
 * Call with the return value rc and a string literal name.
 * Terminates with stderr message if rc != 0. */
static inline void dw_pthread_check(int rc, const char *name) {
    if (rc != 0) {
        fprintf(stderr, "pthread error in %s: %s\n", name, strerror(rc));
        exit(EXIT_FAILURE);
    }
}
#define DW_PTHREAD_CHECK(rc, name) dw_pthread_check((rc), (name))
```

**Usage pattern**:
```c
/* Before (no check): */
pthread_mutex_lock(&my_mutex);

/* After (POS54-C compliant): */
DW_PTHREAD_CHECK(pthread_mutex_lock(&my_mutex), "pthread_mutex_lock(&my_mutex)");
```

**Always-succeeds annotation pattern** (for unlock and broadcast — per POSIX they cannot fail on owned initialized mutexes):
```c
pthread_mutex_unlock(&my_mutex);  /* always succeeds: unlocking owned mutex */
pthread_cond_broadcast(&my_cv);   /* always succeeds */
```

### POS52-C — Lock release/reacquire around `glp_simplex` in `dw_subprob.c`

**Before** (violation — lines 260–343):
```c
pthread_mutex_lock(&sub_data_mutex[id]);
/* ... set GLPK objective coefficients ... */
ret = glp_simplex(lp, simplex_control_params);   /* VIOLATION: holding sub_data_mutex */
my_data->obj = glp_get_obj_val(lp);
if (my_data->globals->enforce_sub_integrality) {
    glp_intopt(lp, int_parm);                    /* VIOLATION: holding sub_data_mutex */
    my_data->obj = glp_mip_obj_val(lp);
}
/* ... copy results to sub_data[id] ... */
pthread_mutex_unlock(&sub_data_mutex[id]);
```

**After** (POS52-C compliant):
```c
DW_PTHREAD_CHECK(pthread_mutex_lock(&sub_data_mutex[id]), "lock sub_data_mutex");
/* ... set GLPK objective coefficients (reads from sub_data[id]) ... */
pthread_mutex_unlock(&sub_data_mutex[id]);  /* always succeeds: release before solve */

/* POS52-C: glp_simplex must not be called while holding any mutex. */
ret = glp_simplex(lp, simplex_control_params);
double obj_val = glp_get_obj_val(lp);
int unbounded = glp_get_status(lp);
if (my_data->globals->enforce_sub_integrality) {
    glp_intopt(lp, int_parm);
    obj_val = glp_mip_obj_val(lp);
}
glp_create_index(lp);
/* ... compute solution values from private lp ... */

DW_PTHREAD_CHECK(pthread_mutex_lock(&sub_data_mutex[id]), "relock sub_data_mutex");
my_data->obj       = obj_val;
my_data->unbounded = unbounded;
/* ... copy results back to sub_data[id] ... */
pthread_mutex_unlock(&sub_data_mutex[id]);  /* always succeeds */
```

### POS53-C / Spurious wakeup — `if` → `while` at `dw_subprob.c:131`

**Before**:
```c
pthread_mutex_lock(&master_lp_ready_mutex);
if (!signals->master_lp_ready) {
    pthread_cond_wait(&master_lp_ready_cv, &master_lp_ready_mutex);
}
pthread_mutex_unlock(&master_lp_ready_mutex);
```

**After**:
```c
DW_PTHREAD_CHECK(pthread_mutex_lock(&master_lp_ready_mutex), "lock master_lp_ready_mutex");
while (!signals->master_lp_ready) {  /* POS53-C: guard against spurious wakeup */
    pthread_cond_wait(&master_lp_ready_cv, &master_lp_ready_mutex);
}
pthread_mutex_unlock(&master_lp_ready_mutex);  /* always succeeds */
```

---

## State Transitions

### Thread lifecycle

```
pthread_create → [subproblem thread starts]
                      ↓
               wait for master_lp_ready_cv (while guard)
                      ↓
               iteration loop:
                 lock sub_data_mutex[id]
                 read master duals → set GLPK objective
                 UNLOCK sub_data_mutex[id]         ← POS52-C fix
                 glp_simplex(lp)                   ← no lock held
                 RELOCK sub_data_mutex[id]          ← POS52-C fix
                 write results back
                 unlock sub_data_mutex[id]
                      ↓
               post semaphore → signal master
                      ↓
               wait on next_iteration_cv (while guard)  ← existing fix
                      ↓
               [repeat or exit on COMMAND_STOP]
                      ↓
               pthread_exit → pthread_join
```
