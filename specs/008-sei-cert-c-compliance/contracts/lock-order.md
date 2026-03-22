# Lock Acquisition Order — POS51-C Compliance

**Feature**: 008-sei-cert-c-compliance  
**Date**: 2026-03-21  
**CERT Rule**: [POS51-C](https://wiki.sei.cmu.edu/confluence/display/c/POS51-C.+Avoid+deadlock+with+POSIX+threads+by+using+correct+locking+order)

---

## Total Acquisition Order

The following total order is established over all synchronization primitives in `src/dw_*.c`. When any thread must acquire more than one primitive, it MUST acquire them in ascending level order. Acquiring a lower-level primitive while holding a higher-level one is **prohibited** (lock-order inversion).

| Level | Primitive | Type | Notes |
|-------|-----------|------|-------|
| 1 | `glpk_mutex` | `pthread_mutex_t` | Innermost guard; held for single GLPK file I/O calls only (`lpx_read_cpxlp`, `lpx_write_cpxlp`). Never held in combination with other mutexes. |
| 2 | `master_lp_ready_mutex` | `pthread_mutex_t` | One-shot startup barrier; also guards `master_lp_ready_cv`. Never nested with other mutexes. |
| 3 | `service_queue_mutex` | `pthread_mutex_t` | Brief critical section for `fg->service_queue[]` head-pointer read/increment. Never nested with other mutexes. |
| 4 | `sub_data_mutex[i]` | `pthread_mutex_t*` (array) | Per-subproblem data guard. May be held when level 5 (`next_iteration_mutex`) is acquired (established pattern). After POS52-C fix, GLPK solves are NOT called while this lock is held. |
| 5 | `next_iteration_mutex` | `pthread_mutex_t` (ERRORCHECK) | Guards `signals->current_iteration` and `next_iteration_cv`. May be acquired while `sub_data_mutex[i]` is held. |
| 6 | `master_mutex` | `pthread_mutex_t` | Guards `md->row_duals[]` dual vector. NEVER acquired while `sub_data_mutex[i]` is held; always used in isolation. |
| — | `reduced_cost_mutex` | `pthread_mutex_t` | Currently unused in main execution path; treated as isolated. |
| — | `fputs_mutex` | `pthread_mutex_t` | Guards thread-safe output; always used in isolation; never nested. |

---

## Multi-Lock Acquisition Site Audit

All code paths that acquire more than one synchronization primitive in proximity are listed below, with their ordering verified against the table above.

### dw_phases.c — phase_1_iteration (lines 87–178)

```
service_queue_mutex (lock :87)
service_queue_mutex (unlock :91)
-- RELEASED before next acquire

sub_data_mutex[index] (lock :97)
next_iteration_mutex (lock :117)        ← Level 4 → Level 5 ✅ CORRECT ORDER
next_iteration_mutex (unlock :124)
sub_data_mutex[index] (unlock :178)
```

**Verdict**: ✅ PASS — acquires sub (4) then next (5), releases in reverse order.

### dw_phases.c — phase_2_iteration (lines 334–533)

```
service_queue_mutex (lock :334)
service_queue_mutex (unlock :338)
-- RELEASED before next acquire

sub_data_mutex[index] (lock :343)
next_iteration_mutex (lock :360)        ← Level 4 → Level 5 ✅ CORRECT ORDER
next_iteration_mutex (unlock :367)
sub_data_mutex[index] (unlock :423)

-- back to top of while loop; above repeated per client ...

for loop over all clients:
  sub_data_mutex[i] (lock :521)
  sub_data_mutex[i] (unlock :523)

master_mutex (lock :511)               ← Level 6, no other lock held ✅
master_mutex (unlock :515)

next_iteration_mutex (lock :530)       ← Level 5, no other lock held ✅
next_iteration_mutex (broadcast :532)
next_iteration_mutex (unlock :533)
```

**Verdict**: ✅ PASS — all nested acquires follow level order; sequential (non-nested) uses of master_mutex and next_iteration_mutex are at correct levels.

### dw_main.c — master thread, per-iteration body (lines ~479–647)

```
sub_data_mutex[j] (lock :479)
sub_data_mutex[j] (unlock :481)         ← sequential per j in loop

sub_data_mutex[j] (lock :541)
sub_data_mutex[j] (unlock :543)         ← sequential per j in loop

sub_data_mutex[i] (lock :640)
sub_data_mutex[i] (unlock :642)         ← RELEASED before next acquire ✅

next_iteration_mutex (lock :644)        ← Level 5, no sub_data lock held ✅
next_iteration_mutex (broadcast :646)
next_iteration_mutex (unlock :647)
```

**Verdict**: ✅ PASS — `sub_data_mutex[i]` is fully released BEFORE `next_iteration_mutex` is acquired (sequential, not nested).

### dw_subprob.c — subproblem thread worker (lines ~127–443)

```
master_lp_ready_mutex (lock :127)       ← Level 2, isolated ✅
pthread_cond_wait (:135)
master_lp_ready_mutex (unlock :139)

master_mutex (lock :179)                ← Level 6, isolated ✅
master_mutex (unlock :234)

sub_data_mutex[id] (lock :258)
  [POS52-C fix: unlock before LP solve]
sub_data_mutex[id] (unlock :301)        ← released before glp_simplex ✅
  glp_simplex / glp_intopt run WITHOUT any lock held ✅
sub_data_mutex[id] (re-lock :311)       ← re-acquired to write results ✅
sub_data_mutex[id] (unlock :348)        ← BEFORE sem_post

service_queue_mutex (lock :373)         ← Level 3, sub_data released ✅
sem_post :379/381
service_queue_mutex (unlock :383)

next_iteration_mutex (lock :386)        ← Level 5, service_queue released ✅
pthread_cond_wait (:388)
next_iteration_mutex (unlock :390)

sub_data_mutex[my_data->my_id] (lock :394)   ← Level 4, alone ✅
sub_data_mutex[my_data->my_id] (unlock :401)

glpk_mutex (lock :440)                 ← Level 1, alone ✅
glpk_mutex (unlock :448)
```

**Verdict**: ✅ PASS — No pairs violate the level order. `master_mutex` and `sub_data_mutex` are used sequentially (never nested together).

---

## Inversion Analysis Summary

No lock-order inversions found in the current codebase.

| Pair | Observed Order | Possible Inversion? |
|------|---------------|-------------------|
| sub_data_mutex → next_iteration_mutex | Level 4 before 5 | ✅ Consistent; no reverse path |
| glpk_mutex → anything | Level 1 first when combined | ✅ Never combined in same scope |
| master_mutex → sub_data_mutex | Never nested | ✅ Always sequential |
| service_queue_mutex → sub_data_mutex | Never nested | ✅ Always sequential |

---

## ThreadSanitizer Verification

After all changes are implemented, the test suite MUST be run with TSan:

```bash
./configure CFLAGS="-fsanitize=thread -g -O1" LDFLAGS="-fsanitize=thread"
make clean && make
./tests/dw-tests.sh
```

Expected result: zero data-race or lock-order warnings for `src/dw_*.c` code.
