# Quickstart: Implementing SEI CERT C Compliance

**Feature**: 008-sei-cert-c-compliance  
**Date**: 2026-03-21

This guide orients a developer starting implementation work on this feature. It describes the key patterns to apply, the order to work in, and how to verify each step.

---

## Overview

Five CERT rule categories need to be addressed, in priority order:

| Priority | Rule | Files Affected | Description |
|----------|------|---------------|-------------|
| P1 | POS54-C | dw_support.c, dw_main.c, dw_phases.c, dw_subprob.c, dw_rounding.c | Add `DW_PTHREAD_CHECK` to all fallible pthread/sem calls |
| P2 | POS52-C | dw_subprob.c | Release mutex before `glp_simplex`; re-acquire after |
| P2 | POS53-C | dw_subprob.c | Change `if (!signals->master_lp_ready)` → `while` (line 131) |
| P2 | POS51-C | None (docs only) | Lock-order already correct; produce `contracts/lock-order.md` ✅ done |
| P3 | MEM32-C | All dw_*.c | Add null-checks post-malloc using `dw_oom_abort` |
| P3 | FIO01/42-C | dw_support.c, dw_main.c, dw_rounding.c | Verify fopen null-checks; verify all opened files are fclose'd |
| P3 | INT30/31-C | dw_*.c | Fix `-Wall -Wextra -Wsign-compare` warnings |

---

## Step 1 — Add `DW_PTHREAD_CHECK` helper to `dw_support.h` (prerequisite for all P1 work)

Open `src/dw_support.h` (or `src/dw_main.c` where headers are included) and add:

```c
#include <string.h>  /* for strerror */

/* POS54-C: Check return value of fallible pthread/sem primitives. */
static inline void dw_pthread_check(int rc, const char *name) {
    if (rc != 0) {
        fprintf(stderr, "FATAL: pthread error in %s: %s\n", name, strerror(rc));
        exit(EXIT_FAILURE);
    }
}
#define DW_PTHREAD_CHECK(rc, name) dw_pthread_check((rc), (name))
```

The existing `dw_oom_abort` in `dw_support.c` is the model for this pattern.

---

## Step 2 — Apply POS54-C checks (P1)

Work through each file in the order listed in `contracts/pthread-call-sites.md`.

**Template — adding a check**:
```c
/* Before: */
pthread_mutex_lock(&my_mutex);

/* After: */
DW_PTHREAD_CHECK(pthread_mutex_lock(&my_mutex), "pthread_mutex_lock(&my_mutex)");
```

**Template — annotating always-succeeds calls**:
```c
pthread_mutex_unlock(&my_mutex);   /* always succeeds: unlocking owned mutex */
pthread_cond_broadcast(&my_cv);    /* always succeeds */
```

After editing each file, compile (`make`) and verify no new warnings.

---

## Step 3 — Fix spurious wakeup at `dw_subprob.c:131` (P2, POS53-C)

Change line 131 from:
```c
if (!signals->master_lp_ready) {
```
to:
```c
while (!signals->master_lp_ready) {  /* POS53-C: guard against spurious wakeup */
```

---

## Step 4 — Fix POS52-C: release mutex before `glp_simplex` (P2)

See `contracts/glpk-lock-state.md` for the full fix specification.

In `dw_subprob.c`, the main iteration loop (around lines 260–400):
1. Find where all `glp_set_obj_coef` calls complete (just before line 303).
2. Insert `pthread_mutex_unlock(&sub_data_mutex[id]); /* always succeeds: release before solve (POS52-C) */`
3. Extract result variables (`obj_val`, `unbounded`) to local variables.
4. After the solve block, insert `DW_PTHREAD_CHECK(pthread_mutex_lock(&sub_data_mutex[id]), "relock sub_data_mutex");`
5. Write local result variables back to `my_data->obj`, `my_data->unbounded`, etc.

**Test**: Run `tests/dw-tests.sh` and confirm all examples still produce matching output.

---

## Step 5 — Add MEM32-C null-checks (P3)

Use `grep -n 'malloc\|calloc' src/dw_*.c` to enumerate sites. For each site:

```c
/* Before: */
my_ptr = malloc(sizeof(my_type));

/* After: */
my_ptr = malloc(sizeof(my_type));
if (my_ptr == NULL) dw_oom_abort();
```

Priority sites (highest risk — dereference immediately follows allocation):
- `dw_support.c:187` — `my_mutex_attr` (used 2 lines later)
- `dw_support.c:193` — `sub_data_mutex` array (iterated in loop immediately after)
- `dw_main.c:119,121` — `md`, `globals` (used throughout `dw_main`)
- `dw_subprob.c:86` — `simplex_control_params` (passed to `glp_init_smcp` immediately)

---

## Step 6 — Verify FIO01-C / FIO42-C (P3)

Run `grep -n 'fopen' src/dw_*.c` and confirm each call:
1. Immediately null-checks the return: `if (p == NULL) { ... return -1; }`
2. Has a corresponding `fclose(p)` reachable on all exit paths within the enclosing function.

Known-good call sites (already have null-checks): `dw_support.c:384`, `513`, `530`, `549`, `dw_main.c:170`, `dw_rounding.c:155`, `dw_support.c:612`, `dw_main.c:781`.

Verify each also has matching `fclose`.

---

## Step 7 — Fix INT30/31-C warnings (P3)

```bash
gcc -Wall -Wextra -Wsign-compare -c src/dw_*.c -I src/ 2>&1 | grep -v 'glp'
```

Work through each warning. Common patterns:
- Loop counter `int i` compared to `glp_get_num_cols()` which returns `int` — safe as-is, add `(int)` cast to suppress if needed.
- `size_t` vs `int` comparisons around malloc sizes: use `(int)` cast with a comment.
- If a GLPK API takes `int` but a local is `size_t`: cast at the call site only, not throughout.

---

## Verification Checklist

After all steps:

- [ ] `grep -n 'pthread_\|sem_' src/dw_*.c | grep -v '//' | grep -v 'DW_PTHREAD_CHECK\|always succeeds\|= 0'` — returns empty (all call sites either checked or annotated)
- [ ] `grep -n 'glp_simplex\|glp_intopt' src/dw_*.c` — manually verify no mutex is held at each site per `contracts/glpk-lock-state.md`
- [ ] `gcc -Wall -Wextra -Wsign-compare -c src/dw_*.c -I src/ 2>&1 | grep -v 'glp'` — zero warnings
- [ ] `tests/dw-tests.sh` — all tests pass
- [ ] TSan build: `./configure CFLAGS="-fsanitize=thread -g -O1" LDFLAGS="-fsanitize=thread" && make clean && make && ./tests/dw-tests.sh` — zero data-race warnings
