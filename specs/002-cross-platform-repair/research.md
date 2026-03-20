# Research: Cross-Platform Build Repair & Safety Hardening

**Feature**: `002-cross-platform-repair`
**Date**: 2026-03-19
**Sources**: Direct source inspection of `src/dw.h` (lines 108-168), `src/Makefile.am` (lines 7-122), `specs/001-as-built/spec.md` (KD-001 through KD-009)

---

## Finding 1 — Root Cause of KD-001 (Linux Link Failure)

**Question**: Why does the build work on macOS/Clang but fail on Linux/GCC?

**Investigation**:

The C99 standard (ISO/IEC 9899:1999 §6.9 ¶5) specifies that a *translation unit*
may have any number of *tentative definitions* of an external-linkage object, but
only one *definition* (with or without an initializer) across all translation units
in a program. A bare declaration such as `pthread_mutex_t master_lp_ready_mutex;`
at file scope in a header is, under strict C99, a tentative *definition* — it
reserves storage for the object.

When `dw.h` is included by five `.c` files, the linker sees five definitions of
each global. GCC, which implements strict C99 by default since GCC 10 (and `-fno-common`
as default since GCC 10), rejects this with:

```
multiple definition of `master_lp_ready_mutex`
```

Clang on macOS applies the traditional "common symbol" merging rule (an extension
from UNIX linkers) which silently deduplicates the definitions. This is why the
code compiled on macOS but not Linux.

**Verified**: The five files that `#include "dw.h"` are exactly:

```
src/dw_main.c      (line 87)
src/dw_support.c   (line 40)
src/dw_phases.c    (line 47)
src/dw_subprob.c   (line 45)
src/dw_rounding.c  (line 54)
```

**Complete symbol list defined (not declared) in `dw.h` at lines 109-167**:

```c
pthread_attr_t   attr;                    // line 109
pthread_mutex_t  master_lp_ready_mutex;   // line 110
pthread_cond_t   master_lp_ready_cv;      // line 111
pthread_mutex_t  service_queue_mutex;     // line 112
pthread_mutex_t  next_iteration_mutex;    // line 113
pthread_cond_t   next_iteration_cv;       // line 114
pthread_mutex_t  master_mutex;            // line 115
pthread_mutex_t  reduced_cost_mutex;      // line 116
pthread_mutex_t  glpk_mutex;              // line 117
pthread_mutex_t  fputs_mutex;             // line 118
pthread_mutex_t* sub_data_mutex;          // line 119
sem_t* customers;   (or sem_t)            // lines 121-125 (#ifdef USE_NAMED_SEMAPHORES)
glp_prob* original_master_lp;            // line 127
glp_prob* master_lp;                     // line 128
glp_iocp* parm;                          // line 129
glp_smcp* simplex_control_params;        // line 130
D_matrix* D;                             // line 144 (after typedef)
signal_data* signals;                    // line 167 (after typedef)
```

Total: **18 symbols** (17 distinct names; `customers` is conditional on one
`#ifdef` branch).

**Decision**: 
- Create `src/dw_globals.c` — defines all 18 symbols.
- Edit `src/dw.h` — convert all 18 to `extern` declarations.
- Edit `src/Makefile.am` — add `dw_globals.c` to `dwsolver_SOURCES`.

**Alternatives rejected**:

| Option | Why rejected |
|--------|-------------|
| `-fcommon` linker flag in `Makefile.am` | Deprecated (removed as GCC default in GCC 10); not available on all linkers; doesn't fix the underlying bug; non-portable |
| `#define DEFINE_GLOBALS` pattern (one TU includes header differently) | Requires include-order discipline; fragile with Makefile parallelism; still non-standard C99 |
| Move all state into a passed struct | Invasive refactor of all function signatures — a rewrite, not a repair; violates Principle IV |

---

## Finding 2 — `sprintf` Call Sites (KD-004)

**All `sprintf` calls in DW-authored files** (verified by grep):

| File | Line | Buffer | Format |
|------|------|--------|--------|
| `dw_main.c` | 275 | `local_buffer[BUFF_SIZE]` | `"sub%d_convexity"` — max 24 chars |
| `dw_main.c` | 496 | `local_buffer[BUFF_SIZE]` | `"phase1_step_0.cpxlp"` — 20 chars |
| `dw_main.c` | 512 | `local_buffer[BUFF_SIZE]` | `"phase1_step_%d.cpxlp"` — max ~26 chars |
| `dw_main.c` | 611 | `local_buffer[BUFF_SIZE]` | `"master_step_%d.cpxlp"` — max ~27 chars |
| `dw_support.c` | 609 | `filename[BUFF_SIZE]` | `"basis_iteration_%d"` — max ~28 chars |

`BUFF_SIZE` = 256 (defined in `dw.h`). All outputs are well within 256 bytes given
plausible iteration counts (≤ 999 → 3-digit int → max 28 chars). No truncation
risk in practice, but the replacement is a zero-risk correctness/portability improvement
required for `-D_FORTIFY_SOURCE=2` builds common on production Linux.

**Decision**: Replace each with `snprintf(buf, BUFF_SIZE, ...)`. No functional change.

---

## Finding 3 — `malloc` Return Value Audit (KD-005)

**Critical allocation sites** (failure → undefined behaviour in same call or thread):

| File | Line | Variable | Consequence of NULL |
|------|------|----------|---------------------|
| `dw_main.c` | 110 | `local_buffer` | Immediate UB (used as printf target line 113) |
| `dw_main.c` | 118 | `md` | Immediate UB (field set line 121) |
| `dw_main.c` | 119 | `globals` | Immediate UB (field set line 124) |
| `dw_main.c` | 154 | `threads` | NULL-deref in thread-create loop |
| `dw_main.c` | 155 | `sub_data` | NULL-deref in same loop |
| `dw_support.c` | 99 | `D` | Immediate UB (field malloc lines 104-107) |
| `dw_support.c` | 173 | `signals` | NULL-deref in all thread signal waits |

Non-critical sites (string buffers, itermediate working storage) are lower priority
and not addressed in this spec.

**Decision**: Introduce a static inline helper `dw_oom_abort(void* ptr, const char* ctx)`
in `dw.h`. For each of the 7 critical sites, add a one-line null-check call immediately
after `malloc`. The exit path prints to `stderr` and calls `exit(EXIT_FAILURE)`.

---

## Finding 4 — Extended Test Coverage (KD-009)

**`four_sea` example**:
- Located at `examples/four_sea/`
- Guide file: `four_sea.dw`; master: `four_sea.cplex`; 4 subproblem LP files
- 1760-line solution output; verified by running `src/dwsolver`:
  objective value = total delay = **12** (deterministic)
- Variables: non-deterministic ordering due to thread scheduling
- Test strategy: run example, read last meaningful line of `relaxed_solution`,
  assert numeric value == 12

**`book_dantzig` example**:
- Located at `examples/book_dantzig/`
- The example is noted as "non-deterministic" in the as-built spec (issue KD comment)
- Need to verify the objective value from a reference run
- Test strategy: same — assert objective value from `relaxed_solution`

**Decision**: Extend `tests/dw-tests.sh` with a new `test_objective` function that:
1. Runs solver with timeout
2. Greps the `relaxed_solution` file for the objective value line
3. Compares to expected constant (pass/fail assertion)

---

## Summary of Decisions

| # | Decision | Risk | Confidence |
|---|----------|------|-----------|
| 1 | Introduce `dw_globals.c`; `extern` in `dw.h` | Very low — pure linker-level fix | High |
| 2 | `sprintf` → `snprintf` with `BUFF_SIZE` | Very low — no format string changes | High |
| 3 | Null-check 7 critical `malloc` sites with exit | Very low — only adds a check | High |
| 4 | Extend test suite with objective-value assertions | Low — read-only solution parsing | High |

All NEEDS CLARIFICATION items from Technical Context are now resolved. The implementation
can proceed directly to task generation.
