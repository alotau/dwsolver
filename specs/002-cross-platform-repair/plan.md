# Implementation Plan: Cross-Platform Build Repair & Safety Hardening

**Branch**: `002-cross-platform-repair` | **Date**: 2026-03-19 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/002-cross-platform-repair/spec.md`

---

## Summary

The DWSOLVER binary currently builds correctly on macOS/Clang but fails to link on
Linux/GCC due to global variables being *defined* (not merely declared) in the shared
header `src/dw.h`. The fix is to introduce a single new translation unit
`src/dw_globals.c` that holds the one authoritative definition of every global, and
to convert the definitions in `dw.h` to `extern` declarations. This is the minimum
change required to make the software buildable on its primary deployment platform.
Secondary work in this spec covers test coverage for the two currently-uncovered
examples and targeted safety hardening (`sprintf`→`snprintf`, `malloc` null-checks).

---

## Technical Context

**Language/Version**: C99 (POSIX; GCC 9+ and Clang 12+ are primary targets)
**Primary Dependencies**: POSIX Threads (pthreads); GLPK 4.44 (embedded, thread-patched)
**Storage**: File I/O — CPLEX LP format input, text output; no database
**Testing**: Shell example runner `tests/dw-tests.sh`; manual objective-value checks for new tests
**Target Platforms**: macOS (Clang, arm64 + x86_64), Linux (GCC, x86_64); Windows deferred
**Project Type**: CLI solver (library-ready internal structure per Principle V)
**Performance Goals**: No performance budget change; this is a correctness repair
**Constraints**: No new runtime dependencies; must compile under C99 `-pedantic` flags
**Scale/Scope**: ~3 000 lines of DW-authored C across 6 source files (excluding embedded GLPK)

---

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | ✅ PASS | Design preserves all synchronization semantics; test suite runs required before merge |
| II. Thread Safety is Mandatory | ✅ PASS | Moving definitions does **not** change mutex/semaphore init/destroy logic; no init order changes introduced |
| III. Cross-Platform Portability | ✅ PASS — this is the fix | After this change, both macOS and Linux build; `#ifdef USE_NAMED_SEMAPHORES` guard preserved |
| IV. Repair Before Extension | ✅ PASS | This IS the repair. No new features are added |
| V. CLI-First, Library-Ready | ✅ PASS | Globals extraction actually improves library-readiness (defined once, extern everywhere) |

**Gate result**: All five principles satisfied. Proceed to Phase 0.

---

## Project Structure

### Documentation (this feature)

```text
specs/002-cross-platform-repair/
├── plan.md              # This file
├── research.md          # Phase 0 output
├── data-model.md        # Phase 1 output (global state inventory)
├── quickstart.md        # Phase 1 output (build + test instructions)
└── tasks.md             # Phase 2 output (/speckit.tasks command)
```

### Source Code (repository root)

```text
src/
├── dw.h                 # MODIFIED: variable definitions → extern declarations
├── dw_globals.c         # NEW: single translation unit owning all global definitions
├── dw_main.c            # unchanged (includes dw.h; gains definitions from dw_globals.c at link)
├── dw_support.c         # MODIFIED: sprintf→snprintf, malloc null-checks
├── dw_phases.c          # unchanged
├── dw_subprob.c         # unchanged
├── dw_rounding.c        # unchanged
├── dw_blas.c            # unchanged
└── Makefile.am          # MODIFIED: add dw_globals.c to dwsolver_SOURCES

tests/
└── dw-tests.sh          # MODIFIED: add two new objective-value test cases
```

**Structure Decision**: This is a minimal-diff repair to an existing single-project C
codebase. Standard `src/` + `tests/` layout is retained as-is.

---

## Phase 0: Research

> *Findings derived from as-built spec analysis and direct source inspection. Full
> detail in [research.md](research.md).*

### Decision 1 — KD-001 fix pattern

**Affected symbols** (defined, not declared, in `src/dw.h` lines 109–167):
`attr`, `master_lp_ready_mutex`, `master_lp_ready_cv`, `service_queue_mutex`,
`next_iteration_mutex`, `next_iteration_cv`, `master_mutex`, `reduced_cost_mutex`,
`glpk_mutex`, `fputs_mutex`, `sub_data_mutex`, `customers` (conditional),
`original_master_lp`, `master_lp`, `parm`, `simplex_control_params`, `D`, `signals`
— **18 symbols across 5 translation units**.

**Decision**: Introduce `src/dw_globals.c` as the sole definition site; convert
`dw.h` to `extern` declarations; add `dw_globals.c` to `dwsolver_SOURCES` in
`src/Makefile.am`.

**Alternatives rejected**:

| Option | Why rejected |
|--------|-------------|
| `-fcommon` GCC flag | Deprecated (removed as GCC 10 default); doesn't fix the bug |
| `#define DEFINE_GLOBALS` include-order trick | Fragile; non-standard C99 |
| Move all state into a passed struct | Invasive rewrite; violates Principle IV |

### Decision 2 — `sprintf` replacement (KD-004)

5 call sites in `dw_main.c` (lines 275, 496, 512, 611) and `dw_support.c` (line 609).
All write to `char[BUFF_SIZE]` (256 bytes). All outputs are bounded integer
substitutions well within 256 bytes. Replace each with `snprintf(buf, BUFF_SIZE, ...)`.

### Decision 3 — `malloc` null-check strategy (KD-005)

7 critical allocation sites where NULL return causes immediate UB in the same scope:
`local_buffer` (dw_main.c:110), `md` (:118), `globals` (:119), `threads` (:154),
`sub_data` (:155), `D` (dw_support.c:99), `signals` (:173).

**Decision**: Add `dw_oom_abort(void*, const char*)` static inline helper to `dw.h`;
apply at each of the 7 critical sites.

### Decision 4 — Extended test coverage (KD-009)

- `four_sea` final objective: `1.200000e+01` (deterministic; confirmed by run)
- `book_dantzig` final objective: `6.357895e+01` (deterministic; confirmed by 2 runs)
- Both values read from stdout (`#### Master objective value = ...` line), not
  from `relaxed_solution` (which does not contain the objective).

---

## Phase 1: Design Artifacts

### 1.1 Global Variable Inventory

See [data-model.md](data-model.md) for the full state inventory including
init/destroy owners and thread-access patterns.

**`src/dw_globals.c`** (new file) contents:

```c
/* src/dw_globals.c — single authoritative definition of all DW globals */
#include "dw.h"

pthread_attr_t   attr;
pthread_mutex_t  master_lp_ready_mutex;
pthread_cond_t   master_lp_ready_cv;
pthread_mutex_t  service_queue_mutex;
pthread_mutex_t  next_iteration_mutex;
pthread_cond_t   next_iteration_cv;
pthread_mutex_t  master_mutex;
pthread_mutex_t  reduced_cost_mutex;
pthread_mutex_t  glpk_mutex;
pthread_mutex_t  fputs_mutex;
pthread_mutex_t* sub_data_mutex;
#ifdef USE_NAMED_SEMAPHORES
sem_t*           customers;
#else
sem_t            customers;
#endif
glp_prob*        original_master_lp;
glp_prob*        master_lp;
glp_iocp*        parm;
glp_smcp*        simplex_control_params;
D_matrix*        D;
signal_data*     signals;
```

**`src/dw.h`** change — same lines, converted to `extern`:

```c
extern pthread_attr_t   attr;
extern pthread_mutex_t  master_lp_ready_mutex;
extern pthread_cond_t   master_lp_ready_cv;
/* ... all 18 symbols with extern prefix ... */
#ifdef USE_NAMED_SEMAPHORES
extern sem_t*  customers;
#else
extern sem_t   customers;
#endif
```

**`src/Makefile.am`** change:

```makefile
dwsolver_SOURCES = \
  ...
  dw_globals.c \
  dw_main.c
```

### 1.2 OOM Helper Design

```c
/* Appended to src/dw.h before #endif */
static inline void dw_oom_abort(void* ptr, const char* ctx) {
    if (!ptr) {
        fprintf(stderr, "dwsolver: out of memory in %s\n", ctx);
        exit(EXIT_FAILURE);
    }
}
```

Requires `<stdlib.h>` and `<stdio.h>` — both already included transitively via GLPK
headers, but should be explicitly guarded.

### 1.3 Quickstart

See [quickstart.md](quickstart.md) for full build and test instructions.

---

## Post-Design Constitution Re-Check

| Principle | Result | Notes |
|-----------|--------|-------|
| I. Correctness First | ✅ PASS | No solver logic changed; full 6-example test suite required before merge |
| II. Thread Safety | ✅ PASS | Moving definitions does not change init/destroy order or mutex semantics |
| III. Cross-Platform Portability | ✅ PASS | `extern` + single definition is the most portable C99 pattern |
| IV. Repair Before Extension | ✅ PASS | Only repairs; no new features |
| V. CLI-First, Library-Ready | ✅ PASS | Fewer globals in header improves eventual library extraction |

**All gates pass. Ready for task generation.**

---

## Implementation Order

1. **T001** — Verify baseline build (macOS)
2. **T002–T004** — KD-001 fix (create `dw_globals.c`, extern `dw.h`, update `Makefile.am`)
3. **T005–T006** — Rebuild and run 4 existing tests (US1 ✅ MVP)
4. **T007–T009** — Extend test suite with `four_sea` and `book_dantzig` (US2)
5. **T010–T015** — Replace 5 `sprintf` → `snprintf` (US3) — parallel with US2
6. **T016–T024** — Add `dw_oom_abort` helper + 7 malloc guards (US4)
7. **T025–T030** — Final audits, commit, push
