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
