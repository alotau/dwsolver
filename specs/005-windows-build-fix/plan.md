# Implementation Plan: Code Quality Improvements — Bug Fixes, Duplication Reduction, and Clarity

**Branch**: `chore/code-quality` *(execute on this branch after `005-windows-build-fix` merges)* | **Date**: 2026-03-20 | **Spec**: [spec.md](./spec.md)
**Input**: Feature specification from `specs/005-windows-build-fix/spec.md`

## Summary

Address a prioritised set of bugs, code clarity issues, and performance
concerns identified by the 2026-03-20 codebase analysis. The most urgent items
are two correctness bugs in `dw_phases.c` and `dw_rounding.c` (undefined
behaviour on `val=NULL`, memory leak via mis-typed `const char*`). Medium-effort
work reduces the near-identical `phase_1_iteration`/`phase_2_iteration` pair to a
single unified function. Lower-priority items improve lock efficiency, remove
dead code, and (longest horizon) move core solver state from true globals into
the existing `faux_globals` / `master_data` structs.

**Non-negotiable constraint**: `tests/dw-tests.sh` must pass on all supported
platforms after every change. Numerical tolerances and convergence behaviour
must not change.

## Technical Context

**Language/Version**: C99 (targeting GCC ≥ 9 and Clang ≥ 12)
**Primary Dependencies**: GLPK 4.44 (thread-patched, embedded in `src/`), pthreads (POSIX or winpthreads on Windows)
**Storage**: File I/O only — CPLEX LP format input (`lpx_read_cpxlp`), text output files
**Testing**: Shell-based runner `tests/dw-tests.sh`; examples in `examples/` provide expected output for diffing
**Target Platform**: macOS, Linux, Windows (MSYS2/MINGW64)
**Project Type**: CLI solver (stand-alone binary `src/dwsolver`)
**Performance Goals**: Correctness is paramount; execution speed is second priority. Changes must not increase wall-clock time on the benchmark examples.
**Constraints**: Thread safety must be preserved or improved. No external dependencies may be added. C99 standard only (no C11 atomics or `_Thread_local` unless guarded).
**Scale/Scope**: ~4,000 lines across 7 source files + 6 headers; 13 test cases in `tests/dw-tests.sh`

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| **I. Correctness First** | ✅ PASS | Every PR in this plan requires `tests/dw-tests.sh` green before merge. P1 bug fixes actively improve correctness. |
| **II. Thread Safety is Mandatory** | ✅ PASS | P1 bug fixes do not touch threading. P2 (phase unification) and P3 (lock reduction) must be analysed with helgrind before merge. P4 (global state) is the highest-risk item and requires TSan sign-off. |
| **III. Cross-Platform Portability** | ✅ PASS | All changes are in portable C99. No platform-specific APIs introduced. CI validates macOS, Linux, Windows. |
| **IV. Repair Before Extension** | ✅ PASS | This spec is entirely repair — no new features are introduced. P1 fixes UB and a memory leak. P2 reduces duplication. |
| **V. CLI-First, Library-Ready** | ✅ PASS | Moving globals into structs (P4) directly advances library-readiness. CLI interface is unchanged. |

**Re-check after Phase 1 design**: Global state migration (P4 / User Story 6)
must be re-evaluated against Principle II after `data-model.md` is finalised.
If the struct layout introduces any shared mutable pointer accessible from
multiple threads without a lock, that is a gate failure.

## Project Structure

### Documentation (this feature)

```text
specs/005-windows-build-fix/
├── plan.md       ← this file
├── spec.md       ← feature specification
├── research.md   ← Phase 0 output (analysis findings)
├── data-model.md ← Phase 1 output (struct extension proposal)
├── quickstart.md ← Phase 1 output (build + verify steps)
└── tasks.md      ← Phase 2 output (/speckit.tasks — not yet created)
```

### Source Code (affected files)

```text
src/
├── dw.h              # P4: add master_lp, D, signals, sync_state to faux_globals
├── dw_globals.c      # P4: remove extern definitions that migrate to structs
├── dw_main.c         # P3: extract master LP build into helper; P4: pass structs
├── dw_phases.c       # P1: fix val=NULL UB; P2: unify phase_1/phase_2; P3: dead code
├── dw_rounding.c     # P1: fix const char* leak; P3: dead code removal
├── dw_subprob.c      # P3: guard glpk_mutex logging; P4: pass structs through
└── dw_support.c      # P3: replace direct printf with dw_printf; P4: struct updates

tests/
└── dw-tests.sh       # unchanged — regression gate throughout
```

**Structure Decision**: Single project, existing `src/` layout. No new files
required for P1–P3. P4 (global state migration) may introduce a
`dw_sync.h`/`dw_sync.c` pair for the synchronisation state struct if the
resulting `faux_globals` becomes unwieldy.

## Complexity Tracking

> No constitution violations that require justification. All items are repairs
> (Principle IV). The architectural change in P4 is a structural improvement,
> not an added abstraction layer.
