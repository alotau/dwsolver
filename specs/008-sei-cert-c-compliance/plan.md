# Implementation Plan: SEI CERT C Standard Compliance

**Branch**: `008-sei-cert-c-compliance` | **Date**: 2026-03-21 | **Spec**: [spec.md](spec.md)  
**Input**: Feature specification from `specs/008-sei-cert-c-compliance/spec.md`

## Summary

Apply SEI CERT C Coding Standard rules to `src/dw_*.c` and `src/dw_*.h`, targeting the five remaining non-conformance categories discovered after the feature-007 architecture review: unchecked `pthread_*`/`sem_*` return values (POS54-C), blocking GLPK solves under mutex (POS52-C), undocumented lock acquisition order (POS51-C), remaining `malloc`/`fopen` null-check gaps (MEM32-C/FIO), and signed/unsigned integer type hazards (INT30/31-C). Implementation proceeds in user-story priority order (P1 → P3) so that the highest-risk threading hazards are addressed before lower-risk improvements.

## Technical Context

**Language/Version**: C99 (GCC/Clang on macOS/Linux; MinGW-w64 on Windows)  
**Primary Dependencies**: POSIX Threads (pthreads); embedded GLPK 4.44 (treated as a library boundary — `src/glp*.c` not modified)  
**Storage**: N/A (solver writes solution files; no database)  
**Testing**: Shell-based example runner (`tests/dw-tests.sh`); ThreadSanitizer build (`-fsanitize=thread`) for threading changes  
**Target Platform**: macOS, Linux, Windows (cross-platform portability required by constitution)  
**Project Type**: CLI tool / solver library  
**Performance Goals**: No solver performance regression; `pthread_*` error checks add negligible overhead (single integer comparison per call)  
**Constraints**: All changes confined to `src/dw_*.c` / `src/dw_*.h`; GLPK files must not be modified; all existing tests must keep passing  
**Scale/Scope**: 7 authored source files; ~50 pthread/sem call sites; ~100 malloc/calloc sites; 12 glp_simplex/glp_intopt calls; ~5 fopen calls

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Verdict | Notes |
|-----------|---------|-------|
| I. Correctness First | **PASS** — Changes do not touch solver math. POS52-C fix (releasing mutex before glp_simplex) requires careful analysis to ensure data availability on re-acquire; covered by existing test suite. |
| II. Thread Safety is Mandatory | **PASS** — This feature *improves* thread safety. POS52-C fix and POS51-C audit directly address threading requirements. TSan verification required for all threading changes per constitution. |
| III. Cross-Platform Portability | **PASS** — `pthread_*` error checks are POSIX-standard; no platform-specific APIs introduced. `DW_PTHREAD_CHECK` macro and `dw_pthread_abort` helper will use only standard `strerror`/`fprintf`. Windows pthreads-win32 returns the same POSIX error codes. |
| IV. Repair Before Extension | **PASS** — This is a pure repair/hardening effort with no new features. |
| V. CLI-First, Library-Ready | **N/A** — No interface changes. |

**Post-design re-check**: No constitution violations identified after Phase 1 design. POS52-C fix (release `sub_data_mutex[id]` before `glp_simplex`) preserves correctness because objective coefficients are fully set into the GLPK LP object before release, and results are read back under re-acquired lock.

## Project Structure

### Documentation (this feature)

```text
specs/008-sei-cert-c-compliance/
├── plan.md              # This file
├── research.md          # Phase 0 output — call-site inventory and findings
├── data-model.md        # Phase 1 output — sync primitive model and error-handling pattern
├── contracts/
│   ├── pthread-call-sites.md   # POS54-C audit table (all ~50 call sites)
│   ├── glpk-lock-state.md      # POS52-C audit table (all 12 glp_simplex/glp_intopt sites)
│   └── lock-order.md           # POS51-C total order table with rationale
├── quickstart.md        # Phase 1 output — developer guide for error-check pattern
└── tasks.md             # Phase 2 output (/speckit.tasks command)
```

### Source Code

```text
src/
├── dw_main.c        # Master thread: pthread_create/join, all mutex lock/unlock/broadcast
├── dw_phases.c      # Phase iteration functions: lock/unlock across 4 primitives
├── dw_subprob.c     # Subproblem thread: POS52-C violation (glp_simplex under lock),
│                    #   additional spurious-wakeup bug (if→while at line 134)
├── dw_support.c     # init_pthread_data(): all 8 pthread_mutex_init + 2 cond_init unchecked
├── dw_rounding.c    # Rounding threads: pthread_create/join (already checked), lock/unlock
├── dw_globals.c     # Global declarations only — no changes expected
└── dw_blas.c        # Math utilities — no threading; minor INT30/31 fixes possible

tests/
└── dw-tests.sh      # Full example suite — must pass before and after
```

**Structure Decision**: Single-project; all changes in `src/`; no new directories in source tree. Audit documents added to `specs/008-sei-cert-c-compliance/contracts/`.

## Complexity Tracking

No constitution violations requiring justification. This is a hardening feature with net-negative complexity (removes latent bugs, adds no new abstractions beyond a helper macro).
