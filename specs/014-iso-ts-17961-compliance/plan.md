# Implementation Plan: ISO/IEC TS 17961:2013 Compliance

**Branch**: `014-iso-ts-17961-compliance` | **Date**: 2026-03-23 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/014-iso-ts-17961-compliance/spec.md`

## Summary

Audit all 22 ISO/IEC TS 17961:2013 analyzable rules against dwsolver's
authored source (`src/dw_*.c`, `src/dw_*.h`). Initial audit finds 20 rules
already PASS or N-A via spec 008/013 prior work. Two rules FAIL and require
targeted remediation in the rounding subsystem: **ga-buffer** (3 `strcpy` calls
from external LP column names into bounded buffers) and **nonnullptr** (5 `strtok`
return values passed to `strcpy` without NULL check). Deliver a formal
compliance matrix, reproducible tooling guide, CI workflow, and acceptance
report suitable for procurement.

## Technical Context

**Language/Version**: C11 (`-std=c11 -pedantic-errors`)
**Primary Dependencies**: GLPK ≥ 4.65, pthreads; cppcheck ≥ 2.12 (CI tool), clang-tidy ≥ 17 (CI tool)
**Storage**: Markdown documents in `specs/014-iso-ts-17961-compliance/audit/`
**Testing**: `tests/dw-tests.sh` (regression), `tests/test_ts17961_enforcement.sh` (new — SC-003 gate)
**Target Platform**: Linux CI runner (`ubuntu-latest`); macOS for local dev
**Project Type**: Compliance audit + targeted source remediation + CI workflow
**Performance Goals**: N/A — correctness audit; zero new build diagnostics
**Constraints**: Changes confined to `src/dw_rounding.c`, `src/dw_support.c`; no functional changes outside the rounding subsystem; `third-party/glpk/` must not be modified
**Scale/Scope**: 22 rules × 1 matrix row each; 8 call sites remediated across 2 files

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-checked after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | ✅ PASS | FR-009: full test suite must pass; no solver logic changes |
| II. Thread Safety | ✅ PASS | No threading changes; datarace audited as PASS (spec 008) |
| III. Cross-Platform Portability | ✅ PASS | CI workflow Linux-only per spec; source changes use standard C11 |
| IV. Repair Before Extension | ✅ PASS | This is a compliance audit of existing repaired code; no new features |
| V. CLI-First, Library-Ready | ✅ PASS | Source changes in rounding subsystem only; no I/O path changes |

No violations. No complexity tracking needed.

## Project Structure

### Documentation (this feature)

```text
specs/014-iso-ts-17961-compliance/
├── plan.md                           # This file
├── research.md                       # Phase 0: per-rule audit findings
├── data-model.md                     # Phase 1: entities and modified files
├── quickstart.md                     # Phase 1: developer quick reference
├── audit/
│   ├── ts17961-compliance-matrix.md  # Primary compliance artifact (22 rules)
│   └── tooling-guide.md              # Reproducible tool invocations
└── tasks.md                          # Phase 2 output (/speckit.tasks — NOT by this command)
```

### Source Code Changes

```text
src/
├── dw_rounding.c    # ga-buffer fix (lines 425, 445) + nonnullptr fix (lines 427-428, 447-448)
└── dw_support.c     # ga-buffer fix (line 668) + nonnullptr fix (line 673)

tests/
└── test_ts17961_enforcement.sh       # SC-003 synthetic violation gate (new)

.github/workflows/
└── ci-compliance.yml                 # Dedicated TS 17961 CI workflow (new)
```

**Structure Decision**: Single-project C library. No contracts directory needed
(this feature produces internal compliance documents, not external API contracts).
