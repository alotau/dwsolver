# Data Model: ISO/IEC TS 17961:2013 Compliance

**Phase 1 — Design**
**Branch**: `014-iso-ts-17961-compliance`
**Date**: 2026-03-23

This spec is a compliance audit and remediation feature, not a data-modelling
feature. The "entities" are the artifacts defined by the spec and the code
constructs being analysed.

---

## 1. Primary Artifact Entities

### ComplianceMatrix (document)

`specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md`

Each row = one TS 17961 rule.

| Field           | Type   | Constraint                              |
|-----------------|--------|-----------------------------------------|
| rule_id         | string | one of the 22 TS 17961 rule mnemonics   |
| rule_name       | string | descriptive name from the TS            |
| cert_mapping    | string | nearest SEI CERT C rule or "none"       |
| prior_work      | string | "spec 008", "spec 013", "none", or "—"  |
| verdict         | enum   | PASS \| FAIL \| N-A                     |
| evidence        | string | citation or grep command + output       |
| remediation_ref | string | task/FR reference if verdict == FAIL    |

**Invariants**:
- Final matrix MUST have zero FAIL verdicts (FR-004).
- Every N-A verdict MUST include a grep-based absence proof (FR-001).
- Every PASS via prior work MUST cite the specific FR from spec 008/013.

### ConformanceRule (enumeration)

The 22 rule mnemonics are the fixed domain:

```
intptrconv, intobjptr, xplicitcomp, boolasgn, trstcmp,
wraparound, argcomp, intnegator, ioilecc, ga-buffer,
strmod, memiograph, accfree, dblfree, nonnullptr,
nullref, diverr, exceptbits, liberr, datarace,
asyncsig, toinit
```

These are immutable — the 2013 edition is the authoritative scope (FR scope
assumption).

### ToolingGuide (document)

`specs/014-iso-ts-17961-compliance/audit/tooling-guide.md`

| Field                | Description                                     |
|----------------------|-------------------------------------------------|
| tool_name            | cppcheck, clang-tidy, Coverity, Klocwork, etc.  |
| version_requirement  | e.g., "≥ 2.12"                                  |
| invocation           | exact command-line string                       |
| exclusion_pattern    | how to exclude third-party/glpk                 |
| checker_names        | tool-specific names for TS 17961 rule set       |
| output_interpretation| how to read findings for this project           |

### AcceptanceReport (document)

`specs/014-iso-ts-17961-compliance/acceptance-report.md`

Format follows `specs/008-sei-cert-c-compliance/acceptance-report.md`:
one section per FR (FR-001 through FR-010), each ending with a ✅ PASS verdict.

---

## 2. Code Entities Modified

### StrncpyCallSite (source change)

Three `strcpy` calls in the rounding/support subsystem replaced with
`strncpy` + explicit NUL-termination (ga-buffer remediation):

| File                | Line  | Description                                         |
|---------------------|-------|-----------------------------------------------------|
| src/dw_rounding.c   | 425   | `strcpy(local_col_name, var_name)` → strncpy        |
| src/dw_rounding.c   | 445   | `strcpy(local_col_name, var_name)` → strncpy        |
| src/dw_support.c    | 668   | `strcpy(local_col_name, var_name)` → strncpy        |

### NonnullptrGuard (source change)

Five `strcpy(dst, strtok(...))` call sites receive a NULL guard and strncpy
replacement (nonnullptr remediation):

| File                | Lines   | Description                                         |
|---------------------|---------|-----------------------------------------------------|
| src/dw_rounding.c   | 427–428 | `curr_flight`, `curr_sector` strtok tokens          |
| src/dw_rounding.c   | 447–448 | `temp_flight`, `prev_sector` strtok tokens          |
| src/dw_support.c    | 673     | `sector_name` strtok token                          |

**State transitions for NonnullptrGuard**:
```
strtok() called
  → token != NULL  → strncpy + NUL-terminate → continue processing
  → token == NULL  → log diagnostic → continue (skip column)
```

Skipping on NULL is safe: the outer loops iterate over LP column ranges;
a malformed column name is non-fatal and already printed as a diagnostic
in the original code.

---

## 3. New Files Created

| Path                                                              | Description                            |
|-------------------------------------------------------------------|----------------------------------------|
| `specs/014-iso-ts-17961-compliance/research.md`                   | Phase 0 audit findings (this spec)     |
| `specs/014-iso-ts-17961-compliance/data-model.md`                 | This file                              |
| `specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md` | Primary compliance artifact        |
| `specs/014-iso-ts-17961-compliance/audit/tooling-guide.md`        | Reproducible tooling instructions      |
| `specs/014-iso-ts-17961-compliance/acceptance-report.md`          | FR sign-off document                   |
| `tests/test_ts17961_enforcement.sh`                               | SC-003 synthetic violation gate        |
| `.github/workflows/ci-compliance.yml`                             | Dedicated compliance CI workflow       |

---

## 4. Validation Rules

- **VR-001**: After remediation, `cppcheck --addon=cert src/dw_*.c` MUST produce
  zero diagnostics mapping to `ga-buffer` or `nonnullptr` (cert-STR31-C,
  cert-EXP34-C).
- **VR-002**: After remediation, `-std=c11 -pedantic-errors -Wall -Wextra` build
  MUST produce zero new diagnostics vs. spec 013 baseline.
- **VR-003**: All 22 rules MUST appear in the final compliance matrix with a
  PASS or N-A verdict.
- **VR-004**: For every N-A verdict, the compliance matrix MUST contain a grep
  command whose output is reproduced inline and shows zero results.
- **VR-005**: `tests/dw-tests.sh` MUST pass after all source changes.
