---
description: Audit dwsolver source against SEI CERT C (spec-008), ISO/IEC TS 17961 (spec-014), and Strict C11/POSIX.1-2008 (spec-013). Run after merging new features to catch regressions and flag stale acceptance reports.
---

## User Input

```text
$ARGUMENTS
```

Consider any user input above (e.g., a specific file or rule to focus on) before proceeding.

## Goal

Perform a non-destructive compliance health check across all three active standards for
the `src/dw_*.c` and `src/dw_*.h` authored source files. Produce a structured audit
report. **Do NOT modify any source or spec files.** All findings are advisory; any
remediation is a separate step requiring explicit user approval.

---

## Phase 1 — Static Analysis

Run all commands from the repository root. If a tool is absent, note it in the report
and skip that check rather than aborting.

### 1a. Clean build — zero-warning gate (SEI CERT C FR-007 / spec-013 FR-001)

```sh
make clean && make 2>&1 | grep -E 'warning:|error:' | grep -v 'glp'
```

**Pass**: zero output lines. Any line is a finding.

### 1b. clang-tidy CERT checks (spec-014 CI gate)

```sh
GLPK_CFLAGS=$(pkg-config --cflags glpk 2>/dev/null || echo "-I/usr/local/include")
clang-tidy -checks='cert-*' \
  src/dw_blas.c src/dw_globals.c src/dw_main.c src/dw_phases.c \
  src/dw_rounding.c src/dw_solver.c src/dw_subprob.c src/dw_support.c \
  -- -std=c11 -D_POSIX_C_SOURCE=200809L -I src/ ${GLPK_CFLAGS} 2>&1
```

Diff output against `specs/008-sei-cert-c-compliance/baseline-warnings.txt`.
**Pass**: zero net-new `cert-*` findings vs baseline.

### 1c. cppcheck (spec-008 FR-007 / spec-013 FR-007)

```sh
cppcheck --enable=all --suppress=missingIncludeSystem \
  -I src/ src/dw_*.c 2>&1
```

Diff against `specs/008-sei-cert-c-compliance/baseline-cppcheck.txt` if present.
**Pass**: zero net-new warnings vs baseline.

### 1d. Strict C11 / POSIX pedantic compile (spec-013 FR-001, FR-002)

```sh
GLPK_CFLAGS=$(pkg-config --cflags glpk 2>/dev/null || echo "-I/usr/local/include")
gcc -std=c11 -pedantic-errors -Wall -Wextra -Wsign-compare \
  -D_POSIX_C_SOURCE=200809L \
  -c src/dw_*.c -I src/ ${GLPK_CFLAGS} 2>&1 | grep -v 'glp'
```

**Pass**: zero errors and zero warnings.

---

## Phase 2 — Source Staleness Detection

Determine when each acceptance report was last committed, then find any authored source
files modified after that date.

### 2a. Get acceptance report commit dates

```sh
git log --format="%ad" --date=short -1 -- \
  specs/008-sei-cert-c-compliance/acceptance-report.md

git log --format="%ad" --date=short -1 -- \
  specs/014-iso-ts-17961-compliance/acceptance-report.md

git log --format="%ad" --date=short -1 -- \
  specs/013-strict-c-posix-compliance/acceptance-report.md
```

### 2b. Find source changes after each date

For each date `<DATE>` returned by 2a, run:

```sh
git log --oneline --after="<DATE>" -- src/dw_*.c src/dw_*.h
```

Record: source file, commit hash, date, and one-line commit message for each hit.

### 2c. Check compliance contract docs for staleness

Run the same `git log --oneline --after="<DATE>"` check against the source files
explicitly cited as evidence in each contract document, using the date the contract
doc was last committed:

```sh
git log --format="%ad" --date=short -1 -- \
  specs/008-sei-cert-c-compliance/contracts/pthread-call-sites.md
git log --format="%ad" --date=short -1 -- \
  specs/008-sei-cert-c-compliance/contracts/lock-order.md
git log --format="%ad" --date=short -1 -- \
  specs/008-sei-cert-c-compliance/contracts/glpk-lock-state.md
git log --format="%ad" --date=short -1 -- \
  specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md
git log --format="%ad" --date=short -1 -- \
  specs/013-strict-c-posix-compliance/contracts/posix-header-audit.md
```

A contract doc that is older than the source files it documents is itself stale —
its line-number citations and evidence may no longer match the current code.

---

## Phase 3 — Rule-to-File Mapping

Read the following documents to determine which compliance rules cite each source file
as evidence:

- `specs/008-sei-cert-c-compliance/contracts/pthread-call-sites.md` — maps CERT
  POS54-C to call sites by file and line number
- `specs/008-sei-cert-c-compliance/contracts/lock-order.md` — maps CON31-C / CON43-C
  to mutex acquisition sites
- `specs/008-sei-cert-c-compliance/contracts/glpk-lock-state.md` — maps GLPK call
  sites to lock-hold state assertions
- `specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md` — maps all
  22 TS 17961 rules to file/function evidence
- `specs/013-strict-c-posix-compliance/contracts/posix-header-audit.md` — maps POSIX
  API usage to source files

For each source file that changed since an acceptance report date (from Phase 2b),
identify every compliance rule that cites that file as evidence. These rules are
**potentially stale** and require human re-audit before the compliance claim can be
considered current.

---

## Phase 4 — Produce Audit Report

Output the following report. Do not truncate any section.

---

### DWSOLVER Compliance Audit Report

**Date**: (today's date)
**Scope**: `src/dw_*.c`, `src/dw_*.h`
**Invocation arguments**: (repeat $ARGUMENTS verbatim, or "none")

#### Overall Status

| Standard | Static Analysis | Source Staleness | Verdict |
|----------|----------------|-----------------|---------|
| SEI CERT C (spec-008) | ✅ PASS / ❌ FAIL / ⚠️ TOOL MISSING | ✅ Current / ⚠️ N file(s) changed | ✅ CLEAN / ⚠️ REVIEW NEEDED |
| ISO/IEC TS 17961 (spec-014) | … | … | … |
| Strict C11/POSIX.1-2008 (spec-013) | … | … | … |

#### Static Analysis Findings

For each tool that ran, list only net-new findings (not already present in the
baseline). Group by tool. If all tools are clean, state:
"No new findings from any tool."

If a tool was absent, state: "SKIPPED — `<toolname>` not found in PATH."

#### Source Files Changed Since Last Audit

| File | Commit | Date | Message | Standards Affected |
|------|--------|------|---------|-------------------|
| … | … | … | … | SEI CERT C / TS 17961 / C11+POSIX |

If no files changed since any acceptance report, state:
"All source files are current with respect to all three acceptance reports."

#### Potentially Stale Compliance Rules

For each rule whose evidence file changed after the corresponding acceptance report:

| Rule ID | Standard | Evidence File(s) Changed | Recommended Action |
|---------|----------|-------------------------|--------------------|
| POS54-C | SEI CERT C | `src/dw_subprob.c` (commit abc1234, 2026-03-25) | Re-audit pthread call sites; update `contracts/pthread-call-sites.md` line citations |
| … | … | … | … |

If none: "No compliance rules are stale."

#### Stale Contract Documents

List any contract or matrix file that is older than the source files it documents
(detected in Phase 2c). These need line-number and evidence updates but do not
necessarily indicate a compliance violation.

#### Recommended Actions

Ordered list, most urgent first:

1. **(If static analysis failures)** Fix net-new `cert-*` or compiler warnings before
   merging further changes to `main`. These are hard regressions.
2. **(If stale rules)** For each entry in the Stale Compliance Rules table: re-audit
   the cited rule against the changed file, update the acceptance report verdict and
   evidence, and refresh the contract doc line-number citations.
3. **(If stale contract docs only)** Update line-number references in the cited
   contract docs to match current source. No acceptance-report change required unless
   a rule's correctness is also affected.
4. **(If all clean)** No action required. Compliance posture is current with respect
   to all three standards.

---

## Operating Constraints

- **Strictly read-only**: Do not edit, create, or delete any file.
- **Graceful tool absence**: If `clang-tidy` or `cppcheck` is not installed, note it
  in the report and continue — do not abort.
- **pkg-config fallback**: If `pkg-config --cflags glpk` fails, substitute
  `-I/usr/local/include` as the GLPK include path.
- **Scope boundary**: Only `src/dw_*.c` and `src/dw_*.h` are in scope. GLPK source
  files (`src/glp*.c`, `src/glp*.h`) are explicitly excluded from all checks.
- **No remediation**: This agent reports findings only. Fixes require explicit user
  approval and must be implemented in a separate session.
- **Focused invocation**: If $ARGUMENTS names a specific file (e.g., `dw_phases.c`)
  or rule (e.g., `POS54-C`), restrict Phase 1 compilation and Phase 3 mapping to that
  file or rule; still run all four phases and produce the full report structure.
