# Tasks: ISO/IEC TS 17961:2013 Compliance

**Input**: Design documents from `/specs/014-iso-ts-17961-compliance/`
**Branch**: `014-iso-ts-17961-compliance`
**Generated**: 2026-03-23
**Prerequisites**: plan.md ✅, spec.md ✅, research.md ✅, data-model.md ✅, quickstart.md ✅

**Scope**: 22 ISO/IEC TS 17961:2013 rules audited; 2 FAIL rules (ga-buffer, nonnullptr)
remediated in `src/dw_rounding.c` and `src/dw_support.c`; compliance matrix, tooling guide,
CI workflow, enforcement script, and acceptance report produced.

**Design note**: Phase 2 (Foundational) is omitted — all research and design work is
complete (research.md and data-model.md authored in Phase 0/1). Baseline capture in
Phase 1 is the sole blocking prerequisite for US2 and US3.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel with other [P]-marked tasks (different files, no interdependency)
- **[Story]**: User story this task belongs to (US1 = matrix, US2 = remediation, US3 = tooling/CI, US4 = procurement)
- Exact file paths are included in every task description

---

## Phase 1: Setup — Pre-Remediation Baseline Capture

**Purpose**: Capture the cppcheck and clang-tidy outputs *before* any source changes.
These outputs are required by US2 (before/after diff for SC-002) and US3 (tooling guide).

**⚠️ CRITICAL**: T001 and T002 MUST complete before any source changes in US2.

- [X] T001 Run `cppcheck --addon=cert --std=c11 --suppress=constParameterPointer -I src src/dw_blas.c src/dw_globals.c src/dw_main.c src/dw_phases.c src/dw_rounding.c src/dw_solver.c src/dw_subprob.c src/dw_support.c` and save full output to `specs/014-iso-ts-17961-compliance/audit/cppcheck-pre-remediation.txt`
- [X] T002 [P] Run `clang-tidy -checks='cert-*' --extra-arg=-std=c11 --extra-arg=-I./src src/dw_blas.c src/dw_globals.c src/dw_main.c src/dw_phases.c src/dw_rounding.c src/dw_solver.c src/dw_subprob.c src/dw_support.c -- -I/usr/include -DHAVE_CONFIG_H` and save full output to `specs/014-iso-ts-17961-compliance/audit/clang-tidy-pre-remediation.txt` (**I-002**: scope widened to all `src/dw_*.c` to match the CI workflow in T014; any cert-* findings in non-rounding files are recorded in the pre-remediation baseline)

**Checkpoint**: Pre-remediation baseline captured — US1, US2, and US3 can now proceed.

---

## Phase 3: User Story 1 — Full Compliance Matrix (Priority: P1) 🎯 MVP

**Goal**: A procurement engineer can read a single document covering all 22 TS 17961 rules
with PASS/FAIL/N-A verdicts, evidence citations, SEI CERT C mappings, and prior-work
references — with no open TBD entries.

**Independent Test**: Inspect `specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md`.
It must: (1) list all 22 rule mnemonics, (2) have a PASS or N-A verdict for every non-gap row,
(3) cite evidence (grep command + zero-result output for N-A; FR citation for prior-work PASS),
(4) mark ga-buffer and nonnullptr as FAIL with a remediation cross-reference, and (5) contain
no TBD or empty evidence cells.

### Implementation for User Story 1

- [X] T003 [US1] Author document header, table schema, and all N-A verdict rows with inline grep absence proofs (intptrconv, xplicitcomp, boolasgn, intnegator, exceptbits, asyncsig — 6 rules) in `specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md`; schema columns: Rule ID | Rule Name | SEI CERT C | Prior Work | Verdict | Evidence
- [X] T004 [P] [US1] Add PASS rows crediting spec 013 prior work for intobjptr (FR-013-pedantic), trstcmp (structural invariant — dwsolver-controlled column names), wraparound (FR-013-wall), argcomp (FR-013-pedantic), ioilecc (FR-013-wformat), and toinit (FR-013-uninitvar) in `specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md`
- [X] T005 [P] [US1] Add PASS rows crediting spec 008 prior work for memiograph (FR-008-mem31), accfree (FR-008-mem31), dblfree (FR-008-mem31), nullref (dw_oom_abort structural invariant), liberr (FR-008-err33), datarace (FR-008-pos51/52/54); and standalone PASS rows for strmod (grep evidence — all strcpy destinations are malloc'd buffers) and diverr (structural invariant — num_clients validated > 0; / 1.0e9 compile-time constant) in `specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md`
- [X] T006 [US1] Add initial FAIL rows for ga-buffer (3 unbounded strcpy call sites: dw_rounding.c:425,445; dw_support.c:668) and nonnullptr (5 strtok-to-strcpy call sites: dw_rounding.c:427-428,447-448; dw_support.c:673) with remediation cross-references to T007 and T008 in `specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md`

**Checkpoint**: Compliance matrix has all 22 rules populated (2 FAIL with cross-references, 20 PASS/N-A). Story 1 independently testable.

---

## Phase 4: User Story 2 — Gap Rule Remediation (Priority: P2)

**Goal**: All 22 TS 17961 rules achieve a clean PASS or N-A verdict in the codebase.
After remediation, cppcheck --addon=cert reports zero cert-STR31-C and cert-EXP34-C
findings in src/dw_*.c, and the full test suite still passes.

**Independent Test**: (a) Run `grep -n "strcpy" src/dw_rounding.c src/dw_support.c` — no
unbounded strcpy copying from glp_get_col_name() return values remains. (b) Run
`cppcheck --addon=cert src/dw_*.c` — zero findings map to ga-buffer or nonnullptr (VR-001).
(c) Run `tests/dw-tests.sh` — all tests pass (FR-009). (d) Build with
`-std=c11 -pedantic-errors -Wall -Wextra` — zero new diagnostics vs. spec 013 baseline (VR-002).

### Implementation for User Story 2

- [X] T007 [US2] Fix ga-buffer rule: **Before implementing**, verify (a) the exact constant symbol used for the destination buffer size in `src/dw_rounding.c` (search for `#define` near the local buffer declaration — confirm name is `BUFF_SIZE` and value is `200`), and (b) confirm GLPK column name limit from the system-installed glpk.h used by the build (e.g., `/usr/include/glpk.h` on Linux or the path reported by `pkg-config --cflags glpk`; search `GLP_MAX_NAME_LEN` or as a comment near the name-length limit — expected 255; note: `third-party/glpk/` contains only docs/patches, not headers). Update the annotation below if either value differs. Then replace `strcpy(local_col_name, var_name)` with `strncpy(local_col_name, var_name, BUFF_SIZE - 1); local_col_name[BUFF_SIZE - 1] = '\0';` at `src/dw_rounding.c` lines 425 and 445, and `src/dw_support.c` line 668; add `/* TS17961-ga-buffer: strncpy caps at BUFF_SIZE-1 (verified=<actual_size>); GLPK name limit is <actual_glpk_limit> > <actual_size> */` annotation at each site (**U-002/I-004**: symbol name and size values must be verified from source before the annotation is written into production code)
- [X] T008 [US2] Fix nonnullptr rule: at each of the 5 call sites extract strtok return to a named `const char *` variable, add NULL guard (branch on NULL: log diagnostic with `dw_printf(IMPORTANCE_DIAG, ...)` then `continue`), replace `strcpy` with `strncpy` + NUL-termination; sites: `src/dw_rounding.c` lines 427–428 (`curr_flight`, `curr_sector`), lines 447–448 (`temp_flight`, `prev_sector`), and `src/dw_support.c` line 673 (`sector_name`)
- [X] T009 [US1|US2] Update ga-buffer and nonnullptr matrix rows from FAIL to PASS in `specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md`; evidence column must name the post-remediation cppcheck clean run and cite this task; verify the matrix now has zero FAIL verdicts (FR-004)
- [X] T010 [US2] Build `src/dw_*.c` with `make` and confirm `gcc -std=c11 -pedantic-errors -Wall -Wextra` produces zero new diagnostics compared to the spec 013 post-implementation baseline captured in `specs/013-strict-c-posix-compliance/warnings-post.txt`; if that file is absent, run a clean build on the unmodified branch tip before applying T007/T008 changes and record the diagnostic count as the local baseline; run `tests/dw-tests.sh` and confirm all tests pass (FR-009, **U-004**: baseline reference is now explicit)

**Checkpoint**: Source is clean for ga-buffer and nonnullptr. Compliance matrix has zero FAIL rows. Story 2 independently testable.

---

## Phase 5: User Story 3 — Tooling and CI Enforcement (Priority: P2)

**Goal**: Any developer with cppcheck ≥ 2.12 and clang-tidy ≥ 17 can reproduce the
compliance verification by following the tooling guide verbatim. The CI pipeline
prevents new TS 17961 violations from being merged. The synthetic violation gate (SC-003)
proves the CI enforcement is functional.

**Independent Test**: (a) Follow `specs/014-iso-ts-17961-compliance/audit/tooling-guide.md`
on a clean checkout and confirm outputs match documented findings. (b) Run
`bash tests/test_ts17961_enforcement.sh` — must exit 0 (gate works: synthetic violation
detected by cppcheck, then reverted). (c) CI workflow file exists at
`.github/workflows/ci-compliance.yml` and lint-validates as valid GitHub Actions YAML.

### Implementation for User Story 3

- [X] T011 [US3] Author mandatory tool sections in `specs/014-iso-ts-17961-compliance/audit/tooling-guide.md`: tool version requirements (cppcheck ≥ 2.12, clang-tidy ≥ 17), exact cppcheck `--addon=cert --std=c11 --suppress=constParameterPointer -I src` invocation against `src/dw_*.c`, GLPK exclusion pattern (`--suppress-dir=third-party`), exact clang-tidy `-checks='cert-*'` invocation against all `src/dw_*.c`, output interpretation guide, and the before/after diff section must include the **exact** `diff` command: `diff specs/014-iso-ts-17961-compliance/audit/cppcheck-pre-remediation.txt specs/014-iso-ts-17961-compliance/audit/cppcheck-post-remediation.txt` with the expected output showing only the `cert-STR31-C` and `cert-EXP34-C` lines removed — demonstrating SC-002 mechanically (**A-001**: diff format is a literal file diff, not a narrative comparison)
- [X] T012 [P] [US3] Add optional commercial tool section to `specs/014-iso-ts-17961-compliance/audit/tooling-guide.md` covering Coverity (checker categories: `BUFFER_SIZE`, `NULL_RETURNS`, `STRING_OVERFLOW`), Klocwork (`SV.STRBO.*`, `NPE.*`), and CodeSonar (`bufferoverrun`, `nullptr`); include GLPK exclusion pattern for each, and note that these tools are not part of the mandatory CI but are sufficient for procurement-grade certification evidence (FR-007)
- [X] T013 [P] [US3] Create `tests/test_ts17961_enforcement.sh` implementing SC-003: (1) shebang `#!/usr/bin/env bash` with `set -e` for general error handling, (2) copy `src/dw_rounding.c` to a temp file `$(mktemp)`, (3) inject a `strcpy(buf, some_ptr);` line without bounds check into the temp file using `sed` or `echo >>`, (4) **disable `set -e` with `set +e` before the cppcheck call**, run `cppcheck --addon=cert` on the temp file and capture exit code with `result=$?`, then **re-enable `set -e`**, (5) assert `[ "$result" -ne 0 ] || { echo "FAIL: cppcheck did not detect synthetic ga-buffer violation" >&2; rm -f "$tmpfile"; exit 1; }`, (6) remove temp file and exit 0; the script must exit 0 when the gate is functional and non-zero when cppcheck fails to catch the violation (**I-001**: bare `set -e` around an intentional non-zero cppcheck invocation causes spurious CI failure; the `set +e` / capture / `set -e` bracket is the correct pattern)
- [X] T019 [US3] Create `specs/014-iso-ts-17961-compliance/audit/cppcheck-na-suppressions.txt` listing the cppcheck checker IDs for the 6 N-A rules: `cert-INT36-C` (intptrconv), `cert-MSC21-C` (xplicitcomp), `cert-INT11-C` (boolasgn), `cert-INT32-C` (intnegator), `cert-FLP30-C` (exceptbits), `cert-SIG30-C` (asyncsig) — one ID per line; this file is the tracked exclusion mechanism used by T014's CI step and must be committed to the repository so auditors can see exactly which rules are excluded and why (**U-001**: concrete exclusion mechanism for N-A rules in CI; file is an auditable artifact)
- [X] T014 [US3] Create `.github/workflows/ci-compliance.yml` with: `runs-on: ubuntu-latest`, tool installation step (`apt-get install -y cppcheck clang-tidy`), tool version verification step (fail with clear message if cppcheck < 2.12 or clang-tidy < 17), cppcheck `--addon=cert --suppressions-list=specs/014-iso-ts-17961-compliance/audit/cppcheck-na-suppressions.txt` step against `src/dw_*.c` (fail on any non-suppressed finding — suppressions file from T019 defines exactly which N-A rule IDs are excluded), clang-tidy `-checks='cert-*'` step against all `src/dw_*.c`, and `bash tests/test_ts17961_enforcement.sh` step; T019 must complete before T014 (FR-006, SC-003, **U-001**: suppressions mechanism now specified)

**Checkpoint**: Tooling guide is reproducible. N-A suppressions file (T019) committed. CI workflow is authored. Enforcement script is executable. Story 3 independently testable.

---

## Phase 6: User Story 4 — Procurement Package (Priority: P3)

**Goal**: A procurement officer can navigate from the compliance matrix to all supporting
evidence without asking clarifying questions. The acceptance report consolidates all 10
FR verdicts in a self-contained sign-off document.

**Independent Test**: Provide `specs/014-iso-ts-17961-compliance/` to a reviewer familiar
with TS 17961 but unfamiliar with dwsolver — they can reach any evidence from the compliance
matrix in two clicks/links without clarifying questions (SC-005).

### Implementation for User Story 4

- [X] T015 [US4] Author `specs/014-iso-ts-17961-compliance/acceptance-report.md` with one section per FR-001 through FR-010 following the format of `specs/008-sei-cert-c-compliance/acceptance-report.md`; each section states the requirement, cites evidence (task IDs, document paths, or grep commands), and ends with `✅ PASS`; summary header gives the consolidated 22/22 rule verdict; zero open FRs or FAIL verdicts (FR-008, SC-001)

**Checkpoint**: Compliance package complete: matrix + tooling guide + acceptance report all authored and cross-referenced. SC-005 reviewer-experience check (T020) recorded in acceptance report.

---

## Polish & Cross-Cutting Verification

**Purpose**: Final validation sweep confirming all success criteria are met before closing the spec.

- [X] T016 [P] Run `cppcheck --addon=cert --std=c11 --suppress=constParameterPointer -I src src/dw_*.c` and confirm zero findings map to ga-buffer (cert-STR31-C) or nonnullptr (cert-EXP34-C); document the clean run output in `specs/014-iso-ts-17961-compliance/audit/cppcheck-post-remediation.txt` (VR-001, SC-002)
- [X] T017 [P] Verify `specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md` has exactly 22 rows, zero FAIL verdicts, every N-A row contains an inline grep command and zero-result output, every PASS row cites a specific FR or grep-based proof, and no cell is empty or contains "TBD" (VR-003, VR-004, FR-001, FR-004)
- [X] T018 Run `bash tests/test_ts17961_enforcement.sh` and confirm exit 0 (synthetic violation caught and reverted), then run `tests/dw-tests.sh` and confirm all tests pass; record both results in `specs/014-iso-ts-17961-compliance/acceptance-report.md` FR-009 section (SC-003, SC-004)
- [X] T020 [US4] Conduct SC-005 reviewer-experience check: provide `specs/014-iso-ts-17961-compliance/` (compliance matrix + tooling guide + acceptance report) to a second reviewer familiar with TS 17961 but unfamiliar with dwsolver; the reviewer must be able to navigate from any compliance matrix row to its supporting evidence without asking clarifying questions; record the outcome (pass/fail and any clarifying questions raised) in the acceptance report FR-008 section as the SC-005 evidence (**U-003**: SC-005 now has an implementation task)

---

## Dependencies & Execution Order

### Phase Dependencies

(Note: Phase 2 is intentionally omitted — research.md and data-model.md were authored as part of the planning workflow; see design note at top of this file.)

- **Setup (Phase 1)**: No dependencies — run T001 and T002 in parallel immediately
- **User Story 1 (Phase 3)**: Depends on Phase 1 (baseline needed for initial FAIL evidence); T003 must complete before T004/T005/T006 (document header first); T004 and T005 can run in parallel
- **User Story 2 (Phase 4)**: Depends on Phase 1 (pre-remediation baseline); T007 and T008 can run in parallel (different call sites but same file — sequence carefully); T009 depends on T007 and T008; T010 depends on T009
- **User Story 3 (Phase 5)**: Depends on Phase 1 (pre-remediation baseline for tooling guide diff); T011 depends on Phase 1 + T010 (post-remediation clean output); T012 and T013 can run in parallel with T011; T019 (suppressions file) must complete before T014
- **User Story 4 (Phase 6)**: Depends on US1 (T006 → T009), US2 (T010), US3 (T014); T020 (SC-005 reviewer check) should be the final task before closing the spec
- **Polish (Final)**: Depends on all story phases complete

### User Story Dependencies

- **US1 (P1)**: Starts after Phase 1 — no dependency on US2 or US3 for initial FAIL rows; T009 (FAIL→PASS update) depends on US2 T007/T008
- **US2 (P2)**: Starts after Phase 1 — independent of US1 document authoring; remediation tasks target source files only
- **US3 (P2)**: Starts after Phase 1 — T011 should start after T010 (post-remediation output needed for before/after diff in tooling guide); T012/T013 are independent of T010
- **US4 (P3)**: Starts after US1 + US2 + US3 are complete

### Parallel Opportunities

- T001 ∥ T002 (Phase 1 setup)
- T004 ∥ T005 (spec 013 rows ∥ spec 008 rows — same document, different sections)
- T007 ∥ T008 are logically parallelizable but edit the same files — sequence T007 then T008 in a single pass
- T012 ∥ T013 ∥ (partial T011) — tooling guide optional section ∥ enforcement script
- T016 ∥ T017 (final verification runs)

---

## Parallel Example: User Story 2 (Remediation)

```bash
# Terminal 1 (or single pass — same files, sequence carefully):
# T007: ga-buffer fix in src/dw_rounding.c:425,445 and src/dw_support.c:668
# T008: nonnullptr fix in src/dw_rounding.c:427-428,447-448 and src/dw_support.c:673

# Then sequentially:
# T009: update compliance matrix FAIL→PASS
# T010: make && tests/dw-tests.sh
```

---

## Implementation Strategy

**MVP (deliver US1 alone)**: T001 → T002 → T003 → T004 ∥ T005 → T006.
Produces a complete compliance matrix with accurate FAIL entries — sufficient to begin a
procurement conversation even before remediation is complete.

**Full compliance (all stories)**: Complete all phases in priority order (US1 → US2 → US3 → US4).
The acceptance report (T015) is the final deliverable and closes the spec.
