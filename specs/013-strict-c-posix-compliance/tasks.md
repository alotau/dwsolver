# Tasks: Strict ISO C and POSIX.1 Compliance

**Input**: Design documents from `specs/013-strict-c-posix-compliance/`
**Prerequisites**: plan.md ✅ spec.md ✅ research.md ✅ data-model.md ✅ contracts/posix-header-audit.md ✅

## Phase Cross-Reference (tasks phases → plan phases)

| Tasks Phase | Plan Phase | Content |
|-------------|------------|---------|
| Phase 1 (Setup) | — (pre-work) | Diagnostic baseline capture (T001) |
| Phase 2 (Foundational) | Plan Phase 1 (partial) | `configure.ac` + `autoheader` (T002–T003) |
| Phase 3 (US1) | Plan Phase 1 | `AM_CFLAGS` hardening; FR-001, FR-003 (T004–T006) |
| Phase 4 (US2) | Plan Phase 1 + Phase 2 | `AM_CPPFLAGS` + `clock_gettime` migration; FR-002, FR-009 (T007–T010) |
| Phase 5 (US3) | Plan Phase 3 (partial) | Extension audit; FR-004 (T011–T012) |
| Phase 6 (US4) | Plan Phase 3 | SEI CERT regression check; FR-007, FR-008 (T013–T015) |
| Final Phase | Plan Phase 4 | Acceptance report + quickstart validation; FR-006, SC-005 (T016–T017) |

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies on incomplete tasks)
- **[Story]**: User story from spec.md ([US1]–[US4])
- Exact file paths included in every task description

---

## Phase 1: Setup

**Purpose**: Capture pre-change diagnostic baseline so Phase 6 diffs have an accurate starting point separate from the feature-008 baseline.

- [ ] T001 Run `gcc -Wall -Wextra -Wsign-compare -c src/dw_*.c -I src/ 2>&1 | grep -v glp` and save output to `specs/013-strict-c-posix-compliance/baseline-pre-warnings.txt` for regression diff reference

---

## Phase 2: Foundational (Blocking Prerequisite)

**Purpose**: Configure `clock_gettime` detection and remove vestigial `gettimeofday` checks in `configure.ac` + `config.h.in`. Must complete before any US1/US2 source changes, because `HAVE_CLOCK_GETTIME` drives the timing guard in `dw_solver.c` and `dw_support.c`.

**⚠️ CRITICAL**: No US2 timing source changes can compile correctly until this phase is complete.

- [ ] T002 Update `configure.ac`: remove the `AC_CHECK_HEADER([sys/time.h])` block (~line 52) and the `AC_CHECK_FUNC([gettimeofday])` block (~line 56); add `AC_SEARCH_LIBS([clock_gettime], [rt])` then `AC_CHECK_FUNC([clock_gettime], AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [Define if clock_gettime is available.]))` after the `AC_CHECK_LIB([pthread])` line in `configure.ac`
- [ ] T003 Regenerate `config.h.in` by running `autoheader` from the project root; verify `HAVE_CLOCK_GETTIME` slot appears and `HAVE_GETTIMEOFDAY` / `HAVE_SYS_TIME_H` slots are removed from `config.h.in`; run `./configure` to regenerate `config.h`

**Checkpoint**: Foundation ready — `HAVE_CLOCK_GETTIME` is defined in `config.h`; US1/US2 source changes may now begin.

---

## Phase 3: User Story 1 — Source Compiles Clean Under ISO C Standard (Priority: P1) 🎯 MVP

**Goal**: Every `src/dw_*.c` compilation is gated by `-std=c11 -pedantic-errors` with zero errors or warnings.

**Independent Test**: `./configure && make` on a clean tree produces zero diagnostics attributable to language-standard violations in any authored `src/dw_*.c` or `src/dw_*.h` file. `make check` passes.

### Implementation for User Story 1

- [ ] T004 [US1] Add `AM_CFLAGS = -std=c11 -pedantic-errors` as the first non-comment line of `src/Makefile.am` (before any target definition); run `./configure && make` to verify zero pedantic errors for all authored source files
- [ ] T005 [P] [US1] Fix any `-pedantic-errors` violations surfaced in `src/dw_*.c` or `src/dw_*.h` by the T004 build (research.md Topic 4 audit found none expected — treat this as the safety-net task; if no violations appear, mark complete after confirming a clean build)
- [ ] T006 [US1] Run `make check` (executes `tests/dw-tests.sh` and `test_lib_api`) to confirm no runtime regression from `-std=c11`; all tests must pass (SC-001, FR-008)

**Checkpoint**: User Story 1 complete — every authored file compiles clean under `-std=c11 -pedantic-errors`. `make check` passes. Delivers FR-001, FR-003.

---

## Phase 4: User Story 2 — POSIX.1-2008 Feature-Test Macro Declared Project-Wide (Priority: P1)

**Goal**: Every translation unit in scope sees `-D_POSIX_C_SOURCE=200809L` before any system header. Timing code migrated to `clock_gettime(CLOCK_MONOTONIC)`.

**Independent Test**: `make` with `_POSIX_C_SOURCE=200809L` active produces zero implicit-declaration diagnostics for any POSIX function. `make check` passes. On a musl-libc host (Alpine Linux), or by manually adding `-D_POSIX_C_SOURCE=200809L -D_XOPEN_SOURCE=700`, zero new diagnostics appear.

### Implementation for User Story 2

- [ ] T007 [US2] Add `AM_CPPFLAGS = -D_POSIX_C_SOURCE=200809L` on the line immediately after `AM_CFLAGS` in `src/Makefile.am`; run `./configure && make` to confirm zero implicit-declaration warnings for POSIX functions across all `src/dw_*.c` (SC-002, FR-002)
- [ ] T008 [US2] Update `src/dw_solver.c` (~lines 118–119): replace `time_t t0 = time(NULL); clock_t c0 = clock();` with an `#ifdef HAVE_CLOCK_GETTIME` / `#else` / `#endif` block — primary path uses `struct timespec dw_ts0; clock_gettime(CLOCK_MONOTONIC, &dw_ts0);`; `#else` retains `time_t t0 = time(NULL); clock_t c0 = clock();`; update any downstream use of `t0` / `c0` within the same function to use `dw_ts0` (or the `#else` variables) under the matching guard (FR-009)
- [ ] T009 [US2] Update `src/dw_support.c::print_timing()` (~lines 828–843): add `#ifdef HAVE_CLOCK_GETTIME` path computing elapsed wall-clock nanoseconds from a `struct timespec` start time passed in or captured inside the function; retain the existing `clock_t`/`time_t` output in the `#else` path; update function declaration in `src/dw_support.h` if the signature changes (FR-009)
- [ ] T010 [US2] Run `./configure && make && make check`; confirm zero implicit-declaration warnings for POSIX functions and all tests pass (SC-002, FR-005, FR-008)

**Checkpoint**: User Story 2 complete — `_POSIX_C_SOURCE=200809L` in effect; `clock_gettime` timing path active where available; `time()`/`clock()` fallback retained. Delivers FR-002, FR-005, FR-009.

---

## Phase 5: User Story 3 — Compiler-Specific Extensions Removed From Authored Source (Priority: P2)

**Goal**: Zero unguarded GCC/Clang extensions remain in `src/dw_*.c` or `src/dw_*.h` outside the existing `#ifdef __GNUC__` and `#if defined(_WIN32)` blocks. All enforced by the `-pedantic-errors` flag added in Phase 3.

**Independent Test**: Build with `-pedantic-errors` (already active from Phase 3) exits zero. A `grep -n '__typeof__\|__builtin_\|__attribute__' src/dw_*.c src/dw_*.h` outside existing guards returns no results.

### Implementation for User Story 3

- [ ] T011 [P] [US3] Audit `src/dw_*.c` and `src/dw_*.h` for any GCC-specific constructs not yet guarded (`__typeof__`, `__builtin_*`, statement expressions, zero-length arrays, VLAs); confirm the existing `-pedantic-errors` build in Phase 3 already rejects all such constructs; update `contracts/posix-header-audit.md` §Extension-Header Disposition with the actual post-Phase-3 audit result
- [ ] T012 [US3] Fix any residual unguarded extension found by T011 in `src/dw_*.c` or `src/dw_*.h` (research.md Topic 4 found none; task exists as a hard gate — if T011 finds zero issues, mark T012 complete immediately)

**Checkpoint**: User Story 3 complete — `-pedantic-errors` build clean; no unguarded extensions in authored source. Confirms FR-003 regression-free after extension removal; Delivers FR-004.

---

## Phase 6: User Story 4 — SEI CERT C Compliance Non-Regression (Priority: P2)

**Goal**: The `cppcheck` and `clang --analyze` output after all changes contains zero additional warnings compared to the feature 008 baselines.

**Independent Test**: `diff specs/008-sei-cert-c-compliance/baseline-cppcheck.txt specs/013-strict-c-posix-compliance/cppcheck-post.txt` shows no new lines (additions). Same for the warning-count diff. `make check` exits 0.

### Implementation for User Story 4

- [ ] T013 [P] [US4] Run `cppcheck --enable=all --suppress=missingIncludeSystem -I src/ src/dw_*.c 2>&1 | tee specs/013-strict-c-posix-compliance/cppcheck-post.txt`; diff against `specs/008-sei-cert-c-compliance/baseline-cppcheck.txt`; zero net-new warnings required to pass (SC-003, FR-007)
- [ ] T014 [P] [US4] Run `gcc -Wall -Wextra -Wsign-compare -std=c11 -D_POSIX_C_SOURCE=200809L -c src/dw_*.c -I src/ 2>&1 | grep -v glp | tee specs/013-strict-c-posix-compliance/warnings-post.txt`; diff against both `specs/013-strict-c-posix-compliance/baseline-pre-warnings.txt` (T001 snapshot — primary gate: zero net-new lines) and `specs/008-sei-cert-c-compliance/baseline-warnings.txt` (FR-007 traceability reference: note any delta attributable to the C11/POSIX flag change vs. substantive code issues); zero net-new diagnostic lines relative to the pre-change snapshot required (SC-003, FR-007)
- [ ] T015 [US4] Run `make check` final verification pass on macOS after all Phase 3–6 changes; confirm all tests pass; confirm the `.github/workflows/ci-linux.yml` GitHub Actions job completes green — satisfying the FR-008 / SC-004 “macOS and Linux” requirement (SC-004, FR-008)

**Checkpoint**: User Story 4 complete — SEI CERT C compliance non-regressed; `make check` passes. Delivers FR-007, FR-008.

---

## Final Phase: Polish & Cross-Cutting Concerns

**Purpose**: Compliance acceptance report and end-to-end validation.

- [x] T016 Write `specs/013-strict-c-posix-compliance/acceptance-report.md` — one section per FR (FR-001 through FR-009) containing: requirement restatement, evidence bullets citing specific files and line ranges, and an explicit ✅ PASS or ❌ FAIL verdict; format must match `specs/008-sei-cert-c-compliance/acceptance-report.md` exactly (FR-006, SC-005)
- [x] T017 [P] Run end-to-end validation per `specs/013-strict-c-posix-compliance/quickstart.md` from a clean build: `./configure && make && make check`; confirm all seven checklist items in quickstart.md produce expected output; fix the quickstart if any step description is inaccurate

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately
- **Phase 2 (Foundational)**: Depends on Phase 1 — BLOCKS Phase 4 timing changes
- **Phase 3 (US1)**: Can begin after Phase 1 (AM_CFLAGS change is independent of configure.ac); best done after Phase 2 for a complete configure run
- **Phase 4 (US2)**: Depends on Phase 2 (needs `HAVE_CLOCK_GETTIME` in `config.h`) and Phase 3 (needs clean pedantic build as starting point)
- **Phase 5 (US3)**: Depends on Phase 3 completion (pedantic errors already enforced)
- **Phase 6 (US4)**: Depends on all previous phases completing
- **Final Phase**: Depends on Phase 6 completion

### User Story Dependencies

- **US1 (Phase 3, P1)**: Depends only on Phase 2 foundation — no dependency on US2
- **US2 (Phase 4, P1)**: Depends on Phase 2 foundation and US1 clean build
- **US3 (Phase 5, P2)**: Depends on US1 (`-pedantic-errors` active from Phase 3)
- **US4 (Phase 6, P2)**: Depends on US1 + US2 + US3 all complete

### Parallel Opportunities Within Phases

- T005 and T006 can proceed in parallel once T004 succeeds
- T011 (audit) can proceed in parallel with T012 (fix gate) since T011 drives T012
- T013 and T014 (different analysis tools) can run in parallel
- T017 (quickstart validation) can run in parallel with T016 (report writing)

---

## Parallel Example: Phase 6 (US4)

```bash
# Launch both analysis tools simultaneously after Phase 5 completes:
Task T013: cppcheck run + diff vs feature-008 baseline
Task T014: gcc -Wall run + diff vs pre-change baseline
# Await both; then run T015 (make check final pass)
```

---

## Implementation Strategy

### MVP (User Stories 1 + 2 Only — Phases 1–4)

1. Phase 1: Capture pre-change baseline
2. Phase 2: Update `configure.ac` and `config.h.in`
3. Phase 3: Add AM_CFLAGS; confirm clean C11 build
4. Phase 4: Add AM_CPPFLAGS + `clock_gettime` timing migration
5. **STOP and VALIDATE**: `make check` passes → FR-001, FR-002, FR-005, FR-008, FR-009 done

### Full Delivery

1. Complete Phases 1–4 (MVP)
2. Phase 5: Audit confirms no unguarded extensions → FR-003, FR-004 done
3. Phase 6: Static analysis non-regression → FR-007 done
4. Final Phase: Write acceptance report → FR-006, SC-005 done

---

## Summary

| Metric | Value |
|--------|-------|
| Total tasks | 17 |
| US1 tasks | 3 (T004–T006) |
| US2 tasks | 4 (T007–T010) |
| US3 tasks | 2 (T011–T012) |
| US4 tasks | 3 (T013–T015) |
| Foundational tasks | 2 (T002–T003) |
| Setup tasks | 1 (T001) |
| Polish tasks | 2 (T016–T017) |
| Parallel opportunities | T005, T011, T013, T014, T017 |
| MVP scope | US1 + US2 (Phases 1–4, tasks T001–T010) |
| Files modified | `configure.ac`, `config.h.in`, `src/Makefile.am`, `src/dw_solver.c`, `src/dw_support.c` (+ `src/dw_support.h` if signature changes) |
