# Tasks: Structured Test Coverage for dwsolver

**Feature**: `006-test-coverage` | **Branch**: `006-test-coverage`
**Input**: [spec.md](spec.md), [plan.md](plan.md), [research.md](research.md), [data-model.md](data-model.md), [quickstart.md](quickstart.md)

## Format: `[ID] [P?] [Story?] Description`

- **[P]**: Can run in parallel (different files, no blocking dependencies)
- **[Story]**: Which user story ([US1]–[US4]) this task belongs to
- Exact file paths are included in every task description

---

## Phase 1: Setup — Pre-condition Fix

**Purpose**: Fix the known pre-existing UB that would make the UBSan CI job fail immediately.
This must land before the sanitizer jobs are added (US2).

- [ ] T001 Fix `printf` format-string UB in `src/dw_support.c`: change `%d` to `%g` (or `%f`) for `DEFAULT_MIP_GAP` in the help/usage output (~line 427); verify with `grep DEFAULT_MIP_GAP src/dw_support.c`

---

## Phase 2: Foundational — Autotools Test Harness

**Purpose**: Wire `tests/Makefile.am` into the Autotools build so `make check` is available.
This phase blocks US1; US2, US3, and US4 can start before it is complete.

**⚠️ Complete T002–T005 before starting Phase 3 (US1)**

- [ ] T002 Add `AC_CONFIG_FILES([tests/Makefile])` to `configure.ac` (after the existing `AC_CONFIG_FILES` line for `src/Makefile`)
- [ ] T003 [P] Add `tests` to `SUBDIRS` in top-level `Makefile.am` (alongside `src`)
- [ ] T004 [P] Create `tests/Makefile.am` with `check_PROGRAMS = test_blas`, `TESTS = test_blas`, `test_blas_SOURCES = test_blas.c`, `test_blas_CFLAGS = -I$(top_srcdir)/src`, `test_blas_LDADD = $(top_builddir)/src/dw_blas.o`
- [ ] T005 Run `autoreconf -fi` from repo root; then `./configure && make` and confirm the build still succeeds with all existing tests passing

**Checkpoint**: `make check` target now exists; `tests/test_blas` will be compiled when `test_blas.c` is created in Phase 3.

---

## Phase 3: User Story 1 — BLAS Unit Tests (Priority: P1) 🎯 MVP

**Goal**: Every function exported from `dw_blas.h` is covered by at least one normal case, one boundary case, and one case with a zero/negative coefficient. `make check` exits 0.

**Independent Test**: `./configure && make check` — passes entirely on its own with no other user stories required.

### Implementation for User Story 1

- [ ] T006 [US1] Create `tests/test_blas.c`: add `#include` block (`stdio.h`, `stdlib.h`, `math.h`, `src/dw_blas.h`), define `EXPECT_APPROX` and `EXPECT_INT` macros exactly as specified in `data-model.md`, and add empty `main()` that prints `ALL TESTS PASSED` and returns 0
- [ ] T007 [US1] Add `dw_daxpy` test cases in `tests/test_blas.c`: normal (len=3, alpha=2.0, verify each element), len=0 (no-op), alpha=0.0 (y unchanged), negative alpha (alpha=-1.0)
- [ ] T008 [US1] Add `dw_ddot` test cases in `tests/test_blas.c`: normal (len=3, known dot product), len=0 (returns 0.0), orthogonal vectors (returns 0.0)
- [ ] T009 [US1] Add `dw_dcopy` test cases in `tests/test_blas.c`: normal (len=3, verify destination matches source), len=0 (no-op, no memory read)
- [ ] T010 [US1] Add `dw_ddoti` test cases in `tests/test_blas.c`: normal (nz=2, known sparse dot product), nz=0 (returns 0.0), single element (nz=1)
- [ ] T011 [US1] Add `dw_dcoogemv` test cases in `tests/test_blas.c`: single non-zero entry (2×2 matrix, verify one output is non-zero and one is zero), full 2×2 dense matrix stored in COO format (verify both output elements)
- [ ] T012 [US1] Run `make check` and confirm exit code 0, `PASS: test_blas` printed, and all assertions pass; introduce a deliberate off-by-one in `dw_daxpy` temporarily to verify the harness catches it and prints the function name and values before reverting

**Checkpoint**: US1 complete — `make check` is a working regression gate for `dw_blas.c`.

---

## Phase 4: User Story 2 — Sanitizer CI Jobs (Priority: P1)

**Goal**: Two new GitHub Actions jobs on `ubuntu-latest` build and run `dw-tests.sh` under ASan+UBSan and TSan respectively.

**Independent Test**: Push branch to GitHub and confirm both new jobs appear in the Actions tab and turn green. Can also be reproduced locally per `quickstart.md` Section 4.

**Depends on**: T001 (UB pre-condition fix must be merged before these jobs are green)

### Implementation for User Story 2

- [ ] T013 [P] [US2] Add `Linux (ASan+UBSan)` job to `.github/workflows/ci-linux.yml`: copy the structure of the existing Linux job (autotools touch steps, `./configure CFLAGS="-fsanitize=address,undefined -g -O1"`, `make`, `export PATH="$PWD/src:$PATH"`, `cd tests && bash dw-tests.sh`); set `continue-on-error: false`
- [ ] T014 [P] [US2] Add `Linux (TSan)` job to `.github/workflows/ci-linux.yml`: same structure as T013 but `CFLAGS="-fsanitize=thread -g -O1"`; TSan and ASan+UBSan must be separate jobs and MUST NOT share a binary
- [ ] T015 [US2] Push branch and verify both sanitizer jobs appear and pass in GitHub Actions; confirm TSan job does not also set ASan flags

**Checkpoint**: US2 complete — any future memory error, UB, or data race in the codebase will block CI.

---

## Phase 5: User Story 3 — New End-to-End Examples (Priority: P2)

**Goal**: Three new example directories under `examples/` each exercise a previously untested code path; all new and existing tests pass in `dw-tests.sh`.

**Independent Test**: `bash tests/dw-tests.sh` (with `dwsolver` on PATH) — each new test block is independently checkable with `pushd examples/<name> && dwsolver guidefile && popd`.

**Depends on**: None — US3 can be worked on in parallel with US1 and US2.

### Implementation for User Story 3 — `examples/single_sub/` (num_clients=1 path)

- [ ] T016 [P] [US3] Design and create `examples/single_sub/sub1.cplex` (single-variable LP: `minimize x; x >= 1`), `examples/single_sub/master.cplex` (coupling bound `x <= 10`), and `examples/single_sub/dw_monolithic.cplex` (union of both); choose variable name consistent across all three files
- [ ] T017 [P] [US3] Create `examples/single_sub/guidefile`: line 1 = `1`, line 2 = `sub1.cplex`, line 3 = `master.cplex`, line 4 = `dw_monolithic.cplex`; run `dwsolver -c examples/single_sub/guidefile | sort > examples/single_sub/ex_relaxed_solution` to capture expected output
- [ ] T018 [US3] Add `single_sub` test block to `tests/dw-tests.sh`: `pushd ../examples/single_sub; dwsolver -c guidefile | sort > /tmp/rs_sorted; diff ex_relaxed_solution /tmp/rs_sorted || exit 1; popd` following the same pattern as existing tests

### Implementation for User Story 3 — `examples/one_iter/` (iteration_count==1 bail path)

- [ ] T019 [P] [US3] Design and create `examples/one_iter/sub1.cplex`, `sub2.cplex`, `master.cplex`, `dw_monolithic.cplex`: craft a 2-subproblem LP whose DW relaxation is already optimal at the initial restricted master (so Phase II terminates on `iteration_count==1`); verify by running `dwsolver -v examples/one_iter/guidefile` and inspecting iteration log
- [ ] T020 [P] [US3] Create `examples/one_iter/guidefile`: line 1 = `2`, lines 2–3 = subproblem filenames, line 4 = master, line 5 = monolithic; note expected obj value from solver output
- [ ] T021 [US3] Add `one_iter` test block to `tests/dw-tests.sh`: run solver and use `grep` to match the known objective value from solver output, following the same grep-based pattern used in existing tests

### Implementation for User Story 3 — `examples/neg_y/` (y-accumulator sign correction path)

- [ ] T022 [P] [US3] Design and create `examples/neg_y/sub1.cplex`, `sub2.cplex`, `master.cplex`, `dw_monolithic.cplex`: craft a 2-subproblem LP with at least one `GLP_LO` (≥) coupling constraint where the y-accumulator sum exceeds the row lower bound during Phase I (triggering `val_local[1] = -1.0`); verify by running solver with debug/verbose output
- [ ] T023 [P] [US3] Create `examples/neg_y/guidefile` and capture `ex_relaxed_solution`: run `dwsolver -c examples/neg_y/guidefile | sort > examples/neg_y/ex_relaxed_solution`
- [ ] T024 [US3] Add `neg_y` test block to `tests/dw-tests.sh` with diff-based check against `ex_relaxed_solution`
- [ ] T025 [US3] Run `bash tests/dw-tests.sh` (with `src/` on PATH) and confirm all 9 tests pass (6 existing + 3 new); fix any example LP that fails to reproduce the correct objective

**Checkpoint**: US3 complete — 3 previously untested code paths each have a named, runnable regression test.

---

## Phase 6: User Story 4 — Guidefile Parsing Tests (Priority: P3)

**Goal**: `tests/test_guidefile.sh` catches guidefile parse errors and edge cases without requiring a working LP solve.

**Independent Test**: `cd tests && bash test_guidefile.sh` — passes entirely on its own once `dwsolver` is built.

**Note**: Per research.md Decision 3, `process_cmdline` cannot be unit-tested without GLPK; these are CLI-level shell tests.

### Implementation for User Story 4

- [ ] T026 [US4] Create `tests/test_guidefile.sh`: add shebang, `set -e` guard disabled (must handle failures), helper `assert_exit` function that checks exit code and prints PASS/FAIL with test name, and `PASS_COUNT`/`FAIL_COUNT` counters; binary assumed to be on PATH
- [ ] T027 [US4] Add "valid 1-subproblem" test case in `tests/test_guidefile.sh`: run `dwsolver examples/single_sub/guidefile`; assert exit code is 0 or 1 (not ≥ 128, which would indicate a crash); assert no `USAGE:` text appears in stderr
- [ ] T028 [US4] Add "valid 4-subproblem" test case in `tests/test_guidefile.sh`: run an existing 4-subproblem example (e.g. the fourth existing example from `dw-tests.sh`); assert exit code ≠ 128+ and no parse-error text
- [ ] T029 [US4] Create `tests/fixtures/bad_count.guidefile` (n=3 declared, only 2 filenames listed) and add test case in `tests/test_guidefile.sh` that runs `dwsolver tests/fixtures/bad_count.guidefile`; assert non-zero exit and that an error message (not silence) appears on stderr
- [ ] T030 [US4] Create `tests/fixtures/missing_file.guidefile` (valid format but references `nonexistent_abc123.cplex`) and add test case in `tests/test_guidefile.sh`; assert non-zero exit and no SIGSEGV (`$?` must not be ≥ 128)
- [ ] T031 [US4] Print test summary at end of `tests/test_guidefile.sh` (`N passed, M failed`); exit 1 if any failures; run `bash tests/test_guidefile.sh` locally and confirm all cases produce expected results

**Checkpoint**: US4 complete — guidefile parse errors and edge-case counts each have a named, runnable check.

---

## Phase 7: Polish & Cross-Cutting Concerns

- [ ] T032 [P] Add `README` to `examples/single_sub/`: describe the problem, known optimal, and which code path it exercises (`num_clients=1` semaphore loop)
- [ ] T033 [P] Add `README` to `examples/one_iter/`: describe the problem, known optimal, and which code path it exercises (`iteration_count==1` bail in `phase_2_iteration`)
- [ ] T034 [P] Add `README` to `examples/neg_y/`: describe the problem, known optimal, and which code path it exercises (y-accumulator sign correction in `phase_1_iteration`)
- [ ] T035 Run the full validation: `make check`, `bash tests/dw-tests.sh`, `bash tests/test_guidefile.sh`; confirm combined run completes in under 2 minutes (SC-006); commit all remaining files

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately
- **Phase 2 (Foundational)**: Depends on Phase 1 completion; blocks Phase 3 (US1) only
- **Phase 3 (US1)**: Depends on Phase 2 (Autotools wiring must be in place)
- **Phase 4 (US2)**: Depends on T001 (UB fix); otherwise independent of Phases 2–3
- **Phase 5 (US3)**: No blocking dependency — can start in parallel with Phases 2–4
- **Phase 6 (US4)**: Depends on T016–T020 (needs `single_sub` and one existing example); otherwise independent
- **Phase 7 (Polish)**: Depends on all story phases being complete

### User Story Dependencies

| Story | Blocks | Blocked by |
|-------|--------|------------|
| US1 (BLAS unit tests) | Nothing | Phase 2 (Autotools) |
| US2 (Sanitizer CI) | Nothing | T001 (UB fix) |
| US3 (New examples) | US4 (needs `single_sub`) | Nothing |
| US4 (Guidefile tests) | Nothing | T017 (single_sub example) |

### Within Each Story

- T006 (skeleton) → T007–T011 (test cases) → T012 (verify): sequential within US1 (same file)
- T013 and T014 ([P]): separate CI job stanzas — write them in the same YAML edit
- T016, T017, T019, T020, T022, T023 ([P]): different example directories — genuinely parallel
- T018, T021, T024 (dw-tests.sh edits): sequential — same file

---

## Parallel Execution Examples

### US3: All Three Example Directories

These three groups can be worked concurrently (different directories):

```bash
# Group A: single_sub
Task T016: Create examples/single_sub/ LP files
Task T017: Capture examples/single_sub/ex_relaxed_solution

# Group B: one_iter (parallel with Group A)
Task T019: Create examples/one_iter/ LP files
Task T020: Create examples/one_iter/guidefile

# Group C: neg_y (parallel with Groups A and B)
Task T022: Create examples/neg_y/ LP files
Task T023: Capture examples/neg_y/ex_relaxed_solution
```

Then sequentially add each to `dw-tests.sh` (T018, T021, T024).

### US2: Both Sanitizer Jobs

```bash
# Both CI job stanzas can be written in one editor pass:
Task T013: Linux (ASan+UBSan) job in ci-linux.yml
Task T014: Linux (TSan) job in ci-linux.yml
```

---

## Implementation Strategy

### MVP First (US1 + US2 — highest value, P1 stories)

1. **Phase 1**: Fix UB pre-condition (T001)
2. **Phase 2**: Wire Autotools (T002–T005)
3. **Phase 3**: Write BLAS unit tests (T006–T012) → `make check` green
4. **Phase 4**: Add sanitizer CI jobs (T013–T015) → CI matrix extended
5. **Stop and validate**: All P1 deliverables complete; branch ready for review

### Incremental Delivery

1. Phases 1–2 → Foundation ready
2. Phase 3 (US1) → `make check` gate live (MVP!)
3. Phase 4 (US2) → Sanitizer gate live
4. Phase 5 (US3) → 3 new code paths covered
5. Phase 6 (US4) → Guidefile errors caught
6. Phase 7 → Polish + full validation

---

## Summary

| Phase | Story | Tasks | Files Changed |
|-------|-------|-------|---------------|
| 1 – Setup | pre-condition | T001 | `src/dw_support.c` |
| 2 – Foundational | build infra | T002–T005 | `configure.ac`, `Makefile.am`, `tests/Makefile.am` |
| 3 – US1 | BLAS tests | T006–T012 | `tests/test_blas.c` |
| 4 – US2 | Sanitizer CI | T013–T015 | `.github/workflows/ci-linux.yml` |
| 5 – US3 | New examples | T016–T025 | `examples/single_sub/`, `examples/one_iter/`, `examples/neg_y/`, `tests/dw-tests.sh` |
| 6 – US4 | Guidefile tests | T026–T031 | `tests/test_guidefile.sh`, `tests/fixtures/` |
| 7 – Polish | cross-cutting | T032–T035 | `examples/*/README` |

**Total tasks**: 35  
**Parallel opportunities**: T003/T004 (Phase 2); T013/T014 (US2); T016/T017/T019/T020/T022/T023 (US3 examples)  
**Suggested MVP scope**: Phases 1–4 (US1 + US2, T001–T015)
