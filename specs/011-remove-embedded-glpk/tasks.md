# Tasks: Remove Embedded GLPK Source (011)

**Input**: Design documents from `specs/011-remove-embedded-glpk/`
**Prerequisites**: plan.md ✅, spec.md ✅, research.md ✅, data-model.md ✅, contracts/build-system.md ✅, quickstart.md ✅

**Tests**: No new tests are written (not requested). User Story 2 validates correctness by running the existing `tests/dw-tests.sh` suite.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies on incomplete tasks)
- **[Story]**: Maps to user story from spec.md (US1, US2, US3, US4)
- Exact file paths are included in every task description

---

## Phase 1: Setup

**Purpose**: Establish a pre-change baseline and confirm local development environment.

- [X] T001 Run `PATH="$PWD/src:$PATH" ./tests/dw-tests.sh 2>&1 | tee specs/011-remove-embedded-glpk/baseline-tests.txt` and record pass/fail counts for each example
- [X] T002 Run `pkg-config --modversion glpk` to confirm GLPK ≥ 4.65 is installed locally; if absent, install via `brew install glpk` (macOS) or `sudo apt-get install libglpk-dev` (Linux)

---

## Phase 2: Foundational — Replace Deprecated `lpx_*` API

**Purpose**: Migrate all `lpx_*` call sites to their modern `glp_*` equivalents **before** the embedded headers are removed. This phase blocks all user story work: the build cannot succeed against system GLPK until these call sites are updated.

**⚠️ CRITICAL**: US1 cannot compile against system GLPK until this phase is complete (modern GLPK 4.65 ships no `lpx_*` declarations).

- [X] T003 [P] [US3] In `src/dw_solver.c`: replace `lpx_read_cpxlp(globals->master_name)` at ~line 184 with `glp_create_prob()` + `glp_read_lp(…, NULL, globals->master_name)` with NULL-check-and-delete-on-failure; and replace all four `lpx_write_cpxlp(master_lp, …)` calls at ~lines 434, 477, 493, 645 with `glp_write_lp(master_lp, NULL, …)` — see data-model.md §B3 for exact before/after
- [X] T004 [P] [US3] In `src/dw_subprob.c`: replace `lp = lpx_read_cpxlp(my_data->infile_name)` at ~line 103 with `lp = glp_create_prob()` + `glp_read_lp(lp, NULL, my_data->infile_name)` with NULL-check-and-delete-on-failure — see data-model.md §B4
- [X] T005 [P] [US3] In `src/dw_phases.c`: replace both occurrences of the `lpx_get_int_parm(master_lp, LPX_K_ITCNT)` pair at ~lines 228–229 and ~lines 460–461 with `glp_get_it_cnt(master_lp)` — see data-model.md §B5

**Checkpoint**: After T003–T005, run `grep -r 'lpx_' src/dw_*.c` — must return no output (SC-004). Commit: `refactor(011): replace deprecated lpx_* API with glp_* equivalents`

---

## Phase 3: User Story 1 — Build Against System GLPK (Priority: P1) 🎯 MVP

**Goal**: Remove embedded GLPK sources and wire the build system to detect and link against a system-installed GLPK ≥ 4.65.

**Independent Test**: A fresh build with no files in `src/glp*`, `src/amd/`, or `src/colamd/` succeeds when GLPK ≥ 4.65 is present, and `configure` fails with a clear error when GLPK is absent.

- [X] T006 [P] [US1] Update `configure.ac`: add `PKG_CHECK_MODULES([GLPK], [glpk >= 4.65], [], [AC_MSG_ERROR([…])])` after the `AC_CHECK_LIB([pthread], …)` line; remove the five `AC_ARG_WITH`/`AC_ARG_ENABLE` blocks for gmp, zlib, dl, odbc, and mysql; remove their corresponding `if test …` application blocks — see contracts/build-system.md §1 for the complete before/after specification
- [X] T007 [P] [US1] Update `src/Makefile.am`: delete the entire `GLPK_SOURCES` variable block (lines 7–113); remove `$(GLPK_SOURCES)` from `libdwsolver_la_SOURCES`; append `$(GLPK_CFLAGS)` to `libdwsolver_la_CPPFLAGS`; add `libdwsolver_la_LIBADD = $(GLPK_LIBS)` after the CPPFLAGS line; replace the `EXTRA_DIST` entry with only the six remaining `dw_*.h` headers — see contracts/build-system.md §2
- [X] T008 [P] [US1] Update `dwsolver.pc.in`: add the line `Requires.private: glpk` between the `Version:` field and the `Cflags:` field — see contracts/build-system.md §3
- [X] T009 [US1] Delete all 141 embedded GLPK source files: run `git rm src/glp*.c src/glp*.h src/amd/amd_*.c src/amd/amd.h src/amd/amd_internal.h src/colamd/colamd.c src/colamd/colamd.h` to stage all deletions (depends on T007 completing first so Makefile.am no longer references them)
- [X] T010 [US1] Run `autoreconf -fi && ./configure && make -j$(nproc || sysctl -n hw.ncpu)` in the project root; confirm build succeeds with zero errors and produces `src/dwsolver` and `src/.libs/libdwsolver.so` (or `.dylib` on macOS) — depends on T006, T007, T008, T009
- [X] T011 [US1] Verify acceptance scenario 2: temporarily unset `PKG_CONFIG_PATH` (or rename `glpk.pc`) and confirm `./configure` exits with a non-zero status and prints the GLPK-missing error message from `PKG_CHECK_MODULES`; then restore the environment (depends on T010)

**Checkpoint**: At this point, User Story 1 is fully functional — a working binary exists, built against system GLPK, with no embedded sources. Commit: `feat(011): remove embedded GLPK; build against system glpk >= 4.65`

---

## Phase 4: User Story 2 — Existing Test Suite Passes (Priority: P1)

**Goal**: Confirm that the refactored binary produces numerically correct results and all previously-passing tests continue to pass.

**Independent Test**: Run `tests/dw-tests.sh`; compare pass/fail counts against `specs/011-remove-embedded-glpk/baseline-tests.txt` (captured in T001); objective values must agree to ≤ 1e-6 relative difference.

- [X] T012 [US2] Run `PATH="$PWD/src:$PATH" ./tests/dw-tests.sh 2>&1 | tee specs/011-remove-embedded-glpk/post-refactor-tests.txt` and compare pass/fail counts to the baseline captured in T001 — all tests that passed before must still pass (SC-003)
- [X] T013 [US2] For any test that numerically compares objective values, confirm the refactored binary's output agrees with the baseline output to within a 1e-6 relative difference; if any value diverges beyond that tolerance, investigate and resolve before proceeding

**Checkpoint**: US2 is complete when all baseline-passing tests pass with the refactored binary. Commit: `test(011): capture post-refactor test results; all tests pass`

---

## Phase 5: User Story 4 — Repository Source Tree Is Clean (Priority: P3)

**Goal**: Confirm the `src/` directory contains only dwsolver's own files with no GLPK remnants.

**Independent Test**: `git ls-files src/ | grep -E 'glp|/amd/|/colamd/'` returns no output.

- [X] T014 [P] [US4] Run `git ls-files src/ | grep -E 'glp|amd|colamd'` and confirm it returns no output (SC-001: all 141 embedded files gone)
- [X] T015 [P] [US4] Run `git ls-files src/ | wc -l` and confirm the count is ≤ 18 (pre-change was 159; 159 − 141 = 18), satisfying SC-001's "shrinks by at least 95 source files"
- [X] T016 [P] [US4] Confirm `third-party/glpk/` is intact: run `git ls-files third-party/glpk/` and verify `glpk-4.44.ThreadReady.patch` and related attribution files are still present (FR-006)

**Checkpoint**: US4 complete — source tree is provably clean. Commit: `chore(011): confirm clean source tree; third-party provenance intact`

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Update CI workflows, Dockerfile, and README to reflect the new external dependency.

- [X] T017 [P] Update `.github/workflows/ci-linux.yml`: add a `- name: Install dependencies` step with `sudo apt-get update && sudo apt-get install -y libglpk-dev` before the `Build` step in **all three** jobs (`linux`, `linux-asan-ubsan`, `linux-tsan`) — see contracts/build-system.md §5a
- [X] T018 [P] Update `.github/workflows/ci-macos.yml`: add a `- name: Install dependencies` step with `brew install glpk` before the `Build` step — see contracts/build-system.md §5b
- [X] T019 [P] Update `Dockerfile`: add `libglpk-dev` to the builder stage's `apt-get install` line; add `RUN apt-get update && apt-get install -y --no-install-recommends libglpk40 && rm -rf /var/lib/apt/lists/*` in the runner stage — see contracts/build-system.md §4
- [X] T020 [P] Update `README.md`: add a **Dependencies** section listing GLPK ≥ 4.65 as required, with install commands for macOS (`brew install glpk`), Ubuntu/Debian (`sudo apt-get install libglpk-dev`), and Fedora/RHEL (`sudo dnf install glpk-devel`); add the ABI note that callers must link against the same GLPK build (FR-007, FR-008a) — see quickstart.md §1 for reference wording
- [X] T021 Run `grep -r 'lpx_' src/dw_*.c` — must return no output (final SC-004 verification)
- [X] T022 Run `docker build -t dwsolver-test .` and confirm it succeeds (SC-007) — depends on T019
- [X] T023 Push the branch and confirm all GitHub Actions CI jobs pass green (SC-006) — depends on T017, T018

**Checkpoint**: All polish tasks complete. Commit: `chore(011): update CI, Dockerfile, README for system GLPK dependency`

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately
- **Phase 2 (Foundational)**: Depends on Phase 1; T003, T004, T005 run in parallel with each other; **BLOCKS** Phase 3
- **Phase 3 (US1)**: T006, T007, T008 run in parallel; T009 depends on T007; T010 depends on T006+T007+T008+T009; T011 depends on T010
- **Phase 4 (US2)**: Depends on Phase 3 (T010 must produce a working binary)
- **Phase 5 (US4)**: Depends on Phase 3 (deletion in T009 must be complete)
- **Phase 6 (Polish)**: T017, T018, T019, T020 run in parallel; T022 depends on T019; T023 depends on T017+T018

### User Story Dependencies

- **US3 (Phase 2)**: Must complete before US1 — provides the code changes that allow the build to link against modern GLPK
- **US1 (Phase 3)**: Must complete before US2 — provides the working binary
- **US2 (Phase 4)**: Independent of US4; both depend on US1
- **US4 (Phase 5)**: Verification-only; depends on US1 for the deletion step

### Parallel Opportunities

```
# Phase 2 — all three files are independent:
T003 src/dw_solver.c   ─┐
T004 src/dw_subprob.c  ─┤─► all can run simultaneously
T005 src/dw_phases.c   ─┘

# Phase 3 — first three files are independent; T009 waits for T007:
T006 configure.ac      ─┐
T007 src/Makefile.am   ─┤─► parallel ─► T009 delete files ─► T010 build
T008 dwsolver.pc.in    ─┘

# Phase 6 — all documentation files are independent:
T017 ci-linux.yml      ─┐
T018 ci-macos.yml      ─┤─► parallel ─► T023 push & CI
T019 Dockerfile        ─┤
T020 README.md         ─┘
```

---

## Implementation Strategy

### MVP First (Phases 1–4, US1 + US2 only)

1. Complete Phase 1: Baseline capture
2. Complete Phase 2: API migration (foundational)
3. Complete Phase 3: Build system changes + file deletion (US1)
4. **STOP and VALIDATE**: Confirm binary builds and tests pass (US2)
5. Branch is functionally complete at this point

### Incremental Delivery

1. Phase 1 → Phase 2 → Phase 3 → Phase 4 = working, tested binary (MVP)
2. Phase 5 = Tree-cleanliness verification (low effort; run after Phase 4)
3. Phase 6 = CI/Docker/README = production-ready (required for merge)

### Suggested Commit Sequence

1. `refactor(011): replace deprecated lpx_* API with glp_* equivalents` — after Phase 2 checkpoint
2. `feat(011): remove embedded GLPK; build against system glpk >= 4.65` — after Phase 3 checkpoint
3. `test(011): capture post-refactor test results; all tests pass` — after Phase 4 checkpoint
4. `chore(011): confirm clean source tree; third-party provenance intact` — after Phase 5
5. `chore(011): update CI, Dockerfile, README for system GLPK dependency` — after Phase 6

---

## Notes

- [P] tasks operate on different files — safe to run simultaneously with no merge conflicts
- [US3] tasks are placed in Phase 2 (Foundational) because they are a prerequisite for US1, even though the spec assigns US3 a P2 priority label
- `glpk_mutex` guards in `dw_globals.c`, `dw_support.c`, `dw_rounding.c` are NOT touched — they protect non-reentrant file-I/O state and remain correct with modern GLPK
- `third-party/glpk/` is NOT touched — preserved for provenance
- Windows CI (`ci-windows.yml`) is NOT touched — out of scope per spec assumption
- `ci-docker.yml` is NOT touched — it triggers `docker build .` which the Dockerfile fix (T019) covers automatically
