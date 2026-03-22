---
description: "Implementation tasks for callable library interface"
---

# Tasks: Callable Library Interface

**Feature**: `010-callable-library`  
**Branch**: `010-callable-library`  
**Input**: Design documents from `/specs/010-callable-library/`  
**Prerequisites**: plan.md ✅ spec.md ✅ research.md ✅ data-model.md ✅ contracts/api.md ✅ quickstart.md ✅

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no shared dependencies on incomplete tasks)
- **[Story]**: Which user story this task belongs to ([US1], [US2], [US3])
- Exact file paths are included in every description

---

## Phase 1: Setup (Infrastructure Bootstrapping)

**Purpose**: Lay the one-line foundation that all other tasks depend on — the visibility
macro that the public header and build system both reference.

**⚠️ CRITICAL**: No user story work can begin until this phase is complete.

- [X] T001 Add `DWSOLVER_API` visibility macro block to `src/dw.h` immediately after the existing include-guard open: conditional `__declspec(dllexport)` / `__declspec(dllimport)` for `_WIN32` / `__CYGWIN__` (guarded by `DWSOLVER_BUILDING_LIB`), `__attribute__((visibility("default")))` for `__GNUC__` / `__clang__`, and empty fallback `#define DWSOLVER_API` — exact text from `contracts/api.md` Visibility Macro section

**Checkpoint**: `src/dw.h` defines `DWSOLVER_API`. Every downstream file that includes `dw.h` will have visibility control available.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Create the public header that is the ABI boundary between the library
and all consumers. Every user story — and the build system — depends on this file
existing and being correct.

**⚠️ CRITICAL**: No user story work can begin until this phase is complete.

- [X] T002 Create `src/dw_solver.h` — full public API header: include guard `DWSOLVER_H_`; `#include "dw.h"` to pull in `DWSOLVER_API`; `dw_status_t` enum (8 codes: `DW_STATUS_OK=0`, `DW_STATUS_ERR_BAD_ARGS=1`, `DW_STATUS_ERR_FILE=2`, `DW_STATUS_ERR_INFEASIBLE=3`, `DW_STATUS_ERR_PHASE1_LIMIT=4`, `DW_STATUS_ERR_PHASE2_LIMIT=5`, `DW_STATUS_ERR_UNBOUNDED=6`, `DW_STATUS_ERR_INTERNAL=99`); `dw_options_t` struct (12 fields per `data-model.md`: `verbosity`, `mip_gap`, `max_phase1_iterations`, `max_phase2_iterations`, `rounding_flag`, `integerize_flag`, `enforce_sub_integrality`, `print_timing_data`, `print_final_master`, `print_relaxed_sol`, `perturb`, `shift`); `dw_result_t` struct (4 fields: `status`, `objective_value`, `num_vars`, `double* x`); thread-safety warning block (verbatim from `contracts/api.md`); declarations for all 4 public functions decorated with `DWSOLVER_API`; `#ifdef __cplusplus extern "C"` guards

**Checkpoint**: `src/dw_solver.h` compiles cleanly with `cc -c src/dw_solver.h` (empty unit test). All types are defined; all 4 function declarations are present.

---

## Phase 3: User Story 1 — Embed Solver in a Host Application (Priority: P1) 🎯 MVP

**Goal**: A C program can `#include "dwsolver.h"`, link `-ldwsolver`, call `dw_solve()` on the `book_bertsimas` example, and get back the same objective value that the `dwsolver` CLI produces on the same input.

**Independent Test**: Compile and run `tests/test_lib_api.c` against the built `libdwsolver.la`. The test program must: complete without crash; return `DW_STATUS_OK` for a valid solve; return `DW_STATUS_ERR_BAD_ARGS` when called with NULL arguments; verify that `dw_result_free()` zeroes all fields. (Requires Phase 4 build system tasks to link, but the source can be written independently.)

### Implementation for User Story 1

- [X] T003 [P] [US1] Create `src/dw_solver.c` — implement `dw_options_init()`: set all 12 `dw_options_t` fields to their documented defaults per `data-model.md` (verbosity=5, mip_gap=0.01, max_phase1_iterations=100, max_phase2_iterations=3000, rounding_flag=0, integerize_flag=0, enforce_sub_integrality=0, print_timing_data=0, print_final_master=0, print_relaxed_sol=0, perturb=0, shift=0.0); include `"dw_solver.h"` and `"dw.h"`

- [X] T004 [US1] Extract solver core from `src/dw_main.c` into a private static function `dw_solver_run(faux_globals *globals)` in `src/dw_solver.c` — move all code from the `init_signals()` call through the `free_globals()` call (~800 lines); move the private `static int hook(glp_tree *T, void *info)` function to `src/dw_solver.c` as well; `src/dw_main.c` retains only `main()` temporarily (it will be overwritten in Phase 5 — T012); do not change any logic or variable names during the move

- [X] T005 [US1] Implement public `dw_solve()` in `src/dw_solver.c`: validate that `master_file`, `subproblem_files`, and `result` are non-NULL (return `DW_STATUS_ERR_BAD_ARGS` immediately); validate `num_subproblems >= 1`; if `opts` is NULL use a local default from `dw_options_init()`; allocate and populate a `faux_globals` struct from the `dw_options_t` fields (translation layer: map all 12 option fields to the corresponding `faux_globals` members used by the original `main()`); call `dw_solver_run()`; on success populate `result->objective_value`, `result->num_vars`, and `result->x` (heap-allocated copy of the x-vector); set `result->status`; return the status code

- [X] T006 [US1] Implement `dw_result_free()` in `src/dw_solver.c`: if `result` is NULL return immediately (no-op); free `result->x`; set `result->x = NULL`, `result->num_vars = 0`, `result->objective_value = 0.0`, `result->status = 0`

- [X] T007 [US1] Implement `dw_version()` in `src/dw_solver.c`: return a pointer to a static string `"dwsolver 1.2.1"` (version matches `PACKAGE_VERSION` in `configure.ac`)

- [X] T008 [P] [US1] Create `tests/test_lib_api.c` — standalone C99 test program (no test framework): test 1: call `dw_options_init()` and assert all 12 fields match documented defaults; test 2: call `dw_solve()` with `examples/book_bertsimas` master and subproblem paths, assert return is `DW_STATUS_OK`, assert `result.num_vars > 0`, assert `result.objective_value` is finite; test 3: call `dw_result_free()` and assert `result.x == NULL` and `result.num_vars == 0`; test 4: call `dw_solve()` with NULL `master_file` and assert return is `DW_STATUS_ERR_BAD_ARGS` and `result.x == NULL`; test 5: call `dw_version()` and assert return value starts with `"dwsolver"`; exit 0 on all assertions passing, exit 1 on first failure; include paths relative to install or build tree (consult quickstart.md)

**Checkpoint**: All code in `src/dw_solver.c` compiles. `tests/test_lib_api.c` compiles (may not yet link until Phase 4 completes the Makefile). Core solver logic is completely extracted; `src/dw_main.c` still contains its old `main()` (to be replaced in T012).

---

## Phase 4: User Story 2 — Build System Produces Both CLI and Library (Priority: P2)

**Goal**: `./configure && make && make install` on a clean checkout produces `bin/dwsolver`, `lib/libdwsolver.a`, `lib/libdwsolver.so` (or `.dylib`), `include/dwsolver.h`, and `lib/pkgconfig/dwsolver.pc` in the install prefix. The existing regression suite passes without change.

**Independent Test**: After `make install` to a local prefix (`./configure --prefix=$(pwd)/tmp/install && make && make install`), verify: `ls tmp/install/bin/dwsolver` exists; `ls tmp/install/lib/libdwsolver*` lists both `.a` and shared variants; `ls tmp/install/include/dwsolver.h` exists; `PKG_CONFIG_PATH=tmp/install/lib/pkgconfig pkg-config --libs dwsolver` outputs `-ldwsolver`.

### Implementation for User Story 2

- [X] T009 [US2] Update `src/Makefile.am` — complete rewrite of the build rules section: define `GLPK_SOURCES` (all `glp*.c` files currently in the `dwsolver_SOURCES` list); define `DW_CORE_SOURCES = dw_blas.c dw_support.c dw_phases.c dw_subprob.c dw_rounding.c dw_globals.c`; add `lib_LTLIBRARIES = libdwsolver.la` with `libdwsolver_la_SOURCES = $(GLPK_SOURCES) $(DW_CORE_SOURCES) dw_solver.c`, `libdwsolver_la_CPPFLAGS = $(AM_CPPFLAGS) -DDWSOLVER_BUILDING_LIB`, `libdwsolver_la_CFLAGS = $(AM_CFLAGS) -fvisibility=hidden`, `libdwsolver_la_LDFLAGS = -version-info 0:0:0`; change `bin_PROGRAMS = dwsolver` to have `dwsolver_SOURCES = dw_main.c`, `dwsolver_CFLAGS = $(AM_CFLAGS)`, `dwsolver_LDADD = libdwsolver.la`; add `nobase_include_HEADERS = dw_solver.h` (installs as `dwsolver.h` in `$(includedir)`)

- [X] T010 [P] [US2] Create `dwsolver.pc.in` at the repository root with the following content: `prefix=@prefix@`, `exec_prefix=@exec_prefix@`, `libdir=@libdir@`, `includedir=@includedir@`, blank line, `Name: dwsolver`, `Description: Dantzig-Wolfe decomposition LP solver library`, `Version: @PACKAGE_VERSION@`, `Cflags: -I${includedir}`, `Libs: -L${libdir} -ldwsolver`

- [X] T011 [US2] Update `configure.ac` — add `dwsolver.pc` to the `AC_CONFIG_FILES` macro call (alongside the existing `Makefile` entries); add `EXTRA_DIST += dwsolver.pc.in` in top-level `Makefile.am` so the template is distributed in tarballs; no new m4 macros or autoconf requirements

- [X] T012 [US2] Update `tests/Makefile.am` — add `check_PROGRAMS = test_lib_api`; add `test_lib_api_SOURCES = test_lib_api.c`; add `test_lib_api_LDADD = $(top_builddir)/src/libdwsolver.la`; add `test_lib_api_CPPFLAGS = -I$(top_srcdir)/src`; add `TESTS = test_lib_api` (so `make check` runs it); preserve all existing `dw-tests.sh` references unchanged

**Checkpoint**: `autoreconf -fi && ./configure && make` succeeds. `ls src/.libs/libdwsolver.*` shows both `.a` and shared library. `src/dwsolver` binary exists and links. `pkg-config --modversion dwsolver` (with appropriate `PKG_CONFIG_PATH`) returns the version string.

---

## Phase 5: User Story 3 — Configure Solver Options Programmatically (Priority: P3)

**Goal**: A host program can set non-default options (silent verbosity, custom MIP gap, custom iteration limit) in `dw_options_t` before calling `dw_solve()` and the solver honors those settings. Unset fields behave as defaults. The CLI (now refactored) exercises the same option pathway as a library consumer.

**Independent Test**: Extend `tests/test_lib_api.c` with option-control scenarios: call `dw_solve()` with `verbosity=0` and assert no output is written to stdout/stderr; call `dw_solve()` with a custom `mip_gap=0.0001` and assert the solve completes; call `dw_solve()` with `opts=NULL` and confirm defaults are used (solve succeeds). Run `tests/dw-tests.sh` and confirm 100% pass rate with the refactored CLI.

### Implementation for User Story 3

- [X] T013 [US3] Rewrite `src/dw_main.c` as a thin CLI wrapper (~70 lines): `#include "dw_solver.h"`; in `main()`: declare `dw_options_t opts`, `dw_result_t result = {0}`, `const char *master_file = NULL`, `const char **subproblem_files = NULL`, `int num_subproblems = 0`; call `dw_options_init(&opts)`; update `process_cmdline` in `dw_support.c`/`dw_support.h` to accept `dw_options_t *opts` plus out-parameters `const char **master_file_out`, `const char ***subproblem_files_out`, `int *num_subproblems_out` — populating the solver option fields **and** the file-path out-parameters from the same argv parse, mapping the same CLI flags (`-g`, `-v`, `--mip-gap`, `--round`, `--iter1`, `--iter2`, etc.) to `dw_options_t` fields; call `process_cmdline(argc, argv, &opts, &master_file, &subproblem_files, &num_subproblems)` and check for error; call `dw_solve(master_file, (const char *const *)subproblem_files, num_subproblems, &opts, &result)`; call `dw_result_free(&result)`; return the status code; remove all solver logic that was previously in `main()` (it now lives in `dw_solver.c`)

- [X] T014 [P] [US3] Add US3 acceptance-scenario tests to `tests/test_lib_api.c`: test 6: set `opts.verbosity = 0`, redirect stdout/stderr to `/dev/null` (or use `dup2`), call `dw_solve()`, restore file descriptors, assert return is `DW_STATUS_OK`; test 7: set `opts.mip_gap = 0.0001`, call `dw_solve()` on `book_bertsimas`, assert return is `DW_STATUS_OK`; test 8: call `dw_solve()` with `opts = NULL` and assert return is `DW_STATUS_OK` (library uses internal defaults)

**Checkpoint**: `make check` runs `test_lib_api` and all 8 tests pass. `tests/dw-tests.sh` passes 100% using the refactored CLI binary. The refactored `dw_main.c` is ≤ 80 lines.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Verify all success criteria from the spec (SC-001 through SC-006) are met. No new code — validation and cleanup only.

- [X] T015 [P] Regenerate autoconf artifacts — run `autoreconf -fi` from the repository root to rebuild `configure`, `Makefile.in`, and `src/Makefile.in`; commit the regenerated files

- [X] T016 Verify full clean build — from the repository root run `./configure && make && make check`; confirm: exit 0, no warnings treated as errors, `src/.libs/libdwsolver.a` exists, `src/.libs/libdwsolver.so` or `.dylib` exists, `src/dwsolver` exists (SC-001)

- [X] T017 [P] Verify symbol visibility — run `nm -gD src/.libs/libdwsolver.so` (Linux) or `nm -gU src/.libs/libdwsolver.dylib` (macOS) and confirm the only exported symbols matching `dw_*` are exactly: `dw_options_init`, `dw_solve`, `dw_result_free`, `dw_version`; confirm zero `glp_*`, `_glp*`, or `faux_*` symbols appear in the export list (SC-004)

- [X] T018 Run full regression suite — `tests/dw-tests.sh` must exit 0 with all previously passing tests still passing; do not modify the test script (SC-003)

- [X] T019 [P] Verify pkg-config integration — install to a local prefix (`make install prefix=$(pwd)/tmp/install`); run `PKG_CONFIG_PATH=$(pwd)/tmp/install/lib/pkgconfig pkg-config --cflags --libs dwsolver`; confirm output contains `-I…/include` and `-ldwsolver` (SC-006)

- [X] T020 [P] Verify memory cleanliness — compile a minimal repeat-solve loop from `quickstart.md` scenario 2 (10 consecutive solves with `dw_result_free()` after each); run under `valgrind --leak-check=full` (Linux) or `leaks` (macOS); confirm zero heap leaks attributed to the library (SC-005)

**Checkpoint**: All 6 success criteria (SC-001 through SC-006) are met. The feature is ready for review.

---

## Dependencies & Execution Order

### Phase Dependencies

```
Phase 1 (Setup: T001)
    └── Phase 2 (Foundational: T002)
            ├── Phase 3 (US1: T003–T008)   ← can begin after T002
            ├── Phase 4 (US2: T009–T012)   ← can begin after T002
            └── Phase 5 (US3: T013–T014)   ← must follow Phase 3 + Phase 4
                    └── Phase 6 (Polish: T015–T020)
```

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2)**: Depends on T001 — BLOCKS all user stories
- **Phase 3 (US1)** and **Phase 4 (US2)**: Both depend on T002; are independent of each other and can proceed in parallel
- **Phase 5 (US3)**: Depends on Phase 3 complete (`dw_solver.c` done) AND Phase 4 complete (build system link works)
- **Phase 6 (Polish)**: Depends on all prior phases complete

### Within-Phase Parallelism

**Phase 3** parallel pairs (different files, no shared incomplete dependencies):
- T003 (dw_solver.c stub), T007 (dw_version), and T008 (test_lib_api.c skeleton) can all start in parallel after T002
- T004 (extract solver core) must precede T005 (dw_solve relies on dw_solver_run)
- T006 and T007 depend on T003 creating `src/dw_solver.c` first; they are NOT independently parallel

**Phase 4** parallel pairs:
- T010 (dwsolver.pc.in) is fully independent and can run alongside T009/T011/T012

**Phase 6** parallel pairs:
- T015, T017, T019, T020 are all read/verify steps that can run in parallel after T016

### User Story Dependencies

- **US1 (P1)**: Depends on Phase 2 only
- **US2 (P2)**: Depends on Phase 2 only (can proceed in parallel with US1)
- **US3 (P3)**: Depends on US1 complete and US2 complete (needs solver logic + working build to test)

---

## Parallel Execution Examples

### Scenario A: Single Developer (recommended sequential order)
```
T001 → T002 → T003 → T004 → T005 → T006 → T007 → T008
     → T009 → T010 → T011 → T012
→ T013 → T014 → T015 → T016 → T017 → T018 → T019 → T020
```

### Scenario B: Two Developers (split US1 and US2)
```
Dev 1: T001 → T002 → T003 → T004 → T005 → T006 → T007 → T008
Dev 2:                T009 → T010 → T011 → T012
Both: → T013 → T014 → T015 → T016 → T017 → T018 → T019 → T020
```

### Scenario C: Fast parallel within US1
```
After T002:
  Thread A: T003 → T004 → T005
  Thread B: T006 (dw_result_free, independent)
  Thread C: T007 (dw_version, independent)
  Thread D: T008 (test_lib_api.c skeleton, independent)
Sync: T005 + T006 + T007 complete → T008 finalized
```

---

## Implementation Strategy

### MVP Scope

**Just Phase 3 + Phase 4** delivers a working, linkable library. A host program can solve
LP problems via the API, the build produces both artifacts, and the CLI regression suite
passes. This is the minimum deliverable that satisfies the feature's primary value proposition.

### Incremental Delivery

1. **T001 + T002** (< 1 hour): Visibility macro + public header. Zero risk — additive only.
2. **T003 + T004 + T005 + T006 + T007** (core of US1): Extract and expose the solver. Highest complexity task; keep diff minimal — move code verbatim, do not refactor.
3. **T008 + T009 + T010 + T011 + T012** (US2): Wire the build system. Mechanical but must be correct.
4. **T013 + T014** (US3): Thin CLI wrapper is the last source change. Tests confirm option routing.
5. **T015–T020** (Polish): Validation only.

### Key Risk: Solver Core Extraction (T004)

The extraction of `dw_solver_run()` from `dw_main.c` is the highest-risk task. Mitigations:
- Move code **verbatim** — no renaming, no logic changes during the move
- Verify intermediate build compiles after T004 before proceeding to T005
- Use `diff` to confirm the moved lines match the original exactly
- Run `tests/dw-tests.sh` after T012 (CLI refactor) as the first validation gate
