# Feature Specification: Structured Test Coverage for dwsolver

**Feature Branch**: `006-test-coverage`
**Created**: 2026-03-20
**Status**: Draft

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Unit Tests for dw_blas Pure Functions (Priority: P1)

A developer making changes to the linear algebra routines in `dw_blas.c` can run a fast, self-contained test suite that exercises every function (`dw_daxpy`, `dw_ddot`, `dw_ddoti`, `dw_dcopy`, `dw_dcoogemv`) with normal inputs, boundary inputs (empty vector, single element, zero alpha), and numerically sensitive inputs. The tests produce a pass/fail result with no external dependencies other than the C standard library.

**Why this priority**: These functions underpin all subproblem solves. They are pure — no globals, no GLPK, no threads — making them the easiest and highest-value unit test target. A bug here silently corrupts solutions.

**Independent Test**: Build and run `tests/test_blas` after any change to `dw_blas.c`; a non-zero exit code means a regression was introduced. Delivers value entirely on its own with no other stories implemented.

**Acceptance Scenarios**:

1. **Given** the repository is checked out and built, **When** `make check` or `tests/test_blas` is run, **Then** all assertions pass and the program exits 0.
2. **Given** a deliberate off-by-one error is introduced into `dw_daxpy`, **When** the test suite is run, **Then** at least one test fails with a clear message identifying which function and what inputs were used.
3. **Given** `dw_ddoti` is called with `nz=0` (empty sparse vector), **When** the test runs, **Then** it returns 0.0 without reading any memory.
4. **Given** `dw_dcoogemv` is called with a matrix containing a single non-zero entry, **When** the test runs, **Then** the output vector has the correct single non-zero product.

---

### User Story 2 — Sanitizer CI Jobs (ASan/UBSan/TSan) (Priority: P1)

CI automatically builds and runs the existing end-to-end example suite under AddressSanitizer + UndefinedBehaviorSanitizer and under ThreadSanitizer, in separate jobs. Any sanitizer finding causes the job to fail and blocks merge.

**Why this priority**: The sanitizer builds require no new test cases — they provide immediate coverage of memory safety and data-race issues across all existing examples with a single CI change. Regressions in `dw_phases.c` or `dw_rounding.c` that pass normal tests will likely be caught here.

**Independent Test**: Add the two jobs to `.github/workflows/`. Each job is independently green or failing; neither depends on Story 1 being implemented.

**Acceptance Scenarios**:

1. **Given** the CI workflow runs on a push to any branch, **When** the ASan+UBSan job builds with `-fsanitize=address,undefined` and runs `dw-tests.sh`, **Then** all 6 tests pass with zero sanitizer errors reported.
2. **Given** the CI workflow runs, **When** the TSan job builds with `-fsanitize=thread` and runs `dw-tests.sh`, **Then** all 6 tests pass with zero data-race reports.
3. **Given** a known UB is reintroduced (e.g., the `val=NULL` dereference from the Phase 3 fix), **When** ASan/UBSan CI runs, **Then** the job fails and the sanitizer output names the offending file and line.
4. **Given** a known data race is introduced in a mutex path, **When** TSan CI runs, **Then** the job fails with a race report identifying the two conflicting accesses.

---

### User Story 3 — Additional End-to-End Examples for Untested Paths (Priority: P2)

Three new example problems are added to `examples/` and integrated into `dw-tests.sh`, each targeting a code path not exercised by the current six examples:

- A problem with a single subproblem (`num_clients = 1`)
- A problem that converges after exactly one Phase II iteration (exercises the `iteration_count == 1` bail path)
- A problem where Phase I requires the y-accumulator sign correction (the `first_run` inner branch in `phase_1_iteration`)

**Why this priority**: The current examples all use 2–4 subproblems and multi-iteration convergence. Single-client and one-iteration paths contain conditionals that have never been validated by a test run.

**Independent Test**: Run `bash tests/dw-tests.sh` after adding each example; each new test is independently verifiable.

**Acceptance Scenarios**:

1. **Given** a single-subproblem example is added, **When** `dw-tests.sh` runs it, **Then** it produces the known-correct objective value and exits 0.
2. **Given** a one-Phase-II-iteration example is added, **When** `dw-tests.sh` runs it, **Then** the solver terminates in one Phase II iteration and produces the correct objective value.
3. **Given** a problem requiring y-accumulator sign correction is added, **When** `dw-tests.sh` runs it, **Then** Phase I produces a correct feasible basis and Phase II converges to the correct optimum.
4. **Given** any of the three new examples is run, **When** the solver exits non-zero, **Then** `dw-tests.sh` reports the failure and exits immediately (same pattern as existing tests).

---

### User Story 4 — Guidefile Parsing Tests (Priority: P3)

A test program directly exercises the guidefile parsing logic with synthetic input files covering valid guidefiles, malformed guidefiles, and edge-case counts (single subproblem, mismatched count vs filenames). Each case asserts the correct struct fields are populated or that an appropriate error is produced.

**Why this priority**: Guidefile parsing is the sole interface between the user and the solver. Errors here are silent and produce bad decompositions. The parsing code is stable, so this has lower urgency than the BLAS and sanitizer work.

**Independent Test**: Build and run `tests/test_guidefile` independently of the full solver binary.

**Acceptance Scenarios**:

1. **Given** a valid guidefile with 4 subproblems, **When** the parsing logic processes it, **Then** `num_clients` is 4 and all filename fields match the guidefile entries exactly.
2. **Given** a guidefile where the declared subproblem count does not match the number of filenames listed, **When** parsing runs, **Then** the program reports an error and does not proceed to solve.
3. **Given** a guidefile with the minimum valid inputs (1 subproblem), **When** parsing runs, **Then** `num_clients` is 1 and all fields are populated correctly.
4. **Given** a guidefile referencing a `.cplex` file that does not exist, **When** the solver is launched, **Then** an informative error is produced at startup rather than a crash or silent wrong answer.

---

### Edge Cases

- What happens when `dw_ddot` or `dw_daxpy` is called with `len = 0`?
- How do the phase functions behave when every subproblem offers no improving column on the first Phase II iteration?
- What does the solver do with a guidefile that has Windows-style line endings (`\r\n`)?
- How does the single-subproblem path interact with the semaphore logic (`customers`)?

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: A `tests/test_blas.c` file MUST be buildable by `make check` (or an explicit `make test_blas`) without requiring any library other than libc and the existing Autotools build system.
- **FR-002**: `tests/test_blas.c` MUST cover all five functions in `dw_blas.c` with at least one normal case, one boundary case (`len=0` or `nz=0`), and one case with a negative or zero coefficient.
- **FR-003**: Each test case MUST produce an informative failure message (function name, inputs, expected vs actual) before aborting, so failures are diagnosable without a debugger.
- **FR-004**: The CI workflow MUST include an ASan+UBSan job that compiles with `-fsanitize=address,undefined -g` and runs the full example test suite.
- **FR-005**: The CI workflow MUST include a TSan job that compiles with `-fsanitize=thread -g` and runs the full example test suite.
- **FR-006**: Sanitizer CI jobs MUST use `./configure` + `make` with the appropriate `CFLAGS` so they exercise the same Autotools build path as the main CI job.
- **FR-007**: Three new example directories MUST be added under `examples/`, each containing a `guidefile`, all required `.cplex` input files, and an `ex_*` expected-output file.
- **FR-008**: The existing `tests/dw-tests.sh` MUST be extended to run the three new examples, following the same pass/fail pattern as the current tests.
- **FR-009**: No existing test MUST be removed or weakened by this work.
- **FR-010**: `tests/test_guidefile` MUST be buildable and runnable independently of the full solver binary; it MUST link only the guidefile parsing translation unit(s) and libc.

### Key Entities

- **test_blas**: A standalone C test executable in `tests/` that exercises `dw_blas.c` in isolation with no GLPK or threading dependency.
- **test_guidefile**: A standalone C test executable that exercises guidefile parsing logic in isolation.
- **Sanitizer CI jobs**: Two new GitHub Actions jobs in `.github/workflows/` that set sanitizer `CFLAGS` and run the existing `dw-tests.sh`.
- **Example problem**: A directory under `examples/` containing a guidefile, CPLEX-format LP files, and an expected-output file consumed by `dw-tests.sh`.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: `make check` completes with exit code 0 on a clean checkout; all unit test assertions pass.
- **SC-002**: All 6 existing end-to-end tests continue to pass after this work with no regressions.
- **SC-003**: The ASan+UBSan CI job and the TSan CI job are both green on `main` after merge.
- **SC-004**: At least 3 code paths in `dw_phases.c` that were previously unexercised by any test (single-client path, one-iteration bail, y-accumulator sign correction) are each covered by a named, runnable test.
- **SC-005**: A developer can identify a failing unit test, its inputs, and the expected vs actual value solely from the test output — without attaching a debugger.
- **SC-006**: Total time to run the full test suite (unit tests + all end-to-end examples) stays under 2 minutes on the CI runner.

## Assumptions

- The three new example LP problems can be constructed as small, hand-crafted CPLEX files with known optimal values; they do not need to represent real-world scheduling problems.
- The guidefile parsing logic can be made available as a linkable translation unit without pulling in the entire solver binary; this may require minor refactoring to extract it from `dw_support.c` or `dw_main.c`.
- The existing Autotools build system can be extended with a `check` target by adding entries to `Makefile.am`; no migration to a different build system is required.
- TSan and ASan/UBSan are available in the GCC/Clang versions used by the Linux CI runner on GitHub Actions.
- Apple's Clang on macOS does not support `detect_leaks` (LSan); sanitizer CI jobs run only on the Linux runner.

