# Feature Specification: Cross-Platform Build Repair & Safety Hardening

**Feature Branch**: `002-cross-platform-repair`
**Created**: 2026-03-19
**Status**: Draft
**Input**: Known Defects KD-001, KD-004, KD-005, KD-009 from `specs/001-as-built/spec.md`; Platform Status gap (Linux broken, Windows unknown)

---

## Summary

The DWSOLVER build is currently broken on Linux due to global variable definitions
placed directly in a shared header (`src/dw.h`). This causes `multiple definition`
linker errors under GCC. The same code compiles on macOS only because Clang applies
a more permissive "common symbol" rule. This repair spec:

1. Fixes the Linux linker failure (KD-001) — the primary blocker.
2. Hardens buffer-handling safety (KD-004: `sprintf` → `snprintf`).
3. Adds `malloc` null-checks at critical allocation sites (KD-005).
4. Extends the test suite to cover the two uncovered examples with objective-value
   assertions (KD-009).

Windows support is tracked separately as a future effort and is **out of scope**
for this spec.

---

## User Scenarios & Testing

### User Story 1 — Build and pass tests on Linux/GCC (Priority: P1)

A developer on a Linux system (Ubuntu or similar) runs the standard
`./configure --enable-named-semaphores && make` and gets a successful build.
Running `tests/dw-tests.sh` passes all four deterministic examples. This is
the direct fix for KD-001 and validates its correctness.

**Why this priority**: The software is currently completely unusable on Linux,
which is a primary target platform. Nothing else matters if this is broken.

**Independent Test**: Run `./configure && make` on Linux (or in a Linux
Docker container / CI runner). Zero linker errors. Then run
`tests/dw-tests.sh` — all four tests pass.

**Acceptance Scenarios**:

1. **Given** a clean Linux build environment with GCC,
   **When** `./configure && make` is run (no `--enable-named-semaphores` needed on Linux; that flag is macOS-only),
   **Then** the build completes with zero errors (warnings from GLPK vendor
   code are acceptable).

2. **Given** a successful Linux build,
   **When** `tests/dw-tests.sh` is run with `dwsolver` on `$PATH`,
   **Then** all four deterministic tests (`book_bertsimas`, `book_lasdon`,
   `web_mitchell`, `web_trick`) produce `PASS`.

3. **Given** the same fix applied to a macOS build,
   **When** `./configure --enable-named-semaphores && make` is run on macOS,
   **Then** the build still succeeds and all tests still pass (no regression).

---

### User Story 2 — Test suite covers all deterministic examples (Priority: P2)

The existing `tests/dw-tests.sh` covers four examples but misses `four_sea`
and `book_dantzig`. Both have deterministic optimal *values* even though
variable assignments are non-deterministic due to thread scheduling. The
test suite should be extended to assert the correct objective value for
these two examples without requiring a specific variable-assignment match.

**Why this priority**: Without this, a regression in either problem would go
undetected by automated testing.

**Independent Test**: Run `tests/dw-tests.sh` — the two new tests run, capture
the final `#### Master objective value = ...` line from solver stdout (not from
the solution file, which does not contain the objective), and assert the
expected values: `four_sea` = `1.200000e+01`, `book_dantzig` = `6.357895e+01`.
The test must also assert the solver exited with code 0 (non-zero exit must
not be mistaken for a correct result).

**Acceptance Scenarios**:

1. **Given** a successfully built `dwsolver`,
   **When** `tests/dw-tests.sh` runs the `four_sea` example,
   **Then** the test reports PASS if the objective value equals 12, FAIL otherwise.

2. **Given** a successfully built `dwsolver`,
   **When** `tests/dw-tests.sh` runs the `book_dantzig` example,
   **Then** the test reports PASS if the objective value equals `6.357895e+01`
   (confirmed by two independent reference runs), FAIL otherwise.

---

### User Story 3 — Buffer safety: sprintf replaced with snprintf (Priority: P3)

All `sprintf` calls in DW-authored source files are replaced with `snprintf`
with the correct buffer length argument. Buffer overflow return codes are
checked. This addresses KD-004.

**Why this priority**: Safety correctness. A buffer overflow in a solver
used on large LP files is a crash risk. Lower priority than the build fix
because it does not affect current functionality on macOS.

**Acceptance Scenarios**:

1. **Given** a code review or grep of DW source files,
   **When** searching for `sprintf(` in `src/dw_*.c` and `src/dw_*.h`,
   **Then** zero results are found (all replaced with `snprintf`).

2. **Given** the modified code,
   **When** the full examples test suite is run,
   **Then** all tests still pass (no functional regression from the change).

---

### User Story 4 — malloc null-checks at critical allocation sites (Priority: P3)

The most critical `malloc` calls in DW-authored source files are checked for
NULL return. On allocation failure the program prints an error message to
`stderr` and exits cleanly rather than segfaulting. This addresses KD-005.

**Why this priority**: Correctness and safety, but a `malloc` failure on
modern systems is rare. Lower priority than the build and test coverage work.

**Acceptance Scenarios**:

1. **Given** the modified code,
   **When** the full examples test suite is run,
   **Then** all tests still pass.

2. **Given** a code review,
   **When** examining `malloc` calls for `globals`, `md`, `threads`,
   `sub_data`, and per-thread structures in `dw_main.c` and `dw_support.c`,
   **Then** each critical allocation checks for NULL and handles failure.

---

### Edge Cases

- Fixing `dw.h` definitions must not change the runtime behaviour on macOS or
  in any way alter the synchronization semantics.
- The test script extensions must not falsely pass when the solver exits with
  a non-zero code (i.e., error exit ≠ correct objective).
- `snprintf` replacements must use the correct buffer sizes and must not
  silently truncate strings that were previously written correctly with `sprintf`.

---

## Requirements

### Functional Requirements

- **FR-001**: `src/dw.h` MUST NOT define (i.e., allocate storage for) any
  variable. All variables currently defined there MUST be declared `extern` in
  the header and defined exactly once in a single `.c` file.
- **FR-002**: After the fix, `./configure && make` MUST succeed on Linux with
  GCC with zero linker errors.
- **FR-003**: The macOS build MUST continue to succeed after the fix.
- **FR-004**: All four existing deterministic `tests/dw-tests.sh` tests MUST
  continue to pass after the fix.
- **FR-005**: `tests/dw-tests.sh` MUST be extended to run `four_sea` and
  `book_dantzig` and assert the correct objective values.
- **FR-006**: All `sprintf(` calls in `src/dw_*.c` files MUST be replaced
  with `snprintf(` using the correct buffer size constant.
- **FR-007**: The primary `malloc` calls for `globals`, `md`, `threads`,
  `sub_data`, and per-iteration solution arrays MUST check for NULL and call
  a clean-exit handler on failure.

### Key Entities

- **`src/dw.h`**: The file being repaired. Currently holds definitions of
  `pthread_attr_t attr`, multiple `pthread_mutex_t` instances, semaphore,
  `glp_prob*` pointers, `D_matrix* D`, and `signal_data* signals`.
  Post-fix: declarations only (`extern`).
- **`src/dw_globals.c`** (new file): Will hold the single authoritative
  definition of all variables moved out of `dw.h`.
- **`src/Makefile.am`**: Must be updated to compile `dw_globals.c`.
- **`tests/dw-tests.sh`**: Extended to add two new objective-value tests.

### Out of Scope

- Windows build support (separate future spec).
- KD-002: Unbounded subproblem handling.
- KD-003: Removing the vestigial monolithic file from guide file format.
- KD-006: Service queue scalability beyond `num_clients`.
- KD-007: Static variable in `phase_1_iteration`.
- KD-008: `purge_nonbasics` function.
