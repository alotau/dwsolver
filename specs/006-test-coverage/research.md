# Research: Structured Test Coverage for dwsolver

**Branch**: `006-test-coverage`
**Status**: Complete — all NEEDS CLARIFICATION items resolved

---

## Decision 1: BLAS Unit Test Harness

**Decision**: Hand-rolled `EXPECT_APPROX` / `EXPECT_INT` macros in `tests/test_blas.c`; no third-party
test framework.

**Rationale**: dwsolver already has zero external runtime dependencies beyond embedded GLPK.
Introducing a test framework (Unity, CUnit, cmocka) would add a build dependency that breaks the
cross-platform CI matrix (Linux/macOS/Windows). A hand-rolled macro harness is ~20 lines and
provides the same failure visibility.

**Alternatives considered**:
- Unity (ThrowTheSwitch): single-header, BSD-licensed — rejected because automake `check_PROGRAMS`
  integration requires adding Unity's .c to LDADD, adding a dependency to the tarball.
- CUnit: requires `pkg-config cunit` — not available on the Windows MSYS2 CI runner without extra
  setup; rejected.

---

## Decision 2: `tests/Makefile.am` structure

**Decision**: Add a new `tests/Makefile.am` with:
```makefile
TESTS = test_blas
check_PROGRAMS = test_blas
test_blas_SOURCES = test_blas.c
test_blas_CFLAGS = -I$(top_srcdir)/src
test_blas_LDADD = $(top_builddir)/src/dw_blas.o
```

Wire it via `AC_CONFIG_FILES([tests/Makefile])` in `configure.ac`.

**Rationale**: `check_PROGRAMS` / `TESTS` is the standard automake mechanism for `make check`. Linking
only `dw_blas.o` (not the full `dwsolver` link set) keeps the unit test binary minimal. `dw_blas.c`
has no includes beyond libc/libm — confirmed by inspection.

**Alternatives considered**:
- Adding `test_blas` as a `check_PROGRAMS` in `src/Makefile.am` directly — works but mixes test code
  with library source; rejected for cleanliness.
- Using a top-level `Makefile.am` TESTS variable — requires more configure.ac changes; rejected.

---

## Decision 3: US4 Scope Revision — CLI-level guidefile tests

**Decision**: US4 `test_guidefile.sh` is a shell-based CLI test (not a C unit test of
`process_cmdline`).

**Rationale**: `process_cmdline` lives in `dw_support.c`, which includes `dw.h`, which includes
`<glpk.h>`. The `faux_globals` struct contains `glp_smcp*` and `glp_iocp*` members (from `dw.h` line
~60). Extracting just `process_cmdline` into a unit-testable module would require either:
(a) splitting `dw.h` to isolate guidefile-related declarations, or
(b) linking the full GLPK (allowed — it's embedded — but then we test `process_cmdline` only by
    running the entire solver anyway).

CLI-level shell tests (`tests/test_guidefile.sh`) are lower-effort, lower-maintenance, and still
satisfy the user story of catching guidefile parsing errors.

**Alternatives considered**:
- Unit test with full GLPK link: technically possible; rejected because it bloats the test binary and
  doesn't add correctness coverage over running `dwsolver` itself.
- Refactoring `process_cmdline` into a GLPK-free module: valid long-term; out of scope for this
  feature branch.

---

## Decision 4: Sanitizer CI targets

**Decision**:
- Job 1: `ubuntu-latest`, `CFLAGS="-fsanitize=address,undefined -g -O1"` — catches memory errors and
  undefined behaviour.
- Job 2: `ubuntu-latest`, `CFLAGS="-fsanitize=thread -g -O1"` — catches data races in the DW master
  loop.

**Rationale**: GCC/Clang on `ubuntu-latest` supports all three sanitizers. macOS CI runner uses
Apple Clang which does not support TSan reliably on macOS 13+ for multi-threaded code using named
POSIX semaphores (pthread + `sem_open` interplay). Windows CI uses MSVC; sanitizers are not
supported in the same way. Therefore sanitizer jobs are Linux-only.

**Pre-condition**: A known pre-existing UB must be fixed before enabling UBSan:
- `dw_support.c` ~line 427: `printf("... Default %d.\n", DEFAULT_MIP_GAP)` — `DEFAULT_MIP_GAP` is
  defined as `0.10` (a `double`). The `%d` format specifier is undefined behaviour
  (`-fsanitize=undefined` will fire with `runtime error: ... calling printf with more arguments than
  format specifies` or a format-type mismatch). Fix: change `%d` to `%g` or `%f`.

**Alternatives considered**:
- MemorySanitizer (MSan): requires recompiling all dependencies (including GLPK) with MSan
  instrumentation; impractical for embedded GLPK without patching the build system; rejected.
- Valgrind in CI: available on Linux but slow; useful for local debugging but not as a CI gate
  requirement; deferred to optional recommendation.

---

## Decision 5: Three new example LP problems

### 5a. `examples/single_sub/` — exercises `num_clients=1` semaphore path

**LP**: 1 coupling constraint, 1 variable in the subproblem. The simplest non-trivial LP that runs
through DW decomposition with a single subproblem.

```
minimize   x
subject to x >= 1        (subproblem constraint, GLP_LO)
           x <= 10       (master coupling, bound)
```

Guidefile: `n=1`, one subproblem file, one master file.

**Code path**: `while(count < num_clients)` loop in DW master thread iterates exactly once (one
`sem_post` / `sem_wait` round-trip). Tests the boundary condition that `num_clients=1` does not
deadlock or skip the semaphore post.

### 5b. `examples/one_iter/` — exercises `iteration_count==1` bail

**LP**: A problem that is solved in Phase I (the relaxed LP is optimal at the initial basis). Phase II
starts, `iteration_count` is 1 on the first call, and the bail condition
`if (iteration_count == 1 && rc <= 0)` triggers immediately.

Guidefile: `n=2`, two small but already-feasible subproblems. The DW master LP should be optimal
after the Phase I restricted master, meaning Phase II adds zero columns and terminates.

**Code path**: `phase_2_iteration` line checking `iteration_count == 1`.

### 5c. `examples/neg_y/` — exercises y-accumulator sign correction

**LP**: A problem with at least one GLP_LO (≥) coupling constraint where, during Phase I, the
y-accumulator sum exceeds the row lower bound, so `glp_get_row_lb - y_accumulator[i] < 0.0`.

This triggers the `val_local[1] = -1.0` branch (negative y-variable coefficient) in
`phase_1_iteration`.

Guidefile: `n=2`, with coupling constraint of type `GLP_LO`.

**Code path**: `phase_1_iteration` inner `if (first_run)` block — the y-variable sign branch.
