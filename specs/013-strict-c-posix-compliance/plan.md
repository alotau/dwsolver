# Implementation Plan: Strict ISO C and POSIX.1 Compliance

**Branch**: `013-strict-c-posix-compliance` | **Date**: 2026-03-23 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/013-strict-c-posix-compliance/spec.md`

## Summary

Add explicit ISO C11 (`-std=c11 -pedantic-errors`) and POSIX.1-2008 (`-D_POSIX_C_SOURCE=200809L`)
enforcement to the dwsolver build system via `AM_CFLAGS` and `AM_CPPFLAGS` in `src/Makefile.am`.
Perform a source audit confirming zero unguarded compiler extensions in `src/dw_*.c` and
`src/dw_*.h`. Introduce a `clock_gettime(CLOCK_MONOTONIC)` timing path (guarded by
`HAVE_CLOCK_GETTIME`) in `dw_solver.c` and `dw_support.c`, replacing the current ISO C
`time()`/`clock()` calls; remove the vestigial `HAVE_GETTIMEOFDAY` configure check.
Verify that the full test suite passes and that all SEI CERT C rules from feature 008 are
non-regressed. Deliver a compliance acceptance report following the feature 008 format.

## Technical Context

**Language/Version**: C11 (GCC/Clang on macOS/Linux; MinGW-w64 on Windows via MSYS2)
**Primary Dependencies**: POSIX.1-2008 (`pthreads`, `clock_gettime`); system GLPK ≥ 4.65;
GNU Autotools (`configure.ac`, `src/Makefile.am`)
**Storage**: N/A (solver writes solution files; no database)
**Testing**: Shell-based example runner (`tests/dw-tests.sh`); C unit tests (`test_lib_api`);
static analysis (`cppcheck`, `clang --analyze`) diffed against feature 008 baselines
**Target Platform**: macOS, Linux, Windows (MSYS2/MinGW-w64, `continue-on-error: true` in CI)
**Project Type**: CLI tool / shared library
**Performance Goals**: No solver performance regression; `clock_gettime` timing path adds
negligible overhead (single syscall vs. two)
**Constraints**: All source changes confined to `src/dw_*.c`, `src/dw_*.h`, `configure.ac`,
`src/Makefile.am`; vendored GLPK (`third-party/glpk/`) and MKL conditional blocks must not
be touched; all existing tests must keep passing
**Scale/Scope**: 8 authored `.c` files, 6 internal headers, 1 public header; 2 timing call
sites (`dw_solver.c` ~lines 118–119, `dw_support.c::print_timing()` ~lines 828–835);
0 actual `gettimeofday` call sites (configure check is vestigial — see `research.md` §Topic 3)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Verdict | Notes |
|-----------|---------|-------|
| I. Correctness First | **PASS** | Changes do not touch solver math. Replacing `time()` with `clock_gettime()` in timing reporting is additive; timing data does not affect solver results. The `clock_gettime` fallback retains `time()`/`clock()` for non-POSIX hosts, ensuring bit-identical behaviour on those platforms. |
| II. Thread Safety is Mandatory | **PASS** | No new threading primitives or lock-order changes introduced. `clock_gettime(2)` is thread-safe per POSIX.1-2008 §6.3.1. Adding `_POSIX_C_SOURCE` does not change pthreads semantics. TSan not required — no synchronization primitives, lock ordering, or shared-state access paths are modified by this feature. |
| III. Cross-Platform Portability | **PASS** | `-std=c11` supported by GCC ≥ 4.7, Clang ≥ 3.1, MinGW-w64 (MSYS2 MINGW64). `_POSIX_C_SOURCE=200809L` recognised by glibc, musl, and MinGW-w64 headers. `clock_gettime` gated by `HAVE_CLOCK_GETTIME`; `time()`/`clock()` fallback retained for non-POSIX hosts. **Note**: Constitution currently lists "C99" as target language; this feature upgrades to C11, which is a strict superset of C99. All existing C99-conforming code is valid C11; no existing functionality is changed. |
| IV. Repair Before Extension | **PASS** | Pure repair/hardening effort. No new features, no new public API surface. |
| V. CLI-First, Library-Ready | **N/A** | No interface changes; public header `dw_solver.h` and `dw_options_t` struct are unmodified. |

**Post-design re-check**: No constitution violations identified after Phase 1 design. Audit
(`contracts/posix-header-audit.md`) confirms all POSIX functions used are in POSIX.1-2008;
audit of compiler extensions finds all non-standard constructs already correctly guarded.

## Project Structure

### Documentation (this feature)

```text
specs/013-strict-c-posix-compliance/
├── plan.md                       # This file
├── research.md                   # Phase 0 output — POSIX header audit, timing call-site inventory
├── data-model.md                 # Phase 1 output — build-flag change model
├── quickstart.md                 # Phase 1 output — "build and verify compliance" developer guide
├── contracts/
│   └── posix-header-audit.md     # Per-file POSIX header table (FR-001/FR-002/FR-005 evidence)
└── tasks.md                      # Phase 2 output (/speckit.tasks command — NOT created by this command)
```

### Source Code

```text
configure.ac             # + AC_SEARCH_LIBS([clock_gettime], [rt])
                         # + AC_CHECK_FUNC([clock_gettime]) → AC_DEFINE([HAVE_CLOCK_GETTIME])
                         # - AC_CHECK_HEADER([sys/time.h]) (vestigial)
                         # - AC_CHECK_FUNC([gettimeofday]) (vestigial)

src/Makefile.am          # + AM_CFLAGS = -std=c11 -pedantic-errors
                         # + AM_CPPFLAGS = -D_POSIX_C_SOURCE=200809L

src/dw_solver.c          # Replace time()/clock() with clock_gettime(CLOCK_MONOTONIC)
                         # under #ifdef HAVE_CLOCK_GETTIME; retain time()/clock() in #else

src/dw_support.c         # Replace time()/clock() in print_timing() under same guard

# All other src/dw_*.c and src/dw_*.h: audit-only; no source changes anticipated
```

**Structure Decision**: Single-project; all source changes in `src/` plus `configure.ac`
and `src/Makefile.am`. No new source directories. Audit documents in
`specs/013-strict-c-posix-compliance/contracts/`.

## Complexity Tracking

No constitution violations requiring justification. This is a hardening feature with direct
precedent in feature 008. C11 upgrade is backwards-compatible (C11 ⊃ C99).

---

## Implementation Phases

### Phase 1 — Build System Hardening (US1, US2: Priority P1)

**Goal**: Every `src/dw_*.c` compilation is gated by `-std=c11 -pedantic-errors
-D_POSIX_C_SOURCE=200809L`.

**Changes**:
1. `src/Makefile.am` — add at the top, before any target:
   ```makefile
   AM_CFLAGS   = -std=c11 -pedantic-errors
   AM_CPPFLAGS = -D_POSIX_C_SOURCE=200809L
   ```
2. `configure.ac` — add after the `AC_CHECK_FUNC([gettimeofday])` block:
   ```autoconf
   AC_SEARCH_LIBS([clock_gettime], [rt])
   AC_CHECK_FUNC([clock_gettime],
     AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [Define if clock_gettime is available.]))
   ```
3. `configure.ac` — remove the vestigial `AC_CHECK_HEADER([sys/time.h])` and
   `AC_CHECK_FUNC([gettimeofday])` blocks (neither macro is consumed by any source).
4. `config.h.in` — run `autoheader` to regenerate: removes `HAVE_GETTIMEOFDAY` /
   `HAVE_SYS_TIME_H` slots; adds `HAVE_CLOCK_GETTIME` slot.

**Verification**: `./configure && make` produces zero warnings or errors for any
`src/dw_*.c` file. Exit criterion: SC-001, SC-002.

---

### Phase 2 — Timing Migration (FR-009)

**Goal**: Replace `time()`/`clock()` calls in timing code with `clock_gettime(CLOCK_MONOTONIC)`
where POSIX is available; retain standard C fallback.

**Changes**:

**`src/dw_solver.c`** (~lines 116–120): Replace
```c
time_t  t0 = time(NULL);
clock_t c0 = clock();
```
with:
```c
#ifdef HAVE_CLOCK_GETTIME
struct timespec dw_ts0;
clock_gettime(CLOCK_MONOTONIC, &dw_ts0);
#else
time_t  t0 = time(NULL);
clock_t c0 = clock();
#endif
```
And at the point where timing is reported (passed to `print_timing` or used directly),
pass `dw_ts0` (or `t0`/`c0`) accordingly.

**`src/dw_support.c::print_timing()`**: Update signature and body to accept either a
`struct timespec` or `time_t`/`clock_t` pair depending on the guard. Alternatively, add
a second `print_timing_monotonic(struct timespec *ts0)` function and call the appropriate
one from `dw_solver.c` — whichever is cleaner given the actual call site.

**Note**: The spec says "preserve `HAVE_GETTIMEOFDAY` as `#else` fallback". Since no
source uses `gettimeofday`, the correct `#else` fallback is the existing `time()`/`clock()`
ISO C path — equally portable and requires no POSIX. See `research.md §Topic 3`.

**Verification**: Build with and without `HAVE_CLOCK_GETTIME` defined (simulate by
commenting out the define in `config.h`); both paths must compile and produce timing output.

---

### Phase 3 — Source Audit and SEI CERT Regression Check (US3, US4)

**Goal**: Confirm zero unguarded extensions; confirm SEI CERT C baseline is not regressed.

**Actions**:
1. Build with `-pedantic-errors` (already enforced by Phase 1). Zero errors = FR-004 PASS.
2. Run `cppcheck --enable=all --suppress=missingIncludeSystem -I src/ src/dw_*.c` and
   diff against `specs/008-sei-cert-c-compliance/baseline-cppcheck.txt`.
3. Run `clang --analyze -std=c11 -D_POSIX_C_SOURCE=200809L -I src/ src/dw_*.c` and
   diff against `specs/008-sei-cert-c-compliance/baseline-warnings.txt`.
4. Run `make check` and confirm all tests pass on macOS and Linux.

**Verification**: Zero new lines in cppcheck diff; zero new lines in clang diff; `make check`
exits 0. Exit criterion: SC-003, SC-004, FR-007, FR-008.

---

### Phase 4 — Compliance Acceptance Report (FR-006, SC-005)

**Goal**: Deliver `acceptance-report.md` following the feature 008 format.

**Actions**:
1. Create `specs/013-strict-c-posix-compliance/acceptance-report.md`.
2. One section per FR (FR-001 through FR-009), each containing:
   - Requirement restatement
   - Evidence bullets (specific files, line ranges, command output)
   - Explicit ✅ PASS or ❌ FAIL verdict
3. All verdicts must be ✅ PASS before this feature is considered done.

**Reference**: `specs/008-sei-cert-c-compliance/acceptance-report.md` format.
