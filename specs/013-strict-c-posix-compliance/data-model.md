# Data Model: Strict ISO C and POSIX.1 Compliance

This feature has no domain data model (no new structs, no database schema, no public API
changes). The "model" for this feature is the set of build-system and source-level changes
and the invariants they establish. This document captures that model.

---

## Build-Flag Change Model

### Before (baseline)

| Mechanism | Value | Where defined |
|-----------|-------|---------------|
| ISO C standard | *(none — compiler default)* | N/A |
| Pedantic mode | *(none)* | N/A |
| POSIX feature-test macro | *(none)* | N/A |
| Timing: wall clock | `time()` — 1-second resolution | `<time.h>` ISO C |
| Timing: CPU time | `clock()` — CLOCKS_PER_SEC ticks | `<time.h>` ISO C |
| `clock_gettime` | *(not checked — not used)* | N/A |
| `gettimeofday` | Checked in configure.ac; never used in source | vestigial |

### After (target)

| Mechanism | Value | Where defined |
|-----------|-------|---------------|
| ISO C standard | `-std=c11` | `AM_CFLAGS` in `src/Makefile.am` |
| Pedantic mode | `-pedantic-errors` | `AM_CFLAGS` in `src/Makefile.am` |
| POSIX feature-test macro | `-D_POSIX_C_SOURCE=200809L` | `AM_CPPFLAGS` in `src/Makefile.am` |
| Timing: wall clock (primary) | `clock_gettime(CLOCK_MONOTONIC)` | `#ifdef HAVE_CLOCK_GETTIME` |
| Timing: wall clock (fallback) | `time()` | `#else` |
| Timing: CPU time | `clock()` where fallback applies | `#else` |
| `clock_gettime` | Checked via `AC_CHECK_FUNC` + `AC_SEARCH_LIBS` | `configure.ac` |
| `gettimeofday` | Removed (vestigial check eliminated) | N/A |

---

## Entity: AM_CFLAGS

**Purpose**: Compiler flags applied to all compilation units in `src/`.
**Owner file**: `src/Makefile.am`
**Current value**: Empty (not set).
**Target value**: `-std=c11 -pedantic-errors`
**Propagation**: Inherited by `libdwblas_internal_la_CFLAGS`, `libdwsolver_la_CFLAGS`,
`dwsolver_CFLAGS` via `$(AM_CFLAGS)` references already present in `src/Makefile.am`.

## Entity: AM_CPPFLAGS

**Purpose**: Preprocessor flags applied to all compilation units.
**Owner file**: `src/Makefile.am`
**Current value**: Empty (not set).
**Target value**: `-D_POSIX_C_SOURCE=200809L`
**Propagation**: Inherited by `libdwsolver_la_CPPFLAGS`, `dwsolver_CPPFLAGS` via `$(AM_CPPFLAGS)`.

## Entity: HAVE_CLOCK_GETTIME

**Purpose**: Config-time detection of `clock_gettime(2)` availability.
**Owner file**: `configure.ac` (check), `config.h.in` (slot), generated `config.h` (value).
**Detection method**:
1. `AC_SEARCH_LIBS([clock_gettime], [rt])` — links `-lrt` on Linux if needed.
2. `AC_CHECK_FUNC([clock_gettime], AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [...]))`.
**Consumer files**: `src/dw_solver.c` (timing start), `src/dw_support.c::print_timing()`.

## Entity: HAVE_GETTIMEOFDAY (removed)

**Purpose**: Was a configure-time check for `gettimeofday(2)`.
**Action**: Remove `AC_CHECK_HEADER([sys/time.h])`, `AC_CHECK_FUNC([gettimeofday])`,
and the `HAVE_SYS_TIME_H` / `HAVE_GETTIMEOFDAY` entries from `configure.ac` and
`config.h.in`. No source file references these macros, so deletion is safe.

---

## Timing Data Flow (after change)

```text
configure.ac
  └── AC_SEARCH_LIBS([clock_gettime], [rt])  →  adds -lrt to LIBS on Linux
  └── AC_CHECK_FUNC([clock_gettime])
      └── AC_DEFINE([HAVE_CLOCK_GETTIME])    →  config.h: #define HAVE_CLOCK_GETTIME 1

src/dw_solver.c  (includes <config.h>)
  #ifdef HAVE_CLOCK_GETTIME
    struct timespec ts0, ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts0);
    ... solver run ...
    clock_gettime(CLOCK_MONOTONIC, &ts1);
  #else
    time_t t0 = time(NULL);
    clock_t c0 = clock();
  #endif

src/dw_support.c::print_timing()
  #ifdef HAVE_CLOCK_GETTIME
    struct timespec ts1;
    clock_gettime(CLOCK_MONOTONIC, &ts1);
    // compute elapsed nanoseconds from stored ts0
  #else
    clock_t c1 = clock(); time_t t1 = time(NULL);
  #endif
```

---

## Invariants Established by This Feature

1. Every `src/dw_*.c` and `src/dw_*.h` compilation is gated by `-std=c11 -pedantic-errors`.
2. Every translation unit sees `_POSIX_C_SOURCE=200809L` before any system header.
3. No unguarded GCC/Clang extension appears in authored source.
4. All timing code uses `CLOCK_MONOTONIC` where available, with `time()`/`clock()` fallback.
5. SEI CERT C rule conformance from feature 008 is not weakened.
