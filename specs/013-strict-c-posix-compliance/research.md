# Research: Strict ISO C and POSIX.1 Compliance — Codebase Audit Findings

## Scope

Audit of all authored files in `src/dw_*.c` and `src/dw_*.h` against the requirements of
FR-001 through FR-009. Vendored files (`third-party/glpk/`) and Intel MKL conditional blocks
(`#ifdef USE_INTEL_MKL`) are excluded.

---

## Topic 1: Current Build System Baseline (FR-001, FR-002, FR-003)

### Decision: What flags are currently active?

- **`-std=` flag**: None. The build system sets no ISO C standard.
  Source: `src/Makefile.am` — `AM_CFLAGS` is never set; per-target `*_CFLAGS` reference
  `$(AM_CFLAGS)` which expands to empty.
- **`_POSIX_C_SOURCE`**: Not defined anywhere. Neither `configure.ac` nor `Makefile.am`
  defines this macro. The macro is not in `config.h.in`.
- **`-pedantic` / `-pedantic-errors`**: Not present in `configure.ac` or any `Makefile.am`.
- **Implied behavior**: Compiler defaults to its own GNU extension mode (GCC: `-std=gnu11`
  or `gnu17`; Clang: `-std=gnu17`). POSIX extensions become available implicitly via
  `_GNU_SOURCE` on glibc, not via a declared interface.

### Decision: Target standard

**C11** (`-std=c11`) via `AM_CFLAGS` in `src/Makefile.am`.
**Rationale**: Per clarification Q1; C11 is a superset of C99 (all existing code is valid C11).
C11 introduced `<stdnoreturn.h>`, `_Static_assert`, and `<threads.h>` — none used by this
codebase, so no new APIs are adopted, but the door is open for `_Static_assert` as a
safe-guard pattern if needed in future.
**Alternatives considered**: C99 (too conservative; C11 has been stable for 15 years);
C17 (identical to C11 except for defect-report integrations; no practical difference here).

### Decision: Pedantic treatment

**`-pedantic-errors`** alongside `-std=c11`.
**Rationale**: Per clarification Q4. Hard errors prevent silent accumulation of extension
usage. Matches the philosophy of feature 008 (zero new violations, not "warnings are OK").

### Decision: POSIX feature-test macro placement

**`AM_CPPFLAGS` in `src/Makefile.am`** with `-D_POSIX_C_SOURCE=200809L`.
**Rationale**: preprocessor defines belong in `CPPFLAGS`, not `CFLAGS`. `AM_CPPFLAGS` is
already referenced by `libdwsolver_la_CPPFLAGS` and `dwsolver_CPPFLAGS` via `$(AM_CPPFLAGS)`,
so a single addition propagates to all targets.
**Alternatives considered**: `AC_DEFINE` in `configure.ac` — valid but less conventional for
feature-test macros, which are typically set before any system header is processed;
AM_CPPFLAGS achieves this without regenerating config.h.in.

---

## Topic 2: POSIX Header Inventory (FR-002, FR-005)

Per-file inventory of all POSIX and implementation-extension headers. "POSIX.1-2008" means
the header is defined in POSIX.1-2008 and its content is gated by `_POSIX_C_SOURCE=200809L`.

| Source File | Direct POSIX Headers | Gets POSIX Via |
|-------------|---------------------|----------------|
| `dw.h` | `<semaphore.h>`, `<fcntl.h>`, `<pthread.h>` | (hub) |
| `dw_main.c` | *(none directly)* | `dw_solver.h` → `dw.h` |
| `dw_globals.c` | *(none directly)* | `dw.h` |
| `dw_phases.c` | `<pthread.h>` | `dw.h` |
| `dw_support.c` | `<unistd.h>`, `<pthread.h>` | `dw.h` |
| `dw_solver.c` | *(none directly)* | `dw.h` |
| `dw_subprob.c` | `<pthread.h>` | `dw.h` |
| `dw_rounding.c` | `<pthread.h>` | `dw.h` |
| `dw_blas.c` | *(none)* | *(none)* |

### POSIX.1-2008 verification for each header

| Header | POSIX.1-2008 Status | Notes |
|--------|---------------------|-------|
| `<semaphore.h>` | ✅ POSIX.1-2008 §2.9 | Named + unnamed semaphores; `_POSIX_C_SOURCE=200809L` exposes both |
| `<fcntl.h>` | ✅ POSIX.1-2008 §6.5 | Used for `O_CREAT`, `O_RDWR` in named semaphore path |
| `<pthread.h>` | ✅ POSIX.1-2008 §17 | pthreads; fully defined under `_POSIX_C_SOURCE=200809L` |
| `<unistd.h>` | ✅ POSIX.1-2008 §2.14 | `sleep()`, `getpid()`; used in `dw_support.c` |
| `<sys/time.h>` | ✅ POSIX.1-2008 (XSI) | `gettimeofday()` — present in `configure.ac` check but **NOT used in any source file** (see Topic 3) |

### include-before-guard pattern in dw.h

`dw.h` places `#include <semaphore.h>`, `#include <fcntl.h>`, and `#include <pthread.h>`
**before** the `#ifndef DECOMPOSE_H_` include guard (lines 35–42 vs. guard at line 43).
This is non-standard practice but not a compliance violation per se — the includes are
idempotent by their own guards. It SHOULD be noted in the acceptance report.
**Decision**: Leave the include order as-is; moving them risks ordering side-effects with
the visibility macro block that follows. The `_POSIX_C_SOURCE` macro (coming from `AM_CPPFLAGS`)
is processed before any `#include` anyway, so its placement is correct.

---

## Topic 3: Timing Code Audit — FR-009 Scope Correction

### Finding: `gettimeofday` is NOT called in any authored source file

The spec's FR-009 references a "`gettimeofday` call in `dw_solver.c`". Investigation reveals:

- `configure.ac` lines 56–58: `AC_CHECK_FUNC([gettimeofday], AC_DEFINE([HAVE_GETTIMEOFDAY], [1], [N/A]))` — this check exists.
- `config.h.in` line 6: `#undef HAVE_GETTIMEOFDAY` — the macro slot exists.
- **No source file** uses `#ifdef HAVE_GETTIMEOFDAY`, `gettimeofday()`, or includes `<sys/time.h>` directly.

The configure check is a **vestigial artifact** from an earlier version of the codebase.

### Actual timing call sites (ISO C only)

| File | Lines | Code | Standard |
|------|-------|------|----------|
| `dw_solver.c` | ~118–119 | `time_t t0 = time(NULL); clock_t c0 = clock();` | ISO C |
| `dw_support.c` | ~828–835 | `print_timing()`: `clock()`, `time(NULL)` | ISO C |

Both call sites use `time()` (`<time.h>`, ISO C) and `clock()` (`<time.h>`, ISO C).

### Decision: FR-009 implementation scope

FR-009 as written calls for migrating from `gettimeofday` to `clock_gettime`. Since
`gettimeofday` is not actually used, the correct implementation is:

1. **Remove** the vestigial `AC_CHECK_FUNC([gettimeofday])` from `configure.ac`.
2. **Remove** the `HAVE_GETTIMEOFDAY` and `HAVE_SYS_TIME_H` checks from `configure.ac`
   (neither is used in any source).
3. **Add** `AC_SEARCH_LIBS([clock_gettime], [rt])` (handles Linux `-lrt` requirement).
4. **Add** `AC_CHECK_FUNC([clock_gettime])` + `AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [Define if clock_gettime is available.])`.
5. **Update `dw_solver.c`**: under `#ifdef HAVE_CLOCK_GETTIME`, use
   `clock_gettime(CLOCK_MONOTONIC, &ts)` for wall-clock timing; `#else` retains `time(NULL)`.
6. **Update `dw_support.c::print_timing()`**: same guard pattern; `#else` retains `clock()`
   and `time(NULL)` for non-POSIX hosts.
7. Update `config.h.in` to add `HAVE_CLOCK_GETTIME` slot; run `autoheader`.

**Rationale**: `clock_gettime(CLOCK_MONOTONIC)` is monotonic (immune to wall-clock
adjustments), has nanosecond precision, and is defined in POSIX.1-2008 §14.2. This aligns
with CERT MSC15-C (use monotonic clocks for elapsed time). `clock()` measures CPU time
and can wrap on 32-bit systems after ~72 minutes; it should also be replaced where
sub-second precision is desired, but it can remain in the `#else` fallback for portability.

**Note**: On macOS, `clock_gettime(CLOCK_MONOTONIC)` is available since macOS 10.12; the
existing CI runs on a recent macOS version so this is not a risk.

---

## Topic 4: Compiler Extension Audit — FR-004

### Audit findings across src/dw_*.c and src/dw_*.h

| Construct | File(s) | Status |
|-----------|---------|--------|
| `__attribute__((visibility("default")))` | `dw.h` | ✅ Correctly guarded by `#ifdef __GNUC__` |
| `__declspec(dllexport)` / `__declspec(dllimport)` | `dw.h` | ✅ Correctly guarded by `#if defined(_WIN32) \|\| defined(__CYGWIN__)` |
| Zero-length arrays | None found | ✅ Clean |
| `__typeof__` | None found | ✅ Clean |
| Statement expressions `({ ... })` | None found | ✅ Clean |
| `__builtin_*` intrinsics | None found | ✅ Clean |
| Variable-length arrays (VLAs) | None found in authored files | ✅ Clean |
| `__attribute__((noreturn))` | None found unguarded | ✅ Clean |

**Decision**: No source changes required for FR-004. The `-pedantic-errors` flag will
enforce this contract in CI going forward.

---

## Topic 5: SEI CERT C Baseline Reference — FR-007

The feature 008 static analysis baselines to use as regression floor:

- `specs/008-sei-cert-c-compliance/baseline-cppcheck.txt` — 8 cppcheck warnings total
  (7 in `dw_support.c`, 1 in `dw_main.c`; all are dead-store / unix.BlockInCriticalSection)
- `specs/008-sei-cert-c-compliance/baseline-warnings.txt` — compiler warning baseline

**Decision**: After applying FR-001–FR-009, re-run `cppcheck --enable=all src/dw_*.c` and
`clang --analyze -I src/ src/dw_*.c` and diff against these baselines. The `clock_gettime`
migration may eliminate the dead-store warnings on `t0`/`c0` in `dw_main.c` if those
variables are now read (improvement, not regression).

---

## Resolution Summary

| Unknown | Resolution |
|---------|------------|
| C standard target | **C11** — clarification Q1 |
| POSIX macro value | **200809L** — POSIX.1-2008 |
| Pedantic level | **`-pedantic-errors`** — clarification Q4 |
| `gettimeofday` call sites | **0 call sites** — configure check is vestigial; FR-009 means adding `clock_gettime` as the new primary path with `time()`/`clock()` fallback |
| `clock_gettime` availability on macOS | ✅ macOS 10.12+; CI runs on recent macOS versions |
| Windows CI | ✅ `ci-windows.yml` (MSYS2 MINGW64) already exists; `continue-on-error: true` |
| Compiler extensions in authored source | **None unguarded** — FR-004 is an audit-only item |
| `AM_CFLAGS` vs per-target `CFLAGS` | Use `AM_CFLAGS` to propagate to all targets uniformly |
