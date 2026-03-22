# Research: Callable Library Interface

**Phase**: 0 — Technology Decisions  
**Feature**: 010-callable-library  
**Date**: 2026-03-22

All NEEDS CLARIFICATION items from the spec and plan have been resolved below.

---

## Topic 1: Autotools Dual-Artifact Build (lib + CLI from same Makefile.am)

**Decision**: Use `lib_LTLIBRARIES` and `bin_PROGRAMS` in the same `src/Makefile.am`.
The binary links the library via `dwsolver_LDADD = libdwsolver.la`.

**Rationale**: This is the canonical autotools/libtool pattern for projects that produce
both a library and a tool using that library. It requires no additional autoconf macros
and no separate build directories. The existing `src/Makefile.am` already lists all
sources; the change is to split them: library sources get everything except `dw_main.c`,
and a new `dw_solver.c`; the binary gets `dw_main.c` and links `libdwsolver.la`.

**How it works**:
```makefile
lib_LTLIBRARIES = libdwsolver.la
libdwsolver_la_SOURCES = ... (all glp*.c and dw_*.c except dw_main.c) dw_solver.c
libdwsolver_la_LDFLAGS = -version-info 0:0:0

bin_PROGRAMS = dwsolver
dwsolver_SOURCES = dw_main.c
dwsolver_LDADD   = libdwsolver.la
```

When `make` runs, libtool compiles library objects with `-fPIC` (or equivalent) and links
them into both `libdwsolver.a` and `libdwsolver.so` / `libdwsolver.dylib`. The binary
compile step is a normal non-PIC compile of `dw_main.c` that links the libtool library.

**Alternatives considered**:
- _Separate top-level Makefile targets_: rejected; would duplicate source lists and
  diverge from the established single-directory build convention.
- _CMake migration_: rejected; out of scope, and autotools is already fully functional
  with libtool configured (`libtool` and `ltmain.sh` are present in the repo).

---

## Topic 2: Shared Library Symbol Visibility

**Decision**: Compile the library with `-fvisibility=hidden` and tag public API functions
with `__attribute__((visibility("default")))` wrapped in a `DWSOLVER_API` macro.

**Rationale**: The embedded GLPK contributes ~100 C files and hundreds of global symbols.
Without visibility control, all of those symbols are exported from the shared library,
which leaks implementation details, causes link conflicts if a host application also uses
a system GLPK, and bloats the dynamic symbol table. `-fvisibility=hidden` defaults all
symbols to hidden; only those explicitly marked `DWSOLVER_API` in `dw_solver.h` are
exported.

**DWSOLVER_API macro definition** (in `dw.h` or `dw_solver.h`):

```c
#if defined(_WIN32) || defined(__CYGWIN__)
#  ifdef DWSOLVER_BUILDING_LIB
#    define DWSOLVER_API __declspec(dllexport)
#  else
#    define DWSOLVER_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define DWSOLVER_API __attribute__((visibility("default")))
#else
#  define DWSOLVER_API
#endif
```

`-DDWSOLVER_BUILDING_LIB` is added to `libdwsolver_la_CFLAGS` so `DWSOLVER_API` expands
to `dllexport` during library compilation and `dllimport` when consumer headers are
included from downstream code.

**Alternatives considered**:
- _GNU export map / linker version script_: provides the same effect but is a separate
  file with its own maintenance burden; the visibility attribute approach is simpler and
  works across GCC, Clang, and MSVC.
- _No visibility control_: rejected; leaks GLPK internals, risks symbol collisions.

---

## Topic 3: pkg-config Integration

**Decision**: Add `dwsolver.pc.in` at the repository root; wire into `configure.ac` via
`AC_CONFIG_FILES([dwsolver.pc])`. Install to `$(libdir)/pkgconfig` via `pkgconfigdir`
and `pkgconfig_DATA` in `Makefile.am`.

**Template (`dwsolver.pc.in`)**:

```
prefix=@prefix@
exec_prefix=@exec_prefix@
libdir=@libdir@
includedir=@includedir@

Name: dwsolver
Description: Dantzig-Wolfe Decomposition Solver Library
Version: @PACKAGE_VERSION@
Libs: -L${libdir} -ldwsolver
Cflags: -I${includedir}
```

**configure.ac addition**:
```
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile dwsolver.pc])
```

**top-level Makefile.am addition**:
```
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = dwsolver.pc
```

**Rationale**: `pkg-config` is the standard mechanism on Linux/macOS for consumers to
discover compiler and linker flags. No additional autoconf macros (`PKG_CHECK_MODULES`,
etc.) are needed for the project itself — those are for consumers. The template uses
only autoconf substitution variables that are already defined.

**Alternatives considered**:
- _cmake-config files_: not a fit for an autotools-based project without a CMake migration.
- _Manual documentation of `-I` and `-L` flags_: insufficient; prevents integration with
  build systems like meson, cmake (via `pkg_check_modules`), and cargo-sys crates.

---

## Topic 4: Public API Struct Isolation vs. Exposing faux_globals

**Decision**: Define new `dw_options_t` and `dw_result_t` structs in `dw_solver.h`.
Keep `faux_globals` strictly internal (not in the public header). The library entry
point converts from `dw_options_t` to `faux_globals` as the first step.

**Rationale**: `faux_globals` contains internal pointers (service queues, thread state,
GLPK problem objects) that are implementation details. Exposing it publicly would lock
the ABI to the internal representation. A narrow `dw_options_t` with only the 11
tunable parameters that callers need is a stable, minimal surface.

**dw_options_t fields** (derived from `faux_globals` + `process_cmdline` analysis):

| Field | Type | Default | Corresponds to |
|-------|------|---------|----------------|
| `verbosity` | `int` | `OUTPUT_NORMAL` (5) | `fg->verbosity` |
| `mip_gap` | `double` | `DEFAULT_MIP_GAP` (0.01) | `fg->mip_gap` |
| `max_phase1_iterations` | `int` | `MAX_PHASE1_ITERATIONS` (100) | `fg->max_phase1_iterations` |
| `max_phase2_iterations` | `int` | `MAX_PHASE2_ITERATIONS` (3000) | `fg->max_phase2_iterations` |
| `rounding_flag` | `int` | 0 | `fg->rounding_flag` |
| `integerize_flag` | `int` | 0 | `fg->integerize_flag` |
| `enforce_sub_integrality` | `int` | 0 | `fg->enforce_sub_integrality` |
| `print_timing_data` | `int` | 0 | `fg->print_timing_data` |
| `print_final_master` | `int` | 0 | `fg->print_final_master` |
| `print_relaxed_sol` | `int` | 0 | `fg->print_relaxed_sol` |
| `perturb` | `int` | 0 | `fg->perturb` |
| `shift` | `double` | 0.0 | `fg->shift` |

**dw_result_t fields**:

| Field | Type | Description |
|-------|------|-------------|
| `status` | `int` | 0 = success; non-zero = error code |
| `objective_value` | `double` | Final LP relaxation objective |
| `num_vars` | `int` | Length of `x` array (number of master LP columns) |
| `x` | `double*` | Primal solution vector (library-allocated; freed by `dw_result_free`) |

**Alternatives considered**:
- _Expose `faux_globals` directly_: rejected; ABI instability, too wide a surface.
- _Opaque handle (void*) pattern_: considered; rejected as unnecessarily complex for a
  single-call API with explicit result return.

---

## Topic 5: Reentrancy and Concurrent Use

**Decision**: Document in `dw_solver.h` (API header) that calling `dw_solve()` concurrently
from multiple threads is **not supported**. The library uses process-wide POSIX global
synchronization objects (`pthread_mutex_t` and `sem_t` globals in `dw_globals.c`) that
are initialized once and reused across phases. Two concurrent calls would race on these objects.

**Rationale**: Fixing reentrancy is out of scope for this feature. Moving the globals
into a per-call context object would require restructuring a large portion of the
threading model — itself a substantial feature. The limitation will be clearly stated
in the header so callers know to serialize calls externally if needed.

**Documentation text** (to appear in `dw_solver.h`):
```
 * THREAD SAFETY: This library is NOT reentrant. Only one call to dw_solve()
 * may be in progress at any time within a process. The library uses process-wide
 * synchronization objects that are not safe for concurrent re-entry. Callers
 * that need to run multiple solves from different threads must serialize calls
 * with an external mutex.
```

---

## Topic 6: test_lib_api.c — Library API Test Program

**Decision**: Add `tests/test_lib_api.c`, a C program that:
1. Calls `dw_options_init()` and verifies defaults match `dw.h` constants.
2. Calls `dw_solve()` on the `book_bertsimas` example and checks the objective value.
3. Calls `dw_result_free()` and verifies no memory is leaked (detectable with valgrind/ASan).
4. Calls `dw_solve()` with a NULL master file and verifies a non-zero error code is returned.

The test is added to `tests/Makefile.am` as a `check_PROGRAMS` target so it runs under
`make check`. It links against `libdwsolver` via the local build tree.

**Rationale**: The existing `dw-tests.sh` shell script only validates the CLI binary.
A dedicated C test program is needed to validate the library API independently.

---

## Summary of Resolved Items

| Was | Now |
|-----|-----|
| `faux_globals` as the only configuration interface | `dw_options_t` public struct + internal `faux_globals` conversion |
| All symbols exported from .so | Only `DWSOLVER_API`-tagged symbols exported; GLPK internals hidden |
| No pkg-config support | `dwsolver.pc` generated and installed |
| CLI-only test coverage | `test_lib_api.c` added for library API test coverage |
| Concurrent call behavior undefined | Documented as unsupported in public header |
