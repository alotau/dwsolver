# Implementation Plan: Callable Library Interface

**Branch**: `010-callable-library` | **Date**: 2026-03-22 | **Spec**: [spec.md](spec.md)  
**Input**: Feature specification from `/specs/010-callable-library/spec.md`

## Summary

Extract the Dantzig-Wolfe solver logic from `dw_main.c` into a dedicated `dw_solver.c`
implementation, expose a minimal public C API via a new `dw_solver.h` header installed
as `dwsolver.h`, and update `src/Makefile.am` to produce `libdwsolver` (static + shared
via libtool) alongside the existing `dwsolver` binary. The CLI becomes a thin wrapper
that calls the public API. A `dwsolver.pc.in` pkg-config template is added and wired
through `configure.ac`.

## Technical Context

**Language/Version**: C99  
**Primary Dependencies**: GLPK 4.44 (thread-patched, embedded in `src/`); POSIX Threads (pthreads); GNU Autotools + libtool  
**Storage**: File-based (LP problem files on disk; no database)  
**Testing**: Shell-based regression suite (`tests/dw-tests.sh`); new C test program linking against the library  
**Target Platform**: Linux, macOS (Windows cross-platform support inherited from existing build)  
**Project Type**: CLI tool + callable C library (dual-artifact build)  
**Performance Goals**: No new performance requirements; solve runtime identical to current CLI  
**Constraints**: All existing `dw-tests.sh` tests must continue to pass; no new solver globals; no external GLPK dependency exposed to consumers  
**Scale/Scope**: 6 core `dw_*.c` source files; ~900-line `dw_main.c` split into solver core (~800 lines) and thin CLI wrapper (~100 lines)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | ✅ PASS | Regression suite (`tests/dw-tests.sh`) is a mandatory acceptance gate. The refactor moves solver logic verbatim — no algorithmic changes. |
| II. Thread Safety is Mandatory | ✅ PASS | No changes to threading model, synchronization primitives, or shared global state. Concurrent API calls with shared context are documented as unsupported rather than introduced or broken. |
| III. Cross-Platform Portability | ✅ PASS | libtool (already present in the autotools setup) handles `.so` / `.dylib` / `.dll` differences transparently. No new platform-specific code is added. |
| IV. Repair Before Extension | ✅ PASS | Repair baseline is complete (specs 001–009). This is confirmed extension work on a clean, passing build. |
| V. CLI-First, Library-Ready | ✅ PASS | This feature directly delivers the "Library-Ready" outcome stated in Principle V. CLI behavior is fully preserved. |

**Gate result**: All five principles satisfied. Proceed to Phase 0.

## Project Structure

### Documentation (this feature)

```text
specs/010-callable-library/
├── plan.md              # This file
├── research.md          # Phase 0: technology decisions
├── data-model.md        # Phase 1: public API types and data structures
├── contracts/
│   └── api.md           # Phase 1: C API contract (functions, parameters, return values)
├── quickstart.md        # Phase 1: consumer integration guide
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code (repository root)

```text
src/
├── dw_solver.c          # NEW: extracted solver core (~800 lines from dw_main.c)
├── dw_solver.h          # NEW: public API header (installed as dwsolver.h)
├── dw_main.c            # MODIFIED: thin CLI wrapper (<100 lines) calling dw_solver API
├── dw.h                 # MODIFIED: add DWSOLVER_API visibility macro
├── dw_support.c / dw_support.h  # MODIFIED: process_cmdline signature updated to accept dw_options_t* + file-path out-parameters
├── dw_blas.c / dw_blas.h
├── dw_globals.c
├── dw_phases.c / dw_phases.h
├── dw_rounding.c / dw_rounding.h
├── dw_subprob.c / dw_subprob.h
├── dw_support.c / dw_support.h
├── glp*.c / glp*.h      # unchanged embedded GLPK
└── Makefile.am          # MODIFIED: add lib_LTLIBRARIES, nobase_include_HEADERS

dwsolver.pc.in           # NEW: pkg-config template (top-level, installed to $(libdir)/pkgconfig)
configure.ac             # MODIFIED: AC_CONFIG_FILES([dwsolver.pc])

tests/
├── dw-tests.sh          # UNCHANGED
└── test_lib_api.c       # NEW: C API regression test (links -ldwsolver, exercises all examples)
```

**Structure Decision**: Single-project layout. The library and CLI share the same `src/`
directory, which is the established convention. New files are additive; no subdirectory
reorganization is required.

## Phases

### Phase 0: Research — Technology Decisions

*Full findings in [research.md](research.md).*

| Topic | Decision | Rationale |
|-------|----------|-----------|
| Dual-artifact build (lib + CLI) | `lib_LTLIBRARIES` + `bin_PROGRAMS` in same `src/Makefile.am`; binary uses `dwsolver_LDADD = libdwsolver.la` | Standard autotools pattern; avoids separate build directories or CMake migration |
| Shared library symbol visibility | `-fvisibility=hidden` on library, `DWSOLVER_API` macro (`__attribute__((visibility("default")))`) on public functions | Hides all GLPK internals and private solver symbols; portable to MSVC via `__declspec(dllexport)` guard |
| pkg-config integration | `dwsolver.pc.in` template substituted by autoconf via `AC_CONFIG_FILES`; no extra m4 macros needed | Standard pattern; zero new autoconf dependencies |
| Public API struct isolation | New `dw_options_t` and `dw_result_t` structs in `dw_solver.h`; internal `faux_globals` stays private | Keeps ABI stable despite future internal refactoring; exposes only what callers need |
| Result memory ownership | Caller allocates `dw_result_t` on the stack; library fills sub-pointer buffers; `dw_result_free()` releases them | Explicit ownership; no opaque handles needed; compatible with stack allocation |
| Concurrency model | Not changed; document in `dw_solver.h` that concurrent calls with shared state are unsupported | Avoids scope creep; global state (GLPK, pthreads primitives) is inherently single-context |

### Phase 1: Design & Contracts

*Full detail in [data-model.md](data-model.md), [contracts/api.md](contracts/api.md),
and [quickstart.md](quickstart.md).*

#### Source Refactor Strategy

`dw_main.c` (891 lines) currently contains one large `main()` function that does:

1. **CLI setup** — `process_cmdline()` + `init_globals()`: ~35 lines → stays in slim
   `dw_main.c` wrapper
2. **Solver core** — everything from `init_signals()` through `free_globals()`: ~800 lines
   → moves to `dw_solver_run()` in `dw_solver.c`
3. **`hook()` static** — GLPK terminal redirect: moves to private static in `dw_solver.c`

The refactored `dw_main.c` becomes a ~70-line file:

```c
int main(int argc, char* argv[]) {
    dw_options_t  opts;
    dw_result_t   result = {0};
    const char   *master_file = NULL;
    const char  **subproblem_files = NULL;
    int           num_subproblems = 0;
    dw_options_init(&opts);
    int rc = process_cmdline(argc, argv, &opts,
                             &master_file,
                             &subproblem_files,
                             &num_subproblems);
    if (rc != 0) return rc;
    rc = dw_solve(master_file,
                  (const char *const *)subproblem_files,
                  num_subproblems,
                  &opts, &result);
    dw_result_free(&result);
    return rc;
}
```

#### Library Build Rules

New `src/Makefile.am` structure (abbreviated):

```makefile
lib_LTLIBRARIES = libdwsolver.la
libdwsolver_la_SOURCES  = $(GLPK_SOURCES) $(DW_CORE_SOURCES) dw_solver.c
libdwsolver_la_CFLAGS   = $(AM_CFLAGS) -fvisibility=hidden
libdwsolver_la_LDFLAGS  = -version-info 0:0:0

bin_PROGRAMS    = dwsolver
dwsolver_SOURCES = dw_main.c
dwsolver_CFLAGS  = $(AM_CFLAGS)
dwsolver_LDADD   = libdwsolver.la

nobase_include_HEADERS = dw_solver.h
```

`DW_CORE_SOURCES` = `dw_blas.c dw_support.c dw_phases.c dw_subprob.c dw_rounding.c dw_globals.c`

#### Library Versioning

Version-info `0:0:0` (current:revision:age) for the initial release. A `dw_version()`
API function returns the version string so consumers can detect it at runtime.

## Post-Phase-1 Constitution Re-check

| Principle | Re-check | Notes |
|-----------|----------|-------|
| I. Correctness First | ✅ | Solver core is moved verbatim; no algorithmic changes. Regression gate mandatory as final acceptance step. |
| II. Thread Safety | ✅ | Global state unchanged. `dw_solver.h` documents concurrent-context limitation explicitly. |
| III. Cross-Platform | ✅ | `DWSOLVER_API` macro degrades gracefully on non-GCC compilers; libtool handles shared-lib suffix differences. |
| IV. Repair Before Extension | ✅ | All existing tests still pass; no existing logic deleted or restructured beyond the extract. |
| V. CLI-First | ✅ | CLI behavior is bit-for-bit identical; it simply delegates through the library. |
