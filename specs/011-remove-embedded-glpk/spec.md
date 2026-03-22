# Feature Specification: Remove Embedded GLPK Source

**Feature Branch**: `011-remove-embedded-glpk`  
**Created**: 2026-03-22  
**Status**: Draft

## Background

The dwsolver project currently carries a full verbatim copy of the GLPK 4.44
source tree (~100 `.c`/`.h` files plus AMD and COLAMD helper libraries) inside
the `src/` directory.  This was necessary in 2010 because GLPK 4.44 stored its
internal environment pointer in a plain C global variable instead of genuine
Thread Local Storage (TLS), breaking multi-threaded use.  A single-file patch
(`third-party/glpk/glpk-4.44.ThreadReady.patch`) replaced the memory-management
routines with raw `malloc`/`free` to work around the race.

Modern GLPK (4.45 onward, currently 5.0) stores the environment pointer in real
per-thread TLS (`__thread` or POSIX `pthread_key_t`), making the workaround
unnecessary.  The project can therefore drop the embedded copy and link against
a system-installed GLPK library instead.

## Clarifications

### Session 2026-03-22

- Q: Should updating the CI workflow and Dockerfile to install system GLPK be in scope for this feature? → A: Yes — in scope; both must be updated in the same change
- Q: How should the ABI compatibility constraint (caller must link the same GLPK build as libdwsolver) be handled? → A: Document-only — state in README and pkg-config; no static linking required
- Q: What numeric tolerance should be used when comparing refactored-binary objective values against the baseline? → A: 1e-6 relative difference
- Q: Should dwsolver's configure propagate GLPK's transitive link dependencies (e.g., GMP) automatically or leave it to the user? → A: Automatic — use `pkg-config --libs glpk` so transitive deps are discovered without user intervention

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Build against system GLPK (Priority: P1)

A developer clones the repository, installs GLPK from their OS package manager
or from source, runs `./configure && make`, and gets a working `dwsolver`
binary and `libdwsolver` shared library — exactly as the existing build does —
without any GLPK source code inside the repository.

**Why this priority**: This is the core goal of the feature.  All other stories
depend on it.

**Independent Test**: A fresh checkout with no files in `src/glp*` or
`src/amd/` or `src/colamd/` must configure, compile, and link cleanly when
GLPK ≥ 4.65 is present on the system.

**Acceptance Scenarios**:

1. **Given** GLPK ≥ 4.65 is installed, **When** the developer runs `./configure && make`, **Then** the build succeeds with zero errors and produces `src/dwsolver` and `src/libdwsolver.la`.
2. **Given** GLPK is absent, **When** the developer runs `./configure`, **Then** configure exits with a clear error message identifying the missing dependency.
3. **Given** GLPK < 4.45 is installed, **When** the developer runs `./configure`, **Then** configure rejects it with a message stating the minimum required version.

---

### User Story 2 - Existing test suite passes (Priority: P1)

A developer runs the existing test suite after the refactoring and every
previously-passing test continues to pass.  The numerical results produced by
`dwsolver` are unchanged.

**Why this priority**: Equal priority with the build story — a build that
produces wrong answers is not a successful refactoring.

**Independent Test**: Run `tests/dw-tests.sh` with the refactored binary;
compare pass/fail counts against the pre-refactoring baseline.

**Acceptance Scenarios**:

1. **Given** the refactored binary, **When** `tests/dw-tests.sh` is run, **Then** all tests that passed before the change still pass.
2. **Given** the refactored binary is run on each example in `examples/`, **When** producing output, **Then** the optimal objective value matches the value produced by the old embedded-GLPK binary to within a relative difference of 1e-6.

---

### User Story 3 - Deprecated `lpx_*` API replaced (Priority: P2)

Three legacy `lpx_*` calls used in the dwsolver source (`lpx_read_cpxlp`,
`lpx_write_cpxlp`, `lpx_get_int_parm`) are replaced with their modern `glp_*`
equivalents so the code compiles cleanly against GLPK ≥ 4.65 without any
compatibility shims.

**Why this priority**: The deprecated API was removed from GLPK after 4.44.
Without this the build fails against any modern GLPK package.

**Independent Test**: The project compiles with zero errors referencing `lpx_`
symbols when built against GLPK 4.65 or later.

**Acceptance Scenarios**:

1. **Given** GLPK ≥ 4.65 installed (no `lpx_*` headers present), **When** the sources are compiled, **Then** no undefined-symbol or implicit-declaration errors for any `lpx_*` function appear.
2. **Given** a test problem read via the new API, **When** solved, **Then** the LP is parsed and solved identically to the old code path.

---

### User Story 4 - Repository source tree is clean (Priority: P3)

After the change, `src/` contains only dwsolver's own source files.  No GLPK
`.c` files, no AMD/COLAMD helpers, and no internal GLPK headers remain.
`git ls-files src/` lists only `dw_*.c`, `dw_*.h`, and `dw.h`.

**Why this priority**: Clean-up goal — verifiable by inspection and `git ls-files`.

**Independent Test**: `git ls-files src/ | grep -E '^src/glp|src/amd|src/colamd'` returns no output.

**Acceptance Scenarios**:

1. **Given** the refactoring is complete and committed, **When** `git ls-files src/` is run, **Then** no `glp*.c`, `glp*.h`, `amd/`, or `colamd/` entries appear.
2. **Given** the historical patch is preserved in `third-party/glpk/`, **When** the repository is browsed, **Then** `glpk-4.44.ThreadReady.patch` and the GLPK attribution files remain for provenance documentation.

---

### Edge Cases

- What happens if the system GLPK was compiled without thread support?  The configure check must detect and reject such installations if they lack TLS-enabled memory management (i.e., GLPK < 4.45).
- What happens if GLPK is installed in a non-standard prefix (e.g., `/opt/homebrew`)?  The configure step must accept `--with-glpk-prefix` or equivalent to locate headers and the library.
- What happens on platforms where GLPK is not available from a package manager?  The `README` must describe how to build GLPK from source and point at the minimum required version.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The build system MUST detect the presence of a system GLPK installation using `pkg-config` and fail with a user-readable error if it is absent or too old.
- **FR-002**: The build system MUST require GLPK version 4.65 or later (the first long-term supported release that guarantees TLS-based environment storage on all target platforms).
- **FR-002a**: The build system MUST obtain GLPK's compiler flags and full linker flags (including transitive dependencies such as GMP and zlib) via `pkg-config --cflags glpk` and `pkg-config --libs glpk`, so that builds against GMP-enabled or zlib-enabled GLPK packages succeed without requiring manual `--with-gmp` or `--with-zlib` flags.
- **FR-003**: The `src/` directory MUST NOT contain any GLPK source files (`glp*.c`, `glp*.h`) or their bundled dependencies (`amd/`, `colamd/`) after the change is applied.
- **FR-004**: All three deprecated `lpx_*` call sites MUST be replaced with their modern `glp_*` equivalents: `lpx_read_cpxlp` → `glp_create_prob` + `glp_read_lp`; `lpx_write_cpxlp` → `glp_write_lp`; `lpx_get_int_parm(lp, LPX_K_ITCNT)` → `glp_get_it_cnt(lp)`.
- **FR-005**: The existing `glpk_mutex` guards around GLPK file-I/O calls MUST be retained unchanged, as they protect sequential file-parsing state that remains non-reentrant even in modern GLPK.
- **FR-006**: The `third-party/glpk/` directory and all files within it MUST be retained in the repository for provenance and attribution purposes.
- **FR-007**: The project's `README.md` MUST document GLPK as an external dependency, including the minimum version, and provide instructions for installing it on all supported platforms.
- **FR-008**: The `dwsolver.pc` pkg-config file MUST list GLPK as a `Requires` dependency so that consumers of `libdwsolver` know they must also link against GLPK.
- **FR-008a**: The `README.md` MUST include a note that callers of `libdwsolver` must link against the same GLPK shared library that dwsolver was built against; mixing GLPK builds produces undefined behaviour. No static-linking enforcement is required.
- **FR-009**: The GitHub Actions CI workflow(s) MUST be updated to install GLPK ≥ 4.65 as a build-dependency step (e.g., `apt-get install libglpk-dev` on Ubuntu runners) before invoking `./configure && make`.
- **FR-010**: The project `Dockerfile` MUST be updated to install GLPK ≥ 4.65 inside the container image (e.g., via the distribution package manager) so that container builds continue to succeed after the embedded sources are removed.

### Key Entities

- **Embedded GLPK sources**: The ~100 `glp*.c`/`.h` files and `amd/`/`colamd/` subdirectories currently in `src/` — to be removed.
- **System GLPK library**: The external `libglpk` installation that the project will link against after the change.
- **Deprecated `lpx_*` API**: Three call sites in `dw_solver.c`, `dw_subprob.c`, and `dw_phases.c` that use removed GLPK 4.44 functions — to be migrated.
- **`glpk_mutex`**: The existing POSIX mutex that serialises GLPK file-I/O across threads — to be preserved.
- **`third-party/glpk/`**: Attribution and patch files — to be preserved.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The repository `src/` directory shrinks by at least 95 source files (all `glp*.c`, `glp*.h`, `amd/`, `colamd/` gone); `git ls-files src/ | wc -l` decreases by at least 95.
- **SC-002**: The project builds successfully against an unmodified system GLPK ≥ 4.65 package with zero linker errors on macOS and Linux.
- **SC-003**: All tests that pass in the pre-refactoring baseline continue to pass after the refactoring, as measured by `tests/dw-tests.sh` pass/fail counts; for any test that compares objective values numerically, the refactored binary must agree with the baseline to within a relative difference of 1e-6.
- **SC-004**: Zero `lpx_*` symbols remain in any dwsolver source file, verifiable by `grep -r 'lpx_' src/dw_*.c` returning no output.
- **SC-005**: A developer who has never built the project can follow the updated `README.md` to install GLPK, build dwsolver, and run the tests successfully without needing any guidance beyond the documented steps.
- **SC-006**: The GitHub Actions CI pipeline passes (green) after the change is merged, confirming that the automated runner can install GLPK and build and test the project without the embedded sources.
- **SC-007**: `docker build .` succeeds on an unmodified checkout after the change, confirming the Dockerfile installs GLPK and produces a working image.

## Assumptions

- GLPK 4.65 (released 2019) is used as the minimum version floor.  It has been available in all major Linux distribution package managers for several years and in Homebrew on macOS.  Earlier versions back to 4.45 also have genuine TLS but lack other API stability guarantees.
- The build target platforms are macOS (Homebrew GLPK) and Linux (apt/yum/dnf GLPK packages).  Windows (MinGW/MSYS2) is out of scope for this feature.
- The `glp_read_lp` / `glp_write_lp` API (CPLEX LP format) present in GLPK 4.65 is functionally equivalent to the old `lpx_read_cpxlp` / `lpx_write_cpxlp` for all LP files used in the dwsolver test suite and examples.
- The `glp_get_it_cnt` function, which replaced `lpx_get_int_parm(lp, LPX_K_ITCNT)`, returns the same iteration counter value.
