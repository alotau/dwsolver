# Research: Remove Embedded GLPK Source

**Feature**: 011-remove-embedded-glpk  
**Phase**: 0 — Outline & Research  
**Date**: 2026-03-22

---

## Decision 1: GLPK Version Floor

**Decision**: Require GLPK ≥ 4.65.

**Rationale**: GLPK 4.45 (2010) introduced genuine thread-local storage for the
internal `ENV *` pointer, making the `glpenv02.c` workaround unnecessary.
However, 4.45–4.64 are all effectively obsolete and unavailable in current
package managers. GLPK 4.65 (released October 2019) is the first version that
is:

- packaged as `libglpk-dev` / `libglpk40` on Ubuntu 20.04+ and Debian 10+
- available via `brew install glpk` on macOS (Homebrew formula tracks 4.65+)
- ships a well-formed `glpk.pc` pkg-config file suitable for `PKG_CHECK_MODULES`

GLPK 5.0 (released 2020) is also available on some platforms but is not yet
universal enough to be a hard floor.

**Alternatives considered**:
- Floor at 4.45: Too old; not in any current package manager. Would require
  users to build from source.
- Floor at 5.0: Too new; Ubuntu 20.04 LTS (supported until 2025) only ships 4.65.

---

## Decision 2: Build System Detection Mechanism

**Decision**: Use `PKG_CHECK_MODULES([GLPK], [glpk >= 4.65])` in `configure.ac`.

**Rationale**: `pkg-config` is the standard autotools mechanism for external
library detection on POSIX. GLPK 4.65 ships `glpk.pc` on all major Linux
distributions and macOS Homebrew. `PKG_CHECK_MODULES` automatically sets
`GLPK_CFLAGS` and `GLPK_LIBS` (including transitive deps such as GMP or zlib
if the installed GLPK was compiled with them), eliminating the manual
`--with-gmp / --with-zlib` flags that previously applied only to the embedded
copy. The macro also produces a clear configure-time error when GLPK is absent
or too old.

**Alternatives considered**:
- `AC_CHECK_LIB([glpk], [glp_create_prob])`: Detects existence but does not
  check version, does not propagate transitive deps, does not set CFLAGS.
- CMake `find_package(GLPK)`: Appropriate only if build system is migrated to
  CMake, which is out of scope.

---

## Decision 3: API Migration for Deprecated `lpx_*` Functions

**Decision**: Replace each deprecated call with the direct modern equivalent.

### Mapping

| Deprecated (GLPK 4.44) | Modern (GLPK 4.65+) | Notes |
|---|---|---|
| `lpx_read_cpxlp(fname)` | `glp_create_prob()` then `glp_read_lp(lp, NULL, fname)` | Returns `glp_prob *`; `NULL` params argument uses default format options |
| `lpx_write_cpxlp(lp, fname)` | `glp_write_lp(lp, NULL, fname)` | `NULL` params argument uses defaults |
| `lpx_get_int_parm(lp, LPX_K_ITCNT)` | `glp_get_it_cnt(lp)` | Returns the simplex iteration counter; same semantics |

**Call site inventory** (from code analysis):

| File | Line(s) | Old Call | New Call |
|------|---------|----------|----------|
| `src/dw_solver.c` | 184 | `lpx_read_cpxlp(globals->master_name)` | `glp_create_prob()` + `glp_read_lp(lp, NULL, globals->master_name)` |
| `src/dw_solver.c` | 434 | `lpx_write_cpxlp(master_lp, "pre_master.cpxlp")` | `glp_write_lp(master_lp, NULL, "pre_master.cpxlp")` |
| `src/dw_solver.c` | 477 | `lpx_write_cpxlp(master_lp, local_buffer)` | `glp_write_lp(master_lp, NULL, local_buffer)` |
| `src/dw_solver.c` | 493 | `lpx_write_cpxlp(master_lp, local_buffer)` | `glp_write_lp(master_lp, NULL, local_buffer)` |
| `src/dw_solver.c` | 645 | `lpx_write_cpxlp(master_lp, local_buffer)` | `glp_write_lp(master_lp, NULL, local_buffer)` |
| `src/dw_subprob.c` | 103 | `lpx_read_cpxlp(my_data->infile_name)` | `glp_create_prob()` + `glp_read_lp(lp, NULL, my_data->infile_name)` |
| `src/dw_phases.c` | 228, 229 | `lpx_get_int_parm(master_lp, LPX_K_ITCNT)` | `glp_get_it_cnt(master_lp)` |
| `src/dw_phases.c` | 460, 461 | `lpx_get_int_parm(master_lp, LPX_K_ITCNT)` | `glp_get_it_cnt(master_lp)` |

**Important**: The two `lpx_read_cpxlp` call sites each return a `glp_prob *`
and immediately use it.  With the new API, `glp_create_prob()` allocates the
struct and `glp_read_lp()` populates it in-place, returning 0 on success or
non-zero on error.  The existing NULL-check on the returned pointer must be
replaced with a return-code check plus an explicit `glp_delete_prob()` on
failure.

**Alternatives considered**:
- Keep the `lpx_*` wrappers by providing a local compatibility shim header:
  Rejected — this would still require the glpk internal headers and defeats
  the purpose of the refactoring.

---

## Decision 4: GLPK Thread Safety in Modern Versions

**Decision**: No change to `glpk_mutex` guards is required for correctness.
The mutex remains appropriate for a different reason than originally intended.

**Rationale**:
- GLPK ≥ 4.45: each thread has its own `ENV *` in genuine TLS.  The original
  thread-safety problem (shared global memory-tracking state) is gone.
- However, `glp_read_lp` and `glp_write_lp` (and the old `lpx_read_cpxlp` /
  `lpx_write_cpxlp`) invoke GLPK's internal file I/O which uses `getenv`,
  locale functions, and format-string functions — all of which are not
  re-entrant under POSIX in all glibc versions.
- The existing `glpk_mutex` was already guarding only file-I/O calls (not
  solver calls), which is the correct granularity.  The guards remain valid
  and should not be removed.

---

## Decision 5: CI Platform Installation Commands

**Decision**: Use standard package manager steps.

| Platform | Step to add | Package |
|---|---|---|
| Ubuntu (GitHub Actions `ubuntu-latest`) | `apt-get install -y libglpk-dev` | `libglpk-dev` (pulls in `libglpk40`) |
| macOS (GitHub Actions `macos-latest`) | `brew install glpk` | `glpk` (Homebrew formula) |
| Docker builder stage | `apt-get install -y libglpk-dev` | `libglpk-dev` |
| Docker runner stage | `apt-get install -y libglpk40` | runtime `.so` only |

Note: All three Linux CI jobs (`linux`, `linux-asan-ubsan`, `linux-tsan`) in
`ci-linux.yml` need the install step added before the configure/build step.

---

## Decision 6: pkg-config File (`dwsolver.pc.in`)

**Decision**: Add `Requires.private: glpk` (not `Requires`).

**Rationale**: `dw_solver.h` exposes zero GLPK types.  Downstream consumers
only need `libglpk.so` at link time (because `libdwsolver.so` links it), not
`glpk.h` at compile time.  `Requires.private` communicates the runtime
dependency to `pkg-config --libs --static dwsolver` without forcing consumers
to install `libglpk-dev`.

---

## Decision 7: Constitution Update

**Decision**: Update the Technology Stack entry in
`.specify/memory/constitution.md` to record GLPK as an external dependency
(not embedded) and reflect the minimum version.

The line currently reads:
```
**LP Backend**: Embedded GLPK 4.44 (thread-patched variant included in `src/`)
```
It should become:
```
**LP Backend**: System GLPK ≥ 4.65 (external shared library, detected via pkg-config)
```
Also: the `**No external runtime dependencies**` line at the end of the
Technology Stack section must be removed or amended to reflect that GLPK is
now an external dependency.
