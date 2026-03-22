# Implementation Plan: Repository Top-Level Cleanup

**Branch**: `009-repo-structure-cleanup` | **Date**: 2026-03-21 | **Spec**: [spec.md](spec.md)

## Summary

Reduce top-level repository clutter by: (1) moving 7 autoconf auxiliary scripts to `build-aux/` via `AC_CONFIG_AUX_DIR`, (2) moving 6 GLPK attribution/patch files to `third-party/glpk/`, and (3) removing two committed editor backup files and adding `*~` to `.gitignore`. All CI workflows, the Dockerfile, and `Makefile.am` must be updated to reference the new paths. Build and tests must remain fully green on all three platforms.

## Technical Context

**Language/Version**: C (C99); shell for CI workflows  
**Primary Dependencies**: GNU Autotools (autoconf ≥ 2.60, automake ≥ 1.11), GLPK 4.44 (embedded)  
**Storage**: N/A  
**Testing**: `tests/dw-tests.sh` (shell, 7 test cases)  
**Target Platform**: macOS, Linux, Windows (MinGW64/MSYS2)  
**Project Type**: CLI tool + autotools build system  
**Performance Goals**: N/A — this is a build/repository structure change  
**Constraints**: Build must produce identical binary; all 7 tests must pass; CI must stay green on all 6 check types  
**Scale/Scope**: ~13 files moved/removed; 4 files updated (configure.ac, Makefile.am, Dockerfile, CI YAMLs)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | ✅ Pass | No solver logic touched; verified by test suite |
| II. Thread Safety | ✅ Pass | No threading code touched |
| III. Cross-Platform Portability | ✅ Pass | All three platform CI configs updated; build verified on each |
| IV. Repair Before Extension | ✅ Pass | This is a chore/cleanup, not an extension; build is currently passing |
| V. CLI-First, Library-Ready | ✅ Pass | CLI interface unchanged |

No gate violations.

## Project Structure

### Documentation (this feature)

```text
specs/009-repo-structure-cleanup/
├── plan.md              ← this file
├── research.md          ← Phase 0 output
├── data-model.md        ← N/A (no data model)
├── quickstart.md        ← Phase 1 output
└── tasks.md             ← Phase 2 output (speckit.tasks)
```

### Repository Changes

```text
# New directories created:
build-aux/               ← 7 autoconf aux scripts (moved from root)
third-party/
└── glpk/                ← 6 GLPK files + new README (moved from root)

# Files deleted from root:
compile  config.guess  config.sub  depcomp  install-sh  ltmain.sh  missing
GLPK_AUTHORS  GLPK_INSTALL  GLPK_NEWS  GLPK_README  GLPK_THANKS
glpk-4.44.ThreadReady.patch
configure~  config.h.in~

# Files modified:
configure.ac              ← add AC_CONFIG_AUX_DIR([build-aux])
Makefile.am               ← update EXTRA_DIST paths for GLPK files
Dockerfile                ← update config.guess/config.sub/ltmain.sh paths
.github/workflows/ci-linux.yml   ← update touch/path references
.github/workflows/ci-macos.yml   ← update touch/path references
.github/workflows/ci-windows.yml ← update touch/cp path references
.gitignore                ← add *~ pattern
README.md                 ← update project structure description
```

---

## Phase 0: Research

*All NEEDS CLARIFICATION items resolved by inspecting the live file system. No external research required.*

### Findings

**Finding 1 — `configure.ac` has no `AC_CONFIG_AUX_DIR`**
- Current: `AM_INIT_AUTOMAKE` only; no `AC_CONFIG_AUX_DIR`
- Decision: Add `AC_CONFIG_AUX_DIR([build-aux])` immediately before `AM_INIT_AUTOMAKE`
- Rationale: Standard autoconf idiom; supported since autoconf 2.60; no compatibility concerns
- Alternatives considered: Keeping aux scripts at root — rejected, defeats the purpose

**Finding 2 — `Makefile.am` lists GLPK files in `EXTRA_DIST` by root-relative name**
- Current (lines 7–8): `EXTRA_DIST = examples glpk-4.44.ThreadReady.patch GLPK_INSTALL GLPK_AUTHORS GLPK_README GLPK_NEWS GLPK_THANKS ...`
- Decision: Update paths to `third-party/glpk/glpk-4.44.ThreadReady.patch` etc.
- Rationale: `make dist` must continue to include the GLPK attribution files in the tarball

**Finding 3 — CI workflows `touch` and copy aux scripts at root**
- `ci-linux.yml` (lines 20, 39, 58): `touch aclocal.m4 configure config.h.in ltmain.sh ...`
- `ci-macos.yml` (line 17): same pattern
- `ci-windows.yml` (lines 34–36): same touch; also `cp ... config.guess` and `cp ... config.sub` to root
- Decision: Update all `touch` and `cp` targets to `build-aux/ltmain.sh`, `build-aux/config.guess`, `build-aux/config.sub`
- Rationale: After the move, the files live under `build-aux/`; touching the old paths would create empty stubs at root, not refresh the real files

**Finding 4 — `Dockerfile` copies `config.guess`/`config.sub` to root**
- Lines 20–23: `touch ... ltmain.sh` and `cp ... config.guess`, `cp ... config.sub`
- Decision: Update same as CI — target `build-aux/config.guess` and `build-aux/config.sub`
- Alternatives considered: Removing the refresh entirely — rejected, the Docker build may run on architectures where the committed config.guess is stale

**Finding 5 — `configure~` and `config.h.in~` are tracked in git**
- Confirmed by `git ls-files | grep '~'`
- Decision: `git rm configure~ config.h.in~`; add `*~` to `.gitignore`

**Finding 6 — GLPK files are not referenced by path in build scripts**
- Confirmed: `Makefile.am` lists them in `EXTRA_DIST` only (no build rule); `Dockerfile` and CI workflows do not reference them
- Decision: Safe to move to `third-party/glpk/` with only the `EXTRA_DIST` update needed

---

## Phase 1: Design & Contracts

### Task Breakdown

The implementation is a sequence of atomic, independently-verifiable steps. No new external interfaces are introduced.

#### Step A — Git hygiene (no build impact; do first)
1. `git rm configure~ config.h.in~`
2. Add `*~` to `.gitignore`
3. Verify `git status` is clean; commit: `chore(009): remove editor backup files`

#### Step B — Move GLPK files
1. `git mv` the 6 GLPK files to `third-party/glpk/`
2. Create `third-party/glpk/README` (brief attribution note)
3. Update `Makefile.am` `EXTRA_DIST` paths
4. Build: `./configure && make` — must succeed unchanged (GLPK files are not compiled)
5. Commit: `chore(009): move GLPK attribution files to third-party/glpk/`

#### Step C — `AC_CONFIG_AUX_DIR` + autoreconf (most impactful step)
1. Add `AC_CONFIG_AUX_DIR([build-aux])` to `configure.ac` before `AM_INIT_AUTOMAKE`
2. Run `autoreconf -fi` — regenerates aux scripts into `build-aux/`
3. `git rm` the old root-level aux scripts that are now in `build-aux/`
4. `git add build-aux/`
5. Build: `./configure && make` — must succeed with scripts in `build-aux/`
6. Commit: `chore(009): move autoconf aux scripts to build-aux/ via AC_CONFIG_AUX_DIR`

#### Step D — Update CI workflows and Dockerfile
1. Update all `touch` commands in `ci-linux.yml`, `ci-macos.yml`, `ci-windows.yml` to reference `build-aux/ltmain.sh`
2. Update `cp` commands in `ci-windows.yml` and `Dockerfile` to target `build-aux/config.guess` and `build-aux/config.sub`
3. Commit: `chore(009): update CI and Dockerfile paths for build-aux/`

#### Step E — Documentation
1. Update `README.md` project structure description to mention `build-aux/` and `third-party/`
2. Commit: `docs(009): update README for new repo structure`

### Interface Contracts

No external interfaces — this feature makes no API, CLI, or file-format changes visible to users. The `./configure && make` invocation is unchanged.

### Quickstart (for reviewer)

After merging:

```sh
# Fresh build — identical to before
./configure && make
cd tests && PATH="$PWD/../src:$PATH" ./dw-tests.sh  # expect 7/7 PASS

# Confirm aux scripts are in build-aux/
ls build-aux/   # compile  config.guess  config.sub  depcomp  install-sh  ltmain.sh  missing

# Confirm GLPK files moved
ls third-party/glpk/   # GLPK_AUTHORS  GLPK_INSTALL  GLPK_NEWS  GLPK_README  GLPK_THANKS  glpk-4.44.ThreadReady.patch  README

# Confirm no backup files
git ls-files | grep '~'    # (no output)
```

---

## Constitution Check (Post-Design)

All five principles remain satisfied. No new complexity introduced. The only risks are:
- CI `touch` commands missing the new path — mitigated by Step D
- `make dist` excluding GLPK files — mitigated by updating `EXTRA_DIST` in Step B
- `autoreconf` leaving stale copies at root — mitigated by explicit `git rm` of root aux scripts in Step C

