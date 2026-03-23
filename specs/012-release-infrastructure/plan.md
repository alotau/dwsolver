# Implementation Plan: Release Infrastructure

**Branch**: `012-release-infrastructure` | **Date**: 2026-03-22 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/012-release-infrastructure/spec.md`

## Summary

Enable reproducible, automated releases for dwsolver by: (1) ensuring `make distcheck` passes cleanly so a complete, self-contained source tarball can always be produced; (2) adding a GitHub Actions release workflow that runs `make distcheck` on Linux and publishes the tarball as a GitHub Release when a `v*` tag is pushed; and (3) documenting the full release process in the README. Implementation is primarily build-system metadata (`EXTRA_DIST`, `DISTCLEANFILES`) and CI configuration rather than C source changes.

## Technical Context

**Language/Version**: C99 (no source changes); GNU Autotools (automake 1.16+, autoconf 2.71+, libtool 2.4+); GitHub Actions YAML  
**Primary Dependencies**: GLPK ≥ 4.65 (external); `softprops/action-gh-release` (SHA-pinned); `actions/checkout` (SHA-pinned)  
**Storage**: N/A — artifacts are source tarballs (`.tar.gz`) and GitHub Release assets  
**Testing**: `make distcheck` (Autotools distcheck gate); `make check` (test_blas, test_lib_api unit tests + dw-tests.sh integration tests run inside the clean distcheck tree)  
**Target Platform**: Release workflow runs on Ubuntu latest; `make distcheck` verified locally on macOS (with named-semaphore workaround); Windows CI is separate and not a release gate  
**Project Type**: Build/release infrastructure (Autotools + GitHub Actions); no new library or CLI code  
**Performance Goals**: N/A  
**Constraints**: Ubuntu's `libglpk-dev` ships no `glpk.pc` — `GLPK_CFLAGS`/`GLPK_LIBS` must be exported to the environment; macOS `sem_init` is not implemented — distcheck on macOS requires `--enable-named-semaphores --enable-recursive-mutex`; `DISTCHECK_CONFIGURE_FLAGS` cannot carry values with spaces (make word-splits them)  
**Scale/Scope**: 4 files modified (`Makefile.am`, `tests/Makefile.am`, `Makefile.in`, `tests/Makefile.in`), 2 files added (`.github/workflows/release.yml`, README section); no new C source files

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | ✅ PASS | No changes to solver logic; `make check` (including dw-tests.sh) runs inside distcheck tree verifying all examples |
| II. Thread Safety | ✅ PASS | No threading changes |
| III. Cross-Platform Portability | ✅ PASS | Release workflow is Linux-only (per spec Assumptions); macOS and Windows continue to be covered by existing CI workflows |
| IV. Repair Before Extension | ✅ PASS | Branch is based on top of the completed repair work (PR #23 / `011-remove-embedded-glpk`); all platforms building correctly |
| V. CLI-First, Library-Ready | ✅ PASS | No changes to CLI or library interface |

**Gate result**: All five principles satisfied. Proceeding to Phase 0.

## Project Structure

### Documentation (this feature)

```text
specs/012-release-infrastructure/
├── plan.md          ← this file
├── research.md      ← Phase 0 output
├── data-model.md    ← Phase 1 output
├── contracts/
│   ├── release-workflow.md   ← Phase 1 output
│   └── tarball-contents.md  ← Phase 1 output
└── tasks.md         ← Phase 2 output (/speckit.tasks — not created here)
```

### Files Changed (repository root)

```text
Makefile.am                          ← EXTRA_DIST expanded
Makefile.in                          ← regenerated from Makefile.am
tests/
├── Makefile.am                      ← EXTRA_DIST + DISTCLEANFILES added
└── Makefile.in                      ← regenerated from tests/Makefile.am
.github/
└── workflows/
    └── release.yml                  ← new: release gate + GitHub Release creation
README.md                            ← "Cutting a Release" section added
```

**Structure Decision**: This feature is pure build/CI infrastructure — no new source directories. All changes are confined to build descriptors, CI YAML, and documentation.
