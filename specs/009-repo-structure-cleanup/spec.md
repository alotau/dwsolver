# Feature Specification: Repository Top-Level Cleanup

**Feature Branch**: `009-repo-structure-cleanup`
**Created**: 2026-03-21
**Status**: Draft

## User Scenarios & Testing *(mandatory)*

### User Story 1 — New contributor clones and builds without confusion (Priority: P1)

A developer visits the repository for the first time, clones it, and runs `./configure && make`. They see a tidy root directory that makes the project structure immediately understandable, without being overwhelmed by autoconf scaffolding scripts, GLPK attribution files, or build artifacts.

**Why this priority**: First impressions matter. A cluttered root is the first thing any contributor sees. The build must still work correctly after restructuring — this is the gate for every other story.

**Independent Test**: Clone a fresh copy of the repo, run `./configure && make` from root, and verify the binary produces correct output on at least one example. Count root-level tracked files before and after — after should be at least 13 fewer.

**Acceptance Scenarios**:

1. **Given** a fresh clone, **When** `./configure && make` is run, **Then** the build succeeds and produces `src/dwsolver` (or `src/dwsolver.exe` on Windows) with no extra steps.
2. **Given** the built repository, **When** the test suite is run, **Then** all existing tests continue to pass.
3. **Given** the root directory after the refactor, **When** a contributor lists its contents, **Then** no autoconf auxiliary scripts (`compile`, `config.guess`, `config.sub`, `depcomp`, `install-sh`, `ltmain.sh`, `missing`) appear at root — they are in `build-aux/`.

---

### User Story 2 — Third-party attribution is findable but not intrusive (Priority: P2)

A contributor who wants to understand the GLPK dependency can easily find its licence and attribution files. Those files do not clutter the root alongside the project's own top-level files.

**Why this priority**: The files must remain in the repository (legal obligation), but their presence at root mixed with project files causes confusion about what is "this project" vs. "vendored dependency".

**Independent Test**: Verify `third-party/glpk/` exists and contains the six GLPK files, while the root no longer contains any `GLPK_*` files or the patch file.

**Acceptance Scenarios**:

1. **Given** the refactored repository, **When** a contributor looks for GLPK attribution, **Then** all GLPK files (`GLPK_AUTHORS`, `GLPK_INSTALL`, `GLPK_NEWS`, `GLPK_README`, `GLPK_THANKS`, `glpk-4.44.ThreadReady.patch`) are under `third-party/glpk/`.
2. **Given** the refactored repository, **When** a contributor lists the root, **Then** no `GLPK_*` files or GLPK patch files appear there.
3. **Given** the `third-party/glpk/` directory, **When** it is inspected, **Then** a brief `README` explains what the files are and why they are present.

---

### User Story 3 — Editor backup and generated files stay out of git (Priority: P3)

A contributor working in the repository never sees `configure~`, `config.h.in~`, or editor `*~` backup files show up in `git status` or `git ls-files`.

**Why this priority**: The `.gitignore` already covers most build artifacts, but `*~` backup files are not currently ignored, and `configure~` / `config.h.in~` are actively committed. These are low-risk, high-visibility fixes.

**Independent Test**: Can be tested in isolation — remove the two committed backup files, add `*~` to `.gitignore`, confirm `git status` is clean after a build and a simulated editor save.

**Acceptance Scenarios**:

1. **Given** the refactored repository, **When** `git ls-files` is run, **Then** no `*~` files appear.
2. **Given** a clean build, **When** `git status` is run, **Then** no untracked `*~` files appear.
3. **Given** the `.gitignore`, **When** it is inspected, **Then** `*~` is listed as an ignored pattern.

---

### Edge Cases

- Does the `Dockerfile` reference any of the files being moved (GLPK files, aux scripts)? It must not break.
- Do any CI workflow files reference moved files by root-relative path? They must be updated if so.
- Does `make dist` / `make distcheck` still work after moving auxiliary scripts to `build-aux/`? This must be verified.
- Does `autoreconf` need to be re-run and its outputs re-committed after the `AC_CONFIG_AUX_DIR` change? Yes — the regenerated files replace the originals.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The repository root MUST NOT contain the autoconf auxiliary scripts `compile`, `config.guess`, `config.sub`, `depcomp`, `install-sh`, `ltmain.sh`, and `missing` — they MUST reside in `build-aux/`.
- **FR-002**: `configure.ac` MUST declare `AC_CONFIG_AUX_DIR([build-aux])` so that `autoreconf` places auxiliary scripts under `build-aux/` automatically on future runs.
- **FR-003**: The repository root MUST NOT contain `GLPK_AUTHORS`, `GLPK_INSTALL`, `GLPK_NEWS`, `GLPK_README`, `GLPK_THANKS`, or `glpk-4.44.ThreadReady.patch` — they MUST reside in `third-party/glpk/`.
- **FR-004**: A `third-party/glpk/README` file MUST exist explaining that the directory contains GLPK attribution and patch files bundled with the original source distribution.
- **FR-005**: The committed files `configure~` and `config.h.in~` MUST be removed from the repository (via `git rm`).
- **FR-006**: The `.gitignore` MUST include a `*~` pattern to prevent editor backup files from appearing as untracked.
- **FR-007**: The project MUST build successfully (`./configure && make`) on macOS, Linux, and Windows (MinGW/MSYS2) after the refactor.
- **FR-008**: All existing tests MUST continue to pass after the refactor.
- **FR-009**: All CI workflows MUST continue to pass after the refactor with no broken file paths.
- **FR-010**: The `README.md` project structure section MUST be updated to reflect the new `build-aux/` and `third-party/` directories.

### Key Entities

- **`build-aux/`**: New subdirectory holding the 7 autoconf-generated auxiliary scripts currently at root.
- **`third-party/glpk/`**: New subdirectory holding the 6 GLPK attribution/patch files currently at root, plus a brief `README`.
- **`configure.ac`**: Modified to add `AC_CONFIG_AUX_DIR([build-aux])`.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The count of tracked, non-hidden, non-directory files at repository root is reduced by at least 13 (7 aux scripts + 6 GLPK files removed from root).
- **SC-002**: `./configure && make` succeeds in a clean working directory on macOS, Linux, and Windows CI — verified by all CI checks passing green.
- **SC-003**: All test suite results remain identical to the pre-refactor baseline (7/7 pass).
- **SC-004**: `git ls-files | grep -E '~$'` returns zero results.
- **SC-005**: All 6 CI check types (macOS, Linux, Linux ASan+UBSan, Linux TSan, Docker, Windows) are green after the refactor merges to `main`.

## Assumptions

- The GLPK files at root are not referenced by any build script path and can be moved without build impact. This will be confirmed by inspecting `Makefile.am`, `configure.ac`, and the `Dockerfile` before moving.
- `AC_CONFIG_AUX_DIR` is supported by the version of autoconf in use (true for autoconf ≥ 2.60; the project's `configure.ac` will be checked to confirm).
- `autoreconf` will be re-run locally after the `configure.ac` change; the regenerated auxiliary scripts in `build-aux/` replace the originals and are committed.
- The `.vscode/`, `.specify/`, and `.github/` directories stay at root — they are intentional tooling directories.
- The project-standard top-level prose files (`AUTHORS`, `COPYING`, `NEWS`, `ChangeLog`, `INSTALL`, `ADDITIONAL_LICENSE_TERMS`) stay at root — GNU coding standards and `make dist` expect them there.

