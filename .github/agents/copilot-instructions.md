# dwsolver-repaired Development Guidelines

Auto-generated from all feature plans. Last updated: 2026-03-22

## Active Technologies
- File I/O — CPLEX LP format input, text output; no database (002-cross-platform-repair)
- YAML (GitHub Actions workflow syntax) + GitHub Actions; `actions/checkout@v4`; `msys2/setup-msys2@v2` (003-github-actions-ci)
- N/A — single workflow YAML file; no persisted state (003-github-actions-ci)
- Dockerfile (no language version constraint); underlying binary is C99 + Docker Engine ≥ 20.10; `ubuntu:24.04` base image; `build-essential` + `automake` (builder stage only) (004-docker-containerized-build)
- None in the image; user data accessed via volume mounts (`-v <host-dir>:/data`) (004-docker-containerized-build)
- C99 (targeting GCC ≥ 9 and Clang ≥ 12) + GLPK 4.44 (thread-patched, embedded in `src/`), pthreads (POSIX or winpthreads on Windows) (005-windows-build-fix)
- File I/O only — CPLEX LP format input (`lpx_read_cpxlp`), text output files (005-windows-build-fix)
- C99 (GCC and Clang; MinGW-w64 on Windows) + Embedded GLPK 4.44, POSIX pthreads, libc/libm (006-test-coverage)
- File-based (CPLEX `.lp` format, shell-script test runner) (006-test-coverage)
- Markdown / Mermaid diagram syntax (GitHub-native rendering, no version dependency) + None — Mermaid renders natively in GitHub Markdown as of 2022 (007-architecture-diagrams)
- Plain `.md` files committed to the repository (007-architecture-diagrams)
- C99 (GCC/Clang on macOS/Linux; MinGW-w64 on Windows) + POSIX Threads (pthreads); embedded GLPK 4.44 (treated as a library boundary — `src/glp*.c` not modified) (008-sei-cert-c-compliance)
- N/A (solver writes solution files; no database) (008-sei-cert-c-compliance)
- C (C99); shell for CI workflows + GNU Autotools (autoconf ≥ 2.60, automake ≥ 1.11), GLPK 4.44 (embedded) (009-repo-structure-cleanup)
- C99 + GLPK 4.44 (thread-patched, embedded in `src/`); POSIX Threads (pthreads); GNU Autotools + libtool (010-callable-library)
- File-based (LP problem files on disk; no database) (010-callable-library)
- C99 + GLPK ≥ 4.65 (external, via `pkg-config`), POSIX pthreads (011-remove-embedded-glpk)
- N/A (solver reads/writes LP files from disk) (011-remove-embedded-glpk)
- C99 (no source changes); GNU Autotools (automake 1.16+, autoconf 2.71+, libtool 2.4+); GitHub Actions YAML + GLPK ≥ 4.65 (external); `softprops/action-gh-release` (SHA-pinned); `actions/checkout` (SHA-pinned) (012-release-infrastructure)
- N/A — artifacts are source tarballs (`.tar.gz`) and GitHub Release assets (012-release-infrastructure)

- C99 (POSIX; GCC 9+ and Clang 12+ are primary targets) + POSIX Threads (pthreads); GLPK 4.44 (embedded, thread-patched) (002-cross-platform-repair)

## Project Structure

```text
src/
tests/
```

## Commands

# Add commands for C99 (POSIX; GCC 9+ and Clang 12+ are primary targets)

## Code Style

C99 (POSIX; GCC 9+ and Clang 12+ are primary targets): Follow standard conventions

## Recent Changes
- 012-release-infrastructure: Added C99 (no source changes); GNU Autotools (automake 1.16+, autoconf 2.71+, libtool 2.4+); GitHub Actions YAML + GLPK ≥ 4.65 (external); `softprops/action-gh-release` (SHA-pinned); `actions/checkout` (SHA-pinned)
- 011-remove-embedded-glpk: Added C99 + GLPK ≥ 4.65 (external, via `pkg-config`), POSIX pthreads
- 010-callable-library: Added C99 + GLPK 4.44 (thread-patched, embedded in `src/`); POSIX Threads (pthreads); GNU Autotools + libtool


<!-- MANUAL ADDITIONS START -->
<!-- MANUAL ADDITIONS END -->
