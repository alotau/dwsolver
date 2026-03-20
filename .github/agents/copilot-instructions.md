# dwsolver-repaired Development Guidelines

Auto-generated from all feature plans. Last updated: 2026-03-20

## Active Technologies
- File I/O — CPLEX LP format input, text output; no database (002-cross-platform-repair)
- YAML (GitHub Actions workflow syntax) + GitHub Actions; `actions/checkout@v4`; `msys2/setup-msys2@v2` (003-github-actions-ci)
- N/A — single workflow YAML file; no persisted state (003-github-actions-ci)
- Dockerfile (no language version constraint); underlying binary is C99 + Docker Engine ≥ 20.10; `ubuntu:24.04` base image; `build-essential` + `automake` (builder stage only) (004-docker-containerized-build)
- None in the image; user data accessed via volume mounts (`-v <host-dir>:/data`) (004-docker-containerized-build)
- C99 (targeting GCC ≥ 9 and Clang ≥ 12) + GLPK 4.44 (thread-patched, embedded in `src/`), pthreads (POSIX or winpthreads on Windows) (005-windows-build-fix)
- File I/O only — CPLEX LP format input (`lpx_read_cpxlp`), text output files (005-windows-build-fix)

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
- 005-windows-build-fix: Added C99 (targeting GCC ≥ 9 and Clang ≥ 12) + GLPK 4.44 (thread-patched, embedded in `src/`), pthreads (POSIX or winpthreads on Windows)
- 004-docker-containerized-build: Added Dockerfile (no language version constraint); underlying binary is C99 + Docker Engine ≥ 20.10; `ubuntu:24.04` base image; `build-essential` + `automake` (builder stage only)
- 003-github-actions-ci: Added YAML (GitHub Actions workflow syntax) + GitHub Actions; `actions/checkout@v4`; `msys2/setup-msys2@v2`


<!-- MANUAL ADDITIONS START -->
<!-- MANUAL ADDITIONS END -->
