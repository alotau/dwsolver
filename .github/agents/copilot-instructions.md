# dwsolver-repaired Development Guidelines

Auto-generated from all feature plans. Last updated: 2026-03-19

## Active Technologies
- File I/O — CPLEX LP format input, text output; no database (002-cross-platform-repair)
- YAML (GitHub Actions workflow syntax) + GitHub Actions; `actions/checkout@v4`; `msys2/setup-msys2@v2` (003-github-actions-ci)
- N/A — single workflow YAML file; no persisted state (003-github-actions-ci)

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
- 003-github-actions-ci: Added YAML (GitHub Actions workflow syntax) + GitHub Actions; `actions/checkout@v4`; `msys2/setup-msys2@v2`
- 002-cross-platform-repair: Added C99 (POSIX; GCC 9+ and Clang 12+ are primary targets) + POSIX Threads (pthreads); GLPK 4.44 (embedded, thread-patched)

- 002-cross-platform-repair: Added C99 (POSIX; GCC 9+ and Clang 12+ are primary targets) + POSIX Threads (pthreads); GLPK 4.44 (embedded, thread-patched)

<!-- MANUAL ADDITIONS START -->
<!-- MANUAL ADDITIONS END -->
