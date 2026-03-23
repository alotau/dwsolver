<!--
SYNC IMPACT REPORT
==================
Version change: 1.0.1 → 1.0.2
Changed: Technology Stack — Language target updated from C99 to C11 (strict superset;
         all existing C99-conforming code is valid C11; no behavioural changes).
Added sections: N/A
Removed sections: N/A
Templates reviewed:
  ✅ .specify/templates/plan-template.md   — generic Constitution Check gate; no changes needed
  ✅ .specify/templates/spec-template.md   — fully generic; no changes needed
  ✅ .specify/templates/tasks-template.md  — fully generic; no changes needed
Follow-up TODOs: none
-->

# DWSOLVER Constitution

## Core Principles

### I. Correctness First (NON-NEGOTIABLE)

Mathematical correctness of Dantzig-Wolfe solutions is the highest priority.
Every change that touches solver logic MUST be verified by running the full
examples suite (`tests/dw-tests.sh`) and confirming all expected optima are
reproduced before the branch is considered ready to merge. Numerical
tolerances MUST NOT be relaxed to paper over failures. Any result deviation,
even seemingly minor floating-point drift, MUST be investigated and resolved,
not suppressed.

### II. Thread Safety is Mandatory

Parallelism using POSIX threads (pthreads) is a core design requirement, not
an optional optimization. All shared state MUST be protected with appropriate
synchronization primitives (mutexes, semaphores). New code MUST NOT introduce
data races, double-frees, or undefined behavior under concurrent access.
Thread-related changes MUST be tested with thread sanitizer (TSan) where
tooling supports it. The principle that "multithreaded C from years ago is the
most dangerous game" applies: changes to threading logic require extra
scrutiny and explicit justification.

### III. Cross-Platform Portability

The codebase MUST build and produce correct results on macOS, Linux, and
Windows. No platform-specific APIs, compiler extensions, or linker flags may
be introduced without providing equivalent implementations or guards for all
three targets. POSIX compatibility headers MUST be conditionally included
where needed. CI MUST validate all supported platforms before a branch is
eligible to merge into main.

### IV. Repair Before Extension

The current state of the codebase is known-broken in measurable ways (Ubuntu
link failures, variable definition issues in `dw.h`). No new features MUST be
added until the existing suite of examples runs correctly on all target
platforms. Every PR MUST state clearly whether it is a repair or an extension.
Extensions require a passing repair baseline first.

### V. CLI-First, Library-Ready

The primary interface is a command-line tool. All user-facing behavior MUST
be reachable from the CLI. Internally, solver logic SHOULD be structured so
that it can be extracted into a callable library in the future without
requiring a rewrite. This means: avoid global mutable state in solver
functions, keep I/O paths separate from computation paths, and use clear
function boundaries.

## Technology Stack

- **Language**: C (targeting C11 for portability across platforms and compilers; C11 is a strict superset of C99 — all existing code remains valid)
- **Parallelism**: POSIX Threads (`pthreads`); on Windows via a pthreads-win32
  compatibility layer
- **LP Backend**: System GLPK ≥ 4.65 (external shared library, detected via
  `pkg-config`; `libglpk-dev` on Debian/Ubuntu, `glpk-devel` on Fedora/RHEL,
  `brew install glpk` on macOS)
- **Build System**: GNU Autotools (`configure`, `Makefile.am`); CMake may be
  introduced later for Windows support
- **Testing**: Shell-based example runner (`tests/dw-tests.sh`); examples live
  in `examples/` organized by source textbook or reference
- **Compiler Requirements**: GCC or Clang on POSIX platforms; MSVC or
  MinGW-w64 on Windows
- **External runtime dependencies**: GLPK ≥ 4.65 shared library; C standard
  library and pthreads

## Git Workflow

- **Never commit directly to `main`**. All work MUST happen on a feature or
  fix branch.
- **Branch naming**: `fix/<short-description>` for repairs, `feat/<short-description>`
  for new features, `chore/<short-description>` for non-functional changes
  (e.g., `fix/ubuntu-linker-errors`, `fix/dw-h-variable-definitions`).
- **Before pushing**: build locally, run `tests/dw-tests.sh`, confirm no
  regressions. Run thread sanitizer if the change touches concurrency.
- **Before opening a PR**: merge (or rebase onto) the current `main` to ensure
  the branch is up to date and conflicts are resolved locally.
- **Commit messages**: imperative mood, present tense, ≤72 chars subject line
  (e.g., `fix: resolve multiple-definition errors in dw.h on Linux`).
- **PRs** MUST describe: what is broken, what was changed, how it was tested,
  and which platforms were verified.

## Governance

This constitution supersedes all other informal practices. Amendments require:

1. A PR against `.specify/memory/constitution.md` with the proposed change.
2. The Sync Impact Report (HTML comment at top of file) updated to reflect the
   version bump, changed principles, and affected templates.
3. Version bump following semantic versioning:
   - **MAJOR**: Removal or redefinition of an existing principle.
   - **MINOR**: New principle or section added.
   - **PATCH**: Clarification, wording fix, or non-semantic refinement.
4. All affected template files updated in the same PR.

All feature specifications and implementation plans MUST include a
"Constitution Check" gate that explicitly verifies compliance with the five
Core Principles before implementation begins.

**Version**: 1.0.2 | **Ratified**: 2026-03-19 | **Last Amended**: 2026-03-23
