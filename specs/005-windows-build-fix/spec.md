# Feature Specification: Code Quality Improvements — Bug Fixes, Duplication Reduction, and Clarity

**Feature Branch**: `chore/code-quality` *(see Note below)*
**Created**: 2026-03-20
**Status**: Draft
**Note**: This spec was scaffolded on branch `005-windows-build-fix` but the
constitution requires a dedicated branch. This work must be executed on a branch
named `chore/code-quality` (or `fix/code-quality-bugs` if scoped to bugs only)
once `005-windows-build-fix` is merged to master.

## Background

A codebase analysis performed 2026-03-20 identified a set of bugs, code clarity
issues, and architectural concerns. This spec captures the actionable items
derived from that analysis, ordered by the priority matrix:

| Issue | Severity | Correctness Risk | Effort |
|-------|----------|-----------------|--------|
| `val = NULL` passed to `glp_get_mat_col` in phase 2 | Bug | High (UB) | Low |
| `const char*` memory leak in `rounding_thread` | Bug | Medium | Low |
| `glpk_mutex` held over logging loop | Perf | None | Low |
| Dead commented-out code blocks | Clarity | None | Low |
| `phase_1` / `phase_2` duplication | Maintainability | Medium | Medium |
| `master_mutex` serializing read-only access | Perf | None | Medium |
| `main()` monolith | Clarity | Low | Medium |
| Pre-allocating 3000 slots per subproblem | Memory | None | Medium |
| Global state (`master_lp`, `D`, `signals`) | Architecture | High | High |
| Mixed `printf`/`dw_printf` | Clarity/correctness | Low | Medium |

**Constraints (non-negotiable):**
- Every change must leave `tests/dw-tests.sh` fully passing on macOS and Linux.
- Numerical tolerances and convergence behaviour MUST NOT change.
- Thread safety MUST be preserved or improved.

---

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Fix Undefined Behaviour in `phase_2_iteration` (Priority: P1)

A developer running the solver with address sanitizer or on any platform where
`glp_get_mat_col(master_lp, col, ind, NULL)` traps (GLPK documents `val` as
required) should not see a crash or silent data corruption.

**Why this priority**: Undefined behaviour in the hot path of Phase II of the
algorithm is a latent correctness hazard. It may silently corrupt memory,
affecting solution values without a visible crash.

**Independent Test**: Run the full test suite under AddressSanitizer
(`CFLAGS="-fsanitize=address" ./configure && make`) and confirm no ASan errors
in the column-addition block of `phase_2_iteration`.

**Acceptance Scenarios**:

1. **Given** the Phase II master iteration is executing, **When** a new column
   is added and `glp_get_mat_col` is called, **Then** both `ind` and `val`
   arrays are valid allocated buffers and the call is safe.
2. **Given** the column-addition diagnostic block in `phase_2_iteration`,
   **When** the block executes, **Then** `ind` and `val` are freed before the
   enclosing mutex scope exits.

---

### User Story 2 — Fix Memory Leak in `rounding_thread` (Priority: P1)

A developer running with a memory analyser (Valgrind, LeakSanitizer) on a run
that uses the `--round` flag should see zero leaks attributable to the rounding
code path.

**Why this priority**: `rounding_thread` is called once per subproblem per
invocation of the rounding pass. Any leak is multiplied by the number of
subproblems.

**Independent Test**: Run an example that exercises rounding (`dwsolver -r ...`)
under Valgrind or LSan and confirm no `rounding_thread` allocation remains
unreleased at process exit.

**Acceptance Scenarios**:

1. **Given** `rounding_thread` is entered, **When** the thread exits,
   **Then** the `malloc` that initialises `local_col_name` is freed before
   `pthread_exit`.
2. **Given** `const char* local_col_name = malloc(...)`, **When** the const
   qualifier is used, **Then** the declaration is changed to `char*` to match
   the allocated type and compiler warnings are resolved.

---

### User Story 3 — Eliminate `phase_1` / `phase_2` Code Duplication (Priority: P2)

A developer fixing a bug in the reduced-cost test (`r - obj > TOLERANCE`) must
only change one place instead of two.

**Why this priority**: The two functions share approximately 60% of their body.
Any bug present in one is almost certainly present in the other, and any fix
must be applied twice — this has already caused divergence (e.g., the
`TOLERANCE` check is `> TOLERANCE` in phase 1 and `> 0.0 + TOLERANCE` in
phase 2).

**Independent Test**: After unification, the shared body exists in exactly one
location; running the full test suite confirms no regression.

**Acceptance Scenarios**:

1. **Given** the column-acceptance test `r - obj > TOLERANCE`, **When** the
   unified function is called for both Phase I and Phase II, **Then** the test
   uses a single consistent expression.
2. **Given** the Phase I `first_run` branch (y_accumulator initialisation, sign
   correction, `obj_names`/`obj_coefs` tracking), **When** the unified function
   runs with `phase == 1` and `first_run == true`, **Then** those branches are
   entered correctly.
3. **Given** the full test suite, **When** run after the refactoring,
   **Then** all tests pass with identical output to pre-refactoring runs.

---

### User Story 4 — Remove Dead Code and Standardise Output Logging (Priority: P3)

A developer reading `dw_phases.c` or `dw_rounding.c` for the first time should
not have to mentally filter out large blocks of commented-out code.  All
diagnostic output should be suppressible by the `--silent` flag.

**Why this priority**: Clarity improvements lower the risk of future bugs; log
standardisation prevents diagnostic messages from bypassing the verbosity gate.

**Independent Test**: After cleanup, `grep -n "//\s*printf\|//\s*glp_" src/dw_phases.c src/dw_rounding.c` returns no substantive commented-out blocks (stray single-line notes are acceptable).

**Acceptance Scenarios**:

1. **Given** `--silent` mode, **When** the solver runs, **Then** no solver
   diagnostic output appears on stdout (GLPK output is already redirected to
   file; the remaining `printf` calls should go through `dw_printf`).
2. **Given** dead commented-out blocks in `dw_phases.c` and `dw_rounding.c`,
   **When** they are removed, **Then** the files have no blocks of more than
   3 consecutive commented-out lines.

---

### User Story 5 — Reduce Lock Contention on `glpk_mutex` and `master_mutex` (Priority: P3)

A developer profiling a run with many subproblems should see reduced lock
wait time on `glpk_mutex` during subproblem logging loops.

**Why this priority**: Both mutex over-use cases are low-risk changes with
measurable performance benefit at scale (many subproblems, verbose logging).

**Independent Test**: Confirm the `glpk_mutex` logging loop in `prepare_column`
is only acquired if the verbosity threshold is met; confirm that `original_master_lp` column reads during subproblem setup can proceed in parallel.

**Acceptance Scenarios**:

1. **Given** `prepare_column` with `verbosity < OUTPUT_ALL`, **When** it is
   called, **Then** `glpk_mutex` is never acquired (the block is skipped
   entirely).
2. **Given** the `original_master_lp` read loop in `subproblem_thread`,
   **When** a `pthread_rwlock_rdlock` (or equivalent) is used, **Then**
   multiple subproblem threads can perform their column mapping in parallel.

---

### User Story 6 — Move Core Solver State Out of True Globals (Priority: P4)

A library consumer (future use case) should be able to create two independent
solver instances in the same process without data collisions.

**Why this priority**: This is a large architectural change. It is scoped to
P4 because it has no immediate user-visible impact and carries the highest
risk of introducing a regression. It should only be attempted after P1–P3
items are complete and the test suite is green.

**Independent Test**: Write a test that calls `main`-equivalent setup/teardown
twice sequentially and confirm both runs produce correct optima independently.
(Full parallel isolation is out of scope for this spec.)

**Acceptance Scenarios**:

1. **Given** `master_lp`, `original_master_lp`, `D`, and `signals` are moved
   into `faux_globals`, **When** the solver runs any example, **Then** all
   tests pass.
2. **Given** the mutexes and semaphore are moved into a new `sync_state` struct
   passed through `faux_globals`, **When** the solver runs, **Then** no
   thread sanitizer errors are reported.
