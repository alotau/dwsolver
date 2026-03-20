# Feature Specification: As-Built — DWSOLVER Core Implementation

**Feature Branch**: `docs/001-as-built-spec`
**Created**: 2026-03-19
**Status**: As-Built (describes existing implementation)
**Version described**: 1.2.1 (from `configure.ac`)

---

## Purpose of This Document

This specification describes what DWSOLVER **already does** as of the initial
repair baseline. It is not a design for new work; it is ground truth for the
state of the code so that repair tasks, regression tests, and future feature
specs can reference agreed-upon behavior.

---

## User Scenarios & Testing

### User Story 1 — Solve a block-angular LP via Dantzig-Wolfe Decomposition (Priority: P1)

A user has a linear program that can be expressed in block-angular form: a
set of coupling constraints (the "master" problem) and one or more
independent subproblems. They provide the master and subproblem LP files in
CPLEX LP format, along with a guide file, and run `dwsolver`. The solver
decomposes the problem, solves the LP relaxation to optimality using the
Dantzig-Wolfe algorithm, and writes the relaxed solution to `relaxed_solution`.

**Why this priority**: This is the entire purpose of the software.

**Independent Test**: Run any example in `examples/` with the correct guidefile
and compare the output `relaxed_solution` (sorted) against the corresponding
`ex_relaxed_solution`. This is exactly what `tests/dw-tests.sh` does.

**Acceptance Scenarios**:

1. **Given** the Bertsimas textbook example files in `examples/book_bertsimas/`,
   **When** `dwsolver -g guidefile` is run from that directory,
   **Then** `relaxed_solution` (sorted) matches `ex_relaxed_solution` with
   `x1=2.0, x2=1.5, x3=2.0` and the objective optimum of -16.5.

2. **Given** the Lasdon textbook example in `examples/book_lasdon/`,
   **When** `dwsolver --no-write-final-master --quiet -g guidefile` is run,
   **Then** `relaxed_solution` (sorted) matches `ex_relaxed_solution`.

3. **Given** the Mitchell web example in `examples/web_mitchell/`,
   **When** `dwsolver --no-write-final-master --quiet -g guidefile` is run,
   **Then** `relaxed_solution` (sorted) matches `ex_relaxed_solution` with
   optimal value -5.

4. **Given** the Trick web example in `examples/web_trick/`,
   **When** `dwsolver --no-write-final-master --quiet -g guidefile` is run,
   **Then** `relaxed_solution` (sorted) matches `ex_relaxed_solution`.

---

### User Story 2 — Solve with multiple parallel subproblems (Priority: P1)

The Bertsimas example can be run with 1, 2, or 3 subproblems using
`guidefile`, `guide_double`, or `guide_triple` respectively. The LP relaxation
optimum MUST be identical regardless of the number of subproblems used; only
the variable assignments across subproblems may differ in multi-subproblem
configurations where variables are not uniquely owned.

**Why this priority**: Parallel subproblem solving is the core value proposition
of the software. Correctness under parallelism must be assured.

**Independent Test**: Run `dwsolver -g guide_double` and `dwsolver -g guide_triple`
in `examples/book_bertsimas/` and verify the objective value matches the
single-subproblem run.

**Acceptance Scenarios**:

1. **Given** `guide_double` (2 subproblems) in `examples/book_bertsimas/`,
   **When** `dwsolver -g guide_double` is run,
   **Then** the optimal objective value equals that of the single-subproblem run.

2. **Given** `guide_triple` (3 subproblems) in `examples/book_bertsimas/`,
   **When** `dwsolver -g guide_triple` is run,
   **Then** the optimal objective value equals that of the single-subproblem run.

---

### User Story 3 — Traffic Flow Management (four_sea) example (Priority: P2)

The `four_sea` example models air traffic flow between airports with capacity
constraints. Four subproblems represent four aircraft/sectors. The optimal
total delay is deterministic (12 minutes), but the assignment of which two of
the eight flights are delayed is non-deterministic due to thread scheduling.

**Why this priority**: This is the domain for which the solver was originally
written and demonstrates non-trivial real-world applicability. The non-
deterministic assignment makes automated pass/fail testing harder, which is
a known gap.

**Independent Test**: Run `dwsolver -g guidefile` in `examples/four_sea/` and
verify the objective value (total delay) equals 12.

**Acceptance Scenarios**:

1. **Given** the four_sea problem files,
   **When** `dwsolver -g guidefile` is run,
   **Then** the optimal objective value equals 12 (total delay = 12 minutes).

2. **Given** multiple runs of the same problem,
   **When** comparing the variable assignments between runs,
   **Then** individual flight delay assignments MAY differ, but the total
   MUST always be 12.

---

### User Story 4 — Dantzig and book_dantzig examples (Priority: P2)

The `book_dantzig` example (from Dantzig & Thapa, Linear Programming, 2003,
example 10.5/10.6) produces multiple valid solutions due to thread ordering.
Upper bounds of 65 are added to subproblem variables because dwsolver does
not support unbounded subproblems.

**Acceptance Scenarios**:

1. **Given** the book_dantzig problem files,
   **When** `dwsolver -g book_test` is run,
   **Then** the solver produces a feasible optimal solution (objective value
   matches the reference, even if variable assignments differ).

---

### Edge Cases

- What happens when fewer than 1 or more than 26,000 subproblems are specified?
  → Solver prints an error and exits with return code 1.
- What happens if a subproblem column name does not appear in the master?
  → Solver prints `"PROBLEM: Column X not found in the master problem"` and
  aborts.
- What happens if a subproblem is unbounded?
  → Currently unsupported. Code prints `"Unbounded subproblems are not supported."`
  Behaviour is undefined and results should not be trusted.
- What happens if the guide file cannot be opened?
  → Solver prints an error and exits with return code 1.

---

## Requirements

### Functional Requirements

- **FR-001**: The solver MUST read a guide file (`-g` / `--guidefile`) that
  specifies the number of subproblems, subproblem LP filenames, the master LP
  filename, and optionally a monolithic LP filename and objective shift.
- **FR-002**: LP files MUST be in CPLEX LP format (read via `lpx_read_cpxlp`).
- **FR-003**: The solver MUST spawn one POSIX thread per subproblem; subproblems
  MUST be solved in parallel.
- **FR-004**: The master MUST coordinate subproblem threads via a semaphore-based
  service queue; newly generated columns MUST be offered to the master in the
  order subproblems complete each iteration.
- **FR-005**: The solver MUST implement a Phase I (Big-M / auxiliary variable
  approach) to find an initial feasible solution to the master when equality
  constraints are present.
- **FR-006**: The solver MUST implement Phase II (standard DW column generation)
  up to `MAX_PHASE2_ITERATIONS` (default: 3000) or until the reduced cost
  criterion is satisfied within `TOLERANCE` (1e-6).
- **FR-007**: The solver MUST write the LP relaxation solution to a file named
  `relaxed_solution` in the current working directory.
- **FR-008**: GLPK terminal output MUST be suppressed/redirected to a file
  named `out_terminal` in the current working directory.
- **FR-009**: The solver MUST free all allocated memory and destroy all
  synchronization primitives before exiting.
- **FR-010**: The solver MUST support verbosity levels: `--quiet` (`-q`),
  `--normal-output` (`-n`), `--verbose`, and `--output-all`.

### Optional / Configurable Behaviors (currently implemented)

- **FR-011**: `--integerize` (`-i`): After LP relaxation, attempt to recover
  an integer solution by re-solving subproblems with integrality constraints
  and a heuristic rounding step.
- **FR-012**: `--round` (`-r`): Enable solution rounding heuristic.
- **FR-013**: `--perturb` (`-p`): Enable perturbation of objective coefficients
  (anti-degeneracy).
- **FR-014**: `--mip-gap <value>`: Set the MIP optimality gap tolerance
  (default: 0.01).
- **FR-015**: `--write-bases`: Write GLPK basis files at each iteration.
- **FR-016**: `--write-int-probs`: Write intermediate MIP LP files.
- **FR-017**: `--no-write-final-master`: Suppress writing the final master LP.
- **FR-018**: `--print-timing`: Print wall-clock and CPU timing data on exit.
- **FR-019**: `--phase1_max <N>`: Override the Phase I iteration limit
  (default: 100).
- **FR-020**: `--phase2_max <N>`: Override the Phase II iteration limit
  (default: 3000).
- **FR-021**: `--skip-monolithic-read`: Skip reading the monolithic file entry
  from the guide file.
- **FR-022**: `--version` (`-v`): Print version string and exit.
- **FR-023**: `--help` (`-h`): Print usage and exit.

### Key Data Structures

- **`faux_globals`**: Central parameter/state struct passed to all functions.
  Contains verbosity, iteration limits, flags, service queue pointers, and
  per-subproblem solution arrays. Avoids true C global variables for solver
  configuration.
- **`D_matrix`**: Compressed sparse row representation of the coupling
  constraint coefficient matrix, built from the master LP file and shared
  read-only across all subproblem threads.
- **`master_data`** (`md`): Dual variables (`row_duals`) for coupling
  constraints; updated each DW iteration by the master and read by subproblems
  to re-price their objective functions.
- **`subprob_struct`**: Per-thread data structure containing the subproblem
  LP (`glp_prob*`), solution vectors, column translation lookup table, reduced
  cost (`r`), objective value (`obj`), and command/signal state.
- **`signal_data`** (`signals`): Iteration counters and ready flags coordinating
  master/subproblem synchronization.

### Synchronization Primitives (globally declared in `dw.h`)

| Primitive | Purpose |
|---|---|
| `master_lp_ready_mutex` / `master_lp_ready_cv` | Subproblems wait until master is set up |
| `next_iteration_mutex` / `next_iteration_cv` | Master signals start of next DW iteration |
| `service_queue_mutex` | Protects the circular service queue (head/tail indices) |
| `master_mutex` | Serializes subproblem access during column translation setup |
| `reduced_cost_mutex` | Protects reduced cost comparisons |
| `glpk_mutex` | Serializes calls to non-reentrant GLPK routines |
| `fputs_mutex` | Serializes output writes across threads |
| `sub_data_mutex[i]` | Per-subproblem data protection |
| `customers` (semaphore) | Counts subproblems ready to be serviced; named or unnamed |

### Source File Map

| File | Responsibility |
|---|---|
| `src/dw_main.c` | Entry point, master thread: setup, phase I/II loop, cleanup |
| `src/dw_phases.c` | `phase_1_iteration()` and `phase_2_iteration()` — master-side DW logic |
| `src/dw_subprob.c` | `subproblem_thread()` — per-thread LP solve, column pricing, signalling |
| `src/dw_support.c` | CLI parsing, initialization, matrix prep, utility functions, `dw_printf()` |
| `src/dw_rounding.c` | Solution rounding heuristic and post-relaxation integer recovery |
| `src/dw_blas.c` | Portable dense/sparse BLAS routines (`dw_daxpy`, `dw_ddot`, `dw_dgemv`, etc.) |
| `src/dw.h` | Shared type definitions, constants, globally-declared synchronization variables |
| `src/dw_subprob.h` | `subprob_struct` definition and subproblem thread interface |
| `src/dw_support.h` | Support function declarations |
| `src/dw_rounding.h` | Rounding function declarations |
| `src/dw_phases.h` | Phase function declarations |
| `src/dw_blas.h` | BLAS function declarations |
| `src/glpk.h` + `src/glp*.c` | Embedded GLPK 4.44 (thread-patched) — not authored by project |

### Build System

- **Build tool**: GNU Autotools (`configure`, `Makefile.am` / `Makefile.in`)
- **Version**: 1.2.1 (declared in `configure.ac`)
- **Configure flags of note**:
  - `--enable-named-semaphores`: Use POSIX named semaphores (`sem_open`/`sem_close`). Required on macOS. Controlled by `USE_NAMED_SEMAPHORES` preprocessor define.
  - `--enable-recursive-mutex`: Use recursive mutexes. May be needed on macOS.
  - `--with-gmp`: Link GNU MP bignum library (via GLPK).
  - `--with-zlib`: Link zlib (via GLPK).
  - `--enable-dl`: Enable shared library / dynamic loading support.
  - `USE_INTEL_MKL`: Optional compile-time flag to substitute Intel MKL sparse BLAS for the built-in `dw_blas` routines. Not currently set by the build system by default.

### Input File Format — Guide File

```
<number_of_subproblems>
<subproblem_1_filename.cplex>
[<subproblem_2_filename.cplex>]
...
<master_filename.cplex>
<monolithic_filename.cplex>        # vestigial; read but not used
[<objective_shift_value>]          # optional floating-point constant
```

All LP files are in CPLEX LP format as consumed by `lpx_read_cpxlp` (GLPK).

### Output Files (written to current working directory)

| File | Contents |
|---|---|
| `relaxed_solution` | Tab-separated variable name / value pairs for the LP relaxation optimum |
| `out_terminal` | Redirected GLPK terminal output (simplex iteration logs) |
| `final_master.cplex` | Master LP at termination (unless `--no-write-final-master`) |
| `basis_<N>` | GLPK basis files per iteration (only with `--write-bases`) |

### Known Defects and Limitations (as of baseline)

- **KD-001 (Critical)**: Multiple definition of global variables in `dw.h` —
  all synchronization primitives and global GLPK objects are *defined* (not
  merely declared) in a header file included by multiple translation units.
  This causes linker errors (`multiple definition of ...`) on Linux with GCC's
  default settings. macOS/Clang is more permissive and links successfully. This
  is the primary known build defect.
- **KD-002**: Unbounded subproblems are not supported. The code detects
  `GLP_UNBND` but prints a warning and does not handle the case correctly.
- **KD-003**: The monolithic file entry in the guide file is vestigial — it is
  read from the guide file but never opened or used. Old guide files must still
  include the line to avoid parse errors unless `--skip-monolithic-read` is passed.
- **KD-004**: `sprintf` is used in some places instead of `snprintf`; buffer
  overflow check is not universally applied. Identified in original author's
  comments.
- **KD-005**: `malloc` usage is inconsistent. Some allocations use
  `sizeof(char)*BUFF_SIZE` patterns that could be simplified; return values
  are not always checked.
- **KD-006**: The service queue is a fixed-size circular array of size
  `num_clients`. With thousands of threads this circular array is re-used;
  correctness depends on threads entering the queue at most once per iteration.
- **KD-007**: Phase II has a static stack-allocated variable `iteration_count`
  in `phase_1_iteration` (shared state via `static` — not thread-safe if ever
  called from multiple threads, though in current design it is only called from
  the master thread).
- **KD-008**: The `purge_nonbasics()` function is implemented but not called
  anywhere. The author noted it increased runtime and caused cycling.
- **KD-009**: The `four_sea` and `book_dantzig` examples have non-deterministic
  variable assignments. `tests/dw-tests.sh` does not cover these two examples.

---

## Test Coverage (as of baseline)

`tests/dw-tests.sh` runs the following examples and diffs sorted output:

| Example | Covered | Deterministic |
|---|---|---|
| `book_bertsimas` | ✅ | ✅ |
| `book_lasdon` | ✅ | ✅ |
| `web_mitchell` | ✅ | ✅ |
| `web_trick` | ✅ | ✅ |
| `four_sea` | ❌ | ❌ (objective deterministic, assignment not) |
| `book_dantzig` | ❌ | ❌ (objective deterministic, assignment not) |

The test script assumes `dwsolver` is on `$PATH`. It does not build the binary.

---

## Platform Status (as of baseline)

| Platform | Build | Tests pass |
|---|---|---|
| macOS (Clang) | ✅ (via GH Actions) | ✅ |
| Linux (GCC) | ❌ (linker: multiple definition in `dw.h`) | ❌ |
| Windows | Unknown — not attempted | Unknown |
