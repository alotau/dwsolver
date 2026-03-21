# Feature Specification: SEI CERT C Standard Compliance

**Feature Branch**: `008-sei-cert-c-compliance`
**Created**: 2026-03-21
**Status**: Draft
**Input**: Apply the SEI CERT C Coding Standard to dwsolver's authored source files (`src/dw_*.c`, `src/dw_*.h`) to eliminate POSIX threading hazards, unchecked error paths, and integer-type issues discovered during the architecture review of feature 007.

## Background

The SEI CERT C Coding Standard (https://wiki.sei.cmu.edu/confluence/display/c) defines rules and recommendations for safe, secure, and reliable C code. A review of the dwsolver source identified several categories of non-conformance after prior safety work (features 002 and 005) was completed:

- **Already fixed (out of scope)**: `STR07-C`/`STR31-C` (`sprintf`→`snprintf`), `MEM32-C` (7 critical `malloc` null-checks), `MEM` leak in `dw_rounding.c`, `POS53-C` (`while` guard on `pthread_cond_wait`).
- **Remaining**: POSIX pthread error-return checking (`POS54-C`), blocking operations under mutex (`POS52-C`), lock acquisition order (`POS51-C`), remaining non-critical `malloc` null-checks (`MEM32-C`), file handle error paths (`FIO01-C`, `FIO42-C`), and signed/unsigned integer type issues (`INT30-C`, `INT31-C`).

**Scope**: `src/dw_*.c` and `src/dw_*.h` only. Vendored GLPK (`src/glp*.c`) is explicitly excluded — it is treated as a third-party library boundary.

## User Scenarios & Testing

### User Story 1 — POSIX threading primitives are error-checked (Priority: P1)

A developer modifying the concurrency code can rely on all `pthread_*` and `sem_*` calls failing loudly and immediately rather than silently returning non-zero error codes that could corrupt program state. Today, none of the 15+ pthread/semaphore calls in `dw_main.c`, `dw_phases.c`, `dw_subprob.c`, and `dw_support.c` check their return values.

**Why this priority**: Threading errors that go unchecked (e.g., a failed `pthread_mutex_lock`) cause silent data races or deadlocks that are extremely difficult to diagnose. This is the highest-risk gap in the current code — CERT rule POS54-C.

**Independent Test**: Running `grep -n 'pthread_\|sem_' src/dw_*.c | grep -v '\/\/' | grep -v '= 0'` returns zero unchecked call sites; the full test suite passes without regression.

**Acceptance Scenarios**:

1. **Given** the modified source, **when** `grep` is run for unchecked `pthread_create`, `pthread_mutex_lock`, `pthread_mutex_unlock`, `pthread_cond_wait`, `pthread_cond_broadcast`, `pthread_join`, `sem_init`, `sem_wait`, and `sem_post` calls in `src/dw_*.c`, **then** every call site either ignores a well-documented always-succeeds operation or assigns and tests the return value.
2. **Given** a threading error at a checked call site, **when** a pthread primitive returns non-zero, **then** the program terminates with a clear diagnostic message to `stderr` identifying the failed primitive and the error code.
3. **Given** the modified code, **when** the full test suite is run, **then** all tests pass with no regression.

---

### User Story 2 — Blocking operations are not called while holding a mutex (Priority: P2)

A developer reviewing the master thread's Phase 1 and Phase 2 iteration loops can confirm that no long-running operations (specifically GLPK LP solves) are called while holding a mutex, which would prevent subproblem threads from making progress and create latent priority-inversion risk.

**Why this priority**: CERT rule POS52-C. Holding a lock across a blocking or long-running call defeats the purpose of the mutex and can cause deadlock if any other thread tries to acquire the same lock. GLPK simplex solves can take seconds on large problems.

**Independent Test**: An audit document in `specs/008-sei-cert-c-compliance/` lists every `glp_simplex` / `glp_intopt` call site in `src/dw_*.c`, its surrounding lock state, and a pass/fail verdict. All verdicts are "pass" (no lock held at call time).

**Acceptance Scenarios**:

1. **Given** the source code, **when** every `glp_simplex(` and `glp_intopt(` call site in `src/dw_*.c` is examined, **then** none are made while any mutex acquired in the same scope is still held.
2. **Given** any code path that reaches a GLPK solve, **when** the sequence of lock acquisitions and releases is traced, **then** all locks acquired before the solve have been released before the call is made.

---

### User Story 3 — Lock acquisition order is documented and enforced (Priority: P2)

A developer can consult a single authoritative lock-order table that defines the permitted acquisition sequence for all five synchronization primitives, and the source code conforms to that order with no inversions.

**Why this priority**: CERT rule POS51-C. Lock-order inversions are the canonical cause of deadlock in multi-threaded C programs. With five primitives (`master_mutex`, `sub_data_mutex[i]`, `next_iteration_mutex`, `service_queue_mutex`, `master_lp_ready_mutex`), the potential inversion pairs require explicit documentation and verification.

**Independent Test**: The lock-order audit document lists all five primitives in their required acquisition order; a manual audit of every multi-lock acquisition in `src/dw_*.c` shows zero inversions; ThreadSanitizer (`-fsanitize=thread`) build reports no lock-order violations when the test suite runs.

**Acceptance Scenarios**:

1. **Given** the specs directory, **when** the lock-order audit document is read, **then** a total order over all five synchronization primitives is defined with rationale for each ordering decision derived from actual acquisition patterns in the source.
2. **Given** the source code, **when** every code path that acquires more than one primitive is examined, **then** no path acquires locks in an order that violates the table.
3. **Given** any inversion found and corrected, **when** the test suite is run with ThreadSanitizer enabled, **then** no data-race or lock-order warnings are emitted.

---

### User Story 4 — All malloc sites and file handles are error-checked (Priority: P3)

A developer can confirm that all `malloc`/`calloc` calls in `src/dw_*.c` null-check their return values (completing the partial coverage from feature 002), and that every `fopen` call null-checks its return and every opened file handle is closed on all exit paths.

**Why this priority**: CERT rules `MEM32-C`, `FIO01-C`, `FIO42-C`. Lower priority than threading issues because these affect only error paths (OOM/missing file), not correctness under normal operation.

**Independent Test**: `grep -n 'malloc\|calloc\|fopen' src/dw_*.c` returns no sites without an immediately following null-check; manual review confirms all `fopen` handles have a matching `fclose` on all exit paths.

**Acceptance Scenarios**:

1. **Given** the modified source, **when** all `malloc` and `calloc` call sites in `src/dw_*.c` are reviewed, **then** each has a null-check that terminates or returns an error before dereferencing the pointer.
2. **Given** the modified source, **when** all `fopen` call sites in `src/dw_*.c` are reviewed, **then** each checks for a `NULL` return and every non-NULL handle has a corresponding `fclose` reachable on all code paths.
3. **Given** the modified code, **when** the full test suite is run, **then** all tests pass with no regression.

---

### User Story 5 — Integer type hazards are eliminated (Priority: P3)

Loop counters, array indices, and size computations in `src/dw_*.c` use types appropriate to their semantics, eliminating signed-overflow risk and signed/unsigned comparison warnings.

**Why this priority**: CERT rules `INT30-C`, `INT31-C`. Lower risk in practice given current problem sizes, but they contribute to a clean `-Wall -Wextra` build and prevent subtle bugs at larger scale.

**Independent Test**: `gcc -Wall -Wextra -Wsign-compare src/dw_*.c` produces zero warnings related to signed/unsigned comparison or integer conversion; the full test suite passes.

**Acceptance Scenarios**:

1. **Given** the modified source, **when** compiled with `-Wall -Wextra -Wsign-compare`, **then** zero warnings about signed/unsigned comparison or implicit integer conversion are emitted for `src/dw_*.c`.
2. **Given** any changed loop or computation, **when** the type change is reviewed, **then** the new type is appropriate to the semantic range of the value.
3. **Given** the modified code, **when** the full test suite is run, **then** all tests pass with no regression.

---

### Edge Cases

- What if a `pthread_mutex_unlock` call fails? POSIX specifies `pthread_mutex_unlock` on a mutex the calling thread owns cannot fail in all defined implementations; documenting these as "always-succeeds" with a comment satisfies the CERT rule without requiring an explicit check.
- What if ThreadSanitizer reports a race on a path that appears correct? TSan can produce false positives with complex condition variable patterns; each report must be triaged individually rather than suppressed wholesale.
- What if fixing a lock-order inversion requires restructuring acquisition order across both master and subproblem threads? Any such change must be covered by an existing test case and reviewed for semantic correctness, not just tool compliance.
- What if an `INT30-C` fix cascades to GLPK API calls that expect `int` parameters? At the GLPK boundary, explicit checked casts with a range assertion are the correct resolution — propagating `size_t` into the GLPK API is not appropriate.

---

## Requirements

### Functional Requirements

- **FR-001**: Every fallible `pthread_*` and `sem_*` call site in `src/dw_*.c` MUST either assign and test the return value, or be annotated with a comment explaining why the return is unconditionally safe to ignore (POS54-C).
- **FR-002**: No call to `glp_simplex` or `glp_intopt` in `src/dw_*.c` MUST be made while any of the five synchronization primitives is held by the calling thread (POS52-C).
- **FR-003**: A total lock acquisition order MUST be defined over all five synchronization primitives (`master_mutex`, `sub_data_mutex[i]`, `next_iteration_mutex`, `service_queue_mutex`, `master_lp_ready_mutex`) and recorded in `specs/008-sei-cert-c-compliance/lock-order.md` (POS51-C).
- **FR-004**: Every code path in `src/dw_*.c` that acquires more than one synchronization primitive MUST acquire them in the order defined by FR-003, with zero inversions (POS51-C).
- **FR-005**: Every `malloc` and `calloc` call in `src/dw_*.c` MUST have a null-check before the returned pointer is first dereferenced (MEM32-C).
- **FR-006**: Every `fopen` call in `src/dw_*.c` MUST null-check its return value; every successfully opened file handle MUST have a corresponding `fclose` reachable on all execution paths within the enclosing scope (FIO01-C, FIO42-C).
- **FR-007**: Signed/unsigned comparison warnings and implicit integer conversion warnings MUST be eliminated from `src/dw_*.c` when compiled with `-Wall -Wextra -Wsign-compare` (INT30-C, INT31-C).
- **FR-008**: A static analysis baseline report MUST be produced by running `cppcheck --enable=all` or `clang --analyze` against `src/dw_*.c` before implementation begins; the same tool MUST report no new findings after all changes are complete.
- **FR-009**: The full test suite MUST pass without regression after all changes.
- **FR-010**: All changes MUST be confined to `src/dw_*.c` and `src/dw_*.h`; vendored GLPK files (`src/glp*.c`, `src/glp*.h`) MUST NOT be modified.

### Key Entities

- **pthread error-check wrapper**: A consistent macro or inline function used to check the return value of every fallible pthread primitive — wraps the call, tests for non-zero, and aborts with a diagnostic message. Analogous to the existing `dw_oom_abort` helper introduced in feature 002.
- **Lock-order table**: A documented total order over the five synchronization primitives stored as `specs/008-sei-cert-c-compliance/lock-order.md`. This is an audit artifact, not source code.
- **Static analysis baseline**: The output of a static analyzer run before any changes begin, used as a before/after comparison to demonstrate net reduction in findings.

---

## Success Criteria

### Measurable Outcomes

- **SC-001**: Zero unchecked `pthread_*` / `sem_*` return values remain in `src/dw_*.c` — verified by grep and code review.
- **SC-002**: Zero `glp_simplex` or `glp_intopt` call sites in `src/dw_*.c` are made while any mutex is held — verified by the lock-state audit document.
- **SC-003**: A complete lock-order table covers all five synchronization primitives and a manual audit of every multi-lock acquisition path in `src/dw_*.c` shows zero inversions — verified by the lock-order audit document.
- **SC-004**: Zero `malloc`/`calloc` or `fopen` return values in `src/dw_*.c` go unchecked — verified by grep and code review.
- **SC-005**: Compilation of `src/dw_*.c` with `-Wall -Wextra -Wsign-compare` produces zero signed/unsigned or integer-conversion warnings — verified by CI build log.
- **SC-006**: The static analysis tool reports no findings in `src/dw_*.c` that were absent in the pre-implementation baseline — verified by diffing before/after reports.
- **SC-007**: The full test suite passes on the CI pipeline after all changes — verified by a green CI run.

---

## Assumptions

- ThreadSanitizer (`-fsanitize=thread`) is available in the CI environment for lock-order verification (GCC ≥ 4.8 or Clang ≥ 3.2, Linux only).
- `cppcheck` or `clang --analyze` is available in the development environment; if neither is present, installation is a prerequisite task, not a blocker for this spec.
- The GLPK simplex and intopt routines are not re-entrant; this is a known constraint that this feature does not change.
- The five-primitive lock set is complete — no additional hidden locks exist in the codebase (verifiable by `grep pthread_mutex_init src/dw_*.c src/dw_*.h`).
- Fixing `INT30-C`/`INT31-C` issues at the GLPK API boundary uses explicit checked casts with range assertions rather than propagating `size_t` into the GLPK API.
