# Tasks: SEI CERT C Standard Compliance

**Input**: Design documents from `specs/008-sei-cert-c-compliance/`
**Prerequisites**: plan.md ✅, spec.md ✅, research.md ✅, data-model.md ✅, contracts/ ✅, quickstart.md ✅

**Organization**: Tasks are grouped by user story. Each story phase is independently testable before the next begins.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no shared dependencies)
- **[Story]**: Which user story this task belongs to
- No test tasks generated (not requested in spec)

---

## Phase 1: Setup

**Purpose**: Add the `DW_PTHREAD_CHECK` error-check helper — the prerequisite macro for all P1 work. No other phase can begin until this is in place.

- [ ] T001 Add `dw_pthread_check()` inline function and `DW_PTHREAD_CHECK` macro to `src/dw_support.h` (see data-model.md pattern; analogous to existing `dw_oom_abort`)

**Checkpoint**: `make` succeeds; macro usable from all `dw_*.c` files.

---

## Phase 2: Foundational

**Purpose**: Establish the static-analysis baseline (FR-008) before making any changes, so regressions are detectable.

**⚠️ Complete before any user story implementation begins**

- [ ] T002 Run `cppcheck --enable=all src/dw_*.c` and capture output to `specs/008-sei-cert-c-compliance/baseline-cppcheck.txt`
- [ ] T003 [P] Compile with `-Wall -Wextra -Wsign-compare` and capture warnings to `specs/008-sei-cert-c-compliance/baseline-warnings.txt` (`gcc -Wall -Wextra -Wsign-compare -c src/dw_*.c -I src/ 2>&1 | grep -v 'glp' > baseline-warnings.txt`)
- [ ] T004 [P] Run full test suite and confirm all tests pass: `./tests/dw-tests.sh`

**Checkpoint**: Baseline captured; all tests green before any code changes.

---

## Phase 3: User Story 1 — POSIX threading primitives error-checked (Priority: P1) 🎯 MVP

**Goal**: Every fallible `pthread_*` and `sem_*` call site in `src/dw_*.c` is either wrapped in `DW_PTHREAD_CHECK` or annotated as always-succeeds. (FR-001, POS54-C)

**Independent Test**: `grep -n 'pthread_\|sem_' src/dw_*.c | grep -v '//' | grep -vE 'DW_PTHREAD_CHECK|always succeeds|pthread_t |pthread_attr_t|pthread_mutex_t|pthread_cond_t|pthread_mutexattr|sem_t'` returns zero unchecked call sites; `make` succeeds; `./tests/dw-tests.sh` passes.

- [ ] T005 [P] [US1] Apply `DW_PTHREAD_CHECK` to all unchecked `pthread_mutex_init`, `pthread_cond_init`, `pthread_attr_init`, `pthread_mutexattr_init`, `sem_init`, `pthread_attr_setdetachstate` calls in `src/dw_support.c` init_pthread_data() (lines ~190–212 per contracts/pthread-call-sites.md); test `rc = pthread_attr_setstacksize` result; annotate `pthread_attr_getstacksize` calls as always-succeeds
- [ ] T006 [P] [US1] Apply `DW_PTHREAD_CHECK` to all unchecked `pthread_mutex_lock` calls in `src/dw_main.c` (lines 221, 258, 458, 479, 541, 640, 644, 826); annotate all `pthread_mutex_unlock` and `pthread_cond_broadcast` calls as always-succeeds; verify `pthread_create` (line 204) and `pthread_join` (line 706) already checked
- [ ] T007 [P] [US1] Apply `DW_PTHREAD_CHECK` to all unchecked `pthread_mutex_lock` calls in `src/dw_phases.c` (lines 87, 97, 117, 271, 280, 334, 343, 360, 511, 521, 530); annotate all `pthread_mutex_unlock` and `pthread_cond_broadcast` as always-succeeds
- [ ] T008 [P] [US1] Apply `DW_PTHREAD_CHECK` to all unchecked `pthread_mutex_lock`, `sem_post`, and `pthread_cond_wait` call sites in `src/dw_subprob.c` (lines 101, 127, 181, 260, 368, 374/376, 381, 389, 435 per contracts/pthread-call-sites.md); annotate all unlocks and broadcasts as always-succeeds
- [ ] T009 [P] [US1] Apply `DW_PTHREAD_CHECK` to unchecked `pthread_mutex_lock` calls in `src/dw_rounding.c` (lines 518, 549, 591); annotate corresponding unlocks as always-succeeds; verify `pthread_create` (lines 175, 181) and `pthread_join` (line 199) already checked
- [ ] T010 [US1] Compile `src/dw_*.c` after all US1 edits (`make`) and fix any build errors introduced by the changes
- [ ] T011 [US1] Run `./tests/dw-tests.sh` and confirm all examples reproduce expected optima; fix any regressions before proceeding

**Checkpoint**: Zero unchecked pthread/sem call sites; all tests pass.

---

## Phase 4: User Story 2 — No blocking operations under mutex (Priority: P2)

**Goal**: `glp_simplex` and `glp_intopt` are never called while holding any synchronization primitive. Fix the 2 violations in `src/dw_subprob.c`. Fix the second spurious-wakeup bug (`if`→`while` at line 131). (FR-002, POS52-C, POS53-C)

**Independent Test**: Re-audit `contracts/glpk-lock-state.md` — all 10 GLPK call sites show "None" in the "Mutexes Held" column; `./tests/dw-tests.sh` passes; all examples reproduce expected optima.

- [ ] T012 [US2] In `src/dw_subprob.c`, fix the spurious-wakeup bug at line ~131: change `if (!signals->master_lp_ready)` to `while (!signals->master_lp_ready)` and add POS53-C comment (distinct from the line-382 fix already on `fix/spurious-wakeup-cond-wait`)
- [ ] T013 [US2] In `src/dw_subprob.c`, fix POS52-C violations at lines ~303 and ~306 (main iteration loop): release `sub_data_mutex[id]` after all `glp_set_obj_coef` calls complete; extract result values (`obj_val`, `result_unbounded`) to local variables; call `glp_simplex` / `glp_intopt` without the lock; re-acquire `sub_data_mutex[id]` before writing `my_data->obj`, `my_data->unbounded`, and solution arrays back (see contracts/glpk-lock-state.md fix specification)
- [ ] T014 [US2] Compile `src/dw_subprob.c` (`make`) and fix any build errors
- [ ] T015 [US2] Run `./tests/dw-tests.sh` and confirm all examples reproduce expected optima; verify the POS52-C fix does not change any output values

**Checkpoint**: GLPK solve sites clean; spurious-wakeup guards in place; all tests pass.

---

## Phase 5: User Story 3 — Lock acquisition order documented (Priority: P2)

**Goal**: A total order over all five synchronization primitives is defined in `contracts/lock-order.md`; source code conforms with zero inversions; audit document produced. (FR-003, FR-004, POS51-C)

**Independent Test**: `contracts/lock-order.md` exists and lists all 5 primitives in required acquisition order with per-site verdicts; manual scan of every multi-lock path in `src/dw_*.c` shows zero inversion rows marked ❌.

- [ ] T016 [US3] Review `contracts/lock-order.md` (already written in plan phase) against current source after US1 and US2 changes; update any line numbers that shifted due to edits in T005–T013; confirm all per-site verdicts remain ✅ PASS
- [ ] T017 [US3] Add a brief lock-order comment block at the top of `src/dw_support.c` (near `init_pthread_data`) documenting the six-level acquisition order (see data-model.md Lock Acquisition Total Order table) so the order is discoverable in source, not just in specs

**Checkpoint**: Lock-order documented in both specs and source; no inversions; no code correctness changes required (research confirmed zero inversions).

---

## Phase 6: User Story 4 — All malloc and fopen sites error-checked (Priority: P3)

**Goal**: Every `malloc`/`calloc` call in `src/dw_*.c` null-checks before first dereference; every `fopen` call has a null-check and a corresponding `fclose`. (FR-005, FR-006, MEM32-C, FIO01-C, FIO42-C)

**Independent Test**: `grep -n 'malloc\|calloc' src/dw_*.c` shows no site without an immediately following null-check; `grep -n 'fopen' src/dw_*.c` shows all sites with null-check; `make` succeeds; `./tests/dw-tests.sh` passes.

- [ ] T018 [P] [US4] Add `dw_oom_abort()` null-checks after every unchecked `malloc`/`calloc` in `src/dw_support.c` (priority sites: lines 99, 105–108, 136–137, 163–164, 173–174, 187, 193, 280, 378, 503–506, 510, 593–594, 609, 670, 709)
- [ ] T019 [P] [US4] Add `dw_oom_abort()` null-checks after every unchecked `malloc`/`calloc` in `src/dw_main.c` (priority sites: lines 110, 119, 121, 156–157, 159, 187, 201, 240, 242, 314, 338, 374, 417, 493, 495, 497, 547, 764)
- [ ] T020 [P] [US4] Add `dw_oom_abort()` null-checks after every unchecked `malloc`/`calloc` in `src/dw_phases.c` (lines 72–75, 114, 166, 322, 389–390, 395–396, 413)
- [ ] T021 [P] [US4] Add `dw_oom_abort()` null-checks after every unchecked `malloc`/`calloc` in `src/dw_subprob.c` (lines 79, 86, 104, 149–160, 173–174) and `src/dw_rounding.c` (lines 76, 78, 87–89, 92, 94, 132, 141, 273–274, 333–335, 380–382)
- [ ] T022 [US4] Verify all `fopen` call sites in `src/dw_*.c` have null-checks (already present per research); verify every successfully-opened file handle has a reachable `fclose` on all exit paths — add missing `fclose` calls if any are found
- [ ] T023 [US4] Compile (`make`) and run `./tests/dw-tests.sh`; fix any regressions

**Checkpoint**: All allocation sites null-checked; all file handles closed; tests pass.

---

## Phase 7: User Story 5 — Integer type hazards eliminated (Priority: P3)

**Goal**: Zero `-Wall -Wextra -Wsign-compare` warnings in `src/dw_*.c` (excluding vendored GLPK). (FR-007, INT30-C, INT31-C)

**Independent Test**: `gcc -Wall -Wextra -Wsign-compare -c src/dw_*.c -I src/ 2>&1 | grep -v 'glp'` produces zero lines of output; `./tests/dw-tests.sh` passes.

- [ ] T024 [P] [US5] In `src/dw_main.c`, fix all signed/unsigned comparison warnings: loop counters compared to `glp_get_num_cols`/`glp_get_num_rows` return (`int`); use explicit `(int)` casts at GLPK boundary call sites where `size_t` → `int` conversions occur
- [ ] T025 [P] [US5] In `src/dw_phases.c` and `src/dw_subprob.c`, fix all signed/unsigned comparison warnings using the same explicit-cast pattern
- [ ] T026 [P] [US5] In `src/dw_support.c` and `src/dw_rounding.c`, fix all signed/unsigned comparison warnings
- [ ] T027 [US5] Compile the full project (`make`) with `-Wall -Wextra -Wsign-compare` added to `CFLAGS` in `Makefile` (or `configure.ac`); confirm zero warnings remain from `dw_*.c` files; run `./tests/dw-tests.sh`

**Checkpoint**: Zero integer-type warnings; tests pass; baseline-warnings.txt new count = 0 for dw_*.c.

---

## Phase 8: Polish & Cross-Cutting Concerns

**Purpose**: Final verification that all CERT requirements are satisfied end-to-end.

- [ ] T028 [P] Re-run `cppcheck --enable=all src/dw_*.c` and confirm no new findings compared to `baseline-cppcheck.txt` (FR-008)
- [ ] T029 [P] Run ThreadSanitizer build and test suite per quickstart.md verification checklist: `./configure CFLAGS="-fsanitize=thread -g -O1" LDFLAGS="-fsanitize=thread" && make clean && make && ./tests/dw-tests.sh`; confirm zero data-race warnings
- [ ] T030 Run final acceptance grep checks from spec.md and quickstart.md: unchecked pthread grep, GLPK lock-state grep, sign-compare compile, full test suite — record results in `specs/008-sei-cert-c-compliance/acceptance-report.md`

**Checkpoint**: All FR-001 through FR-010 satisfied; TSan clean; acceptance report written.

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 (`DW_PTHREAD_CHECK` macro must exist before it can be used in any phase)
- **US1 (Phase 3)**: Depends on Phase 1 + Phase 2 (baseline captured before changes begin)
- **US2 (Phase 4)**: Depends on Phase 1 (needs `DW_PTHREAD_CHECK`); can proceed independently of US1 but benefits from US1 completing first to avoid merge conflicts in the same files
- **US3 (Phase 5)**: Depends on US1 + US2 completion (line numbers in audit may shift after earlier edits)
- **US4 (Phase 6)**: Depends on Phase 2 (baseline captured); independent of US1/US2/US3 (different change type, different lines); can proceed in parallel once baseline is done
- **US5 (Phase 7)**: Depends on Phase 2 (baseline captured); independent of other user stories; can proceed in parallel once baseline is done
- **Polish (Phase 8)**: Depends on all user stories complete

### User Story Dependencies

| Story | Depends On | Can Parallelize With |
|-------|-----------|---------------------|
| US1 (POS54-C) | Phase 1 + 2 | US4, US5 (different change type) |
| US2 (POS52-C) | Phase 1 | US4, US5 |
| US3 (POS51-C, docs) | US1 + US2 | US4, US5 |
| US4 (MEM32-C/FIO) | Phase 2 | US1, US2, US5 |
| US5 (INT30/31-C) | Phase 2 | US1, US2, US4 |

### Within Each User Story

- Within US1: Tasks T005–T009 can all be run in parallel (each targets a different file); T010 (build) and T011 (test) must follow
- Within US4: T018–T021 can all run in parallel (each targets different files); T022–T023 must follow

### Parallel Opportunities Per Story

**US1** — parallel file edits:
```
T005 (dw_support.c) ─┐
T006 (dw_main.c)    ──┤─→ T010 (build) → T011 (test)
T007 (dw_phases.c)  ──┤
T008 (dw_subprob.c) ──┤
T009 (dw_rounding.c)─┘
```

**US4** — parallel file edits:
```
T018 (dw_support.c)  ─┐
T019 (dw_main.c)     ──┤─→ T022 (fopen verify) → T023 (build + test)
T020 (dw_phases.c)   ──┤
T021 (dw_subprob.c + rounding) ─┘
```

**US5** — parallel file edits:
```
T024 (dw_main.c)    ─┐
T025 (phases+subprob)─┤─→ T027 (build with Wextra)
T026 (support+round) ─┘
```

---

## Implementation Strategy

**MVP scope**: US1 (Phase 3) alone. Completing POS54-C — adding `DW_PTHREAD_CHECK` to all ~24 call sites — delivers the highest-risk fix with the smallest blast radius. It is entirely additive (no behavior change, only error detection on previously-silent failures).

**Recommended delivery order**:
1. Phase 1 + Phase 2 (setup + baseline) — 30 minutes
2. US1 tasks T005–T011 (POS54-C, highest risk) — ~1 day
3. US2 tasks T012–T015 (POS52-C, requires restructuring dw_subprob.c iteration loop) — ~1 day
4. US3 tasks T016–T017 (POS51-C docs, no code changes) — ~30 minutes
5. US4 tasks T018–T023 (MEM32-C/FIO, mechanical null-checks) — ~1 day
6. US5 tasks T024–T027 (INT30/31-C warnings) — ~0.5–1 day
7. Phase 8 (TSan + final acceptance) — ~1 hour
