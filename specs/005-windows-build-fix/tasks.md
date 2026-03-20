# Tasks: Code Quality Improvements — Bug Fixes, Duplication Reduction, and Clarity

**Input**: Design documents from `specs/005-windows-build-fix/`
**Branch**: `chore/code-quality` (create from master after `005-windows-build-fix` merges)
**Prerequisites**: plan.md ✅, spec.md ✅, research.md ✅, data-model.md ✅, quickstart.md ✅

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies between these tasks)
- **[Story]**: Which user story this task belongs to (US1–US6)
- Exact file paths are included in every description

---

## Phase 1: Setup

**Purpose**: Create the working branch and confirm the baseline is green before any code changes.

- [ ] T001 Create branch `chore/code-quality` from master after `005-windows-build-fix` merges
- [ ] T002 Run `bash tests/dw-tests.sh` on the new branch and confirm all 13 tests pass; record baseline warning count (`make clean && make 2>&1 | grep -c 'warning:'`)

---

## Phase 2: Foundational (Baseline Verification)

**Purpose**: Reproduce the known bugs locally so each fix can be verified against a concrete failure mode.

**⚠️ CRITICAL**: Confirm each bug is reproducible before fixing it.

- [ ] T003 Build with ASan (`CFLAGS="-fsanitize=address -g" ./configure && make`) and run any example to confirm the `val=NULL` UB is caught in `src/dw_phases.c` `phase_2_iteration`
- [ ] T004 Build with LSan or run under Valgrind (`valgrind ./dwsolver -r examples/example1.lp`) to confirm the `local_col_name` malloc leak appears under `src/dw_rounding.c:rounding_thread`

**Checkpoint**: Both bugs are reproducible — story-by-story fixes can now begin

---

## Phase 3: User Story 1 — Fix Undefined Behaviour in `phase_2_iteration` (Priority: P1) 🎯 MVP

**Goal**: Eliminate the `val=NULL` UB in the Phase II diagnostic block so the solver is safe under ASan, Valgrind, and any strict GLPK build.

**Independent Test**: Build with `-fsanitize=address`; run all 13 test cases; confirm zero ASan errors in `phase_2_iteration`.

### Implementation for User Story 1

- [ ] T005 [US1] In `src/dw_phases.c` `phase_2_iteration`, locate the diagnostic block that calls `glp_get_mat_col(master_lp, col, ind, val)` with `val=NULL`; replace with `ind2`/`val2` allocations sized `glp_get_num_rows(master_lp)+1` as specified in `specs/005-windows-build-fix/data-model.md`
- [ ] T006 [US1] In the same block in `src/dw_phases.c`, add `free(ind2); free(val2);` after the diagnostic loop so neither allocation leaks when the block exits
- [ ] T007 [US1] Rebuild with `-fsanitize=address` and run `bash tests/dw-tests.sh` to confirm all 13 tests pass with no ASan errors; rebuild without sanitiser flags and re-run to confirm clean release build
- [ ] T008 [US1] Commit: `fix: allocate ind2/val2 in phase_2_iteration diagnostic block (was UB)`

**Checkpoint**: User Story 1 complete — Phase II UB is eliminated and test suite is green.

---

## Phase 4: User Story 2 — Fix Memory Leak in `rounding_thread` (Priority: P1)

**Goal**: Replace the `const char*` malloc with a stack buffer so no per-subproblem allocation leaks during rounding passes.

**Independent Test**: Run `dwsolver -r examples/example1.lp` under Valgrind (`--leak-check=full`); confirm zero leaks attributed to `rounding_thread`.

### Implementation for User Story 2

- [ ] T009 [US2] In `src/dw_rounding.c` `rounding_thread`, change `const char* local_col_name = malloc(sizeof(char)*BUFF_SIZE);` to `char local_col_name[BUFF_SIZE];` as specified in `specs/005-windows-build-fix/data-model.md`
- [ ] T010 [US2] In `src/dw_rounding.c`, replace the direct pointer assignment `local_col_name = glp_get_col_name(...)` with `strncpy(local_col_name, glp_get_col_name(...), BUFF_SIZE-1); local_col_name[BUFF_SIZE-1] = '\0';`
- [ ] T011 [US2] Remove any `free(local_col_name)` call in the same function that was paired with the old `malloc` (to avoid a double-use-after-stack bug)
- [ ] T012 [US2] Run Valgrind (`valgrind --leak-check=full ./dwsolver -r examples/example1.lp`) and confirm zero leaks in `rounding_thread`; run `bash tests/dw-tests.sh` to confirm all 13 tests pass
- [ ] T013 [US2] Commit: `fix: use stack buffer for local_col_name in rounding_thread (was malloc leak)`

**Checkpoint**: User Story 2 complete — rounding memory leak is gone; both P1 bugs are fixed.

---

## Phase 5: User Story 3 — Eliminate `phase_1`/`phase_2` Code Duplication (Priority: P2)

**Goal**: Replace the near-identical `phase_1_iteration` and `phase_2_iteration` functions with a single `dw_iteration` function that accepts `phase` and `first_run` parameters.

**Independent Test**: Run `bash tests/dw-tests.sh` after each sub-step; all 13 tests must pass with identical stdout as the pre-refactoring baseline.

### Implementation for User Story 3

- [ ] T014 [US3] In `src/dw_phases.h`, add the `dw_iteration` declaration as specified in `specs/005-windows-build-fix/data-model.md` (parameters: `subprob_struct*`, `faux_globals*`, `master_data*`, `int phase`, `int first_run`, `char**`, `double*`, `int*`)
- [ ] T015 [US3] In `src/dw_phases.c`, implement `dw_iteration` by diff-merging the two existing functions: shared body in one place, Phase I–only branches (`first_run`, `obj_names`/`obj_coefs`/`obj_count`, y-accumulator setup, sign correction) guarded by `if (phase == 1 && first_run)`
- [ ] T016 [US3] In `src/dw_phases.c`, replace `phase_1_iteration` and `phase_2_iteration` bodies with thin wrappers that delegate to `dw_iteration` (keeps the diff reviewable without breaking any external call site)
- [ ] T017 [US3] In `src/dw_main.c`, update call sites to call `dw_iteration` directly instead of the wrapper names
- [ ] T018 [US3] Remove the wrapper functions from `src/dw_phases.c` and their declarations from `src/dw_phases.h` now that no call sites remain
- [ ] T019 [US3] Run `bash tests/dw-tests.sh` and confirm all 13 tests pass; `grep -n 'phase_1_iteration\|phase_2_iteration' src/*.c src/*.h` must return no matches
- [ ] T020 [US3] Commit: `refactor: unify phase_1_iteration/phase_2_iteration into dw_iteration`

**Checkpoint**: User Story 3 complete — single iteration function; the TOLERANCE divergence between phases is resolved.

---

## Phase 6: User Story 4 — Remove Dead Code and Standardise Output Logging (Priority: P3)

**Goal**: Delete commented-out code blocks; route all solver `printf` calls through `dw_printf` so `--silent` suppresses them.

**Independent Test**: `grep -n '//\s*printf\|//\s*glp_' src/dw_phases.c src/dw_rounding.c` returns no substantive commented-out blocks; run `./dwsolver --silent examples/example1.lp` and confirm no solver diagnostic output on stdout.

### Implementation for User Story 4

- [ ] T021 [P] [US4] In `src/dw_phases.c`, delete all commented-out blocks of more than 3 consecutive lines (verify via `grep -n` before and after); no logic may be deleted — only already-commented code
- [ ] T022 [P] [US4] In `src/dw_rounding.c`, delete dead commented-out blocks including the `print_zeros` dead function and any commented-out diagnostic `printf` calls
- [ ] T023 [US4] In `src/dw_support.c`, identify every `printf(` call that is not the CLI usage/help output; replace with the appropriate `dw_printf(globals, ...)` or `fprintf(stderr, ...)` call
- [ ] T024 [US4] Run `bash tests/dw-tests.sh` to ensure all 13 tests pass; run with `--silent` flag and confirm only GLPK file output is written (no unexpected stdout)
- [ ] T025 [US4] Commit: `refactor: remove dead code; route printf through dw_printf in dw_support.c`

**Checkpoint**: User Story 4 complete — clean files, suppressible output.

---

## Phase 7: User Story 5 — Reduce Lock Contention on `glpk_mutex` and `master_mutex` (Priority: P3)

**Goal**: Guard the `glpk_mutex` logging block with a verbosity check; replace `master_mutex` in the subproblem read-only loop with a `pthread_rwlock_t`.

**Independent Test**: Confirm `glpk_mutex` is not acquired when `verbosity < OUTPUT_ALL`; add `pthread_rwlock_t master_data_rwlock` and confirm `bash tests/dw-tests.sh` passes; optionally verify with Helgrind that no new lock-order violations appear.

### Implementation for User Story 5

- [ ] T026 [US5] In `src/dw_subprob.c` `prepare_column`, wrap the `pthread_mutex_lock(&glpk_mutex)` / `pthread_mutex_unlock(&glpk_mutex)` logging block with `if (globals->verbose >= OUTPUT_ALL) { ... }` so the lock is never acquired at lower verbosity levels
- [ ] T027 [US5] In `src/dw.h`, add `extern pthread_rwlock_t master_data_rwlock;` alongside the existing mutex externs; in `src/dw_globals.c`, add `pthread_rwlock_t master_data_rwlock = PTHREAD_RWLOCK_INITIALIZER;`
- [ ] T028 [US5] In `src/dw_support.c` `init_pthread_data`, add `pthread_rwlock_init(&master_data_rwlock, NULL);`; in the corresponding teardown / `free_globals`, add `pthread_rwlock_destroy(&master_data_rwlock);`
- [ ] T029 [US5] In `src/dw_subprob.c` `subproblem_thread`, replace `pthread_mutex_lock(&master_mutex)` / `pthread_mutex_unlock(&master_mutex)` around the column-mapping read loop with `pthread_rwlock_rdlock(&master_data_rwlock)` / `pthread_rwlock_rdunlock(&master_data_rwlock)`
- [ ] T030 [US5] Run `bash tests/dw-tests.sh` to confirm all 13 tests pass; optionally run TSan build (`CFLAGS="-fsanitize=thread -g" ./configure && make && bash tests/dw-tests.sh`) to confirm no new data race reports
- [ ] T031 [US5] Commit: `perf: guard glpk_mutex logging with verbosity check; use rwlock for original_master_lp reads`

**Checkpoint**: User Story 5 complete — reduced lock contention at scale.

---

## Phase 8: User Story 6 — Move Core Solver State Out of True Globals (Priority: P4)

**Goal**: Migrate `master_lp`, `original_master_lp`, `D`, `signals`, and the sync primitives into `faux_globals` / a new `sync_state` struct so the solver can be instantiated cleanly without global side effects.

**Independent Test**: All 13 tests pass; TSan reports no new races; `extern` count in `dw.h` is reduced to zero for the migrated symbols.

### Implementation for User Story 6

- [ ] T032 [US6] In `src/dw.h`, define `sync_state` struct containing all `pthread_mutex_t`, `pthread_cond_t`, `sem_t`, `pthread_rwlock_t`, and `pthread_attr_t` members currently declared as extern — follow the layout in `specs/005-windows-build-fix/data-model.md`
- [ ] T033 [US6] Add `sync_state sync` as an embedded (not pointer) member to `faux_globals` in `src/dw.h`; keep the extern declarations intact for now (compile-time only change at this step — no behaviour change)
- [ ] T034 [US6] Update `init_pthread_data` in `src/dw_support.c` to initialise all fields through `fg->sync.*` instead of the top-level externs; update `free_globals` similarly; confirm `bash tests/dw-tests.sh` still passes
- [ ] T035 [US6] Remove the now-redundant top-level extern mutex/cond/sem declarations from `src/dw.h` and their definitions from `src/dw_globals.c`; update all call sites that referenced `master_mutex`, `glpk_mutex`, etc., to use `fg->sync.master_mutex`, `fg->sync.glpk_mutex`, etc.
- [ ] T036 [US6] Add `D_matrix* D` to `faux_globals` in `src/dw.h`; update all call sites (search `extern D_matrix\* D` and every `D` reference in `src/`) to use `fg->D`; remove `extern D_matrix* D` from `src/dw.h` and its definition from `src/dw_globals.c`; run tests
- [ ] T037 [US6] Add `signal_data* signals` to `faux_globals` in `src/dw.h`; migrate all `signals` references to `fg->signals`; remove extern; run tests
- [ ] T038 [US6] Add `glp_prob* master_lp` and `glp_prob* original_master_lp` to `faux_globals` in `src/dw.h`; migrate all references across `src/dw_main.c`, `src/dw_phases.c`, `src/dw_subprob.c`, `src/dw_rounding.c`, `src/dw_support.c`; remove externs; run tests after each file is updated
- [ ] T039 [US6] Confirm `src/dw_globals.c` now contains only the `faux_globals` allocator / initialiser (no remaining extern definitions); if the file is empty, remove it and update `Makefile.am`
- [ ] T040 [US6] Run full TSan build and `bash tests/dw-tests.sh`; confirm zero new data race reports
- [ ] T041 [US6] Commit: `refactor: migrate master_lp, D, signals, sync primitives into faux_globals / sync_state`

**Checkpoint**: User Story 6 complete — no true globals remain for solver state.

---

## Final Phase: Polish & Cross-Cutting Concerns

- [ ] T042 [P] Verify warning count did not increase: `make clean && make 2>&1 | grep -c 'warning:'` must be ≤ baseline from T002
- [ ] T043 [P] Run `bash tests/dw-tests.sh` one final time on both macOS and Linux (or CI) to confirm all 13 tests pass with the full change set
- [ ] T044 Update `specs/005-windows-build-fix/spec.md` status from `Draft` to `Implemented`

---

## Dependencies

```
T001 → T002 → T003–T004 (parallel) → T005–T008 (US1) → T009–T013 (US2)
                                       ↓ (both P1 bugs fixed)
                                      T014–T020 (US3)
                                       ↓
                           T021–T025 (US4, P3) ←parallel→ T026–T031 (US5, P3)
                                       ↓
                                      T032–T041 (US6, P4)
                                       ↓
                                      T042–T044 (Polish)
```

**US4 and US5 are independent of each other** — they touch different files (`dw_phases.c`/`dw_rounding.c`/`dw_support.c` vs `dw_subprob.c`/`dw_globals.c`) and can be implemented in parallel branches if desired.

---

## Parallel Execution Opportunities

| Stories | Can run in parallel? | Notes |
|---------|---------------------|-------|
| US1 and US2 | ✅ Yes | Different files — `dw_phases.c` vs `dw_rounding.c` |
| US4 and US5 | ✅ Yes | Different primary files |
| T021 and T022 (within US4) | ✅ Yes | Different files |
| US3 and any P3 story | ❌ No | US3 restructures `dw_phases.c` — wait for it to merge first |
| US6 and anything | ❌ No | Touches every file — do last |

---

## Implementation Strategy

**MVP Scope (deliver first)**: US1 + US2 only — two targeted bug fixes, minimal diff, immediately verifiable via sanitisers. This is safe to PR and merge without waiting for the refactoring work.

**Increment 2**: US3 (phase unification) — medium-effort, builds confidence in the shared path, no behaviour change.

**Increment 3**: US4 + US5 in parallel — clarity and lock efficiency, low risk.

**Increment 4**: US6 (global state migration) — highest risk, should have its own dedicated PR with TSan verification.

---

## Summary

| Metric | Value |
|--------|-------|
| Total tasks | 44 |
| US1 (P1 UB fix) | 4 tasks |
| US2 (P1 leak fix) | 5 tasks |
| US3 (P2 duplication) | 7 tasks |
| US4 (P3 dead code / logging) | 5 tasks |
| US5 (P3 lock contention) | 6 tasks |
| US6 (P4 global state) | 10 tasks |
| Setup / Polish | 7 tasks |
| Parallelisable tasks [P] | 5 tasks |
| MVP scope | US1 + US2 (9 tasks) |
