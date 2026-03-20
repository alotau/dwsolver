# Tasks: Cross-Platform Build Repair & Safety Hardening

**Input**: Design documents from `specs/002-cross-platform-repair/`
**Prerequisites**: plan.md ✅, spec.md ✅, research.md ✅, data-model.md ✅, quickstart.md ✅

**Tests**: No test-first tasks — this is a repair spec. Verification steps confirm each story's acceptance criteria.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story?] Description`

- **[P]**: Can run in parallel (different files, no incomplete dependencies)
- **[Story]**: Which user story this task belongs to (US1–US4)
- Exact file paths are included in each description

---

## Phase 1: Setup

**Purpose**: Confirm the pre-repair baseline so we have a known-good starting point.

- [X] T001 Verify baseline macOS build: `make clean && ./configure --enable-named-semaphores && make` in `/Users/joey/Development/dwsolver-repaired` — confirm zero errors and `src/dwsolver` produced

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: The KD-001 fix — introducing `src/dw_globals.c` and converting `src/dw.h` to `extern` declarations. This is a prerequisite for all user story verifications and all subsequent code changes.

**⚠️ CRITICAL**: No user story verification can run until this phase is complete.

- [X] T002 Create `src/dw_globals.c` — new file with `#include "dw.h"` and the single authoritative definition of all 18 globals: `attr`, `master_lp_ready_mutex`, `master_lp_ready_cv`, `service_queue_mutex`, `next_iteration_mutex`, `next_iteration_cv`, `master_mutex`, `reduced_cost_mutex`, `glpk_mutex`, `fputs_mutex`, `sub_data_mutex`, `customers` (with `#ifdef USE_NAMED_SEMAPHORES` pointer/value guard), `original_master_lp`, `master_lp`, `parm`, `simplex_control_params`, `D`, `signals`

- [X] T003 [P] Edit `src/dw.h` lines 109–167: convert all 18 variable definitions to `extern` declarations (prefix each with `extern`; preserve the `#ifdef USE_NAMED_SEMAPHORES` guard on `customers`). The `typedef` blocks for `D_matrix`, `signal_data`, `master_data`, `new_column`, `hook_struct`, `faux_globals`, `subprob_struct` are NOT globals and must remain unchanged.

- [X] T004 [P] Edit `src/Makefile.am` line 122: add `dw_globals.c \` to `dwsolver_SOURCES` immediately before `dw_main.c`

**Checkpoint**: Foundational changes complete — user story work can now begin.

---

## Phase 3: User Story 1 — Build and Pass Tests on Linux/GCC (Priority: P1) 🎯 MVP

**Goal**: Confirm the KD-001 fix compiles correctly on macOS (regression check) and documents the verification path for Linux.

**Independent Test**: `make clean && ./configure --enable-named-semaphores && make` completes with zero errors; `tests/dw-tests.sh` reports PASS for all 4 existing examples.

- [X] T005 [US1] Rebuild after KD-001 fix: from `/Users/joey/Development/dwsolver-repaired` run `make clean && ./configure --enable-named-semaphores && make` — confirm zero linker errors, `src/dwsolver` produced

- [X] T006 [US1] Run existing tests: `export PATH="$PWD/src:$PATH" && bash tests/dw-tests.sh` — confirm all 4 tests (`book_bertsimas`, `book_lasdon`, `web_mitchell`, `web_trick`) report PASS

**Checkpoint**: User Story 1 complete — macOS build and existing tests verified. Linux verification requires a GCC environment (Docker/CI); the fix is structurally correct per C99 §6.9 ¶5.

---

## Phase 4: User Story 2 — Test Suite Covers All Deterministic Examples (Priority: P2)

**Goal**: Extend `tests/dw-tests.sh` with two new objective-value tests covering `four_sea` and `book_dantzig`, the two examples currently lacking automated coverage.

**Independent Test**: Run `bash tests/dw-tests.sh` — six tests execute, including the two new ones which parse the final `#### Master objective value = ...` line from solver output and assert the expected value.

**Objective values confirmed by reference runs**:
- `four_sea`: final `#### Master objective value = 1.200000e+01`
- `book_dantzig`: final `#### Master objective value = 6.357895e+01`

- [X] T007 [P] [US2] Add `four_sea` objective-value test to `tests/dw-tests.sh`: `pushd` to `../examples/four_sea`, run `dwsolver -g guidefile > out_obj.txt`, check exit code is 0 (non-zero → FAIL immediately), extract the last `#### Master objective value = ...` line from `out_obj.txt` via `grep`, string-match it against `1.200000e+01` exactly, report PASS/FAIL, `popd`

- [X] T008 [P] [US2] Add `book_dantzig` objective-value test to `tests/dw-tests.sh`: same pattern as T007 — run `dwsolver -g guidefile > out_obj.txt`, check exit code, grep last `Master objective value` line, string-match against `6.357895e+01` exactly, report PASS/FAIL

- [X] T009 [US2] Run full extended test suite: `bash tests/dw-tests.sh` — confirm all 6 tests (4 original + 2 new) report PASS

**Checkpoint**: User Story 2 complete — all six deterministic examples now have automated coverage.

---

## Phase 5: User Story 3 — sprintf Replaced with snprintf (Priority: P3)

**Goal**: Replace all 5 `sprintf` calls in DW-authored source files with `snprintf(buf, BUFF_SIZE, ...)` to satisfy FR-006 and harden against buffer overflows.

**Independent Test**: `grep 'sprintf(' src/dw_*.c` returns zero results; full test suite passes.

- [X] T010 [P] [US3] Edit `src/dw_main.c` line 275: replace `sprintf(local_buffer, "sub%d_convexity", ...)` with `snprintf(local_buffer, BUFF_SIZE, "sub%d_convexity", ...)`

- [X] T011 [P] [US3] Edit `src/dw_main.c` line 496: replace `sprintf(local_buffer, "phase1_step_0.cpxlp")` with `snprintf(local_buffer, BUFF_SIZE, "phase1_step_0.cpxlp")`

- [X] T012 [P] [US3] Edit `src/dw_main.c` line 512: replace `sprintf(local_buffer, "phase1_step_%d.cpxlp", ...)` with `snprintf(local_buffer, BUFF_SIZE, "phase1_step_%d.cpxlp", ...)`

- [X] T013 [P] [US3] Edit `src/dw_main.c` line 611: replace `sprintf(local_buffer, "master_step_%d.cpxlp", ...)` with `snprintf(local_buffer, BUFF_SIZE, "master_step_%d.cpxlp", ...)`

- [X] T014 [P] [US3] Edit `src/dw_support.c` line 609: replace `sprintf(filename, "basis_iteration_%d", ...)` with `snprintf(filename, BUFF_SIZE, "basis_iteration_%d", ...)`

- [X] T015 [US3] Rebuild and run full test suite: `make && bash tests/dw-tests.sh` — confirm all 6 tests still pass

**Checkpoint**: User Story 3 complete — zero `sprintf` calls remain in DW source files.

---

## Phase 6: User Story 4 — malloc Null-Checks at Critical Sites (Priority: P3)

**Goal**: Add a `dw_oom_abort()` helper to `src/dw.h` and apply it at the 7 critical `malloc` call sites where a NULL return would cause immediate undefined behaviour (dereference in same function scope).

**Independent Test**: Full test suite passes; code review shows each critical `malloc` is followed by `dw_oom_abort(ptr, "context")`.

- [X] T016 [US4] Add `dw_oom_abort` static inline helper to `src/dw.h` immediately before the `#endif /*DECOMPOSE_H_*/` closing guard:
  ```c
  static inline void dw_oom_abort(void* ptr, const char* ctx) {
      if (!ptr) { fprintf(stderr, "dwsolver: out of memory in %s\n", ctx); exit(EXIT_FAILURE); }
  }
  ```
  Add `#include <stdlib.h>` near the top of `dw.h` if not already present.

- [X] T017 [P] [US4] Edit `src/dw_main.c` line 110: add `dw_oom_abort(local_buffer, "local_buffer");` immediately after the `malloc` for `local_buffer`

- [X] T018 [P] [US4] Edit `src/dw_main.c` line 118: add `dw_oom_abort(md, "md");` immediately after the `malloc` for `md`

- [X] T019 [P] [US4] Edit `src/dw_main.c` line 119: add `dw_oom_abort(globals, "globals");` immediately after the `malloc` for `globals`

- [X] T020 [P] [US4] Edit `src/dw_main.c` line 154: add `dw_oom_abort(threads, "threads");` immediately after the `malloc` for `threads`

- [X] T021 [P] [US4] Edit `src/dw_main.c` line 155: add `dw_oom_abort(sub_data, "sub_data");` immediately after the `malloc` for `sub_data`

- [X] T022 [P] [US4] Edit `src/dw_support.c` line 99: add `dw_oom_abort(D, "D");` immediately after the `malloc` for `D`

- [X] T023 [P] [US4] Edit `src/dw_support.c` line 173: add `dw_oom_abort(signals, "signals");` immediately after the `malloc` for `signals`

- [X] T024 [US4] Rebuild and run full test suite: `make && bash tests/dw-tests.sh` — confirm all 6 tests still pass

**Checkpoint**: User Story 4 complete — 7 critical allocation sites are guarded against OOM crashes.

---

## Final Phase: Polish & Cross-Cutting Concerns

- [X] T025 [P] Audit: `grep 'sprintf(' src/dw_*.c` — confirm zero results (FR-006 complete)
- [X] T026 [P] Audit: `grep -n 'malloc(' src/dw_main.c src/dw_support.c` — confirm each of the 7 critical sites listed in research.md is followed by `dw_oom_abort`
- [X] T027 Run `make clean && ./configure --enable-named-semaphores && make` — final clean build on macOS, zero errors
- [X] T028 Run `bash tests/dw-tests.sh` — final run, all 6 tests PASS
- [X] T029 Stage and commit all modified files (`src/dw.h`, `src/dw_globals.c`, `src/Makefile.am`, `tests/dw-tests.sh`) with message: `fix: resolve KD-001 Linux linker errors, extend tests, snprintf, malloc guards`
- [X] T030 Push branch: `git push origin 002-cross-platform-repair`

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 (baseline confirmed) — **BLOCKS all user story verifications**
- **US1 (Phase 3)**: Depends on Foundational — must rebuild after T002–T004
- **US2 (Phase 4)**: Depends on Foundational (fixed build)
- **US3 (Phase 5)**: Depends on Foundational (fixed build); independent of US2
- **US4 (Phase 6)**: Depends on Foundational (fixed build); independent of US2, US3
- **Polish (Final)**: Depends on all desired user stories complete

### User Story Dependencies

- **US1 (P1)**: Directly verifies the Foundational phase output — no other story dependency
- **US2 (P2)**: Needs working build; independent of US1's Linux verification
- **US3 (P3)**: Needs working build; independent of US2, US4
- **US4 (P3)**: Needs `dw_oom_abort` helper in `dw.h`; independent of US2, US3

### Within Each Phase

- T002 (create file) before T003/T004 (edits referencing it) — though file edits can be done in same session
- T010–T014 (all [P]) can be done simultaneously — each touches a different location
- T017–T023 (all [P]) can be done simultaneously — each touches a different line

### Parallel Execution Examples

**Foundational phase** (fastest path):
```
T002 (create dw_globals.c)
  → T003 [P] (edit dw.h)    } can overlap as file edits
  → T004 [P] (edit Makefile) }
→ T005 (make clean && make)
→ T006 (run tests)
```

**US3 — all 5 sprintf replacements in parallel**:
```
T010 [P] dw_main.c:275
T011 [P] dw_main.c:496   } all simultaneously
T012 [P] dw_main.c:512
T013 [P] dw_main.c:611
T014 [P] dw_support.c:609
  → T015 (rebuild + test)
```

**US4 — all 7 malloc guards in parallel**:
```
T016 (add helper to dw.h first)
  → T017–T023 [P] all simultaneously
  → T024 (rebuild + test)
```

---

## Implementation Strategy

**MVP = Phase 1 + Phase 2 + Phase 3 (US1)**: The KD-001 fix alone makes the software
buildable on Linux for the first time. All other phases are independent improvements.

**Incremental delivery order**: US1 → US2 → US3 + US4 (US3 and US4 are equivalent
priority and fully independent — either can go first).

**Total tasks**: 30
**Parallelizable tasks**: 17 (marked [P])
**Per user story**: US1 = 2 tasks, US2 = 3 tasks, US3 = 6 tasks, US4 = 9 tasks
