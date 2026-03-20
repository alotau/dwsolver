# Tasks: Containerized Build for dwsolver

**Input**: Design documents from `specs/004-docker-containerized-build/`
**Branch**: `004-docker-containerized-build`
**Date**: 2026-03-20
**Prerequisites**: plan.md ✅, spec.md ✅, research.md ✅, data-model.md ✅, contracts/docker-interface.md ✅, quickstart.md ✅

**Tests**: Not requested — no test tasks generated.

**Organization**: Tasks grouped by user story (US1 → US2 → US3 in priority order).

## Format: `[ID] [P?] [Story?] Description`

- **[P]**: Parallelizable (different files / no dependency on an incomplete sibling task)
- **[Story]**: User story this task belongs to

---

## Phase 1: Setup

**Purpose**: Confirm local prerequisites before any file creation.

- [X] T001 Confirm branch `004-docker-containerized-build` is checked out and `docker info` succeeds (Docker Engine ≥ 20.10 available locally)

**Checkpoint**: Branch active, Docker running — proceed to Phase 2.

---

## Phase 2: Foundational (Blocking Prerequisite)

**Purpose**: Create `.dockerignore` — required before `docker build` can correctly filter the build context. This is a prerequisite for both US1 (the build) and US3 (artifact exclusion).

**⚠️ CRITICAL**: `.dockerignore` MUST exist before any `docker build` attempt.

- [X] T002 Create `.dockerignore` at repository root excluding `.git/`, `*.o`, `*.a`, `*.Po`, `src/dwsolver`, `config.log`, `config.status`, `config.h`, `stamp-h1`, `libtool`, `autom4te.cache/`, `**/out_terminal`, `**/relaxed_solution`, `**/rs_sorted`, `**/done.cpxlp`, `**/pre_master.cpxlp`, `tmp/`, `.vscode/`, `.specify/` (per data-model.md build context exclusion table and contracts/docker-interface.md)
  - **Note**: Used `**/*.o`, `**/*.Po` (recursive glob) — Docker's `*.o` only matches root-level files.

**Checkpoint**: `.dockerignore` in place — US1, US2, US3 implementation can now begin.

---

## Phase 3: User Story 1 — Run Solver Against Local Data Files (Priority: P1) 🎯 MVP

**Goal**: Any user with Docker can build the image and run the solver against their own data using a single `docker run` command with a volume mount — no build toolchain required.

**Independent Test**: `docker build -t dwsolver .` succeeds, image is ≤ 200 MB, and `docker run --rm -v "$(pwd)/examples/book_bertsimas:/data" dwsolver --no-write-final-master --quiet -g /data/guidefile` exits 0.

### Implementation for User Story 1

- [X] T003 [US1] Create `Dockerfile` at repository root with two-stage build: builder stage (`FROM ubuntu:24.04 AS builder`, `apt-get install -y --no-install-recommends build-essential automake`, `WORKDIR /build`, `COPY . .`, `touch aclocal.m4 configure config.h.in ltmain.sh Makefile.in src/Makefile.in`, `./configure && make`) and runner stage (`FROM ubuntu:24.04 AS runner`, `COPY --from=builder /build/src/dwsolver /usr/local/bin/dwsolver`, `ENTRYPOINT ["/usr/local/bin/dwsolver"]`) — per quickstart.md Step 2 and contracts/docker-interface.md
- [X] T004 [US1] Run `docker build -t dwsolver .` from repository root and confirm it completes with exit code 0 and no error output
  - **Note**: Also refreshed `config.guess`/`config.sub` from system automake to support aarch64 (Apple Silicon Docker).
- [X] T005 [P] [US1] Run `docker run --rm dwsolver 2>&1 | grep -q "Usage:"` and confirm it exits 0 (smoke test: binary executes and prints usage when invoked with no arguments)
- [X] T006 [P] [US1] Run `docker image inspect dwsolver --format '{{.Size}}'` and confirm the reported byte count is ≤ 209715200 (200 MB), satisfying SC-002 and FR-006
  - **Result**: 99.2 MB ✅
- [X] T007 [US1] Run `docker run --rm -v "$(pwd)/examples/book_bertsimas:/data" -w /data dwsolver --no-write-final-master --quiet -g /data/guidefile` and confirm exit code 0, satisfying FR-001/FR-004/SC-003 (solver reads from mounted volume and solves correctly)

**Checkpoint**: User Story 1 fully functional. A user with only Docker can solve the Bertsimas example. ✅ MVP complete.

---

## Phase 4: User Story 2 — Build the Image in CI (Priority: P2)

**Goal**: A `docker-build` job in `ci.yml` detects Dockerfile regressions on every PR without blocking the required `build-linux` / `build-macos` checks.

**Independent Test**: Inspect `.github/workflows/ci.yml` — a `docker-build` job exists, `continue-on-error: true` is set, and the smoke-test step (`grep -q "Usage:"`) is present. Open a PR and confirm the job appears in the checks list.

### Implementation for User Story 2

- [X] T008 [US2] Add `docker-build` job to `.github/workflows/ci.yml` after the `build-windows` job: `runs-on: ubuntu-latest`, `timeout-minutes: 30`, `continue-on-error: true`, with steps `actions/checkout@v4` → `docker build -t dwsolver .` → `docker run --rm dwsolver 2>&1 | grep -q "Usage:"` (per contracts/docker-interface.md CI contract and research.md R7)

**Checkpoint**: CI job defined. Will become active when the PR is opened in Phase Final. ✅

---

## Phase 5: User Story 3 — Exclude Build Artifacts from the Image (Priority: P3)

**Goal**: Confirm the final runnable image contains no C source files, object files, or `.git` directory — only the binary. Enforced structurally by the multi-stage Dockerfile and the `.dockerignore`.

**Independent Test**: `.dockerignore` contains all critical patterns; `docker run --rm --entrypoint /bin/sh dwsolver -c 'find / -name "*.c" 2>/dev/null | wc -l'` outputs `0`.

### Implementation for User Story 3

- [X] T009 [P] [US3] Verify `.dockerignore` completeness: run `grep -cE '`.git|`.o|src/dwsolver' .dockerignore` and confirm all three patterns are present, satisfying FR-005 and US3 acceptance scenario 1
- [X] T010 [P] [US3] Verify runtime image contains no C source files: run `docker run --rm --entrypoint /bin/sh dwsolver -c 'find / -name "*.c" 2>/dev/null | wc -l'` and confirm output is `0`, satisfying FR-006 (multi-stage isolation) and US3 acceptance scenario 2

**Checkpoint**: Image is confirmed lean — no source files, no artifacts, binary only. ✅

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Commit all deliverables to the branch and open a PR.

- [ ] T011 Commit `Dockerfile`, `.dockerignore`, and `.github/workflows/ci.yml` to branch `004-docker-containerized-build` with message `feat: add Docker containerized build (Spec 004)` (see quickstart.md Step 5 for full commit message)
- [ ] T012 Push branch and open a PR against `master`; confirm the `docker-build` job appears in the CI checks list and reports a result (pass or informational fail) without blocking `build-linux` or `build-macos`

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately
- **Phase 2 (Foundational)**: Depends on Phase 1. Blocks T003–T010.
- **Phase 3 (US1)**: Depends on Phase 2. T003 → T004 → {T005, T006, T007} (T005 and T006 are parallel after T004)
- **Phase 4 (US2)**: Depends on Phase 2 only. Can start after T002 (independent of US1 build verification)
- **Phase 5 (US3)**: Depends on T004 (image must be built). T009 and T010 are parallelizable.
- **Phase 6 (Polish)**: Depends on T003–T010 all complete

### User Story Dependencies

- **US1 (P1)**: After T002 — no dependency on US2 or US3
- **US2 (P2)**: After T002 — no dependency on US1 (editing `ci.yml` is independent of building the image)
- **US3 (P3)**: After T004 — verifications depend on the image being built

### Parallel Opportunities

```
T001 → T002 → T003 → T004 ─┬─ T005 (parallel)
                             ├─ T006 (parallel)
                             └─ T007

          T002 → T008        (independent of T003–T007)

          T004 ─┬─ T009 (parallel)
                └─ T010 (parallel)

T003 + T007 + T008 + T009 + T010 all complete → T011 → T012
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Confirm prerequisites (T001)
2. Complete Phase 2: Create `.dockerignore` (T002) — **blocks everything**
3. Complete Phase 3: Create Dockerfile and verify all three acceptance scenarios (T003–T007)
4. **STOP and VALIDATE**: `docker run --rm -v "$(pwd)/examples/book_bertsimas:/data" dwsolver --no-write-final-master --quiet -g /data/guidefile` exits 0
5. If MVP is sufficient, skip Phase 4–5 and go straight to Phase 6 (commit + PR)

### Full Delivery (All Three User Stories)

6. Complete Phase 4: Add CI job (T008) — quick edit to `ci.yml`
7. Complete Phase 5: Verify artifact exclusion (T009–T010)
8. Complete Phase 6: Commit + PR (T011–T012)

---

## Summary

| Phase | Tasks | User Story | Parallelizable |
|-------|-------|----------|----------------|
| 1 — Setup | T001 | — | No |
| 2 — Foundational | T002 | — | No |
| 3 — US1 (P1) 🎯 | T003–T007 | US1 | T005, T006 after T004 |
| 4 — US2 (P2) | T008 | US2 | Parallel with US1 from T002 |
| 5 — US3 (P3) | T009–T010 | US3 | Both parallel after T004 |
| 6 — Polish | T011–T012 | — | No |
| **Total** | **12 tasks** | | |

**Suggested MVP scope**: Phases 1–3 + Phase 6 (T001–T007, T011–T012 = 9 tasks, delivers US1 completely)
