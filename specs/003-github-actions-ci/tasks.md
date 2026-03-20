# Tasks: GitHub Actions CI — Multi-Platform Build & Test

**Input**: Design documents from `specs/003-github-actions-ci/`
**Prerequisites**: plan.md ✅ spec.md ✅ research.md ✅ data-model.md ✅ contracts/ ✅

**Organization**: Tasks are grouped by user story. Each story phase is independently
deliverable. No tests are generated (not requested in spec).

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to

---

## Phase 1: Setup

**Purpose**: Ensure the feature branch is rebased on master and ready for implementation.
This is a CI spec; there is no project initialization — the repository already exists.

- [ ] T001 Rebase `003-github-actions-ci` onto `master` to pick up any recent changes

**Checkpoint**: Branch is up to date. Implementation can begin.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Remove the broken legacy workflow before creating the replacement. This
prevents duplicate CI runs and eliminates the incorrect macOS configuration that omits
`--enable-named-semaphores`.

**⚠️ CRITICAL**: Must complete before creating `ci.yml` — having both files active
simultaneously on `master` would double-trigger CI on every push.

- [ ] T002 Delete `.github/workflows/build.yaml` (superseded; legacy macOS-only, missing `--enable-named-semaphores`, outdated action pins)

**Checkpoint**: `build.yaml` removed. Clean slate for the new workflow file.

---

## Phase 3: User Story 1 — PR and Push Gating on macOS and Linux (Priority: P1) 🎯 MVP

**Goal**: Every PR against `master` and every push to `master` automatically triggers
build-and-test on macOS and Linux. Failures are visible on the PR checks tab.

**Independent Test**: Open any PR against `master` (or push directly) — the GitHub
Actions "Checks" tab shows `CI / build-macos` and `CI / build-linux` appearing within
60 seconds, completing within 15 minutes, and reporting pass/fail correctly.

### Implementation for User Story 1

- [ ] T003 [US1] Create `.github/workflows/ci.yml` with top-level fields:
  `name: CI`, triggers `on.push.branches: [master]` and
  `on.pull_request.branches: [master]`
- [ ] T004 [P] [US1] Add `build-macos` job to `.github/workflows/ci.yml`:
  `runs-on: macos-latest`, `timeout-minutes: 30`, steps: checkout@v4,
  `./configure --enable-named-semaphores && make`,
  `echo "$GITHUB_WORKSPACE/src" >> $GITHUB_PATH`,
  `cd tests && bash dw-tests.sh`
- [ ] T005 [P] [US1] Add `build-linux` job to `.github/workflows/ci.yml`:
  `runs-on: ubuntu-latest`, `timeout-minutes: 30`, steps: checkout@v4,
  `./configure && make`,
  `echo "$GITHUB_WORKSPACE/src" >> $GITHUB_PATH`,
  `cd tests && bash dw-tests.sh`

**Checkpoint**: `ci.yml` contains two required jobs. Push the branch, open a PR, confirm
`CI / build-macos` and `CI / build-linux` both appear and pass in the Checks tab.

---

## Phase 4: User Story 2 — Windows CI Job (Priority: P2)

**Goal**: A Windows job appears on every PR and every push to `master`. It attempts to
build and test on MSYS2 MinGW-w64. Failures are visible but non-blocking via
`continue-on-error: true`.

**Independent Test**: With `ci.yml` already triggering on PRs (from Phase 3), a Windows
job (`CI / build-windows`) appears alongside macOS/Linux. If it fails, the PR remains
mergeable; the failed job is visible in the job log.

### Implementation for User Story 2

- [ ] T006 [US2] Add `build-windows` job to `.github/workflows/ci.yml`:
  `runs-on: windows-latest`, `timeout-minutes: 30`, `continue-on-error: true`,
  `defaults.run.shell: msys2 {0}`, steps: checkout@v4,
  `msys2/setup-msys2@v2` with `msystem: MINGW64 update: true cache: true install: mingw-w64-x86_64-gcc make`,
  `./configure && make`,
  `echo "$GITHUB_WORKSPACE/src" >> $GITHUB_PATH`,
  `cd tests && bash dw-tests.sh`

**Checkpoint**: All three jobs (`build-macos`, `build-linux`, `build-windows`) appear on
the PR checks tab. Windows job failure does not block merge.

---

## Phase 5: User Story 3 — Branch Protection on `master` (Priority: P1)

**Goal**: Configure branch protection rules on `master` so that `CI / build-macos` and
`CI / build-linux` are required status checks. No PR can be merged to `master` until
both pass. This is enforced at the repository level.

**Independent Test**: With branch protection active, attempt to merge a PR with a failing
macOS or Linux CI job — the "Merge pull request" button is greyed out with the message
"Required status checks have not passed."

**Prerequisite**: At least one CI run must have completed on the branch (so the check
names `CI / build-macos` and `CI / build-linux` appear in GitHub's status check search).

### Implementation for User Story 3

- [ ] T007 [US3] Open a pull request from `003-github-actions-ci` to `master` to trigger
  the first CI run and register the status check names in GitHub's system
- [ ] T008 [US3] In GitHub repository settings → Branches → `master` branch protection rule:
  enable "Require status checks to pass before merging", add required checks
  `CI / build-macos` and `CI / build-linux` (do NOT add `CI / build-windows`)
- [ ] T009 [US3] Verify enforcement: confirm that a PR with a failing `build-macos` or
  `build-linux` check cannot be merged (Merge button disabled)

**Checkpoint**: Branch protection is active. `master` is protected. The CI spec is
fully operational.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Final validation of the complete implementation against spec success
criteria and quickstart.

- [ ] T010 [P] Validate `ci.yml` YAML syntax locally: `python3 -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml'))"` — must parse without error
- [ ] T011 [P] Verify `build.yaml` is absent from `.github/workflows/` — confirm no duplicate workflow exists
- [ ] T012 Run through the quickstart.md checklist in `specs/003-github-actions-ci/quickstart.md` to confirm all implementation steps are complete
- [ ] T013 Confirm status check names in branch protection exactly match `CI / build-macos` and `CI / build-linux` (case-sensitive, space around `/`)

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately
- **Phase 2 (Foundational)**: After Phase 1 — BLOCKS Phases 3, 4, 5
- **Phase 3 (US1 — macOS + Linux)**: After Phase 2 — MVP deliverable; required before Phase 5
- **Phase 4 (US2 — Windows)**: After Phase 2 — independent of Phase 3; can run in parallel with Phase 3
- **Phase 5 (US3 — Branch Protection)**: After Phase 3 (CI must have run once in GH)
- **Phase 6 (Polish)**: After Phases 3, 4, 5

### User Story Dependencies

- **US1 (P1)**: Depends on Phase 2 (legacy file removed) only
- **US2 (P2)**: Depends on Phase 2 only — adds a job to the same file as US1
- **US3 (P1)**: Depends on US1 being live (CI must have run to register check names)

### Task Dependencies Within Phases

- T003 before T004 and T005 (file must exist before jobs are added)
- T004 and T005 are [P] — same file, but non-overlapping sections (separate jobs)
- T006 is [P] relative to T004/T005 — adds a third non-overlapping section to `ci.yml`
- T007 before T008 (PR must exist and CI must have run before check names are searchable)
- T010 and T011 are [P] — independent validation steps

### Note on T004/T005/T006 Parallelism

All three jobs can be written in a single edit pass to `ci.yml` if working alone:
- Write `build-macos` (T004), `build-linux` (T005), and `build-windows` (T006) in one commit
- This satisfies all three tasks simultaneously

---

## Parallel Execution Example: Phases 3 & 4 Together

Since both Phase 3 and Phase 4 add jobs to the same file (`ci.yml`), and the file is
created in T003, the fastest path is to write all three jobs in one editing session:

```
T001 → T002 → T003 → [T004 + T005 + T006 in one pass] → T007 → T008 → T009 → [T010 + T011] → T012 → T013
```

**Total task count**: 13 tasks
**Parallelizable**: T004+T005+T006 (write together); T010+T011

---

## Implementation Strategy

**MVP (Phases 1–3, Tasks T001–T005)**:
- Removes broken `build.yaml`
- Creates `ci.yml` with macOS and Linux jobs
- Provides immediate CI coverage on all PR

**Full delivery adds (Phase 4–6)**:
- Windows visibility job
- Branch protection enforcement
- Final validation

The MVP alone satisfies US1 (P1) and partially US3 (P1 — CI visible but not yet enforced
until branch protection rules are set in T008).
