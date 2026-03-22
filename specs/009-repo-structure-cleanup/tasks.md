# Tasks: Repository Top-Level Cleanup

**Feature**: 009-repo-structure-cleanup  
**Branch**: `009-repo-structure-cleanup`  
**Input**: [spec.md](spec.md) · [plan.md](plan.md)

---

## Phase 1: Setup

**Purpose**: Verify baseline build and test state before making any changes, so regressions are immediately attributable.

- [x] T001 Confirm clean baseline: run `./configure && make` and `cd tests && PATH="$PWD/../src:$PATH" ./dw-tests.sh`; record 7/7 PASS in a tmp note
- [x] T002 Count current tracked root-level files: `git ls-files | grep -v '/' | wc -l`; record for SC-001 comparison

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Git hygiene tasks with zero build impact. Must be done before moving files to keep commits atomic.

**⚠️ Each task here is independently committable.**

- [x] T003 [US3] Remove committed editor backup files from git: `git rm configure~ config.h.in~`
- [x] T004 [US3] Add `*~` pattern to `.gitignore` (append after existing `*.plist` entry)
- [x] T005 [US3] Commit: `chore(009): remove editor backup files and ignore *~ pattern`
- [x] T006 [US3] Verify: `git ls-files | grep '~'` returns no output; `git status` is clean

**Checkpoint**: No backup files in git. US3 (editor clutter) is now complete and independently verifiable.

---

## Phase 3: User Story 2 — Move GLPK Files to `third-party/glpk/` (Priority: P2)

**Goal**: All GLPK attribution and patch files live in `third-party/glpk/`; root contains none of them.

**Independent Test**: `ls third-party/glpk/` shows 7 files (6 GLPK + README); `git ls-files | grep 'GLPK_\|glpk-4.44'` shows only `third-party/glpk/` paths; `./configure && make` still succeeds.

- [x] T007 [US2] Create `third-party/glpk/README` explaining the directory contains GLPK 4.44 attribution and patch files bundled with the original dwsolver source distribution
- [x] T008 [US2] Move files with `git mv`: `GLPK_AUTHORS`, `GLPK_INSTALL`, `GLPK_NEWS`, `GLPK_README`, `GLPK_THANKS`, `glpk-4.44.ThreadReady.patch` → `third-party/glpk/`
- [x] T009 [US2] Update `Makefile.am` `EXTRA_DIST` (lines 7–8): replace root-relative GLPK filenames with `third-party/glpk/` prefixed paths so `make dist` still includes them
- [x] T010 [US2] Build verify: `./configure && make` succeeds; `cd tests && PATH="$PWD/../src:$PATH" ./dw-tests.sh` shows 7/7 PASS
- [x] T011 [US2] Commit: `chore(009): move GLPK attribution files to third-party/glpk/`

**Checkpoint**: GLPK files are in `third-party/glpk/`. US2 is fully verifiable without any further tasks.

---

## Phase 4: User Story 1 — Move Autoconf Aux Scripts to `build-aux/` (Priority: P1)

**Goal**: The 7 autoconf auxiliary scripts are under `build-aux/`; root contains none of them. Build and tests pass unchanged.

**Independent Test**: `ls build-aux/` shows all 7 scripts; root `ls` shows none of them; `./configure && make && cd tests && PATH="$PWD/../src:$PATH" ./dw-tests.sh` all pass; CI is green on all 6 check types.

### configure.ac change

- [x] T012 [US1] In `configure.ac`, add `AC_CONFIG_AUX_DIR([build-aux])` on a new line immediately before `AM_INIT_AUTOMAKE` (line 9)

### Regenerate with autoreconf

- [x] T013 [US1] Run `autoreconf -fi` from repo root; confirm `build-aux/` is created and contains: `compile`, `config.guess`, `config.sub`, `depcomp`, `install-sh`, `ltmain.sh`, `missing`
- [x] T014 [US1] Remove the old root-level copies: `git rm compile config.guess config.sub depcomp install-sh ltmain.sh missing`
- [x] T015 [US1] Stage new files: `git add build-aux/ configure Makefile.in aclocal.m4 src/Makefile.in tests/Makefile.in`
- [x] T016 [US1] Build verify (macOS): `./configure && make` succeeds; `cd tests && PATH="$PWD/../src:$PATH" ./dw-tests.sh` shows 7/7 PASS
- [x] T017 [US1] Commit: `chore(009): move autoconf aux scripts to build-aux/ via AC_CONFIG_AUX_DIR`

### Update CI workflows

- [x] T018 [US1] In `.github/workflows/ci-linux.yml`: update all three `touch` steps — change `ltmain.sh` → `build-aux/ltmain.sh` in each `touch aclocal.m4 configure config.h.in ltmain.sh ...` line (lines 20, 39, 58)
- [x] T019 [US1] In `.github/workflows/ci-macos.yml`: update the `touch` step — change `ltmain.sh` → `build-aux/ltmain.sh`
- [x] T020 [US1] In `.github/workflows/ci-windows.yml`: update the `touch` step (line 34) — `build-aux/ltmain.sh`; update the two `cp` commands (lines 35–36) to copy to `build-aux/config.guess` and `build-aux/config.sub` instead of root

### Update Dockerfile

- [x] T021 [US1] In `Dockerfile` (lines 20–23): update `touch ... ltmain.sh` → `build-aux/ltmain.sh`; update the two `cp` commands to target `build-aux/config.guess` and `build-aux/config.sub`

### Final verification

- [x] T022 [US1] Commit: `chore(009): update CI workflows and Dockerfile for build-aux/ paths`

**Checkpoint**: Aux scripts are in `build-aux/`. CI should go green on all 6 check types after push.

---

## Phase 5: Polish & Cross-Cutting Concerns

**Purpose**: Documentation update and final validation.

- [x] T023 Update `README.md` project structure description: add a mention of `build-aux/` (autoconf auxiliary scripts) and `third-party/` (GLPK attribution files) alongside the existing directory descriptions
- [x] T024 Final root count: `git ls-files | grep -v '/' | wc -l`; confirm reduction of ≥ 13 files vs baseline recorded in T002
- [x] T025 Commit: `docs(009): update README for new build-aux/ and third-party/ directories`
- [x] T026 Push branch and verify all 6 CI checks pass (macOS, Linux, Linux ASan+UBSan, Linux TSan, Docker, Windows)

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2 / US3)**: No build dependencies — can be done first to get it out of the way
- **Phase 3 (US2 / GLPK)**: Optionally after Phase 2; no dependency on Phase 2
- **Phase 4 (US1 / build-aux)**: Independent of Phase 3; most impactful, requires autoreconf
- **Phase 5 (Polish)**: After Phases 3 and 4 are both committed

### User Story Dependencies

- **US3 (P3 — git hygiene)**: Fully independent; complete first for quick win
- **US2 (P2 — GLPK move)**: Independent of US3 and US1; Makefile.am update is the only file interaction
- **US1 (P1 — build-aux)**: Independent of US2 and US3; blocking on autoreconf + CI update

### Parallel Opportunities

```
[T001, T002]                     ← Setup, run together

[T003→T006]                      ← US3 (git hygiene, sequential internal steps)

[T007→T011] can overlap with     ← US2 and US1 are independent of each other;
[T012→T022]                         can be done in any order or simultaneously
                                    by different contributors

[T023→T026]                      ← Polish, after both US1 and US2 are committed
```

### Within US1 (most complex)

```
T012 (configure.ac edit)
  → T013 (autoreconf -fi)
    → T014 (git rm old aux scripts)
    → T015 (git add build-aux/)
      → T016 (build verify)
        → T017 (commit)
          → [T018, T019, T020, T021] (CI + Dockerfile, parallel)
            → T022 (commit)
```

---

## Implementation Strategy

**MVP**: US3 (T003–T006) delivers immediate, zero-risk value in one commit.  
**Then US2** (T007–T011) removes 6 GLPK files from root in one commit.  
**Then US1** (T012–T022) is the highest-value, highest-effort step — requires local autoreconf run and CI validation.  
**Finally Polish** (T023–T026) wraps up docs and confirms success criteria.

Total tasks: 26 across 5 phases.
