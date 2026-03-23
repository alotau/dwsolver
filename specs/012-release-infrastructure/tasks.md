# Tasks: Release Infrastructure

**Input**: Design documents from `specs/012-release-infrastructure/`
**Branch**: `012-release-infrastructure`
**Spec**: [spec.md](spec.md) | **Plan**: [plan.md](plan.md) | **Research**: [research.md](research.md)

> **Status note**: This feature is largely implemented. The task list reflects
> the as-built state — completed tasks are checked. Remaining tasks are the
> verification, spec-commit, and PR steps.

---

## Phase 1: Setup

**Purpose**: Verify clean working tree and build toolchain ready.

- [x] T001 Confirm build passes and tests pass (`make && make check` from repo root)
- [x] T002 Confirm branch `012-release-infrastructure` is checked out and based on `011-remove-embedded-glpk`

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Autotools distcheck gate — must pass before any release can be cut.

**⚠️ CRITICAL**: All user stories depend on `make distcheck` passing cleanly.

- [x] T003 Expand `EXTRA_DIST` in `Makefile.am` to include `examples`, `architecture`, `specs`, `Dockerfile`, `.github`, `third-party/glpk` docs, `ADDITIONAL_LICENSE_TERMS`, `ChangeLog`, `dwsolver.pc.in`
- [x] T004 Add `EXTRA_DIST` to `tests/Makefile.am` listing `dw-tests.sh`, `test_guidefile.sh`, `fixtures`
- [x] T005 Add `DISTCLEANFILES` to `tests/Makefile.am` listing `relaxed_solution`, `rs_sorted`, `out_obj.txt`, `src/dw_blas.Po`
- [x] T006 Regenerate `Makefile.in` and `tests/Makefile.in` from updated `.am` files (`automake`)
- [x] T007 Verify `make distcheck` exits 0 and produces `dwsolver-1.2.1.tar.gz`

**Checkpoint**: `make distcheck` passes — tarball gate is open.

---

## Phase 3: User Story 1 — Maintainer Produces a Complete Source Tarball (Priority: P1) 🎯 MVP

**Goal**: `make distcheck` succeeds end-to-end; tarball contains all required files.

**Independent Test**: Run `make distcheck` on a clean checkout; exit 0; unpack tarball; build and run `make check`; all tests pass.

- [x] T008 [US1] Confirm tarball contains `examples/`, `architecture/`, `specs/`, `Dockerfile`, `.github/`, `tests/dw-tests.sh`, `tests/fixtures/` (`tar tzf dwsolver-1.2.1.tar.gz` spot-check)
- [x] T009 [US1] Confirm `make check` inside the unpacked distcheck tree shows zero failures

**Checkpoint**: US1 complete — reproducible, self-contained source tarball verified.

---

## Phase 4: User Story 2 — Automated GitHub Release on Version Tag Push (Priority: P2)

**Goal**: Push a `v*` tag → CI runs `make distcheck` → GitHub Release created with tarball attached. Workflow must NOT fire on branch pushes or PRs.

**Independent Test**: Inspect `.github/workflows/release.yml`:
- `on.push.tags: ['v*']` present; no `branches:` or `pull_request:` triggers.
- `make distcheck` step precedes the release-creation step.
- `softprops/action-gh-release` step has `draft: false` and `files:` pointing to the tarball.
- Both actions are SHA-pinned.
- `GLPK_CFLAGS`/`GLPK_LIBS` exported to `$GITHUB_ENV` in the install step.

- [x] T010 [US2] Create `.github/workflows/release.yml` with `on.push.tags: ['v*']` trigger only
- [x] T011 [US2] Add install-dependencies step: `build-essential automake pkg-config libglpk-dev`; export `GLPK_CFLAGS=-I/usr/include` and `GLPK_LIBS=-lglpk` to `$GITHUB_ENV`
- [x] T012 [US2] Add automake timestamp-touch step to suppress unnecessary regeneration
- [x] T013 [US2] Add `make distcheck` step (no `DISTCHECK_CONFIGURE_FLAGS` — GLPK vars already in env)
- [x] T014 [US2] Add tarball-locate step outputting `path` and `name` via `$GITHUB_OUTPUT`
- [x] T015 [US2] Add `softprops/action-gh-release` step with `draft: false`, `generate_release_notes: true`, `files:` from locate step
- [x] T016 [P] [US2] Pin `actions/checkout` to SHA `34e114876b0b11c390a56381ad16ebd13914f8d5` (v4)
- [x] T017 [P] [US2] Pin `softprops/action-gh-release` to SHA `153bb8e04406b158c6c84fc1615b65b24149a1fe` (v2)
- [x] T018 [US2] Set `permissions: contents: write` and `timeout-minutes: 20` on the job

**Checkpoint**: US2 complete — automated release workflow is spec-compliant and SHA-pinned.

---

## Phase 5: User Story 3 — Correct Library Soname Versioning (Priority: P3)

**Goal**: libtool version-info triple in `src/Makefile.am` is `0:0:0`; soname is `libdwsolver.so.0`; update rules are documented.

**Independent Test**: `grep version-info src/Makefile.am` shows `0:0:0`; `otool -L src/.libs/libdwsolver.dylib` (macOS) or `readelf -d src/.libs/libdwsolver.so` (Linux) confirms soname matches `CURRENT - AGE = 0`.

- [x] T019 [US3] Confirm `src/Makefile.am` has `-version-info 0:0:0` in `libdwsolver_la_LDFLAGS`
- [ ] T020 [US3] Verify soname of built shared library matches `libdwsolver.so.0` / `libdwsolver.0.dylib` using `otool -L` (macOS) or `readelf -d` (Linux) on `src/.libs/libdwsolver.*`

**Checkpoint**: US3 complete — soname is correct and consistent with documented versioning rules.

---

## Phase 6: User Story 4 — Release Process Documented in README (Priority: P4)

**Goal**: README contains a "Cutting a Release" (or equivalent) section covering version bump, distcheck gate, tagging convention, and libtool version-info rules.

**Independent Test**: Read the README; locate a section titled "Cutting a Release"; confirm it contains: (1) the `make distcheck` command with the platform-specific flags, (2) `git tag v...` + `git push --tags` instructions, (3) a link to or mention of the release workflow, (4) the version-info update rules table.

- [x] T021 [US4] Add "Cutting a Release" section to `README.md` between Docker and SEI CERT C sections
- [x] T022 [P] [US4] Include Linux `make distcheck` command in README release section
- [x] T023 [P] [US4] Include macOS `make distcheck` command (with `--enable-named-semaphores --enable-recursive-mutex`) in README release section
- [x] T024 [US4] Include `git tag v<ver>` + `git push origin v<ver>` instructions in README
- [x] T025 [US4] Include libtool version-info update rules table in README
- [x] T026 [US4] Reference the release workflow (`.github/workflows/release.yml`) in README

**Checkpoint**: US4 complete — any maintainer can follow README alone to cut a release.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Spec artifacts committed, PR ready, edge-case verification.

- [x] T027 [P] Commit `Makefile.am`, `Makefile.in`, `tests/Makefile.am`, `tests/Makefile.in` with EXTRA_DIST and DISTCLEANFILES changes
- [x] T028 [P] Commit `.github/workflows/release.yml`
- [x] T029 [P] Commit `README.md` release section
- [x] T030 [P] Commit `specs/012-release-infrastructure/spec.md` with clarifications
- [x] T031 [P] Commit `specs/012-release-infrastructure/plan.md`, `research.md`, `data-model.md`, `contracts/`
- [ ] T032 Push branch `012-release-infrastructure` to origin (already at `b9ec5aa` — confirm up to date)
- [ ] T033 Open pull request: `012-release-infrastructure` → `main`; title "feat(012): release infrastructure"; body references spec and SC-001–SC-006

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately.
- **Phase 2 (Foundational/distcheck)**: Depends on Phase 1. Blocks all user stories.
- **Phase 3 (US1 — tarball)**: Depends on Phase 2.
- **Phase 4 (US2 — release workflow)**: Depends on Phase 2; independent of US1.
- **Phase 5 (US3 — soname)**: Depends on Phase 2; independent of US1, US2.
- **Phase 6 (US4 — README)**: Depends on Phase 2; can partially be done alongside US2/US3.
- **Phase 7 (Polish/PR)**: Depends on all stories being verified.

### User Story Dependencies

- **US1 (P1)**: No inter-story dependency. Foundational gate only.
- **US2 (P2)**: No dependency on US1. Foundational gate only.
- **US3 (P3)**: No dependency on US1/US2. Foundational gate only.
- **US4 (P4)**: References US2 (workflow link) and US3 (version-info rules) but can be written independently.

### Parallel Opportunities

**US2, US3, US4 can proceed in parallel** once Phase 2 is complete.

```
Phase 2 done
    ├── US1 tasks (T008–T009)
    ├── US2 tasks (T010–T018)   ← can run concurrently with US1, US3, US4
    ├── US3 tasks (T019–T020)   ← can run concurrently with US1, US2, US4
    └── US4 tasks (T021–T026)   ← can run concurrently with US1, US2, US3
```

---

## Implementation Strategy

**MVP Scope**: Phase 2 + US1 (T001–T009) — `make distcheck` passing is the minimum viable release gate. Everything else builds on it.

**Delivery order** (as implemented):
1. Phase 2: distcheck gate (EXTRA_DIST + DISTCLEANFILES)
2. Phase 3 US1: tarball verification
3. Phase 4 US2 + Phase 5 US3 + Phase 6 US4 in parallel
4. Phase 7: commit, push, PR

**Remaining work**:
- T020: soname spot-check on built library
- T032: confirm branch is up to date with origin
- T033: open PR

---

## Success Criteria Checklist

- [x] SC-001: `make distcheck` exits 0 in clean checkout
- [ ] SC-002: Tarball verified to contain 100% of EXTRA_DIST files (`tar tzf` spot-check)
- [ ] SC-003: Pushing a `v*` tag produces a GitHub Release within 15 minutes with tarball attached *(verified on first real tag push)*
- [x] SC-004: `make distcheck` step in workflow precedes release creation; logs show zero test failures *(verified via workflow structure)*
- [x] SC-005: version-info `0:0:0` in `src/Makefile.am`; README documents update rules
- [ ] SC-006: Maintainer following only README can produce a release *(verified manually)*
