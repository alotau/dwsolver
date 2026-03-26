# Tasks: OpenSSF Best Practices Badge (Passing Level)

**Input**: Design documents from `/specs/015-openssf-passing-badge/`  
**Prerequisites**: plan.md ✅ | spec.md ✅ | research.md ✅ | data-model.md ✅ | contracts/ ✅ | quickstart.md ✅  
**Branch**: `015-openssf-passing-badge`  
**Total tasks**: 33  
**Tests**: Not requested — no test tasks generated.

## Format: `[ID] [P?] [Story?] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[US1/US2/US3]**: Owning user story

---

## Phase 1: Setup

**Purpose**: Confirm branch state and verify no CI changes are required (FR-007 already satisfied).  
**Independent test**: `grep "linux-asan-ubsan" .github/workflows/ci-linux.yml` confirms the job exists and runs `dw-tests.sh` under `-fsanitize=address,undefined`.

- [X] T001 Verify `linux-asan-ubsan` job exists in `.github/workflows/ci-linux.yml`, runs `dw-tests.sh` under `-fsanitize=address,undefined`, and fails the build on any sanitizer error — confirming FR-007 is fully met with no CI changes required

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: All artifacts that MUST be committed before badge registration (US1) can proceed: SPDX identifiers in every source file, `SECURITY.md`, and `CONTRIBUTING.md`. All SPDX tasks are mutually independent and may run in parallel.

**⚠️ CRITICAL**: Badge registration (T021) cannot proceed until T002–T019 are all committed to the branch.

### SPDX License Identifiers (FR-004)

Add `SPDX-License-Identifier: GPL-3.0-or-later` inside each file's existing `/* ... */` license header, on a new line immediately before the closing ` * **...***/` delimiter, formatted as ` *   SPDX-License-Identifier: GPL-3.0-or-later`. Do not alter any other line in the header.

- [X] T002 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw.h` inside the existing license header block
- [X] T003 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_blas.c` inside the existing license header block
- [X] T004 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_blas.h` inside the existing license header block
- [X] T005 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_globals.c` inside the existing license header block
- [X] T006 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_main.c` inside the existing license header block
- [X] T007 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_phases.c` inside the existing license header block
- [X] T008 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_phases.h` inside the existing license header block
- [X] T009 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_rounding.c` inside the existing license header block
- [X] T010 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_rounding.h` inside the existing license header block
- [X] T011 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_solver.c` inside the existing license header block
- [X] T012 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_solver.h` inside the existing license header block
- [X] T013 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_subprob.c` inside the existing license header block
- [X] T014 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_subprob.h` inside the existing license header block
- [X] T015 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_support.c` inside the existing license header block
- [X] T016 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `src/dw_support.h` inside the existing license header block
- [X] T017 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `tests/test_blas.c` inside the existing license header block
- [X] T018 [P] Add SPDX-License-Identifier: GPL-3.0-or-later to `tests/test_lib_api.c` inside the existing license header block

### Security Policy (FR-001, FR-012)

- [X] T019 Create `SECURITY.md` at repository root with all six required sections per `contracts/SECURITY-contract.md`: (1) Supported Versions table showing 2.x supported / 1.x unsupported, (2) Reporting — direct to `https://github.com/alotau/dwsolver/security/advisories/new` with email fallback; explicit "do NOT open a public issue" warning, (3) Response Timeline — 14 business day acknowledgment, 60-day resolution target, (4) Scope — memory safety, integer overflow via crafted LP input, data races, GLPK dependency vulns; out of scope: DoS via pathological inputs, build-time tool vulns, (5) Security Design Notes — no crypto (N/A to all crypto OpenSSF criteria), LP file input only, pthreads mutex protection, `dw_solve()` not thread-safe across concurrent calls on same context, no network access, (6) Disclosure Policy — private advisory, coordinated disclosure, CVE via GitHub Advisory

### Contribution Guidelines (FR-002)

- [X] T020 [P] Create `CONTRIBUTING.md` at repository root with all seven required sections per `contracts/CONTRIBUTING-contract.md`: (1) Welcome paragraph, (2) Bug Reports — GitHub Issues at `https://github.com/alotau/dwsolver/issues`; required info: OS, GLPK version, build config, input file, full output, (3) Enhancement Requests — GitHub Issues with `enhancement` label; describe use case first, (4) Pull Requests — fork → branch (`fix/`, `feat/`, or `chore/` prefix) → PR against `main`; all CI must pass, (5) Coding Standards — C11; SEI CERT C (spec-008); ISO/IEC TS 17961 (spec-013); POSIX.1-2008; static analysis clean, (6) DCO Sign-Off — `git commit -s` required; include DCO 1.1 text or link to `https://developercertificate.org/`; explain what signing certifies, (7) CI Requirements — all five CI workflows must pass before merge (macOS, Linux matrix, Windows, Docker, Compliance)

**Checkpoint**: T002–T020 committed → badge registration (T021) can proceed.

---

## Phase 3: User Story 1 — Maintainer Earns Passing Badge (Priority: P1) 🎯 MVP

**Goal**: All prerequisite gap items committed; project registered on the OpenSSF platform; self-assessment completed with 100% mandatory criteria Met; badge status shows "Passing".

**Independent Test**: Visit the badge assessment URL and confirm badge status reads "Passing" with no "Unmet" mandatory criteria.

### Release Notes Infrastructure (FR-005, FR-006)

- [X] T021 [P] [US1] Create `CHANGELOG.md` at repository root using Keep a Changelog 1.1.0 format (`https://keepachangelog.com`): top-level `## [Unreleased]` section (empty), then `## [2.0.0] - 2026-xx-xx` entry with subsections `### Added` (summary of all work from specs 001–015 in bullet form), `### Changed` (N/A), `### Fixed` (N/A), `### Deprecated` (N/A), `### Removed` (N/A), `### Security` ("Initial public release. No known vulnerabilities."); add preamble explaining the format and that every release must include a Security subsection
- [X] T022 [P] [US1] Update `.github/workflows/release.yml` so the GitHub Release body template includes a `## Security` section (default text: "No security changes in this release.") and a link to `CHANGELOG.md`; this establishes the release-notes format required by the `release_notes_vulns` criterion (FR-006)

### README Badge Placeholder (FR-003, FR-011)

- [X] T023 [US1] Add OpenSSF Passing badge placeholder to `README.md` in the existing CI badge row (after the Docker badge line) using the format: `[![OpenSSF Best Practices](https://bestpractices.coreinfrastructure.org/projects/NNNNN/badge)](https://bestpractices.coreinfrastructure.org/projects/NNNNN)` with `NNNNN` as a literal placeholder; add a comment `<!-- Replace NNNNN with real project ID from badge registration (task T024) -->`

### Badge Registration (FR-009, FR-010, FR-011) — Manual Action Required

- [ ] T024 [US1] Register the project on the OpenSSF Best Practices Badge platform: (a) visit `https://bestpractices.coreinfrastructure.org/en/projects/new`, (b) authenticate with the GitHub account that owns `alotau/dwsolver`, (c) enter repository URL `https://github.com/alotau/dwsolver`, (d) confirm auto-detected fields (license: GPL-3.0, language: C), (e) record the assigned Project ID integer from the resulting URL `.../projects/NNNNN`

- [ ] T025 [US1] Complete the OpenSSF self-assessment using the evidence table in `specs/015-openssf-passing-badge/quickstart.md`: answer every mandatory Passing criterion with "Met" or "N/A" (crypto section); for each "Met" response provide the evidence URL or file path per the quickstart evidence table; confirm zero criteria remain "Unmet"; submit the assessment

### README Finalization (FR-003, FR-011)

- [ ] T026 [US1] Replace `NNNNN` placeholder with the real Project ID from T024 in `README.md` (both the badge image URL and the badge link URL); add a prose line in the README Compliance section: "OpenSSF Best Practices assessment: `https://bestpractices.coreinfrastructure.org/projects/<ID>`"

**Checkpoint**: At this point US1 is complete — the badge platform shows "Passing" and the README displays the live badge widget.

---

## Phase 4: User Story 2 — Downstream Consumer Verifies Project Quality (Priority: P2)

**Goal**: The Passing badge widget is visible in the README and the assessment page is publicly accessible and shows "Passing".

**Independent Test**: Visit the repository README on GitHub.com; confirm the OpenSSF badge image renders in the badge row; click it and confirm the assessment page loads with "Passing" status.

- [ ] T027 [P] [US2] Verify the OpenSSF badge image renders correctly in `README.md` on GitHub.com — navigate to `https://github.com/alotau/dwsolver`; confirm the badge appears in the CI badge row and is not a broken image
- [ ] T028 [P] [US2] Verify the badge link navigates to the live assessment page at `https://bestpractices.coreinfrastructure.org/projects/<ID>` and the page shows badge status "Passing" (green) with the project name "dwsolver"

**Checkpoint**: US2 complete — downstream consumers can verify project quality via the badge.

---

## Phase 5: User Story 3 — New Contributor Knows How to Engage (Priority: P3)

**Goal**: `CONTRIBUTING.md` and `SECURITY.md` are present, contain all required sections, are discoverable from the README, and surfaced by the GitHub UI on new issue/PR creation.

**Independent Test**: Confirm `CONTRIBUTING.md` and `SECURITY.md` exist at the repository root; open a new issue draft on GitHub and confirm the contributing guide link appears; verify the SECURITY.md advisory link resolves.

- [X] T029 [US3] Validate `SECURITY.md` content against `contracts/SECURITY-contract.md`: confirm all six sections are present; visit `https://github.com/alotau/dwsolver/security/advisories/new` (while authenticated) to confirm the GitHub Security Advisory form is accessible; confirm `grep -c "14" SECURITY.md` returns ≥ 1 (14-day acknowledgment documented)
- [X] T030 [US3] Validate `CONTRIBUTING.md` content against `contracts/CONTRIBUTING-contract.md`: confirm all seven sections are present; confirm `grep -c "Signed-off-by\|DCO" CONTRIBUTING.md` returns ≥ 2; confirm `grep -c "GitHub Issues" CONTRIBUTING.md` returns ≥ 2
- [ ] T031 [US3] Verify GitHub surfaces `CONTRIBUTING.md` automatically: navigate to `https://github.com/alotau/dwsolver/issues/new` and confirm GitHub shows the contributing guide link in the issue creation UI (GitHub does this automatically when CONTRIBUTING.md exists at the root)

**Checkpoint**: US3 complete — new contributors and security researchers have clear, discoverable onboarding paths.

---

## Final Phase: Polish & Cross-Cutting Validation

**Purpose**: Zero-defect validation across all produced artifacts before the PR is opened against `main`.

- [X] T032 [P] Validate SPDX coverage: run `grep -rL "SPDX-License-Identifier" src/ tests/test_blas.c tests/test_lib_api.c` from the repository root; confirm the command returns no output (zero files missing the identifier)
- [X] T033 [P] Validate no CI regressions: run `make check` locally and confirm the test suite passes; confirm all GitHub Actions CI jobs are green on the feature branch before opening the PR

---

## Dependencies

```
T001 (verify ASAN CI)
  └─ (no code change; informational only)

T002–T018 (SPDX headers) — all parallel, no cross-dependencies
T019 (SECURITY.md) — parallel with T020
T020 (CONTRIBUTING.md) — parallel with T019

T021 (CHANGELOG.md) — parallel with T022, T023; no dependency on T002–T020
T022 (release.yml) — parallel with T021, T023
T023 (README placeholder) — parallel with T021, T022; needed before T024

T024 (badge registration) ← requires T002–T023 committed
T025 (self-assessment) ← requires T024
T026 (README real ID) ← requires T024

T027, T028 (US2 verification) ← requires T025, T026
T029, T030, T031 (US3 verification) ← requires T019, T020 committed

T032 (SPDX validation) ← requires T002–T018
T033 (CI validation) ← requires all prior tasks
```

## Parallel Execution Examples

### Batch A (no sequential deps — start immediately)
Run together: T001, T002–T018 (all SPDX files), T019 (SECURITY.md), T020 (CONTRIBUTING.md)

### Batch B (after Batch A is committed)
Run together: T021 (CHANGELOG.md), T022 (release.yml), T023 (README placeholder)

### Batch C (after Batch B is committed) — manual
T024 → T025 → T026 (sequential: register, assess, finalize README)

### Batch D (after Batch C)
Run together: T027, T028 (US2), T029, T030, T031 (US3), T032 (SPDX grep)

### Final
T033 (CI check before PR open)

## Implementation Strategy

**MVP**: User Story 1 (Passing badge earned). Delivers SC-001 and SC-002 immediately.

**Suggested delivery order**:
1. Commit Batch A (17 SPDX edits + SECURITY.md + CONTRIBUTING.md) — can be done in one commit pass
2. Commit Batch B (CHANGELOG.md + release.yml + README placeholder) — one commit
3. Manual registration (T024) + self-assessment (T025) — requires web browser
4. Commit T026 (real project ID in README) — one-line edit
5. Validate (T027–T033) — verify and confirm
6. Open PR against `main`

**Note on T024–T025 (manual web steps)**: These cannot be automated — they require authenticated browser login to the OpenSSF platform. Complete all code commits first so the self-assessment evidence links resolve to committed files.
