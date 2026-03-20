# Feature Specification: GitHub Actions CI — Multi-Platform Build & Test

**Feature Branch**: `003-github-actions-ci`
**Created**: 2026-03-19
**Status**: Draft
**Input**: CI must run `tests/dw-tests.sh` on every pull request and every push to `main`, across macOS, Linux, and Windows.

---

## Summary

DWSOLVER currently has no automated CI. The constitution (Principle III) mandates that
"CI MUST validate all supported platforms before a branch is eligible to merge into main."
This spec defines a GitHub Actions workflow that builds the solver and runs the test
suite on macOS, Linux, and Windows for every pull request targeting `main` and every
direct push to `main`.

Windows presents unique build challenges: it requires a POSIX-compatible build
environment (MinGW-w64 or Cygwin) and a pthreads compatibility layer. The CI workflow
must handle the platform-specific dependencies transparently.

---

## User Scenarios & Testing

### User Story 1 — PRs are automatically gated on macOS and Linux (Priority: P1)

A developer opens a pull request against `main`. GitHub Actions automatically triggers
a build-and-test run on both macOS (Clang) and Linux (GCC). If either platform's
tests fail, the PR is blocked. The developer sees per-platform pass/fail status
on the PR page within a reasonable time.

**Why this priority**: macOS and Linux are the two platforms with confirmed working
builds. Gating PRs on these immediately protects `main` from regressions and
actively enforces Principle III.

**Independent Test**: Open a PR against `main` (or push a branch), observe the
"Checks" tab — the macOS and Linux jobs appear, run, and report green for the
current working state of the code.

**Acceptance Scenarios**:

1. **Given** a PR is opened against `main`,
   **When** the CI workflow triggers,
   **Then** a build-and-test job runs on macOS (latest) with
   `./configure --enable-named-semaphores && make` and `tests/dw-tests.sh`,
   and reports a pass/fail status on the PR within 15 minutes.

2. **Given** a PR is opened against `main`,
   **When** the CI workflow triggers,
   **Then** a build-and-test job runs on Ubuntu (latest) with
   `./configure && make` and `tests/dw-tests.sh`,
   and reports a pass/fail status on the PR within 15 minutes.

3. **Given** a direct push to `main`,
   **When** the CI workflow triggers,
   **Then** both macOS and Linux jobs run and results are recorded.

4. **Given** a branch where `tests/dw-tests.sh` fails on Linux,
   **When** a PR from that branch is opened,
   **Then** the Linux CI job reports failure and the PR check shows a red ✗,
   preventing merge until resolved.

---

### User Story 2 — Windows build is attempted in CI (Priority: P2)

The CI workflow includes a Windows job using a MinGW-w64 / MSYS2 environment that
attempts to build the solver and run the test suite. Because Windows support is a
known work-in-progress (separate spec 004), the Windows job is allowed to fail
without blocking the PR — but its result is always visible so regressions can be
tracked over time.

**Why this priority**: Visibility is valuable even before full support is achieved.
The constitution requires Windows support; CI must at minimum attempt it and surface
the current status. Non-blocking until spec 004 is complete.

**Independent Test**: Open any PR — a "Windows" CI job appears and runs. If the
build or tests fail, the PR is still mergeable (job is marked as non-required),
but the failure is visible.

**Acceptance Scenarios**:

1. **Given** a PR is opened,
   **When** the CI workflow triggers,
   **Then** a Windows job runs using MSYS2 with MinGW-w64 toolchain, attempts
   `./configure && make`, captures the result, and reports it on the PR.

2. **Given** the Windows job fails due to known platform gaps,
   **When** the PR check results are displayed,
   **Then** the macOS and Linux required checks are green, and the Windows
   job is visible but not a merge blocker.

3. **Given** the Windows job eventually passes (after spec 004 is implemented),
   **When** the Windows job is promoted to a required check,
   **Then** PRs with Windows failures are blocked from merging.

---

### User Story 3 — CI is configured as a required status check on main (Priority: P1)

The GitHub repository is configured so that the macOS and Linux CI jobs are
*required status checks* on the `main` branch. No pull request can be merged to
`main` until both pass. This is enforced at the repository level (branch protection
rules), not just by convention.

**Why this priority**: Without branch protection, the CI workflow is advisory only.
Mandatory status checks are what give CI its enforcement power. This is the
highest-leverage protection for `main`.

**Independent Test**: Attempt to merge a PR with a failing macOS or Linux job —
the GitHub UI prevents the merge with "Required status checks have not passed."

**Acceptance Scenarios**:

1. **Given** the branch protection rule is active on `main`,
   **When** a PR has a failing macOS or Linux CI job,
   **Then** the "Merge pull request" button is disabled until the jobs pass.

2. **Given** all required CI jobs pass,
   **When** the PR author clicks "Merge pull request",
   **Then** the merge proceeds normally.

---

### Edge Cases

- If a CI job hangs (e.g., `four_sea` on a slow runner), it must time out and
  report failure rather than blocking indefinitely. Per-job timeout must be set.
- The Windows MSYS2 environment setup (toolchain install, autotools) must be
  cached to keep CI times reasonable. Cache must be invalidated when the
  MSYS2 package list changes.
- `tests/dw-tests.sh` uses `pushd`/`popd` (bash built-ins) — the shebang is
  `#!/bin/sh`. On Windows under MSYS2, `sh` must resolve to a bash-compatible
  shell.
- The CI workflow file must not embed secrets or tokens. The `GITHUB_TOKEN`
  provided automatically by Actions is sufficient for status reporting.
- The `four_sea` example is the most computationally intensive test (4 subproblem
  threads). CI runners must have at least 2 vCPUs for it to complete reliably.

---

## Requirements

### Functional Requirements

- **FR-001**: A GitHub Actions workflow file MUST exist at `.github/workflows/ci.yml`.
- **FR-002**: The workflow MUST trigger on `pull_request` events targeting `main`
  and on `push` events to `main`.
- **FR-003**: The workflow MUST include a macOS job that installs autotools, runs
  `./configure --enable-named-semaphores && make`, and runs `tests/dw-tests.sh`
  with `dwsolver` on `$PATH`.
- **FR-004**: The workflow MUST include a Linux (Ubuntu) job that installs autotools
  and GCC, runs `./configure && make`, and runs `tests/dw-tests.sh` with `dwsolver`
  on `$PATH`.
- **FR-005**: The workflow MUST include a Windows job using MSYS2 (MinGW-w64) that
  attempts `./configure && make` and `tests/dw-tests.sh`. The Windows job MUST be
  marked `continue-on-error: true` until Windows support is confirmed complete.
- **FR-006**: Each job MUST have a `timeout-minutes` limit (suggested: 30 minutes)
  to prevent runaway jobs.
- **FR-007**: The macOS and Linux jobs MUST be registered as required status checks
  on the `main` branch via GitHub branch protection rules.
- **FR-008**: The workflow MUST NOT require any secrets beyond the automatically
  provided `GITHUB_TOKEN`.
- **FR-009**: The MSYS2 package cache MUST be used via `actions/cache` or the
  `msys2/setup-msys2` action's built-in caching to keep Windows job setup time
  under 5 minutes.

### Key Entities

- **`.github/workflows/ci.yml`**: The workflow definition file. Single workflow,
  three jobs (`build-macos`, `build-linux`, `build-windows`), triggered by
  `pull_request` and `push` to `main`.
- **Branch protection rule on `main`**: Requires `build-macos` and `build-linux`
  to pass before merge. Configured in GitHub repository settings (not in code).
- **`tests/dw-tests.sh`**: The existing test script, unmodified by this spec.
  Extended by spec 002 (KD-009 fix) to 6 tests; all 6 must pass in CI.

### Out of Scope

- Self-hosted runners — use GitHub-hosted runners only.
- Code coverage reporting.
- Artifact upload (binaries, logs) beyond CI job logs.
- Performance benchmarking jobs.
- Windows branch protection (deferred until spec 004 completes Windows support).

---

## Success Criteria

### Measurable Outcomes

- **SC-001**: Every PR against `main` triggers CI within 60 seconds of opening
  and completes macOS + Linux jobs within 15 minutes in the normal case.
- **SC-002**: A PR that breaks `tests/dw-tests.sh` on either macOS or Linux
  cannot be merged to `main` without resolving the failure — enforced by branch
  protection, not convention.
- **SC-003**: The Windows CI job runs on every PR and its result is visible,
  even if non-blocking.
- **SC-004**: Zero CI secrets required beyond `GITHUB_TOKEN` (automatic).
- **SC-005**: The CI configuration can be understood and modified by a developer
  familiar with GitHub Actions without reading additional documentation.

### Assumptions

- The GitHub repository at `github.com/alotau/dwsolver` allows GitHub Actions
  workflows and has at least one maintainer with permission to set branch
  protection rules.
- GitHub-hosted `macos-latest` runners provide Homebrew; `ubuntu-latest` provides
  `apt-get`; `windows-latest` provides MSYS2 via `msys2/setup-msys2`.
- The `tests/dw-tests.sh` script works correctly after the spec 002 repairs land.
  CI for this spec assumes spec 002 is merged first (or this branch is based on it).

## User Scenarios & Testing *(mandatory)*

<!--
  IMPORTANT: User stories should be PRIORITIZED as user journeys ordered by importance.
  Each user story/journey must be INDEPENDENTLY TESTABLE - meaning if you implement just ONE of them,
  you should still have a viable MVP (Minimum Viable Product) that delivers value.
  
  Assign priorities (P1, P2, P3, etc.) to each story, where P1 is the most critical.
  Think of each story as a standalone slice of functionality that can be:
  - Developed independently
  - Tested independently
  - Deployed independently
  - Demonstrated to users independently
-->

### User Story 1 - [Brief Title] (Priority: P1)

[Describe this user journey in plain language]

**Why this priority**: [Explain the value and why it has this priority level]

**Independent Test**: [Describe how this can be tested independently - e.g., "Can be fully tested by [specific action] and delivers [specific value]"]

**Acceptance Scenarios**:

1. **Given** [initial state], **When** [action], **Then** [expected outcome]
2. **Given** [initial state], **When** [action], **Then** [expected outcome]

---

### User Story 2 - [Brief Title] (Priority: P2)

[Describe this user journey in plain language]

**Why this priority**: [Explain the value and why it has this priority level]

**Independent Test**: [Describe how this can be tested independently]

**Acceptance Scenarios**:

1. **Given** [initial state], **When** [action], **Then** [expected outcome]

---

### User Story 3 - [Brief Title] (Priority: P3)

[Describe this user journey in plain language]

**Why this priority**: [Explain the value and why it has this priority level]

**Independent Test**: [Describe how this can be tested independently]

**Acceptance Scenarios**:

1. **Given** [initial state], **When** [action], **Then** [expected outcome]

---

[Add more user stories as needed, each with an assigned priority]

### Edge Cases

<!--
  ACTION REQUIRED: The content in this section represents placeholders.
  Fill them out with the right edge cases.
-->

- What happens when [boundary condition]?
- How does system handle [error scenario]?

## Requirements *(mandatory)*

<!--
  ACTION REQUIRED: The content in this section represents placeholders.
  Fill them out with the right functional requirements.
-->

### Functional Requirements

- **FR-001**: System MUST [specific capability, e.g., "allow users to create accounts"]
- **FR-002**: System MUST [specific capability, e.g., "validate email addresses"]  
- **FR-003**: Users MUST be able to [key interaction, e.g., "reset their password"]
- **FR-004**: System MUST [data requirement, e.g., "persist user preferences"]
- **FR-005**: System MUST [behavior, e.g., "log all security events"]

*Example of marking unclear requirements:*

- **FR-006**: System MUST authenticate users via [NEEDS CLARIFICATION: auth method not specified - email/password, SSO, OAuth?]
- **FR-007**: System MUST retain user data for [NEEDS CLARIFICATION: retention period not specified]

### Key Entities *(include if feature involves data)*

- **[Entity 1]**: [What it represents, key attributes without implementation]
- **[Entity 2]**: [What it represents, relationships to other entities]

## Success Criteria *(mandatory)*

<!--
  ACTION REQUIRED: Define measurable success criteria.
  These must be technology-agnostic and measurable.
-->

### Measurable Outcomes

- **SC-001**: [Measurable metric, e.g., "Users can complete account creation in under 2 minutes"]
- **SC-002**: [Measurable metric, e.g., "System handles 1000 concurrent users without degradation"]
- **SC-003**: [User satisfaction metric, e.g., "90% of users successfully complete primary task on first attempt"]
- **SC-004**: [Business metric, e.g., "Reduce support tickets related to [X] by 50%"]
