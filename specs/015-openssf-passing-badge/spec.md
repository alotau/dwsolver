# Feature Specification: OpenSSF Best Practices Badge (Passing Level)

**Feature Branch**: `015-openssf-passing-badge`  
**Created**: 2026-03-24  
**Status**: Draft  
**Input**: Implement OpenSSF Best Practices Badge (CII Badge Program, Linux Foundation) — Passing level

## Background

The [OpenSSF Best Practices Badge Program](https://bestpractices.coreinfrastructure.org/) (formerly CII Best Practices) is a self-assessed, tiered badge program (Passing → Silver → Gold) administered by the Linux Foundation's Open Source Security Foundation. It evaluates open-source projects across security, testing, CI, documentation, and supply-chain hygiene.

Achieving the **Passing** badge signals to downstream consumers of `libdwsolver` that the project follows recognized open-source quality and security practices. The project already satisfies many Passing criteria (GPL-3.0 license, CI matrix, Docker build, test suite, static/dynamic analysis, distributed VCS). This feature closes the remaining gaps.

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Maintainer Earns Passing Badge (Priority: P1)

The project maintainer registers dwsolver on the OpenSSF Best Practices Badge platform, fills in the self-assessment, links evidence for each criterion, and achieves the Passing badge.

**Why this priority**: This is the primary deliverable. All other stories depend on it. A visible badge immediately communicates project health to anyone landing on the repository.

**Independent Test**: Can be fully tested by visiting the badge assessment URL and confirming all Passing criteria are marked "Met" with no "Unmet" entries.

**Acceptance Scenarios**:

1. **Given** the project is registered on the OpenSSF badge platform, **When** the maintainer submits the self-assessment, **Then** every mandatory Passing criterion is marked "Met" and the badge status shows "Passing".
2. **Given** the badge assessment is complete, **When** the README is viewed, **Then** the Passing badge image is visible and links to the live assessment page.
3. **Given** the badge is awarded, **When** a new release tag is pushed, **Then** the release workflow completes without badge-related CI failures.

---

### User Story 2 — Downstream Consumer Verifies Project Quality (Priority: P2)

An engineer evaluating `libdwsolver` as a dependency visits the repository, sees the Passing badge, and can follow the badge link to verify the security and quality practices in detail before committing to adoption.

**Why this priority**: This is the primary value for external users — the badge provides an independently verifiable signal of project credibility without requiring deep code review.

**Independent Test**: Can be fully tested by visiting the repository README and confirming the badge is present and links to a "Passing" assessment page showing all criteria met.

**Acceptance Scenarios**:

1. **Given** the README is displayed on the repository homepage, **When** a visitor views it, **Then** they see the OpenSSF Passing badge image within the badge row at the top of the document.
2. **Given** the visitor clicks the badge image, **When** the browser follows the link, **Then** they land on the OpenSSF assessment page confirming "Passing" status for this project.
3. **Given** the assessment page is loaded, **When** the visitor reviews the security section, **Then** they find a documented vulnerability reporting process with a private channel and a 14-day response commitment.

---

### User Story 3 — New Contributor Knows How to Engage (Priority: P3)

A developer who wants to submit a patch or report a bug finds clear, documented processes for both, without needing to search issue history or ask in forums.

**Why this priority**: CONTRIBUTING.md and SECURITY.md are required Passing criteria and also fulfill the OpenSSF requirement to document contribution and vulnerability-reporting processes. They benefit contributors independently of the badge.

**Independent Test**: Can be fully tested by confirming CONTRIBUTING.md and SECURITY.md exist at the repository root with required content, and that both are referenced or discoverable from the README.

**Acceptance Scenarios**:

1. **Given** a contributor wants to submit a patch, **When** they open CONTRIBUTING.md, **Then** they find: how to open an issue or PR, the coding standards reference, and the Developer Certificate of Origin (DCO) sign-off requirement.
2. **Given** a security researcher discovers a vulnerability, **When** they open SECURITY.md, **Then** they find a private reporting channel (email or GitHub private advisory) and a documented policy stating reports will be acknowledged within 14 days.
3. **Given** SECURITY.md exists, **When** the project receives a vulnerability report, **Then** the maintainer acknowledges it within 14 days and tracks resolution through a private advisory prior to public disclosure.

---

### Edge Cases

- What happens if a Passing criterion is disputed by a badge reviewer? The self-assessment allows free-text justification; the maintainer must provide links to verifiable evidence (CI logs, file paths, commit SHAs).
- What if a vulnerability is reported and no maintainer is reachable within 14 days? SECURITY.md must document a secondary contact or escalation path so the 14-day commitment can still be honored.
- What if a new release omits the security-advisory section in release notes? The release notes template must include a "Security" section (even if it reads "No security changes in this release") to satisfy the `release_notes_vulns` criterion on an ongoing basis.

## Requirements *(mandatory)*

### Functional Requirements

#### Documentation Gaps

- **FR-001**: The project MUST have a `SECURITY.md` file at the repository root documenting: (a) the private vulnerability reporting channel, (b) the scope of what constitutes a security issue, and (c) a commitment to acknowledge reports within 14 days.
- **FR-002**: The project MUST have a `CONTRIBUTING.md` file at the repository root documenting: (a) how to open bug reports and enhancement requests via the issue tracker, (b) how to submit patches via pull requests, (c) coding standards (reference to existing compliance specs), and (d) the Developer Certificate of Origin (DCO) sign-off requirement.
- **FR-003**: The README MUST include the OpenSSF Passing badge image in the CI badge row, linking to the live assessment page on bestpractices.coreinfrastructure.org.
- **FR-004**: Source files MUST carry SPDX license identifier comments (`SPDX-License-Identifier: GPL-3.0-or-later`) so automated license scanners can verify the project's FLOSS license without parsing COPYING manually.

#### Release Notes

- **FR-005**: The `ChangeLog` (or a new `CHANGELOG.md`) MUST be updated to use a consistent per-release format that includes a "Security" subsection. Releases with no security changes MUST state "No security changes" rather than omitting the section.
- **FR-006**: Each tagged release MUST have a corresponding GitHub Release entry whose body includes: a summary of changes, a "Security" section, and a link to the full ChangeLog.

#### Continuous Integration

- **FR-007**: The dynamic analysis (AddressSanitizer) build MUST be executed as a CI job (not only as a local artifact), and the job MUST fail the build on any AddressSanitizer error. This documents the `dynamic_analysis` and `dynamic_analysis_unsafe` Passing criteria.
- **FR-008**: The existing static-analysis CI job MUST be explicitly cross-referenced in the OpenSSF badge assessment as evidence for the `static_analysis` and `static_analysis_common_vulnerabilities` criteria.

#### Badge Registration & Self-Assessment

- **FR-009**: The project MUST be registered on the OpenSSF Best Practices Badge platform (https://bestpractices.coreinfrastructure.org) with the correct repository URL.
- **FR-010**: All mandatory Passing-level criteria MUST be answered and evidence links provided for each "Met" response. No criterion may remain "Unmet" at submission.
- **FR-011**: The badge assessment URL MUST be stored in the repository (e.g., in README or a `.openssf` config file) so it remains discoverable if the badge image ever breaks.

#### Security Documentation

- **FR-012**: `SECURITY.md` MUST include a brief security design overview covering the main risk surface: numerical precision and integer overflow in LP computations, thread-safety guarantees of the public API, and absence of cryptographic operations (not applicable / N/A to crypto criteria).

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The project achieves "Passing" status on the OpenSSF Best Practices Badge platform, with 100% of mandatory Passing criteria marked "Met".
- **SC-002**: The Passing badge widget is visible in the README within the existing CI badge row, and clicking it loads the live assessment page confirming "Passing" status.
- **SC-003**: `SECURITY.md` and `CONTRIBUTING.md` are present at the repository root and contain all required content sections before the branch is merged to main.
- **SC-004**: All source files under `src/` carry an SPDX license identifier, verifiable by running a license-scanning tool with zero "no-license" findings.
- **SC-005**: A CI job for dynamic analysis (AddressSanitizer) runs on every push to main and pull request, completing with zero errors against the full test suite.
- **SC-006**: At least one GitHub Release entry exists with the structured format (summary + Security section), establishing the template for future releases.

## Assumptions

- The project's GitHub repository is publicly accessible; all OpenSSF badge criteria requiring a public repo are already satisfied.
- The GPL-3.0 license is OSI-approved and satisfies the `floss_license` and `floss_license_osi` Passing criteria without any changes.
- Semantic versioning (MAJOR.MINOR.PATCH, currently `2.0.0`) is already in use via `configure.ac` and `v*` release tags; no version scheme changes are required.
- GitHub Issues serves as the public bug and enhancement tracker; the `report_tracker` criterion is met once CONTRIBUTING.md documents this.
- The existing clang-tidy + cppcheck CI job satisfies `static_analysis` and `static_analysis_common_vulnerabilities`; only documentation and cross-referencing are needed.
- The project does not perform any cryptographic operations; all OpenSSF crypto criteria will be answered "Not applicable" in the self-assessment.
- "Passing" is the target badge level. Silver and Gold are out of scope for this feature.

## Dependencies

- Spec 003 (GitHub Actions CI) — AddressSanitizer CI job will be added to the existing CI matrix.
- Spec 008 (SEI CERT C Compliance) — Static analysis CI job provides evidence for `static_analysis` criteria.
- Spec 013 (ISO/IEC TS 17961 Compliance) — `ci-compliance.yml` provides additional static-analysis evidence.
- Spec 012 (Release Infrastructure) — Release workflow and tagging practices provide evidence for `version_tags` and `version_semver` criteria.
