# Feature Specification: ISO/IEC TS 17961:2013 Compliance

**Feature Branch**: `014-iso-ts-17961-compliance`
**Created**: 2026-03-23
**Status**: Draft
**Input**: Add a specification to bring dwsolver's authored source files (`src/dw_*.c`, `src/dw_*.h`) into compliance with ISO/IEC TS 17961:2013 ("C Secure Coding Rules"). This is the ISO formalization of the secure C coding rules that SEI CERT C expresses as guidelines. Tool vendors (Coverity, Klocwork, CodeSonar) certify against TS 17961 rather than the wiki. The project already has SEI CERT C compliance (spec 008) and strict ISO C11 / POSIX.1-2008 compliance (spec 013). This spec leverages that prior work, audits the source against the 22 TS 17961 rules, documents any gaps not already covered by spec 008, and produces an auditable compliance artifact suitable for formal procurement contexts.

## Background

ISO/IEC TS 17961:2013 ("C Secure Coding Rules") is the ISO Technical Specification that formalizes secure C coding requirements as **normative rules** with defined severity levels and diagnostic requirements. Commercial static analysis vendors — Coverity, Klocwork, CodeSonar, and LDRA — certify their tools against TS 17961 rather than the SEI CERT C wiki, making TS 17961 the standard referenced in defense, aerospace, automotive (ISO 26262), and medical device (IEC 62443) procurement requirements.

The specification defines **22 analyzable rules** covering:
- Initialization (intptrconv, intobjptr, xplicitcomp, boolasgn)
- Integer safety (trstcmp, wraparound, argcomp, intnegator, ioilecc, ga-buffer)
- String/memory operations (strmod, memiograph, accfree, dblfree, nonnullptr, nullref)
- Arithmetic and I/O (diverr, exceptbits, liberr)
- Concurrency (datarace, asyncsig, toinit)

**Relationship to prior work**:
- **Spec 008 (SEI CERT C)**: Addressed the primary threading hazards (POS51/52/54-C), memory checks (MEM32-C), file handle errors (FIO01/42-C), and integer type issues (INT30/31-C). The SCI CERT C rules have approximate correspondences to TS 17961 rules, but the mapping is not 1-to-1; TS 17961 includes rules with no CERT equivalent and vice versa.
- **Spec 013 (Strict C11/POSIX.1-2008)**: Enforced `-std=c11 -pedantic-errors`, set `_POSIX_C_SOURCE=200809L`, and eliminated GNU extensions. These address the TS 17961 conformance baseline (rules require conforming C programs as their starting point).

This spec's job is to **close the gap**: identify which of the 22 TS 17961 rules are already satisfied by spec 008/013 work, identify any remaining gaps, remediate those gaps, and produce a formal compliance matrix.

**Scope**: `src/dw_*.c` and `src/dw_*.h` only. Vendored GLPK (`third-party/glpk/`) is excluded — it is a third-party library boundary.

---

## Clarifications

### Session 2026-03-23

- Q: For TS 17961 rules that have no automated checker in `cppcheck --addon=cert` or `clang-tidy cert-*`, what form of manual audit evidence is required in the compliance matrix? → A: Grep-based absence proof only — a command and its zero-result output showing the triggering construct is absent from `src/dw_*.c`
- Q: Which CI platform/workflow file should the TS 17961 enforcement job be added to? → A: New dedicated workflow `ci-compliance.yml` (runs on Linux only, separate from build/test matrix)
- Q: What is the expected distribution of verdicts in the initial audit (PASS via prior work / N-A / FAIL requiring remediation)? → A: Majority PASS via spec 008/013 prior work; a small number of rules require fresh verification; N-A for signal/setjmp rules whose triggering constructs are absent from the codebase
- Q: How should the synthetic-violation CI enforcement test (SC-003) be implemented? → A: A shell script in `tests/` that injects a known-bad construct, runs `cppcheck --addon=cert`, asserts non-zero exit, then reverts — executed as a dedicated step in `ci-compliance.yml`
- Q: What is the canonical term for the 22-rule evidence table artifact — "compliance matrix" or "audit matrix"? → A: Compliance matrix — used everywhere; "audit matrix" retired (formerly used interchangeably in draft)

---

## User Scenarios & Testing

### User Story 1 — Full compliance matrix is produced (Priority: P1)

A procurement engineer or auditor evaluating dwsolver for use in a regulated environment can read a single document that maps all 22 ISO/IEC TS 17961:2013 rules to the dwsolver codebase, records which rules are satisfied by prior work (specs 008 and 013), documents any newly remediated gaps, and provides evidence for each verdict.

**Why this priority**: The compliance matrix is the primary deliverable of this feature. All other stories either feed into it or depend on it. Without it, no formal compliance claim can be made regardless of how clean the code is. This is the artifact that procurement and certification processes actually require.

**Independent Test**: Can be fully tested in isolation by inspecting `specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md`. The document must exist, list all 22 rules, provide a PASS/FAIL/N-A verdict for each with a citation to evidence (spec 008/013 acceptance reports, static analysis findings, or grep-based absence proofs for rules with no automated checker), and be free of open "TBD" entries. No code changes are required to deliver this artifact for rules already covered.

**Acceptance Scenarios**:

1. **Given** the compliance matrix, **when** it is reviewed by a reader unfamiliar with the codebase, **then** every one of the 22 TS 17961 rules appears with: (a) a plain-language description, (b) a verdict (PASS / FAIL / N-A), (c) a rationale citing the specific prior work or new evidence, and (d) the approximate SEI CERT C rule correspondence where one exists.
2. **Given** a rule marked PASS by reference to spec 008 or spec 013, **when** the cited acceptance report or task is consulted, **then** the evidence described in the compliance matrix is actually present and accurate.
3. **Given** a rule marked N-A, **when** the rationale is read, **then** it identifies why the rule cannot be violated by this codebase (e.g., the code never performs the triggering operation).
4. **Given** a rule marked FAIL, **when** the compliance matrix is read, **then** it references an open task or FR in this spec that addresses the gap.

---

### User Story 2 — Gap rules are remediated (Priority: P2)

A developer implementing this feature can identify which subset of the 22 TS 17961 rules are not already satisfied by spec 008/013 work, and apply targeted fixes so that the codebase achieves clean verdicts for those rules.

**Why this priority**: The audit in Story 1 may reveal gaps not captured by the CERT-to-TS-17961 mapping. The expected scope is small: the majority of the 22 rules are anticipated to resolve as PASS via spec 008/013 prior work, with N-A verdicts for rules whose triggering constructs (e.g., `signal()`, `setjmp`/`longjmp`) are absent from this codebase. Story 2 targets only the residual rules requiring fresh verification or new code changes. It depends on Story 1 (you must audit before you fix) but is independently testable once the gap list is known.

**Independent Test**: Can be fully tested in isolation by: (a) identifying the gap rules from the compliance matrix, (b) applying fixes confined to `src/dw_*.c` / `src/dw_*.h`, (c) re-running the static analysis suite, and (d) verifying the full test suite still passes. Delivers provably gap-free code for the identified rules.

**Acceptance Scenarios**:

1. **Given** the gap rules identified in Story 1, **when** each gap rule's specific triggering construct is searched for in `src/dw_*.c` after remediation, **then** no triggering construct remains.
2. **Given** the remediated source, **when** `cppcheck --addon=cert` and `clang-tidy -checks=cert-*` are run against `src/dw_*.c`, **then** no diagnostics map to a TS 17961 gap rule (some CERT rule names do not map to TS 17961; only TS 17961-mapped diagnostics are required to be clean).
3. **Given** the remediated source, **when** the full test suite (`tests/dw-tests.sh`) is run, **then** all tests pass with no regression.
4. **Given** the remediated source, **when** it is compiled with `-std=c11 -pedantic-errors -Wall -Wextra`, **then** zero new diagnostics are introduced relative to the spec 013 baseline.

---

### User Story 3 — Static analysis tooling is configured and baselining is documented (Priority: P2)

A developer or CI administrator can run the exact tool commands listed in this specification to reproduce the compliance verification results, and the CI pipeline enforces that no new TS 17961 violations are introduced in subsequent commits.

**Why this priority**: Without reproducible tooling instructions and CI enforcement, the compliance audit is a point-in-time artifact that degrades immediately as the codebase evolves. This story makes compliance continuous.

**Independent Test**: Can be fully tested in isolation by following the tool commands in `specs/014-iso-ts-17961-compliance/audit/tooling-guide.md` on a clean checkout and confirming the outputs match the documented findings. Independently delivers a reproducible compliance-checking workflow.

**Acceptance Scenarios**:

1. **Given** the tooling guide, **when** a developer with `cppcheck` (≥ 2.12) and `clang-tidy` (≥ 17) installed follows it verbatim on a clean checkout of this branch, **then** the tool outputs match the documented baseline findings.
2. **Given** the CI pipeline, **when** a pull request introduces a new violation of a TS 17961-mapped `clang-tidy` or `cppcheck --addon=cert` rule in `src/dw_*.c`, **then** the `.github/workflows/ci-compliance.yml` job (Linux runner) fails and the violation is surfaced in the job output.
3. **Given** the tooling guide, **when** the optional Coverity or CodeSonar invocation instructions are followed on a system where those tools are licensed, **then** the checker categories to enable, the exclusion pattern for vendored GLPK, and the expected certificate-relevant findings are documented (even if those tools are not part of the automated CI).
4. **Given** `tests/test_ts17961_enforcement.sh`, **when** it is run in the `ci-compliance.yml` pipeline, **then** it injects a known TS 17961-violating construct into a temporary copy of a source file, confirms `cppcheck --addon=cert` exits non-zero and reports the expected finding, then reverts the file — proving the enforcement gate is functional, not merely present.

---

### User Story 4 — Compliance claim is suitable for formal procurement (Priority: P3)

A project maintainer responding to a procurement request that cites ISO/IEC TS 17961:2013 can provide the compliance matrix, tooling guide, and acceptance report as a self-contained compliance package without requiring additional interpretation by the evaluator.

**Why this priority**: This is a presentation and packaging concern — the underlying content is produced by Stories 1–3. It is lower priority because it only matters once the code and audit are complete, but it determines whether the compliance artifacts are actually usable in a formal context.

**Independent Test**: Can be tested by providing the compliance package to a reviewer who is familiar with TS 17961 but unfamiliar with dwsolver — if that reviewer can navigate from the compliance matrix to the evidence without asking clarifying questions, the packaging is complete.

**Acceptance Scenarios**:

1. **Given** the compliance package (compliance matrix + tooling guide + acceptance report), **when** an evaluator reads the compliance matrix, **then** each row is self-contained: rule ID, rule name, verdict, evidence reference, and remediation reference (if any).
2. **Given** the acceptance report (`specs/014-iso-ts-17961-compliance/acceptance-report.md`), **when** it is reviewed, **then** it follows the same format as `specs/008-sei-cert-c-compliance/acceptance-report.md` for continuity and presents a consolidated PASS verdict for all 22 rules.
3. **Given** the entire spec directory, **when** a reader looks for any open FRs or unresolved FAIL verdicts, **then** none exist in the final acceptance report.

---

### Edge Cases

- What if a TS 17961 rule does not have a corresponding `cppcheck --addon=cert` or `clang-tidy` checker? The compliance matrix must note "no automated checker available" and provide a grep-based absence proof: the exact `grep` command run against `src/dw_*.c` and its zero-result output confirming the triggering construct is absent from the codebase. Code review notes or narrative rationale alone are not sufficient.
- What if a rule is mechanically violated (triggering construct is present) but cannot actually trigger the defect due to program structure (e.g., a pointer is always non-null at a particular call site due to an invariant)? The compliance matrix must record the structural invariant as the evidence for PASS and verify it with a grep or manual code review — annotating source with a comment is the correct signal.
- What if spec 008 or spec 013 acceptance reports contain an error in their evidence citations? This spec's compliance matrix must cite its own independently verified evidence rather than copying citations from prior specs without re-checking.
- What if `cppcheck --addon=cert` is not available in the CI environment? The CI step must fail the build with a clear message identifying the missing tool, not silently skip the check.

---

## Requirements

### Functional Requirements

- **FR-001**: A compliance matrix MUST be produced at `specs/014-iso-ts-17961-compliance/audit/ts17961-compliance-matrix.md` covering all 22 rules defined in ISO/IEC TS 17961:2013, with a PASS / FAIL / N-A verdict and evidence citation for each rule. For rules with no automated `cppcheck --addon=cert` or `clang-tidy cert-*` checker, the evidence MUST be a grep-based absence proof: the exact `grep` command and its zero-result output; narrative rationale alone is not acceptable evidence.
- **FR-002**: The compliance matrix MUST map each TS 17961 rule to its nearest SEI CERT C rule correspondence (or note "no direct mapping") and indicate whether that correspondence was addressed by spec 008 or spec 013.
- **FR-003**: Any TS 17961 rule receiving a FAIL verdict in the initial audit MUST be remediated by a targeted code change confined to `src/dw_*.c` and `src/dw_*.h`, and the verdict updated to PASS with new evidence before the spec is closed.
- **FR-004**: All 22 rules MUST achieve a PASS or documented N-A verdict in the final compliance matrix; zero FAIL verdicts may remain in the acceptance report.
- **FR-005**: A tooling guide MUST be produced at `specs/014-iso-ts-17961-compliance/audit/tooling-guide.md` specifying the exact command-line invocations for `cppcheck --addon=cert` and `clang-tidy -checks=cert-*` used to verify compliance, including version requirements, exclusion patterns for vendored GLPK, and instructions for interpreting the output.
- **FR-006**: A new dedicated GitHub Actions workflow `.github/workflows/ci-compliance.yml` MUST be created. It runs on a Linux runner only (separate from the build/test matrix), executes `cppcheck --addon=cert --suppressions-list=specs/014-iso-ts-17961-compliance/audit/cppcheck-na-suppressions.txt` and `clang-tidy -checks=cert-*` against `src/dw_*.c`, and fails the job if any non-suppressed finding maps to a TS 17961 rule. The suppressions file lists the cppcheck checker IDs for the 6 N-A rules and is committed to the repository as an auditable artifact; it is the sole mechanism by which N-A rules are excluded from CI failure.
- **FR-007**: Instructions for optional Coverity/CodeSonar/Klocwork verification MUST be documented in the tooling guide, including the checker category names relevant to TS 17961 certification for each tool, even though those tools are not part of the mandatory CI.
- **FR-008**: An acceptance report MUST be produced at `specs/014-iso-ts-17961-compliance/acceptance-report.md` following the format of `specs/008-sei-cert-c-compliance/acceptance-report.md`, with one section per FR and an explicit PASS verdict for each.
- **FR-009**: All code changes MUST leave the full test suite passing and MUST NOT introduce any new `-std=c11 -pedantic-errors -Wall -Wextra` diagnostics relative to the spec 013 baseline.
- **FR-010**: All code changes MUST be confined to `src/dw_*.c` and `src/dw_*.h`; vendored GLPK files in `third-party/glpk/` MUST NOT be modified.

### Key Entities

- **TS 17961 Compliance Matrix** (canonical term; "audit matrix" retired): A structured document (`audit/ts17961-compliance-matrix.md`) that is the primary compliance artifact. Each row covers one of the 22 TS 17961 rules with columns for: Rule ID, Rule Name, SEI CERT C correspondence, Prior Work (spec 008/013 reference or "none"), Verdict (PASS/FAIL/N-A), and Evidence citation.
- **Tooling Guide**: A reproducibility document (`audit/tooling-guide.md`) that lets any developer or auditor re-run the compliance verification checks and obtain the same results. It covers mandatory tools (cppcheck, clang-tidy) and optional commercial tools (Coverity, Klocwork, CodeSonar).
- **Acceptance Report**: A final sign-off document (`acceptance-report.md`) that consolidates all FR verdicts into a single PASS/FAIL summary, following the format established by spec 008.

---

## Success Criteria

### Measurable Outcomes

- **SC-001**: All 22 ISO/IEC TS 17961:2013 rules are covered in the compliance matrix with a PASS or N-A verdict — verified by counting non-FAIL rows in `ts17961-compliance-matrix.md`.
- **SC-002**: Zero new findings attributable to a TS 17961 rule are introduced by this feature's code changes — verified by before/after `cppcheck --addon=cert` diff.
- **SC-003**: The CI pipeline catches at least one synthetic TS 17961 violation and reports failure — verified by a shell script (`tests/test_ts17961_enforcement.sh`) that injects a known-bad construct, runs `cppcheck --addon=cert`, asserts a non-zero exit code, then reverts the construct; the script is executed as a dedicated step in `ci-compliance.yml`.
- **SC-004**: The full test suite passes after all code changes — verified by the CI run on this branch.
- **SC-005**: The compliance package (matrix + tooling guide + acceptance report) is navigable by a reader familiar with TS 17961 but unfamiliar with dwsolver, requiring no additional interpretation — verified by a maintainer review pass.
- **SC-006**: The tooling guide is reproducible: a developer following it on a clean checkout obtains tool outputs consistent with the documented findings — verified by a second developer executing the guide independently.

---

## Assumptions

- The expected distribution of verdicts from the initial audit is: majority PASS (credited to spec 008 or spec 013 prior work), a small number requiring fresh automated-tool or grep verification, and N-A for rules whose triggering constructs (`signal()`, `setjmp`/`longjmp`, etc.) are entirely absent from `src/dw_*.c`. The actual distribution is determined by the audit; this expectation is a scoping assumption, not a requirement.
- The 22 TS 17961 rules enumerated in the 2013 edition of the Technical Specification are the authoritative scope; any corrigenda or draft updates are out of scope.
- `cppcheck` ≥ 2.12 with the CERT addon (`--addon=cert`) is available or installable in the CI environment; this is a prerequisite, not an optional capability.
- `clang-tidy` ≥ 17 is available in the CI environment; the `cert-*` check family covers the subset of TS 17961 rules that have automated checkers.
- The mapping from TS 17961 rule IDs to `cppcheck --addon=cert` and `clang-tidy cert-*` checker names is documented in the tooling guide based on the tool documentation available at the time of authoring; the mapping may not be exhaustive.
- Spec 008 and spec 013 acceptance reports are accurate and their evidence citations have been verified; where this spec cross-references them, it does so trusting their verdicts. Any discrepancy found during this audit is noted in the matrix but does not block this spec's progress.
- Coverity, Klocwork, and CodeSonar are not available in the project's CI environment; their inclusion in the tooling guide is for documentation purposes only and does not gate any CI pass/fail decision.
- The TS 17961 `datarace` rule (concurrent data access) is addressed by spec 008 (POS51/52/54-C); this spec verifies that mapping is complete and documents it, but does not re-implement threading fixes.
- Any TS 17961 rule that requires a construct entirely absent from `src/dw_*.c` (e.g., `signal()` usage for `asyncsig`, or `setjmp`/`longjmp` usage) receives an N-A verdict with a grep-based absence proof.
