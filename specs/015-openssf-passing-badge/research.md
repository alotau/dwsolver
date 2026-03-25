# Research: OpenSSF Best Practices Badge (Passing Level)

**Feature**: 015-openssf-passing-badge  
**Date**: 2026-03-24

---

## 1. Criteria Gap Analysis — What dwsolver Already Satisfies

### Decision
Many Passing criteria are already met. Implementation scope is narrow.

### Evidence of met criteria
| Criterion | Status | Evidence |
|-----------|--------|----------|
| `floss_license` | ✅ Met | GPL-3.0; OSI-approved; `COPYING` at repo root |
| `floss_license_osi` | ✅ Met | GPL-3.0 is OSI-approved |
| `repo_public` | ✅ Met | Public GitHub repository |
| `repo_vcs` | ✅ Met | Git (distributed VCS) |
| `repo_distributed` | ✅ Met | GitHub + every clone is a full history copy |
| `build` | ✅ Met | GNU Autotools (`./configure && make`); documented in README |
| `build_common_tools` | ✅ Met | GCC/Clang + standard POSIX tools; MSYS2/MinGW on Windows |
| `test` | ✅ Met | `tests/dw-tests.sh`; `make check` via `make distcheck` in CI |
| `test_invocation` | ✅ Met | `make check` or `bash dw-tests.sh` documented in README |
| `automated_integration_testing` | ✅ Met | `dw-tests.sh` exercises full solver I/O end-to-end in CI |
| `static_analysis` | ✅ Met | `ci-compliance.yml`: clang-tidy (cert-* rules) + cppcheck |
| `static_analysis_common_vulnerabilities` | ✅ Met | cert-* rules in clang-tidy cover CWE/CERT vulnerability classes |
| `dynamic_analysis` | ✅ Met | `linux-asan-ubsan` job in `ci-linux.yml` |
| `dynamic_analysis_unsafe` | ✅ Met | ASan+UBSan flags; job fails on any detected error |
| `dynamic_analysis_enable_assertions` | ✅ Met | NDEBUG not set in sanitizer build; assertions active |
| `version_controlled` | ✅ Met | Git tags; `v*` tag pattern triggers release workflow |
| `version_semver` | ✅ Met | `AC_INIT([dwsolver], [2.0.0], ...)` in `configure.ac`; SemVer pattern |
| `version_tags` | ✅ Met | `release.yml` triggered by `v*` tags; GitHub Releases documented |
| `crypto_call_network` | ✅ N/A | No cryptographic operations performed |
| `crypto_used` | ✅ N/A | Solver is purely numerical; no crypto |
| `delivery_mitm` | ✅ Met | HTTPS delivery via GitHub; release tarballs verifiable via `make distcheck` |

### What Remains (Gaps)
| Criterion | Gap | Required Action |
|-----------|-----|----------------|
| `vulnerability_report_process` | No `SECURITY.md` | Create `SECURITY.md` with private reporting channel |
| `vulnerability_report_response` | Undocumented 14-day SLA | Add to `SECURITY.md` |
| `vulnerabilities_fixed_60_days` | No policy stated | Document in `SECURITY.md` — commit to 60-day resolution target |
| `contribution_requirements` | No `CONTRIBUTING.md` | Create with DCO sign-off requirement and issue/PR process |
| `report_tracker` | Not cross-referenced | Point to GitHub Issues in `CONTRIBUTING.md` |
| `enhancement_tracker` | Not documented | Same — GitHub Issues (feature requests) |
| `license_per_file` | No SPDX identifiers in source | Add `SPDX-License-Identifier: GPL-3.0-or-later` to 17 files |
| `release_notes` | ChangeLog exists but lacks consistent per-release security subsection | Create `CHANGELOG.md` with structured format |
| `release_notes_vulns` | No "Security" section in existing changelog | New `CHANGELOG.md` template enforces this |

---

## 2. SPDX License Identifier Placement

### Decision
Add `SPDX-License-Identifier: GPL-3.0-or-later` as the **final line inside the existing `/* ... */` license header block** in each file, as a `/* */` comment on its own line before the closing `****/`.

### Rationale
- SPDX best practice places the identifier in a comment near the top of the file, inside the copyright/license block.
- Preserving the existing verbose GPL notice maintains human readability; the SPDX tag adds machine-parseability for tools like `licensee`, `scancode`, and `FOSSology`.
- `GPL-3.0-or-later` is the correct SPDX expression for "GPL v3 or (at your option) any later version" (used throughout the existing headers).

### Files requiring SPDX identifiers
`src/`: `dw.h`, `dw_blas.c`, `dw_blas.h`, `dw_globals.c`, `dw_main.c`, `dw_phases.c`, `dw_phases.h`, `dw_rounding.c`, `dw_rounding.h`, `dw_solver.c`, `dw_solver.h`, `dw_subprob.c`, `dw_subprob.h`, `dw_support.c`, `dw_support.h`  
`tests/`: `test_blas.c`, `test_lib_api.c`

### Alternatives considered
- Separate `.license` sidecar files per REUSE spec: rejected — more complex, requires REUSE tooling; SPDX in-file is simpler and widely accepted.
- A `LICENSES/` directory: complementary to in-file SPDX, not a replacement; out of scope for Passing level.

---

## 3. SECURITY.md Content Requirements

### Decision
`SECURITY.md` must cover: (1) supported versions, (2) reporting channel — GitHub private Security Advisory, (3) 14-day acknowledgment SLA, (4) 60-day resolution target, (5) brief security design overview scoped to this project.

### Rationale
- OpenSSF `vulnerability_report_process` requires a documented *private* reporting mechanism. GitHub's Security Advisory feature is purpose-built for this and requires zero additional infrastructure.
- `vulnerability_report_response` requires documented response time. 14-day acknowledgment is the OpenSSF-stated threshold.
- The security design note (FR-012) covers: numerical overflow risk in LP operations, pthreads data-race invariants, and explicit N/A for crypto.

### Alternatives considered
- Email-only reporting: acceptable but less structured; GitHub Security Advisories preferred because they support private CVE assignment and are integrated with the repository.
- HackerOne: overkill for a research tool; not warranted at current adoption level.

---

## 4. CONTRIBUTING.md Content Requirements

### Decision
`CONTRIBUTING.md` must cover: (1) bug reporting via GitHub Issues, (2) enhancement requests via GitHub Issues with the `enhancement` label, (3) PR process (fork → branch → PR against `main`), (4) branch naming convention from the constitution, (5) coding standards reference (C11; CERT C; ISO/IEC TS 17961), (6) DCO sign-off requirement (`git commit -s`), (7) CI requirements before merge.

### Rationale
- `contribution_requirements` explicitly requires a document telling users how contributions are accepted — it is one of the most commonly missing Passing criteria found in assessments of similar OSS projects.
- DCO is the lightweight alternative to a CLA appropriate for a NASA-originated research tool. It documents that contributors affirm they have the rights to submit the contribution without requiring a separate legal agreement.
- Constitution Git workflow section already defines branch naming; CONTRIBUTING.md simply references it rather than duplicating it.

### Alternatives considered
- CLA (Contributor License Agreement): heavier process; not warranted for a research tool; DCO is sufficient.
- Inline in README: less discoverable; GitHub surfaces CONTRIBUTING.md automatically in the new-issue and new-PR UI.

---

## 5. CHANGELOG.md Format

### Decision
Adopt [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format with sections: `Added`, `Changed`, `Fixed`, `Security`. Create `CHANGELOG.md` as the going-forward changelog. The existing `ChangeLog` (GNU format, sparse) is kept for historical provenance but marked as superseded in its header.

### Rationale
- OpenSSF `release_notes_vulns` requires that each release's notes address known vulnerabilities. Keep a Changelog mandates a `Security` subsection, enforcing this naturally.
- The format is widely understood, machine-parseable by changelog tools, and compatible with SemVer.
- Existing `ChangeLog` covers versions up to 1.2 (2010); `CHANGELOG.md` will start at `2.0.0` (current release).

### Alternatives considered
- Continue using GNU `ChangeLog` format: rejected because it has no standard `Security` section and is harder for downstream tooling to parse.
- GitHub Releases only: releases already exist in `release.yml` but they don't satisfy the `release_notes` criterion without a committed file; a committed `CHANGELOG.md` is more durable.

---

## 6. Badge Registration Process

### Decision
Register via `https://bestpractices.coreinfrastructure.org/en/projects/new`, authenticate with the project's GitHub account, link `https://github.com/alotau/dwsolver`. Record the badge project ID and assessment URL in `README.md` (both in the badge widget and as a prose link in the Compliance section).

### Rationale
- The OpenSSF platform auto-detects public GitHub repos and pre-fills many criteria from repository metadata (license file, CI configuration, etc.), reducing the self-assessment burden.
- Storing the assessment URL in the README (FR-011) ensures it remains discoverable even if the badge image CDN changes.

### Evidence links to prepare for the self-assessment
| Criterion | Evidence URL pattern |
|-----------|---------------------|
| `floss_license` | `COPYING` file in repo |
| `build` | `README.md` → Build section |
| `test` | `tests/dw-tests.sh` |
| `static_analysis` | `.github/workflows/ci-compliance.yml` |
| `dynamic_analysis` | `.github/workflows/ci-linux.yml` → `linux-asan-ubsan` job |
| `contribution_requirements` | `CONTRIBUTING.md` (to be created) |
| `vulnerability_report_process` | `SECURITY.md` (to be created) |
| `license_per_file` | `src/*.{c,h}` after SPDX tagging |
| `release_notes_vulns` | `CHANGELOG.md` (to be created) |
