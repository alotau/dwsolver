# Data Model: OpenSSF Best Practices Badge (Passing Level)

**Feature**: 015-openssf-passing-badge  
**Date**: 2026-03-24

This feature introduces no database or structured data. "Data model" captures the
inventory of artifacts to be created/modified and their required content structure.

---

## 1. New Repository Files

### 1.1 `SECURITY.md`
**Location**: repository root  
**Required sections**:
- **Supported Versions** — table of versions receiving security updates (currently: `2.x`)
- **Reporting a Vulnerability** — GitHub Security Advisory URL; do NOT open a public issue
- **Response Timeline** — acknowledgment within 14 business days; resolution target 60 days
- **Scope** — what counts as a security issue for this project
- **Security Design Notes** — numerical overflow risks in LP operations; pthreads safety model; N/A for crypto
- **Disclosure Policy** — coordinated disclosure; CVE assignment via GitHub Advisory; embargo period

### 1.2 `CONTRIBUTING.md`
**Location**: repository root  
**Required sections**:
- **Code of Conduct** — brief statement (contributor expectations)
- **Bug Reports** — open a GitHub Issue, use issue template; provide: OS, GLPK version, input file, command output
- **Enhancement Requests** — GitHub Issue with `enhancement` label; describe use case before implementation
- **Pull Requests** — fork → feature branch (see constitution naming) → PR against `main`; all CI jobs must pass
- **Coding Standards** — C11; SEI CERT C (see spec-008); ISO/IEC TS 17961 (see spec-013); POSIX.1-2008 (see spec-013)
- **DCO Sign-Off** — `git commit -s` required on every commit; by signing off the contributor certifies DCO 1.1
- **CI Requirements** — new PRs must not break: macOS CI, Linux CI, Windows CI, Docker CI, compliance CI; check all pass before requesting review

### 1.3 `CHANGELOG.md`
**Location**: repository root  
**Format**: Keep a Changelog 1.1.0 (https://keepachangelog.com)  
**Required top-level sections**: `Unreleased`, then per-release sections in descending order  
**Per-release subsections**: `Added`, `Changed`, `Fixed`, `Deprecated`, `Removed`, `Security`  
**Rule**: `Security` section is mandatory in every release entry; if empty, must state "No security changes in this release."  
**Initial content**: Entry for `[2.0.0] - 2026-xx-xx` covering all work from specs 001–015; `Security` section: "Initial public release. No known vulnerabilities."

---

## 2. Modified Repository Files

### 2.1 `README.md` — Badge Row Update
**Location**: repository root, badge row at top  
**Change**: Add OpenSSF Passing badge image (Markdown), linking to the live assessment page  
**Position**: After the last existing CI badge (Docker badge), before any blank line  
**Placeholder URL**: `https://bestpractices.coreinfrastructure.org/projects/NNNNN` (replace with real project ID after registration)

### 2.2 Source File SPDX Headers (17 files)
**Location**: `src/*.{c,h}` (15 files) + `tests/test_blas.c` + `tests/test_lib_api.c`  
**Change**: Add `SPDX-License-Identifier: GPL-3.0-or-later` inside the existing `/* ... */` license block, on the line immediately before the closing ` * ******************/` delimiter  
**Identifier**: `GPL-3.0-or-later` (correct SPDX expression for the "version 3 or later" wording used in the existing headers)  
**Existing header structure** (no change to other lines):
```c
/* ****...****
 *
 *   DWSOLVER - ...
 *   Copyright 2010 ...
 *   ...
 *   SPDX-License-Identifier: GPL-3.0-or-later        ← ADD THIS LINE
 *
 **...***/ 
```

---

## 3. External Artifact: OpenSSF Badge Assessment

**Platform**: https://bestpractices.coreinfrastructure.org  
**Required fields at registration**:
- Repository URL: `https://github.com/alotau/dwsolver`
- Project license: GPL-3.0-or-later
- Programming language: C

**Self-assessment answer map** (all Passing criteria):
- ~45 criteria: answer "Met" with one-line evidence link
- ~10 criteria: answer "N/A" (crypto section, ≥ 95% statement coverage not applicable at this maturity)
- ≥ 0 criteria: answer "Unmet" — goal is zero Unmet at submission

**Output**: Project badge ID (integer), assessment URL  
**Must be stored in**: `README.md` badge widget href and prose link in Compliance section

---

## 4. Validation Rules

| Artifact | Validation | Pass Condition |
|----------|-----------|----------------|
| `SECURITY.md` | Manual review + GitHub renders it | All required sections present; private reporting channel documented |
| `CONTRIBUTING.md` | Manual review + GitHub renders it | DCO section present; issue/PR process described; coding standards referenced |
| `CHANGELOG.md` | Manual review | Keep a Changelog format; v2.0.0 entry present; Security subsection included |
| SPDX headers | `grep -rL SPDX-License-Identifier src/ tests/*.c` returns empty | Zero files missing the identifier |
| README badge | Visual inspection; link resolves | Badge renders; link returns HTTP 200 to the project's assessment page |
| OpenSSF assessment | Platform shows "Passing" | 100% mandatory criteria marked Met; badge status = Passing |
