# Contract: SECURITY.md

**Feature**: 015-openssf-passing-badge  
**Satisfies**: FR-001, FR-012 (OpenSSF `vulnerability_report_process`, `vulnerability_report_response`)

This document defines the required content contract for the `SECURITY.md` file
that will be created at the repository root.

---

## Required Sections and Content

### Section 1 — Supported Versions
A Markdown table listing which versions receive security updates.

| Version | Supported |
|---------|-----------|
| 2.x (current) | ✅ Yes |
| 1.x and earlier | ❌ No |

### Section 2 — Reporting a Vulnerability
Instructions must:
- Explicitly state: do NOT open a public GitHub issue for security vulnerabilities
- Direct reporters to GitHub's private Security Advisory feature:
  - URL: `https://github.com/alotau/dwsolver/security/advisories/new`
- Include a contact email as a fallback for reporters who cannot use GitHub

### Section 3 — Response Timeline
Must state:
- Acknowledgment: within 14 business days of report receipt
- Resolution target: 60 days from acknowledgment (or coordinated disclosure timeline if a CVE is involved)
- If a fix requires longer, the reporter will be notified with a status update

### Section 4 — Scope
Define what counts as a security vulnerability for this project:
- Memory safety issues in solver or library code (buffer overflows, use-after-free, etc.)
- Integer overflow or undefined behavior that could be exploited via crafted LP input files
- Data races in multi-threaded execution paths
- Dependency vulnerabilities in GLPK that affect dwsolver's security posture

Out of scope (document explicitly):
- Denial of service via pathological LP inputs (inherent to NP-hard computation — not a vulnerability)
- Vulnerabilities in third-party tools used only at build time (automake, libtool, etc.)
- Theoretical issues with no known exploit path

### Section 5 — Security Design Notes (FR-012)
Brief (3–5 bullet) security overview:
- **Numerical operations only**: dwsolver performs no cryptographic operations; all crypto-related OpenSSF criteria are N/A
- **Input validation**: LP input is read from `.cpxlp` / GLPK-format files; malformed input triggers GLPK error paths; no shell escaping or injection surface
- **Integer overflow risk**: loop indices and matrix dimension arithmetic use `int`; very large LP instances could theoretically overflow; this is documented as a known limitation
- **Thread safety**: shared state is protected by `pthread_mutex_t`; the public library API is documented as not thread-safe across multiple concurrent calls to `dw_solve()` with the same context object
- **No network access**: dwsolver makes no network calls; attack surface is limited to local file input and GLPK shared library

### Section 6 — Disclosure Policy
- Private advisory opened; triaged and acknowledged within 14 business days
- Fix developed on private branch; coordinated with reporter
- CVE requested via GitHub Advisory if severity warrants
- Public disclosure after fix is released and tagged; reporter credited (with permission)

---

## Validation
The implementation is complete when:
- `SECURITY.md` exists at the repository root
- All six sections above are present
- The GitHub Security Advisory link is live (test: visit the URL while authenticated)
- Running `cat SECURITY.md | grep -c "14"` returns ≥ 1 (14-day acknowledgment documented)
