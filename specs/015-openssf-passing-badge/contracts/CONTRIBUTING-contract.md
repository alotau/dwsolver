# Contract: CONTRIBUTING.md

**Feature**: 015-openssf-passing-badge  
**Satisfies**: FR-002 (OpenSSF `contribution_requirements`, `report_tracker`, `enhancement_tracker`)

This document defines the required content contract for the `CONTRIBUTING.md` file
that will be created at the repository root.

---

## Required Sections and Content

### Section 1 — Welcome
One paragraph: who can contribute, spirit of the project (research tool, NASA origins), openness to external contributions.

### Section 2 — Bug Reports
Must state:
- Use GitHub Issues: `https://github.com/alotau/dwsolver/issues`
- Include in the report:
  - Operating system and version
  - GLPK version (`glpsol --version` or `pkg-config --modversion glpk`)
  - Build configuration (output of `./configure ...`)
  - Minimal reproducing input file (`.cpxlp` or equivalent)
  - Full command invocation and output
- Check for duplicate issues before filing

### Section 3 — Enhancement Requests
Must state:
- Use GitHub Issues with the `enhancement` label
- Describe the use case and motivation *before* proposing an implementation
- If the enhancement touches solver math, reference the relevant paper or textbook section

### Section 4 — Pull Requests
Step-by-step process:
1. Fork the repository and clone locally
2. Create a branch following the project naming convention: `fix/<description>`, `feat/<description>`, or `chore/<description>`
3. Make changes; ensure all CI jobs pass locally where possible
4. Sign off every commit with `git commit -s` (DCO sign-off — see below)
5. Push the branch and open a PR against `main`
6. Fill in the PR description: what is broken/missing, what was changed, how it was tested, platforms verified
7. Address review comments; do not force-push after review begins unless asked

### Section 5 — Coding Standards
Reference (do not duplicate) the compliance specifications:
- Language: C11 (GCC or Clang; MinGW-w64 on Windows)
- Memory: no heap leaks; AddressSanitizer clean (`make check` with ASAN)
- Thread safety: no data races; TSan clean
- Compliance: SEI CERT C (spec-008); ISO/IEC TS 17961 (spec-013); POSIX.1-2008 (spec-013)
- Static analysis must pass: clang-tidy (cert-* rules) + cppcheck

### Section 6 — Developer Certificate of Origin (DCO)
Must state:
- All contributions require a DCO sign-off
- Add `Signed-off-by: Your Name <email@example.com>` using `git commit -s`
- By signing off, the contributor certifies they have the right to submit the contribution under the project's GPL-3.0-or-later license, per the [Developer Certificate of Origin v1.1](https://developercertificate.org/)
- Include a verbatim or linked copy of the DCO 1.1 text

### Section 7 — CI Requirements
State explicitly:
- All five CI workflows must pass before a PR can be merged: macOS, Linux (standard + ASAN/UBSan + TSan), Windows, Docker, Compliance
- PRs that break any CI job will not be merged until fixed

---

## Validation
The implementation is complete when:
- `CONTRIBUTING.md` exists at the repository root
- All seven sections above are present
- `grep -c "Signed-off-by\|DCO" CONTRIBUTING.md` returns ≥ 2
- `grep -c "GitHub Issues" CONTRIBUTING.md` returns ≥ 2 (bug reports + enhancements)
- GitHub renders the file and shows it in the new-issue UI dropdown
