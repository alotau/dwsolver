# Quickstart: OpenSSF Best Practices Badge — Registration & Self-Assessment

**Feature**: 015-openssf-passing-badge  
**Audience**: Project maintainer performing badge registration

---

## Prerequisites

Before beginning the self-assessment, confirm these artifacts are committed and pushed to `main`:

- [ ] `SECURITY.md` at repository root
- [ ] `CONTRIBUTING.md` at repository root
- [ ] `CHANGELOG.md` at repository root
- [ ] SPDX identifiers in all 17 source files (`grep -rL SPDX-License-Identifier src/ tests/test_blas.c tests/test_lib_api.c` returns empty)
- [ ] README badge row updated (placeholder `NNNNN` will be replaced after registration)
- [ ] All CI jobs green on `main`

---

## Step 1 — Register the Project

1. Visit: `https://bestpractices.coreinfrastructure.org/en/projects/new`
2. Click **"Log in with GitHub"** and authorize with the `alotau` GitHub account (or the account with write access to `alotau/dwsolver`)
3. Enter the repository URL: `https://github.com/alotau/dwsolver`
4. The platform auto-detects: license (GPL-3.0), language (C), name (dwsolver)
5. Note the assigned **Project ID** (an integer in the URL: `.../projects/NNNNN`)
6. Update `README.md` badge widget: replace `NNNNN` with the real project ID

---

## Step 2 — Complete the Self-Assessment

The assessment has ~75 questions across categories. Use the evidence table below.

### Basics

| Question | Answer | Evidence |
|----------|--------|----------|
| What is the project's license? | GPL-3.0-or-later | `COPYING` |
| Is the license OSI-approved? | Met | GPL-3.0 is OSI-approved |
| Does the project have a public repo? | Met | `https://github.com/alotau/dwsolver` |
| Is the VCS distributed? | Met | Git |

### Documentation

| Question | Answer | Evidence |
|----------|--------|----------|
| Does the project have basic documentation? | Met | `README.md` |
| Does the project have reference documentation for the external interface? | Met | `README.md` → Usage + Library API section; `src/dw_solver.h` header |
| Is there a CONTRIBUTING guide? | Met | `CONTRIBUTING.md` |

### Change Control

| Question | Answer | Evidence |
|----------|--------|----------|
| Do releases have version identifiers? | Met | `configure.ac` `AC_INIT([dwsolver],[2.0.0])` |
| Are version identifiers unique and semver? | Met | SemVer; `v*` tags in `release.yml` |
| Are release notes provided? | Met | `CHANGELOG.md` |
| Do release notes address known vulnerabilities? | Met | `CHANGELOG.md` Security subsection |

### Reporting

| Question | Answer | Evidence |
|----------|--------|----------|
| Is there a bug-reporting process? | Met | GitHub Issues; documented in `CONTRIBUTING.md` |
| Is the bug tracker public? | Met | GitHub Issues — public |
| Is there an enhancement-request tracker? | Met | GitHub Issues with `enhancement` label |

### Quality

| Question | Answer | Evidence |
|----------|--------|----------|
| Is there a build system? | Met | GNU Autotools (`./configure && make`) |
| Is the build system common/standard? | Met | GNU Autotools is standard for C projects |
| Does the project have a test suite? | Met | `tests/dw-tests.sh`; `make check` |
| Is the test suite invocable? | Met | `make check` or `bash tests/dw-tests.sh` |
| Is there automated CI? | Met | GitHub Actions: 5 platform/configuration workflows |
| Does the project have static analysis? | Met | `ci-compliance.yml`: clang-tidy cert-* + cppcheck |
| Does static analysis catch common vulnerabilities? | Met | cert-* rules cover CWE classes |
| Is dynamic analysis (e.g., sanitizers) used? | Met | `ci-linux.yml` `linux-asan-ubsan` job |
| Does dynamic analysis use memory-unsafe options disabled (i.e., are unsafe constructs caught)? | Met | `-fsanitize=address,undefined` fails on detection |

### Security

| Question | Answer | Evidence |
|----------|--------|----------|
| Is there a vulnerability reporting process? | Met | `SECURITY.md` → GitHub Security Advisory |
| Does the project respond within 14 days? | Met | `SECURITY.md` states 14-day acknowledgment SLA |
| Does the project fix known vulnerabilities within 60 days? | Met | `SECURITY.md` states 60-day resolution target |
| Is cryptography used? | N/A | dwsolver performs no cryptographic operations |
| (All crypto sub-criteria) | N/A | No crypto |

### License

| Question | Answer | Evidence |
|----------|--------|----------|
| Is the license a FLOSS license? | Met | GPL-3.0-or-later |
| Is the license in the repo? | Met | `COPYING` |
| Is the license on each source file? | Met | SPDX `GPL-3.0-or-later` in all `src/*.{c,h}` |

---

## Step 3 — Submit and Record

1. Review all answers; ensure no mandatory criterion is marked "Unmet"
2. Click **"Submit"**
3. The platform shows the badge status: **Passing** ✅
4. Copy the badge Markdown from the platform and confirm the `README.md` badge widget matches
5. Record the assessment URL (`.../projects/NNNNN`) in the README Compliance section prose

---

## Step 4 — Post-Registration Maintenance

- **Each release**: Update `CHANGELOG.md` with a `Security` subsection
- **Each vulnerability report**: Follow the process in `SECURITY.md`; update the OpenSSF assessment if response/resolution timelines change
- **Re-assessment**: OpenSSF badges do not expire, but the platform sends annual reminders to re-confirm; update evidence links if CI or file paths change
