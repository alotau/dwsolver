# Contributing to dwsolver

Thank you for your interest in contributing to **dwsolver** — a Dantzig-Wolfe decomposition LP solver written in C11. Contributions of all kinds are welcome: bug reports, code fixes, documentation improvements, and example problems.

Please read this guide before opening an issue or pull request. It keeps the process smooth for everyone.

---

## Table of Contents

1. [Bug Reports](#bug-reports)
2. [Enhancement Requests](#enhancement-requests)
3. [Pull Requests](#pull-requests)
4. [Coding Standards](#coding-standards)
5. [DCO Sign-Off](#dco-sign-off)
6. [CI Requirements](#ci-requirements)
7. [Security Issues](#security-issues)

---

## Bug Reports

Open a new ticket on [GitHub Issues](https://github.com/alotau/dwsolver/issues) with the label **`bug`**.

To help us reproduce and fix the problem quickly, please include **all** of the following:

| Item | Details |
|------|---------|
| **Operating system and version** | e.g., Ubuntu 24.04 LTS, macOS 14.5, Windows 11 |
| **Compiler and version** | e.g., `gcc --version` → `gcc (Ubuntu 13.3.0) 13.3.0` |
| **GLPK version** | e.g., `glpsol --version` → `GLPK 5.0` |
| **Build configuration** | The exact `./configure` invocation you used |
| **Minimal input file** | A `.cpxlp` or MPS file that reproduces the problem (attach it to the issue) |
| **Full terminal output** | Everything from `./configure` (or `make`) through the error — no truncation |
| **Expected vs. actual behaviour** | What you expected to happen and what actually happened |

Before filing, please search existing issues to avoid duplicates.

---

## Enhancement Requests

Open a new ticket on [GitHub Issues](https://github.com/alotau/dwsolver/issues) and apply the **`enhancement`** label.

Describe:
- The problem or limitation you are trying to address
- Your proposed solution or the behaviour you would like to see
- Any relevant references (papers, algorithms, other solver implementations)

Enhancement requests will be discussed in the issue thread before any code is written.

---

## Pull Requests

### Workflow

1. **Fork** the repository on GitHub.
2. Create a branch off `main` with a descriptive name, using one of these prefixes:
   - `fix/` — bug fix (e.g., `fix/negative-y-bounds-crash`)
   - `feat/` — new feature (e.g., `feat/mps-format-output`)
   - `chore/` — maintenance, docs, CI (e.g., `chore/update-changelog`)
3. Make your changes, following the [Coding Standards](#coding-standards) below.
4. Ensure **all CI checks pass** (see [CI Requirements](#ci-requirements)).
5. Add or update tests for the changed behaviour in `tests/`.
6. Commit with a [DCO sign-off](#dco-sign-off) (`git commit -s`).
7. Open a pull request against the `main` branch. Fill in the PR template.

### What makes a good PR

- **Focused**: one logical change per PR. Large refactors should be discussed in an issue first.
- **Well-described**: explain *why* the change is needed, not just *what* it does.
- **Tested**: add or modify tests that demonstrate the fix or feature.
- **Clean history**: squash WIP commits before requesting review.

---

## Coding Standards

All C code in this project must conform to the following standards. CI enforces these automatically.

| Standard | Rationale |
|----------|-----------|
| **C11** (`-std=c11`) | Minimum required language version |
| **POSIX.1-2008** | Portable system interfaces |
| **SEI CERT C Coding Standard** | Memory safety, integer safety, undefined-behaviour avoidance |
| **ISO/IEC TS 17961:2013** | C secure coding rules — static-analysis verifiable |

### Practical guidelines

- No compiler warnings. Build must succeed with `-Wall -Wextra -Wpedantic`.
- No undefined behaviour. All CI jobs run with `-fsanitize=address,undefined`.
- Functions that allocate memory must document and handle allocation failure.
- Use fixed-width integer types (`int32_t`, etc.) when overflow behaviour matters.
- Do not use `gets()`, `sprintf()`, or other deprecated interfaces.
- String buffers must be null-terminated; use `snprintf()` with explicit size limits.
- All new public API symbols must be documented in the relevant header file.

---

## DCO Sign-Off

This project uses the [Developer Certificate of Origin (DCO) 1.1](https://developercertificate.org/) instead of a Contributor Licence Agreement. By signing off your commits you certify that you have the right to submit the code under the project's licence (GPL-3.0-or-later).

### How to sign off

Add the `-s` flag to every commit:

```sh
git commit -s -m "fix: correct off-by-one in phase boundary check"
```

This appends the following trailer to the commit message:

```
Signed-off-by: Your Name <your.email@example.com>
```

The name and email must match your Git identity (`git config user.name` / `git config user.email`).

### What signing certifies

By adding your `Signed-off-by` line you certify, to the best of your knowledge, that:

> (a) The contribution was created in whole or in part by me and I have the right to submit it under the open source licence indicated in the file; or
>
> (b) The contribution is based upon previous work that, to the best of my knowledge, is covered under an appropriate open source licence and I have the right to submit that work with modifications under the same licence; or
>
> (c) The contribution was provided directly to me by some other person who certified (a), (b), or (c) and I have not modified it.

See [https://developercertificate.org/](https://developercertificate.org/) for the full text.

Pull requests with unsigned commits will not be merged.

---

## CI Requirements

All five CI workflows must report **green** before a pull request can be merged:

| Workflow | What it checks |
|----------|---------------|
| `ci-linux.yml` | Standard build + tests on Linux (GCC and Clang); ASAN+UBSan; TSan |
| `ci-macos.yml` | Build + tests on macOS (Apple Clang) |
| `ci-windows.yml` | Cross-compilation via MinGW-w64 |
| `ci-docker.yml` | Reproducible build in the project Docker image |
| `ci-compliance.yml` | Static analysis and SEI CERT C / ISO TS 17961 conformance |

If CI fails on your PR, fix the issues before requesting a review. If you cannot determine why a CI job is failing, describe the failure in the PR comments and ask for help — do not ignore it.

---

## Security Issues

**Please do not report security vulnerabilities in a public GitHub issue.**  
See [SECURITY.md](SECURITY.md) for the private disclosure process.
