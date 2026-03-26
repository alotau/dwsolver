# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

## [2.0.0] - 2025-07-01

### Added

- Callable library API (`libdwsolver`) with public header `dw.h`, enabling
  callers to embed Dantzig-Wolfe decomposition into their own applications.
- `pkg-config` support (`dwsolver.pc`) for downstream build-system integration.
- POSIX thread-based parallel subproblem solving.
- GitHub Actions CI across Linux (GCC and Clang, ASAN+UBSan, TSan), macOS,
  Windows (MinGW-w64), and Docker.
- SEI CERT C and ISO/IEC TS 17961 static-analysis compliance CI job.
- SPDX-License-Identifier comment (`GPL-3.0-or-later`) in every source file.
- `SECURITY.md` with vulnerability reporting process and security design notes.
- `CONTRIBUTING.md` with DCO sign-off requirement, coding standards, and CI
  expectations for contributors.

### Security

Initial public release. No known vulnerabilities. The solver accepts LP input
files only; it performs no network access, no shell command construction, and no
cryptographic operations.

---

[Unreleased]: https://github.com/alotau/dwsolver/compare/v2.0.0...HEAD
[2.0.0]: https://github.com/alotau/dwsolver/releases/tag/v2.0.0
