# Feature Specification: Strict ISO C and POSIX.1 Compliance

**Feature Branch**: `013-strict-c-posix-compliance`
**Created**: 2026-03-23
**Status**: Draft
**Input**: Add a spec to make this project Strict ISO C and POSIX.1 compliant. Ensure that SEI CERT C compliance is not compromised in the process.

## Background

The dwsolver source files (`src/dw_*.c`, `src/dw_*.h`) currently build without any enforced C language standard. No `-std=` flag is passed by the build system, `configure.ac` defines no `_POSIX_C_SOURCE` or `_GNU_SOURCE` feature-test macro, and `Makefile.am` sets no `-pedantic` flag. The code reaches into POSIX-only headers — `<unistd.h>`, `<pthread.h>`, `<semaphore.h>`, `<fcntl.h>`, `<sys/time.h>` — without declaring which standard it requires from the system. One file (`dw_subprob.c`, `dw_rounding.c`) conditionally includes Intel MKL headers under `#ifdef USE_INTEL_MKL`; those guards already exist and are out of scope.

This spec defines the work needed to:

1. Explicitly target a single ISO C standard revision across all authored files.
2. Add the appropriate POSIX feature-test macro so the system exposes exactly the POSIX.1-2008 API surface and nothing beyond it.
3. Eliminate any compiler-specific language extensions from the authored source.
4. Verify that the SEI CERT C rules addressed in feature 008 are not weakened by these changes.

**Scope**: `src/dw_*.c` and `src/dw_*.h`, `configure.ac`, and `src/Makefile.am`.  
Vendored GLPK (`third-party/glpk/`) and Intel MKL conditional blocks are explicitly excluded.

---

## Clarifications

### Session 2026-03-23

- Q: Which ISO C standard revision should be targeted (C99, C11, or C17)? → A: C11 (`-std=c11`)
- Q: Should `gettimeofday` (POSIX.1-2008 obsolescent) be migrated to `clock_gettime`, or retained under the existing configure guard? → A: Migrate to `clock_gettime(CLOCK_MONOTONIC, ...)` guarded by a new `HAVE_CLOCK_GETTIME` autoconf check; retain `gettimeofday` as a `#else` fallback
- Q: Does Windows CI verification need a new FR/SC, or is the existing `ci-windows.yml` job sufficient? → A: Existing `ci-windows.yml` (MSYS2/MinGW-w64, `continue-on-error: true`) is sufficient; no new requirement needed
- Q: Should pedantic diagnostics be hard build errors (`-pedantic-errors`) or warnings only (`-pedantic`)? → A: `-pedantic-errors` — violations are hard errors, enforced in CI
- Q: What format should the compliance audit document follow? → A: Align with the feature 008 acceptance-report format (`specs/008-sei-cert-c-compliance/acceptance-report.md`): one section per FR with evidence bullets and an explicit ✅ PASS / ❌ FAIL verdict

---

## User Scenarios & Testing

### User Story 1 — Source compiles clean under an explicit ISO C standard (Priority: P1)

A developer building dwsolver can pass `-std=c11 -pedantic` to any conforming C compiler and receive zero warnings or errors for the authored source files `src/dw_*.c` and `src/dw_*.h`.

**Why this priority**: Without a declared standard, the compiler silently accepts GNU extensions, VLAs used as goto targets, and other non-portable constructs. Locking in a standard is the prerequisite for every other item in this feature and signals portability intent to downstream users.

**Independent Test**: Can be fully tested in isolation by configuring the build system with the new flags and running `make` on a clean tree; any diagnostic produced for an authored file is a failure. Delivers a portable, reproducible baseline for contributors.

**Acceptance Scenarios**:

1. **Given** a clean build tree, **when** `./configure && make` is run on a system with GCC or Clang, **then** compilation of every `src/dw_*.c` file produces zero warnings or errors attributable to language-standard or pedantic violations.
2. **Given** the updated `configure.ac`, **when** `AC_PROG_CC` runs with strict standard enforcement enabled, **then** the configure script correctly detects whether the compiler supports the chosen standard and aborts with a clear error if it does not.
3. **Given** the compiled artifacts, **when** the full test suite (`tests/dw-tests.sh`) is run, **then** all tests pass — confirming that enforcing the standard did not silently change any runtime behaviour.

---

### User Story 2 — POSIX.1-2008 feature-test macro is declared project-wide (Priority: P1)

A developer reading any `src/dw_*.c` file can see immediately which POSIX API surface the file depends on. The build system sets `_POSIX_C_SOURCE=200809L` globally so that every translation unit uses a consistent, documented API contract — no file relies on the implicit GNU extension namespace.

**Why this priority**: Without `_POSIX_C_SOURCE`, glibc on Linux silently augments the POSIX namespace with GNU extensions. Code written under that implicit contract will fail to compile on stricter POSIX-only platforms (musl libc, FreeBSD default mode, MSYS2 with strict headers) and may silently use deprecated or removed POSIX interfaces.

**Independent Test**: Can be fully tested in isolation by building on a system with musl libc (e.g., Alpine Linux) or by passing `-D_POSIX_C_SOURCE=200809L -D_XOPEN_SOURCE=700` manually and confirming zero diagnostics; the test suite passes.

**Acceptance Scenarios**:

1. **Given** the updated `configure.ac` or `src/Makefile.am`, **when** the build system prepends `-D_POSIX_C_SOURCE=200809L` to `AM_CPPFLAGS` (or equivalent), **then** every translation unit in scope receives that macro before any system header is included.
2. **Given** the macro in effect, **when** GCC or Clang compiles the authored source, **then** zero implicit-declaration warnings are emitted for any POSIX function used by the code (e.g., `pthread_create`, `sem_init`, `clock_gettime`, `nanosleep`).
3. **Given** the macro in effect, **when** all function calls in `src/dw_*.c` are reviewed, **then** every function used is defined in POSIX.1-2008 or ISO C99 — no call relies on a GNU-only or POSIX.1-2024-only extension.

---

### User Story 3 — Compiler-specific extensions are removed from authored source (Priority: P2)

A developer porting dwsolver to a non-GCC/non-Clang toolchain encounters no unguarded GCC-specific language constructs in `src/dw_*.c` or `src/dw_*.h`. Any remaining compiler-specific syntax is isolated behind a feature-detection guard consistent with existing patterns in `dw.h`.

**Why this priority**: The `dw.h` header already correctly guards `__attribute__((visibility))` and `__declspec` behind `#ifdef __GNUC__` blocks. This story ensures no other compiler extensions have crept in and that `-pedantic` enforces that contract continuously.

**Independent Test**: Can be tested independently by compiling the authored files with a pedantic-mode Clang scan (`-Wpedantic -Wno-language-extension-token`) and confirming zero extension warnings for unguarded constructs; no test-suite regression.

**Acceptance Scenarios**:

1. **Given** the source files, **when** an audit is run for unguarded GCC-specific constructs (zero-length arrays, `__typeof__`, statement expressions, `__builtin_*` calls, `__attribute__` outside the existing visibility guard), **then** none are found in `src/dw_*.c` or `src/dw_*.h`.
2. **Given** `-pedantic-errors` is added to the compiler flags, **when** the project builds, **then** compilation succeeds with zero new errors relative to the baseline (i.e., no previously-latent extension errors surface beyond those in the pre-existing vendor code).
3. **Given** the modified source, **when** the full test suite is run, **then** all tests pass with no regression.

---

### User Story 4 — SEI CERT C compliance is verified to be non-regressed (Priority: P2)

After all ISO C and POSIX.1 compliance changes are applied, the same static analysis tools used in feature 008 (`cppcheck`, `clang --analyze`) are re-run and their output is compared against the feature 008 baselines. The compliance count does not increase.

**Why this priority**: Mechanical changes — especially type changes to fix signed/unsigned warnings, removal of implicit declarations, or adjusting header include order — can inadvertently introduce new code paths or hide existing error-check logic. An explicit regression gate prevents that.

**Independent Test**: Can be tested independently by running `cppcheck --enable=all src/dw_*.c` and `clang --analyze -I src/ src/dw_*.c` and diffing against `specs/008-sei-cert-c-compliance/baseline-cppcheck.txt` and `specs/008-sei-cert-c-compliance/baseline-warnings.txt`; zero new issues is the pass criterion.

**Acceptance Scenarios**:

1. **Given** the modified source, **when** `cppcheck --enable=all src/dw_*.c` is run, **then** no new warnings appear that were not present in the feature 008 baseline.
2. **Given** the modified source, **when** `clang --analyze -I src/ src/dw_*.c` is run, **then** no new analyzer findings appear beyond those recorded in the feature 008 baseline.
3. **Given** the modified source is compiled with `-Wall -Wextra -Wsign-compare`, **then** the warning count for `src/dw_*.c` does not exceed the count recorded in `specs/008-sei-cert-c-compliance/baseline-warnings.txt`.
4. **Given** the modified code, **when** the full test suite is run, **then** all tests pass with no regression.

---

### Edge Cases

- `gettimeofday` is **not present** in any authored source file (zero call sites confirmed by Phase 0 research). The `HAVE_GETTIMEOFDAY` and `HAVE_SYS_TIME_H` configure checks in `configure.ac` are vestigial and MUST be removed (see FR-009, T002). The `#else` fallback for platforms lacking `clock_gettime` is the ISO C `time()`/`clock()` path — equally portable and requiring no POSIX. This migration adds a monotonic clock path aligned with CERT MSC15-C while removing a misleading vestigial configure check.
- `<sys/time.h>` is not ISO C but is POSIX. It is already conditionally included via `HAVE_SYS_TIME_H`; that guard must be preserved.
- `dw.h` includes `<semaphore.h>`, `<fcntl.h>`, and `<pthread.h>` unconditionally at file scope before the header guard. When `_POSIX_C_SOURCE` is active, these headers expose their POSIX.1-2008 interfaces as expected. No change is required, but the include-before-guard pattern should be noted in the audit.
- Windows compatibility (MinGW/MSYS2) is covered by the existing `ci-windows.yml` GitHub Actions job (MSYS2 MINGW64, `continue-on-error: true`). `_POSIX_C_SOURCE=200809L` is recognised by MinGW-w64 headers. The existing `#if defined(_WIN32)` visibility guards in `dw.h` must not be broken; their correctness will be confirmed by that CI job.
- The Intel MKL optional path (`#ifdef USE_INTEL_MKL`) introduces non-standard headers. Because those blocks are already fully guarded, they are explicitly out of scope for this feature and must not be touched.

---

## Requirements

### Functional Requirements

- **FR-001**: The build system MUST pass `-std=c11` to all `src/dw_*.c` compilations via `AM_CFLAGS` or an equivalent autoconf mechanism.
- **FR-002**: The build system MUST define `_POSIX_C_SOURCE=200809L` globally via `AM_CPPFLAGS` or an `AC_DEFINE` call in `configure.ac` so that every translation unit in scope receives it.
- **FR-003**: Compilation of all `src/dw_*.c` and `src/dw_*.h` files MUST produce zero errors when `-pedantic-errors` is active alongside `-std=c11`. Pedantic violations are treated as hard build errors, not warnings.
- **FR-004**: No unguarded compiler-specific language extension (other than those already guarded by the existing `#ifdef __GNUC__` / `#if defined(_WIN32)` blocks in `dw.h`) MUST appear in any authored source file.
- **FR-005**: Every POSIX function used in `src/dw_*.c` MUST be defined in POSIX.1-2008 (`_POSIX_C_SOURCE=200809L`) or ISO C99/C11.
- **FR-006**: A compliance acceptance report MUST be written in `specs/013-strict-c-posix-compliance/acceptance-report.md` following the same format as `specs/008-sei-cert-c-compliance/acceptance-report.md`: one section per FR (FR-001 through FR-009) containing a requirement restatement, evidence bullets citing specific files and line ranges, and an explicit ✅ PASS or ❌ FAIL verdict.
- **FR-007**: The SEI CERT C static analysis results MUST be re-run after all changes and the delta against the feature 008 baselines MUST be captured in the audit document; net new violations MUST be zero.
- **FR-008**: The full test suite (`tests/dw-tests.sh`) MUST pass without modification on macOS and Linux after all changes are applied.
- **FR-009**: The timing code in `dw_solver.c` and `dw_support.c` MUST be migrated to `clock_gettime(CLOCK_MONOTONIC, ...)`, guarded by a new `HAVE_CLOCK_GETTIME` autoconf feature check. The vestigial `HAVE_GETTIMEOFDAY` and `HAVE_SYS_TIME_H` configure checks MUST be removed from `configure.ac` (zero call sites exist in any authored source file — confirmed by Phase 0 research). The `#else` fallback for platforms lacking `clock_gettime` MUST use the existing ISO C `time()`/`clock()` path.

### Key Entities

- **Authored Source Files**: `src/dw_*.c` and `src/dw_*.h` — the files subject to this feature.
- **Build Configuration**: `configure.ac`, `src/Makefile.am`, `Makefile.am` — where compiler flags and feature-test macros are declared.
- **Static Analysis Baselines**: `specs/008-sei-cert-c-compliance/baseline-cppcheck.txt` and `specs/008-sei-cert-c-compliance/baseline-warnings.txt` — the reference against which regressions are measured.
- **Compliance Audit Document**: A new file in `specs/013-strict-c-posix-compliance/` that records per-file findings and the pass/fail verdict for each requirement.

---

## Success Criteria

### Measurable Outcomes

- **SC-001**: 100% of `src/dw_*.c` files compile without warnings or errors under `-std=c11 -pedantic-errors` on both GCC and Clang.
- **SC-002**: 100% of `src/dw_*.c` compilations succeed with `-D_POSIX_C_SOURCE=200809L` active and zero implicit-declaration diagnostics for POSIX functions.
- **SC-003**: Zero new findings appear in `cppcheck` or `clang --analyze` output compared to the feature 008 baselines — SEI CERT C compliance does not regress.
- **SC-004**: 100% of tests in `tests/dw-tests.sh` pass after all changes on macOS and Linux.
- **SC-005**: A compliance acceptance report is delivered in `specs/013-strict-c-posix-compliance/acceptance-report.md`, formatted consistently with the feature 008 report, with a ✅ PASS verdict for every FR-001 through FR-009.
