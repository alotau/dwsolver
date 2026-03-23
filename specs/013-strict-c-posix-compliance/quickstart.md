# Quickstart: Build and Verify ISO C / POSIX.1 Compliance

This guide walks a developer through building dwsolver with strict standard enforcement
and verifying that the compliance requirements of feature 013 are met.

---

## Prerequisites

- GCC ≥ 4.7 or Clang ≥ 3.1
- GNU Autotools (autoconf, automake, libtool)
- System GLPK ≥ 4.65 (`libglpk-dev` on apt, `glpk-devel` on dnf, `brew install glpk` on macOS)
- `cppcheck` (optional, for static analysis verification)

---

## Step 1: Configure a clean build

```sh
./configure
```

After this feature is implemented, `config.h` will contain:

```c
#define HAVE_CLOCK_GETTIME 1
```

on any POSIX.1-2008 host (macOS 10.12+, Linux with glibc ≥ 2.17, or musl).

---

## Step 2: Build with strict flags

The updated `src/Makefile.am` sets `AM_CFLAGS = -std=c11 -pedantic-errors` and
`AM_CPPFLAGS = -D_POSIX_C_SOURCE=200809L`, so a plain `make` is sufficient:

```sh
make
```

A successful build with zero warnings or errors is the acceptance criterion for FR-001,
FR-002, and FR-003.

### Expected compiler invocation (illustrative)

```sh
gcc -std=c11 -pedantic-errors -D_POSIX_C_SOURCE=200809L \
    -DDWSOLVER_BUILDING_LIB -I. -I/usr/include \
    -c dw_solver.c -o dw_solver.o
```

### Common errors and fixes

| Error | Cause | Fix |
|-------|-------|-----|
| `implicit declaration of function 'clock_gettime'` | `_POSIX_C_SOURCE` not set or host predates POSIX.1-2008 | Ensure `AM_CPPFLAGS` contains `-D_POSIX_C_SOURCE=200809L`; check `config.h` |
| `error: ISO C forbids ...` | Unguarded GCC extension | Add `#ifdef __GNUC__` guard or remove the extension |
| `undefined reference to 'clock_gettime'` | Linux needs `-lrt` | `AC_SEARCH_LIBS([clock_gettime], [rt])` in `configure.ac` handles this automatically |

---

## Step 3: Run the test suite

```sh
make check
```

The BLAS unit test (`test_blas`) and all suites in `tests/dw-tests.sh` must pass. This is
the acceptance criterion for FR-008 and SC-004.

> **macOS note**: `test_lib_api` fails on macOS with `sem_init: Function not implemented`
> (`ENOSYS`). This is a pre-existing macOS limitation (unnamed POSIX semaphores are not
> supported on macOS; named semaphores via `sem_open` work correctly). This failure predates
> this feature and is unrelated to POSIX compliance. Linux CI is the authoritative pass
> criterion for FR-008.

---

## Step 4: Verify POSIX compliance with a strict build (optional)

To manually simulate a POSIX-strict build without GNU extensions:

```sh
make CFLAGS="-std=c11 -pedantic-errors -D_POSIX_C_SOURCE=200809L -Wall -Wextra"
```

Zero warnings or errors for `src/dw_*.c` is the criterion for SC-001 and SC-002.

---

## Step 5: Run static analysis (FR-007, SC-003)

### cppcheck

```sh
cppcheck --enable=all --suppress=missingIncludeSystem \
         -I src/ src/dw_*.c 2>&1 | tee /tmp/cppcheck-013.txt
diff specs/008-sei-cert-c-compliance/baseline-cppcheck.txt /tmp/cppcheck-013.txt
```

No new lines in the diff output (additions) is the pass criterion.

### clang static analyser

```sh
clang --analyze -std=c11 -D_POSIX_C_SOURCE=200809L \
      -I src/ src/dw_*.c 2>&1 | tee /tmp/clang-analyze-013.txt
diff specs/008-sei-cert-c-compliance/baseline-warnings.txt /tmp/clang-analyze-013.txt
```

---

## Step 6: Windows build (CI verification)

The Windows CI job (`ci-windows.yml`) runs on MSYS2 MINGW64. Push the branch to trigger it,
or inspect it locally with `act` if configured. The job uses `continue-on-error: true`,
so a Windows failure is non-blocking but should be investigated.

`_POSIX_C_SOURCE=200809L` is recognised by MinGW-w64 headers. The existing
`#if defined(_WIN32)` guards in `dw.h` must continue to compile cleanly.

---

## Summary checklist

- [x] `./configure && make` — zero warnings/errors for `src/dw_*.c` ⟹ FR-001, FR-002, FR-003
- [x] `make check`: `test_blas` PASS; `test_lib_api` FAIL is pre-existing macOS limitation (see Step 3 note) ⟹ FR-008
- [x] `cppcheck` diff has zero new additions (`specs/013-strict-c-posix-compliance/cppcheck-post.txt`) ⟹ FR-007, SC-003
- [x] `clang --analyze` diff has zero new additions (`specs/013-strict-c-posix-compliance/clang-analyze-post.txt`) ⟹ FR-007, SC-003
- [x] `acceptance-report.md` written with ✅ PASS per FR ⟹ FR-006, SC-005
