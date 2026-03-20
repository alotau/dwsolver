# Quickstart: Build, Test, and Repair DWSOLVER

**Feature**: `002-cross-platform-repair`
**Date**: 2026-03-19

---

## Prerequisites

### macOS

```bash
# Xcode Command Line Tools (provides Clang, make, autotools)
xcode-select --install
# Or via Homebrew:
brew install autoconf automake libtool
```

### Linux (Ubuntu / Debian)

```bash
sudo apt-get install build-essential autoconf automake libtool
```

### Windows

Windows support is out of scope for this repair spec. See future spec for MinGW / WSL approach.

---

## Build

### macOS (named semaphores required)

```bash
cd /path/to/dwsolver-repaired
./configure --enable-named-semaphores
make
# Binary: src/dwsolver
```

### Linux (unnamed semaphores — default)

```bash
cd /path/to/dwsolver-repaired
./configure
make
# Binary: src/dwsolver
```

Optional configure flags:

| Flag | Purpose |
|------|---------|
| `--enable-named-semaphores` | Required on macOS; safe to use on Linux too |
| `--enable-recursive-mutex` | Use `PTHREAD_MUTEX_RECURSIVE` on macOS |
| `--with-gmp` | Enable GMP integer arithmetic in GLPK |
| `--with-zlib` | Enable compressed MPS file support in GLPK |

---

## Run All Tests

After a successful build:

```bash
export PATH="$PWD/src:$PATH"   # put dwsolver on PATH
bash tests/dw-tests.sh
```

Expected output (after this repair is complete):

```
PASS: book_bertsimas
PASS: book_lasdon
PASS: web_mitchell
PASS: web_trick
PASS: four_sea
PASS: book_dantzig
```

Any `FAIL` line indicates a regression. Investigate before merging.

---

## Run a Single Example Manually

Examples must be run from their own directory because guide files reference LP files
by relative path.

```bash
cd examples/book_bertsimas
../../src/dwsolver book_bertsimas.dw
cat relaxed_solution
```

Expected result for `book_bertsimas`:
```
x1 = 2.0
x2 = 1.5
x3 = 2.0
```

---

## Verify the KD-001 Fix (Linux)

If you have a Linux environment (Docker, VM, or CI):

```bash
# On Linux, standard configure (no named-semaphores flag)
./configure && make 2>&1 | grep -E 'error:|warning:|^make'
```

Before the fix, you will see:
```
multiple definition of `master_lp_ready_mutex`; dw_support.o: first defined here
```

After the fix, the build completes cleanly (zero linker errors).

---

## Verify the KD-001 Fix Does Not Regress macOS

```bash
# macOS: use existing flags
./configure --enable-named-semaphores && make && bash tests/dw-tests.sh
```

All 6 tests should pass.

---

## Development Tips

- **Rebuild after header changes**: `make clean && make` — touching `dw.h` does not
  always trigger recompilation of all dependent `.c` files under `make`.
- **Check for remaining `sprintf`**: `grep -n 'sprintf(' src/dw_*.c` — should return
  zero results after the KD-004 fix.
- **Check for unchecked `malloc`**: After the KD-005 fix, each critical `malloc` call
  is immediately followed by `dw_oom_abort(ptr, "context")`.
- **Thread sanitizer** (optional, when tooling supports it):
  ```bash
  CFLAGS="-fsanitize=thread -g" ./configure --enable-named-semaphores && make
  cd examples/book_bertsimas && ../../src/dwsolver book_bertsimas.dw
  ```
