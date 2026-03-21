# Quickstart: Building and Running the New Tests

**Branch**: `006-test-coverage`

This guide covers the four new test-related workflows introduced by this feature:
1. Running the BLAS unit tests locally (`make check`)
2. Running the guidefile parser shell tests
3. Running a full sanitizer build manually
4. Running the three new example problems

All commands are run from the repository root unless noted otherwise.

---

## Prerequisites

```sh
# Standard configure + make (first time or after clean)
autoreconf -fi        # only needed if configure.ac or Makefile.am changed
./configure
make
```

---

## 1. BLAS Unit Tests (`make check`)

Once `tests/Makefile.am` is wired into the build:

```sh
./configure
make check
```

`make check` compiles `tests/test_blas.c` against `src/dw_blas.o` and runs the resulting binary.
Output on success:

```
PASS: test_blas
==================
1 of 1 tests passed
==================
```

On failure the failed assertion is printed to stderr:

```
FAIL [dw_daxpy: len=3, alpha=2.0]: expected 7.000000000000000, got 6.999999999999999
FAIL: test_blas
```

To run the test binary directly (useful for debugging):

```sh
cd tests
./test_blas
```

---

## 2. Guidefile Parser Shell Tests

```sh
cd tests
bash test_guidefile.sh
```

The script runs `dwsolver` with several synthetic guidefiles and checks exit codes and output.
A passing run prints nothing. A failing test prints the test name and reason.

---

## 3. Full Integration Tests (existing)

```sh
cd tests
bash dw-tests.sh
```

This runs all integration tests including the three new examples (`single_sub`, `one_iter`, `neg_y`).
The `dwsolver` binary must be on `PATH` (or the script must be updated to use a relative path).

To add `dwsolver` to PATH from the repo root:

```sh
export PATH="$PWD/src:$PATH"
cd tests
bash dw-tests.sh
```

---

## 4. Manual Sanitizer Build

### ASan + UBSan (address + undefined behaviour)

```sh
# Fix the pre-condition UB first (dw_support.c ~line 427: %d -> %g for DEFAULT_MIP_GAP)
./configure CFLAGS="-fsanitize=address,undefined -g -O1"
make
export PATH="$PWD/src:$PATH"
cd tests && bash dw-tests.sh
```

Any ASan/UBSan violation will print a report to stderr and exit non-zero.

### TSan (thread sanitizer)

```sh
# Must be a separate build — do NOT combine with ASan
make distclean   # or ./configure in a fresh build directory
./configure CFLAGS="-fsanitize=thread -g -O1"
make
export PATH="$PWD/src:$PATH"
cd tests && bash dw-tests.sh
```

### Clean up after sanitizer build

```sh
make distclean
./configure      # back to normal build
make
```

---

## 5. Running the New Example Problems by Hand

Each example is a subdirectory of `examples/`. To run dwsolver on a new example:

```sh
export PATH="$PWD/src:$PATH"

# single-subproblem example
dwsolver examples/single_sub/guidefile

# one-Phase-II-iteration example
dwsolver examples/one_iter/guidefile

# y-accumulator sign-correction example
dwsolver examples/neg_y/guidefile
```

To check the relaxed LP solution (flag `-c` prints sorted column values):

```sh
dwsolver -c examples/single_sub/guidefile | sort > /tmp/rs_sorted
diff examples/single_sub/ex_relaxed_solution /tmp/rs_sorted
```

---

## 6. CI Reference

After this feature is merged, the GitHub Actions CI matrix will include four Linux jobs:

| Job | Trigger |
|-----|---------|
| `Linux (GCC)` | push / PR — existing |
| `Linux (ASan+UBSan)` | push / PR — new |
| `Linux (TSan)` | push / PR — new |
| `macOS`, `Windows`, `Docker` | unchanged |

All jobs must pass before any PR can be merged to `main`.
