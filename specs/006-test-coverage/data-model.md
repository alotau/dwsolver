# Data Model: Structured Test Coverage for dwsolver

**Branch**: `006-test-coverage`
**Date**: 2026-03-20

This feature does not introduce persistent data structures or storage schemas.
The "data model" here describes the structural entities that must be created and
their invariants, so implementers have a precise checklist.

---

## Entity: BLAS Test Case (`tests/test_blas.c`)

| Field | Type | Description |
|-------|------|-------------|
| `name` | `const char *` | Human-readable ID printed on failure, e.g. `"dw_daxpy: len=3, alpha=2.0"` |
| inputs | typed C values | Passed directly to the function under test |
| `expected` | `double` or `double[]` | Known-correct output |
| `tolerance` | `double` | Absolute comparison epsilon — use `1e-12` for all BLAS tests |

**Invariants**:
- Every function exported from `dw_blas.h` must appear in at least one test case.
- The test `name` must contain the function name so failures are diagnosable without source lookup.
- On failure: print `FAIL [name]: expected <x>, got <y>` to stderr, then call `exit(1)`.
- On success: print nothing (or a single trailing `ALL TESTS PASSED` line) and exit 0.

**Macro interface** (defined at the top of `tests/test_blas.c`):

```c
/* Passes if |got - expected| <= tol */
#define EXPECT_APPROX(name, expected, got, tol) \
    do { \
        double _e = (expected), _g = (got), _d = _e - _g; \
        if (_d < -(tol) || _d > (tol)) { \
            fprintf(stderr, "FAIL [%s]: expected %.15g, got %.15g\n", \
                    (name), _e, _g); \
            exit(1); \
        } \
    } while(0)

/* Passes if got == expected (integer comparison) */
#define EXPECT_INT(name, expected, got) \
    do { \
        if ((expected) != (got)) { \
            fprintf(stderr, "FAIL [%s]: expected %d, got %d\n", \
                    (name), (int)(expected), (int)(got)); \
            exit(1); \
        } \
    } while(0)
```

**Build structure** (`tests/Makefile.am`):

```makefile
TESTS           = test_blas
check_PROGRAMS  = test_blas
test_blas_SOURCES = test_blas.c
test_blas_CFLAGS  = -I$(top_srcdir)/src
test_blas_LDADD   = $(top_builddir)/src/dw_blas.o
```

`configure.ac` change required: add `tests/Makefile` to `AC_CONFIG_FILES`.

---

## Entity: Example Problem Directory (`examples/<name>/`)

| File | Required | Description |
|------|----------|-------------|
| `guidefile` | yes | dwsolver guidefile: line 1 = `n`, lines 2..n+1 = subproblem .cplex paths, line n+2 = master path, line n+3 = monolithic path, optional line n+4 = objective constant |
| `sub<k>.cplex` (k = 1..n) | yes | Subproblem LP in CPLEX LP format |
| `master.cplex` | yes | Master LP in CPLEX LP format |
| `dw_monolithic.cplex` | yes | Union LP — same optimal as DW result |
| `ex_relaxed_solution` | for diff-based tests | Sorted expected output of `dwsolver -c` |
| `README` | recommended | Problem description, known optimal, code path exercised |

**Invariants**:
- Variable names must be globally consistent across all `.cplex` files in the same directory.
- The monolithic LP optimal must match the Dantzig-Wolfe relaxation optimal.
- The `guidefile` integer `n` must equal the exact number of subproblem files listed.
- Each new example must be registered in `tests/dw-tests.sh` with a `pushd`/`popd` block.

**New examples summary**:

| Directory | n | Code path exercised | Test method |
|-----------|---|---------------------|-------------|
| `single_sub` | 1 | `num_clients=1` semaphore round-trip | `diff ex_relaxed_solution` |
| `one_iter` | 2 | `iteration_count==1` bail in Phase II | `grep` objective value |
| `neg_y` | 2 | y-accumulator sign correction (`val=-1.0`) in Phase I | `diff ex_relaxed_solution` |

---

## Entity: Sanitizer CI Job (`.github/workflows/`)

| Field | ASan+UBSan job | TSan job |
|-------|---------------|----------|
| `jobs.<id>.name` | `Linux (ASan+UBSan)` | `Linux (TSan)` |
| `runs-on` | `ubuntu-latest` | `ubuntu-latest` |
| `CFLAGS` env var | `-fsanitize=address,undefined -g -O1` | `-fsanitize=thread -g -O1` |
| autotools touch step | yes (copy from `ci-linux.yml`) | yes (copy from `ci-linux.yml`) |
| configure step | `./configure CFLAGS="$CFLAGS"` | `./configure CFLAGS="$CFLAGS"` |
| make step | `make` | `make` |
| test step | `cd tests && bash dw-tests.sh` | `cd tests && bash dw-tests.sh` |

**Invariants**:
- ASan and TSan MUST be in separate jobs — combining them in one binary is unsupported by GCC/Clang.
- The pre-condition UB fix (`%d` → `%g` for `DEFAULT_MIP_GAP` in `dw_support.c`) MUST land before
  these jobs are added or UBSan will fire immediately.
- macOS and Windows CI jobs are NOT extended with sanitizers (Apple Clang TSan + named POSIX
  semaphores is unreliable; MSVC has no `-fsanitize` support).

---

## Entity: Guidefile Parser Test Script (`tests/test_guidefile.sh`)

| Test case | Mechanism | Expected result |
|-----------|-----------|-----------------|
| Valid 1-subproblem guidefile | Run `dwsolver guidefile` against `single_sub` | Exit 0 or 1; no `USAGE:` on stderr |
| Valid 4-subproblem guidefile | Run against an existing example | Exit 0 or 1; no parse-error text |
| Guidefile declares n=3, lists 2 files | Synthetic bad guidefile | Non-zero exit; informative error message |
| Missing `.cplex` file | Guidefile lists nonexistent path | Non-zero exit; no crash / no SIGSEGV |

**Note**: The solver may exit non-zero for LP infeasibility even on "valid" synthetic inputs. The
test checks exit code relative to parse errors only — and the absence of a crash (non-zero from
`SIGSEGV`/`SIGABRT` is distinguishable in bash via `$?` ≥ 128).
