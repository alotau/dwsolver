# Quickstart: Using libdwsolver

**Feature**: 010-callable-library  
**Date**: 2026-03-22

This guide shows how to integrate `libdwsolver` into a C project after it has been
installed via `make install` (or `make install DESTDIR=...`).

---

## Prerequisites

- `libdwsolver` installed (see Build & Install below)
- GCC or Clang (C99 or later)
- `pkg-config` (recommended; manual flags are an alternative)

---

## Build & Install

```sh
# From the dwsolver source root:
./configure --prefix=/usr/local
make
make install
```

This installs:
- `/usr/local/bin/dwsolver` — the CLI (unchanged behavior)
- `/usr/local/lib/libdwsolver.a` — static library
- `/usr/local/lib/libdwsolver.so` (Linux) or `libdwsolver.dylib` (macOS) — shared library
- `/usr/local/include/dwsolver.h` — public API header
- `/usr/local/lib/pkgconfig/dwsolver.pc` — pkg-config metadata

---

## Minimal Example

```c
/* solve_example.c */
#include <stdio.h>
#include "dwsolver.h"

int main(void) {
    dw_options_t opts;
    dw_result_t  result = {0};

    /* Initialize options to defaults */
    dw_options_init(&opts);

    /* Optionally override defaults */
    opts.verbosity = DW_OUTPUT_QUIET;

    /* Subproblem file list */
    const char *subproblems[] = {
        "examples/book_bertsimas/sub1.lp",
        "examples/book_bertsimas/sub2.lp"
    };

    /* Run the solver */
    int rc = dw_solve(
        "examples/book_bertsimas/master.lp",
        subproblems,
        2,      /* num_subproblems */
        &opts,
        &result
    );

    if (rc == DW_STATUS_OK) {
        printf("Optimal objective: %.6f\n", result.objective_value);
        printf("Solution has %d variables.\n", result.num_vars);
        for (int i = 0; i < result.num_vars && i < 5; i++)
            printf("  x[%d] = %.6f\n", i, result.x[i]);
    } else {
        fprintf(stderr, "Solver failed with status %d\n", rc);
    }

    /* Always free result resources */
    dw_result_free(&result);
    return rc;
}
```

---

## Compiling with pkg-config (recommended)

```sh
cc $(pkg-config --cflags dwsolver) -o solve_example solve_example.c \
   $(pkg-config --libs dwsolver)
```

---

## Compiling without pkg-config

```sh
cc -I/usr/local/include -o solve_example solve_example.c \
   -L/usr/local/lib -ldwsolver -lpthread -lm
```

---

## Verifying the Install

```sh
# Check CLI still works:
dwsolver --version

# Check pkg-config:
pkg-config --modversion dwsolver
pkg-config --cflags --libs dwsolver

# Check exported symbols (should show only dw_options_init, dw_solve, dw_result_free, dw_version):
nm -D /usr/local/lib/libdwsolver.so | grep ' T dw_'
```

---

## Common Option Overrides

```c
/* Run silently */
opts.verbosity = DW_OUTPUT_SILENT;

/* Tighten MIP gap to 0.1% */
opts.mip_gap = 0.001;

/* Increase phase 2 iteration budget */
opts.max_phase2_iterations = 10000;

/* Request rounding after LP relaxation */
opts.rounding_flag = 1;
```

---

## Thread Safety Note

`dw_solve()` is **not reentrant**. If multiple threads need to run solves, serialize
calls with a mutex:

```c
pthread_mutex_lock(&my_solver_mutex);
int rc = dw_solve(master, subs, n, &opts, &result);
pthread_mutex_unlock(&my_solver_mutex);
```

---

## Error Handling

```c
int rc = dw_solve(master, subs, n, &opts, &result);
switch (rc) {
    case DW_STATUS_OK:            /* success */           break;
    case DW_STATUS_ERR_FILE:      /* bad file path */     break;
    case DW_STATUS_ERR_INFEASIBLE:/* infeasible problem */break;
    case DW_STATUS_ERR_PHASE1_LIMIT:
    case DW_STATUS_ERR_PHASE2_LIMIT:
        /* iteration limit hit; partial result may exist */
        break;
    default:
        fprintf(stderr, "Unexpected solver error: %d\n", rc);
}
dw_result_free(&result);  /* always safe to call, even on error */
```

---

## Linking in a Makefile.am (Autotools project)

```makefile
bin_PROGRAMS = my_app
my_app_SOURCES = main.c
my_app_CFLAGS  = $(shell pkg-config --cflags dwsolver)
my_app_LDADD   = $(shell pkg-config --libs dwsolver)
```

Or using `PKG_CHECK_MODULES` in your `configure.ac`:

```
PKG_CHECK_MODULES([DWSOLVER], [dwsolver >= 1.2])
```

Then in `Makefile.am`:
```makefile
my_app_CFLAGS = $(DWSOLVER_CFLAGS)
my_app_LDADD  = $(DWSOLVER_LIBS)
```

---

## Linking in a CMakeLists.txt

```cmake
find_package(PkgConfig REQUIRED)
pkg_check_modules(DWSOLVER REQUIRED dwsolver)

target_include_directories(my_app PRIVATE ${DWSOLVER_INCLUDE_DIRS})
target_link_libraries(my_app PRIVATE ${DWSOLVER_LIBRARIES})
```
