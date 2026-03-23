# DWSOLVER (Dantzig-Wolfe Solver)

[![macOS](https://github.com/alotau/dwsolver/actions/workflows/ci-macos.yml/badge.svg)](https://github.com/alotau/dwsolver/actions/workflows/ci-macos.yml)
[![Linux](https://github.com/alotau/dwsolver/actions/workflows/ci-linux.yml/badge.svg)](https://github.com/alotau/dwsolver/actions/workflows/ci-linux.yml)
[![Windows](https://github.com/alotau/dwsolver/actions/workflows/ci-windows.yml/badge.svg)](https://github.com/alotau/dwsolver/actions/workflows/ci-windows.yml)
[![Docker](https://github.com/alotau/dwsolver/actions/workflows/ci-docker.yml/badge.svg)](https://github.com/alotau/dwsolver/actions/workflows/ci-docker.yml)

## Overview
The Dantzig-Wolfe Solver program is a stand-alone implementation of the
Dantzig-Wolfe Decomposition algorithm.  The GNU Linear Programming Kit provides
the functions for all of the necessary linear programming (reading in problems,
performing the simplex algorithm, querying various LP data strucures, etc.).
This software also uses POSIX Threads (pthreads) for parallel solving of 
subproblems.  The implementation is general in that any linear program that can
be presented to the software in block-angular form (with some current 
limitations) can be solved.

## Dependencies

dwsolver requires **GLPK ≥ 4.65** (the GNU Linear Programming Kit) at build
time and at runtime.  Install it from your OS package manager before running
`./configure`:

| Platform | Command |
|----------|---------|
| macOS (Homebrew) | `brew install glpk` |
| Ubuntu / Debian | `sudo apt-get install libglpk-dev` |
| Fedora / RHEL | `sudo dnf install glpk-devel` |
| Windows (MSYS2/MINGW64) | `pacman -S mingw-w64-x86_64-glpk` |

> **ABI note for libdwsolver callers**: programs that link against
> `libdwsolver` must also link — directly or transitively — against a GLPK
> shared library with the same **soname** as the one used at build time.
> GLPK ≥ 4.60 through the current 5.x series all ship as `libglpk.so.40`
> (Linux) / `libglpk.40.dylib` (macOS), so any GLPK ≥ 4.65 installation is
> ABI-compatible with any other at runtime; an exact version match is **not**
> required.  When dynamically linking, `libdwsolver` records this dependency
> via `DT_NEEDED`/`LC_LOAD_DYLIB`, so you normally do not need to list
> `-lglpk` explicitly.  For static linking, use
> `pkg-config --libs --static dwsolver` to obtain all required libraries
> including GLPK.

## Building

### macOS and Linux

Prerequisites: `gcc`, `make`, `autoconf`, `automake`, and GLPK (via your package
manager, e.g. `brew install glpk` or `apt-get install libglpk-dev`).

```sh
./configure && make
```

The compiled binary is `src/dwsolver`.

### Windows (MSYS2 / MINGW64)

1. Install [MSYS2](https://www.msys2.org/) and open an **MINGW64** shell.
2. Install the required packages:
   ```sh
   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-winpthreads-git make automake diffutils
   ```
3. Build:
   ```sh
   ./configure && make
   ```

The compiled binary is `src/dwsolver.exe`.

> **Note:** Named POSIX semaphores are not supported on Windows. The solver uses
> unnamed semaphores automatically on that platform; no extra configure flag is
> needed.

### Docker

A pre-built image is available, or you can build one locally:

```sh
# Build the image
docker build -t dwsolver .

# Run the solver (mount a directory containing your .lp files)
docker run --rm -v "$PWD/examples/book_dantzig:/data" dwsolver book_dantzig.lp
```

The container's working directory is `/data`.  Any `.lp` input files and output
files live there; mount the directory that holds your problem files.

## SEI CERT C Compliance

The solver's C source code has been audited and remediated against the [SEI CERT C Coding Standard](https://wiki.sei.cmu.edu/confluence/display/c/SEI+CERT+C+Coding+Standard), focusing on the concurrency and POSIX rules most relevant to a multi-threaded solver:

| Rule | Description | Status |
|------|-------------|--------|
| POS54-C | Detect and handle POSIX errors | ✅ All `pthread_*` and `sem_*` call sites checked via `DW_PTHREAD_CHECK` / `DW_SEM_CHECK` |
| MEM35-C | Allocate sufficient memory | ✅ `malloc` operator-precedence bugs corrected |
| CON31-C / CON43-C | Thread and lock safety | ✅ Mutex acquisition order documented and enforced; no lock inversions |
| EXP12-C | Do not ignore operands in expressions | ✅ `errno` snapshotted before `strerror()` / `fprintf()` calls |
| MSC24-C | Do not use deprecated functions | ✅ `sem_open` return value validated against `SEM_FAILED` |

**Artifacts**:
- [Acceptance report](specs/008-sei-cert-c-compliance/acceptance-report.md) — all 10 functional requirements with evidence
- [Compliance spec](specs/008-sei-cert-c-compliance/spec.md) — full rule set and rationale
- [Static analysis baseline](specs/008-sei-cert-c-compliance/baseline-warnings.txt) — zero new warnings from remediation

## Architecture

The [`architecture/`](./architecture/) directory contains Mermaid diagrams covering the threading model, module dependencies, and input/output data flow. These render natively in GitHub without any plugins.

The [`build-aux/`](./build-aux/) directory contains Autoconf/Automake auxiliary scripts (`config.guess`, `config.sub`, `ltmain.sh`, etc.) generated by `autoreconf`. These are committed so the build works without requiring a specific Automake version to be installed.

The [`third-party/glpk/`](./third-party/glpk/) directory contains GLPK 4.44
attribution and patch files bundled with the original dwsolver source
distribution.  dwsolver now links against a system-installed GLPK (≥ 4.65)
rather than carrying its own copy.

## Callable Library (libdwsolver)

In addition to the `dwsolver` CLI, the build also produces **`libdwsolver`** — a
C library that lets you call the solver directly from your own code.

### What gets installed

```
/usr/local/bin/dwsolver          # CLI (unchanged)
/usr/local/lib/libdwsolver.a     # static library
/usr/local/include/dw_solver.h   # public API header
/usr/local/lib/pkgconfig/dwsolver.pc  # pkg-config metadata
```

### Quick example

```c
#include <stdio.h>
#include "dw_solver.h"

int main(void) {
    const char *subs[] = { "sub1.lp", "sub2.lp" };
    dw_options_t opts;
    dw_result_t  result = {0};

    dw_options_init(&opts);          /* set all options to defaults */

    int rc = dw_solve("master.lp", subs, 2, &opts, &result);
    if (rc == DW_STATUS_OK)
        printf("Objective: %g\n", result.objective_value);

    dw_result_free(&result);
    return rc;
}
```

Compile with pkg-config:

```sh
gcc -o myapp myapp.c $(pkg-config --cflags --libs dwsolver)
```

### Public API

| Symbol | Description |
|--------|-------------|
| `dw_options_init(opts)` | Populate `dw_options_t` with documented defaults |
| `dw_solve(master, subs, n, opts, result)` | Run the solver; returns a `dw_status_t` code |
| `dw_result_free(result)` | Release heap memory owned by a `dw_result_t` |
| `dw_version()` | Return a static version string |

Status codes: `DW_STATUS_OK`, `DW_STATUS_ERR_BAD_ARGS`, `DW_STATUS_ERR_FILE`,
`DW_STATUS_ERR_INFEASIBLE`, `DW_STATUS_ERR_PHASE1_LIMIT`,
`DW_STATUS_ERR_PHASE2_LIMIT`, `DW_STATUS_ERR_UNBOUNDED`, `DW_STATUS_ERR_INTERNAL`.

See [`specs/010-callable-library/quickstart.md`](specs/010-callable-library/quickstart.md) for a fuller walkthrough
and [`specs/010-callable-library/contracts/api.md`](specs/010-callable-library/contracts/api.md) for the full API contract.

## Examples and Tests
There is a an [examples](./examples/) directory with several problems plucked from popular textbooks and some websites. There is also a toy version of the problem that initiated the writing of this code related to air traffic management. Each example has a README and should run successfully with a properly built dwsolver executable. There is also a new [tests](./tests/) directory that contains a [test script](./tests/dw-tests.sh) that runs most of the examples.  A couple of the examples have non-deterministic solutions (but deterministic optimum values), and I didn't take the time to write tests that check for the correct optimum for those problems. Again, happy for any PR that runs those examples, parses the output files to pluck out the optimum, and shows it is the expected value.

## Copyright
Copyright  2010 United States Government National Aeronautics and Space 
Administration (NASA).  No copyright is claimed in the United States under 
Title 17, U.S. Code. All Other Rights Reserved.

## Legal
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

In accordance with the GNU General Public License Version 3, 29 June 2007
(GPL V3) Section 7. Additional Terms, Additional Permissions are added as
exceptions to the terms of the GPL V3 for this program.  These additional
terms should have been received with this program in a file entitled
["ADDITIONAL_LICENSE_TERMS"](./ADDITIONAL_LICENSE_TERMS).  If a copy was not provided, you may request
one from the contact author listed below.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Additional Info
See the file [COPYING](COPYING) for the GNU Public License.

All Dantzig-Wolfe source files created by the author begins with the prefix
"dw" while all other source files in the src/ directory are GLPK files.  All
GLPK files are provided as published by the author(s) of GLPK except for those
specified in the .patch file provided with this release.  The changes to GLPK
were made to implement a thread-friendly version of GLPK necessary for 
implementation of a parallel Dantzig-Wolfe algorithm.

~~Please report bugs/comments/suggestions/patches to Joseph.L.Rios@nasa.gov.~~
I'm not longer with NASA, so best to just submit through GitHub as an issue.
