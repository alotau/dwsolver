# Feature Specification: Callable Library Interface

**Feature Branch**: `010-callable-library`  
**Created**: 2026-03-21  
**Status**: Draft  
**Input**: User description: "Create a callable library from this project. We will keep the command line interface, but add appropriate code and artifacts to have the codebase produce a callable library as well."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Embed Solver in a Host Application (Priority: P1)

A developer building a larger optimization or decision-support system wants to call the Dantzig-Wolfe solver directly from their own C or C++ application without spawning a subprocess or parsing files round-trip through the command line. They link against `libdwsolver` and call a documented API to configure, run, and retrieve the solution.

**Why this priority**: This is the core value of the feature — without a linkable library API, no embedding is possible. All other stories build on this capability.

**Independent Test**: Link a minimal C test program against `libdwsolver`, call the solver API with one of the existing example problem files, and verify the returned solution values match those produced by the `dwsolver` CLI on the same inputs.

**Acceptance Scenarios**:

1. **Given** `libdwsolver` is installed on the system, **When** a C program `#include`s `dwsolver.h` and links with `-ldwsolver`, **Then** the program compiles and links without errors.
2. **Given** a valid master LP file and one or more subproblem LP files, **When** the host program calls the library solve function with those file paths and default options, **Then** the function returns a status code indicating success and populates a solution buffer with primal variable values.
3. **Given** the host program runs the solver on the `book_bertsimas` example, **When** the returned solution is compared to the `dwsolver` CLI output for the same example, **Then** the objective values agree within the configured MIP gap tolerance.
4. **Given** the host program calls the library with invalid file paths, **When** the solve function is invoked, **Then** it returns a non-zero error code and does not crash or produce undefined behavior.

---

### User Story 2 - Build System Produces Both CLI and Library (Priority: P2)

A developer building or packaging the project runs the standard `./configure && make && make install` workflow and receives both the `dwsolver` executable and the `libdwsolver` static and shared libraries, along with the public header, in the expected install locations.

**Why this priority**: Delivery of both artifacts from one build is what guarantees the CLI is preserved while the library is added. Without this, one or the other is missing.

**Independent Test**: Run the build and install steps from a clean checkout; verify that `bin/dwsolver`, `lib/libdwsolver.a`, `lib/libdwsolver.so` (or `.dylib`), and `include/dwsolver.h` are all present in the install prefix.

**Acceptance Scenarios**:

1. **Given** a clean source checkout, **When** the standard autotools configure-make-install sequence completes, **Then** the install prefix contains the `dwsolver` binary, both library variants (`libdwsolver.a` and a shared library), and the `dwsolver.h` public header.
2. **Given** the installed shared library, **When** a consumer links against it at runtime, **Then** the shared library's exported symbols are limited to those declared in `dwsolver.h` (no internal GLPK or solver implementation symbols are inadvertently exported).
3. **Given** the `dwsolver` executable is built after the library refactor, **When** the existing test suite (`tests/dw-tests.sh`) is run, **Then** all tests that previously passed continue to pass without modification.

---

### User Story 3 - Configure Solver Options Programmatically (Priority: P3)

A developer wants to tune solver behavior — verbosity, MIP gap tolerance, iteration limits, rounding — without constructing command-line argument strings or writing configuration files. They use an options struct provided by the library API to set parameters before calling the solver.

**Why this priority**: Programmatic option control is what distinguishes a real library API from a thin wrapper that still parses argv. It is lower priority because file-path-based invocation with defaults already delivers usable value.

**Independent Test**: Write a host program that sets non-default options (e.g., silent verbosity, custom MIP gap) via the options struct, runs the solver on a known example, and verifies the output matches expectations for those settings (e.g., no diagnostic output is produced when verbosity is silent).

**Acceptance Scenarios**:

1. **Given** an initialized options struct, **When** the caller sets verbosity to silent and invokes the solver, **Then** the library produces no output to stdout or stderr during the solve.
2. **Given** an options struct with a custom MIP gap tolerance, **When** the solver runs, **Then** the stopping criterion reflects the caller-provided gap rather than the compiled-in default.
3. **Given** a caller who initializes the options struct via `dw_options_init()` and then overrides only selected fields, **When** the solver is invoked, **Then** fields not explicitly overridden retain the values set by `dw_options_init()` and the solver produces correct results.

---

### Edge Cases

- What happens when the caller invokes the library solver function concurrently from multiple threads without external serialization? The library uses shared global state inherited from the original implementation; the spec should document whether concurrent calls are supported or explicitly prohibited.
- What happens when the host application has already initialized GLPK independently before calling the library? The library must document any GLPK state requirements and avoid double-initialization crashes.
- What happens when the solve is called with zero subproblems or an empty master LP file? The library must return a well-defined error code rather than crashing.
- What happens on platforms where named POSIX semaphores behave differently (e.g., macOS vs. Linux)? The feature inherits existing platform-conditional code; the library artifact must build and run correctly on both platforms.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The build system MUST produce both `libdwsolver` (static archive and shared library) and the `dwsolver` binary from a single build invocation without requiring separate configuration flags.
- **FR-002**: A public C header file (`dwsolver.h`) MUST be defined and installed, exposing the complete callable API and all types needed to use it; no internal headers need be included by the caller.
- **FR-003**: The public API MUST provide a function to initialize a solver options struct to well-documented default values, mirroring the defaults currently applied by the CLI's `init_globals()`.
- **FR-004**: The public API MUST provide a solve function accepting the master LP file path, an array of subproblem LP file paths, a subproblem count, and a pointer to a solver options struct; the function MUST return a numeric status code indicating success or a specific category of failure.
- **FR-005**: The public API MUST provide a mechanism for the caller to retrieve the primal solution values (the final `x` vector) and the optimal objective value after a successful solve.
- **FR-006**: The public API MUST provide a cleanup/free function for any resources allocated by the library on behalf of the caller, preventing memory leaks in long-running host applications.
- **FR-007**: The `dw_main.c` CLI entry point MUST be refactored to call the library API rather than duplicating solver logic, so that the CLI and library share a single code path.
- **FR-008**: All symbols in the compiled library that are not declared in `dwsolver.h` MUST be hidden from the shared library's exported symbol table (using visibility controls or an export map).
- **FR-009**: The existing test suite MUST continue to pass without modification after the refactor.
- **FR-010**: The library MUST be documented — in `dwsolver.h` comments — with whether concurrent calls from multiple threads using the same solver context are supported; if not supported, callers MUST be warned at the API level.
- **FR-011**: A `pkg-config` metadata file (`dwsolver.pc`) MUST be generated and installed so that consumers can query compile and link flags using standard toolchain integration.

### Key Entities

- **Solver Options**: A struct grouping all tunable parameters (verbosity, MIP gap, iteration limits, rounding flags, output preferences); populated to defaults by the init function and passed by pointer to the solve function.
- **Solver Result**: A struct or out-parameter set populated by the solve function containing the objective value, the primal solution vector, and solution status.
- **Public Header (`dwsolver.h`)**: The single include file that defines the API; the boundary between library internals and consumers.
- **`libdwsolver`**: The compiled library artifact (both `.a` and shared variants) containing all solver logic except `main()`.
- **`dwsolver` Binary**: The existing CLI, now acting as a thin wrapper that parses command-line arguments and delegates to the library API.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The full build (configure, make) completes without errors on all platforms currently supported by the project (Linux and macOS), and produces all three artifacts: the binary, the static library, and the shared library.
- **SC-002**: A host program using only `#include "dwsolver.h"` and linking `-ldwsolver` can solve any of the existing example problems and obtain results that match the CLI output to within the MIP gap tolerance (default 1%).
- **SC-003**: The existing regression test suite passes at 100% after the refactor — no previously passing tests regress.
- **SC-004**: The shared library's exported symbol list contains only symbols declared in `dwsolver.h` — zero unintended symbol leakage from GLPK internals or private solver functions.
- **SC-005**: A host program that calls the library's cleanup function after each solve does not leak memory across repeated solve calls, as verified by a memory-checking tool run on a simple loop of 10 consecutive solves.
- **SC-006**: `pkg-config --cflags --libs dwsolver` produces correct compiler and linker flags on a system where the library is installed, enabling downstream build systems to integrate without manual flag management.

## Assumptions

- The library will initially support C callers; future language bindings (Python, etc.) are out of scope for this feature but the API design should not preclude them.
- The current GLPK fork (patched for thread readiness) is considered a private dependency bundled with the library; consumers do not need a separately installed GLPK.
- Concurrent calls to the library from multiple threads sharing the same solver context are not required to be safe; the limitation will be documented rather than resolved in this feature.
- Standard autotools (`libtool`) will be used to manage shared library versioning and build portability, consistent with the existing build infrastructure.
- The install convention will follow GNU standards: binaries to `$(bindir)`, libraries to `$(libdir)`, the header to `$(includedir)`, and the pkg-config file to `$(libdir)/pkgconfig`.
