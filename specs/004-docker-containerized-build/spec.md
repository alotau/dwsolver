# Feature Specification: Containerized Build for dwsolver

**Feature Branch**: `004-docker-containerized-build`
**Created**: 2025-07-22
**Status**: Draft
**Supersedes**: `set-up-docker` (abandoned WIP branch — broken Dockerfile, incomplete source changes, never merged)

## Overview

dwsolver requires a C build toolchain (autoconf, automake, libtool, gcc) that is non-trivial to install on macOS, Windows, or unfamiliar Linux distributions. This feature delivers a `Dockerfile` and `.dockerignore` at the repository root so that any user with Docker can build and run dwsolver without installing any native build tools.

The container image must support the primary user workflow: solving a Dantzig–Wolfe decomposition problem by pointing the solver at a guide file and associated data files that live on the user's local filesystem.

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Run Solver Against Local Data Files (Priority: P1)

A researcher or engineer has a set of input files (guide file + data files) on their local machine and wants to solve a DW problem without configuring a build toolchain.

**Why this priority**: This is the only reason to have a container image. Everything else (build, CI check) is in service of this outcome.

**Independent Test**: Can be fully tested by building the image, mounting the `examples/book_bertsimas` directory, and running the solver against the bundled example guide file. The solver must exit successfully and produce the same output as a native run.

**Acceptance Scenarios**:

1. **Given** a user has Docker installed and has cloned the repository, **When** they run `docker build -t dwsolver .` from the repo root, **Then** the image builds without errors and the resulting image is under 200 MB.
2. **Given** a successfully built image and a directory containing a dwsolver guide file and its data files, **When** the user runs `docker run --rm -w /data -v /path/to/data:/data dwsolver -g guidefile`, **Then** the solver runs, reads all data correctly from the mounted path (including any relative paths referenced from the guide file), and exits with code 0 when the problem is feasible.
3. **Given** the same image and the bundled `examples/book_bertsimas` example, **When** the user runs the container against the example guide file, **Then** the output matches the expected result from a native build.

---

### User Story 2 — Build the Image in CI (Priority: P2)

A contributor wants confidence that the Dockerfile is not broken after source code changes, without gating the primary CI checks on Docker availability.

**Why this priority**: CI coverage prevents the Dockerfile from silently breaking (as the `set-up-docker` branch demonstrates), but Docker image builds are slower than the native build jobs and should not block merges.

**Independent Test**: A CI job named `docker-build` can be added that runs `docker build` and validates the resulting image with a smoke test. The job is marked `continue-on-error: true` so it does not block the required check gates.

**Acceptance Scenarios**:

1. **Given** a pull request that modifies `Dockerfile`, `.dockerignore`, or any source file, **When** the CI runs, **Then** a `docker-build` job executes and reports pass or fail without blocking the merge.
2. **Given** a successful Docker build in CI, **When** the job runs a smoke test (`docker run --rm dwsolver --help` or a quick example solve), **Then** the job exits with code 0.
3. **Given** a broken Dockerfile, **When** the CI runs, **Then** the `docker-build` job fails visibly in the PR checks but the required `build-linux` and `build-macos` checks are not affected.

---

### User Story 3 — Exclude Build Artifacts from the Image (Priority: P3)

A user pulling the image should receive only the solver binary and its runtime dependencies — not the full source tree, intermediate object files, or any files that are irrelevant to running the solver.

**Why this priority**: A bloated image containing `.o` files, `Makefile` fragments, or `.git` history increases pull time and storage cost for every user. Keeping final images lean is a baseline best practice.

**Independent Test**: After building, inspect the image layers. The final runnable image MUST NOT contain `.o` files, the `.git` directory, or C source files. A `.dockerignore` file at the repo root enforces this at build time.

**Acceptance Scenarios**:

1. **Given** a `.dockerignore` file at the repo root, **When** `docker build` runs, **Then** the `.git` directory, any `*.o` object files, and any previously compiled `src/dwsolver` binary are excluded from the build context.
2. **Given** the final image, **When** a user inspects its filesystem, **Then** C source files (`.c`, `.h`) are not present in any layer visible to the running container.

---

### Edge Cases

- What happens when the mounted data directory is empty or the guide file path is wrong? The solver must exit with a non-zero code and print a user-readable error message to stderr (existing solver behavior — unchanged by this feature).
- What happens when the host architecture is ARM64 (Apple Silicon)? The image must build and run correctly on both `linux/amd64` and `linux/arm64`.
- What happens when the user builds the image on a machine with a newer or older Docker version? The Dockerfile must use only stable, widely-supported Dockerfile instructions with no experimental syntax.
- What happens if someone runs the container with no arguments? The solver must print its usage help and exit with a non-zero code (existing solver behavior — unchanged by this feature).

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The repository root MUST contain a `Dockerfile` that, when built with `docker build`, produces a runnable image containing the `dwsolver` binary.
- **FR-002**: The `Dockerfile` MUST use a pinned, versioned base image (e.g., `ubuntu:24.04`) so that builds are reproducible and do not silently break when upstream images are updated.
- **FR-003**: The build inside the container MUST use the pre-generated autotools files already committed to the repository. The build MUST `touch` the six pre-generated files (`aclocal.m4`, `configure`, `config.h.in`, `ltmain.sh`, `Makefile.in`, `src/Makefile.in`) before invoking `./configure && make` — consistent with the approach already proven in the GitHub Actions CI workflow (Spec 003).
- **FR-004**: The `Dockerfile` MUST produce an image where the `dwsolver` binary is on the default `PATH` so users can invoke it without specifying a full path inside the container.
- **FR-005**: The repository root MUST contain a `.dockerignore` file that excludes at minimum: the `.git` directory, any compiled binary (`src/dwsolver`), and object files (`*.o`, `*.a`).
- **FR-006**: The final runnable container image MUST be no larger than 200 MB. A multi-stage build (separate build stage and runtime stage, copying only the compiled binary and required shared libraries to the runtime stage) is the preferred mechanism to achieve this.
- **FR-007**: The `ENTRYPOINT` or `CMD` in the `Dockerfile` MUST be structured so that users can pass dwsolver command-line arguments directly after the image name (e.g., `docker run --rm dwsolver --help`).
- **FR-008**: The `Dockerfile` MUST NOT modify any source files (`src/*.c`, `src/*.h`) during the build. All necessary source fixes are already present on the `master` branch (Spec 002).
- **FR-009**: The CI workflow (`ci.yml`) MUST be extended with a `docker-build` job that builds the image and runs a smoke test against the bundled example. This job MUST be marked `continue-on-error: true` so it never blocks the required merge checks.
- **FR-010**: The `set-up-docker` branch is abandoned. No code or configuration from that branch should be ported into this feature. Its broken partial source changes are superseded by Spec 002.

### Key Entities

- **Docker image**: The built artifact. Consists of a build stage (full toolchain + sources) and a runtime stage (minimal OS + binary). Identified by a local tag (e.g., `dwsolver`).
- **Build context**: The set of files sent to the Docker daemon when `docker build` runs. Controlled by `.dockerignore` to exclude irrelevant or large files.
- **Volume mount**: A host directory bound into the container at a known path (e.g., `/data`) so the solver can read guide files and data files from the user's local filesystem.
- **Guide file**: A dwsolver-specific input file that references other data files by relative path. When using the container, all referenced files must exist within the mounted directory.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Users with only Docker installed can build the image from a fresh `git clone` with a single command and begin solving problems within 5 minutes on a standard broadband connection.
- **SC-002**: The final runnable container image is no larger than 200 MB, ensuring fast pull times and low storage cost.
- **SC-003**: The solver produces identical output inside the container versus a native Linux build when run against the bundled `examples/book_bertsimas` example.
- **SC-004**: The CI `docker-build` job detects a broken Dockerfile on every pull request that modifies the Dockerfile or source files, without adding to the required check gates.
- **SC-005**: Users can run the container against their own data using a single `docker run` command with a volume mount — no additional configuration files or environment variables required.

## Assumptions

- Docker Engine (or Docker Desktop) is the only prerequisite on the user's machine. No other tools (make, gcc, autoconf, etc.) are required.
- The bundled example data in `examples/book_bertsimas` remains in the repository and serves as the smoke-test input for both local and CI validation.
- The `master` branch already has all source-level fixes in place (Spec 002: `dw_globals.c`, `extern` declarations, `snprintf` replacements, OOM guards). The Dockerfile does not need to apply any source patches.
- The CI runner for the `docker-build` job has Docker available (GitHub-hosted `ubuntu-latest` runners include Docker by default).
- Named POSIX semaphores (required for the macOS build) are not relevant inside the Linux container; the Linux build path in the autotools configuration handles this automatically.
