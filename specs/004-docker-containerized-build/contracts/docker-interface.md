# Docker Interface Contract: dwsolver Container

**Feature**: 004 — Containerized Build
**Date**: 2026-03-19
**Status**: Authoritative

This document defines the stable interface contract for the `dwsolver` Docker
image. All implementation decisions MUST honour these semantics. Any change to
the contract requires a corresponding spec update.

---

## Build Contract

### Command

```sh
docker build -t <image-tag> <build-context>
```

**Canonical invocation** (from repository root):
```sh
docker build -t dwsolver .
```

### Prerequisites

| Requirement | Detail |
|-------------|--------|
| Docker Engine | ≥ 20.10 (BuildKit enabled by default) |
| Build context | Repository root (`Dockerfile` + `.dockerignore`) |
| Internet access | Required during build to pull base image (`ubuntu:24.04`) and run `apt-get` |

### Build Stages

| Stage name | Base | Purpose |
|-----------|------|---------|
| `builder` | `ubuntu:24.04` | Installs toolchain, builds binary |
| `runner`  | `ubuntu:24.04` | Copies binary only; no toolchain |

### Expected Outcomes

- Image builds without errors
- Final image tag is `dwsolver` (or the tag supplied by `-t`)
- Final image size ≤ 200 MB
- Final image contains `/usr/local/bin/dwsolver` on `PATH`
- Final image contains NO `.c`/`.h` source files, `.o` object files, or build tools

---

## Run Contract

### Command

```sh
docker run [--rm] [-v <host-dir>:<container-dir>] <image-tag> [dwsolver-args]
```

### ENTRYPOINT

```dockerfile
ENTRYPOINT ["/usr/local/bin/dwsolver"]
```

- All arguments after the image name are forwarded directly to `dwsolver`
- No `CMD` default is defined; running the image with no arguments prints
  dwsolver's usage message and exits non-zero (identical to native behaviour)
- Shell wrapping is NOT used; signals are delivered directly to the process

### Volume Mounts

```sh
-v <host-absolute-path>:/data
```

| Item | Detail |
|----|---|
| Convention | Mount the directory containing the guide file and all data files to `/data` inside the container |
| Access mode | Read-write (dwsolver writes output files — e.g., `relaxed_solution`, `out_terminal` — to the same directory as the guide file) |
| Multiple mounts | Supported; use any container path, then reference it in the `-g` argument |
| No-mount usage | Valid for smoke tests (`docker run --rm dwsolver` → usage message) |

### Argument Forwarding

All standard dwsolver flags are passed through unchanged:

| Flag | Description |
|------|-------------|
| `-g <path>` | Guide file path (required for a solve run) |
| `--quiet` | Suppress verbose output |
| `--no-write-final-master` | Do not write the final master LP file |
| `--help` / no args | Print usage and exit |

### Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Solver completed successfully |
| Non-zero | Solver error, bad arguments, or infeasible problem |

Exit codes are inherited directly from the `dwsolver` binary. No container-level
exit code mapping is performed.

---

## Canonical Usage Examples

### Solve a problem with local data

```sh
# Build once
docker build -t dwsolver .

# Run against local data directory
docker run --rm \
  -v "$(pwd)/examples/book_bertsimas:/data" \
  dwsolver \
  --no-write-final-master --quiet -g /data/guidefile
```

### Print usage

```sh
docker run --rm dwsolver
# Prints: Usage: dwsolver -g <guide_file_name> [options]
# Exits non-zero
```

### Run against an arbitrary directory

```sh
docker run --rm \
  -v "/absolute/path/to/my/problem:/data" \
  dwsolver -g /data/my_guidefile
```

---

## CI Contract

### Job name

`docker-build`

### Location in `ci.yml`

Added as a fourth job alongside `build-macos`, `build-linux`, `build-windows`.

### Blocking behaviour

```yaml
continue-on-error: true
```

The `docker-build` job MUST NOT be a required status check. It MUST NOT block
merges. Its sole purpose is visibility — catching Dockerfile regressions early.

### Smoke test command

```sh
docker run --rm dwsolver 2>&1 | grep -q "Usage:"
```

This verifies:
1. The image was built successfully (prior `docker build` step)
2. The ENTRYPOINT is correctly configured
3. The binary executes and produces expected output

### Trigger

Same as all other CI jobs: push to `master` and pull requests targeting `master`.

---

## Invariants

The following invariants MUST hold after every merge to `master`:

1. `docker build -t dwsolver .` completes without error from the repository root
2. `docker run --rm dwsolver 2>&1 | grep -q "Usage:"` exits with code 0
3. The final image layers contain no `.c`/`.h` source files (verified by multi-stage build)
4. The final image size is ≤ 200 MB
