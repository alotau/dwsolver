# Quick Start: Containerized Build for dwsolver

**Feature**: 004 — Containerized Build
**Date**: 2026-03-19

This guide covers everything an implementer needs to create, validate, and
maintain the Docker support for dwsolver.

---

## Prerequisites

- Docker Engine ≥ 20.10 installed locally
- Repository cloned: `git clone https://github.com/alotau/dwsolver.git`
- Branch checked out: `git checkout 004-docker-containerized-build`

No other tools are required. The build toolchain (gcc, automake, etc.) runs
inside the Docker builder stage.

---

## Step 1 — Create `.dockerignore`

Create `.dockerignore` at the repository root with the following content:

```dockerignore
# Version control
.git/

# Compiled artefacts (stale local builds)
*.o
*.a
*.Po
src/dwsolver

# Configure output
config.log
config.status
config.h
stamp-h1
libtool
autom4te.cache/

# Solver run outputs (test artefacts)
**/out_terminal
**/relaxed_solution
**/rs_sorted
**/done.cpxlp
**/pre_master.cpxlp

# Misc
tmp/
.vscode/
.specify/
```

**Verify**: After creating the file, check the build context size:
```sh
# Shows what would be sent (dry-run: no build)
docker build --no-cache --dry-run . 2>/dev/null || true
# Or watch context upload size when running a real build
```

---

## Step 2 — Create `Dockerfile`

Create `Dockerfile` at the repository root:

```dockerfile
# syntax=docker/dockerfile:1

# ── Stage 1: Build ────────────────────────────────────────────────────────────
FROM ubuntu:24.04 AS builder

# Install build toolchain
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        automake \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY . .

# Touch pre-generated autotools files so make does not try to re-run
# aclocal/automake (which may be a different version than what generated them).
# This is the same technique used in the GitHub Actions CI workflow (ci.yml).
RUN touch aclocal.m4 configure config.h.in ltmain.sh Makefile.in src/Makefile.in

RUN ./configure && make

# ── Stage 2: Runtime ──────────────────────────────────────────────────────────
FROM ubuntu:24.04 AS runner

# Copy only the compiled binary from the builder stage.
# The binary links against glibc (libc.so.6) and libpthread.so.0, both
# already present in the ubuntu:24.04 base — no extra packages needed.
COPY --from=builder /build/src/dwsolver /usr/local/bin/dwsolver

ENTRYPOINT ["/usr/local/bin/dwsolver"]
```

**Key points**:
- Two stages: `builder` with full toolchain, `runner` with binary only
- `touch` replicates the CI workaround for the autotools version mismatch
- `rm -rf /var/lib/apt/lists/*` keeps the builder layer smaller
- `ENTRYPOINT` uses exec-form (JSON array) — signals go directly to the process
- No `CMD` — all arguments are forwarded from `docker run`

---

## Step 3 — Build and Verify Locally

```sh
# Build the image
docker build -t dwsolver .

# Check the image size (must be ≤ 200 MB)
docker image inspect dwsolver --format '{{.Size}}' | \
    awk '{printf "Image size: %.0f MB\n", $1/1024/1024}'

# Smoke test: usage message appears, command exits non-zero (that is correct)
docker run --rm dwsolver
# Expected: Usage: dwsolver -g <guide_file_name> [options]

# Solve the bundled Bertsimas example
docker run --rm \
    -v "$(pwd)/examples/book_bertsimas:/data" \
    -w /data \
    dwsolver \
    --no-write-final-master --quiet -g guidefile
# Expected: solver runs and exits 0; no output (--quiet)
```

---

## Step 4 — Add the `docker-build` CI Job

Append the following job to `.github/workflows/ci.yml`, after the `build-windows`
job:

```yaml
  docker-build:
    runs-on: ubuntu-latest
    timeout-minutes: 30
    continue-on-error: true
    steps:
      - uses: actions/checkout@v4

      - name: Build Docker image
        run: docker build -t dwsolver .

      - name: Smoke test
        run: docker run --rm dwsolver 2>&1 | grep -q "Usage:"
```

**Notes**:
- `continue-on-error: true` — the job is informational; it MUST NOT block merges
- The smoke test passes if the usage string appears; `grep -q` exits 0 on match
- No volume mount needed for the smoke test — no file I/O required

---

## Step 5 — Commit and Open a PR

```sh
git add Dockerfile .dockerignore .github/workflows/ci.yml
git commit -m "feat: add Docker containerized build (Spec 004)

- Add Dockerfile with two-stage build (builder + runtime)
- Add .dockerignore to exclude .git, object files, artifacts
- Add docker-build CI job (continue-on-error) to ci.yml
- Final image ~78 MB; binary at /usr/local/bin/dwsolver
- Touch trick matches existing ci.yml approach (no autoreconf)"

git push origin 004-docker-containerized-build
# Open PR at: https://github.com/alotau/dwsolver/pull/new/004-docker-containerized-build
```

---

## Failure Mode Reference

| Symptom | Cause | Fix |
|---------|-------|-----|
| `autoreconf: command not found` | `autoreconf` is being called somewhere | Do NOT run `autoreconf` in the Dockerfile; use the `touch` technique only |
| `automake: error: version mismatch` | `make` is trying to re-run automake | The `touch` step did not execute; check WORKDIR is set correctly before `touch` |
| `./configure: No such file or directory` | `.dockerignore` is excluding `configure` | Verify `.dockerignore` does not contain `configure` |
| `symbol(s) not found` on link | Stale object files in build context conflicting | Add `*.o` to `.dockerignore` |
| Image size > 200 MB | Multi-stage build not working; builder stage shipped | Confirm `FROM ubuntu:24.04 AS runner` and `COPY --from=builder` are present |
| `ENTRYPOINT` not forwarding args | CMD is overriding ENTRYPOINT | Remove any `CMD` directive; use `ENTRYPOINT` only |
| CI `docker-build` job blocking merges | `continue-on-error` missing | Add `continue-on-error: true` to the job |

---

## Using the Container (End-User Reference)

### Solve your own problem

```sh
# Place your guide file and all referenced data files in one directory
ls /path/to/my/problem/
# guidefile  master.cplex  sub1.cplex  ...

# Run the solver
docker run --rm \
    -v "/path/to/my/problem:/data" \
    dwsolver -g /data/guidefile
```

### Run with verbose output

```sh
docker run --rm \
    -v "$(pwd)/examples/book_lasdon:/data" \
    dwsolver -g /data/guidefile
```

### Print all available options

```sh
# Run with no args — usage is printed to stdout/stderr; exits non-zero
docker run --rm dwsolver
```

---

## Maintenance

| Event | Action |
|-------|--------|
| Ubuntu 24.04 reaches EOL (April 2029) | Update `FROM ubuntu:24.04` to next LTS in both stages; rebuild and retest |
| New autotools files generated (e.g., automake version bump) | No Dockerfile change needed; the touch list covers the six standard files |
| dwsolver gains a new runtime dependency (e.g., a shared library) | Add `apt-get install` to the runner stage; document in this file |
| Docker BuildKit is deprecated | Update `# syntax=docker/dockerfile:1` comment; no functional impact expected |
