# Research: Containerized Build for dwsolver

**Phase**: 0 — Unknowns & Best Practices
**Branch**: `004-docker-containerized-build`
**Date**: 2026-03-19

All unknowns from the Technical Context have been resolved. This document
captures each key decision, the selected approach, its rationale, and the
alternatives that were evaluated and rejected.

---

## R1 — Base Image Selection

**Decision**: `ubuntu:24.04` (Noble Numbat, LTS)

**Rationale**: This is the same platform the `build-linux` CI job runs on
(`ubuntu-latest` → 24.04 as of 2026). Pinning to an exact version ensures
that the Dockerfile is reproducible and does not silently break when the
upstream `ubuntu:latest` tag rolls forward. Ubuntu 24.04 LTS is supported
until April 2029, giving ample runway before the next refresh.

**Alternatives Considered**:
- `ubuntu:latest` — rejected: unpinned; rolls forward silently, breaks reproducibility
- `ubuntu:22.04` — viable but older; 24.04 is current LTS and matches the CI environment
- `debian:bookworm-slim` — smaller base, but diverges from the known-working CI platform
- `alpine` — uses musl libc; the pre-compiled binary links against glibc (`libc.so.6`);
  running a glibc binary on musl requires `gcompat` shims and is fragile

---

## R2 — Build Technique: Touch Trick vs. Autoreconf

**Decision**: `touch aclocal.m4 configure config.h.in ltmain.sh Makefile.in src/Makefile.in`
before `./configure && make` — identical to what `ci.yml` already does.

**Rationale**: The pre-generated autotools files committed to the repository were
produced with `automake-1.17` on macOS. When `make` runs, it compares the timestamps
of `Makefile.am` (the input) against `Makefile.in` (the generated output). If the
input appears newer, `make` invokes `automake` to regenerate. `ubuntu:24.04`'s
`automake` package provides version 1.16.5 — not 1.17 — which produces different
(incompatible) output and causes the build to fail. `touch`-ing the pre-generated files
advances their timestamps so `make` treats them as up-to-date and skips regeneration.
This technique is already proven in Spec 003 (`ci.yml`) on every CI run.

**Alternatives Considered**:
- `autoreconf -fi` inside the container — rejected: requires automake-1.17 (not
  available in Ubuntu package repos), and was the exact failure mode of the
  abandoned `set-up-docker` WIP branch
- `apt-get install automake-1.17` — no such package exists in Ubuntu 24.04 repos
- Pinning `automake=1.17` via backport — no backport available; would require
  source builds of automake, significantly complicating the Dockerfile

---

## R3 — Multi-Stage Build

**Decision**: Two-stage Dockerfile: a `builder` stage with full toolchain, and a
`runner` stage with only the binary.

**Rationale**: The builder stage needs `build-essential` (~200 MB of gcc, make, etc.)
and `automake` on top of the base Ubuntu image, plus all source files. The runtime
stage needs none of that — glibc (`libc.so.6`) and `libpthread.so.0` are already
present in the Ubuntu base image. A single-stage build would produce an image of
600+ MB. The multi-stage approach yields a final image of approximately 78 MB (77 MB
Ubuntu base + 1 MB binary), comfortably within the 200 MB requirement (SC-002).

**Verified runtime dependencies**: Running `otool -L src/dwsolver` on macOS confirms
the binary links only against the system library. On Linux the equivalent is confirmed
by the `build-linux` CI job succeeding on a stock `ubuntu-latest` runner with no extra
runtime package installs. The runtime stage requires no `apt-get install` commands.

**Alternatives Considered**:
- Single-stage build — simple but violates SC-002 (<200 MB) and User Story 3
- Distroless base for runtime — minimal, but requires additional investigation of
  glibc variant compatibility; not worth the complexity for a small CLI binary

---

## R4 — ENTRYPOINT and Argument Forwarding

**Decision**: `ENTRYPOINT ["/usr/local/bin/dwsolver"]` with no `CMD`.

**Rationale**: With `ENTRYPOINT` set and no `CMD`, Docker appends any arguments the
user passes on `docker run` directly to the entrypoint binary. This means
`docker run --rm dwsolver -g /data/guidefile` works intuitively — `-g /data/guidefile`
is passed to `dwsolver` as-is. The exec-form (`[]`) is used to avoid shell wrapping,
which ensures signals are delivered directly to the process (Linux container best
practice).

**Alternatives Considered**:
- `CMD ["dwsolver"]` without `ENTRYPOINT` — using CMD means the entire command is
  replaced when the user passes arguments; `docker run img -g foo` would run `-g foo`
  as a command, not pass it to dwsolver
- `ENTRYPOINT ["/bin/sh", "-c", "dwsolver \"$@\""]` — shell form; unnecessary complexity
- `ENTRYPOINT ["dwsolver"]` with `CMD ["--help"]` — adds an unintuitive default;
  complicates the CI smoke test since it requires an override

---

## R5 — Binary Install Location

**Decision**: Copy `src/dwsolver` from the builder stage directly into
`/usr/local/bin/dwsolver` in the runtime stage, without running `make install`.

**Rationale**: `make install` installs to `$(prefix)` (default `/usr/local`), which
means invoking `make install DESTDIR=/staging` and copying would work — but is
unnecessary complexity. The build output is a single binary at `src/dwsolver`. A
direct `COPY --from=builder /build/src/dwsolver /usr/local/bin/dwsolver` is cleaner.
`/usr/local/bin` is on the default `PATH` in Ubuntu containers.

**Alternatives Considered**:
- `make install` in builder stage, then `COPY --from=builder /usr/local/bin/dwsolver ...`
  — equivalent outcome, but requires the builder's `/usr/local/bin` to be the install
  target, adding a confusing layer boundary
- Copy to `/usr/bin/dwsolver` — works, but `/usr/local/bin` is the conventional home
  for locally compiled software

---

## R6 — .dockerignore Contents

**Decision**: Exclude `.git`, `*.o`, `*.a`, `src/dwsolver`, `config.log`,
`config.status`, `autom4te.cache/`, `libtool`, `stamp-h1`, `tmp/`, and example
output files (`**/out_terminal`, `**/relaxed_solution`, `**/rs_sorted`,
`**/done.cpxlp`, `**/pre_master.cpxlp`). Include source files, pre-generated
autotools files, and the `examples/` directory (needed for CI smoke tests).

**Rationale**: The `.git` directory is the single largest exclusion — it can be 100+
MB on a repository with history. Object files (`*.o`) and the compiled binary
(`src/dwsolver`) should not be sent as build context; they are irrelevant and could
cause confusion if the builder picks up a stale local binary. Example output files are
test artifacts that would be written during previous local runs; excluding them keeps
the build context clean. The entire list mirrors the project's `.gitignore`.

**Alternatives Considered**:
- Not including a `.dockerignore` — Docker would send the entire working tree including
  `.git` as build context, wasting bandwidth and potentially including stale artifacts
- Excluding `examples/` — the examples are needed for the CI docker-build smoke test
  (volume-mount approach requires the directory to exist on the runner)

---

## R7 — CI Docker Build Job Design

**Decision**: Add a `docker-build` job to `ci.yml` that:
1. Checks out the repository (`actions/checkout@v4`)
2. Runs `docker build -t dwsolver .`
3. Runs a smoke test: `docker run --rm dwsolver 2>&1 | grep -q "Usage:"`
4. Is marked `continue-on-error: true`

**Rationale**: The smoke test verifies that the ENTRYPOINT works and the binary
executes correctly inside the container. Invoking `dwsolver` with no arguments
causes it to print a usage message to stderr and exit non-zero — `2>&1 | grep -q "Usage:"`
confirms the binary is functional without requiring a file solve (which would need
writable volume mounts and reference data comparison). Marking the job
`continue-on-error: true` means it never blocks the required `build-linux` and
`build-macos` checks.

**Alternatives Considered**:
- Full solve smoke test with volume mount: `docker run --rm -v $WORKSPACE/examples/book_bertsimas:/data dwsolver --no-write-final-master --quiet -g /data/guidefile`
  — more thorough but requires write access to the mounted directory for output files;
  complicates teardown and is slower than a usage check
- No CI job — rejected; the entire purpose of the `set-up-docker` branch being broken
  for years is that no CI caught it
- Separate workflow file — unnecessary; the `docker-build` job belongs logically with
  the other build jobs in `ci.yml`
