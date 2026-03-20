# Research: GitHub Actions CI — Multi-Platform Build & Test

**Feature**: 003-github-actions-ci
**Date**: 2026-03-19

---

## R1: Default Branch Is `master`, Not `main`

**Discovery**: `git remote show origin` confirms `HEAD branch: master`. The GitHub
repository at `github.com/alotau/dwsolver` uses the old Git default `master` branch
name.

**Decision**: All workflow trigger references and branch protection documentation MUST
use `master`. The `spec.md` was updated in place to reflect this.

**Impact**: `ci.yml` must have `push: branches: [master]` and the branch protection
rule must be applied to `master`, not `main`.

---

## R2: Pre-Generated `configure` Script (No Bootstrap Step Needed)

**Discovery**: Both `configure` and `configure.ac` are committed to the repository.
The `configure` script is already generated and ready to run.

**Decision**: CI jobs run `./configure [flags] && make` directly. No `autoreconf`
or `libtoolize` step is necessary. Adding one would be harmless but wasteful.

**Alternatives considered**: Running `autoreconf -i` before configure (rejected: unnecessary,
would also require `autotools` to be installed on all platforms). 

---

## R3: Existing `build.yaml` Partial Implementation

**Discovery**: `.github/workflows/build.yaml` already exists with:
- Single platform: `macos-13`
- Trigger: `push` to `master` only (no PR trigger)
- Build steps: `./configure && make && make install && cd tests && ./dw-tests.sh`
- Missing `--enable-named-semaphores` (broken on macOS)
- Uses `actions/checkout@v3` (outdated)
- No Linux, no Windows, no per-job timeout

**Decision**: Create a new `ci.yml` that supersedes `build.yaml`. Delete `build.yaml`
as part of this implementation to avoid having two conflicting CI definitions. The new
`ci.yml` is strictly a superset of `build.yaml`'s intent.

**Rationale**: Keeping both would create confusion about which is authoritative, and
`build.yaml` is already incorrect (missing `--enable-named-semaphores` on macOS).

---

## R4: Test Runner Must Be Invoked via `bash`, Not `sh`

**Discovery**: `tests/dw-tests.sh` has shebang `#!/bin/sh` but uses `pushd`/`popd`.
The script wraps them as shell functions:
```sh
pushd () { command pushd "$@" > /dev/null; }
popd ()  { command popd  "$@" > /dev/null; }
```
`command pushd` calls the built-in `pushd` directly. In bash, `pushd` is a built-in, so
this works. On Ubuntu, `/bin/sh` is dash; dash has no `pushd` built-in, so `command pushd`
would look for an external `pushd` binary (which doesn't exist) and fail.

**Decision**: Invoke as `bash tests/dw-tests.sh` in all CI jobs. This also works on
macOS (bash pre-installed) and Windows MSYS2 (bash provided). Using `bash` explicitly
is safer and more portable than relying on the shebang.

**Alternatives considered**: Fixing the shebang to `#!/bin/bash` (rejected: modifies
the test script outside this spec's scope; spec 002 owns the test script). Wrapping
with `env bash` (equivalent to just calling `bash tests/dw-tests.sh`).

---

## R5: PATH Setup — Prefer Injecting `src/` Over `make install`

**Discovery**: The existing `build.yaml` uses `make install` to put `dwsolver` on PATH
(default prefix `/usr/local`). On Ubuntu GitHub Actions runners, writing to `/usr/local/bin`
requires `sudo` — the existing workflow would fail on Ubuntu without `sudo make install`.
On macOS GitHub Actions runners, the user has admin rights and `make install` works without
`sudo`.

**Decision**: After `make`, inject `src/` into `$GITHUB_PATH`:
```yaml
- name: Add dwsolver to PATH
  run: echo "$GITHUB_WORKSPACE/src" >> $GITHUB_PATH
```
This is idiomatic for GitHub Actions, requires no `sudo`, is cross-platform, and
avoids modifying the system install directories.

**Alternatives considered**: `sudo make install` on Linux only (rejected: platform-specific
step divergence). `make install` with custom `--prefix=$HOME/.local` (rejected: more
configuration complexity than necessary).

**Windows note**: `$GITHUB_WORKSPACE` is available in MSYS2 shells as `/d/a/repo/repo`
or via `$GITHUB_WORKSPACE` env var. The equivalent in MSYS2: use `$GITHUB_WORKSPACE/src`
with `shell: msys2 {0}` — MSYS2 environment exposes the Windows-style path correctly.

---

## R6: No Pre-Install Step Needed for macOS or Linux

**Discovery**: GitHub-hosted `macos-latest` runners have Xcode CLT (provides `cc`,
`make`, `ld`) pre-installed. GitHub-hosted `ubuntu-latest` runners have `gcc`, `g++`,
and `make` pre-installed. Since `configure` is already generated, neither platform
requires any `apt-get install` or `brew install` step before `./configure && make`.

**Decision**: macOS and Linux CI jobs have zero setup steps beyond `actions/checkout`.
This keeps CI extremely lean and fast.

**Validated by**: The existing `build.yaml` runs on `macos-13` with no install step and
it completes (albeit with the wrong configure flags). 

---

## R7: Windows — MSYS2 MinGW-w64 Toolchain

**Discovery**: Windows GitHub-hosted runners do NOT have a C toolchain pre-installed
that supports POSIX-compatible configure scripts. MSYS2 is the recommended approach for
building Autotools-based projects on Windows in CI.

**Decision**:
- Use `msys2/setup-msys2@v2` action with `msystem: MINGW64`
- Install packages: `mingw-w64-x86_64-gcc` and `make`
  - `autoconf`/`automake` are NOT needed since `configure` is pre-generated
- Enable built-in caching: `cache: true` (satisfies FR-009)
- Run all build+test steps with `shell: msys2 {0}`
- Use `continue-on-error: true` on the job (FR-005) — Windows is best-effort until spec 004

**pthreads on MinGW-w64**: The `mingw-w64-x86_64-gcc` toolchain bundles `winpthreads`
(`libwinpthread-1.dll`). When `configure` detects pthreads, it should link against
`-lpthread` which resolves to winpthreads. However, the project's `configure.ac` may
not correctly probe for Windows pthreads — this is why the job is non-blocking.

**Known risk**: The `dw.h` multiple-definition bug (KD-001, fixed in spec 002) must be
merged before Windows CI can be expected to pass. Even with KD-001 fixed, there may be
additional Windows-specific issues that spec 004 will address.

---

## R8: GitHub Actions Versions to Pin

| Action | Version | Notes |
|--------|---------|-------|
| `actions/checkout` | `@v4` | Current stable; v3 used in build.yaml is outdated |
| `msys2/setup-msys2` | `@v2` | Current stable; built-in caching support |

No other third-party actions required.

---

## R9: Workflow Trigger — `push` + `pull_request`

**Decision**:
```yaml
on:
  push:
    branches: [master]
  pull_request:
    branches: [master]
```
`pull_request` targeting `master` covers all feature branches opening PRs.
`push` to `master` covers direct merges and catches any post-merge regressions.

**Note**: GitHub Actions runs `pull_request` workflows in the context of the merge
commit (target + PR branch), so test failures reflect what would land in `master`
if the PR were merged — exactly what we want.

---

## R10: `four_sea` Test Timing

**Context**: `four_sea` (per spec 002 T007) launches 4 subproblem threads and is the
most CPU-intensive test. GitHub-hosted runners have 2-4 vCPUs. The test has historically
completed in under 60 seconds on a MacBook Pro.

**Decision**: Per-job `timeout-minutes: 30` (FR-006) is generous enough to handle slow
runners and the full 6-test suite (bertsimas, lasdon, mitchell, trick, four_sea, book_dantzig).

**Dependency**: The `four_sea` and `book_dantzig` tests are added by spec 002 (T007/T008).
If spec 002 is NOT merged before this CI spec, the test suite currently has 4 tests and
`four_sea`/`book_dantzig` are not present. CI will still pass with 4 tests.

---

## Summary: All NEEDS CLARIFICATION Items Resolved

| ID | Question | Answer |
|----|----------|--------|
| NC-1 | Default branch name | `master` |
| NC-2 | Is `configure` pre-generated? | Yes — committed to repo |
| NC-3 | Existing CI state | `build.yaml` is partial; supersede with `ci.yml` |
| NC-4 | `pushd`/`popd` compatibility | Must invoke via `bash tests/dw-tests.sh` |
| NC-5 | PATH setup for `dwsolver` | Inject `src/` into `$GITHUB_PATH` |
| NC-6 | macOS/Linux pre-install steps needed? | No — toolchains pre-installed on runners |
| NC-7 | Windows approach | MSYS2 + MinGW-w64 via `msys2/setup-msys2@v2` |
| NC-8 | Actions versions | `checkout@v4`, `msys2/setup-msys2@v2` |
