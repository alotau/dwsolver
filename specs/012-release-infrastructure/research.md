# Research: Release Infrastructure

**Feature**: 012-release-infrastructure  
**Date**: 2026-03-22  
**Status**: Complete — all NEEDS CLARIFICATION resolved

---

## Decision 1: How to pass GLPK build variables into `make distcheck`

**Question**: Ubuntu's `libglpk-dev` does not ship a `glpk.pc` pkg-config file.
`configure` detects GLPK via `PKG_CHECK_MODULES`, which fails unless
`GLPK_CFLAGS` and `GLPK_LIBS` are pre-supplied.  How should these be passed
when `make distcheck` unpacks into a clean tree and re-runs `./configure`?

**Decision**: Export `GLPK_CFLAGS` and `GLPK_LIBS` as shell environment
variables before invoking `make distcheck`.  Do NOT pass them inside
`DISTCHECK_CONFIGURE_FLAGS`.

**Rationale**: `DISTCHECK_CONFIGURE_FLAGS` is space-split by make.  A value
like `GLPK_LIBS="-L/path -lglpk"` (two tokens after the `=`) causes the
second token (`-lglpk`) to be passed as a bare `./configure` argument, which
`./configure` rejects with "unrecognized option".  Environment variables are
inherited by every subprocess including the re-invoked `./configure` inside the
distcheck tree, so they "just work" with no quoting games.

**Alternatives considered**:
- Pass `GLPK_CFLAGS=-I/usr/include GLPK_LIBS=-lglpk` inside
  `DISTCHECK_CONFIGURE_FLAGS` — works only when `GLPK_LIBS` has no spaces
  (true on Ubuntu where `-lglpk` is a single token), but is fragile and
  inconsistent with the existing CI workflows.
- Add a `pkg-config` wrapper that generates `glpk.pc` — unnecessary complexity;
  the env-var approach is simpler.

**Verified**: On macOS via three distcheck runs; the env-var approach passes
cleanly on both the outer `./configure` and the inner distcheck `./configure`.

---

## Decision 2: macOS `sem_init` not implemented — distcheck workaround

**Question**: Running `make distcheck` on macOS fails at `make check` inside
the clean tree with `sem_init: Function not implemented`.  macOS does not
implement unnamed POSIX semaphores.

**Decision**: On macOS, add `--enable-named-semaphores --enable-recursive-mutex`
to `DISTCHECK_CONFIGURE_FLAGS`.  The release workflow does not run on macOS
(Linux only), so this workaround applies only to local developer use and is
documented in the README.

**Rationale**: The project has configure-time options specifically to work
around macOS semaphore limitations.  Using them is the intended path.

**Alternatives considered**:
- Skip `make distcheck` on macOS — defeats the purpose of local verification.
- Change the release workflow to also run on macOS — the spec explicitly limits
  the release gate to Linux.

---

## Decision 3: Files missing from `EXTRA_DIST` / `DISTCLEANFILES`

**Question**: Which files were missing from `EXTRA_DIST` (causing distcheck to
fail with file-not-found in the clean tree) and which generated files were not
cleaned up by `make distclean` (causing `distcleancheck` to fail)?

**Findings**:

*Missing from `EXTRA_DIST` (top-level `Makefile.am`)*:
- `examples/` — example LP inputs used by `dw-tests.sh`
- `architecture/` — architecture diagrams
- `specs/` — specification documents
- `Dockerfile` — container build file
- `.github/` — CI workflow definitions
- `third-party/glpk/` patch and docs — needed to document what was replaced

*Missing from `EXTRA_DIST` (`tests/Makefile.am`)*:
- `dw-tests.sh` — main integration test driver
- `test_guidefile.sh` — helper script
- `fixtures/` — test fixture files

*Missing from `DISTCLEANFILES` (`tests/Makefile.am`)*:
- `relaxed_solution` — written by `dw-tests.sh` during test execution
- `rs_sorted` — written by `dw-tests.sh` during test execution
- `out_obj.txt` — written by `dw-tests.sh` during test execution
- `src/dw_blas.Po` — GCC dependency tracking file in `tests/src/`; not inside
  a standard automake `_SOURCES` directory so automake's own distclean does not
  remove it

**Decision**: Add all of the above to the appropriate `EXTRA_DIST` and
`DISTCLEANFILES` variables.

**Verified**: After adding `DISTCLEANFILES`, `make distcheck` passes end-to-end
on macOS and produces `dwsolver-1.2.1.tar.gz`.

---

## Decision 4: GitHub Actions supply-chain pinning

**Question**: Should GitHub Actions be referenced by floating version tags
(e.g., `@v4`, `@v2`) or pinned to specific commit SHAs?

**Decision**: Pin all actions — both first-party (`actions/checkout`) and
third-party (`softprops/action-gh-release`) — to their full commit SHA.
Annotate each pinned reference with a comment showing the corresponding version
tag.

**Rationale**: Floating tags are mutable; a compromised or accidentally
force-pushed tag would silently alter the workflow without a code review.
SHA-pinning is required by OSSF Scorecard and GitHub's own Actions hardening
guide.  The cost is a manual update step when upgrading the action; the benefit
is a tamper-evident, auditable workflow.

**Resolved SHAs** (as of 2026-03-22):
- `actions/checkout@v4` → `34e114876b0b11c390a56381ad16ebd13914f8d5`
- `softprops/action-gh-release@v2` → `153bb8e04406b158c6c84fc1615b65b24149a1fe`

---

## Decision 5: Duplicate tag push behaviour

**Question**: If a maintainer pushes the same `v*` tag twice (e.g., to fix a
botched release), what should the release workflow do?

**Decision**: The workflow fails — `softprops/action-gh-release` will error
because the GitHub Release already exists.  No asset is overwritten.  The
maintainer must manually delete the GitHub Release and the tag, then re-push
the corrected tag.

**Rationale**: Silently overwriting a published release is dangerous (consumers
may have already downloaded and cached the previous tarball).  Failing loudly
forces a deliberate, auditable remediation step.

---

## Decision 6: GitHub Release draft vs. publish status

**Question**: Should the workflow create the release as a draft (requiring
manual promotion) or publish it immediately?

**Decision**: Publish immediately (`draft: false`).

**Rationale**: `make distcheck` already gates the upload — the tarball is only
attached after the full build+test cycle passes inside a clean tree.  A
separate "promote from draft" step adds friction with no corresponding safety
benefit.

---

## All NEEDS CLARIFICATION resolved

No open research items remain.  All six decisions above are captured in the
spec's Clarifications section and reflected in FR-005, FR-007a, FR-010, and
FR-011.
