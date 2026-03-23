# Data Model: Release Infrastructure

**Feature**: 012-release-infrastructure  
**Date**: 2026-03-22

This feature does not introduce persistent data storage.  The "entities" are
build-time and CI-time artifacts with defined structures, relationships, and
validation rules.

---

## Entity: Source Tarball

**Produced by**: `make dist` / `make distcheck`  
**Filename pattern**: `dwsolver-<major>.<minor>.<patch>.tar.gz`  
**Location**: Repository root after a successful `make distcheck`

### Fields

| Field | Source | Example | Validation rule |
|-------|--------|---------|-----------------|
| `version` | `AC_INIT` in `configure.ac` | `1.2.1` | Must match `<major>.<minor>.<patch>` (semver, numeric components only) |
| `filename` | Derived from version | `dwsolver-1.2.1.tar.gz` | Must match `dwsolver-${version}.tar.gz` |
| `contents` | Autotools `make dist` | See tarball contents below | All items in `EXTRA_DIST` across all `Makefile.am` files must be present |

### Required Contents

The tarball must contain (verified by `tar tzf dwsolver-<ver>.tar.gz`):

```
dwsolver-<ver>/
├── src/           — all C source and header files
├── tests/         — dw-tests.sh, test_guidefile.sh, fixtures/, C unit test sources
├── examples/      — LP input files for all reference problems
├── specs/         — feature specification documents
├── architecture/  — architecture diagrams
├── Dockerfile     — container build definition
├── .github/       — CI and release workflow definitions
├── third-party/   — GLPK patch and license documentation
├── README.md
├── COPYING
├── AUTHORS
├── ChangeLog
├── NEWS
├── INSTALL
├── configure.ac
├── Makefile.am
└── dwsolver.pc.in
```

### Validation Rules

- Must be produced only from a clean working tree (no generated files present
  before `make check`).
- `make distcheck` must exit 0 — this implies: the tarball unpacks, configures,
  builds, passes `make check`, installs, uninstalls, and leaves no files behind
  after `make distclean`.
- The version in the filename must match `AC_INIT(dwsolver, VERSION, ...)` in
  `configure.ac`.

---

## Entity: Version-Info Triple

**Location**: `src/Makefile.am`, `libdwsolver_la_LDFLAGS` variable  
**Format**: `CURRENT:REVISION:AGE` (libtool convention)  
**Current value**: `0:0:0`

### Fields

| Field | Type | Current | Meaning |
|-------|------|---------|---------|
| `CURRENT` | non-negative integer | `0` | Number of the current interface (incremented when the public API changes incompatibly or gains new symbols) |
| `REVISION` | non-negative integer | `0` | Implementation revision of the current interface (incremented on bug-fix releases; reset to 0 when CURRENT changes) |
| `AGE` | non-negative integer | `0` | Number of previous interfaces still supported (incremented when adding new symbols; reset to 0 when removing or changing symbols) |

### Derived soname

The installed shared library soname is `libdwsolver.so.(CURRENT - AGE)`.  With
`0:0:0` the soname is `libdwsolver.so.0`.

### Update Rules (from spec FR-008)

| Release type | CURRENT | REVISION | AGE |
|-------------|---------|----------|-----|
| Bug fix only (no API change) | unchanged | +1 | unchanged |
| New public API symbols added (backward compatible) | +1 | reset 0 | +1 |
| Public API symbol removed or changed (breaking) | +1 | reset 0 | reset 0 |

---

## Entity: Release Tag

**Format**: `v<major>.<minor>.<patch>`  
**Example**: `v1.2.1`  
**Location**: Git repository (pushed to `origin`)

### Fields

| Field | Constraint |
|-------|-----------|
| `name` | Must match regex `^v[0-9]+\.[0-9]+\.[0-9]+$` |
| `target` | Must point to a commit on `main` |
| `version` | Numeric portion must match `AC_INIT` version in `configure.ac` at the tagged commit |

### State Transitions

```
[no tag] → push v* tag → [tag exists on origin]
                               ↓
                    release workflow triggered
                               ↓
                    make distcheck passes → [GitHub Release created + tarball attached]
                               ↓ (if re-tag needed)
                    maintainer deletes release + tag → [no tag] → re-push
```

---

## Entity: GitHub Release

**Created by**: `.github/workflows/release.yml` via `softprops/action-gh-release`  
**Trigger**: Successful `make distcheck` on a `v*` tag push

### Fields

| Field | Value | Notes |
|-------|-------|-------|
| `tag_name` | The pushed `v*` tag | Set automatically by the action |
| `name` | Same as `tag_name` (default) | Can be overridden in workflow `with:` |
| `draft` | `false` | Published immediately; no manual promotion |
| `prerelease` | `false` (default) | Can be set manually for `-rc` tags |
| `body` | Auto-generated from commits | `generate_release_notes: true` in workflow |
| `assets` | `dwsolver-<ver>.tar.gz` | Uploaded at release creation time |

### Constraints

- One GitHub Release per tag: pushing the same tag twice will fail the workflow
  (release already exists).
- The tarball asset name must match `dwsolver-*.tar.gz`.
- The release is only created after `make distcheck` exits 0.
