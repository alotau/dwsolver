# Contract: Release Workflow

**File**: `.github/workflows/release.yml`  
**Feature**: 012-release-infrastructure

Documents the inputs, outputs, and behavioural guarantees of the automated
release workflow.

---

## Trigger Contract

| Property | Value |
|----------|-------|
| Event | `push` |
| Filter | `tags: ['v*']` |
| NOT triggered by | branch pushes, pull requests, manual `workflow_dispatch` |

**Guarantee**: The workflow fires exactly once per `v*` tag push.  Re-pushing
the same tag does not re-trigger a separate run (GitHub deduplicates by ref
SHA only if the SHA changes; a force-push will trigger a new run but the
release-creation step will fail because the release already exists).

---

## Environment Contract

### Required Secrets / Permissions

| Resource | Value | Purpose |
|----------|-------|---------|
| `GITHUB_TOKEN` | Automatically provided | Create GitHub Release and upload assets |
| `permissions.contents` | `write` | Required by `softprops/action-gh-release` |

No additional repository secrets are needed.

### Required System Dependencies (installed in workflow)

| Package | Installed via | Purpose |
|---------|--------------|---------|
| `build-essential` | `apt-get` | GCC compilers for the build |
| `automake` | `apt-get` | Regenerate `Makefile` from `Makefile.in` |
| `pkg-config` | `apt-get` | GLPK detection (partially — see env workaround) |
| `libglpk-dev` | `apt-get` | GLPK headers and shared library |

### Required Environment Variables (set via `$GITHUB_ENV`)

| Variable | Value | Reason |
|----------|-------|--------|
| `GLPK_CFLAGS` | `-I/usr/include` | `libglpk-dev` on Ubuntu ships no `glpk.pc`; bypass pkg-config |
| `GLPK_LIBS` | `-lglpk` | Same reason; inherited by `make distcheck`'s inner `./configure` |

---

## Step Sequence and Exit Conditions

```
checkout
    ↓
install dependencies + set GLPK env vars
    ↓
touch autotools timestamps (suppress regeneration)
    ↓
make distcheck                      ← EXIT non-zero → workflow fails; no release created
    ↓
locate dwsolver-*.tar.gz            ← EXIT non-zero if glob matches nothing
    ↓
softprops/action-gh-release         ← EXIT non-zero if release already exists
    ↓
workflow success
```

**Guarantee**: No GitHub Release is created unless `make distcheck` exits 0.

---

## Output Contract

On success the workflow produces:

| Output | Description |
|--------|-------------|
| GitHub Release | Public release associated with the pushed `v*` tag |
| Release name | Same as tag name (e.g., `v1.2.1`) |
| Release body | Auto-generated from commits since last tag |
| Draft status | `false` — immediately public |
| Asset | `dwsolver-<ver>.tar.gz` attached as a downloadable file |

---

## Failure Modes

| Failure | Cause | Recovery |
|---------|-------|---------|
| `make distcheck` fails | Missing file in `EXTRA_DIST`, test failure, or build error | Fix the issue, delete the tag, push a corrected tag |
| Release already exists | Same tag pushed twice | Delete the GitHub Release and the tag manually, then re-push |
| `GITHUB_TOKEN` lacks `contents: write` | Repository permission settings | Enable "Allow GitHub Actions to create and approve pull requests" and release permissions in repository Settings |

---

## Action Version Pins

| Action | SHA | Tag |
|--------|-----|-----|
| `actions/checkout` | `34e114876b0b11c390a56381ad16ebd13914f8d5` | `v4` |
| `softprops/action-gh-release` | `153bb8e04406b158c6c84fc1615b65b24149a1fe` | `v2` |

SHAs were verified against the GitHub API on 2026-03-22.  Update by running:
```sh
gh api repos/actions/checkout/git/ref/tags/v4 --jq '.object.sha'
gh api repos/softprops/action-gh-release/git/ref/tags/v2 --jq '.object.sha'
```
