# Quickstart: GitHub Actions CI — Implementation & Maintenance

**Feature**: 003-github-actions-ci
**Date**: 2026-03-19
**Audience**: Developer implementing or maintaining the CI workflow

---

## Prerequisites

- Push access to `github.com/alotau/dwsolver`
- Admin access to repository settings (for branch protection rules)
- Working on branch `003-github-actions-ci` (branched from `master`)

---

## Implementation Steps

### Step 1: Create `.github/workflows/ci.yml`

Create the workflow file at `.github/workflows/ci.yml`. The complete workflow
is defined in [contracts/ci-workflow.md](contracts/ci-workflow.md). It contains:

- Triggers: `push` and `pull_request` on `master`
- Three jobs: `build-macos`, `build-linux`, `build-windows`
- macOS: `./configure --enable-named-semaphores && make`
- Linux: `./configure && make`
- Windows: MSYS2 MinGW-w64, `continue-on-error: true`
- All jobs: `cd tests && bash dw-tests.sh`, timeout 30 minutes

### Step 2: Delete `.github/workflows/build.yaml`

Remove the existing partial CI workflow to avoid duplicate runs and confusion.
`build.yaml` is superseded by `ci.yml`.

```bash
git rm .github/workflows/build.yaml
```

### Step 3: Commit and Push

```bash
git add .github/workflows/ci.yml
git commit -m "ci: add multi-platform GitHub Actions workflow"
git push origin 003-github-actions-ci
```

### Step 4: Verify CI Runs on the PR

Open a pull request from `003-github-actions-ci` to `master`. Go to the "Checks"
tab of the PR. All three jobs (`build-macos`, `build-linux`, `build-windows`)
should appear within ~60 seconds of opening the PR.

Expected initial results:
- `build-macos`: PASS (macOS build was verified in the project setup phase)
- `build-linux`: Depends on spec 002 merge status. Pre-fix: FAIL (KD-001). Post-fix: PASS.
- `build-windows`: FAIL acceptable (continue-on-error, non-blocking)

### Step 5: Set Branch Protection Rules

Once the PR is open and CI is running, configure branch protection on `master`:

1. Go to `github.com/alotau/dwsolver` → **Settings** → **Branches**
2. Click **Add rule** (or **Edit** if a rule already exists for `master`)
3. Branch name pattern: `master`
4. Check: **Require status checks to pass before merging**
5. Search for and add required checks:
   - `CI / build-macos`
   - `CI / build-linux`
6. Do NOT add `CI / build-windows` as a required check
7. (Optional) Check **Require branches to be up to date before merging**
8. Save

> **Note**: The status check names `CI / build-macos` and `CI / build-linux` only appear
> in the search box after at least one CI run has completed. Complete step 4 before step 5.

---

## Testing Locally Before Push

To simulate what CI runs:

**macOS:**
```bash
cd /path/to/dwsolver-repaired
./configure --enable-named-semaphores
make
export PATH="$PWD/src:$PATH"
cd tests && bash dw-tests.sh
```

**Linux:**
```bash
cd /path/to/dwsolver-repaired
./configure
make
export PATH="$PWD/src:$PATH"
cd tests && bash dw-tests.sh
```

All tests should report `PASS` and exit 0.

---

## Common Failure Modes

| Symptom | Likely Cause | Fix |
|---------|-------------|-----|
| macOS job fails at configure | Missing `--enable-named-semaphores` | Check `ci.yml` macOS configure step |
| Linux job fails at link step | KD-001 not fixed (spec 002) | Merge spec 002 into this branch first |
| All jobs fail at test step | `dwsolver` not on PATH | Check the `echo ... >> $GITHUB_PATH` step runs before the test step |
| Tests fail with `pushd: not found` | Tests invoked via `sh` not `bash` | Ensure test step uses `bash dw-tests.sh`, not `./dw-tests.sh` |
| Windows job fails | Expected (KD-001 + platform gaps) | Confirm `continue-on-error: true` is set |
| Branch protection check name not found | CI hasn't completed a run yet | Run CI once (push the branch), then search for check names |

---

## Maintenance

**When adding a new test to `dw-tests.sh`**: No CI changes needed. The test script
runs automatically; new tests are included on next run.

**When promoting Windows to required**: After spec 004 confirms Windows builds cleanly:
1. Remove `continue-on-error: true` from `build-windows` in `ci.yml`
2. Add `CI / build-windows` to branch protection required checks

**When changing the workflow `name`**: Update branch protection required check names
to `{new name} / build-macos` and `{new name} / build-linux`.

**To force re-run a failed CI job**: Go to the PR → Checks → click the failed job →
"Re-run failed jobs". No code change needed.
