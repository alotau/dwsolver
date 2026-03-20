# Contract: GitHub Actions Workflow — `ci.yml`

**Feature**: 003-github-actions-ci
**Date**: 2026-03-19
**File**: `.github/workflows/ci.yml`

---

## Purpose

This contract defines the observable interface of the CI workflow: what triggers it,
what status checks it produces, and what success/failure signals downstream systems
(branch protection, PR UI) consume.

---

## Trigger Contract

| Event | Condition | Effect |
|-------|-----------|--------|
| `push` | branch = `master` | All 3 jobs run |
| `pull_request` | base branch = `master` | All 3 jobs run (in merge-commit context) |

No other triggers. Workflow does NOT run on:
- Direct pushes to feature branches (only PRs targeting master)
- Tag pushes
- Manual dispatch (not configured)
- Schedule

---

## Output Contract: Status Checks

The workflow produces the following status check identifiers (as they appear in GitHub
PR checks and branch protection configuration):

| Check Name | Job ID | Required to Merge? |
|------------|--------|--------------------|
| `CI / build-macos` | `build-macos` | Yes — add to branch protection |
| `CI / build-linux` | `build-linux` | Yes — add to branch protection |
| `CI / build-windows` | `build-windows` | No — informational only |

**IMPORTANT**: If the top-level workflow `name` field changes from `CI`, the status
check names change accordingly. Branch protection rules must be updated to match.

---

## Failure Semantics

A job fails (`failure` exit status) if and only if:

1. Any step exits with a non-zero code, OR
2. The job exceeds `timeout-minutes: 30`

For `build-macos` and `build-linux`: failure blocks PRs via branch protection.

For `build-windows` (`continue-on-error: true`): GitHub records the workflow run
as "succeeded" regardless of step failures, so the `CI / build-windows` status check
always reports green at the workflow level. The actual step-level failure IS visible
in the job log. This is why `build-windows` MUST NOT be added as a required check.

---

## Annotated Workflow Skeleton

```yaml
name: CI

on:
  push:
    branches: [master]
  pull_request:
    branches: [master]

jobs:

  build-macos:
    runs-on: macos-latest
    timeout-minutes: 30
    steps:
      - uses: actions/checkout@v4
      - name: Build
        run: ./configure --enable-named-semaphores && make
      - name: Add to PATH
        run: echo "$GITHUB_WORKSPACE/src" >> $GITHUB_PATH
      - name: Test
        run: cd tests && bash dw-tests.sh

  build-linux:
    runs-on: ubuntu-latest
    timeout-minutes: 30
    steps:
      - uses: actions/checkout@v4
      - name: Build
        run: ./configure && make
      - name: Add to PATH
        run: echo "$GITHUB_WORKSPACE/src" >> $GITHUB_PATH
      - name: Test
        run: cd tests && bash dw-tests.sh

  build-windows:
    runs-on: windows-latest
    timeout-minutes: 30
    continue-on-error: true
    defaults:
      run:
        shell: msys2 {0}
    steps:
      - uses: actions/checkout@v4
      - uses: msys2/setup-msys2@v2
        with:
          msystem: MINGW64
          update: true
          cache: true
          install: >-
            mingw-w64-x86_64-gcc
            make
      - name: Build
        run: ./configure && make
      - name: Add to PATH
        run: echo "$GITHUB_WORKSPACE/src" >> $GITHUB_PATH
      - name: Test
        run: cd tests && bash dw-tests.sh
```

---

## Invariants

1. **Job names are stable**: `build-macos`, `build-linux`, `build-windows` MUST NOT
   be renamed without updating branch protection rules.
2. **Workflow name is stable**: Top-level `name: CI` MUST NOT be changed without
   updating branch protection rules.
3. **No secrets**: The workflow MUST NOT reference any secrets (not even `GITHUB_TOKEN`
   explicitly — it is injected automatically).
4. **No artifact upload**: This workflow does not upload binaries or test logs. CI
   job logs are the sole output.
