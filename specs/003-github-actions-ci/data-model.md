# Data Model: GitHub Actions CI — Multi-Platform Build & Test

**Feature**: 003-github-actions-ci
**Date**: 2026-03-19
**Derived from**: research.md, spec.md FR-001 through FR-009

---

## Overview

This feature introduces a single YAML workflow file. It has no database entities,
no persisted state, and no runtime data model. The "entities" here are the structural
components of the workflow definition and the external configuration objects they
interact with.

---

## Entity 1: `ci.yml` — GitHub Actions Workflow File

**Path**: `.github/workflows/ci.yml`
**Purpose**: Defines when CI runs and what it does on each platform.
**Replaces**: `.github/workflows/build.yaml` (deleted as part of this implementation)

### Workflow-Level Fields

| Field | Value | Notes |
|-------|-------|-------|
| `name` | `CI` | Displayed in GitHub Actions tab |
| `on.push.branches` | `[master]` | Triggers on direct push to master |
| `on.pull_request.branches` | `[master]` | Triggers on PRs targeting master |

### Job: `build-macos`

| Field | Value |
|-------|-------|
| `runs-on` | `macos-latest` |
| `timeout-minutes` | `30` |
| `continue-on-error` | `false` (required check) |

**Steps**:
1. `actions/checkout@v4`
2. `./configure --enable-named-semaphores && make`
3. Add `$GITHUB_WORKSPACE/src` to `$GITHUB_PATH`
4. `cd tests && bash dw-tests.sh`

### Job: `build-linux`

| Field | Value |
|-------|-------|
| `runs-on` | `ubuntu-latest` |
| `timeout-minutes` | `30` |
| `continue-on-error` | `false` (required check) |

**Steps**:
1. `actions/checkout@v4`
2. `./configure && make`
3. Add `$GITHUB_WORKSPACE/src` to `$GITHUB_PATH`
4. `cd tests && bash dw-tests.sh`

**Note**: No pre-install step needed; `gcc` and `make` are pre-installed on Ubuntu runners.

### Job: `build-windows`

| Field | Value |
|-------|-------|
| `runs-on` | `windows-latest` |
| `timeout-minutes` | `30` |
| `continue-on-error` | `true` (non-blocking, FR-005) |

**Steps**:
1. `actions/checkout@v4`
2. `msys2/setup-msys2@v2` with `msystem: MINGW64`, packages: `mingw-w64-x86_64-gcc make`, `cache: true`
3. (MSYS2 shell) `./configure && make`
4. (MSYS2 shell) Add `src` to PATH
5. (MSYS2 shell) `cd tests && bash dw-tests.sh`

---

## Entity 2: Branch Protection Rule (External to Codebase)

**Location**: GitHub repository settings → Branches → `master`
**Purpose**: Enforces that required CI checks pass before any PR can be merged.
**Configuration**: Set by a repository maintainer in the GitHub web UI.

| Required Status Check | Job Name | Must Pass? |
|-----------------------|----------|-----------|
| CI / build-macos | `build-macos` | Yes |
| CI / build-linux | `build-linux` | Yes |
| CI / build-windows | `build-windows` | No (informational) |

**Note**: The status check name format is `{workflow name} / {job id}`, so if the
workflow `name` is `CI` and the job is `build-macos`, the check name is `CI / build-macos`.
This exact name must be entered in branch protection settings.

---

## Entity 3: `build.yaml` — To Be Deleted

**Path**: `.github/workflows/build.yaml`
**Status**: Superseded by `ci.yml`
**Action**: Delete this file as part of the implementation.

**Why delete rather than keep**: Having two CI workflow files with overlapping triggers
on `master` would run both on every push to master. `build.yaml` is strictly worse
(macOS-only, missing `--enable-named-semaphores`, outdated action versions). Keeping
it creates confusion and wastes CI minutes.

---

## Dependency Map

```
Feature Branch PR
        │
        ▼
GitHub Actions Trigger (pull_request → master)
        │
   ┌────┴─────────────┐
   │                  │                  │
build-macos     build-linux       build-windows
 (required)      (required)      (non-blocking)
   │                  │
   └────────┬─────────┘
            │
     Branch Protection Rule
            │
     "Merge" enabled only if both pass
```

---

## State Transitions

The workflow has three outcome states per job:

| State | Meaning | Effect on PR |
|-------|---------|-------------|
| `success` | Build and all tests passed | Contributes to merge eligibility |
| `failure` | Build failed or ≥1 test failed | Blocks merge (required checks only) |
| `cancelled` | Job was cancelled manually | Treated as not-yet-run; does not block |

For `build-windows` with `continue-on-error: true`: GitHub marks the job as
"success" from the perspective of the workflow run even if the underlying steps
fail, so branch protection rules treating the Windows check as required would
always see it as passing. For this reason, the Windows job is informational only
and should NOT be added as a required status check.
