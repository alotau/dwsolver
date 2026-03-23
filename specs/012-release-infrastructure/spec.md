# Feature Specification: Release Infrastructure

**Feature Branch**: `012-release-infrastructure`  
**Created**: 2026-03-22  
**Status**: Draft  
**Input**: User description: "Release infrastructure for dwsolver C/Autotools project"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Maintainer Produces a Complete Source Tarball (Priority: P1)

A project maintainer runs `make distcheck` from the repository root after completing a change.  The command builds a source tarball, unpacks it into a clean temporary directory, configures, builds, and runs all tests from that clean tree — with zero file-not-found errors.  Every file that any downstream user needs to build and verify the project (source, tests, specs, example inputs, CI configuration, documentation) is present in the tarball.

**Why this priority**: Without a passing `make distcheck` there is no releasable artifact.  All other stories depend on this gate passing cleanly.

**Independent Test**: Run `make distcheck` on a clean checkout.  The command exits 0 and places `dwsolver-<version>.tar.gz` in the working directory.  Unpacking and building from the tarball alone produces a working binary and all tests pass.

**Acceptance Scenarios**:

1. **Given** a clean working tree at version 1.2.1, **When** a maintainer runs `make distcheck`, **Then** the command exits successfully and produces `dwsolver-1.2.1.tar.gz`.
2. **Given** the produced tarball, **When** it is unpacked and built on a system where only GLPK is installed (no other dwsolver files present), **Then** `make check` within the unpacked tree passes with no failures.
3. **Given** the produced tarball, **When** it is inspected, **Then** files from `tests/`, `examples/`, `specs/`, `architecture/`, `Dockerfile`, `.github/workflows/`, and `README.md` are all present.

---

### User Story 2 - Automated GitHub Release on Version Tag Push (Priority: P2)

A maintainer pushes a version tag (e.g., `v1.2.1`) to the repository.  Within minutes, a GitHub Actions job runs `make distcheck` on Linux, then creates a GitHub Release for that tag and attaches the resulting source tarball as a release asset.  The maintainer does not need to manually build or upload anything.

**Why this priority**: Automating the release eliminates manual upload steps and ensures every published release has a verifiably tested tarball attached.

**Independent Test**: Push a `v*` tag to the repository; verify a GitHub Release appears with the tarball attached and a `make distcheck` run in the workflow log that exited 0.

**Acceptance Scenarios**:

1. **Given** a version tag is pushed to the `main` branch, **When** the release workflow runs, **Then** `make distcheck` completes without errors and the tarball is uploaded to a new GitHub Release.
2. **Given** a non-tag push or a PR, **When** CI runs, **Then** the release workflow does NOT trigger (only the existing CI workflows run).
3. **Given** `make distcheck` fails on the tag push (e.g., due to a missing file), **When** the workflow completes, **Then** no GitHub Release is created and the workflow exits with a non-zero status.

---

### User Story 3 - Correct Library Soname Versioning (Priority: P3)

A maintainer preparing a release updates the libtool version-info triple in `src/Makefile.am` according to documented rules before cutting the tag.  The resulting shared library soname reflects the interface compatibility level, enabling package maintainers and downstream users to install multiple ABI-incompatible versions side-by-side if needed.

**Why this priority**: Correct soname versioning is required for safe library packaging across Linux distributions and is a prerequisite for downstream binary compatibility guarantees.  It can be addressed independently of the tarball and release workflow.

**Independent Test**: Build `libdwsolver` and verify the soname encoded in the shared library matches the expected value for the current version-info triple (e.g., `libdwsolver.so.0` when `CURRENT=0, AGE=0`).

**Acceptance Scenarios**:

1. **Given** the current version-info `0:0:0`, **When** `libdwsolver.so` is built, **Then** the soname is `libdwsolver.so.0` and the install name includes the CURRENT interface number.
2. **Given** a new release that adds public API functions without removing any, **When** the maintainer follows the documented update procedure, **Then** CURRENT and AGE are each incremented by 1 and REVISION is reset to 0.
3. **Given** a release that changes or removes an existing public API function, **When** the maintainer follows the documented update procedure, **Then** CURRENT is incremented, AGE is reset to 0, and REVISION is reset to 0.

---

### User Story 4 - Release Process Documented in README (Priority: P4)

A contributor or new maintainer reads the project README and finds a concise "How to cut a release" section that walks them through the complete release workflow: updating version numbers, running `make distcheck`, tagging and pushing to trigger the automated release, and updating the version-info triple when the API changes.

**Why this priority**: Documentation ensures the release process is reproducible by any authorized maintainer without tribal knowledge.

**Independent Test**: A person unfamiliar with the project follows only the README instructions and successfully produces a release from start to finish.

**Acceptance Scenarios**:

1. **Given** the README, **When** a maintainer reads the release section, **Then** they can identify all four steps: version bump, distcheck gate, tagging convention, and library version rules.
2. **Given** the README, **When** a maintainer follows the documented steps, **Then** `make distcheck` produces a passing tarball and a pushed `v*` tag triggers the release workflow.

---

### Edge Cases

- What happens when `make distcheck` is run but a required file is missing from `EXTRA_DIST`?  The build inside the clean tree fails with a clear file-not-found error; `make distcheck` exits non-zero and no tarball is produced.
- What happens when a tag is pushed that does not match `v*` (e.g., `release-1.2.1`)?  The release workflow does not trigger; only the regular CI runs.
- What happens when the release workflow runs but the GitHub token lacks permission to create releases?  The workflow fails at the release-creation step with an authentication error; the job exit code is non-zero.
- What happens when `make distcheck` is run on a system that does not have GLPK installed?  The configure step inside the clean tree fails immediately with a descriptive error; `make distcheck` exits non-zero.
- What happens when the version-info triple is not updated before a release that changes the public API?  The soname of the installed library is unchanged, causing silent ABI conflicts for programs that depend on the old interface.
- What happens when a maintainer pushes the same `v*` tag a second time (e.g., to fix a botched release)?  The release workflow fails at the release-creation step because the GitHub Release already exists; no asset is overwritten.  The maintainer must manually delete the tag and the GitHub Release, then re-push the corrected tag.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: `make distcheck` MUST exit 0 and produce a `dwsolver-<version>.tar.gz` tarball that contains all source files, test files, example input files, specification documents, CI workflow definitions, architecture documents, the Dockerfile, and README.
- **FR-002**: The source tarball MUST be self-contained: a user who unpacks it on a system with only GLPK available MUST be able to build the project and pass `make check` without any additional files.
- **FR-003**: `EXTRA_DIST` in the top-level `Makefile.am` MUST include: `tests/`, `specs/`, `architecture/`, `Dockerfile`, `.github/workflows/`, and `README.md`.
- **FR-004**: `EXTRA_DIST` in `src/Makefile.am` MUST include all non-installed header files that are part of the source distribution.
- **FR-005**: A GitHub Actions release workflow MUST trigger on any tag push matching the pattern `v*` and MUST NOT trigger on branch pushes or pull requests.  If the tag already has an associated GitHub Release, the workflow MUST fail (exit non-zero) rather than overwriting or silently skipping the release; no asset upload occurs in this case.
- **FR-006**: The release workflow MUST run `make distcheck` on Linux (Ubuntu latest) before creating the GitHub Release; if `make distcheck` fails, no release is created.
- **FR-007**: The release workflow MUST upload the `dwsolver-<version>.tar.gz` tarball as a release asset on the automatically created GitHub Release.
- **FR-007a**: The GitHub Release MUST be published immediately (not as a draft) when the release workflow succeeds; no manual promotion step is required.
- **FR-008**: The libtool version-info triple in `src/Makefile.am` MUST follow the standard update rules: increment REVISION on bug-fix releases, increment CURRENT and AGE (reset REVISION) when adding new public APIs, increment CURRENT and reset both AGE and REVISION when removing or changing public APIs.
- **FR-009**: The README MUST contain a "How to cut a release" section that describes, at minimum: updating the version in `configure.ac`, updating the version-info triple when the API changes, running `make distcheck` as the release gate, and pushing a `v<major>.<minor>.<patch>` tag to trigger the automated release.
- **FR-010**: The release workflow MUST install GLPK and any other build dependencies on the CI runner before executing `make distcheck`, and MUST export `GLPK_CFLAGS` and `GLPK_LIBS` as environment variables (not via `DISTCHECK_CONFIGURE_FLAGS`) so the configure invocation inside the clean distcheck tree picks them up correctly — consistent with the existing CI workflows.
- **FR-011**: All GitHub Actions used in the release workflow (both first-party `actions/*` and third-party) MUST be pinned to a specific commit SHA (not a floating tag); each pinned reference MUST include a comment identifying the corresponding version tag.

### Key Entities

- **Source tarball**: The versioned `.tar.gz` archive produced by `make dist` / `make distcheck`, containing everything needed to build and test dwsolver from scratch.
- **Version-info triple**: The `CURRENT:REVISION:AGE` libtool parameter that controls the shared library soname and ABI compatibility guarantees.
- **Release tag**: A git tag matching `v<major>.<minor>.<patch>` that serves as the trigger for the automated release workflow.
- **GitHub Release**: The release record created on GitHub, associated with a tag, containing the release notes and the attached source tarball asset.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: `make distcheck` exits 0 in a clean checkout of the repository with no manual preparation beyond installing GLPK.
- **SC-002**: The source tarball produced by `make distcheck` contains 100% of the files listed in `EXTRA_DIST` across all `Makefile.am` files, verified by `tar tzf dwsolver-<version>.tar.gz`.
- **SC-003**: Pushing a `v*` tag to the repository results in a GitHub Release being created within 15 minutes, with the tarball attached as a downloadable asset.
- **SC-004**: The `make distcheck` step within the release workflow exits 0 before the release is created; workflow logs show zero test failures.
- **SC-005**: The libtool version-info triple in `src/Makefile.am` is set to `0:0:0` for the current 1.2.1 release, and the README documents the rules for updating it in future releases.
- **SC-006**: A maintainer following only the README's release section can produce a tagged GitHub Release without consulting any other documentation.

## Assumptions

- The project remains at version 1.2.1 for this feature; the version number in `configure.ac` is not changed as part of this work.
- The libtool version-info stays at `0:0:0` for the current release because no public API has changed yet; the infrastructure documents when and how to change it.
- The release workflow runs only on Linux; macOS and Windows builds are covered by existing CI workflows and are not required for the release gate.
- GitHub Actions permissions (`contents: write`) are granted to the `GITHUB_TOKEN` for the release workflow to create and upload releases; the repository settings permit this.
- `make distcheck` may require the same `GLPK_CFLAGS`/`GLPK_LIBS` workaround used in the existing Linux CI workflow (since Ubuntu's libglpk-dev does not ship a `.pc` file).
- `.github/` is a hidden directory; automake supports listing it in `EXTRA_DIST` and it will be included in the tarball when referenced explicitly.

## Clarifications

### Session 2026-03-22

- Q: When a maintainer pushes the same version tag a second time, what should the release workflow do? → A: Fail the workflow — the release already exists; maintainer must manually delete it first
- Q: Should the `softprops/action-gh-release` action version be pinned to a commit SHA or kept as a floating tag? → A: Pin to a commit SHA (supply-chain security); annotate with version tag comment
- Q: Should action SHA-pinning apply to all actions (including first-party `actions/*`) or only third-party actions? → A: All actions including `actions/checkout`; consistent pinning eliminates the whole category of supply-chain risk
- Q: Should the GitHub Release be published immediately or created as a draft for manual promotion? → A: Publish immediately — `make distcheck` is already the gate; no manual promotion step needed
- Q: Should `GLPK_CFLAGS`/`GLPK_LIBS` be passed inside `DISTCHECK_CONFIGURE_FLAGS` or exported as environment variables? → A: Export as environment variables (consistent with `ci-linux.yml`; avoids word-splitting risk)
