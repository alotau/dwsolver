# Data Model: Containerized Build for dwsolver

**Feature**: 004 — Containerized Build
**Date**: 2026-03-19

This feature introduces no new data types, database entities, or persistent
application state. The "data model" for a containerized build is the filesystem
structure of the Docker image layers and the mount conventions through which the
running container accesses user data.

---

## Image Layer Model

### Builder Stage (`builder`)

| Layer | Content | Estimated Size |
|-------|---------|----------------|
| Base OS | `ubuntu:24.04` | ~77 MB |
| Toolchain | `build-essential`, `automake` (via `apt-get`) | ~200 MB |
| Source tree | All files from build context (after `.dockerignore` filtering) | ~5 MB |
| Build artifacts | `.o` files, binary `src/dwsolver` (generated during `make`) | ~10 MB |
| **Builder total** | (discarded — not in final image) | ~292 MB |

### Runner Stage (`runner`)

| Layer | Content | Estimated Size |
|-------|---------|----------------|
| Base OS | `ubuntu:24.04` (fresh copy — no toolchain) | ~77 MB |
| Binary | `/usr/local/bin/dwsolver` (copied from builder) | ~1 MB |
| **Runner total** | **Final image size** | **~78 MB** ✅ (under 200 MB requirement) |

---

## Build Context Model

The build context is the set of files sent to the Docker daemon at `docker build`
time. It is governed by `.dockerignore`.

### Included (sent to daemon)

| Path | Reason |
|------|--------|
| `src/*.c`, `src/*.h` | Compiled in the builder stage |
| `src/Makefile.am`, `src/Makefile.in` | Required for autotools build |
| `aclocal.m4`, `configure`, `config.h.in`, `ltmain.sh`, `Makefile.in` | Pre-generated autotools files; `touch`-ed before `make` |
| `configure.ac`, `Makefile.am` | Top-level autotools inputs |
| `examples/` | Needed for CI smoke test volume mount |
| `tests/` | Available for extended smoke tests |

### Excluded (`.dockerignore`)

| Pattern | Reason |
|---------|--------|
| `.git/` | Can be 100+ MB of history; irrelevant to build |
| `*.o`, `*.a` | Stale local build artifacts |
| `src/dwsolver` | Stale locally compiled binary |
| `config.log`, `config.status` | Local configure artifacts |
| `autom4te.cache/` | Autoconf cache; irrelevant |
| `libtool`, `stamp-h1` | Local build state |
| `tmp/` | Temporary files |
| `**/out_terminal` | Test output |
| `**/relaxed_solution` | Test output |
| `**/rs_sorted` | Test output |
| `**/done.cpxlp`, `**/pre_master.cpxlp` | Solver intermediate files |
| `*.Po` | GCC dependency tracking files |

---

## Volume Mount Convention

The running container has no persistent storage. User data is provided exclusively
via Docker volume mounts.

```
Host filesystem                     Container filesystem
─────────────────                   ────────────────────
/my/problem/                   →    /data/
  ├── guidefile                         ├── guidefile
  ├── master.cplex                      ├── master.cplex
  ├── sub1.cplex                        ├── sub1.cplex
  └── ...                               └── ...
```

### Mount Properties

| Property | Value |
|----------|-------|
| Host path | Any absolute path containing the guide file and referenced data files |
| Container path | `/data` (conventional; user may choose any path) |
| Access mode | Read-write — dwsolver writes output files (e.g., `relaxed_solution`) to the working directory |
| Requirement | All files referenced by the guide file MUST be within the mounted directory |

### Guide File Path Resolution

dwsolver resolves data file paths **relative to the guide file's location**.
When using the container, all problem files must reside in the same directory
as the guide file, which must itself be within the mounted volume.

**Correct**: All files in the mounted directory:
```sh
-v /my/problem:/data   →   -g /data/guidefile
```

**Incorrect**: Guide file in mounted directory, data files referencing absolute
host paths — data files are not accessible from inside the container.

---

## File Delivery Model

### Files Created by This Feature

| File | Location | Content |
|------|----|------|
| `Dockerfile` | Repository root | Two-stage build definition |
| `.dockerignore` | Repository root | Build context exclusion list |

### Files Modified by This Feature

| File | Location | Change |
|------|----------|--------|
| `.github/workflows/ci.yml` | Repository root | Add `docker-build` job |

### Files NOT Changed

- No `src/*.c` or `src/*.h` files — all source fixes are on `master` (Spec 002)
- No `tests/dw-tests.sh` — test script unchanged
- No `Makefile.am` or autotools files
