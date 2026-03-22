# Build System Contract: GLPK Dependency

**Feature**: 011-remove-embedded-glpk  
**Contract type**: Build system (configure.ac, Makefile.am, pkg-config)  
**Date**: 2026-03-22

This document captures the exact before/after changes to the build-system
interface.  It is the normative specification that the implementation tasks
must reproduce faithfully.

---

## 1. `configure.ac` — GLPK detection block

### 1a. Additions

Insert after the existing `AC_CHECK_LIB([pthread], ...)` line:

```m4
dnl Require system GLPK >= 4.65
PKG_CHECK_MODULES([GLPK], [glpk >= 4.65],
  [],
  [AC_MSG_ERROR([GLPK >= 4.65 is required. Install libglpk-dev (Debian/Ubuntu), glpk-devel (Fedora/RHEL), or brew install glpk (macOS).])])
```

`PKG_CHECK_MODULES` sets `GLPK_CFLAGS` and `GLPK_LIBS` automatically.
These are substituted into `Makefile.am` (see §2).

### 1b. Removals from `configure.ac`

Remove all of the following blocks (lines 32–177 in the current file).
They exclusively configured optional features of the embedded GLPK build;
none of them affect dwsolver's own code:

| Block | Lines (approx) | Reason for removal |
|-------|--------|--------------------|
| `AC_ARG_WITH(gmp, ...)` | 32–39 | Configured embedded GLPK's optional GNU MP support |
| `AC_ARG_WITH(zlib, ...)` | 41–48 | Configured embedded GLPK's optional zlib support |
| `AC_ARG_ENABLE(dl, ...)` | 50–56 | Configured embedded GLPK's dynamic-load (MathProg ODBC/MySQL) |
| `AC_ARG_ENABLE(odbc, ...)` | 58–65 | Configured embedded GLPK's MathProg ODBC support |
| `AC_ARG_ENABLE(mysql, ...)` | 67–74 | Configured embedded GLPK's MathProg MySQL support |
| `if test "$with_gmp" ...` | 118–122 | Applied `HAVE_GMP` and `-lgmp` for embedded GLPK |
| `if test "$with_zlib" ...` | 126–132 | Applied `HAVE_ZLIB` and `-lz` for embedded GLPK |
| `if test "$enable_dl" ...` | 134–145 | Applied `HAVE_LTDL`/`HAVE_DLFCN` for embedded GLPK |
| `case $host_os in darwin*` ... `esac` | 147–157 | Set `LIBIODBC`/`LIBODBC`/`LIBMYSQL` dylib names |
| `if test "$enable_odbc" ...` | 159–171 | Defined `ODBC_DLNAME` for MathProg ODBC |
| `if test "$enable_mysql" ...` | 173–184 | Defined `MYSQL_DLNAME` for MathProg MySQL |

> **Note**: The `darwin*` `case` block also added `-DUSE_NAMED_SEMAPHORES`
> unconditionally on macOS.  That behaviour is intentionally preserved by
> keeping the existing `AC_ARG_ENABLE(named_semaphores, ...)` code path
> (dwsolver-level option, unrelated to GLPK).

### 1c. Retained blocks

These lines remain unchanged:

- `AC_ARG_ENABLE(named_semaphores, ...)`
- `AC_ARG_ENABLE(recursive_mutex, ...)`
- `AC_CHECK_LIB([m], [exp])`
- `AC_CHECK_LIB([pthread], [sem_wait])`
- `AC_CHECK_HEADER([sys/time.h], ...)`
- `AC_CHECK_FUNC([gettimeofday], ...)`
- All `AC_MSG_CHECKING` / `AC_MSG_RESULT` for named semaphores and recursive mutex

---

## 2. `src/Makefile.am` — build rule changes

### 2a. Remove: `GLPK_SOURCES` variable

Delete the entire block from `GLPK_SOURCES = \` through `colamd/colamd.c`
(currently lines 7–113, 107 entries including the trailing `colamd/colamd.c`
with no backslash).

### 2b. Modify: `libdwsolver_la_SOURCES`

Before:
```makefile
libdwsolver_la_SOURCES  = $(GLPK_SOURCES) $(DW_CORE_SOURCES) dw_solver.c
```

After:
```makefile
libdwsolver_la_SOURCES  = $(DW_CORE_SOURCES) dw_solver.c
```

### 2c. Modify: `libdwsolver_la_CPPFLAGS`

Before:
```makefile
libdwsolver_la_CPPFLAGS = $(AM_CPPFLAGS) -DDWSOLVER_BUILDING_LIB
```

After:
```makefile
libdwsolver_la_CPPFLAGS = $(AM_CPPFLAGS) -DDWSOLVER_BUILDING_LIB $(GLPK_CFLAGS)
```

### 2d. Add: `libdwsolver_la_LIBADD`

Insert after the `libdwsolver_la_CPPFLAGS` line:

```makefile
libdwsolver_la_LIBADD   = $(GLPK_LIBS)
```

### 2e. Modify: `EXTRA_DIST`

Before (excerpt):
```makefile
EXTRA_DIST = amd/amd.h amd/amd_internal.h colamd/colamd.h \
dw_blas.h      glpapi.h  glpfhv.h  glplib.h  glpnet.h  glpspm.h \
dw.h           glpavl.h  glpgmp.h  glplpf.h  glpnpp.h  glpspx.h \
dw_phases.h    glpbfd.h  glphbm.h  glpluf.h  glpqmd.h  glpsql.h \
dw_rounding.h  glpbfx.h  glpios.h  glplux.h  glprgr.h  glpssx.h \
dw_subprob.h   glpdmp.h  glpipm.h  glpmat.h  glprng.h  glpstd.h \
dw_support.h   glpenv.h  glpmpl.h  glpscf.h  glptsp.h  glpk.h
```

After (remove all `glp*.h`, `amd/`, and `colamd/` entries):
```makefile
EXTRA_DIST = \
dw_blas.h \
dw.h \
dw_phases.h \
dw_rounding.h \
dw_subprob.h \
dw_support.h
```

---

## 3. `dwsolver.pc.in` — pkg-config metadata

### 3a. Before

```
prefix=@prefix@
exec_prefix=@exec_prefix@
libdir=@libdir@
includedir=@includedir@

Name: dwsolver
Description: Dantzig-Wolfe decomposition LP solver library
Version: @PACKAGE_VERSION@
Cflags: -I${includedir}
Libs: -L${libdir} -ldwsolver
```

### 3b. After

```
prefix=@prefix@
exec_prefix=@exec_prefix@
libdir=@libdir@
includedir=@includedir@

Name: dwsolver
Description: Dantzig-Wolfe decomposition LP solver library
Version: @PACKAGE_VERSION@
Requires.private: glpk
Cflags: -I${includedir}
Libs: -L${libdir} -ldwsolver
```

**Rationale**: `Requires.private` (not `Requires`) because `dw_solver.h` exposes
zero GLPK types.  Consumers of `libdwsolver` do not need `glpk.h` in their
own includes.  `Requires.private` ensures transitive link flags appear in
`pkg-config --libs --static dwsolver` output but not in the default
`pkg-config --libs dwsolver` output.

---

## 4. `Dockerfile` — system dependency additions

### 4a. Builder stage (`ubuntu:24.04`)

Add `libglpk-dev` to the builder's `apt-get install` line so the build can
find `glpk.h` and `libglpk.so` at configure time and link time.

Before (excerpt):
```dockerfile
RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    libtool \
    pkg-config
```

After:
```dockerfile
RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    automake \
    libtool \
    pkg-config \
    libglpk-dev
```

### 4b. Runner stage

Add a step that installs only the GLPK runtime shared library into the
minimal runner image.  The exact package name depends on the runner base:

| Base image | Package to install |
|------------|--------------------|
| `ubuntu:24.04` or `ubuntu:22.04` | `libglpk40` |
| `ubuntu:20.04` | `libglpk40` |
| `debian:bookworm` | `libglpk40` |
| `alpine:3.x` | `glpk` |

```dockerfile
RUN apt-get update && apt-get install -y --no-install-recommends libglpk40 \
 && rm -rf /var/lib/apt/lists/*
```

---

## 5. CI workflow changes

### 5a. `.github/workflows/ci-linux.yml`

Three jobs require identical treatment: `linux`, `linux-asan-ubsan`,
`linux-tsan`.  Before the existing `Build` step in each:

```yaml
- name: Install dependencies
  run: |
    sudo apt-get update
    sudo apt-get install -y libglpk-dev
```

### 5b. `.github/workflows/ci-macos.yml`

Before the existing `Build` step:

```yaml
- name: Install dependencies
  run: brew install glpk
```

### 5c. `.github/workflows/ci-windows.yml`

**No change** — Windows cross-build support is out of scope for this feature.

### 5d. `.github/workflows/ci-docker.yml`

**No change** — Docker CI triggers `docker build .`; the Dockerfile change
(§4) covers it.

---

## 6. Compatibility matrix

| Platform | GLPK package | Min version in distro |
|----------|--------------|-----------------------|
| Ubuntu 20.04 LTS | `libglpk-dev` | 4.65 ✓ |
| Ubuntu 22.04 LTS | `libglpk-dev` | 5.0 ✓ |
| Ubuntu 24.04 LTS | `libglpk-dev` | 5.0 ✓ |
| macOS (Homebrew) | `glpk` | 5.0 ✓ |
| Fedora 38+ | `glpk-devel` | 5.0 ✓ |
| Debian Bookworm | `libglpk-dev` | 5.0 ✓ |

All supported platforms provide ≥ 4.65; the `PKG_CHECK_MODULES` floor of
`glpk >= 4.65` will always pass on these systems.
