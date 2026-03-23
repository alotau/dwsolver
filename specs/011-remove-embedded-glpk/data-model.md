# Data Model: Change Inventory

**Feature**: 011-remove-embedded-glpk  
**Phase**: 1 — Design  
**Date**: 2026-03-22

This document is the complete inventory of every file that will be **deleted**,
**modified**, or **created** as part of this feature.  It serves as the
authoritative reference for task generation.

---

## File Counts (pre-change)

| Category | Count |
|---|---|
| Total tracked files in `src/` | 159 |
| `glp*.c` / `glp*.h` files to delete | 125 |
| `src/amd/` files to delete | 14 |
| `src/colamd/` files to delete | 2 |
| **Total deletions** | **141** |
| **Files remaining in `src/`** | **18** |

---

## A. Files to Delete

### A1. Embedded GLPK source and headers (`src/glp*`)

All 125 files matching `src/glp[^k]*` — i.e., every file whose name begins
with `glp` except `glpk.h` (which is the public header that will be provided
by the system installation instead).  Specifically:

- `glpapi01.c` … `glpapi19.c` (19 files)
- `glpavl.c`, `glpavl.h`
- `glpbfd.c`, `glpbfd.h`
- `glpbfx.c`, `glpbfx.h`
- `glpcpx.c`
- `glpdmp.c`, `glpdmp.h`
- `glpdmx.c`
- `glpenv.h`, `glpenv01.c` … `glpenv08.c` (1 header + 8 sources)
- `glpfhv.c`, `glpfhv.h`
- `glpgmp.c`, `glpgmp.h`
- `glphbm.c`, `glphbm.h`
- `glpini01.c`, `glpini02.c`
- `glpios.h`, `glpios01.c` … `glpios12.c` (1 header + 12 sources)
- `glpipm.c`, `glpipm.h`
- `glplib.h`, `glplib01.c`, `glplib02.c`, `glplib03.c`
- `glplpf.c`, `glplpf.h`
- `glplpx01.c`, `glplpx02.c`, `glplpx03.c`
- `glpluf.c`, `glpluf.h`
- `glplux.c`, `glplux.h`
- `glpmat.c`, `glpmat.h`
- `glpmpl.h`, `glpmpl01.c` … `glpmpl06.c` (1 header + 6 sources)
- `glpmps.c`
- `glpnet.h`, `glpnet01.c` … `glpnet09.c` (1 header + 9 sources)
- `glpnpp.h`, `glpnpp01.c` … `glpnpp05.c` (1 header + 5 sources)
- `glpqmd.c`, `glpqmd.h`
- `glprgr.c`, `glprgr.h`
- `glprng.h`, `glprng01.c`, `glprng02.c`
- `glpscf.c`, `glpscf.h`
- `glpscl.c`
- `glpsdf.c`
- `glpspm.c`, `glpspm.h`
- `glpspx.h`, `glpspx01.c`, `glpspx02.c`
- `glpsql.c`, `glpsql.h`
- `glpssx.h`, `glpssx01.c`, `glpssx02.c`
- `glpstd.h`
- `glptsp.c`, `glptsp.h`
- `glpapi.h` (internal GLPK umbrella header)
- **`glpk.h`**: this is the public GLPK API header — it will be **deleted from
  `src/`** because the system-installed GLPK provides it; `#include <glpk.h>`
  (angle-bracket form) in all DW source files continues to work once the
  system header path is passed via `$(GLPK_CFLAGS)`.

### A2. AMD library (`src/amd/`)

All 14 files: `amd_1.c`, `amd_2.c`, `amd_aat.c`, `amd_control.c`,
`amd_defaults.c`, `amd_dump.c`, `amd_info.c`, `amd_order.c`,
`amd_post_tree.c`, `amd_postorder.c`, `amd_preprocess.c`, `amd_valid.c`,
`amd.h`, `amd_internal.h`.

### A3. COLAMD library (`src/colamd/`)

All 2 files: `colamd.c`, `colamd.h`.

---

## B. Files to Modify

### B1. `configure.ac`

| Change | Detail |
|--------|--------|
| Remove | `AC_CHECK_LIB([m], [exp])` remains (math lib still needed) |
| Remove | Manual `AC_DEFINE` blocks for `HAVE_GMP`, `HAVE_ZLIB`, `HAVE_LTDL`, `HAVE_DLFCN` that only existed to enable optional GLPK built-in features. These are no longer needed because the system GLPK already encodes those choices. |
| Add | `PKG_CHECK_MODULES([GLPK], [glpk >= 4.65])` — sets `GLPK_CFLAGS` and `GLPK_LIBS` |
| Remove | `AC_ARG_WITH(gmp, ...)`, `AC_ARG_WITH(zlib, ...)`, `AC_ARG_ENABLE(dl, ...)`, `AC_ARG_ENABLE(odbc, ...)`, `AC_ARG_ENABLE(mysql, ...)` — all only configured the embedded GLPK optional features |
| Retain | `AC_ARG_ENABLE(named_semaphores, ...)` and `AC_ARG_ENABLE(recursive_mutex, ...)` — these configure dwsolver's own threading behaviour, unrelated to GLPK |
| Retain | `AC_CHECK_LIB([pthread], ...)` |
| Retain | `AC_CHECK_HEADER([sys/time.h], ...)`, `AC_CHECK_FUNC([gettimeofday], ...)` |

### B2. `src/Makefile.am`

| Change | Detail |
|--------|--------|
| Remove | Entire `GLPK_SOURCES` variable block (all `glpapi*.c`, `glp*.c`, `amd/*.c`, `colamd/*.c`) |
| Modify | `libdwsolver_la_SOURCES` — remove `$(GLPK_SOURCES)` |
| Add | `libdwsolver_la_CPPFLAGS = $(AM_CPPFLAGS) -DDWSOLVER_BUILDING_LIB $(GLPK_CFLAGS)` |
| Add | `libdwsolver_la_LIBADD = $(GLPK_LIBS)` |
| Modify | `EXTRA_DIST` — remove all `glp*.h`, `amd/amd.h`, `amd/amd_internal.h`, `colamd/colamd.h` entries |
| Retain | All `DW_CORE_SOURCES`, `dwsolver_SOURCES`, `nobase_include_HEADERS` |

### B3. `src/dw_solver.c` — 5 `lpx_*` call sites

| Location | Old | New |
|----------|-----|-----|
| ~line 184 (read master LP) | `original_master_lp = lpx_read_cpxlp(globals->master_name);` | `original_master_lp = glp_create_prob(); if (glp_read_lp(original_master_lp, NULL, globals->master_name) != 0) { glp_delete_prob(original_master_lp); original_master_lp = NULL; }` |
| ~line 434 | `lpx_write_cpxlp(master_lp, "pre_master.cpxlp");` | `glp_write_lp(master_lp, NULL, "pre_master.cpxlp");` |
| ~line 477 | `lpx_write_cpxlp(master_lp, local_buffer);` | `glp_write_lp(master_lp, NULL, local_buffer);` |
| ~line 493 | `lpx_write_cpxlp(master_lp, local_buffer);` | `glp_write_lp(master_lp, NULL, local_buffer);` |
| ~line 645 | `lpx_write_cpxlp(master_lp, local_buffer);` | `glp_write_lp(master_lp, NULL, local_buffer);` |

### B4. `src/dw_subprob.c` — 1 `lpx_*` call site

| Location | Old | New |
|----------|-----|-----|
| ~line 103 (read subproblem LP) | `lp = lpx_read_cpxlp(my_data->infile_name);` | `lp = glp_create_prob(); if (glp_read_lp(lp, NULL, my_data->infile_name) != 0) { glp_delete_prob(lp); lp = NULL; }` |

### B5. `src/dw_phases.c` — 4 `lpx_*` call sites (2 logical, 4 lines)

| Location | Old | New |
|----------|-----|-----|
| ~lines 228–229 | `if (simplex_iterations < lpx_get_int_parm(master_lp, LPX_K_ITCNT)) { simplex_iterations = lpx_get_int_parm(master_lp, LPX_K_ITCNT); }` | `if (simplex_iterations < glp_get_it_cnt(master_lp)) { simplex_iterations = glp_get_it_cnt(master_lp); }` |
| ~lines 460–461 | same pattern | same replacement |

### B6. `dwsolver.pc.in`

Add one line:
```
Requires.private: glpk
```

### B7. `Dockerfile`

| Stage | Change |
|-------|--------|
| `builder` | Add `libglpk-dev` to the `apt-get install` line |
| `runner` | Add `apt-get install -y libglpk40` step (runtime `.so` only) |

### B8. `.github/workflows/ci-linux.yml`

Add a `Install dependencies` step before `Build` in **all three** jobs
(`linux`, `linux-asan-ubsan`, `linux-tsan`):
```yaml
- name: Install dependencies
  run: apt-get update && apt-get install -y libglpk-dev
  # or: sudo apt-get update && sudo apt-get install -y libglpk-dev
```

### B9. `.github/workflows/ci-macos.yml`

Add a `Install dependencies` step before `Build`:
```yaml
- name: Install dependencies
  run: brew install glpk
```

### B10. `README.md`

Add a **Dependencies** section documenting:
- GLPK ≥ 4.65 is required
- macOS install: `brew install glpk`
- Ubuntu/Debian install: `sudo apt-get install libglpk-dev`
- Fedora/RHEL install: `sudo dnf install glpk-devel`
- Build from source note with minimum version

Add ABI note: callers of `libdwsolver` must link against the same GLPK
shared library that dwsolver was compiled against.

### B11. `.specify/memory/constitution.md`

Update the Technology Stack section:
- `**LP Backend**`: change from `Embedded GLPK 4.44 (thread-patched variant included in `src/`)` to `System GLPK ≥ 4.65 (external shared library, detected via pkg-config)`
- Remove / amend the `**No external runtime dependencies**` line

---

## C. Files to Create

None. All changes are modifications or deletions.

---

## D. Files Explicitly Not Changed

| File | Reason |
|------|--------|
| `src/dw.h` | Internal header; `#include <glpk.h>` resolves to system header via `$(GLPK_CFLAGS)` |
| `src/dw_solver.h` | Public API; already exposes zero GLPK types |
| `src/dw_globals.c` | Contains `glpk_mutex` definition — retained unchanged |
| `src/dw_support.c` | Contains `glpk_mutex` guards — retained unchanged |
| `src/dw_rounding.c` | Contains `glpk_mutex` guards — retained unchanged |
| `third-party/glpk/` | Preserved for provenance (all attribution files + patch) |
| `.github/workflows/ci-windows.yml` | Windows is out of scope for this feature |
| `.github/workflows/ci-docker.yml` | Triggers `docker build .`; Dockerfile fix (B7) covers it |
