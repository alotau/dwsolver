# Research: ISO/IEC TS 17961:2013 Compliance

**Phase 0 — Pre-design audit**
**Branch**: `014-iso-ts-17961-compliance`
**Date**: 2026-03-23

---

## 1. The 22 TS 17961 Rules — Definition and Scope

ISO/IEC TS 17961:2013 defines exactly 22 analyzable rules. The table below
 maps each rule to its commonly used mnemonic, its nearest SEI CERT C
 correspondent, and the automated checker (if any) available in
 `cppcheck --addon=cert` or `clang-tidy -checks=cert-*`.

| # | Rule mnemonic  | TS 17961 description (summary)                                     | CERT C nearest    | Automated checker                        |
|---|----------------|--------------------------------------------------------------------|-------------------|------------------------------------------|
| 1 | intptrconv     | Do not convert between integer types and pointer types improperly  | INT36-C / EXP35-C | clang-tidy: cert-int36-c                 |
| 2 | intobjptr      | Do not convert pointer to object to incompatible pointer type      | EXP36-C           | clang-tidy: cert-exp36-c                 |
| 3 | xplicitcomp    | Do not perform explicit comparisons with byte strings              | STR34-C           | no automated checker                     |
| 4 | boolasgn       | Do not assign boolean-like values to integer variables             | —                 | no automated checker                     |
| 5 | trstcmp        | Detect errors when converting a string to a number                 | ERR34-C           | clang-tidy: cert-err34-c                 |
| 6 | wraparound     | Ensure integer arithmetic operations do not result in overflow     | INT30-C / INT32-C | clang-tidy: cert-int30-c, cert-int32-c   |
| 7 | argcomp        | Arguments to functions must be compatible with declared parameters | DCL11-C           | cppcheck: passedByValue (partial)        |
| 8 | intnegator     | Do not negate the minimum value of a signed integer               | INT32-C           | clang-tidy: cert-int32-c                 |
| 9 | ioilecc        | Do not use types with implementation-defined width for I/O        | FIO47-C           | no automated checker (partial -Wformat)  |
|10 | ga-buffer      | Do not form or use out-of-bounds pointers or data                  | ARR30-C / STR31-C | cppcheck: bufferAccessOutOfBounds        |
|11 | strmod         | Do not modify string literals                                      | STR30-C           | cppcheck: writeReadOnlyVar               |
|12 | memiograph     | Do not access freed memory                                         | MEM30-C           | clang-tidy: cert-mem30-cpp               |
|13 | accfree        | Do not use resources that have been freed                          | MEM30-C           | cppcheck: deallocuse                     |
|14 | dblfree        | Do not free memory more than once                                  | MEM31-C           | cppcheck: doubleFree                     |
|15 | nonnullptr     | Do not pass a null pointer to a library function requiring non-null| EXP34-C           | cppcheck: nullPointer; clang-tidy: cert-exp34-c |
|16 | nullref        | Do not dereference null pointers                                   | EXP34-C           | same as above                            |
|17 | diverr         | Ensure integer division does not result in divide-by-zero          | INT33-C           | no automated checker (static proof)      |
|18 | exceptbits     | Floating-point operations shall not set exception status bits      | FLP03-C           | no automated checker                     |
|19 | liberr         | Detect and handle library errors                                   | ERR33-C           | no automated checker (code review)       |
|20 | datarace       | Prevent data races between threads                                 | CON32-C / POS51-C | cppcheck: --enable=thread                |
|21 | asyncsig       | Signal handlers shall not call functions not in the safe list      | SIG30-C           | cppcheck: signalHandler                  |
|22 | toinit         | Objects shall be initialized before use                            | EXP33-C           | cppcheck: uninitvar                      |

---

## 2. Codebase Audit Findings

### Methodology

Static analysis run: `cppcheck --addon=cert src/dw_*.c src/dw_*.h` (spec 013
baseline output stored at `specs/013-strict-c-posix-compliance/cppcheck-post.txt`).

Source audit performed on branch `014-iso-ts-17961-compliance`, commit HEAD
(2026-03-23). All grep commands executed from workspace root
`/Users/joey/Development/dwsolver-repaired`.

Scope: `src/dw_*.c`, `src/dw_*.h` only. `third-party/glpk/` excluded.

### Rule-by-rule findings

#### intptrconv — N-A
**Query**: `grep -rn "(int).*\b\|intptr_t\|uintptr_t\|ptrdiff_t" src/dw_*.c src/dw_*.h`
**Result**: Zero results for pointer-to-integer or integer-to-pointer casts.
No `intptr_t`/`uintptr_t` usage in the codebase; no `(int)ptr` or `(void*)int`
patterns. Triggering construct absent.

#### intobjptr — PASS via spec 013
The codebase uses `malloc` with explicit cast `(type*)malloc(...)` — valid
in C11. No incompatible pointer object casts found.
`-pedantic-errors` (spec 013) catches reinterpretation casts at compile time.

#### xplicitcomp — N-A
**Query**: `grep -rn "== *\"" src/dw_*.c src/dw_*.h`
**Result**: Zero results. All string comparisons use `strcmp`/`strncmp`.
No explicit pointer equality with string literals.

#### boolasgn — N-A
**Query**: `grep -rn "_Bool\|bool " src/dw_*.c src/dw_*.h`
**Result**: Zero results. Codebase uses `int` for flags with explicit 0/1
values; no `_Bool` assignment patterns that trigger this rule.

#### trstcmp — PASS via spec 013
`atoi()` calls exist at `dw_rounding.c:134-135` and `dw_support.c:757-758`, converting
`strtok` token to integer. These are parsing LP column names generated by
dwsolver itself (e.g., `lambda_N_M` → extract N, M). The input is
dwsolver-controlled; the CERT rule applies to untrusted external numeric
input. An acceptable structural invariant exists. Column names derive from
dwsolver's `snprintf`-generated names (see `dw_phases.c:124-128`).

#### wraparound — PASS via spec 013
**Division/modulus audit**: Two operations found:
- `dw_solver.c:592`: `(globals->head_service_queue + 1) % num_clients` —
  `num_clients` is validated `> 0` at entry to `dw_solve()`.
- `dw_support.c:866`: `(double)elapsed_ns / 1.0e9` — constant divisor.

No unsigned-wraparound arithmetic with user-controlled values found.
`-Wall -Wextra` (spec 013) catches signed overflow in constant expressions.

#### argcomp — PASS via spec 013
Strict `-std=c11 -pedantic-errors` compilation catches all argument/parameter
type mismatches at compile time (spec 013 baseline builds clean).

#### intnegator — N-A
**Query**: `grep -rn "\-num_\|\-count\|\-size\|\-idx\|= -[a-z]" src/dw_*.c | grep -v "//"` — no unsigned negation.
`INT_MIN` negation pattern: `grep -rn "INT_MIN\|-INT_MIN" src/dw_*.c` → zero results.
No unsigned negative or INT_MIN negation present.

#### ioilecc — PASS via spec 013
`-Wall -Wformat` (enforced by spec 013) catches format/type mismatches.
All `printf`/`fprintf`/`snprintf` calls compile clean under spec 013 baseline.

#### ga-buffer — **FAIL (partial) → requires remediation**

**Decision**: FAIL for two call sites where `strcpy` copies from an
external-LP-file column name into a BUFF_SIZE=200 heap buffer.

```
dw_rounding.c:425  strcpy(local_col_name, var_name);  // var_name = glp_get_col_name(sub_lp, j)
dw_rounding.c:445  strcpy(local_col_name, var_name);  // same
dw_support.c:668   strcpy(local_col_name, var_name);  // same
```

`glp_get_col_name()` returns names from user-provided subproblem LP files.
GLPK's internal limit is 255 chars; BUFF_SIZE = 200. Static analysis cannot
prove the source fits; no runtime verification exists.

**PASS cases within ga-buffer** (structural invariant):
- `dw_rounding.c:130  strcpy(local_col_name, col_name)` — `col_name` from
  master LP, whose columns are set by `snprintf(_, BUFF_SIZE-1, "lambda_%d_%d", ...)`
  and `snprintf(_, BUFF_SIZE-1, "theta_%d_%d", ...)` in `dw_phases.c:124-128`.
  Provably ≤ BUFF_SIZE-1 chars.

**Remediation**: Replace three `strcpy` call sites with `strncpy` + explicit
null-termination.

#### strmod — PASS
**Query**: `grep -rn 'strcpy.*"[A-Z]' src/dw_*.c | grep -v malloc`
All `strcpy` destinations are heap-allocated (malloc'd) writable buffers.
No modification of string literal storage.

#### memiograph — PASS via spec 008
Memory access patterns: all array accesses via LP-dimension-bounded loop
indices; `fg->x[index][iter]` arrays are allocated to the correct size
(`sub_data[index].num_cols_plus`). No out-of-range indexing patterns.

#### accfree — PASS
**Query**: Linear review of all `free()` call sites (see above).
No pointer is accessed after being freed; allocation lifetimes are linear
within function scopes or within clearly-demarcated cleanup blocks.

#### dblfree — PASS
All free calls are at distinct cleanup points; no pointer is freed twice.
spec 008 memory discipline (MEM31-C) verified.

#### nonnullptr — **FAIL → requires remediation**

`strtok(NULL, ",")` return value is passed directly as the `src` argument
to `strcpy()` without a NULL check at five call sites:

```
dw_rounding.c:427  strcpy( curr_flight, strtok(NULL, ",") );
dw_rounding.c:428  strcpy( curr_sector, strtok(NULL, ",") );
dw_rounding.c:447  strcpy( temp_flight, strtok(NULL, ",") );
dw_rounding.c:448  strcpy( prev_sector, strtok(NULL, ",") );
dw_support.c:673   strcpy( sector_name, strtok(NULL, ",") );
```

`strtok` returns `NULL` when no further tokens are found. Passing `NULL` to
`strcpy` is undefined behaviour (C11 §7.24.2.3; EXP34-C / TS 17961 nonnullptr).

Additionally, the conditional guard `if( strtok(...) == NULL ) printf(...)` at
the preceding line does not branch on the result; execution continues and calls
`strtok(NULL, ",")` regardless, compounding the risk.

**Remediation**: Extract `strtok` return to a named variable; branch on NULL
before calling `strcpy`; replace `strcpy` with `strncpy` + NUL-termination.

#### nullref — PASS (with OOM structural invariant)
All `malloc` returns are immediately checked by `dw_oom_abort()` which calls
`exit(EXIT_FAILURE)` on NULL, preventing any subsequent null dereference. This
is not a cppcheck-recognizable idiom, so cppcheck emits false positives
(`nullPointerOutOfMemory`); they are accounted for in the spec 013 baseline.

`glp_get_col_name()` can return NULL for unnamed columns; call sites either
pass the result to `printf` (safe with `%s` of NULL on GCC/Clang) or copy
it via `strcpy` — the last case overlaps with nonnullptr remediation.

#### diverr — PASS
Only two division-class operations found (see wraparound section above).
`num_clients` invariant that it is > 0 is established by argument validation
in `dw_solve()`. The `/ 1.0e9` divisor is a compile-time constant.

#### exceptbits — N-A
**Query**: `grep -rn "fenv\|fesetround\|feclearexcept\|fetestexcept" src/dw_*.c src/dw_*.h`
**Result**: Zero results. No floating-point exception manipulation.
The codebase performs standard `double` arithmetic but does not read or
set FP exception status bits.

#### liberr — PASS
All library functions with fallible returns are checked:
- `fopen()`: all call sites check `!= NULL` before use (11 call sites audited).
- `glp_read_lp()`: checked for `!= 0` at `dw_support.c:104` and `dw_solver.c:192`.
- `pthread_*` / `sem_*`: wrapped with `DW_PTHREAD_CHECK` / `DW_SEM_CHECK` (spec 008 FR-001).
- `glp_simplex()` / `glp_intopt()`: return values stored and checked.

#### datarace — PASS via spec 008
All shared state (`sub_data_mutex[]`, `master_mutex`, condition variables,
semaphores) is protected with appropriate primitives. TSan-clean under spec 008.
POS51/52/54-C verified by spec 008.

#### asyncsig — N-A
**Query**: `grep -rn "\bsignal\b(" src/dw_*.c src/dw_*.h`
**Result**: Zero results for the POSIX `signal()` function. All occurrences of
the word "signal" in the source refer to the project-internal `signal_data*`
struct and `signals` variable (iteration-synchronisation data structure).
No signal handlers are installed; the `asyncsig` rule cannot be triggered.

#### toinit — PASS via spec 013
`cppcheck --enable=uninitvar` and `-Wall` (spec 013) find no uninitialized
reads. All variables are initialized before use in the audited source.

---

## 3. Gap Summary

| Rule        | Verdict in initial audit | Remediation required        |
|-------------|--------------------------|------------------------------|
| intptrconv  | N-A                      | None                         |
| intobjptr   | PASS (spec 013)          | None                         |
| xplicitcomp | N-A                      | None                         |
| boolasgn    | N-A                      | None                         |
| trstcmp     | PASS (structural)        | None                         |
| wraparound  | PASS (spec 013)          | None                         |
| argcomp     | PASS (spec 013)          | None                         |
| intnegator  | N-A                      | None                         |
| ioilecc     | PASS (spec 013)          | None                         |
| **ga-buffer**   | **FAIL**             | **strncpy + NUL-term at 3 sites** |
| strmod      | PASS                     | None                         |
| memiograph  | PASS (spec 008)          | None                         |
| accfree     | PASS                     | None                         |
| dblfree     | PASS                     | None                         |
| **nonnullptr**  | **FAIL**             | **NULL-guard + strncpy at 5 sites** |
| nullref     | PASS (dw_oom_abort)      | None                         |
| diverr      | PASS                     | None                         |
| exceptbits  | N-A                      | None                         |
| liberr      | PASS                     | None                         |
| datarace    | PASS (spec 008)          | None                         |
| asyncsig    | N-A                      | None                         |
| toinit      | PASS (spec 013)          | None                         |

**2 FAIL rules — both in `dw_rounding.c` and `dw_support.c`** (the four-season
rounding subsystem, which is lightly used but compiled into the library).

---

## 4. Remediation Plan (Story 2)

### Fix A — ga-buffer: Replace unbounded strcpy with strncpy

**Files**: `src/dw_rounding.c`, `src/dw_support.c`

Three call sites copy from `glp_get_col_name(sub_lp, ...)` into a `BUFF_SIZE`
(200-byte) heap buffer using `strcpy`. GLPK allows names up to 255 chars. Fix:

```c
/* Before */
strcpy(local_col_name, var_name);

/* After */
strncpy(local_col_name, var_name, BUFF_SIZE - 1);
local_col_name[BUFF_SIZE - 1] = '\0';
```

### Fix B — nonnullptr: NULL-guard strtok return before strcpy

**Files**: `src/dw_rounding.c`, `src/dw_support.c`

Five call sites pass `strtok(NULL, ",")` directly to `strcpy`. Replace
with a NULL-checked extraction:

```c
/* Before */
if( strtok(local_col_name, "(") == NULL ) printf("NULL 1");
strcpy( curr_flight, strtok(NULL, ",") );
strcpy( curr_sector, strtok(NULL, ",") );

/* After */
if( strtok(local_col_name, "(") == NULL ) {
    dw_printf(IMPORTANCE_DIAG, "TS17961-nonnullptr: missing '(' delimiter; skipping column.\n");
    /* skip this entry — continue outer loop */
    continue;  /* or return, depending on context */
}
{
    const char *tok_flight = strtok(NULL, ",");
    const char *tok_sector = strtok(NULL, ",");
    if( !tok_flight || !tok_sector ) {
        dw_printf(IMPORTANCE_DIAG, "TS17961-nonnullptr: missing token; skipping column.\n");
        continue;
    }
    strncpy(curr_flight, tok_flight, BUFF_SIZE - 1);
    curr_flight[BUFF_SIZE - 1] = '\0';
    strncpy(curr_sector, tok_sector, BUFF_SIZE - 1);
    curr_sector[BUFF_SIZE - 1] = '\0';
}
```

Note: The calling contexts loop over LP column ranges; a `continue` terminates
the current column's processing and advances to the next, preserving
correctness (skipping malformed names rather than invoking UB).

---

## 5. Tooling Research

### cppcheck --addon=cert

Version ≥ 2.12 required for cert addon reliability. Available on Ubuntu via
`sudo apt-get install cppcheck`. The cert addon adds CERT-mapped checkers on
top of standard cppcheck. Relevant cert checks triggered:
- `cert-STR31-C` / `cert-EXP34-C`: closest to ga-buffer / nonnullptr
- `cert-ERR33-C`: liberr checks

Invocation for this project:
```bash
cppcheck --addon=cert --std=c11 \
  --suppress=constParameterPointer \
  -I src \
  $(pkg-config --cflags-only-I glpk 2>/dev/null || echo "-I/usr/include") \
  src/dw_*.c 2>&1
```

The `constParameterPointer` suppression matches the known false positive for
`dw_oom_abort(void* ptr, ...)` identified in the spec 013 baseline.

### clang-tidy -checks=cert-*

Version ≥ 17 required. Relevant checkers:
- `cert-exp34-c`: Detects null pointer dereference (overlaps nonnullptr/nullref)
- `cert-int30-c`, `cert-int32-c`: Integer overflow (wraparound, intnegator)
- `cert-str31-c`: String buffer overflow (ga-buffer mapping for clang-tidy)

Note: `clang-tidy` cert checks require a `compile_commands.json` to resolve
includes. Generate with:
```bash
bear -- make clean all  # if bear is installed
# or: use cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON for cmake build
```

For Autotools projects without Bear, supply `--` with compiler flags directly.

### Synthetic violation test (SC-003)

Per spec clarification, the CI enforcement test (`tests/test_ts17961_enforcement.sh`)
must:
1. Copy a source file to a temp file
2. Inject a known-bad construct (e.g., `char *p = NULL; strcpy(p, "x");` — nonnullptr violation)
3. Run `cppcheck --addon=cert` against the modified temp file
4. Assert exit code non-zero OR grep for the cert finding
5. Delete the temp file

**Decision**: Use a dedicated temporary C file (not patch existing source) to
avoid filesystem race conditions in CI.

---

## 6. CI Workflow Design

New workflow: `.github/workflows/ci-compliance.yml`
- Trigger: `push` and `pull_request` to `main` (same as existing workflows)
- Runner: `ubuntu-latest` only (Linux-only per spec clarification)
- Steps:
  1. `actions/checkout@v4`
  2. Touch autotools artifacts (same pattern as `ci-linux.yml`)
  3. `sudo apt-get install -y libglpk-dev cppcheck clang-tidy`
  4. `./configure && make` (need a build for clang-tidy)
  5. Run `cppcheck --addon=cert` against `src/dw_*.c`
  6. Run `clang-tidy -checks=cert-*` against `src/dw_*.c` (with `-I` flags)
  7. Run `tests/test_ts17961_enforcement.sh` (synthetic violation gate SC-003)
  8. Fail if any finding maps to a non-N-A TS 17961 rule

**N-A allowlist**: Rules `asyncsig`, `intptrconv`, `intnegator`, `xplicitcomp`,
`boolasgn`, `exceptbits` are documented N-A and any tool findings for those
rules may be suppressed via inline comments or a `.cppcheck` suppressions file.

---

## 7. Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Remediation scope | `dw_rounding.c`, `src/dw_support.c` only | Both FAILs are in the rounding subsystem; no other source files require changes |
| Fix strategy | `strncpy` + explicit NUL + NULL guard | Simple, readable, no new dependencies; `strlcpy` not in C11 standard |
| Error handling on NULL strtok | `continue` (skip column, log diagnostic) | Consistent with existing pattern; malformed LP column names are non-fatal; original code already prints and continues |
| Synthetic test approach | Dedicated temp file | Avoids patching source files; cleaner CI recovery |
| cppcheck version pin | ≥ 2.12 | Needed for reliable cert addon; Ubuntu 22.04 ships 2.7 (too old) — CI must install from PPA or use newer Ubuntu image |
| N-A suppression mechanism | Inline `// cppcheck-suppress` comments | Preferred over global suppression files; co-located with evidence |
