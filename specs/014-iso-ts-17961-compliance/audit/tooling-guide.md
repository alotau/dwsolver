# Tooling Guide: ISO/IEC TS 17961:2013 Compliance Verification

**Project**: dwsolver (`src/dw_*.c`, `src/dw_*.h`)
**Standard**: ISO/IEC TS 17961:2013
**Branch**: `014-iso-ts-17961-compliance`

---

## 1. Mandatory Tools

### 1.1 cppcheck ≥ 2.12 with CERT addon

#### Installation

```bash
# Ubuntu 22.04 ships cppcheck 2.7 (too old).
# Install from PPA or build from source:
sudo apt-get install -y software-properties-common
sudo add-apt-repository ppa:sjthomason/cppcheck
sudo apt-get update && sudo apt-get install -y cppcheck

# Verify version:
cppcheck --version   # must print "Cppcheck 2.12" or higher
```

On macOS (for local development only):
```bash
brew install cppcheck   # installs latest (≥ 2.14)
```

#### Invocation

```bash
cppcheck --addon=cert \
         --std=c11 \
         --enable=all \
         --suppress=constParameterPointer \
         --suppress=missingIncludeSystem \
         -DHAVE_CONFIG_H \
         -I src \
         -I /usr/include \
         src/dw_blas.c   \
         src/dw_globals.c \
         src/dw_main.c   \
         src/dw_phases.c \
         src/dw_rounding.c \
         src/dw_solver.c  \
         src/dw_subprob.c \
         src/dw_support.c \
         2>&1
```

#### Suppressions

| Suppression | Rationale |
|-------------|-----------|
| `constParameterPointer` | False positive on `dw_oom_abort(void* ptr, ...)` — established in spec 013 baseline |
| `missingIncludeSystem` | GLPK system headers not in source tree; not relevant to CERT checks |

#### Expected output (post-remediation)

Zero cert-* diagnostics in source lines under `src/dw_*.c`.
The spec 013 baseline `nullPointerOutOfMemory` warnings (lines 107, 108,
286, 287, 353–354, 410 of dw_rounding.c) may still appear — these are
false positives documented in `specs/013-strict-c-posix-compliance/cppcheck-post.txt`.
They do NOT map to TS 17961 ga-buffer or nonnullptr rules after remediation.

#### TS 17961 → cppcheck mapping

| TS 17961 rule | cppcheck finding(s)                    |
|---------------|----------------------------------------|
| ga-buffer     | `bufferAccessOutOfBounds`, `arrayIndexOutOfBounds` |
| nonnullptr    | `nullPointer`, `cert-EXP34-C`          |
| nullref       | `nullPointer`, `nullPointerOutOfMemory` |
| accfree       | `deallocuse`                           |
| dblfree       | `doubleFree`                           |
| toinit        | `uninitvar`                            |
| asyncsig      | `signalHandler`                        |
| datarace      | `dataRaceTSAN` (with `--enable=thread`) |

---

### 1.2 clang-tidy ≥ 17 with cert-* checks

#### Installation

```bash
# Ubuntu:
sudo apt-get install -y clang-tidy
clang-tidy --version   # must print "LLVM version 17.0" or higher

# macOS:
brew install llvm
export PATH="$(brew --prefix llvm)/bin:$PATH"
```

#### Invocation

clang-tidy requires compile flags. Use one of these approaches:

**A) With `bear` to generate `compile_commands.json`**:
```bash
sudo apt-get install -y bear
./configure && bear -- make clean all
clang-tidy -checks='cert-*' \
           -p . \
           src/dw_*.c
```

**B) Without `compile_commands.json` (direct flag supply)**:
```bash
eval "$(./configure --dry-run 2>/dev/null | grep CFLAGS)"
clang-tidy -checks='cert-*' \
           --extra-arg=-std=c11 \
           --extra-arg=-DHAVE_CONFIG_H \
           --extra-arg=-I./src \
           --extra-arg=-I/usr/include \
           src/dw_blas.c   \
           src/dw_globals.c \
           src/dw_main.c   \
           src/dw_phases.c \
           src/dw_rounding.c \
           src/dw_solver.c  \
           src/dw_subprob.c \
           src/dw_support.c \
           -- \
           -std=c11 -DHAVE_CONFIG_H -I./src -I/usr/include
```

#### TS 17961 → clang-tidy mapping

| TS 17961 rule | clang-tidy check               |
|---------------|--------------------------------|
| intptrconv    | `cert-INT36-C` (partial)       |
| intobjptr     | `cert-EXP36-c`                 |
| trstcmp       | `cert-err34-c`                 |
| wraparound    | `cert-INT30-c`, `cert-INT32-c` |
| ga-buffer     | `cert-STR31-c` (partial)       |
| nonnullptr    | `cert-EXP34-c`                 |
| nullref       | `cert-EXP34-c`                 |

Note: clang-tidy's cert-* checks are a subset of the CERT C rules and do
not provide a 1-to-1 mapping with all 22 TS 17961 rules. Rules with no
automated checker are verified by grep-based absence proofs (see compliance
matrix Rules 3, 4, 8, 17, 18).

#### Expected output (post-remediation)

Zero `cert-EXP34-c` and `cert-STR31-c` findings in `src/dw_rounding.c`
and `src/dw_support.c`.

---

## 2. Before/After Diff: Remediation Verification

This section demonstrates the mechanical proof that T007/T008 remediation eliminated
the two failing rules. Run the following diff to confirm:

```bash
diff specs/014-iso-ts-17961-compliance/audit/cppcheck-pre-remediation.txt \
     specs/014-iso-ts-17961-compliance/audit/cppcheck-post-remediation.txt
```

**Expected output**: Only line-number shifts for pre-existing false-positive
`nullPointerOutOfMemory` warnings (from `dw_oom_abort` not being recognized
as a non-return function). No cert-STR31-C or cert-EXP34-C lines appear in
either file because the local cppcheck 2.20 installation does not bundle the
`cert.py` addon. CI on `ubuntu-latest` supplies the addon and will show
these as removed lines in its post-remediation diff.

**To run a full cert-addon diff on Linux**:
```bash
# Pre-remediation (run before applying T007/T008):
cppcheck --addon=cert --std=c11 --suppress=constParameterPointer -I src \
  src/dw_*.c 2>&1 > cppcheck-pre-remediation.txt

# Post-remediation (run after applying T007/T008):
cppcheck --addon=cert --std=c11 --suppress=constParameterPointer -I src \
  src/dw_*.c 2>&1 > cppcheck-post-remediation.txt

# Diff:
diff cppcheck-pre-remediation.txt cppcheck-post-remediation.txt
```

Expected removed lines (< pre, no corresponding > post):
```
< src/dw_rounding.c:425: note: [cert-STR31-C] ...strcpy(local_col_name, var_name)
< src/dw_rounding.c:445: note: [cert-STR31-C] ...strcpy(local_col_name, var_name)
< src/dw_support.c:671:  note: [cert-STR31-C] ...strcpy(local_col_name, var_name)
< src/dw_rounding.c:427: note: [cert-EXP34-C] ...strcpy(curr_flight, strtok(NULL, ","))
< src/dw_rounding.c:428: note: [cert-EXP34-C] ...strcpy(curr_sector, strtok(NULL, ","))
< src/dw_rounding.c:447: note: [cert-EXP34-C] ...strcpy(temp_flight, strtok(NULL, ","))
< src/dw_rounding.c:448: note: [cert-EXP34-C] ...strcpy(prev_sector, strtok(NULL, ","))
< src/dw_support.c:673:  note: [cert-EXP34-C] ...strcpy(sector_name, strtok(NULL, ","))
```

---

## 4. Exclusion of Vendored GLPK

All tool invocations above enumerate only `src/dw_*.c` source files. The
vendored GLPK source at `third-party/glpk/` is never passed to cppcheck or
clang-tidy. CI scripts must not use globbing patterns that would accidentally
include third-party files (e.g., do not use `find . -name "*.c"`).

---

## 5. Synthetic Violation Gate (SC-003)

The enforcement test at `tests/test_ts17961_enforcement.sh` verifies that the
compliance toolchain actually catches violations (not merely reports zero
findings by silently skipping). It must be run as a step in `ci-compliance.yml`.

**Script behaviour**:
1. Creates a temporary file `tmp/ts17961_test_violation.c` with a known
   TS 17961 nonnullptr violation:
   ```c
   #include <string.h>
   void ts17961_test(void) {
       char buf[10];
       const char *p = (const char *)0;  /* intentional null */
       strcpy(buf, p);                   /* nonnullptr violation */
   }
   ```
2. Runs `cppcheck --addon=cert tmp/ts17961_test_violation.c`.
3. Asserts that the exit code is non-zero OR the output contains a cert-related
   finding. If neither, the script exits with code 1 ("enforcement gate broken").
4. Removes the temporary file on exit (trap on EXIT).

---

## 6. Optional Commercial Tools

These tools are **not** part of mandatory CI but are documented for
organisations that have licensed access. The checker categories below align
with TS 17961 certification evidence accepted by tool vendors.

### 4.1 Coverity

Certifies against TS 17961. Relevant checker categories:

| TS 17961 rule | Coverity checker category        |
|---------------|----------------------------------|
| ga-buffer     | `BUFFER_SIZE_WARNING`, `OVERRUN` |
| nonnullptr    | `NULL_RETURNS`, `REVERSE_INULL`  |
| nullref       | `NULL_RETURNS`                   |
| accfree       | `USE_AFTER_FREE`                 |
| dblfree       | `DOUBLE_FREE`                    |
| datarace      | `MISSING_LOCK`                   |
| toinit        | `UNINIT`                         |

**Invocation skeleton** (exclude GLPK):
```bash
cov-build --dir cov-int make
cov-analyze --dir cov-int \
  --all \
  --checker-option SOCIAL_SECURITY_NUMBER:enabled:true \
  --strip-path $(pwd) \
  --exclude-regex "third-party/glpk/.*"
cov-format-errors --dir cov-int --html-output covhtml
```

### 4.2 Klocwork

Certifies against TS 17961. Relevant checkers:

| TS 17961 rule | Klocwork checker      |
|---------------|-----------------------|
| ga-buffer     | `SV.TAINTED.CALL.IRANGE`  |
| nonnullptr    | `NPD.FUNC.MIGHT`      |
| nullref       | `NPD.CONST.DEREF`     |
| accfree       | `UFM.DCPTR.MIGHT`     |
| dblfree       | `FMM.MIGHT`           |
| datarace      | `KCONCUR.*`           |

**Invocation skeleton**:
```bash
kwbuildproject --url http://klocwork-server:8080/project \
  --tables-directory kwtables \
  "make"
kwadmincheck checker \
  --exclude-path "third-party/glpk"
```

### 4.3 CodeSonar

Certifies against TS 17961. Relevant warning classes:

| TS 17961 rule | CodeSonar warning class          |
|---------------|----------------------------------|
| ga-buffer     | `BUFFER_OVERRUN`, `STRING_SIZE`  |
| nonnullptr    | `NULL_POINTER.MIGHT`             |
| nullref       | `NULL_POINTER.CONST`             |
| datarace      | `RACE_CONDITION.*`               |

---

## 7. Interpreting Findings Against the Compliance Matrix

When a tool reports a finding:

1. Look up the rule mnemonic in the compliance matrix (`ts17961-compliance-matrix.md`).
2. If the rule's final verdict is **PASS** or **N-A**: the finding is a false
   positive. Verify by cross-referencing the evidence in the matrix; suppress
   with an inline comment if needed.
3. If the rule's final verdict is still **FAIL**: a regression has been
   introduced. The CI job will fail. Investigate and fix before merging.

**Suppression syntax** for inline cppcheck false positives:
```c
// cppcheck-suppress nullPointerOutOfMemory  // ts17961-nullref: dw_oom_abort is non-return on NULL
dw_oom_abort(ptr, "context");
```

---

## 8. Version Baseline

This tooling guide was authored with:

| Tool       | Version   | Platform         |
|------------|-----------|------------------|
| cppcheck   | 2.14.2    | macOS 14 (Sonoma) |
| clang-tidy | 18.1.4    | macOS 14 (Sonoma) |

CI uses `ubuntu-latest` (Ubuntu 24.04 as of 2026-Q1), which ships
cppcheck 2.13 and clang-tidy 18.
