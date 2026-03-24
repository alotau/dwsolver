# Quickstart: ISO/IEC TS 17961:2013 Compliance

**Branch**: `014-iso-ts-17961-compliance`
**For**: Developers and auditors working on or reviewing this feature

---

## What is this feature?

This feature audits dwsolver's authored source files against all 22 analyzable
rules in ISO/IEC TS 17961:2013 ("C Secure Coding Rules"), remediates the two
gap rules found (`ga-buffer` and `nonnullptr` in the rounding subsystem), and
produces a formal compliance matrix and tooling guide suitable for procurement.

---

## Quick compliance check

```bash
# From the repository root — Linux with cppcheck ≥ 2.12 installed:
cppcheck --addon=cert --std=c11 \
  --suppress=constParameterPointer \
  -I src -I /usr/include \
  src/dw_blas.c src/dw_globals.c src/dw_main.c src/dw_phases.c \
  src/dw_rounding.c src/dw_solver.c src/dw_subprob.c src/dw_support.c

# Expected: zero cert-* findings (after remediation)
```

---

## Understanding the audit artifacts

```text
specs/014-iso-ts-17961-compliance/
├── spec.md                            Feature specification
├── research.md                        Phase 0: per-rule audit findings
├── data-model.md                      Phase 1: entities, modified files, VR
├── quickstart.md                      This file
├── audit/
│   ├── ts17961-compliance-matrix.md   Primary compliance artifact (22 rules)
│   └── tooling-guide.md               Reproducible tool invocations
└── acceptance-report.md               FR sign-off (produced by /speckit.tasks)
```

---

## Gap rules and remediations

Two rules failed initial audit — both in the four-season rounding subsystem:

| Rule        | Files                              | Fix                                         |
|-------------|------------------------------------|---------------------------------------------|
| `ga-buffer` | `src/dw_rounding.c`, `src/dw_support.c` | Replace `strcpy` with `strncpy` + NUL-term |
| `nonnullptr`| `src/dw_rounding.c`, `src/dw_support.c` | NULL-guard `strtok` return before `strcpy` |

See `research.md §4` for exact before/after code patterns.

---

## Running the CI compliance job locally

```bash
# Install tools (Ubuntu):
sudo apt-get install -y cppcheck clang-tidy

# cppcheck cert check:
cppcheck --addon=cert --std=c11 --suppress=constParameterPointer \
  -I src src/dw_*.c

# clang-tidy cert checks (requires compile_commands.json or explicit flags):
clang-tidy -checks='cert-*' \
  --extra-arg=-std=c11 \
  --extra-arg=-I./src \
  src/dw_*.c -- \
  -I/usr/include -DHAVE_CONFIG_H

# Synthetic violation gate (SC-003):
bash tests/test_ts17961_enforcement.sh
```

---

## Key TS 17961 → CERT C rule mappings

| TS 17961 rule | Nearest CERT C | cppcheck checker         |
|---------------|----------------|--------------------------|
| ga-buffer     | STR31-C        | bufferAccessOutOfBounds  |
| nonnullptr    | EXP34-C        | nullPointer              |
| nullref       | EXP34-C        | nullPointer              |
| datarace      | CON32-C        | (TSan; spec 008)         |
| asyncsig      | SIG30-C        | signalHandler            |

---

## Prior work references

- **Spec 008** (`specs/008-sei-cert-c-compliance/acceptance-report.md`):
  threading safety (POS51/52/54-C), OOM guards (MEM32-C), file handle
  errors (FIO01/42-C), integer types (INT30/31-C). Covers: `datarace`,
  `liberr` (pthread), `memiograph`.

- **Spec 013** (`specs/013-strict-c-posix-compliance/acceptance-report.md`):
  strict `-std=c11 -pedantic-errors -Wall -Wextra`, `_POSIX_C_SOURCE=200809L`.
  Covers: `argcomp`, `ioilecc`, `wraparound`, `intobjptr`, `toinit`.
