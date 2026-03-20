# Quickstart: Code Quality Improvement Workflow

**Branch**: `chore/code-quality` (create from master after `005-windows-build-fix` merges)
**Spec**: `specs/005-windows-build-fix/spec.md`

This guide walks a developer through building, verifying, and iterating on the
code quality changes described in the spec.

---

## Prerequisites

| Tool | Minimum version | Install hint |
|------|-----------------|--------------|
| GCC or Clang | GCC ≥ 9 / Clang ≥ 12 | `brew install gcc` / `apt install gcc` |
| GNU Autotools | autoconf ≥ 2.69, automake ≥ 1.16 | `brew install autoconf automake` |
| GLPK headers | 4.44 (bundled) | already in `src/` |
| pthreads | POSIX or winpthreads | system-provided |
| bash | ≥ 3.2 | system-provided |

---

## 1. Initial build (fresh checkout)

```bash
# From repo root
autoreconf -fiv
./configure
make -j$(nproc || sysctl -n hw.ncpu)
```

Expected: no warnings on `src/dw*.c` files (there are existing warnings; note
their current count before making changes so regressions can be spotted).

Capture baseline warning count:

```bash
make clean && make 2>&1 | grep -c 'warning:'
```

Record that number. After your changes the warning count should be equal or
lower — never higher.

---

## 2. Run the test suite

```bash
cd tests
bash dw-tests.sh
```

Expected: all 13 tests pass. The suite runs each example in `examples/` and
compares solver output.

If a test hangs, `Ctrl-C` and inspect with:

```bash
bash dw-tests.sh 2>&1 | head -40
```

---

## 3. Applying the P1 fixes

### Fix `val = NULL` (dw_phases.c)

1. Locate the diagnostic block in `phase_2_iteration` — search for
   `glp_get_mat_col` inside that function.
2. Follow the before/after pattern in `data-model.md`.
3. Recompile and run tests.
4. Confirm no `valgrind` errors (optional but recommended):
   ```bash
   valgrind --error-exitcode=1 ./dwsolver examples/example1.lp
   ```

### Fix `const char*` leak (dw_rounding.c)

1. Locate `local_col_name` declaration in `rounding_thread`.
2. Replace `malloc` with a stack buffer as shown in `data-model.md`.
3. Remove the matching `free(local_col_name)` call (if one exists) below.
4. Recompile and run tests.
5. Valgrind will report the before-fix leak under
   `dw_rounding.c: rounding_thread` — confirm it is gone after.

---

## 4. Applying the P2 refactor (`dw_iteration` unification)

1. Create the new `dw_iteration` function signature in `dw_phases.h` (see
   `data-model.md`).
2. Implement the unified body in `dw_phases.c`. Diff both original functions
   side-by-side to verify every line is covered.
3. Update `dw_main.c` call sites.
4. Leave `phase_1_iteration` / `phase_2_iteration` as thin wrappers calling
   `dw_iteration` until the PR is merged — this makes the diff reviewable
   without a massive single-commit rewrite.
5. Run tests after each incremental change.

---

## 5. Applying P3 changes

### Verbosity guard around `glpk_mutex` (dw_subprob.c)

Locate `prepare_column` and wrap the `pthread_mutex_lock(&glpk_mutex)` block
in a verbosity check:

```c
if (globals->verbose >= OUTPUT_ALL) {
    pthread_mutex_lock(&glpk_mutex);
    /* ... logging ... */
    pthread_mutex_unlock(&glpk_mutex);
}
```

Verify no output is lost when verbosity is set to `OUTPUT_ALL`:
```bash
./dwsolver --verbose=3 examples/example1.lp   # adjust flag as needed
```

### Replace `master_mutex` with `pthread_rwlock_t`

1. Add `pthread_rwlock_t master_data_rwlock` to `dw.h` externs and
   `dw_globals.c` definitions.
2. In `subproblem_thread`, replace `pthread_mutex_lock(&master_mutex)` around
   the column-mapping read loop with `pthread_rwlock_rdlock».
3. Ensure `init_pthread_data()` initialises the rwlock and `free_globals()`
   destroys it.
4. Run tests with `TSAN` if available:
   ```bash
   make clean && CFLAGS="-fsanitize=thread -g" LDFLAGS="-fsanitize=thread" \
     ./configure && make && bash tests/dw-tests.sh
   ```

---

## 6. Verification checklist before opening a PR

- [ ] `make` completes with no new warnings
- [ ] `bash tests/dw-tests.sh` — all 13 tests pass
- [ ] `valgrind` reports no leaks in `rounding_thread` (P1 fix)
- [ ] `valgrind` reports no invalid read/write in `phase_2_iteration` (P1 fix)
- [ ] Each P3 lock change verified by running with both `verbose=0` and max verbosity
- [ ] `git diff --stat` shows no unintended files changed

---

## 7. Creating the branch

After `005-windows-build-fix` is merged to master:

```bash
git checkout master
git pull origin master
git checkout -b chore/code-quality
```

Commit each P1 fix separately (smallest diff = easiest review):

```bash
git add src/dw_phases.c
git commit -m "fix: allocate val/ind2 in phase_2_iteration diagnostic block"

git add src/dw_rounding.c
git commit -m "fix: use stack buffer for local_col_name in rounding_thread"
```

Then P2:
```bash
git add src/dw_phases.c src/dw_phases.h
git commit -m "refactor: unify phase_1/phase_2 iteration into dw_iteration"
```

Then P3 in one or two commits.
