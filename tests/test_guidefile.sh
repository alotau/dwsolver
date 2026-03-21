#! /bin/sh

# test_guidefile.sh: CLI-level tests for dwsolver guidefile parsing.
# Run from the tests/ directory with dwsolver on PATH.
# Usage: cd tests && bash test_guidefile.sh

PASS_COUNT=0
FAIL_COUNT=0

assert_exit() {
    local name="$1"
    local expected_check="$2"   # "zero", "nonzero", or "no_crash"
    local actual_rc="$3"
    local stderr_text="$4"

    case "$expected_check" in
        zero)
            if [ "$actual_rc" -eq 0 ]; then
                echo "  PASS: $name (exit $actual_rc)"
                PASS_COUNT=$((PASS_COUNT + 1))
            else
                echo "  FAIL: $name (expected exit 0, got $actual_rc)"
                FAIL_COUNT=$((FAIL_COUNT + 1))
            fi
            ;;
        nonzero)
            if [ "$actual_rc" -ne 0 ]; then
                echo "  PASS: $name (exit $actual_rc, non-zero as expected)"
                PASS_COUNT=$((PASS_COUNT + 1))
            else
                echo "  FAIL: $name (expected non-zero exit, got 0)"
                FAIL_COUNT=$((FAIL_COUNT + 1))
            fi
            ;;
        no_crash)
            if [ "$actual_rc" -lt 128 ]; then
                echo "  PASS: $name (exit $actual_rc, no crash)"
                PASS_COUNT=$((PASS_COUNT + 1))
            else
                echo "  FAIL: $name (crash/signal, exit $actual_rc)"
                FAIL_COUNT=$((FAIL_COUNT + 1))
            fi
            ;;
    esac

    # Check for unexpected USAGE output in stderr
    if [ -n "$5" ] && echo "$stderr_text" | grep -q "USAGE:"; then
        echo "  FAIL: $name (unexpected USAGE text in output)"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
}

echo "=== test_guidefile.sh ==="
echo ""

# T027: valid 1-subproblem guidefile (single_sub)
echo "Test: valid 1-subproblem guidefile"
pushd ../examples/single_sub > /dev/null
_out=$(dwsolver --no-write-final-master --quiet -g guidefile 2>&1)
_rc=$?
popd > /dev/null
assert_exit "valid single_sub" zero "$_rc"
if echo "$_out" | grep -q "USAGE:"; then
    echo "  FAIL: valid single_sub printed USAGE text"
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# T028: valid 4-subproblem guidefile (four_sea)
echo "Test: valid 4-subproblem guidefile"
pushd ../examples/four_sea > /dev/null
_out=$(dwsolver --quiet -g guidefile 2>&1)
_rc=$?
popd > /dev/null
assert_exit "valid four_sea (no crash)" no_crash "$_rc"
if echo "$_out" | grep -qi "parse error\|USAGE:\|invalid guidefile"; then
    echo "  FAIL: four_sea printed unexpected parse-error text"
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# T029: bad_count guidefile (n=3 declared, only 2 files listed)
echo "Test: bad_count guidefile (too few filenames)"
_out=$(dwsolver -g fixtures/bad_count.guidefile 2>&1)
_rc=$?
assert_exit "bad_count non-zero exit" nonzero "$_rc"
if echo "$_out" | grep -q ".\|Problem"; then
    : # any output is acceptable (just checking non-zero exit)
fi
# Verify something was printed to stderr/stdout (not silent failure)
if [ -z "$_out" ]; then
    echo "  FAIL: bad_count produced no output (expected an error message)"
    FAIL_COUNT=$((FAIL_COUNT + 1))
else
    echo "  PASS: bad_count produced error output"
    PASS_COUNT=$((PASS_COUNT + 1))
fi

# T030: missing_file guidefile (references nonexistent file)
echo "Test: missing_file guidefile (file not found)"
_out=$(dwsolver -g fixtures/missing_file.guidefile 2>&1)
_rc=$?
assert_exit "missing_file no crash" no_crash "$_rc"
assert_exit "missing_file non-zero exit" nonzero "$_rc"

echo ""
echo "=== Results: $PASS_COUNT passed, $FAIL_COUNT failed ==="
if [ "$FAIL_COUNT" -ne 0 ]; then
    exit 1
fi
exit 0
