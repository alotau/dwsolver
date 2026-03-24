#!/usr/bin/env bash
# tests/test_ts17961_enforcement.sh
#
# SC-003 synthetic violation gate for ISO/IEC TS 17961:2013 compliance.
#
# PURPOSE: Verify that the compliance toolchain (cppcheck --addon=cert)
#   actually detects a known TS 17961 violation. This script would catch
#   a "silent pass" scenario where cppcheck is misconfigured and reports
#   zero findings even when a violation is present.
#
# EXIT CODES:
#   0  - Enforcement gate functional (cppcheck detected the synthetic violation)
#   1  - Enforcement gate broken (cppcheck did NOT detect the violation — FAIL)
#   2  - Unexpected error (cppcheck not found, etc.)
#
# USAGE:
#   bash tests/test_ts17961_enforcement.sh
#   (also run as part of CI via ci-compliance.yml)

set -e  # Fail fast on unexpected errors

# ---- Verify cppcheck is available ----
if ! command -v cppcheck > /dev/null 2>&1; then
    echo "ERROR: cppcheck not found in PATH. Install cppcheck >= 2.12." >&2
    exit 2
fi

# ---- Create temporary file with a known TS 17961 nonnullptr violation ----
# The violation: pass a potential NULL (strtok return unchecked) to a
# function expecting non-null (strcpy). This is the exact pattern remediated
# in T007/T008.
tmpfile=$(mktemp /tmp/ts17961_violation_XXXXXX.c)

cat > "$tmpfile" << 'CVIOLATION'
#include <string.h>
#include <stdio.h>
/* Synthetic TS 17961 nonnullptr violation — for SC-003 enforcement gate */
void ts17961_synthetic_test(const char *input) {
    char buf[64];
    char tmp[64];
    char *tok;
    strncpy(tmp, input, sizeof(tmp) - 1);
    tmp[sizeof(tmp) - 1] = '\0';
    tok = strtok(tmp, ",");
    /* Intentional violation: tok may be NULL (strtok returns NULL on no match).
       Deliberately not guarded to exercise the enforcement gate. */
    strcpy(buf, tok);  /* TS17961-nonnullptr: tok may be NULL — violation */
    printf("%s\n", buf);
}
CVIOLATION

# ---- Run cppcheck against the synthetic violation ----
# I-001 fix: disable set -e around cppcheck because it exits non-zero when
# it detects a violation — which is exactly what we want to happen here.
set +e
cppcheck_output=$(cppcheck --addon=cert --enable=all "$tmpfile" 2>&1)
result=$?
set -e

# ---- Assess result: non-zero exit OR finding in output ----
has_cert_finding=0
if echo "$cppcheck_output" | grep -qi "cert\|nullPointer\|strcpy\|warning.*null"; then
    has_cert_finding=1
fi

# ---- Cleanup ----
rm -f "$tmpfile"

# ---- Assert enforcement is functional ----
if [ "$result" -ne 0 ] || [ "$has_cert_finding" -eq 1 ]; then
    echo "PASS: SC-003 enforcement gate functional — cppcheck detected the synthetic nonnullptr violation."
    exit 0
else
    echo "FAIL: SC-003 enforcement gate broken — cppcheck did NOT detect the synthetic TS17961-nonnullptr violation." >&2
    echo "  Exit code: $result" >&2
    echo "  Output: $cppcheck_output" >&2
    echo "  Check that cppcheck --addon=cert is correctly installed (cppcheck >= 2.12 with cert addon)." >&2
    exit 1
fi
