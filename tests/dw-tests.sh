#! /bin/sh

# Quiet down pushd and popd
pushd () {
    command pushd "$@" > /dev/null
}
popd () {
    command popd "$@" > /dev/null
}

echo "Test: Bertsimas textbook example"
pushd ../examples/book_bertsimas
dwsolver --no-write-final-master --quiet -g guidefile
sort relaxed_solution > rs_sorted
if diff -w -B ex_relaxed_solution rs_sorted ; then
   echo "  PASS: Got expected solution."
   echo ""
else
   echo "  FAIL: Unexpected solution. Test failed. Exiting."
   exit 1
fi
popd

echo "Test: Lasdon textbook example"
pushd ../examples/book_lasdon
dwsolver --no-write-final-master --quiet -g guidefile
sort relaxed_solution > rs_sorted
if diff -w -B ex_relaxed_solution rs_sorted ; then
   echo "  PASS: Got expected solution."
   echo ""
else
   echo "  FAIL: Unexpected solution. Test failed. Exiting."
   exit 1
fi
popd

echo "Test: Mitchell web example"
pushd ../examples/web_mitchell
dwsolver --no-write-final-master --quiet -g guidefile
sort relaxed_solution > rs_sorted
if diff -w -B ex_relaxed_solution rs_sorted ; then
   echo "  PASS: Got expected solution."
   echo ""
else
   echo "  FAIL: Unexpected solution. Test failed. Exiting."
   exit 1
fi
popd

echo "Test: Trick web example"
pushd ../examples/web_trick
# Verify by objective value (-40.0) rather than variable assignment, since
# multiple optimal bases may exist.
dwsolver -g guidefile > out_obj.txt
_rc=$?
if [ "$_rc" -ne 0 ]; then
   echo "  FAIL: dwsolver exited with code $_rc. Test failed. Exiting."
   popd
   exit 1
fi
_obj=$(grep "Master objective value" out_obj.txt | tail -1)
if echo "$_obj" | grep -qF -- "-4.000000e+01"; then
   echo "  PASS: Got expected objective value."
   echo ""
else
   echo "  FAIL: Unexpected objective value. Got: $_obj. Test failed. Exiting."
   popd
   exit 1
fi
popd

echo "Test: Four seas example"
pushd ../examples/four_sea
dwsolver -g guidefile > out_obj.txt
_rc=$?
if [ "$_rc" -ne 0 ]; then
   echo "  FAIL: dwsolver exited with code $_rc. Test failed. Exiting."
   popd
   exit 1
fi
_obj=$(grep "Master objective value" out_obj.txt | tail -1)
if echo "$_obj" | grep -qF "1.200000e+01"; then
   echo "  PASS: Got expected objective value."
   echo ""
else
   echo "  FAIL: Unexpected objective value. Got: $_obj. Test failed. Exiting."
   popd
   exit 1
fi
popd

echo "Test: Dantzig textbook example"
pushd ../examples/book_dantzig
dwsolver -g guidefile > out_obj.txt
_rc=$?
if [ "$_rc" -ne 0 ]; then
   echo "  FAIL: dwsolver exited with code $_rc. Test failed. Exiting."
   popd
   exit 1
fi
_obj=$(grep "Master objective value" out_obj.txt | tail -1)
if echo "$_obj" | grep -qF "6.357895e+01"; then
   echo "  PASS: Got expected objective value."
   echo ""
else
   echo "  FAIL: Unexpected objective value. Got: $_obj. Test failed. Exiting."
   popd
   exit 1
fi
popd

echo "Test: single_sub (num_clients=1 semaphore path)"
pushd ../examples/single_sub
dwsolver --no-write-final-master --quiet -g guidefile
sort relaxed_solution > rs_sorted
if diff -w -B ex_relaxed_solution rs_sorted ; then
   echo "  PASS: Got expected solution."
   echo ""
else
   echo "  FAIL: Unexpected solution. Test failed. Exiting."
   exit 1
fi
popd

echo "Test: one_iter (fast Phase II convergence path)"
pushd ../examples/one_iter
dwsolver -g guidefile > out_obj.txt
_rc=$?
if [ "$_rc" -ne 0 ]; then
   echo "  FAIL: dwsolver exited with code $_rc. Test failed. Exiting."
   popd
   exit 1
fi
_obj=$(grep "Master objective value" out_obj.txt | tail -1)
if echo "$_obj" | grep -qF -- "-6.000000e+00"; then
   echo "  PASS: Got expected objective value."
   echo ""
else
   echo "  FAIL: Unexpected objective value. Got: $_obj. Test failed. Exiting."
   popd
   exit 1
fi
popd

echo "Test: neg_y (GLP_LO coupling / y-accumulator sign path)"
pushd ../examples/neg_y
dwsolver --no-write-final-master --quiet -g guidefile
sort relaxed_solution > rs_sorted
if diff -w -B ex_relaxed_solution rs_sorted ; then
   echo "  PASS: Got expected solution."
   echo ""
else
   echo "  FAIL: Unexpected solution. Test failed. Exiting."
   exit 1
fi
popd

