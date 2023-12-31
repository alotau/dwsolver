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

