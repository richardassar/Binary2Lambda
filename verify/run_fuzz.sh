#!/usr/bin/env bash
# Cross-language differential fuzz for the compression decoder.
#
# C++ and Rust decode the full blob set (both compiled, fast) and must agree
# byte for byte. The Python reference decodes a representative subset and must
# agree with them there; valid blobs additionally check the round-trip digest.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VDIR="$ROOT/verify/work"
CPP="$ROOT/cpp/lambda_bijection_cpp"
RUST="$ROOT/rust/lambda_bijection_rs"
PY="python3 $ROOT/verify/fuzz_harness.py"
SCALE="${1:-1}"     # multiplies random-sample counts
PYSUB="${2:-2500}"  # garbage blobs the Python reference re-checks

mkdir -p "$VDIR"
echo "== generating blob files (scale=$SCALE) =="
$PY gen "$VDIR" "$SCALE" || exit 1

status=0
fail() { echo "  *** $1 ***"; status=1; }

for set in garbage valid; do
  blob="$VDIR/$set.txt"
  echo "== differential fuzz: $set =="
  "$CPP"  --fuzz "$blob" > "$VDIR/out_cpp_$set.txt"
  "$RUST" --fuzz "$blob" > "$VDIR/out_rust_$set.txt"
  n=$(wc -l < "$VDIR/out_cpp_$set.txt"); ok=$(grep -c '^OK' "$VDIR/out_cpp_$set.txt")
  if diff -q "$VDIR/out_cpp_$set.txt" "$VDIR/out_rust_$set.txt" >/dev/null; then
    echo "  C++ vs Rust: identical over all $n blobs ($ok OK, $((n-ok)) ERR)"
  else
    fail "C++ vs Rust DIVERGENCE in $set"
    diff "$VDIR/out_cpp_$set.txt" "$VDIR/out_rust_$set.txt" | head -20
  fi
done

echo "== Python reference: $PYSUB-blob garbage subset =="
head -n "$PYSUB" "$VDIR/garbage.txt" > "$VDIR/garbage_sub.txt"
$PY run "$VDIR/garbage_sub.txt" > "$VDIR/out_py_garbage_sub.txt"
head -n "$PYSUB" "$VDIR/out_cpp_garbage.txt" > "$VDIR/out_cpp_garbage_sub.txt"
if diff -q "$VDIR/out_py_garbage_sub.txt" "$VDIR/out_cpp_garbage_sub.txt" >/dev/null; then
  echo "  Python agrees with C++/Rust on all $PYSUB"
else
  fail "Python vs C++ DIVERGENCE on garbage subset"
  diff "$VDIR/out_py_garbage_sub.txt" "$VDIR/out_cpp_garbage_sub.txt" | head -20
fi

echo "== Python reference: all valid blobs + round-trip digests =="
$PY run "$VDIR/valid.txt" > "$VDIR/out_py_valid.txt"
if diff -q "$VDIR/out_py_valid.txt" "$VDIR/out_cpp_valid.txt" >/dev/null \
&& diff -q "$VDIR/out_py_valid.txt" "$VDIR/valid_expected.txt" >/dev/null; then
  echo "  all valid blobs round-trip to the original term in every language"
else
  fail "valid round-trip mismatch"
  diff "$VDIR/valid_expected.txt" "$VDIR/out_py_valid.txt" | head -20
fi

[ $status -eq 0 ] && echo "== PASS: all decoders agree ==" || echo "== FAIL =="
exit $status
