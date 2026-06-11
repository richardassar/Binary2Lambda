# Contributing

This project is four independent implementations of one specification — a
Python reference, a C++ port, a Rust port and a Wolfram Language port. The
binding requirement is that they agree: identical decode, encode and compress
output for every input, and interchangeable table files. A change is correct
only when every implementation still agrees.

## Building and testing

No build is needed for Python. C++ uses CMake, Rust uses Cargo; both are
single dependency-free files.

```sh
# Python reference + compression self-tests
python3 python/lambda_bijection.py
python3 python/lambda_compress.py

# C++ (CMake)
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure

# Rust (Cargo)
cargo test --manifest-path rust/Cargo.toml
cargo run  --release --manifest-path rust/Cargo.toml          # self-test + demo

# Wolfram Language
wolframscript -code 'Get["wolfram/LambdaBinarization.wl"]; Print[LambdaBijectionSelfTest[]]'
```

Each self-test validates counts against exhaustive enumeration, round-trips
both directions across several caps, exercises table save/load with
corruption detection, and round-trips the compression coder. The C++ and Rust
binaries also build with a single compiler command (see the README) if you
prefer not to use CMake/Cargo.

## Cross-language checks

Three checks guard the agreement invariant. Run them after any change to the
bijection or the coder:

```sh
# 1. decode vectors identical (cap inf and cap 3)
diff <(cpp/lambda_bijection_cpp --vectors) <(rust/lambda_bijection_rs --vectors)

# 2. compression bytes identical on the shared vector set
diff <(cpp/lambda_bijection_cpp --compress-vectors) \
     <(rust/lambda_bijection_rs --compress-vectors)

# 3. differential fuzz: same blobs -> same decode in all three languages
bash verify/run_fuzz.sh           # see verify/README.md
```

The Python vectors must match too: `python3 python/lambda_compress.py --vectors`
prints the compression vectors, and `python3 python/lambda_bijection.py` cross-
checks the decode vectors embedded in the self-test.

## Conventions

- **Keep the implementations dependency-free.** The bignum, the coder and the
  table format are written out in each language on purpose. Do not introduce a
  bignum crate, a test framework, or a serialization library.
- **Change all ports together.** A new feature or a fixed bug lands in Python,
  C++, Rust and (where applicable) Wolfram in the same change, with the cross-
  language checks passing.
- **Comments document the present.** A comment states a current invariant or
  constraint a reader would otherwise miss. It does not narrate history or
  contrast with code that is not there; git history records the past.
- **Match the surrounding style** in each language: the existing naming,
  comment density and idiom. The four ports deliberately mirror each other's
  structure.
- **Guard internal invariants with `assert` / `debug_assert!`** (compiled out
  of release), and report bad *input* with exceptions / `Result`, never with a
  crash, hang, or silent wrong answer.

## Adding a test vector

The compression vectors live in `compressionVectorTerms` (C++),
`compression_vector_terms` (Rust) and `_vector_terms` (Python). Add the same
term to all three; the bytes it produces must match across them.
