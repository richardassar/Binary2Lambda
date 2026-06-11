# Cross-language differential fuzz

`compress`/`decompress` are implemented three times (Python, C++, Rust) as
independent ports of one deterministic specification. For any input bytes,
`decompress` is therefore a pure function whose result must be identical in
every language — not only on valid `compress` outputs but on arbitrary and
adversarial bytes. This harness checks that.

## What it does

`fuzz_harness.py gen` writes blob files; the C++/Rust `--fuzz` mode and
`fuzz_harness.py run` decode each blob and emit, per line, either

```
OK<TAB><digest><TAB><nodes>
```

or `ERR`. The digest is a structural FNV-1a-64 over the decoded term,
computed identically in all three languages, so two outputs are equal exactly
when the decoded terms are equal. `run_fuzz.sh` feeds one blob file to all
three and diffs the outputs.

Two blob sets:

- **garbage** — every one-byte blob, header bytes spanning the valid and
  rejected ranges, long constant-fill runs that drive deep interval-narrowing
  chains, uniform random blobs across many lengths, and truncations of valid
  blobs. The three decoders must produce byte-identical results (same decoded
  term or all `ERR`); any difference is a cross-implementation bug.
- **valid** — `compress` outputs of random and structured terms. Each must
  decode back to the original term (the harness checks the per-term digests
  against `valid_expected.txt`) and do so identically in all three languages.

This exercises the paths static review cannot: end-of-stream handling, the
node/depth termination guards, and the 32-bit range-coder arithmetic and byte
renormalization under out-of-range inputs.

## Running

```sh
bash verify/run_fuzz.sh            # default sample size
bash verify/run_fuzz.sh 4          # 4x the random samples (thorough)
bash verify/run_fuzz.sh 0.2        # smaller, faster (CI smoke)
```

Build the C++ and Rust binaries first (the script expects
`cpp/lambda_bijection_cpp` and `rust/lambda_bijection_rs`). Scratch files are
written under `verify/work/` (git-ignored). A non-zero exit means a
divergence or a round-trip mismatch, with the differing lines printed.

## Files

- `fuzz_harness.py` — blob generator (`gen`) and Python reference oracle
  (`run`); defines the shared structural digest.
- `run_fuzz.sh` — builds the blob files, runs all three decoders, diffs.
