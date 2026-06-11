#!/usr/bin/env python3
"""Cross-language differential fuzz harness for the compression decoder.

Two modes:

  gen <dir>          write garbage.txt, valid.txt and valid_expected.txt
  run <blobfile>     decompress every hex blob, print "OK\\t<digest>\\t<nodes>"
                     or "ERR" per line (the Python reference oracle)

The structural digest matches the C++ and Rust `--fuzz` modes byte for byte,
so feeding the same blob file to all three and diffing the outputs proves the
three decoders agree on every input (valid and adversarial), not just on the
fixed compression vectors.
"""

import os
import random
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "python"))

from lambda_bijection import App, Lam, Table, Var, decode  # noqa: E402
import lambda_compress as lc  # noqa: E402

_FNV_PRIME = 0x100000001B3
_FNV_OFFSET = 0xCBF29CE484222325
_MASK64 = (1 << 64) - 1


def digest(root):
    """Iterative pre-order FNV-1a 64-bit digest; identical across languages."""
    h = _FNV_OFFSET
    nodes = 0
    stack = [root]
    while stack:
        term = stack.pop()
        nodes += 1
        if isinstance(term, Var):
            h = ((h ^ 0x56) * _FNV_PRIME) & _MASK64
            index = term.index & 0xFFFFFFFF
            for shift in (0, 8, 16, 24):
                h = ((h ^ ((index >> shift) & 0xFF)) * _FNV_PRIME) & _MASK64
        elif isinstance(term, Lam):
            h = ((h ^ 0x4C) * _FNV_PRIME) & _MASK64
            stack.append(term.body)
        else:
            h = ((h ^ 0x41) * _FNV_PRIME) & _MASK64
            stack.append(term.arg)
            stack.append(term.fun)
    return h, nodes


def _result_line(blob):
    try:
        term = lc.decompress(blob)
    except Exception:                  # any malformed input -> ERR, never crash
        return "ERR"
    h, nodes = digest(term)
    return f"OK\t{h:016x}\t{nodes}"


def _lam_chain(n):
    body = Var(1)
    for _ in range(n):
        body = Lam(body)
    return body


def _church(n):
    body = Var(1)
    for _ in range(n):
        body = App(Var(2), body)
    return Lam(Lam(body))


def gen(out_dir, scale=1):
    """Write the blob files. `scale` multiplies the random-sample counts so a
    smoke run (scale << 1) and a thorough sweep (scale >> 1) share one path;
    the exhaustive structured blobs are emitted regardless of scale."""
    os.makedirs(out_dir, exist_ok=True)
    rng = random.Random(20260612)
    per_len = max(1, round(40 * scale))

    garbage = []
    # every one-byte blob
    garbage += [bytes([b]).hex() for b in range(256)]
    # two-byte blobs across a spread of leading and trailing bytes
    for first in list(range(16)) + [0x7F, 0x80, 0xFF]:
        for second in (0x00, 0x01, 0x7F, 0x80, 0xAA, 0x55, 0xFF):
            garbage.append(bytes([first, second]).hex())
    # constant-fill runs of varied leading byte and length, which drive long
    # decode chains
    for lead in range(8):
        for fill in (0x00, 0xFF, 0xAA, 0x55):
            for length in range(1, 25):
                garbage.append(bytes([lead] + [fill] * length).hex())
    # uniform random blobs across a range of lengths
    for length in range(1, 41):
        for _ in range(per_len):
            garbage.append(bytes(rng.randrange(256)
                                 for _ in range(length)).hex())

    valid_terms = []
    table = Table()
    for length in range(0, 200):
        for _ in range(max(1, round(6 * scale))):
            bits = "".join(rng.choice("01") for _ in range(length))
            valid_terms.append(decode(table, bits))
    valid_terms += [_lam_chain(48), _lam_chain(1000), _church(2000),
                    _church(50), _lam_chain(300)]

    valid = []
    expected = []
    for term in valid_terms:
        blob = lc.compress(term)
        if not blob:
            continue  # an all-zero stream (e.g. Lam(Var(1))) has no line form
        valid.append(blob.hex())
        h, nodes = digest(term)
        expected.append(f"OK\t{h:016x}\t{nodes}")

    # truncations of valid blobs are adversarial-but-near-valid garbage
    for hexblob in valid[:max(1, round(400 * scale))]:
        raw = bytes.fromhex(hexblob)
        for drop in (1, 2, 3):
            if len(raw) > drop:
                garbage.append(raw[:-drop].hex())

    with open(os.path.join(out_dir, "garbage.txt"), "w") as handle:
        handle.write("\n".join(garbage) + "\n")
    with open(os.path.join(out_dir, "valid.txt"), "w") as handle:
        handle.write("\n".join(valid) + "\n")
    with open(os.path.join(out_dir, "valid_expected.txt"), "w") as handle:
        handle.write("\n".join(expected) + "\n")
    print(f"garbage blobs: {len(garbage)}")
    print(f"valid blobs:   {len(valid)}")


def run(path):
    out = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            try:
                blob = bytes.fromhex(line)
            except ValueError:
                out.append("ERR")
                continue
            out.append(_result_line(blob))
    sys.stdout.write("\n".join(out) + ("\n" if out else ""))


def main():
    if len(sys.argv) >= 3 and sys.argv[1] == "gen":
        scale = float(sys.argv[3]) if len(sys.argv) >= 4 else 1.0
        gen(sys.argv[2], scale)
    elif len(sys.argv) >= 3 and sys.argv[1] == "run":
        run(sys.argv[2])
    else:
        sys.exit("usage: fuzz_harness.py gen <dir> [scale] | run <blobfile>")


if __name__ == "__main__":
    main()
