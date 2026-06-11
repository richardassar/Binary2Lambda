#!/usr/bin/env python3
"""Research harness: measure compressor variants against the corpus and the
general-purpose compressors. Separate from the shipped coder. Variants:

  rank    the bijective rank code (entropy floor for structureless terms)
  base    the adaptive grammar coder (lambda_compress.compress)
  prior   base + a measured counting-measure prior on the kind/index model
  share   prior + subterm sharing (grammar-native repeat references)
"""

from __future__ import annotations

import bz2
import gzip
import lzma

from lambda_bijection import Table, Term, Var, Lam, decode, encode, term_size
from lambda_compress import (_KIND_APP, _KIND_LAM, _KIND_VAR, _RangeEncoder,
                             _encode_symbol, compress)
import compression_benchmark as cb

TABLE = Table()
MIN_SHARE = 6  # min BLC size for a subterm to be worth a back-reference


def rank_bytes_len(term: Term) -> int:
    bits = len(encode(TABLE, term))
    return (bits + 7) // 8


# ----------------------------------------------- measured kind/index prior

def measure_prior(samples: int = 1500):
    """Asymptotic constructor and index frequencies over uniform-random
    terms; seeds the model so short terms start near the rank-code optimum."""
    kinds = [[0, 0, 0] for _ in range(4)]
    indices = [0] * 8
    table = Table()
    for i in range(samples):
        term = decode(table, bin(((i * 2654435761) % (1 << 40)) + 1)[3:])
        stack = [(term, 0)]
        while stack:
            t, m = stack.pop()
            b = min(m, 3)
            if isinstance(t, Var):
                kinds[b][_KIND_VAR] += 1
                indices[min(t.index - 1, 7)] += 1
            elif isinstance(t, Lam):
                kinds[b][_KIND_LAM] += 1
                stack.append((t.body, m + 1))
            else:
                kinds[b][_KIND_APP] += 1
                stack.append((t.fun, m))
                stack.append((t.arg, m))

    def scale(row, total):
        s = sum(row) or 1
        return [max(1, round(total * v / s)) for v in row]

    return [scale(row, 24) for row in kinds], scale(indices, 32)


# --------------------------------------------- coder with prior + sharing

class _Model:
    def __init__(self, kind_prior, index_prior, sharing):
        self._kinds = [list(row) for row in kind_prior]
        self._indices = list(index_prior)
        self.sharing = sharing
        self.repeat = [6, 1]        # adaptive [literal, repeat]
        self.refs: list[int] = []   # adaptive frequency over seen ids
        self.seen: dict = {}        # subterm -> id

    def kind_weights(self, m):
        v, lam, app = self._kinds[min(m, 3)]
        return [v if m >= 1 else 0, lam, app]

    def saw_kind(self, m, kind):
        self._kinds[min(m, 3)][kind] += 16

    def index_weights(self, m):
        return [self._indices[min(i, 7)] for i in range(m)]

    def saw_index(self, index):
        self._indices[min(index - 1, 7)] += 16

    def register(self, term):
        if self.sharing and term_size(term) >= MIN_SHARE and term not in self.seen:
            self.seen[term] = len(self.refs)
            self.refs.append(1)


def _enc(term, m, coder, model):
    if model.sharing:
        if term_size(term) >= MIN_SHARE and term in model.seen:
            _encode_symbol(coder, model.repeat, 1)
            model.repeat[1] += 8
            idx = model.seen[term]
            _encode_symbol(coder, model.refs, idx)
            model.refs[idx] += 16
            return
        _encode_symbol(coder, model.repeat, 0)
        model.repeat[0] += 8
    if isinstance(term, Var):
        _encode_symbol(coder, model.kind_weights(m), _KIND_VAR)
        model.saw_kind(m, _KIND_VAR)
        _encode_symbol(coder, model.index_weights(m), term.index - 1)
        model.saw_index(term.index)
    elif isinstance(term, Lam):
        _encode_symbol(coder, model.kind_weights(m), _KIND_LAM)
        model.saw_kind(m, _KIND_LAM)
        _enc(term.body, m + 1, coder, model)
    else:
        _encode_symbol(coder, model.kind_weights(m), _KIND_APP)
        model.saw_kind(m, _KIND_APP)
        _enc(term.fun, m, coder, model)
        _enc(term.arg, m, coder, model)
    model.register(term)


def variant_len(term, kind_prior, index_prior, sharing):
    coder = _RangeEncoder()
    _enc(term, 0, coder, _Model(kind_prior, index_prior, sharing))
    return len(coder.finish())


# -------------------------------------------------------------- measure

def main():
    kind_prior, index_prior = measure_prior()
    print("kind prior (depth 0..3):", kind_prior)
    print("index prior:", index_prior)

    # bounded corpora so the rank-code table build stays fast (the point is to
    # compare coders, not to stress huge terms)
    structured = [t for name, t in cb.structured_corpus()
                  if not (name.startswith("S^") and int(name[2:]) > 16)]
    uniform = []
    seen_table = Table()
    for length in range(8, 60, 2):
        for j in range(4):
            idx = (1 << length) + j * ((1 << length) // 4 + 1)
            uniform.append(decode(seen_table, bin(idx + 1)[3:]))
    corpora = {"structured": structured, "uniform": uniform}
    general = {"gzip": lambda b: gzip.compress(b, 9),
               "bzip2": lambda b: bz2.compress(b, 9),
               "lzma": lambda b: lzma.compress(b, preset=9)}

    for name, terms in corpora.items():
        blc_per = sum(len(cb.blc_bytes(t)) for t in terms)
        blob = b"".join(cb.blc_bytes(t) for t in terms)
        print(f"\n=== {name} ({len(terms)} terms), total bytes ===")
        print(f"  raw BLC (per-term)         {blc_per:>7}  100%")
        for label, fn in general.items():
            per = sum(len(fn(cb.blc_bytes(t))) for t in terms)
            arch = len(fn(blob))
            print(f"  BLC+{label:<5} per-term {per:>7}  arch {arch:>7}")
        for label, sizes in (
                ("rank code", [rank_bytes_len(t) for t in terms]),
                ("rank base", [len(compress(t)) for t in terms]),
                ("rank prior",
                 [variant_len(t, kind_prior, index_prior, False) for t in terms]),
                ("rank share",
                 [variant_len(t, kind_prior, index_prior, True) for t in terms])):
            total = sum(sizes)
            print(f"  {label:<24}   {total:>7}  {100*total/blc_per:>3.0f}%")


if __name__ == "__main__":
    main()
