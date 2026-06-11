#!/usr/bin/env python3
"""Binary2Lambda compression comparison.

Two ways to write a closed lambda term as raw bytes:

  BLC   - Tromp's binary lambda calculus code (var i = 1^i 0, abstraction
          = 00 <body>, application = 01 <fun> <arg>), packed into bytes.
          Structure-preserving, so an LZ-family compressor can match on it.
  rank  - the bijective rank code: the term's index in the canonical
          enumeration, in binary (what encode() returns).  Maximally dense
          (no wasted strings), but it scrambles structure into one integer,
          so general compressors find nothing to exploit.

For each representation the raw size and the size after gzip,
bzip2 and LZMA/xz (all Python standard library — no external dependency).

Separately, this project's own compressor (lambda_compress) is NOT a general
bitstream compressor: it consumes the *term* through the grammar with an
adaptive model.  It cannot be "run on BLC bytes"; it is its own row.

Two regimes: per term (store one term, independently decodable - the genetic
programming genome unit) and whole corpus (serialize every term, concatenate,
compress once - where a general compressor amortizes its container and can
match repetition across terms).

Run: python3 compression_benchmark.py          (prints the tables)
     python3 compression_benchmark.py --plot    (also writes plots/*.png)
"""

from __future__ import annotations

import bz2
import gzip
import lzma

from lambda_bijection import App, Lam, Table, Term, Var, decode, encode
from lambda_compress import compress

TABLE = Table()


# ---------------------------------------------------- representations

def node_count(term: Term) -> int:
    """Number of constructors (Var + Lam + App) in the term's syntax tree.
    A representation-neutral size measure (the 'natural size' of a lambda
    term) - unlike a BLC bit count, it does not privilege any encoding or
    charge for de Bruijn index magnitude."""
    if isinstance(term, Var):
        return 1
    if isinstance(term, Lam):
        return 1 + node_count(term.body)
    return 1 + node_count(term.fun) + node_count(term.arg)


def _pack(bits: str) -> bytes:
    bits += "0" * (-len(bits) % 8)
    return int(bits, 2).to_bytes(len(bits) // 8, "big") if bits else b""


def blc_bits(term: Term) -> str:
    if isinstance(term, Var):
        return "1" * term.index + "0"
    if isinstance(term, Lam):
        return "00" + blc_bits(term.body)
    return "01" + blc_bits(term.fun) + blc_bits(term.arg)


def blc_bytes(term: Term) -> bytes:
    return _pack(blc_bits(term))


def rank_bytes(term: Term) -> bytes:
    return _pack(encode(TABLE, term))


GENERAL = {
    "gzip": lambda data: gzip.compress(data, 9),
    "bzip2": lambda data: bz2.compress(data, 9),
    "lzma": lambda data: lzma.compress(data, preset=9),
}

# the comparison matrix: representation x {raw, gzip, bzip2, lzma}, plus the
# grammar-aware coder as its own column.
COLUMNS = (["BLC", "BLC+gzip", "BLC+bzip2", "BLC+lzma",
            "rank", "rank+gzip", "rank+bzip2", "rank+lzma", "coder"])


def per_term_bytes(term: Term) -> dict[str, int]:
    blc, rank = blc_bytes(term), rank_bytes(term)
    sizes = {"BLC": len(blc), "rank": len(rank), "coder": len(compress(term))}
    for name, fn in GENERAL.items():
        sizes[f"BLC+{name}"] = len(fn(blc))
        sizes[f"rank+{name}"] = len(fn(rank))
    return sizes


# ------------------------------------------------------------ corpora

def church(n: int) -> Term:
    body: Term = Var(1)
    for _ in range(n):
        body = App(Var(2), body)
    return Lam(Lam(body))


def combinators() -> dict[str, Term]:
    """Standard combinators and Church-encoded data/operations, written in
    de Bruijn form (indices from 1)."""
    return {
        "I": Lam(Var(1)),
        "K": Lam(Lam(Var(2))),
        "S": Lam(Lam(Lam(App(App(Var(3), Var(1)), App(Var(2), Var(1)))))),
        "B": Lam(Lam(Lam(App(Var(3), App(Var(2), Var(1)))))),
        "C": Lam(Lam(Lam(App(App(Var(3), Var(1)), Var(2))))),
        "W": Lam(Lam(App(App(Var(2), Var(1)), Var(1)))),
        "true": Lam(Lam(Var(2))),
        "false": Lam(Lam(Var(1))),
        "pair": Lam(Lam(Lam(App(App(Var(1), Var(3)), Var(2))))),
        "not": Lam(App(App(Var(1), Lam(Lam(Var(1)))), Lam(Lam(Var(2))))),
        "add": Lam(Lam(Lam(Lam(App(App(Var(4), Var(2)),
                                   App(App(Var(3), Var(2)), Var(1))))))),
        "mul": Lam(Lam(Lam(App(Var(3), App(Var(2), Var(1)))))),
        "Y": Lam(App(Lam(App(Var(2), App(Var(1), Var(1)))),
                     Lam(App(Var(2), App(Var(1), Var(1)))))),
        "omega": App(Lam(App(Var(1), Var(1))), Lam(App(Var(1), Var(1)))),
    }


def structured_corpus() -> list[tuple[str, Term]]:
    """Combinators, Church numerals, Church-arithmetic applications, and
    deliberately repetitive terms - data with real intra/cross-term
    regularity, like a library of small programs."""
    combs = combinators()
    add, mul = combs["add"], combs["mul"]
    terms = [(name, t) for name, t in combs.items()]
    terms += [(f"church{n}", church(n)) for n in (0, 1, 2, 3, 5, 8, 13, 21, 34)]
    # applied programs (stored unreduced)
    terms += [("add 2 3", App(App(add, church(2)), church(3))),
              ("mul 3 4", App(App(mul, church(3)), church(4))),
              ("add(add 2 2)1", App(App(add, App(App(add, church(2)),
                                                 church(2))), church(1)))]
    repeated: Term = combs["S"]
    for i in range(1, 5):  # S^2 .. S^16; larger repeats blow up the rank table
        repeated = App(repeated, repeated)
        terms.append((f"S^{2**i}", repeated))
    return terms


def uniform_corpus() -> list[tuple[str, Term]]:
    """A deterministic, reproducible population of *diverse* uniformly-random
    closed terms - a stand-in for a genetic-programming population, and the
    honest worst case for any compressor (random terms are near-incompressible).

    Indices are drawn by a fixed-seed LCG spread across the full range of each
    target bit length, so the terms are varied (not near-duplicate
    low-rank terms, which would let a general compressor cheat)."""
    terms = []
    state = 0x2545F4914F6CDD1D
    for length in range(16, 150, 6):
        for _ in range(6):
            state = (state * 6364136223846793005 + 1442695040888963407) \
                & ((1 << 64) - 1)
            index = (state % (1 << (length - 1))) | (1 << (length - 1))
            terms.append((f"L{length}", decode(TABLE, bin(index + 1)[3:])))
    return terms


# --------------------------------------------------------- reporting

def corpus_totals(terms: list[Term]) -> dict[str, dict[str, int]]:
    per_term = {c: 0 for c in COLUMNS}
    for term in terms:
        for col, size in per_term_bytes(term).items():
            per_term[col] += size

    blc_blob = b"".join(blc_bytes(t) for t in terms)
    rank_blob = b"".join(rank_bytes(t) for t in terms)
    archive = {"BLC": len(blc_blob), "rank": len(rank_blob),
               "coder": sum(len(compress(t)) for t in terms)}
    for name, fn in GENERAL.items():
        archive[f"BLC+{name}"] = len(fn(blc_blob))
        archive[f"rank+{name}"] = len(fn(rank_blob))
    return {"per_term": per_term, "archive": archive}


def print_matrix(title: str, terms: list[Term]) -> None:
    totals = corpus_totals(terms)
    base = totals["per_term"]["BLC"]
    arch_base = totals["archive"]["BLC"]
    print(f"\n{title}  ({len(terms)} terms, totals in bytes; "
          f"% is of raw BLC in that regime)")
    print(f"  {'method':<12} {'per-term store':>18} {'one archive':>16}")
    for col in COLUMNS:
        pt, ar = totals["per_term"][col], totals["archive"][col]
        print(f"  {col:<12} {pt:>9} {100*pt/base:>5.0f}%   "
              f"{ar:>9} {100*ar/arch_base:>5.0f}%")


# ------------------------------------------------------------- plots

def make_plots(corpora: dict[str, list[Term]]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pathlib import Path

    plots = Path(__file__).resolve().parents[1] / "plots"
    series = {  # (column, color, style, linewidth)
        "raw BLC": ("BLC", "tab:gray", "--", 1.2),
        "BLC + gzip": ("BLC+gzip", "tab:blue", "--", 1.2),
        "BLC + lzma": ("BLC+lzma", "tab:cyan", "--", 1.2),
        "rank code": ("rank", "tab:red", "-", 2.6),
        "grammar coder": ("coder", "tab:orange", ":", 1.6),
    }

    # 1. per-term bytes vs term size, over diverse random terms (the GP unit)
    sweep = []
    state = 0x9E3779B97F4A7C15
    for length in range(8, 256, 3):
        state = (state * 6364136223846793005 + 1442695040888963407) \
            & ((1 << 64) - 1)
        idx = (state % (1 << (length - 1))) | (1 << (length - 1))
        sweep.append(decode(TABLE, bin(idx + 1)[3:]))
    sweep.sort(key=node_count)
    xs = [node_count(t) for t in sweep]
    rows = [per_term_bytes(t) for t in sweep]
    fig, ax = plt.subplots(figsize=(8, 5))
    for label, (col, color, style, lw) in series.items():
        ax.plot(xs, [r[col] for r in rows], style, color=color,
                linewidth=lw, label=label)
    ax.set_xlabel("term size (nodes: variables + abstractions + applications)")
    ax.set_ylabel("bytes to store one term")
    ax.set_title("Per-term storage: rank code vs general compressors\n"
                 "(one term, independently decodable - the GP genome unit)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(plots / "compression_per_term.png", dpi=150)
    plt.close(fig)

    # 2. whole-corpus archive totals (grouped bars over the full matrix)
    bar_cols = ["BLC", "BLC+gzip", "BLC+lzma", "rank", "rank+lzma", "coder"]
    labels = ["BLC", "BLC\n+gzip", "BLC\n+lzma", "rank", "rank\n+lzma",
              "grammar\ncoder"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, (name, terms) in zip(axes, corpora.items()):
        totals = corpus_totals(terms)
        per = [totals["per_term"][c] / 1024 for c in bar_cols]
        arch = [totals["archive"][c] / 1024 for c in bar_cols]
        x = range(len(bar_cols))
        ax.bar([i - 0.2 for i in x], per, width=0.4, label="per-term store",
               color="tab:red", edgecolor="black", linewidth=0.4)
        ax.bar([i + 0.2 for i in x], arch, width=0.4, label="one archive",
               color="tab:blue", edgecolor="black", linewidth=0.4)
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("KiB")
        ax.set_title(f"{name} corpus ({len(terms)} terms)")
        ax.grid(True, axis="y", alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("Corpus storage by method (lower is better)")
    fig.tight_layout()
    fig.savefig(plots / "compression_corpus.png", dpi=150)
    plt.close(fig)
    print(f"\nwrote {plots}/compression_per_term.png and "
          f"{plots}/compression_corpus.png")


# -------------------------------------------------------------- main

def main() -> None:
    import sys

    print("Single terms (bytes to store one, independently decodable):")
    print(f"  {'term':<16} {'BLC':>5} {'rank':>5} {'coder':>6} "
          f"{'BLC+gz':>7} {'BLC+xz':>7}")
    combs = combinators()
    samples = [("I = λx.x", combs["I"]), ("S", combs["S"]),
               ("add 2 3", App(App(combs["add"], church(2)), church(3))),
               ("church 20", church(20)), ("church 100", church(100)),
               ("S^16 (repeat)", structured_corpus()[-1][1])]
    for label, term in samples:
        s = per_term_bytes(term)
        print(f"  {label:<16} {s['BLC']:>5} {s['rank']:>5} {s['coder']:>6} "
              f"{s['BLC+gzip']:>7} {s['BLC+lzma']:>7}")

    corpora = {"structured": [t for _, t in structured_corpus()],
               "uniform-random": [t for _, t in uniform_corpus()]}
    for name, terms in corpora.items():
        print_matrix(f"{name} corpus", terms)

    if "--plot" in sys.argv:
        make_plots(corpora)


if __name__ == "__main__":
    main()
