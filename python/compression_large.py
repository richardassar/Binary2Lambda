#!/usr/bin/env python3
"""Large-term compression behaviour: are there discontinuities in gzip / lzma
as terms grow by one or two orders of magnitude in size?

Big random closed terms are built directly (no counting table needed), then
each is BLC-encoded and compressed.  Plotted vs node count: raw BLC, gzip,
bzip2, lzma, and a table-free grammar coder.  (The rank code is omitted at
this scale - it needs the full counting table built to the term's size.)

Run: python3 compression_large.py        (writes plots/compression_large.png)
"""

from __future__ import annotations

import bz2
import gzip
import lzma

from lambda_bijection import App, Lam, Term, Var
from lambda_compress import compress
from compression_benchmark import blc_bytes, node_count


def gen_term(budget: int, depth: int, state: list[int]) -> tuple[Term, int]:
    """A deterministic random closed term of roughly `budget` nodes; every
    variable is bound (index <= enclosing lambdas). The term reliably grows to
    its budget: a leaf (Var) is emitted only when the budget is spent."""
    state[0] = (state[0] * 6364136223846793005 + 1442695040888963407) \
        & ((1 << 64) - 1)
    r = state[0]
    if budget <= 1 and depth >= 1:
        return Var(1 + (r >> 2) % depth), 1
    if depth == 0 or r % 3 == 0:                 # abstraction (always at root)
        body, used = gen_term(budget - 1, depth + 1, state)
        return Lam(body), used + 1
    left, lu = gen_term(max(1, (budget - 1) // 2), depth, state)
    right, ru = gen_term(max(1, budget - 1 - lu), depth, state)
    return App(left, right), lu + ru + 1


def main() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pathlib import Path

    general = {"BLC + gzip": (lambda b: gzip.compress(b, 9), "tab:blue", "--"),
               "BLC + bzip2": (lambda b: bz2.compress(b, 9), "tab:green", "--"),
               "BLC + lzma": (lambda b: lzma.compress(b, preset=9), "tab:cyan",
                              "--")}

    state = [0x12345678ABCDEF01]
    records = []  # (node_count, {method: bytes}) per term, sorted by size below
    # node counts from ~30 up to ~70000 (BLC bytes past gzip's 32 KB window)
    for target in [30, 60, 120, 250, 500, 1000, 2000, 4000, 8000, 16000,
                   24000, 32000, 48000, 64000]:
        term, _ = gen_term(target, 0, state)
        blc = blc_bytes(term)
        row = {"raw BLC": len(blc), "grammar coder": len(compress(term))}
        for name, (fn, _c, _s) in general.items():
            row[name] = len(fn(blc))
        records.append((node_count(term), row))
    records.sort(key=lambda rec: rec[0])
    sizes = [n for n, _ in records]
    data = {k: [row[k] for _, row in records]
            for k in ("raw BLC", *general, "grammar coder")}

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    ax.plot(sizes, data["raw BLC"], "-", color="tab:gray", label="raw BLC")
    for name, (_fn, color, style) in general.items():
        ax.plot(sizes, data[name], style, color=color, marker="o",
                markersize=3, label=name)
    ax.plot(sizes, data["grammar coder"], "-", color="tab:red", linewidth=2,
            label="grammar coder")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("term size (nodes: variables + abstractions + applications)")
    ax.set_ylabel("bytes to store one term")
    ax.set_title("Large-term compression: any gzip / lzma discontinuities?\n"
                 "(single random closed term, log-log)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out = Path(__file__).resolve().parents[1] / "plots" / "compression_large.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)

    print(f"{'nodes':>7} {'BLC':>8} {'gzip':>8} {'bzip2':>8} {'lzma':>8} "
          f"{'coder':>8}")
    for i, n in enumerate(sizes):
        print(f"{n:>7} {data['raw BLC'][i]:>8} {data['BLC + gzip'][i]:>8} "
              f"{data['BLC + bzip2'][i]:>8} {data['BLC + lzma'][i]:>8} "
              f"{data['grammar coder'][i]:>8}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
