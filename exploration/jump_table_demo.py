#!/usr/bin/env python3
"""Anatomy of one size class and a full node-by-node trace of both
directions of the map through it.

Class (n=18, m=0): closed terms of BLC size 18. Shows the split prefix-sum
row (SPLIT), the jump table built over it (JUMP), then a traced decode and
the traced encode of the same term, listing every table cell touched.
"""

import sys

sys.setrecursionlimit(100000)

from bijective_lambda import CUM, T, extend_tables, size, show

extend_tables(40)

_rows = {}


def split_row(n, m):
    m = min(m, n - 1)
    if (n, m) not in _rows:
        pre, acc, anatomy = [], 0, []
        for k in range(2, n - 3):
            left, right = T(k, m), T(n - 2 - k, m)
            acc += left * right
            pre.append(acc)
            anatomy.append((k, left, right, left * right, acc))
        B = len(pre)
        W = pre[-1] // B + 1 if B else 0
        J, kk = [], 0
        for j in range(B):
            while kk < B - 1 and pre[kk] <= j * W:
                kk += 1
            J.append(kk)
        _rows[(n, m)] = (pre, anatomy, W, J)
    return _rows[(n, m)]


# ------------------------------------------------- A. table contents

N0, M0 = 18, 0
pre, anatomy, W, J = split_row(N0, M0)
print(f"A. CLASS (n={N0}, m={M0}): T = {T(N0, M0)} closed terms of size {N0}")
print(f"   layout of [0, {T(N0, M0)}): vars 0 | lams [0, {T(N0 - 1, M0 + 1) if False else T(N0 - 2, M0 + 1)}) | apps [{T(N0 - 2, M0 + 1)}, {T(N0 - 2, M0 + 1) + pre[-1]})")
print()
print("   SPLIT row (apps by left-subterm size k):")
print(f"   {'k':>4} {'T(k,0)':>7} {'T(16-k,0)':>10} {'block':>6} {'pre (cum)':>10}")
for k, left, right, block, cum in anatomy:
    print(f"   {k:>4} {left:>7} {right:>10} {block:>6} {cum:>10}")
print()
print(f"   JUMP over that row: {len(J)} buckets of rank-mass W = {W}")
print(f"   {'bucket j':>9} {'covers ranks':>13} {'J[j]':>5} {'-> first candidate k':>21}")
for j, idx in enumerate(J):
    lo, hi = j * W, min((j + 1) * W - 1, pre[-1] - 1)
    if lo > hi:
        rng = "(empty)"
    else:
        rng = f"{lo}..{hi}"
    print(f"   {j:>9} {rng:>13} {idx:>5} {anatomy[idx][0]:>21}")

# ------------------------------------------------- B. traced decode


def unrank_traced(r, n, m, depth):
    pad = "      " + "    " * depth
    nvar = 1 if n - 1 <= m else 0
    if nvar:
        if r == 0:
            print(f"{pad}(n={n},m={m}) r=0 -> Var({n - 1})    [1 compare]")
            return ('var', n - 1)
        r -= 1
    lams = T(n - 2, m + 1)
    if r < lams:
        print(f"{pad}(n={n},m={m}) r={r} < lams T({n - 2},{m + 1})={lams} -> Lam, recurse"
              f"    [{1 + nvar} compares]")
        return ('lam', unrank_traced(r, n - 2, m + 1, depth + 1))
    r0, r = r, r - lams
    p, anat, w, jj = split_row(n, m)
    j = min(r // w, len(jj) - 1)
    idx = jj[j]
    probes = 1
    while p[idx] <= r:
        idx += 1
        probes += 1
    k = anat[idx][0]
    base = p[idx - 1] if idx > 0 else 0
    right = T(n - 2 - k, m)
    q, rr = divmod(r - base, right)
    print(f"{pad}(n={n},m={m}) r={r0} >= lams {lams} -> App, app-rank {r}:"
          f" j={r}//{w}={j}, J[{j}]={idx if probes == 1 else jj[j]},"
          f" {probes} probe{'s' if probes > 1 else ''} -> k={k};"
          f" divmod({r}-{base}, T({n - 2 - k},{m})={right}) = ({q},{rr})")
    return ('app', unrank_traced(q, k, m, depth + 1),
            unrank_traced(rr, n - 2 - k, m, depth + 1))


lams18 = T(N0 - 2, M0 + 1)
r_app = 37
N = CUM[N0 - 1] + lams18 + r_app
s = bin(N + 1)[3:]
print()
print(f"B. DECODE the {len(s)}-bit string {s}  (N = {N})")
print(f"      size class: CUM[{N0 - 1}]={CUM[N0 - 1]} <= N < CUM[{N0}]={CUM[N0]} -> n={N0};"
      f" in-class rank {N - CUM[N0 - 1]}")
t = unrank_traced(N - CUM[N0 - 1], N0, 0, 0)
print(f"      term: {show(t)}")

# ------------------------------------------------- C. traced encode


def rank_traced(t, m, depth):
    pad = "      " + "    " * depth
    if t[0] == 'var':
        print(f"{pad}Var({t[1]}) -> 0")
        return 0
    n = size(t)
    nvar = 1 if n - 1 <= m else 0
    if t[0] == 'lam':
        print(f"{pad}Lam (n={n},m={m}) -> {nvar} + rank(body)")
        return nvar + rank_traced(t[1], m + 1, depth + 1)
    f, a = t[1], t[2]
    fs = size(f)
    p, anat, _, _ = split_row(n, m)
    idx = fs - 2
    base = p[idx - 1] if idx > 0 else 0
    lams = T(n - 2, m + 1)
    rb = T(size(a), m)
    print(f"{pad}App (n={n},m={m}), left size {fs} KNOWN -> no search:"
          f" {nvar} + lams T({n - 2},{m + 1})={lams} + pre[{idx - 1}]={base}"
          f" + rank(left)*T({size(a)},{m})={rb} + rank(right)")
    return (nvar + lams + base
            + rank_traced(f, m, depth + 1) * rb + rank_traced(a, m, depth + 1))


print()
print(f"C. ENCODE the same term back")
r = rank_traced(t, 0, 0)
N2 = CUM[size(t) - 1] + r
print(f"      in-class rank {r}; N = CUM[{size(t) - 1}]={CUM[size(t) - 1]} + {r} = {N2};"
      f" string = {bin(N2 + 1)[3:]}")
assert N2 == N
