#!/usr/bin/env python3
"""Two measurements for the bounded bijection:

1. Incremental table extension: the recurrence T(n,m) references only rows
   with smaller n, so the table is a cache of a prefix of one infinite
   triangular array. Extending the length cap appends rows; nothing already
   computed is touched, and the map on short strings is unchanged.
   Marginal cost of extending cap a -> b scales like (b^3 - a^3)/3 word ops.

2. Jump tables for the app-split search in unrank (the num->lam hot spot):
   bucket the rank mass of each (n,m) split distribution into equal-width
   cells storing the first candidate split; lookup = one division + a short
   scan, expected O(1) probes per app node, replacing the O(log n) bisect.
"""

import sys
import time
from bisect import bisect_right

sys.setrecursionlimit(100000)

from bijective_lambda import TT, CUM, T, extend_tables, size, show, decode

OUT = []


def report(line=""):
    OUT.append(line)


# ---------------------------------------------------------------- part 1

report("1. INCREMENTAL EXTENSION of the count table")
s_probe = "01100101"
probe_before = show(decode(s_probe))     # forces tables to exist for short strings

stamps = []
caps = [100, 200, 400]
for cap in caps:
    t0 = time.perf_counter()
    extend_tables(cap)
    stamps.append(time.perf_counter() - t0)

probe_after = show(decode(s_probe))
assert probe_after == probe_before
report(f"   map stability: decode({s_probe!r}) = {probe_before}  "
       f"(identical before and after extending cap to 400) ✓")

pred = [100 ** 3, 200 ** 3 - 100 ** 3, 400 ** 3 - 200 ** 3]
report(f"   {'cap':>12} {'marginal time':>14} {'measured ratio':>15} {'n^3 model':>10}")
for i, cap in enumerate(caps):
    ratio = stamps[i] / stamps[0]
    report(f"   {'-> ' + str(cap):>12} {stamps[i]:>12.3f} s {ratio:>14.1f}x "
           f"{pred[i] / pred[0]:>9.1f}x")
report(f"   extending 400 -> 401 (one row): "
       f"~{401 * 401 // 2:,} multiply-adds, microseconds")

# table footprint at two caps
for cap in (76, 200):
    entries = sum(len(TT[n]) for n in range(2, cap + 1))
    bits = sum(v.bit_length() for n in range(2, cap + 1) for v in TT[n])
    report(f"   T-table footprint at cap n={cap}: {entries:,} entries, "
           f"{bits / 8 / 1024:.0f} KiB of count data")

# ---------------------------------------------------------------- part 2

L1_MAX = 60
N_MAX = (1 << (L1_MAX + 1)) - 2
n_max = 4
while CUM[n_max] <= N_MAX:
    n_max += 1

SPLIT = {}
for n in range(6, n_max + 1):
    for m in range(n):
        acc = 0
        pre = []
        for k in range(2, n - 3):
            acc += T(k, m) * T(n - 2 - k, m)
            pre.append(acc)
        SPLIT[(n, m)] = pre

# jump tables: equal rank-mass buckets over each split distribution
t0 = time.perf_counter()
JUMP = {}
for key, pre in SPLIT.items():
    if not pre or pre[-1] == 0:
        continue
    B = len(pre)
    W = pre[-1] // B + 1
    J = []
    k = 0
    for j in range(B):
        target = j * W
        while k < B - 1 and pre[k] <= target:
            k += 1
        J.append(k)
    JUMP[key] = (W, J)
jt_build = time.perf_counter() - t0

OPS = 0
PROBES = 0


def unrank_jump(r, n, m):
    global OPS, PROBES
    OPS += 1
    if n - 1 <= m:
        if r == 0:
            return ('var', n - 1)
        r -= 1
    lam_count = T(n - 2, m + 1)
    OPS += 1
    if r < lam_count:
        return ('lam', unrank_jump(r, n - 2, m + 1))
    r -= lam_count
    key = (n, min(m, n - 1))
    pre = SPLIT[key]
    W, J = JUMP[key]
    k = J[min(r // W, len(J) - 1)]
    while pre[k] <= r:
        k += 1
        PROBES += 1
    OPS += 2
    PROBES += 1
    base = pre[k - 1] if k > 0 else 0
    right = T(n - 2 - (k + 2), m)
    a, b = divmod(r - base, right)
    OPS += 1
    return ('app', unrank_jump(a, k + 2, m), unrank_jump(b, n - 2 - (k + 2), m))


def unrank_bisect(r, n, m):
    global OPS
    OPS += 1
    if n - 1 <= m:
        if r == 0:
            return ('var', n - 1)
        r -= 1
    lam_count = T(n - 2, m + 1)
    OPS += 1
    if r < lam_count:
        return ('lam', unrank_bisect(r, n - 2, m + 1))
    r -= lam_count
    pre = SPLIT[(n, min(m, n - 1))]
    i = bisect_right(pre, r)
    OPS += max(1, len(pre).bit_length())
    base = pre[i - 1] if i > 0 else 0
    right = T(n - 2 - (i + 2), m)
    a, b = divmod(r - base, right)
    OPS += 1
    return ('app', unrank_bisect(a, i + 2, m), unrank_bisect(b, n - 2 - (i + 2), m))


def decode_with(unranker, s):
    N = int('1' + s, 2) - 1
    n = bisect_right(CUM, N, 4, n_max + 1)
    return unranker(N - CUM[n - 1], n, 0)


report("")
report(f"2. JUMP TABLES vs BINARY SEARCH in num->lam (cap {L1_MAX} bits)")
report(f"   jump-table build: {1000 * jt_build:.0f} ms, "
       f"{sum(len(j) for _, j in JUMP.values()):,} extra words "
       f"(rows append-only under cap extension, same as T)")
report(f"   {'L':>4} {'ops bisect':>11} {'ops jump':>9} {'probes/app node':>16}")
for L in (10, 20, 30, 40, 50, 60):
    Ns = [(1 << L) - 1 + (j * ((1 << L) // 200 + 1)) % (1 << L) for j in range(200)]
    ss = [bin(N + 1)[3:] for N in Ns]
    OPS = 0
    terms_b = [decode_with(unrank_bisect, s) for s in ss]
    ops_b = OPS / len(ss)
    OPS = 0
    PROBES = 0
    terms_j = [decode_with(unrank_jump, s) for s in ss]
    ops_j = OPS / len(ss)
    assert terms_b == terms_j
    apps = sum(1 for t in terms_j for _ in [0] if t) or 1
    n_apps = sum(str(t).count("'app'") for t in terms_j)
    report(f"   {L:>4} {ops_b:>11.1f} {ops_j:>9.1f} "
           f"{PROBES / max(1, n_apps):>16.2f}")

print("\n".join(OUT))
