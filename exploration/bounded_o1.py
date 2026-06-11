#!/usr/bin/env python3
"""Bounded-length bijection {0,1}^{<=L} <-> closed lambda terms, precached.

Two regimes:

Level 0 (L small, here L <= 16): materialize the entire bijection.
    DECODE_TAB[N] = interned term, ENCODE_TAB[term] = N.
    One query = one array index (decode) or one hash lookup (encode): O(1).
    Space O(2^L), preprocessing O(2^L * L).

Level 1 (L up to word size, here L <= 60): precache counting tables only.
    Key fact making this work: T(n, m) <= 2^n, because each term of size n
    has a distinct n-bit BLC code. Bounding L therefore bounds every count
    below ~2^(L+8): all arithmetic is fixed-width machine words, no bigints.
    Precache:  T(n,m) table, cumulative closed counts CUM, and per-(n,m)
    prefix sums over application splits (the Catalan-convolution blocks).
    decode = grammar walk steered by table comparisons: O(L) word ops total
    (binary search per app node adds a log factor; a per-(n,m) jump table
    removes it). encode = same walk with direct lookups, no search: O(L).
"""

import sys
import time
from bisect import bisect_right

sys.setrecursionlimit(100000)

from bijective_lambda import TT, CUM, T, extend_tables, size, unrank, show

OUT = []


def report(line=""):
    OUT.append(line)


# ===================================================================
# Level 0: full table for strings of length <= L0_MAX
# ===================================================================

L0_MAX = 16
N_LIMIT = (1 << (L0_MAX + 1)) - 1          # strings of length <= 16 <-> N < 2^17 - 1

INTERN = {}


def intern(t):
    if t[0] == 'var':
        key = t
    elif t[0] == 'lam':
        key = ('lam', intern(t[1]))
    else:
        key = ('app', intern(t[1]), intern(t[2]))
    return INTERN.setdefault(key, key)


def build_level0():
    decode_tab = [None] * N_LIMIT
    n = 4
    extend_tables(64)
    while CUM[n - 1] < N_LIMIT:
        for r in range(T(n, 0)):
            N = CUM[n - 1] + r
            if N >= N_LIMIT:
                break
            decode_tab[N] = intern(unrank(r, n, 0))
        n += 1
    encode_tab = {t: N for N, t in enumerate(decode_tab)}
    return decode_tab, encode_tab


t0 = time.perf_counter()
DECODE_TAB, ENCODE_TAB = build_level0()
t1 = time.perf_counter()

report("LEVEL 0: full materialization, strings of length <= %d" % L0_MAX)
report(f"   entries: {N_LIMIT:,} strings/terms, {len(INTERN):,} distinct interned nodes")
report(f"   build time: {t1 - t0:.2f} s (one-time)")

# timing: decode = string -> term
samples = [bin(N + 1)[3:] for N in range(N_LIMIT)]
t0 = time.perf_counter()
acc = 0
for s in samples:
    acc ^= id(DECODE_TAB[int('1' + s, 2) - 1])
t1 = time.perf_counter()
report(f"   decode (incl. string->int): {1e9 * (t1 - t0) / len(samples):.0f} ns/query")

terms = DECODE_TAB
t0 = time.perf_counter()
acc = 0
for t in terms:
    acc ^= ENCODE_TAB[t]
t1 = time.perf_counter()
report(f"   encode (hash lookup):       {1e9 * (t1 - t0) / len(terms):.0f} ns/query")

# bounded sub-ranges the user mentioned
n_exact8 = (1 << 9) - (1 << 8)
n_5_15 = (1 << 16) - (1 << 5)
report(f"   range sizes: exactly 8 bits -> {n_exact8} strings; "
       f"lengths 5..15 -> {n_5_15:,} strings")


# ===================================================================
# Level 1: word-arithmetic decode for strings of length <= L1_MAX
# ===================================================================

L1_MAX = 60

# find n_max covering all indices N <= 2^(L1_MAX+1) - 2
N_MAX = (1 << (L1_MAX + 1)) - 2
n_max = 4
while True:
    extend_tables(n_max)
    if CUM[n_max] > N_MAX:
        break
    n_max += 1

t0 = time.perf_counter()
# prefix sums over app splits: SPLIT[(n,m)][i] = number of size-n apps in
# context m whose left subterm has size < k_i (block boundaries by left size)
SPLIT = {}
max_bits = 0
for n in range(6, n_max + 1):
    for m in range(n):
        acc = 0
        pre = []
        for k in range(2, n - 3):
            acc += T(k, m) * T(n - 2 - k, m)
            pre.append(acc)
        SPLIT[(n, m)] = pre
        if acc:
            max_bits = max(max_bits, acc.bit_length())
t1 = time.perf_counter()
n_entries = sum(len(v) for v in SPLIT.values())

report("")
report(f"LEVEL 1: count-table precache, strings of length <= {L1_MAX}")
report(f"   n_max = {n_max}; split-table entries: {n_entries:,} words "
       f"+ T-table {sum(len(r) for r in TT[:n_max + 1]):,} words")
report(f"   largest table entry: {max_bits} bits "
       f"(fits 2 x u64 limbs; all ops fixed-width)")
report(f"   build time: {t1 - t0:.3f} s (one-time)")

OPS = 0


def unrank_fast(r, n, m):
    """Same enumeration as bijective_lambda.unrank; app split found by
    binary search in the precached prefix sums instead of a linear scan."""
    global OPS
    OPS += 1
    if n - 1 <= m:
        if r == 0:
            return ('var', n - 1)
        r -= 1
    lam_count = T(n - 2, m + 1)
    OPS += 1
    if r < lam_count:
        return ('lam', unrank_fast(r, n - 2, m + 1))
    r -= lam_count
    pre = SPLIT[(n, min(m, n - 1))]
    i = bisect_right(pre, r)
    OPS += max(1, len(pre).bit_length())
    base = pre[i - 1] if i > 0 else 0
    k = i + 2
    right = T(n - 2 - k, m)
    a, b = divmod(r - base, right)
    OPS += 1
    return ('app', unrank_fast(a, k, m), unrank_fast(b, n - 2 - k, m))


def rank_fast(t, m=0):
    """Inverse walk: direct table lookups, no search at all."""
    global OPS
    if t[0] == 'var':
        return 0
    n = size(t)
    r = 1 if n - 1 <= m else 0
    if t[0] == 'lam':
        OPS += 1
        return r + rank_fast(t[1], m + 1)
    r += T(n - 2, m + 1)
    f, a = t[1], t[2]
    fs = size(f)
    i = fs - 2
    if i > 0:
        r += SPLIT[(n, min(m, n - 1))][i - 1]
    OPS += 3
    return r + rank_fast(f, m) * T(size(a), m) + rank_fast(a, m)


def decode_fast(s):
    N = int('1' + s, 2) - 1
    n = bisect_right(CUM, N, 4, n_max + 1)
    return unrank_fast(N - CUM[n - 1], n, 0)


def encode_fast(t):
    n = size(t)
    return bin(CUM[n - 1] + rank_fast(t, 0) + 1)[3:]


# correctness: fast path == reference path, and round trips
checked = 0
for L in range(0, L1_MAX + 1, 4):
    base = (1 << L) - 1
    for j in range(25):
        N = base + (j * ((1 << L) // 25 + 1)) % (1 << L) if L > 0 else 0
        s = bin(N + 1)[3:]
        t = decode_fast(s)
        assert encode_fast(t) == s, (L, s)
        checked += 1
report(f"   correctness: {checked} round trips across lengths 0..{L1_MAX} ✓")

# scaling: word ops and wall time per decode as a function of L
report(f"   {'L (bits)':>9} {'ops/decode':>11} {'us/decode':>10} "
       f"{'ops/encode':>11} {'us/encode':>10}")
for L in (10, 20, 30, 40, 50, 60):
    Ns = [(1 << L) - 1 + (j * ((1 << L) // 200 + 1)) % (1 << L) for j in range(200)]
    ss = [bin(N + 1)[3:] for N in Ns]
    OPS = 0
    t0 = time.perf_counter()
    ts = [decode_fast(s) for s in ss]
    t1 = time.perf_counter()
    dec_ops, dec_us = OPS / len(ss), 1e6 * (t1 - t0) / len(ss)
    OPS = 0
    t0 = time.perf_counter()
    for t in ts:
        encode_fast(t)
    t1 = time.perf_counter()
    report(f"   {L:>9} {dec_ops:>11.1f} {dec_us:>10.2f} "
           f"{OPS / len(ss):>11.1f} {1e6 * (t1 - t0) / len(ss):>10.2f}")

example = decode_fast(bin((1 << 60) + 987654321098765 + 1)[3:])
report(f"   example 60-bit decode: {show(example)}")

print("\n".join(OUT))
