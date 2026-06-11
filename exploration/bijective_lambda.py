#!/usr/bin/env python3
"""Bijection between binary strings {0,1}* and closed lambda terms.

Construction: order closed terms by (size, canonical rank within size),
then identify {0,1}* with N via bijective binary numeration
    s  <->  N = int('1' + s, 2) - 1        ("" -> 0, "0" -> 1, "1" -> 2, ...)
so every binary string of every length (all 2^L strings per length L,
including the empty string) denotes exactly one closed lambda term, and
every closed lambda term has exactly one binary representation.

Size measure = Tromp's binary lambda calculus (BLC) bit size:
    |var i|   = i + 1        (de Bruijn index i >= 1, code 1^i 0)
    |lam M|   = |M| + 2      (code 00 <M>)
    |app M N| = |M| + |N| + 2 (code 01 <M> <N>)

T(n, m) = number of terms of size n whose free de Bruijn indices are all <= m.
Closed terms of size n = T(n, 0).  This recurrence is what makes the bijection
polynomial-time: decoding walks the grammar, using precomputed counts to skip
whole subclasses at once, instead of enumerating terms one by one.

Terms are tuples: ('var', i), ('lam', body), ('app', fun, arg).
"""

import sys
import time

sys.setrecursionlimit(100000)


# ---------------------------------------------------------------- counting

# TT[n] is a list indexed by m in [0, n-1]; a term of size n cannot contain a
# variable with index > n-1, so T(n, m) saturates at m = n-1 and lookups clamp.
TT = [[], []]
CUM = [0, 0]  # CUM[n] = number of closed terms of size <= n


def T(n, m):
    if n < 2:
        return 0
    if m >= n:
        m = n - 1
    return TT[n][m]


def extend_tables(n_max):
    for n in range(len(TT), n_max + 1):
        row = []
        for m in range(n):
            v = 1 if m == n - 1 else 0          # var (n-1), needs index <= m
            v += T(n - 2, m + 1)                # lam: body in extended context
            for k in range(2, n - 3):           # app: left size k, right n-2-k
                v += T(k, m) * T(n - 2 - k, m)
            row.append(v)
        TT.append(row)
        CUM.append(CUM[-1] + (row[0] if row else 0))


# ------------------------------------------------------- rank / unrank

def size(t):
    if t[0] == 'var':
        return t[1] + 1
    if t[0] == 'lam':
        return size(t[1]) + 2
    return size(t[1]) + size(t[2]) + 2


def rank(t, m=0):
    """Rank of t within the class (size(t), m). Order: var < lam < app,
    apps by (left size asc, left rank, right rank)."""
    if t[0] == 'var':
        return 0
    n = size(t)
    r = 1 if n - 1 <= m else 0
    if t[0] == 'lam':
        return r + rank(t[1], m + 1)
    r += T(n - 2, m + 1)
    f, a = t[1], t[2]
    fs, as_ = size(f), size(a)
    for k in range(2, fs):
        r += T(k, m) * T(n - 2 - k, m)
    return r + rank(f, m) * T(as_, m) + rank(a, m)


def unrank(r, n, m=0):
    """Inverse of rank: the r-th term of size n in context m."""
    if n - 1 <= m:
        if r == 0:
            return ('var', n - 1)
        r -= 1
    lam_count = T(n - 2, m + 1)
    if r < lam_count:
        return ('lam', unrank(r, n - 2, m + 1))
    r -= lam_count
    for k in range(2, n - 3):
        right = T(n - 2 - k, m)
        block = T(k, m) * right
        if r < block:
            i, j = divmod(r, right)
            return ('app', unrank(i, k, m), unrank(j, n - 2 - k, m))
        r -= block
    raise ValueError(f"rank {r} out of range for (n={n}, m={m})")


# ------------------------------------------- the bijection {0,1}* <-> terms

def encode(t):
    """Closed lambda term -> binary string."""
    n = size(t)
    extend_tables(n)
    N = CUM[n - 1] + rank(t, 0)
    return bin(N + 1)[3:]


def decode(s):
    """Binary string -> closed lambda term."""
    N = int('1' + s, 2) - 1
    n = 4
    while True:
        extend_tables(n)
        if CUM[n] > N:
            break
        n += 1
    return unrank(N - CUM[n - 1], n, 0)


# --------------------------------------------------- Tromp BLC, for compare

def blc(t):
    if t[0] == 'var':
        return '1' * t[1] + '0'
    if t[0] == 'lam':
        return '00' + blc(t[1])
    return '01' + blc(t[1]) + blc(t[2])


def blc_parse(s, i=0):
    """Parse one term starting at i; returns (term, next_index) or None."""
    if i >= len(s):
        return None
    if s[i] == '0':
        if i + 1 >= len(s):
            return None
        if s[i + 1] == '0':
            r = blc_parse(s, i + 2)
            if r is None:
                return None
            t, j = r
            return ('lam', t), j
        r = blc_parse(s, i + 2)
        if r is None:
            return None
        f, j = r
        r = blc_parse(s, j)
        if r is None:
            return None
        a, k = r
        return ('app', f, a), k
    j = i
    while j < len(s) and s[j] == '1':
        j += 1
    if j >= len(s):
        return None
    return ('var', j - i), j + 1


def max_free(t, depth=0):
    if t[0] == 'var':
        return t[1] - depth
    if t[0] == 'lam':
        return max_free(t[1], depth + 1)
    return max(max_free(t[1], depth), max_free(t[2], depth))


def blc_decode_closed(s):
    """Full-string BLC decode to a closed term, or None."""
    r = blc_parse(s)
    if r is None:
        return None
    t, j = r
    if j != len(s) or max_free(t) > 0:
        return None
    return t


# ------------------------------------------------------------ pretty print

def show(t, ctx=0):
    if t[0] == 'var':
        return str(t[1])
    if t[0] == 'lam':
        s = 'λ' + show(t[1])
        return '(' + s + ')' if ctx > 0 else s
    s = show(t[1], 1) + ' ' + show(t[2], 2)
    return '(' + s + ')' if ctx == 2 else s


# ------------------------------------------------------------------- report

def main():
    out = []

    # --- 1. validate the counting recurrence against brute force over BLC ---
    out.append("1. Counting recurrence vs brute force (all 2^n strings, BLC-decoded):")
    out.append(f"   {'n':>3} {'2^n':>8} {'valid closed (brute)':>21} {'T(n,0)':>8} {'density':>10}")
    terms_by_size = {}
    extend_tables(20)
    for n in range(4, 18):
        found = []
        for bits in range(1 << n):
            s = format(bits, f'0{n}b')
            t = blc_decode_closed(s)
            if t is not None:
                found.append(t)
        terms_by_size[n] = found
        ok = "✓" if len(found) == T(n, 0) else "✗ MISMATCH"
        out.append(f"   {n:>3} {1 << n:>8} {len(found):>21} {T(n, 0):>8} "
                   f"{len(found) / (1 << n):>10.6f} {ok}")

    # --- 2. the bijection on the first 16 strings ---
    out.append("")
    out.append("2. First 16 binary strings and their lambda terms:")
    for N in range(16):
        s = bin(N + 1)[3:]
        t = decode(s)
        out.append(f"   {('ε' if s == '' else s):>5}  ->  {show(t)}")

    # --- 3. round-trip checks, both directions ---
    out.append("")
    fails = 0
    for N in range(50000):
        s = bin(N + 1)[3:]
        t = decode(s)
        if encode(t) != s:
            fails += 1
    out.append(f"3. Round trip string -> term -> string: first 50000 strings "
               f"(lengths 0..{len(bin(50000)) - 3}), {fails} failures")
    fails = 0
    total = 0
    for n, terms in terms_by_size.items():
        for t in terms:
            total += 1
            if decode(encode(t)) != t:
                fails += 1
    out.append(f"   Round trip term -> string -> term: all {total} closed terms "
               f"of BLC size <= 17, {fails} failures")

    # --- 4. code length: bijective vs BLC ---
    out.append("")
    out.append("4. Code length, bijective vs Tromp BLC:")
    omega_half = ('lam', ('app', ('var', 1), ('var', 1)))
    omega = ('app', omega_half, omega_half)
    y_half = ('lam', ('app', ('var', 2), ('app', ('var', 1), ('var', 1))))
    ycomb = ('lam', ('app', y_half, y_half))
    s2 = ('lam', ('lam', ('lam', ('app', ('app', ('var', 3), ('var', 1)),
                                  ('app', ('var', 2), ('var', 1))))))  # S
    for name, t in [("λx.x", ('lam', ('var', 1))),
                    ("S = λλλ(3 1 (2 1))", s2),
                    ("Ω = (λx.x x)(λx.x x)", omega),
                    ("Y = λf.(λx.f(x x))(λx.f(x x))", ycomb)]:
        b = encode(t)
        out.append(f"   {name:<32} BLC {len(blc(t)):>3} bits   "
                   f"bijective {len(b):>3} bits   ({b if len(b) <= 32 else b[:29] + '...'})")

    # --- 5. asymptotic growth rate of closed-term counts ---
    extend_tables(220)
    ratios = [TT[n][0] / TT[n - 1][0] for n in range(216, 221)]
    out.append("")
    out.append(f"5. Growth of closed-term counts T(n,0): ratios at n=216..220: "
               + ", ".join(f"{r:.6f}" for r in ratios))
    import math
    rho = ratios[-1]
    out.append(f"   => T(n,0) ~ rho^n with rho ≈ {rho:.4f}; "
               f"bijective code length ≈ {math.log2(rho):.4f} * BLC length "
               f"(saves ≈ {100 * (1 - math.log2(rho)):.2f}% asymptotically)")

    # --- 6. decoding a 202-bit string: poly(L), not O(2^L) ---
    s_big = bin(3 ** 127)[2:]
    t0 = time.perf_counter()
    t_big = decode(s_big)
    t1 = time.perf_counter()
    assert encode(t_big) == s_big
    t2 = time.perf_counter()
    N_big = int('1' + s_big, 2) - 1
    out.append("")
    out.append(f"6. Decode a {len(s_big)}-bit string (index N ≈ {float(N_big):.3e}):")
    out.append(f"   decoded to closed term of BLC size {size(t_big)} "
               f"in {t1 - t0:.3f} s (re-encode round trip ✓, {t2 - t1:.3f} s)")
    out.append(f"   naive 'enumerate up to N' would need ~{float(N_big):.0e} steps.")
    out.append(f"   term: {show(t_big)[:120]}...")

    print("\n".join(out))


if __name__ == '__main__':
    main()
