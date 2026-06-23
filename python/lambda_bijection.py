#!/usr/bin/env python3
"""Binary2Lambda: a total bijection between binary strings and closed untyped
lambda terms.

Every binary string (including the empty string) denotes exactly one closed
lambda term, and every closed lambda term has exactly one binary string.
Both directions share one incrementally-built counting table.

Specification
-------------
Terms are de Bruijn-indexed:  Var(i) with i >= 1,  Lam(body),  App(fun, arg).

Size measure (bits of Tromp's binary lambda calculus):
    |Var(i)| = i + 1      |Lam(b)| = |b| + 2      |App(f, a)| = |f| + |a| + 2

Canonical enumeration of closed terms: ascending size, and within one size
class the order  Var < Lam < App,  abstractions ordered by body rank,
applications ordered by (left subterm size, left rank, right rank).

A string s corresponds to the natural number N = int('1' + s, 2) - 1
(bijective binary numeration: the leading 1 is implicit), and N indexes the
canonical enumeration.

De Bruijn index cap
-------------------
A Table built with index_cap=K enumerates the sublanguage of closed terms
whose de Bruijn indices never exceed K.  Each value of K (including None,
meaning unbounded) defines a DIFFERENT bijection; encode and decode must use
the same cap.  The capped and unbounded bijections agree on all terms of
size <= K + 1.  A finite cap shrinks the counting table from Theta(n^2) to
Theta(K n) entries, which is what makes very long strings affordable.

Incrementality
--------------
Both table axes grow without recomputation where mathematically possible:
size extension appends rows and touches nothing; changing the cap reuses
every row of size <= min(old, new) + 1 (those are cap-independent) and
rebuilds only the rest.

Recursion depth and untrusted input
-----------------------------------
encode/decode and the structural helpers recurse to a depth proportional to
a term's nesting (Theta(size) for a degenerate spine), so the interpreter
recursion limit is raised below.  A decoded term of size n comes from a
string of about n bits, so the default raise comfortably covers multi-kilobit
strings; for adversarial, deeply degenerate inputs far beyond that, raise the
limit further or reimplement the walks iteratively.  decode also grows the
counting table to a size proportional to the input length (Theta(L^2) bignum
entries unbounded, Theta(K*L) capped); a caller handling untrusted input of
unbounded length should bound the input length itself.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import Optional, Union

# Terms recurse to a depth proportional to their nesting; the default limit of
# 1000 overflows around 1 Kbit degenerate strings.  This covers ~14 Kbit.
sys.setrecursionlimit(max(sys.getrecursionlimit(), 15000))

# --------------------------------------------------------------------- terms


@dataclass(frozen=True)
class Var:
    index: int


@dataclass(frozen=True)
class Lam:
    body: "Term"


@dataclass(frozen=True)
class App:
    fun: "Term"
    arg: "Term"


Term = Union[Var, Lam, App]


def term_size(term: Term) -> int:
    if isinstance(term, Var):
        return term.index + 1
    if isinstance(term, Lam):
        return term_size(term.body) + 2
    return term_size(term.fun) + term_size(term.arg) + 2


def max_de_bruijn_index(term: Term) -> int:
    """Largest index occurring anywhere in the term (0 if none occur)."""
    if isinstance(term, Var):
        return term.index
    if isinstance(term, Lam):
        return max_de_bruijn_index(term.body)
    return max(max_de_bruijn_index(term.fun), max_de_bruijn_index(term.arg))


def check_closed(term: Term, depth: int = 0) -> None:
    """Raise ValueError unless term is a well-formed closed term: every de
    Bruijn index is >= 1 and bound by an enclosing lambda (index <= depth).

    encode requires this; without it a non-closed term would silently encode
    to some other closed term's string, breaking the bijection."""
    if isinstance(term, Var):
        if term.index < 1:
            raise ValueError(f"de Bruijn index must be >= 1, got {term.index}")
        if term.index > depth:
            raise ValueError("term is not closed (has a free variable)")
    elif isinstance(term, Lam):
        check_closed(term.body, depth + 1)
    else:
        check_closed(term.fun, depth)
        check_closed(term.arg, depth)


def show_term(term: Term, context: int = 0) -> str:
    """Render a term: λ-bodies extend right, application is left-associative.

    context: 0 = top level / lambda body, 1 = left of application,
    2 = right of application.
    """
    if isinstance(term, Var):
        return str(term.index)
    if isinstance(term, Lam):
        rendered = "λ" + show_term(term.body, 0)
        return f"({rendered})" if context > 0 else rendered
    rendered = show_term(term.fun, 1) + " " + show_term(term.arg, 2)
    return f"({rendered})" if context == 2 else rendered


def show_bits(bits: str) -> str:
    return bits if bits else "ε"


_FNV_OFFSET = 0xCBF29CE484222325
_FNV_PRIME = 0x100000001B3
_FNV_MASK = 0xFFFFFFFFFFFFFFFF

# Magic for the portable binary table format (.lamtab); the last byte is the
# format version. Shared byte-for-byte with the C++, Rust and Wolfram ports.
_TABLE_MAGIC = b"LAMTAB\x01"


def _fnv1a64_update(h: int, data: bytes) -> int:
    """Fold data into a running FNV-1a 64-bit hash, used as a table-file
    integrity checksum. Identical integer arithmetic in the C++ and Rust
    implementations, so a table file's checksum is the same in every language.
    Folding incrementally lets save and load stream the file without holding
    the whole body in memory."""
    for byte in data:
        h = ((h ^ byte) * _FNV_PRIME) & _FNV_MASK
    return h


# --------------------------------------------------------------------- table


class Table:
    """Counting table T(n, m): closed-context term counts, both axes growable.

    T(n, m) = number of terms of size n whose free indices are all <= m and
    whose indices nowhere exceed the cap.  Stored rows are normalised: the
    effective context is min(m, n - 1, cap), since a size-n term cannot
    contain an index above n - 1 and the cap masks everything beyond it.
    """

    def __init__(self, index_cap: Optional[int] = None):
        if index_cap is not None and index_cap < 1:
            raise ValueError("index_cap must be at least 1 (or None)")
        self.index_cap = index_cap
        self._rows: list[list[int]] = [[], []]  # rows for n = 0, 1 are empty
        self._cum: list[int] = [0, 0]           # closed terms of size <= n

    # -- public interface ---------------------------------------------------

    @property
    def built_size(self) -> int:
        return len(self._rows) - 1

    def extend(self, size_limit: int) -> None:
        """Append rows up to size_limit; existing rows are never touched."""
        for n in range(len(self._rows), size_limit + 1):
            row = [self._row_value(n, m) for m in range(self._width(n))]
            self._rows.append(row)
            self._cum.append(self._cum[-1] + (row[0] if row else 0))

    def set_index_cap(self, new_cap: Optional[int]) -> None:
        """Change the de Bruijn cap, reusing all cap-independent rows.

        Rows of size n are identical under caps K and K' whenever
        n - 1 <= min(K, K'); those are kept, the rest are rebuilt.
        """
        if new_cap is not None and new_cap < 1:
            raise ValueError("index_cap must be at least 1 (or None)")
        if new_cap == self.index_cap:
            return
        finite = [c for c in (self.index_cap, new_cap) if c is not None]
        keep = min(finite) + 1
        target = self.built_size
        self.index_cap = new_cap
        del self._rows[keep + 1:]
        del self._cum[keep + 1:]
        self.extend(target)

    def save(self, path: str) -> None:
        """Write the table to disk (portable binary format, shared with C++ and
        Rust; the bytes are identical in every language).

        Layout: the 7-byte magic `LAMTAB\\x01` (the last byte is the version),
        a one-byte cap kind (0 unbounded, 1 finite) and 4-byte little-endian
        cap value, the 4-byte little-endian built size, then for each size from
        2 upward its row of counts, each count a 4-byte little-endian length
        and that many big-endian magnitude bytes, then an 8-byte little-endian
        FNV-1a-64 of every preceding byte. Cumulative counts are derivable, so
        they are not stored. The body is streamed, so memory stays flat."""
        checksum = _FNV_OFFSET

        def emit(handle, chunk: bytes) -> None:
            nonlocal checksum
            handle.write(chunk)
            checksum = _fnv1a64_update(checksum, chunk)

        cap = self.index_cap
        with open(path, "wb") as handle:
            emit(handle, _TABLE_MAGIC)
            emit(handle, b"\x00" if cap is None else b"\x01")
            emit(handle, (0 if cap is None else cap).to_bytes(4, "little"))
            emit(handle, self.built_size.to_bytes(4, "little"))
            for n in range(2, self.built_size + 1):
                for value in self._rows[n]:
                    width = (value.bit_length() + 7) // 8
                    emit(handle, width.to_bytes(4, "little"))
                    emit(handle, value.to_bytes(width, "big"))
            handle.write(checksum.to_bytes(8, "little"))

    @classmethod
    def load(cls, path: str) -> "Table":
        """Read a table written by save (any implementation).  The checksum is
        verified, so any corruption or truncation is rejected before the table
        is returned."""
        with open(path, "rb") as handle:
            data = handle.read()
        head = len(_TABLE_MAGIC) + 1 + 4 + 4
        if len(data) < head + 8 or data[:len(_TABLE_MAGIC)] != _TABLE_MAGIC:
            raise ValueError("not a lambda-binarization table file")
        if _fnv1a64_update(_FNV_OFFSET, memoryview(data)[:-8]) \
                != int.from_bytes(data[-8:], "little"):
            raise ValueError("table file failed checksum (corrupt or truncated)")
        pos = len(_TABLE_MAGIC)
        cap_kind = data[pos]
        cap_value = int.from_bytes(data[pos + 1:pos + 5], "little")
        built = int.from_bytes(data[pos + 5:pos + 9], "little")
        pos = head
        if cap_kind not in (0, 1) or built < 1:
            raise ValueError("table file has a malformed header")
        table = cls(index_cap=None if cap_kind == 0 else cap_value)
        end = len(data) - 8
        for n in range(2, built + 1):
            row = []
            for _ in range(table._width(n)):
                if pos + 4 > end:
                    raise ValueError("table file is truncated")
                width = int.from_bytes(data[pos:pos + 4], "little")
                pos += 4
                if pos + width > end:
                    raise ValueError("table file is truncated")
                row.append(int.from_bytes(data[pos:pos + width], "big"))
                pos += width
            table._rows.append(row)
            table._cum.append(table._cum[-1] + (row[0] if row else 0))
        if pos != end:
            raise ValueError("table file size does not match its rows")
        return table

    def entry_count(self) -> int:
        """Total number of stored count entries (a table-size metric)."""
        return sum(len(self._rows[n]) for n in range(2, self.built_size + 1))

    def total_bit_length(self) -> int:
        """Total bit length of all stored counts (an in-memory size metric)."""
        return sum(v.bit_length()
                   for n in range(2, self.built_size + 1)
                   for v in self._rows[n])

    def count(self, n: int, m: int) -> int:
        """T(n, m); extends the table if size n is not yet built."""
        if n < 2:
            return 0
        if n > self.built_size:
            self.extend(n)
        return self._rows[n][self._effective_context(n, m)]

    def closed_cumulative(self, n: int) -> int:
        """Number of closed terms of size <= n."""
        if n < 2:
            return 0
        if n > self.built_size:
            self.extend(n)
        return self._cum[n]

    # -- internals ----------------------------------------------------------

    def _effective_context(self, n: int, m: int) -> int:
        m = min(m, n - 1)
        if self.index_cap is not None:
            m = min(m, self.index_cap)
        return m

    def _width(self, n: int) -> int:
        return self._effective_context(n, n - 1) + 1

    def _row_value(self, n: int, m: int) -> int:
        """T(n, m) computed from already-built smaller rows; m is effective."""
        value = 1 if m == n - 1 else 0                      # the variable n-1
        value += self._built_count(n - 2, m + 1)            # abstractions
        for k in range(2, n - 3):                           # applications
            value += self._built_count(k, m) * self._built_count(n - 2 - k, m)
        return value

    def _built_count(self, n: int, m: int) -> int:
        if n < 2:
            return 0
        return self._rows[n][self._effective_context(n, m)]


# --------------------------------------------------------------------- codec


def encode(table: Table, term: Term) -> str:
    """Closed lambda term -> binary string (inverse of decode).

    Raises ValueError if term is not closed or uses an index above the cap;
    encode is only defined on the closed terms the table enumerates."""
    check_closed(term)
    cap = table.index_cap
    if cap is not None and max_de_bruijn_index(term) > cap:
        raise ValueError(f"term uses indices above the table cap {cap}")
    n = term_size(term)
    number = table.closed_cumulative(n - 1) + _rank(table, term, 0)
    return bin(number + 1)[3:]


def decode(table: Table, bits: str) -> Term:
    """Binary string -> closed lambda term (total: never fails on 0/1 input).

    The table grows as needed; for untrusted input see the module docstring."""
    if any(c not in "01" for c in bits):
        raise ValueError("input must consist of 0s and 1s")
    number = int("1" + bits, 2) - 1
    n = 4  # the smallest closed term, λ1, has size 4
    while table.closed_cumulative(n) <= number:
        n += 1
    return _unrank(table, number - table.closed_cumulative(n - 1), n, 0)


def decode_index(table: Table, index: int) -> Term:
    """Closed term at enumeration index `index` (`index >= 0`): the integer
    view of the bijection. `decode_index(0)` is the smallest term, λ1, and
    `decode_index(N)` equals `decode(table, b)` for `b` the bijective-binary
    form of N. Use this when the data are integers (fixed-width genomes,
    Gödel-style IDs, uniform sampling) rather than variable-length bit
    strings."""
    if index < 0:
        raise ValueError("index must be >= 0")
    return decode(table, bin(index + 1)[3:])


def encode_index(table: Table, term: Term) -> int:
    """Enumeration index of a closed term (inverse of decode_index)."""
    return int("1" + encode(table, term), 2) - 1


def _var_count(table: Table, n: int, m: int) -> int:
    if n < 2 or n - 1 > m:
        return 0
    if table.index_cap is not None and n - 1 > table.index_cap:
        return 0
    return 1


def _rank(table: Table, term: Term, m: int) -> int:
    if isinstance(term, Var):
        return 0
    n = term_size(term)
    rank = _var_count(table, n, m)
    if isinstance(term, Lam):
        return rank + _rank(table, term.body, m + 1)
    rank += table.count(n - 2, m + 1)
    left_size = term_size(term.fun)
    for k in range(2, left_size):
        rank += table.count(k, m) * table.count(n - 2 - k, m)
    return (rank
            + _rank(table, term.fun, m) * table.count(term_size(term.arg), m)
            + _rank(table, term.arg, m))


def _unrank(table: Table, rank: int, n: int, m: int) -> Term:
    # Precondition: 0 <= rank < T(n, m).  The divmod split below is exact only
    # because rank < count(k,m)*right, so the quotient stays below count(k,m).
    if _var_count(table, n, m):
        if rank == 0:
            return Var(n - 1)
        rank -= 1
    lam_block = table.count(n - 2, m + 1)
    if rank < lam_block:
        return Lam(_unrank(table, rank, n - 2, m + 1))
    rank -= lam_block
    for k in range(2, n - 3):
        right = table.count(n - 2 - k, m)
        block = table.count(k, m) * right
        if rank < block:
            left_rank, right_rank = divmod(rank, right)
            return App(_unrank(table, left_rank, k, m),
                       _unrank(table, right_rank, n - 2 - k, m))
        rank -= block
    raise AssertionError("rank out of range for size class")


# ----------------------------------------------------------------- self-test


def _blc_parse(bits: str, pos: int = 0):
    """Parse one BLC-coded term; returns (term, next_pos) or None."""
    if pos >= len(bits):
        return None
    if bits[pos] == "0":
        if pos + 1 >= len(bits):
            return None
        tag, rest = bits[pos + 1], pos + 2
        first = _blc_parse(bits, rest)
        if first is None:
            return None
        if tag == "0":
            return Lam(first[0]), first[1]
        second = _blc_parse(bits, first[1])
        if second is None:
            return None
        return App(first[0], second[0]), second[1]
    end = pos
    while end < len(bits) and bits[end] == "1":
        end += 1
    if end >= len(bits):
        return None
    return Var(end - pos), end + 1


def _all_blc_terms(n: int):
    """Every BLC-parseable term of size exactly n (closed or open)."""
    terms = []
    for value in range(1 << n):
        parsed = _blc_parse(format(value, f"0{n}b"))
        if parsed is not None and parsed[1] == n:
            terms.append(parsed[0])
    return terms


def _brute_force_terms(n: int, cap: Optional[int]):
    """All closed terms of size n (indices <= cap) by exhaustive BLC decode."""
    return [term for term in _all_blc_terms(n)
            if _max_free(term) <= 0
            and (cap is None or max_de_bruijn_index(term) <= cap)]


def _max_free(term: Term, depth: int = 0) -> int:
    if isinstance(term, Var):
        return term.index - depth
    if isinstance(term, Lam):
        return _max_free(term.body, depth + 1)
    return max(_max_free(term.fun, depth), _max_free(term.arg, depth))


def _self_test() -> None:
    for cap in (None, 1, 2, 5):
        table = Table(index_cap=cap)
        for n in range(4, 15):  # counts match exhaustive enumeration
            assert table.count(n, 0) == len(_brute_force_terms(n, cap)), (cap, n)
        for number in range(3000):  # string -> term -> string
            bits = bin(number + 1)[3:]
            assert encode(table, decode(table, bits)) == bits, (cap, number)
    idx = Table()  # integer view: index -> term -> index
    for number in range(300):
        assert encode_index(idx, decode_index(idx, number)) == number, number
    assert decode_index(idx, 0) == Lam(Var(1)), "decode_index(0) must be λ1"
    table = Table()
    capped = Table(index_cap=8)
    for number in range(table.closed_cumulative(9)):  # capped agrees on small sizes
        bits = bin(number + 1)[3:]
        assert decode(table, bits) == decode(capped, bits)
    # counts at m > 0 (open-term contexts) match exhaustive enumeration
    for cap in (None, 3):
        ctx_table = Table(index_cap=cap)
        for n in range(2, 13):
            terms_n = _all_blc_terms(n)
            for m in range(5):
                brute = sum(1 for term in terms_n
                            if _max_free(term) <= m
                            and (cap is None or max_de_bruijn_index(term) <= cap))
                assert ctx_table.count(n, m) == brute, (cap, n, m)

    # set_index_cap, every direction, agrees with a freshly built table
    raised = Table(index_cap=2)
    raised.extend(30)
    raised.set_index_cap(7)  # cap bump reuses rows of size <= 3 only
    lowered = Table(index_cap=7)
    lowered.extend(30)
    lowered.set_index_cap(2)
    uncapped = Table(index_cap=3)
    uncapped.extend(25)
    uncapped.set_index_cap(None)
    for built, observed in ((30, raised), (30, lowered), (25, uncapped)):
        reference = Table(index_cap=observed.index_cap)
        reference.extend(built)
        for n in range(4, built + 1):
            assert observed.count(n, 0) == reference.count(n, 0), n

    # encode rejects everything that is not a closed, in-cap term
    table = Table()
    table.extend(14)
    for bad_bits in ("2", "01x", "10 1"):
        try:
            decode(table, bad_bits)
            raise AssertionError(f"decode accepted {bad_bits!r}")
        except ValueError:
            pass
    for bad_term in (Var(1), App(Var(1), Var(1)), Lam(Var(2)), Lam(Var(0))):
        try:
            encode(table, bad_term)
            raise AssertionError(f"encode accepted {bad_term}")
        except ValueError:
            pass
    try:                                    # index 3 exceeds cap 2
        encode(Table(index_cap=2), Lam(Lam(Lam(Var(3)))))
        raise AssertionError("encode accepted an over-cap term")
    except ValueError:
        pass

    # save / load round trip (capped and unbounded) and corruption rejection
    import os
    import tempfile
    for cap, built in ((5, 40), (None, 30)):
        saved = Table(index_cap=cap)
        saved.extend(built)
        with tempfile.NamedTemporaryFile(suffix=".lamtab", delete=False) as f:
            path = f.name
        try:
            saved.save(path)
            loaded = Table.load(path)
            assert loaded.index_cap == cap and loaded.built_size == built
            for n in range(4, built + 1):
                assert loaded.count(n, 0) == saved.count(n, 0), n
            bits = bin(1235 + 1)[3:]
            assert encode(loaded, decode(loaded, bits)) == bits
            good = bytearray(open(path, "rb").read())
            mid = len(good) // 2
            flipped = bytearray(good)
            flipped[mid] ^= 1
            corruptions = [
                (bytes(flipped), "flipped byte"),
                (bytes(good[:-3]), "truncation"),
                (bytes([good[0] ^ 1]) + bytes(good[1:]), "bad magic"),
            ]
            for blob, label in corruptions:
                with open(path, "wb") as handle:
                    handle.write(blob)
                try:
                    Table.load(path)
                    raise AssertionError(f"load accepted {label}")
                except ValueError:
                    pass
        finally:
            os.unlink(path)
    print("self-test passed")


def _print_vectors() -> None:
    for cap, limit in ((None, 500), (3, 300)):
        table = Table(index_cap=cap)
        print(f"# cap={'inf' if cap is None else cap}")
        for number in range(limit):
            bits = bin(number + 1)[3:]
            term = decode(table, bits)
            print(f"{number}\t{show_bits(bits)}\t{show_term(term)}")


if __name__ == "__main__":
    import sys

    if "--vectors" in sys.argv:
        _print_vectors()
    else:
        _self_test()
        table = Table()
        print("first strings of the canonical bijection:")
        for number in range(8):
            bits = bin(number + 1)[3:]
            print(f"  {show_bits(bits):>4}  ->  {show_term(decode(table, bits))}")
