#!/usr/bin/env python3
"""Binary2Lambda compression layer: lambda-specific compression of closed terms.

This is the compression axis of the project, deliberately separate from the
bijection: a compressor needs compress + decompress, but it does NOT need
"every bit string is a term" (it would not be used where that property
matters, e.g. genetic programming). Use lambda_bijection for the total
bijective representation; use this module to store or transmit terms
compactly.

Method: a renormalizing range coder over the term grammar with an adaptive
model. The coder walks the term in pre-order; at each node the choice set is
{Var(1..m), Lam, App} where m is the number of enclosing lambdas, so

  - zero probability mass is spent on syntactically impossible
    continuations (a byte-level compressor cannot know that a variable
    cannot exceed its lambda depth, or that the walk ends exactly when the
    tree closes);
  - the adaptive counts learn the corpus statistics of constructors
    (conditioned on depth) and of de Bruijn indices (which decay
    geometrically in practice).

Termination is structural (the walk ends when all subtrees are complete), so
no size header is needed; the output is a bare byte stream. The range coder
keeps its state in 32-bit registers and renormalizes a byte at a time, so
both directions run in time linear in the node count, and the model counts
stay bounded by periodic rescaling.

On uniform-random terms the bijective rank code is the Shannon-optimal code
and this coder's output is larger by its adaptive-model overhead; on
structured terms (real programs, repetitive corpora) this coder's output is
smaller. The demo at the bottom reports both.
"""

from __future__ import annotations

import sys

from lambda_bijection import (App, Lam, Table, Term, Var, check_closed, decode,
                              encode)

_KIND_VAR, _KIND_LAM, _KIND_APP = 0, 1, 2

# 32-bit Subbotin range coder. All registers are unsigned 32-bit; arithmetic
# wraps modulo 2^32 (emulated here by masking), which is what makes the byte
# stream identical across the Python, C++ and Rust implementations. Every
# symbol's total frequency must stay below _BOT; the model guarantees that by
# rescaling, and indices above the bucket count are coded bit by bit.
_MASK = 0xFFFFFFFF
_TOP = 1 << 24
_BOT = 1 << 16


# --------------------------------------------------------- range coding


class _RangeEncoder:
    """Narrows [low, low+range) per symbol, emitting settled high bytes."""

    def __init__(self) -> None:
        self.low = 0
        self.range = _MASK
        self.out = bytearray()

    def encode(self, c_low: int, freq: int, total: int) -> None:
        self.range //= total
        self.low = (self.low + ((c_low * self.range) & _MASK)) & _MASK
        self.range = (self.range * freq) & _MASK
        while (((self.low ^ ((self.low + self.range) & _MASK)) & _MASK) < _TOP
               or (self.range < _BOT
                   and self._shrink_range())):
            self.out.append((self.low >> 24) & 0xFF)
            self.low = (self.low << 8) & _MASK
            self.range = (self.range << 8) & _MASK

    def _shrink_range(self) -> bool:
        self.range = (-self.low) & (_BOT - 1)
        return True

    def finish(self) -> bytes:
        for _ in range(4):
            self.out.append((self.low >> 24) & 0xFF)
            self.low = (self.low << 8) & _MASK
        return bytes(self.out)


class _RangeDecoder:
    """Mirrors the encoder; the code register steers symbol choice."""

    def __init__(self, data: bytes) -> None:
        self.data = data
        self.pos = 0
        self.low = 0
        self.range = _MASK
        self.code = 0
        for _ in range(4):
            self.code = ((self.code << 8) | self._next_byte()) & _MASK

    def _next_byte(self) -> int:
        if self.pos < len(self.data):
            byte = self.data[self.pos]
            self.pos += 1
            return byte
        return 0  # past the stream the decoder reads zero bytes

    def target(self, total: int) -> int:
        self.range //= total
        value = ((self.code - self.low) & _MASK) // self.range
        return total - 1 if value >= total else value

    def consume(self, c_low: int, freq: int) -> None:
        self.low = (self.low + ((c_low * self.range) & _MASK)) & _MASK
        self.range = (self.range * freq) & _MASK
        while (((self.low ^ ((self.low + self.range) & _MASK)) & _MASK) < _TOP
               or (self.range < _BOT
                   and self._shrink_range())):
            self.code = ((self.code << 8) | self._next_byte()) & _MASK
            self.low = (self.low << 8) & _MASK
            self.range = (self.range << 8) & _MASK

    def _shrink_range(self) -> bool:
        self.range = (-self.low) & (_BOT - 1)
        return True


def _encode_symbol(coder: _RangeEncoder, weights: list[int],
                   symbol: int) -> None:
    c_low = sum(weights[:symbol])
    coder.encode(c_low, weights[symbol], sum(weights))


def _decode_symbol(coder: _RangeDecoder, weights: list[int]) -> int:
    total = sum(weights)
    target = coder.target(total)
    c_low = 0
    for symbol, weight in enumerate(weights):
        if target < c_low + weight:
            coder.consume(c_low, weight)
            return symbol
        c_low += weight
    raise ValueError("malformed compressed data (code point out of range)")


def _encode_bit(coder: _RangeEncoder, bit: int) -> None:
    coder.encode(bit, 1, 2)


def _decode_bit(coder: _RangeDecoder) -> int:
    bit = 1 if coder.target(2) >= 1 else 0
    coder.consume(bit, 1)
    return bit


# Elias-gamma over an equiprobable bit model, used for the rare de Bruijn
# index that exceeds the bucket count.  _GAMMA_MAX_BITS caps the unary length
# so malformed input cannot spin in the prefix.
_GAMMA_MAX_BITS = 40


def _encode_gamma(coder: _RangeEncoder, value: int) -> None:
    v = value + 1
    n = v.bit_length() - 1
    for _ in range(n):
        _encode_bit(coder, 0)
    _encode_bit(coder, 1)
    for i in range(n - 1, -1, -1):
        _encode_bit(coder, (v >> i) & 1)


def _decode_gamma(coder: _RangeDecoder) -> int:
    n = 0
    while _decode_bit(coder) == 0:
        n += 1
        if n > _GAMMA_MAX_BITS:
            raise ValueError("malformed compressed data (index code too long)")
    v = 1
    for _ in range(n):
        v = (v << 1) | _decode_bit(coder)
    return v - 1


# --------------------------------------------------------- adaptive model


class _AdaptiveModel:
    """Counts shared by compress and decompress, updated identically.

    Kinds are conditioned on the lambda depth bucket min(m, 3). Variable
    indices share an 8-entry bucketed frequency table (bucket min(i, 8) - real
    index usage decays geometrically, so the tail shares one bucket); an index
    in the tail bucket is followed by its offset, coded bit by bit. Counts are
    halved when a context total reaches _RESCALE_LIMIT, keeping every total
    below the coder's _BOT bound.
    """

    _INDEX_BUCKETS = 8
    _INCREMENT = 16
    _RESCALE_LIMIT = 1 << 14

    def __init__(self) -> None:
        self._kinds = [[1, 1, 1] for _ in range(4)]
        self._indices = [1] * self._INDEX_BUCKETS

    def kind_weights(self, m: int) -> list[int]:
        var_weight, lam_weight, app_weight = self._kinds[min(m, 3)]
        return [var_weight if m >= 1 else 0, lam_weight, app_weight]

    def saw_kind(self, m: int, kind: int) -> None:
        context = self._kinds[min(m, 3)]
        context[kind] += self._INCREMENT
        if sum(context) >= self._RESCALE_LIMIT:
            for i in range(3):
                context[i] = (context[i] >> 1) or 1

    def index_bucket_weights(self, alphabet: int) -> list[int]:
        return self._indices[:alphabet]

    def saw_index(self, index: int) -> None:
        self._indices[min(index - 1, self._INDEX_BUCKETS - 1)] += self._INCREMENT
        if sum(self._indices) >= self._RESCALE_LIMIT:
            for i in range(self._INDEX_BUCKETS):
                self._indices[i] = (self._indices[i] >> 1) or 1


def _encode_index(coder: _RangeEncoder, model: _AdaptiveModel, m: int,
                  index: int) -> None:
    alphabet = min(m, model._INDEX_BUCKETS)
    bucket = min(index - 1, model._INDEX_BUCKETS - 1)
    _encode_symbol(coder, model.index_bucket_weights(alphabet), bucket)
    if bucket == model._INDEX_BUCKETS - 1:
        _encode_gamma(coder, index - model._INDEX_BUCKETS)
    model.saw_index(index)


def _decode_index(coder: _RangeDecoder, model: _AdaptiveModel, m: int) -> int:
    alphabet = min(m, model._INDEX_BUCKETS)
    bucket = _decode_symbol(coder, model.index_bucket_weights(alphabet))
    if bucket < model._INDEX_BUCKETS - 1:
        index = bucket + 1
    else:
        index = model._INDEX_BUCKETS + _decode_gamma(coder)
        if index > m:
            raise ValueError("malformed compressed data (index exceeds depth)")
    model.saw_index(index)
    return index


# ------------------------------------------------------------ public API

# Termination guards. The node ceiling caps a stream that never closes the
# tree; the depth cap (below the interpreter recursion limit) caps a stream
# that keeps nesting. Both directions are linear in the node count, so either
# limit is reached cheaply.
_MAX_DECODE_NODES = 1 << 20
_MAX_DECODE_DEPTH = 12000


def compress(term: Term) -> bytes:
    """Closed lambda term -> compact bytes (inverse of decompress).

    The output is the range coder's byte stream (a four-byte flush tail makes
    it self-delimiting against the structural end of the walk). Raises
    ValueError on a non-closed or malformed term, like encode."""
    check_closed(term)
    coder = _RangeEncoder()
    model = _AdaptiveModel()
    _walk_encode(term, 0, model, coder)
    # The decoder reads zero bytes past the stream, so trailing zero bytes are
    # redundant; drop them (an all-zero stream becomes empty).
    return coder.finish().rstrip(b"\x00")


def decompress(data: bytes) -> Term:
    """Bytes produced by compress -> the original closed lambda term.

    Any byte string terminates: it yields a term or raises ValueError."""
    coder = _RangeDecoder(data)
    model = _AdaptiveModel()
    return _walk_decode(0, 0, model, coder, [_MAX_DECODE_NODES])


def compressed_bits(data: bytes) -> int:
    """Size of a compress() result in bits."""
    return 8 * len(data)


def _walk_encode(term: Term, m: int, model: _AdaptiveModel,
                 coder: _RangeEncoder) -> None:
    weights = model.kind_weights(m)
    if isinstance(term, Var):
        _encode_symbol(coder, weights, _KIND_VAR)
        model.saw_kind(m, _KIND_VAR)
        _encode_index(coder, model, m, term.index)
    elif isinstance(term, Lam):
        _encode_symbol(coder, weights, _KIND_LAM)
        model.saw_kind(m, _KIND_LAM)
        _walk_encode(term.body, m + 1, model, coder)
    else:
        _encode_symbol(coder, weights, _KIND_APP)
        model.saw_kind(m, _KIND_APP)
        _walk_encode(term.fun, m, model, coder)
        _walk_encode(term.arg, m, model, coder)


def _walk_decode(m: int, depth: int, model: _AdaptiveModel,
                 coder: _RangeDecoder, budget: list[int]) -> Term:
    budget[0] -= 1
    if budget[0] < 0 or depth > _MAX_DECODE_DEPTH:
        raise ValueError("malformed compressed data (term does not terminate)")
    kind = _decode_symbol(coder, model.kind_weights(m))
    model.saw_kind(m, kind)
    if kind == _KIND_VAR:
        return Var(_decode_index(coder, model, m))
    if kind == _KIND_LAM:
        return Lam(_walk_decode(m + 1, depth + 1, model, coder, budget))
    fun = _walk_decode(m, depth + 1, model, coder, budget)
    return App(fun, _walk_decode(m, depth + 1, model, coder, budget))


# ------------------------------------------------------------- self-test


def _church(n: int) -> Term:
    body: Term = Var(1)
    for _ in range(n):
        body = App(Var(2), body)
    return Lam(Lam(body))


def _lam_chain(n: int) -> Term:
    """A nest of n lambdas over Var(1): n+1 nodes whose compressed size is a
    few bytes, so its node count far exceeds the byte count."""
    body: Term = Var(1)
    for _ in range(n):
        body = Lam(body)
    return body


def _vector_terms() -> list[tuple[str, Term]]:
    """Fixed term set shared with the C++ and Rust implementations; their
    compressed bytes must match Python's exactly."""
    s_comb = Lam(Lam(Lam(App(App(Var(3), Var(1)), App(Var(2), Var(1))))))
    y_half = Lam(App(Var(2), App(Var(1), Var(1))))
    repetitive: Term = s_comb
    for _ in range(5):
        repetitive = App(repetitive, repetitive)
    uniform = decode(Table(), format(987654322, "0192b"))
    return [("S", s_comb), ("Y", Lam(App(y_half, y_half))),
            ("church10", _church(10)), ("church100", _church(100)),
            ("rep32S", repetitive), ("uniform192", uniform)]


def _print_vectors() -> None:
    for name, term in _vector_terms():
        print(f"{name}\t{compress(term).hex()}")


def _self_test() -> None:
    from lambda_bijection import _brute_force_terms

    cases = []
    for n in range(4, 14):
        cases.extend(_brute_force_terms(n, None))
    table = Table()
    cases.extend(decode(table, bin((1 << length) + 12345 + 1)[3:])
                 for length in (32, 64, 128))
    cases.extend([_church(5), _church(50)])
    # Highly compressible terms whose node count far exceeds their byte count,
    # exercising the node ceiling that bounds decompression.
    cases.extend([_lam_chain(48), _lam_chain(1000), _church(2000)])
    cases.extend(term for _, term in _vector_terms())
    for term in cases:
        assert decompress(compress(term)) == term, term

    # compress rejects non-closed / malformed terms, like encode
    for bad_term in (Var(1), Lam(Var(2)), App(Var(1), Var(1)), Lam(Var(0))):
        try:
            compress(bad_term)
            raise AssertionError(f"compress accepted {bad_term}")
        except ValueError:
            pass

    # decompress is always well-behaved on malformed bytes: it returns some
    # term or raises ValueError, never looping, overflowing, or crashing.  Its
    # domain is compress outputs, so (lacking an integrity check) it need not
    # detect every corrupt blob -- only stay well-behaved on all of them.  The
    # constant-fill blobs once drove the decoder to grind; a linear coder makes
    # them resolve to a term or an error quickly.
    bad_blobs = [b"", b"\x00", b"\xff",
                 bytes([0x00]) + b"\xff" * 8,
                 bytes([0xff]) + b"\xff" * 16,
                 b"\xaa" * 64,
                 compress(cases[0])[:-1]]
    for blob in bad_blobs:
        try:
            assert isinstance(decompress(blob), (Var, Lam, App))
        except ValueError:
            pass
    print(f"self-test passed ({len(cases)} round trips + error paths)")


def _demo() -> None:
    import gzip

    from lambda_bijection import term_size

    table = Table()

    def gzip_bits(term: Term) -> int:
        rank_bits = encode(table, term)
        packed = int("1" + rank_bits, 2).to_bytes(
            (len(rank_bits) + 8) // 8, "big")
        return 8 * len(gzip.compress(packed, 9))

    labels = {"S": "S combinator", "Y": "Y combinator",
              "church10": "Church numeral 10", "church100": "Church numeral 100",
              "rep32S": "32 x S combinator (repetitive)",
              "uniform192": "uniform 192-bit term"}
    print(f"\n{'term':<32} {'BLC':>6} {'rank code':>10} "
          f"{'compressed':>11} {'gzip(rank)':>11}")
    for key, term in _vector_terms():
        blc_bits = term_size(term)
        rank_bits = len(encode(table, term))
        comp_bits = compressed_bits(compress(term))
        print(f"{labels[key]:<32} {blc_bits:>6} {rank_bits:>10} "
              f"{comp_bits:>11} {gzip_bits(term):>11}")
    print("\nbits shown; compressed = range-coder byte stream including its "
          "four-byte flush tail; gzip column includes its ~18-byte format "
          "overhead")


if __name__ == "__main__":
    if "--vectors" in sys.argv:
        _print_vectors()
    else:
        _self_test()
        _demo()
