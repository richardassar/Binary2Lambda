#!/usr/bin/env python3
"""Binary2Lambda performance profile: scaling of the bijection (Python reference).

For each regime (unbounded, and de Bruijn cap 32) and each bit length L:
  - table construction: wall time and tracemalloc high-water mark,
  - decode (number -> lambda): mean / std / min / max over fixed samples,
  - encode (lambda -> number): same statistics,
  - table footprint: entry count and total bits of stored counts.

Sample strings are deterministic: 24 indices spread uniformly across the
length-L block of the numeration. Empirical scaling exponents are reported
between successive (doubling) values of L. C++ follows the same curves with
a much smaller constant; the exponents are language-independent.
"""

import math
import os
import statistics
import tempfile
import time
import tracemalloc

from lambda_bijection import Table, decode, encode

SAMPLES_PER_LENGTH = 24


def profile_block(cap, length):
    tracemalloc.start()
    table = Table(index_cap=cap)
    top = (1 << (length + 1)) - 2  # largest index with bit length `length`

    t0 = time.perf_counter()
    n_max = 4
    while table.closed_cumulative(n_max) <= top:
        n_max += 1
    build_seconds = time.perf_counter() - t0

    numbers = [(1 << length) - 1 + j * ((1 << length) // SAMPLES_PER_LENGTH)
               for j in range(SAMPLES_PER_LENGTH)]
    decode_times, encode_times, terms = [], [], []
    for number in numbers:
        bits = bin(number + 1)[3:]
        t0 = time.perf_counter()
        term = decode(table, bits)
        decode_times.append(time.perf_counter() - t0)
        terms.append((term, bits))
    for term, bits in terms:
        t0 = time.perf_counter()
        assert encode(table, term) == bits
        encode_times.append(time.perf_counter() - t0)

    peak_bytes = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()

    entries = table.entry_count()
    table_bits = table.total_bit_length()
    with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as handle:
        disk_path = handle.name
    table.save(disk_path)
    disk_kib = os.path.getsize(disk_path) / 1024
    os.unlink(disk_path)
    return {
        "L": length, "n_max": n_max, "build_s": build_seconds,
        "peak_mb": peak_bytes / 2**20,
        "entries": entries, "table_kib": table_bits / 8 / 1024,
        "disk_kib": disk_kib,
        "dec": decode_times, "enc": encode_times,
    }


def stats_ms(times):
    mean = statistics.mean(times) * 1e3
    std = (statistics.stdev(times) if len(times) > 1 else 0.0) * 1e3
    return mean, std, min(times) * 1e3, max(times) * 1e3


def report(regime_name, lengths, cap):
    print(f"\n=== regime: {regime_name} ===")
    header = (f"{'L':>5} {'n_max':>6} {'build s':>8} {'peak MB':>8} "
              f"{'entries':>8} {'tbl KiB':>8} "
              f"{'dec ms':>8} {'±std':>7} {'min':>7} {'max':>7} "
              f"{'enc ms':>8} {'±std':>7}")
    print(header)
    rows = []
    for length in lengths:
        row = profile_block(cap, length)
        rows.append(row)
        d_mean, d_std, d_min, d_max = stats_ms(row["dec"])
        e_mean, e_std, _, _ = stats_ms(row["enc"])
        print(f"{row['L']:>5} {row['n_max']:>6} {row['build_s']:>8.3f} "
              f"{row['peak_mb']:>8.2f} {row['entries']:>8} "
              f"{row['table_kib']:>8.1f} "
              f"{d_mean:>8.3f} {d_std:>7.3f} {d_min:>7.3f} {d_max:>7.3f} "
              f"{e_mean:>8.3f} {e_std:>7.3f}")
    print("empirical exponents between successive L "
          "(value ~ L^e; L doubles each step):")
    for prev, cur in zip(rows, rows[1:]):
        ratio = math.log2(cur["L"] / prev["L"])

        def exponent(prev_value, cur_value):
            return math.log2(max(cur_value, 1e-12)
                             / max(prev_value, 1e-12)) / ratio

        print(f"  L {prev['L']:>4} -> {cur['L']:>4}:"
              f"  build {exponent(prev['build_s'], cur['build_s']):>5.2f}"
              f"  peak-mem {exponent(prev['peak_mb'], cur['peak_mb']):>5.2f}"
              f"  table-bits {exponent(prev['table_kib'], cur['table_kib']):>5.2f}"
              f"  decode {exponent(statistics.mean(prev['dec']), statistics.mean(cur['dec'])):>5.2f}"
              f"  encode {exponent(statistics.mean(prev['enc']), statistics.mean(cur['enc'])):>5.2f}")


if __name__ == "__main__":
    report("unbounded indices (canonical bijection)",
           [16, 32, 64, 128, 256], cap=None)
    report("de Bruijn cap K = 32",
           [64, 128, 256, 512, 1024], cap=32)
