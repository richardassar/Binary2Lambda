#!/usr/bin/env python3
"""Binary2Lambda cross-language table interop check.

Usage:
  python3 check_cross_tables.py load <path>   # load a table saved by any
                                              # implementation, verify it
                                              # against a fresh build
  python3 check_cross_tables.py save <path>   # save a cap-5 size-40 table
                                              # (the C++ --load-table
                                              # counterpart reads it)
"""

import sys

from lambda_bijection import Table, decode, encode


def main() -> None:
    if len(sys.argv) != 3 or sys.argv[1] not in ("save", "load"):
        sys.exit(__doc__.strip())
    mode, path = sys.argv[1], sys.argv[2]
    if mode == "save":
        table = Table(index_cap=5)
        table.extend(40)
        table.save(path)
        print(f"saved cap-5 size-40 table to {path}")
        return
    loaded = Table.load(path)
    fresh = Table(index_cap=loaded.index_cap)
    fresh.extend(loaded.built_size)
    for n in range(2, loaded.built_size + 1):
        assert loaded.count(n, 0) == fresh.count(n, 0), n
    bits = bin(1235 + 1)[3:]
    assert encode(loaded, decode(loaded, bits)) == bits
    cap = "inf" if loaded.index_cap is None else loaded.index_cap
    print(f"table file OK (cap {cap}, size {loaded.built_size})")


if __name__ == "__main__":
    main()
