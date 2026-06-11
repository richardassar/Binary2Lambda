#!/usr/bin/env python3
"""Binary2Lambda: emit validated constants in WL syntax for the .wl self-test."""

from lambda_bijection import Lam, Table, Var, decode

def wl_term(t):
    if isinstance(t, Var):
        return f"LambdaVar[{t.index}]"
    if isinstance(t, Lam):
        return f"LambdaAbs[{wl_term(t.body)}]"
    return f"LambdaApp[{wl_term(t.fun)}, {wl_term(t.arg)}]"

unbounded = Table()
capped = Table(index_cap=2)
print("uncapped counts n=4..14:",
      "{" + ", ".join(str(unbounded.count(n, 0)) for n in range(4, 15)) + "}")
print("cap-2 counts n=4..14:   ",
      "{" + ", ".join(str(capped.count(n, 0)) for n in range(4, 15)) + "}")
terms = [wl_term(decode(unbounded, bin(number + 1)[3:])) for number in range(10)]
print("first 10 terms: {" + ",\n   ".join(terms) + "}")
