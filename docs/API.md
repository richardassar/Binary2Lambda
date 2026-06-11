# Binary2Lambda — API reference

The programmer-facing surface of all four implementations. For the theory and
the design rationale see [../NOTES.md](../NOTES.md); for a guided tour see
[../README.md](../README.md).

Two layers are documented here:

- **Bijection** (`lambda_bijection`): the total map between bit strings and
  closed de Bruijn lambda terms — `decode`, `encode`, and the counting
  `Table` they share. Every bit string is a valid term; every closed term has
  exactly one bit string.
- **Compression** (`lambda_compress`, Python/C++/Rust): a renormalizing range
  coder over the term grammar with an adaptive model — `compress`,
  `decompress`. A
  separate axis: it stores a term compactly but does not give the
  every-string-is-a-term property.

Conventions used below: a *de Bruijn index* is `>= 1` (1 = the innermost
binder). A *closed* term has every index bound (no index exceeds its binder
depth). A *bit string* is a string over `{'0','1'}`; the empty string is
valid and denotes index 0 of the enumeration.

---

## 1. Term representation

| | construct a term | inspect |
|---|---|---|
| **Python** | `Var(index)`, `Lam(body)`, `App(fun, arg)` (dataclasses) | `.index`, `.body`, `.fun`, `.arg`; `isinstance` |
| **C++** | `var(int)`, `lam(TermPtr)`, `app(TermPtr, TermPtr)` → `TermPtr` (`shared_ptr<const Term>`) | `term->kind` (`Term::Kind::{Var,Lam,App}`), `term->index`, `term->left`, `term->right` |
| **Rust** | `Term::Var(u32)`, `Term::Lam(Box<Term>)`, `Term::App(Box<Term>, Box<Term>)` | `match` on the enum |
| **Wolfram** | `LambdaVar[i]`, `LambdaAbs[body]`, `LambdaApp[f, x]` | pattern matching |

Terms are immutable values; share them freely (C++ uses `shared_ptr`).

---

## 2. Bijection API

### 2.1 Term utilities

| operation | Python | C++ | Rust |
|---|---|---|---|
| BLC size | `term_size(t) -> int` | `termSize(t) -> int` | `term_size(&t) -> usize` |
| largest free index | `max_de_bruijn_index(t) -> int` | `maxDeBruijnIndex(t) -> int` | `max_de_bruijn_index(&t) -> u32` |
| validate closed | `check_closed(t, depth=0)` | `checkClosed(t, depth=0)` | `check_closed(&t, depth) -> Result<(),String>` |
| structural equality | `==` | `sameTerm(a, b) -> bool` | `==` |
| render term | `show_term(t, context=0) -> str` | `showTerm(t, context=0) -> string` | `show_term(&t, context) -> String` |
| render bit string | `show_bits(bits) -> str` | `showBits(bits) -> string` | `show_bits(&bits) -> &str` |

- **`term_size`** — BLC size: `|Var i| = i+1`, `|Lam b| = |b|+2`, `|App f a| = |f|+|a|+2`.
- **`max_de_bruijn_index`** — the largest index appearing anywhere; `0` for a
  term with no variables. Use it to choose a cap that admits a given term
  (`cap >= max_de_bruijn_index(t)`).
- **`check_closed`** — returns normally if the term is closed; otherwise
  signals an error (Python/C++ raise `ValueError`/`std::invalid_argument`;
  Rust returns `Err`). `depth` is the number of enclosing binders assumed
  around `t` (0 for a top-level closed term).
- **`show_term`** — human-readable, e.g. `(λ1 (λλ1)) (λ1)`. `context` is an
  internal precedence level; pass `0`.

### 2.2 The counting table

The `Table` holds the counting recurrence `T(n, m)` (number of closed-under-m
terms of BLC size `n` with free indices `<= m`), plus an optional de Bruijn
index cap. `decode`/`encode` consume it and grow it on demand.

| operation | Python | C++ | Rust |
|---|---|---|---|
| construct | `Table(index_cap=None)` | `Table(optional<int> cap = {})` | `Table::new(cap: Option<u32>)` |
| grow to size | `extend(size_limit)` | `extend(int)` | `extend(usize)` |
| change cap | `set_index_cap(new_cap)` | `setIndexCap(optional<int>)` | `set_index_cap(Option<u32>)` |
| built size | `built_size() -> int` | `builtSize() -> int` | `built_size() -> usize` |
| count entry | `count(n, m) -> int` | `count(n, m) -> BigNat` | `count_built(n, m) -> &BigNat` |
| persist | `save(path)` | `saveToFile(path)` | `save_to_file(path) -> io::Result` |
| load | `Table.load(path)` (classmethod) | `Table::loadFromFile(path)` (static) | `Table::load_from_file(path) -> Result` |

- **`index_cap` / `cap`** — `None`/`nullopt`/`None` means unbounded indices.
  A finite cap `K` restricts the enumeration to terms whose indices never
  exceed `K`. **Each cap value (unbounded included) is a different
  bijection** — fix it before encoding anything and use the same value on
  both sides. Capped and unbounded bijections agree on all terms of size
  `<= K+1`.
- **`extend`** — builds counts up to the given BLC size. Idempotent and
  append-only: existing rows never change. `decode`/`encode` call it
  automatically when an input needs a larger size, so explicit `extend` is an
  optimization (prebuild once, then no growth at run time).
- **`set_index_cap`** — switches the cap in place, reusing every row of size
  `<= min(old, new) + 1`.
- **`count(n, m)`** — `T(n, m)` as an arbitrary-precision integer.

### 2.3 The bijection

| operation | Python | C++ | Rust |
|---|---|---|---|
| decode (string → term) | `decode(table, bits) -> Term` | `decode(table, bits) -> TermPtr` | `decode(&mut table, &bits) -> Result<Term,String>` |
| encode (term → string) | `encode(table, term) -> str` | `encode(table, term) -> string` | `encode(&mut table, &term) -> Result<String,String>` |

- **`decode`** — maps a bit string to its closed term. Total over the chosen
  cap: every `{0,1}*` string is valid. The integer index is
  `N = int("1" + bits, 2) - 1` (the implicit leading 1 makes all `2^L`
  strings of every length distinct); `N` indexes the size-ordered
  enumeration. Raises/`Err` only on a non-`{0,1}` character.
- **`encode`** — the exact inverse: maps a closed term to the unique bit
  string that decodes to it. Raises/`Err` on a non-closed term, or on a term
  whose largest index exceeds the table's cap (it is not in that bijection).
  `encode(table, decode(table, bits)) == bits` for every `bits`;
  `decode(table, encode(table, t)) == t` for every closed in-cap `t`.

```python
from lambda_bijection import Table, decode, encode, show_term
table = Table()                       # unbounded
t = decode(table, "010001100")
assert show_term(t) == "(λ1 (λλ1)) (λ1)"
assert encode(table, t) == "010001100"
```

---

## 3. Compression API

`lambda_compress` (Python), and the `lambda_compress` namespace / module-level
functions in the C++ and Rust files. Byte-compatible across the three: the
same term compresses to the same bytes everywhere.

| operation | Python | C++ | Rust |
|---|---|---|---|
| compress | `compress(term) -> bytes` | `lambda_compress::compress(term) -> vector<uint8_t>` | `compress(&term) -> Result<Vec<u8>,String>` |
| decompress | `decompress(data) -> Term` | `lambda_compress::decompress(data) -> TermPtr` | `decompress(&data) -> Result<Term,String>` |
| size in bits | `compressed_bits(data) -> int` | — | — |

- **`compress`** — a renormalizing range coder over the term grammar with an
  adaptive model. The output is the coder's byte stream; its four-byte flush
  tail makes it self-delimiting against the structural end of the walk, so
  there is no size header. Raises/`Err` on a non-closed term (like `encode`).
- **`decompress`** — the inverse. On any malformed input it stays
  well-behaved: it returns some term or raises/`Err`, never looping,
  overflowing, or crashing, and runs in time linear in the node count. Its
  domain is `compress` outputs; lacking an integrity check it need not detect
  every corrupt blob, only stay well-behaved. Termination is structural (the
  grammar walk ends when the tree closes), guarded against adversarial input
  by a node ceiling (`2^20`) and a nesting-depth cap (`12000`).
- **`compressed_bits`** — the size of a blob in bits (`8 * len`).

`compress`/`decompress` consume and produce terms, not bit strings, so they
do not give the every-string-is-a-term property; use the bijection for that.

---

## 4. On-disk table format

One portable text format, byte-identical across Python, C++ and Rust and
cross-loadable between them. Lines, in order:

```
lambda-binarization-table v1
cap <K|inf>
size <built_size>
<hex counts for size 0>
<hex counts for size 1>
...
checksum <fnv1a64-hex>
```

Each count line holds the row `T(n, 0..width(n)-1)` as space-separated
lowercase hex. The checksum is FNV-1a-64 over the body (everything before the
checksum line), streamed so memory stays flat. Load verifies the checksum and
rejects any corruption or truncation (flipped digit, dropped row, missing
checksum). The Wolfram port persists its memoized store as WL data via
`Put`/`Get` (a different, non-interchangeable format).

---

## 5. Command-line interface (C++ and Rust binaries)

The C++ (`lambda_bijection_cpp`) and Rust (`lambda_bijection_rs`) binaries
take the same flags:

| invocation | effect |
|---|---|
| *(no arguments)* | run the full self-test, then print the first strings of the bijection |
| `--vectors` | print decode vectors (`number`, bits, term) for cap `inf` and cap `3` |
| `--compress-vectors` | print `name<TAB>hex` for the shared compression vector terms |
| `--fuzz <file>` | decompress each hex-blob line; print `OK<TAB><digest><TAB><nodes>` or `ERR` (drives the differential fuzz) |
| `--save-table <path>` | write a sample cap-5 size-40 table |
| `--load-table <path>` | load a table file and verify its counts |
| `--bench <cap> <length>` | run one benchmark block, emitting CSV (`cap` is `inf` or an integer; `length >= 5`) |

The binary exits non-zero with a one-line `error:` message on any failure
(self-test assertion, bad argument, corrupt table). The `--fuzz` digest is a
structural FNV-1a-64 over the decoded term, computed identically in all three
languages, so feeding one blob file to each and diffing the outputs proves
the decoders agree on every input. See [../verify/README.md](../verify/README.md).

Python entry points:

- `python3 python/lambda_bijection.py` — bijection self-test.
- `python3 python/lambda_compress.py` — compression self-test and a demo
  table; `--vectors` prints the compression vectors.

---

## 6. Error model

| condition | Python | C++ | Rust |
|---|---|---|---|
| non-`{0,1}` char in `decode` | `ValueError` | `std::invalid_argument` | `Err(String)` |
| non-closed term in `encode`/`compress` | `ValueError` | `std::invalid_argument` | `Err(String)` |
| index over cap in `encode` | `ValueError` | exception | `Err(String)` |
| malformed `decompress` input | `ValueError` | `std::runtime_error` | `Err(String)` |
| corrupt/truncated table file | `ValueError` | exception | `Err(String)` |

The Rust API reports bad input as `Err(String)`; bad input never panics.
Debug builds additionally check internal invariants with `debug_assert!`
(C++: `assert`), e.g. big-integer subtraction never underflows and the
decoder target always fits in 64 bits; these are compiled out of release
builds.
