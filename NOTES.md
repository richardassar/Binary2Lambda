# Binary2Lambda — Notes

Bijection between binary strings and closed untyped lambda terms: every
binary string (including the empty string) denotes exactly one closed term,
every closed term has exactly one string, no error states. This document is
the design record: theory, complexity, the two-axis incremental table,
design rationale with literature, the measured scaling curves, and the
positioning against prior work. The README covers usage; this file covers
*why*.

---

## 1. Problem statement and constraints

**Firm constraints**
- Total: every bit string of every length decodes; no error states.
- Bijective: every closed term is covered exactly once.
- Implicit leading digit: string `s` corresponds to `N = int('1' + s, 2) − 1`
  (bijective binary numeration), so all 2^L strings of each length L are
  distinct values and the empty string is N = 0.

**Relaxed constraint**
- number→lambda (decode) is the performance-critical direction;
  lambda→number (encode) is secondary. (This costs nothing: encode turns
  out to be the cheaper direction — see §4.)

**Out of scope but adjacent**
- Compression is an orthogonal axis on top of the raw representation (§9).

## 2. The fundamental picture

A de Bruijn lambda term is a unary-binary (Motzkin) tree with one
decoration and one constraint:

```
Term(m) → Var(i)   1 ≤ i ≤ m     (m = number of enclosing lambdas)
        | Lam Term(m+1)
        | App Term(m) Term(m)
closed term = Term(0)
```

Size measure (= bit length of Tromp's binary lambda calculus code):
`|Var i| = i+1`, `|Lam b| = |b|+2`, `|App f a| = |f|+|a|+2`.

**Structure theorem.** Any bijection {0,1}* ≅ closed terms factors as
numeration ∘ mixed-radix positional system whose place values are the
grammar's completion counts. The decoder's invariant is an exact partition
of `[0, T(n,m))` into var ⊎ lam ⊎ app blocks; the app block is a
Catalan-style convolution `Σ_k T(k,m)·T(n−2−k,m)`. Consequences:

- If the counts were powers of two, raw bit-peeling would already be a
  bijection. They are not (growth ρ ≈ 1.963 ≠ 2), and the deficit is
  exactly BLC's redundancy. Any bijection therefore needs exact counts.
- "Exhaustive combinatorics" means exhaustive **counting** (a polynomial
  recurrence), never enumeration of terms. Decoding a ~200-bit string takes
  microseconds (C++) to milliseconds (Python), not the ~10⁶⁰ steps of
  enumerate-to-index.
- Prefix-free + bijective is impossible (both `0` and `00` are codewords;
  Kraft sum of {0,1}* diverges), so no streaming O(1)/bit bijection exists.
  This is the property Tromp sacrificed bijectivity to keep.

## 3. The counting table

`T(n, m)` = number of terms of size n, free indices ≤ m, all indices ≤ cap K
(K = ∞ for the canonical bijection):

```
T(n, m) = [n−1 ≤ min(m, K)]                      # the variable n−1
        + T(n−2, m+1)                            # abstractions
        + Σ_{k=2..n−4} T(k, m) · T(n−2−k, m)     # applications
closed terms of size n = T(n, 0);   CUM(n) = Σ_{j≤n} T(j, 0)
```

**Saturation lemmas** (these drive everything):
- A size-n term cannot contain an index above n−1, so
  `T(n, m) = T(n, min(m, n−1, K))` — rows are stored under that effective
  context, and the effective cap is `min(K, n−1)`.
- `T(n, m) ≤ 2^n` (each term has a distinct n-bit BLC code), so with a
  string bound L every entry fits in ~L+8 bits — at L ≤ ~60, machine words.
- Rows of size n ≤ K+1 are **cap-independent**: identical for every cap
  ≥ n−1, including ∞.

**No shortcut exists**: closed terms have no finite admissible
specification (the m-indexed system is irreducibly infinite — Bendkowski
2020 survey, §2.5), so no holonomic/closed-form trick can replace the
table; Θ(L²) exact numbers is the memory floor of the canonical bijection.

## 4. The bijection

Canonical order: ascending size; within a size class Var < Lam < App;
abstractions by body rank; applications by (left size, left rank, right
rank). Global index N = CUM(size−1) + in-class rank; bits = binary(N+1)
with the leading 1 dropped.

- **decode** (number→lambda): locate the size class in CUM (binary
  search), then walk the grammar splitting the rank against the partition;
  app nodes use one divmod. This is the direction with searches.
- **encode** (lambda→number): the same walk in reverse with **no searches**
  — the term itself names its blocks, so offsets are direct lookups.

The two directions traverse the identical partition in opposite
directions; they are inverses by construction, share all cached state, and
encode is the cheaper one (so the relaxed constraint costs nothing).
Per-query transient memory is Θ(L): divmod conserves bit length, so the
pending remainders on the stack never total more than L bits.

**De Bruijn cap semantics.** Each cap K (including ∞) defines a DIFFERENT
bijection — onto closed terms with indices ≤ K. Encode and decode must use
the same K (treat K as a codec version). Capped and unbounded bijections
agree on all terms of size ≤ K+1. `MaxDeBruijnIndex` supports the
workflow "measure the term, then pick any table with cap ≥ that".

**If you know the bit bound but not the index bound** (the common case):
no decision is needed. A length-L string decodes to a term of size
n(L) ≈ 1.028·L, and a size-n term cannot contain an index above n−1 — so
the bit bound already implies an index bound, and the canonical (∞) table
truncated at n(L) is automatically the right object. The cap axis is purely
an optional memory optimization for *large* L, never a requirement.

## 5. Complexity

Let L = string bit length, n_max ≈ 1.028·L, w = word size.

**Resident memory (table)**

| configuration | entries | bits | decode cost |
|---|---|---|---|
| T only, cap ∞ (the four shipped implementations) | Θ(L²) | Θ(L³) | O(L²) mults (splits re-scanned in-walk) |
| T only, cap K | Θ(K·L) | Θ(K·L²) | O(L²) mults |
| + SPLIT prefix sums | Θ(L³) words | Θ(L⁴) | O(L log L) |
| + JUMP tables | +Θ(L³) small words | — | Θ(L) expected (measured 1.4 probes/app node) |

Concrete anchors (measured): the uncapped T-table holds 47 MiB of count
*data* in memory at L=1024 and writes a 43 MiB file (the L³/4-bit
formula is asymptotic and under-predicts ~30% here); the capped K=32 table
is 2.2 MiB in memory / 2.1 MiB on disk at L=1024. SPLIT+JUMP ≈ 1 MB at
L=60. Save and load stream the file through an incremental checksum, so
peak process RSS stays near the in-memory table size (~80 MB at L=1024
uncapped) rather than doubling for the file buffer.

**Time**
- Build: Θ(L³) multiplications uncapped, Θ(K·L²) capped; append-only
  marginal cost `(b³−a³)/3` for an extension a→b (measured ratios match).
- Decode/encode after build: Θ(L) word ops with the full accelerator stack
  in the word regime (~1.9/~0.7 ops per bit measured); O(L²) bigint-ops
  with the T-only core (measured exponent ~1.8–2.1 in Python).
- Lower bound Ω(L) both directions (read every bit, emit Θ(L) nodes);
  the accelerated word-regime decoder matches that bound.

**Bounded-table ladder** (what you can cut and what it costs)
1. Drop JUMP: decode Θ(L) → O(L log L). Negligible.
2. Drop SPLIT (T-core only): memory Θ(L⁴)→Θ(L³) bits; ops O(L²).
3. Cap the table below the string's needs: hard wall — a capped table is a
   bounded universe (strings ≤ ~0.97·C bits); longer strings *require*
   growing the core. The Θ(L²)-entries floor is fundamental (no finite
   specification), not implementational.
4. Cap de Bruijn indices at K: memory falls to Θ(K·L) entries; the term
   side becomes "closed terms with indices ≤ K" (deep indices are
   exponentially rare — measured index frequencies decay geometrically).
5. Research direction: chunked bijections with a structural combiner
   `List(Term) ≅ Term` would give memory poly(C) independent of L; no such
   combiner is known for binding syntax (LCRS works only for free algebras).

**Practical ceilings**: exact canonical bijection to L ≈ low thousands of
bits; capped variant comfortably to L ≈ 10⁴–10⁵.

## 6. Incrementality (both axes) — and what it costs

- **Size bump**: append-only. Row n depends only on rows < n; incremental
  total cost = one-shot cost; the *map itself is cap- and bound-independent*
  (tables are a materialized prefix of one infinite array), verified by
  map-stability tests.
- **Cap bump K→K′**: rows of size ≤ min(K,K′)+1 are reused verbatim
  (cap-independence lemma); the rest are recomputed. This is not a defect
  of incrementalization — different caps are different counting
  functions. The WL implementation keys its memo by normalized cap, so
  cross-cap sharing happens automatically.

**Does incrementalization lose anything?** Two things, both recoverable:
1. Row-at-a-time construction forecloses *batch convolution*: with a known
   bound, all app-sums for a whole m-stripe form polynomial products that
   fast multiplication could batch (Õ(L²) instead of Θ(L³) multiplications).
   Hybrid recovers it: batch-build to the known bound, extend
   incrementally past it.
2. Open-ended growth forces a general bigint representation; a known bit
   bound ≤ ~60 allows fixed-width words (large constant-factor reduction).
   Same hybrid: specialize for the known bound, promote on overflow.
With known bounds you can have both; with unknown bounds incrementality
costs only these missed specializations, never redundant work.

## 6b. Design rationale — why each piece is the way it is, with literature

**Why a counting table at all.** Any bijection is forced to be an exact
mixed-radix numeral system over the grammar (§2), so the radixes — the
counts — must be available exactly. Computing them by dynamic programming
over the recurrence is the *recursive method* of combinatorial generation
(Nijenhuis & Wilf, *Combinatorial Algorithms*, 1978; systematised by
Flajolet–Zimmermann–Van Cutsem 1994; Knuth TAOCP 4A on generating trees),
and ranking/unranking through such counts is the standard bijective
machinery (Goldberg & Sipser 1991 frame ranking as optimal compression).

**Why exact big integers, not floats.** Denise & Zimmermann (1999) show
certified floating-point cuts *sampling* to near-linear time — but a
bijection needs exact comparisons at every block boundary: a one-ulp error
decodes a different term. Floats can accelerate the search for the right
block; they cannot replace the exact offsets that must then be subtracted.
Equivalently: the bijection is *exact* arithmetic coding (Rissanen 1976,
Pasco 1976), specialised to the counting measure — which is precisely
Cover's enumerative source coding (1973).

**Why jump tables (the `exploration/` accelerator layer).** Decoding must
locate a rank inside the cumulative app-split distribution of each size
class — a *static predecessor search* over a CDF known entirely at
preprocessing time. Binary search costs O(log n) probes per application
node. Because the distribution is fixed, one can bucket its rank mass into
equal cells, each storing its first candidate block: lookup becomes one
division plus expected O(1) probes (measured 1.4, flat in L). This is the
search-side analogue of Walker's alias method for O(1) static discrete
*sampling* (Walker 1977; Vose 1991) and a degenerate, distribution-aware
case of interpolation search (Peterson 1957; Perl–Itai–Avni 1978). In the
word-RAM model, deterministic O(1) static predecessor structures
(fusion-tree packing, Fredman & Willard 1993) could remove even the
worst case — unnecessary in practice at ≤ L split points per node. The
clean four-language implementations omit this layer deliberately (binary
search inside the walk keeps them readable); the measured prototypes live
in `exploration/`.

**Why prefix-freeness was given up.** A bijection cannot be prefix-free
(`0` and `00` are both codewords; equivalently the Kraft–McMillan sum of
{0,1}* diverges). Self-delimiting codes — Levenshtein's universal code
(1968), Tromp's BLC — buy streaming decode and concatenable programs at
~2.7% redundancy; totality + bijectivity buys exhaustive coverage at the
price of needing the string's length as side information. These are the
plain-vs-prefix (KS vs KP) poles of algorithmic information theory, both
defined concretely in Tromp's *Functional Bits*.

**Why no closed form can replace the table.** Closed λ-terms have no
finite admissible combinatorial specification — the context parameter
nests indefinitely (Bendkowski 2020, §2.5; Bodini–Gittenberger–Gołębiewski
2018) — so no holonomic recurrence in n alone exists, unlike Catalan
numbers for plain binary trees. Capping the de Bruijn index at K makes the
language unambiguous context-free, the generating function algebraic, and
the table Θ(K·n); that is all the cap does.

**Why local/linear-time bijections don't exist here.** Rémy's algorithm
(1985) gives linear-time bijective generation for plain binary trees, and
Bacher–Bodini–Jacquot extended it to Motzkin trees — but the extension to
closed λ-terms fails because a leaf must know how many binders sit above
it (Bendkowski 2020, Remark 8). The binding constraint is exactly what
makes the counting two-dimensional and the local tricks unavailable.

**Why the compressor shares the bijection's walk.** Coding grammar decisions
under learned weights is the Krichevsky–Trofimov/PPM construction restricted
to valid syntax; swapping the counting measure for adaptive weights converts
the zero-redundancy uniform code into a corpus-adaptive one. The same
pre-order walk and depth conditioning drive both; they differ only in the
arithmetic backend, the bijection using exact counting arithmetic and the
compressor a renormalizing range coder. Bijection = uniform measure;
compression = any other measure; one walk serves both.

## 7. Measured scaling curves (Python reference; C++ ≈ 50× faster constants)

Unbounded (canonical) — `profile_scaling.py`, 24 samples/L:

| L | n_max | build s | peak MB | entries | table KiB | decode ms (±std, min–max) | encode ms |
|---|---|---|---|---|---|---|---|
| 16 | 28 | 0.006 | 0.05 | 405 | 0.6 | 0.042 ±0.015, 0.019–0.079 | 0.032 |
| 64 | 81 | 0.167 | 0.21 | 3,320 | 18 | 0.265 ±0.085, 0.138–0.445 | 0.259 |
| 256 | 282 | 7.94 | 2.6 | 39,902 | 841 | 3.53 ±0.94, 1.93–5.60 | 3.51 |

Capped K = 32:

| L | n_max | build s | peak MB | entries | table KiB | decode ms | encode ms |
|---|---|---|---|---|---|---|---|
| 128 | 148 | 0.48 | 0.33 | 4,355 | 38 | 1.12 ±0.26 | 1.15 |
| 512 | 547 | 8.99 | 1.7 | 17,522 | 562 | 15.1 ±3.8 | 14.9 |
| 1024 | 1074 | 40.1 | 4.6 | 34,913 | 2,207 | 55.0 ±18.4, 20.3–102.0 | 53.7 |

Empirical exponents (value ~ L^e) converge to the derived ones: unbounded
build → 2.9 (theory 3), table bits → 2.8 (theory 3); capped build → 2.16
(theory 2 + bigint growth), entries → 1.0·K·L, table bits → 1.97 (theory 2);
decode/encode ≈ 1.8–2.1 (theory: T-only core O(L²)). The min–max spread
(~5×) is term-shape variance: degenerate spines walk more split blocks.
Earlier accelerator measurements: full-table level-0 lookup 197 ns/query
(L ≤ 16); SPLIT+JUMP word-regime decode ≈ 1.9 word-ops/bit with 1.4
probes per app node (L ≤ 60); 202-bit canonical decode 39 ms in Python.

**Compiled-language benchmarks** (C++ and Rust `--bench`, one process per
block; rendered in `plots/`, raw data `plots/bench_results.csv`): at
L = 64 unbounded, C++ builds in ~3 ms and decodes in 10 µs; at L = 1024
unbounded, build 21 s, decode 1.1 ms, encode 0.6 ms, table 43 MiB on
disk, ~80 MB peak RSS; the K = 32 cap at the same length: build 0.6 s,
table 2.1 MiB on disk, ~9 MB peak. Rust tracks C++ closely — decode a touch
faster, build ~40% slower (allocator-heavy bignums), peak RSS comparable.
Table metrics (entries, bits, disk bytes) are bit-identical between the
two — an additional cross-validation.

## 8. AIT positioning (vs Tromp's BLC)

- BLC is injective, prefix-free, streaming O(1)/bit; valid closed programs
  of length n number ~ρⁿ, ρ = 1.963447954 (Grygiel–Lescanne), so its
  asymptotic redundancy is log₂(2/ρ) ≈ 2.7% plus per-term constants.
- The bijection is zero-redundancy by construction (it is exact
  enumerative/arithmetic coding under the counting measure): Ω 18→9 bits,
  S 23→12, Y 30→19, λx.x → empty string.
- What BLC's redundancy buys: Kraft ⇒ Chaitin's Ω converges, universal
  prior, program∥data stream composition — structurally impossible for any
  bijection. The two codes are the plain-vs-prefix (KS vs KP) poles of the
  same theory; under the invariance theorem they define the same K up to
  an additive constant.

## 9. Compression

The rank code (`encode(table, term)`, the bijection's own output) is a
compressor. Per closed λ-term it produces fewer bytes than general-purpose
compressors, for two reasons (both AIT):

1. The rank code is the Shannon-optimal code for the uniform distribution
   over closed terms: on a uniformly random term, no code is shorter. A
   general-purpose compressor does not model the grammar, so it leaves BLC's
   ~2.7% redundancy (rank ≈ log₂ρ·BLC ≈ 0.973·BLC asymptotically; 84% on the
   finite corpus after byte packing) and adds a per-stream container.
2. Per term, gzip/bzip2/lzma add 20–90 bytes of framing; the rank code adds
   none. Per term it is 3–10× smaller than these (measured: uniform-random
   corpus, rank 84% of BLC vs gzip 259%, lzma 557%; structured corpus 83% vs
   343% / 858%).

Counting argument: any bijection {0,1}* ↔ closed terms has the same multiset
of codeword lengths — 2^L strings of length L — hence the same mean length
over all terms. No bijection has a smaller mean than the rank code; it is the
bijection that assigns the smallest terms the shortest codes. BLC is not a
bijection (valid closed programs are an exponentially-vanishing fraction of
strings: 0.084% at length 20), which is its redundancy and why the rank code
is smaller. A bijection smaller than the rank code on a *specific*
non-uniform source exists: reorder the enumeration so that source's frequent
terms get shorter codes — a bijective compressor (a probability-sorted rank
code / bijective arithmetic coding), shorter on the source and longer
elsewhere. The rank code is a bijection that compresses; a source-tuned
bijection compresses further while staying total and invertible.

The structure-aware coder (`lambda_compress`, Python/C++/Rust,
byte-compatible via `--compress-vectors`; no WL port) is an adaptive
renormalizing range coder over the grammar walk: constructor model
conditioned on depth, bucketed index model with bit-coded high indices,
structural termination, a bare byte stream (a four-byte flush tail, no size
header). Both directions are linear in the node count. It consumes a term,
not bytes, so it cannot run on BLC. On small random terms
its adaptive startup costs more than the rank code; as terms grow the model
converges and its output is smaller than raw BLC and than gzip/bzip2/lzma at
every size into the 10⁴–10⁵-node range (the general-purpose compressors show
no discontinuities, near 85% of BLC). The rank code is the smaller codec for
small/medium terms, the grammar coder for large terms.

Grammar-native subterm sharing (`python/compress_research.py`, measured)
adds a repeat-reference layer — the α-invariant analogue of LZ matching. On
the structured corpus it cuts the grammar coder from 112% to 95% of raw BLC
per term, below the LZ-family compressors but, on these small terms, still
above the rank code's 83%; the constant per-term coder overhead amortizes on
larger repetitive terms. (A research script, separate from the shipped coder.)

Plots: `plots/compression_per_term.png` (per term), `…_corpus.png`
(both regimes, both corpora), `…_large.png` (the large-term sweep). Bench:
`python/compression_benchmark.py` (the full BLC/rank × {raw, gzip, bzip2,
lzma} matrix; the x-axis is *node count*, a representation-neutral size, not
BLC bits — which would be circular, since raw BLC is one of the methods).

## 10. Literature map and novelty

**Closest published prior art — Grygiel & Lescanne 2014 (arXiv:1401.0379),
§7 "Unrankings".** It uses the same counting recurrence (their S_{m,n} =
the T(n,m) here, same Tromp size model, same OEIS A114852) and unranking
functions s_{m,n} — bijections from {1,…,S_{m,n}} to terms *of a single
size class n*, with (m, n, k) supplied by the caller, as a half-page Haskell
program that recomputes counts on the fly. What it does **not** have, and
this project adds: the *ranking* direction (term→number) at all (their use
case, random generation, needs only unranking); composition across size
classes into a *single* map; the bijective-numeration step that makes
**every bit string** a term; a global de Bruijn cap (their m is the
free-variable context, the index ceiling K here is a different axis); table
persistence; incrementality; and the engineering. Even the within-class
order differs (their variable block comes last, size-first here), so this
enumeration is a distinct canonical object built with a shared technique.

**Other points:**
- Tromp, *Functional Bits* — BLC, plain (KS) vs prefix (KP) complexity, the
  n ↔ string numeration table (§4.4), prefix codes & Kraft (§4.5).
- Tarau — size-proportionate bijective Gödel numberings of λ-terms
  (compressed de Bruijn ↔ Catalan objects ↔ naturals); the closest work on
  the *both-directions* axis, but in research-grade Prolog over tree-based
  naturals, with no totality-over-bitstrings framing and no engineering.
- Bendkowski, *How to generate random lambda terms?* (arXiv:2005.08856) —
  recursive method (Nijenhuis–Wilf), O(n^{3+o(1)}) bit complexity,
  Denise–Zimmermann certified floats, no finite admissible specification,
  Rémy-style local bijections do **not** extend to closed terms (Remark 8).
- Goldberg & Sipser, *Compression and Ranking* (SICOMP 1991) — ranking as
  optimal compression; the complexity-theoretic frame.
- Code search: `BinaryLambdaCalculus.jl` reimplements GL's counting/unranking
  (per-size-class only, no ranking, caps, persistence, or benchmarks — its
  README says as much); `tromp/AIT` is interpreters, not ranking. No
  grammar-aware λ-term *compressor* exists in any repo or paper I found.
- Barker's **Jot** — prior art for "every bitstring is a program", achieved
  by changing the *language* (a combinator left-fold over the bits), not a
  bijection onto closed λ-terms; the genetic-programming "closure property"
  (every symbol tree is a valid λ-expression) is the same idea at the tree
  level, never as a bitstring↔term bijection.

**Ranking algorithms are also mature prior art.** Ranking/unranking of
*decomposable* combinatorial structures is solved generically: Martínez &
Molinero's generic unranking (Random Structures & Algorithms, 2001), the
Flajolet–Zimmermann–Van Cutsem recursive method, and tree-specific
linear-time results (Catalan-cipher-vector ranking; cool-lex k-ary trees;
B-tree ranking in O(n) after O(n²) preprocessing — the same complexity
profile used here). The jump-table acceleration is predecessor search on a
known CDF (Walker's alias method / interpolation-search territory); it
matches these established bounds; it does not improve on them. So the word-RAM
treatment is competent application, not a new algorithm.

**Assembled here, not found in a single prior source:** implicit-1
numeration ∘ this specific BLC-size-ordered closed-term unranking, with both
directions, as a single total bijection onto closed untyped λ-terms; the
BLC-vs-bijective redundancy measurements; the two-axis incremental table with
cap-independence reuse and the capped-memory regime; and the empirical
comparison of λ-term binary codings against general-purpose compressors
(gzip/bzip2/lzma). The mathematical and algorithmic components are published
and cited above; this project is their engineering synthesis, with
attribution. (Caveat: the ranking-algorithms literature is large; the "not
found in one source" claims rest on a focused, non-exhaustive sweep.)

## 11. Implementations

Four self-contained single-file implementations, one specification,
byte-identical cross-validated output (800 test vectors: 500 unbounded +
300 capped K=3, identical across Python ≡ C++ ≡ Rust ≡ WL; table files
byte-identical and cross-loadable between Python, C++ and Rust):

Each self-test covers: brute-force counts vs exhaustive BLC enumeration
(n≤14 × 4 caps, including open contexts m>0), 12,000 round trips, cap change
in every direction, save/load with checksum-corruption rejection, the
compression round trips, and error paths (bad bits, non-closed/over-cap
terms, malformed compressed bytes).

| | file | status |
|---|---|---|
| Python (reference) | `python/lambda_bijection.py` | self-test ✓; pyflakes-clean |
| C++17 (no deps, own BigNat; GMP swap noted in header) | `cpp/lambda_bijection.cpp` | self-test ✓; vectors ≡ Python; `-Wall -Wextra` clean; `--bench` mode |
| Rust (no deps, own BigNat) | `rust/lambda_bijection.rs` | self-test ✓; vectors ≡ Python; clippy-clean; `--bench` mode |
| Wolfram Language (pure WL, notebook-style .wl) | `wolfram/LambdaBinarization.wl` | executed via WolframScript: `LambdaBijectionSelfTest[]` all True, vectors ≡ Python |

**Robustness.** encode/compress reject non-closed and over-cap terms (so an
out-of-domain term never silently takes another's code); decode rejects
non-0/1 input; table load verifies an FNV-1a checksum and rejects any
corruption; decompress bounds total nodes and nesting depth, so any input
terminates, yielding a term or an error. The
structural recursions are depth-proportional to term nesting; Python raises
its recursion limit to cover multi-kilobit strings, and the C++/Rust
decompressors cap nesting below the native stack. A caller decoding
untrusted, unbounded-length input should bound the input length itself
(decode grows the table with L). `Table` is single-threaded (it extends
lazily during decode); share one per thread.

**Uniform interface** (names adapted per language):
- objects: `Var/Lam/App` (WL: `LambdaVar/LambdaAbs/LambdaApp`), the table
  (`Table(index_cap)` / `Table` class / memoised `TermCount`).
- table: `extend` / `BuildLambdaTable` (size axis), `set_index_cap` /
  `setIndexCap` (cap axis; WL shares across caps automatically via
  normalised memo keys), `count`, `closed_cumulative` /
  `ClosedTermCumulative`.
- codec: `encode` / `EncodeLambdaTerm`, `decode` / `DecodeBitString`.
- printing: `show_term`/`showTerm`/`LambdaTermForm` (λ-rendering),
  `show_bits`/`BitStringForm` (ε for empty), WL also `LambdaTermTree`
  (expression trees).
- workflow helper: `max_de_bruijn_index` / `MaxDeBruijnIndex`.
- persistence: `table.save`/`Table.load` (Python), `saveToFile`/
  `loadFromFile` (C++), `save_to_file`/`load_from_file` (Rust),
  `SaveLambdaTable`/`LoadLambdaTable` (WL) — one shared binary format
  (`.lamtab`: magic, cap, size, length-prefixed big-endian count rows,
  FNV-1a-64 trailer), byte-identical and cross-loadable across all four; WL
  also reads and writes native `.wxf` and `.mx` for the same table.

`--vectors` emits the shared test vectors (Python/C++/Rust);
`--save-table`/`--load-table` exercise cross-language table interop;
`python/check_cross_tables.py` is the Python side of that check.

## 12. File inventory

- `README.md` — usage; `NOTES.md` — this document; `LICENSE` — MIT.
- `python/lambda_bijection.py`, `cpp/lambda_bijection.cpp`,
  `rust/lambda_bijection.rs`, `wolfram/LambdaBinarization.wl` — the four
  implementations.
- `python/lambda_compress.py` — the shipped compression coder (§9);
  `python/compression_benchmark.py` — the BLC/rank × {raw,gzip,bzip2,lzma}
  matrix + plots; `python/compression_large.py` — the large-term sweep;
  `python/compress_research.py` — the subterm-sharing / prior research
  harness (§9).
- `python/profile_scaling.py` — Python-reference scaling (§7);
  `python/plot_performance.py` — C++/Rust benchmark driver and plots
  (`--from-csv` re-renders without re-running); `plots/` — rendered curves +
  `bench_results.csv`.
- `python/check_cross_tables.py`, `python/gen_wl_constants.py` — test
  utilities.
- `exploration/` — research record from the design sessions:
  `bijective_lambda.py` (first validated implementation + AIT
  comparisons), `bounded_o1.py` (level-0 full materialization and
  word-regime SPLIT accelerators), `extension_scaling.py` (incremental
  extension law + jump tables), `jump_table_demo.py` (worked class
  anatomy and traces).

## 13. References

The mathematical components of this project are published; the citations
below are the prior work it builds on and is positioned against (see §6b,
§8–§10). To the best of a non-exhaustive sweep, the *assembled* artifact —
the size-ordered total bijection {0,1}* ↔ closed de Bruijn λ-terms with
both directions, the de Bruijn-cap memory regime, the cross-language
engineering, and the empirical comparison of λ-term binary codings against
general-purpose compressors — has no single published source; the bijection
(Tarau), λ-calculus-as-compression (Kobayashi), and bijective arithmetic
coding (Timmermans) are individually prior art.

**λ-term enumeration, counting, ranking/unranking**

1. K. Grygiel and P. Lescanne. *Counting and generating lambda terms.*
   Journal of Functional Programming 23(5):594–628, 2013.
   doi:10.1017/S0956796813000178.
2. K. Grygiel and P. Lescanne. *Counting and generating terms in the binary
   lambda calculus.* Journal of Functional Programming 25:e24, 2015.
   arXiv:1401.0379; extended version arXiv:1511.05334.
3. P. Tarau. *Ranking/Unranking of Lambda Terms with Compressed de Bruijn
   Indices.* In CICM 2015, LNCS 9150, pp. 118–133. Springer.
   doi:10.1007/978-3-319-20615-8_8.
4. P. Tarau. *A Size-Proportionate Bijective Encoding of Lambda Terms as
   Catalan Objects Endowed with Arithmetic Operations.* In PADL 2016,
   LNCS 9585, pp. 99–116. Springer. doi:10.1007/978-3-319-28228-2_7.
5. P. Tarau. *A Logic Programming Playground for Lambda Terms, Combinators,
   Types and Tree-based Arithmetic Computations.* arXiv:1507.06944, 2015.
6. M. Bendkowski. *How to generate random lambda terms?* arXiv:2005.08856,
   2020.
7. M. Bendkowski, K. Grygiel, P. Lescanne, M. Zaionc. *A Natural Counting of
   Lambda Terms.* In SOFSEM 2016, LNCS 9587, pp. 183–194. Springer.
   arXiv:1506.02367.
8. B. Gittenberger and Z. Gołębiewski. *On the Number of Lambda Terms With
   Prescribed Size of Their De Bruijn Representation.* In STACS 2016, LIPIcs
   47, 40:1–40:13. arXiv:1509.06139.
9. O. Bodini, B. Gittenberger, Z. Gołębiewski. *Enumerating lambda terms by
   weighted length of their De Bruijn representation.* Discrete Applied
   Mathematics 239:45–61, 2018. arXiv:1707.02101.
10. A. Nijenhuis and H. S. Wilf. *Combinatorial Algorithms.* 2nd ed.,
    Academic Press, 1978. (the recursive method)
11. P. Flajolet, P. Zimmermann, B. Van Cutsem. *A calculus for the random
    generation of labelled combinatorial structures.* Theoretical Computer
    Science 132(1–2):1–35, 1994.
12. A. Denise and P. Zimmermann. *Uniform random generation of decomposable
    structures using floating-point arithmetic.* Theoretical Computer
    Science 218(2):233–248, 1999.
13. J.-L. Rémy. *Un procédé itératif de dénombrement d'arbres binaires et son
    application à leur génération aléatoire.* RAIRO Theoretical Informatics
    and Applications 19(2):179–195, 1985.
14. P. Duchon, P. Flajolet, G. Louchard, G. Schaeffer. *Boltzmann Samplers
    for the Random Generation of Combinatorial Structures.* Combinatorics,
    Probability and Computing 13(4–5):577–625, 2004.
15. N. G. de Bruijn. *Lambda calculus notation with nameless dummies.*
    Indagationes Mathematicae 34:381–392, 1972.

**Algorithmic information theory, binary lambda calculus, compression**

16. J. Tromp. *Binary Lambda Calculus and Combinatory Logic.* In Kolmogorov
    Complexity and Applications, Dagstuhl Seminar Proceedings 06051, 2006;
    revised as *Functional Bits*, 2023. https://tromp.github.io/cl/LC.pdf
17. M. Li and P. Vitányi. *An Introduction to Kolmogorov Complexity and Its
    Applications.* 3rd ed., Springer, 2008.
18. N. Kobayashi, K. Matsuda, A. Shinohara. *Functional Programs as
    Compressed Data.* In PEPM 2012; Higher-Order and Symbolic Computation,
    2013. doi:10.1007/s10990-013-9093-z.
19. A. V. Goldberg and M. Sipser. *Compression and Ranking.* SIAM Journal on
    Computing 20(3):524–536, 1991. (STOC 1985)
20. T. M. Cover. *Enumerative source encoding.* IEEE Transactions on
    Information Theory 19(1):73–77, 1973.
21. J. Rissanen. *Generalized Kraft inequality and arithmetic coding.* IBM
    Journal of Research and Development 20(3):198–203, 1976.
22. M. Timmermans. *Bijective Arithmetic Encoding with Really Good End
    Treatment.* 1999. http://www3.sympatico.ca/mt0000/biacode/biacode.html
23. V. I. Levenshtein. *On the redundancy and delay of decodable coding of
    natural numbers.* Problems of Cybernetics 20:173–179, 1968.
24. C. Barker. *Iota and Jot.* 2001. https://esolangs.org/wiki/Iota_and_Jot

**Discrete sampling and word-RAM search (the accelerator layer, §6b)**

25. A. J. Walker. *An Efficient Method for Generating Discrete Random
    Variables with General Distributions.* ACM Transactions on Mathematical
    Software 3(3):253–256, 1977. (the alias method)
26. M. D. Vose. *A Linear Algorithm for Generating Random Numbers with a
    Given Distribution.* IEEE Transactions on Software Engineering
    17(9):972–975, 1991.
27. M. L. Fredman and D. E. Willard. *Surpassing the Information Theoretic
    Bound with Fusion Trees.* Journal of Computer and System Sciences
    47(3):424–436, 1993.
28. Y. Perl, A. Itai, H. Avni. *Interpolation search—a log log N search.*
    Communications of the ACM 21(7):550–553, 1978.
29. D. E. Knuth. *The Art of Computer Programming, Vol. 4A: Combinatorial
    Algorithms, Part 1.* Addison-Wesley, 2011. (generating all trees)
30. A. Bacher, O. Bodini, A. Jacquot. *Exact-size Sampling for Motzkin Trees
    in Linear Time via Boltzmann Samplers and Holonomic Specification.* In
    ANALCO 2013, pp. 52–61. SIAM.

**Ranking/unranking algorithms (the maturity of the area, §10)**

31. C. Martínez and X. Molinero. *A generic approach for the unranking of
    labeled combinatorial classes.* Random Structures & Algorithms
    19(3–4):472–497, 2001.
32. W. Myrvold and F. Ruskey. *Ranking and unranking permutations in linear
    time.* Information Processing Letters 79(6):281–284, 2001.
33. M. C. Er. *Enumerating ordered trees lexicographically.* The Computer
    Journal 28(5):538–542, 1985. (and the broader Catalan-object
    ranking/unranking line, incl. cool-lex order for k-ary trees)

**Coding and search (named in §6b)**

34. R. C. Pasco. *Source Coding Algorithms for Fast Data Compression.* PhD
    thesis, Stanford University, 1976.
35. R. E. Krichevsky and V. K. Trofimov. *The performance of universal
    encoding.* IEEE Transactions on Information Theory 27(2):199–207, 1981.
36. W. W. Peterson. *Addressing for random-access storage.* IBM Journal of
    Research and Development 1(2):130–146, 1957.
