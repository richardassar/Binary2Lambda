// Binary2Lambda — lambda_bijection.rs: self-contained single-file Rust
// implementation.
//
// Bijection between binary strings and closed untyped lambda terms.
// Every binary string (including the empty string) denotes exactly one
// closed lambda term and vice versa; encode and decode share one
// incrementally-built counting table.
//
// Specification (identical to the Python, C++ and Wolfram Language
// versions):
//   terms        Var(i), i >= 1 (de Bruijn) | Lam(body) | App(fun, arg)
//   size         |Var i| = i+1,  |Lam b| = |b|+2,  |App f a| = |f|+|a|+2
//   order        ascending size; within a size class Var < Lam < App,
//                abstractions by body rank, applications by
//                (left size, left rank, right rank)
//   numeration   string s  <->  N = value of "1"+s in binary, minus 1
//
// A table built with a finite de Bruijn index cap K enumerates the
// sublanguage of closed terms whose indices never exceed K. Each cap value
// (including "unbounded") defines a DIFFERENT bijection; encode and decode
// must use the same cap. Capped and unbounded bijections agree on all terms
// of size <= K+1. Tables can be saved to and loaded from disk in a portable
// text format shared with the Python and C++ implementations.
//
// Arbitrary-precision arithmetic is provided by the self-contained BigNat
// type below, keeping this file dependency-free; for heavy workloads its
// operations map directly onto a bignum crate such as `rug` or `num-bigint`.
//
// The file also contains the project's compression layer: a renormalizing
// range coder over the term grammar with an adaptive model, byte-compatible
// with python/lambda_compress.py. Compression is a separate axis from the
// bijection - it needs compress + decompress, not "every bit string is a
// term".
//
// Build:  rustc -O -o lambda_bijection_rs lambda_bijection.rs
// Run:    ./lambda_bijection_rs                     (self-test and demo)
//         ./lambda_bijection_rs --vectors           (cross-language vectors)
//         ./lambda_bijection_rs --compress-vectors  (compression vectors)
//         ./lambda_bijection_rs --save-table <path> (write a sample table)
//         ./lambda_bijection_rs --load-table <path> (load and verify a table)
//         ./lambda_bijection_rs --bench <cap> <L>   (one benchmark block, CSV)

use std::cmp::Ordering;
use std::fmt::Write as _;
use std::fs;

// ------------------------------------------------------------------ numbers

/// Arbitrary-precision unsigned integer, little-endian base-2^64 limbs,
/// trimmed (no trailing zero limbs). Only the operations the bijection
/// needs: addition, subtraction (no underflow), multiplication, divmod,
/// comparison, and binary/hex string conversion.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct BigNat {
    limbs: Vec<u64>,
}

impl BigNat {
    fn from_u64(value: u64) -> BigNat {
        if value == 0 {
            BigNat::default()
        } else {
            BigNat { limbs: vec![value] }
        }
    }

    fn is_zero(&self) -> bool {
        self.limbs.is_empty()
    }

    fn trim(&mut self) {
        while self.limbs.last() == Some(&0) {
            self.limbs.pop();
        }
    }

    fn add_assign(&mut self, other: &BigNat) {
        if self.limbs.len() < other.limbs.len() {
            self.limbs.resize(other.limbs.len(), 0);
        }
        let mut carry: u128 = 0;
        for i in 0..self.limbs.len() {
            let sum = carry
                + self.limbs[i] as u128
                + *other.limbs.get(i).unwrap_or(&0) as u128;
            self.limbs[i] = sum as u64;
            carry = sum >> 64;
        }
        if carry != 0 {
            self.limbs.push(carry as u64);
        }
    }

    /// Requires self >= other.
    fn sub(&self, other: &BigNat) -> BigNat {
        debug_assert!(*other <= *self, "BigNat subtraction would underflow");
        let mut result = self.clone();
        let mut borrow: u128 = 0;
        for i in 0..result.limbs.len() {
            let take = borrow + *other.limbs.get(i).unwrap_or(&0) as u128;
            let current = result.limbs[i] as u128;
            if take <= current {
                result.limbs[i] = (current - take) as u64;
                borrow = 0;
            } else {
                result.limbs[i] = ((1u128 << 64) + current - take) as u64;
                borrow = 1;
            }
        }
        result.trim();
        result
    }

    fn mul(&self, other: &BigNat) -> BigNat {
        if self.is_zero() || other.is_zero() {
            return BigNat::default();
        }
        let mut out = vec![0u64; self.limbs.len() + other.limbs.len()];
        for i in 0..self.limbs.len() {
            let mut carry: u128 = 0;
            for j in 0..other.limbs.len() {
                let current = self.limbs[i] as u128 * other.limbs[j] as u128
                    + out[i + j] as u128
                    + carry;
                out[i + j] = current as u64;
                carry = current >> 64;
            }
            out[i + other.limbs.len()] = carry as u64;
        }
        let mut result = BigNat { limbs: out };
        result.trim();
        result
    }

    fn divmod(num: &BigNat, den: &BigNat) -> (BigNat, BigNat) {
        assert!(!den.is_zero());
        let mut quotient = BigNat::default();
        let mut remainder = BigNat::default();
        for i in (0..num.bit_length()).rev() {
            remainder.shift_left_one();
            if num.bit(i) {
                remainder.set_bit(0);
            }
            if *den <= remainder {
                remainder = remainder.sub(den);
                quotient.set_bit(i);
            }
        }
        (quotient, remainder)
    }

    fn bit_length(&self) -> usize {
        match self.limbs.last() {
            None => 0,
            Some(top) => self.limbs.len() * 64 - top.leading_zeros() as usize,
        }
    }

    fn bit(&self, i: usize) -> bool {
        self.limbs
            .get(i / 64)
            .is_some_and(|limb| (limb >> (i % 64)) & 1 != 0)
    }

    fn set_bit(&mut self, i: usize) {
        if i / 64 >= self.limbs.len() {
            self.limbs.resize(i / 64 + 1, 0);
        }
        self.limbs[i / 64] |= 1u64 << (i % 64);
    }

    fn shift_left_one(&mut self) {
        let mut carry = 0u64;
        for limb in &mut self.limbs {
            let next = *limb >> 63;
            *limb = (*limb << 1) | carry;
            carry = next;
        }
        if carry != 0 {
            self.limbs.push(carry);
        }
    }

    /// Most significant bit first; "0" for zero.
    fn to_binary_string(&self) -> String {
        if self.is_zero() {
            return "0".to_string();
        }
        (0..self.bit_length())
            .rev()
            .map(|i| if self.bit(i) { '1' } else { '0' })
            .collect()
    }

    fn from_binary_string(bits: &str) -> BigNat {
        let mut result = BigNat::default();
        for c in bits.chars() {
            result.shift_left_one();
            if c == '1' {
                result.set_bit(0);
            }
        }
        result
    }

    /// Lowercase, no leading zeros.
    fn to_hex_string(&self) -> String {
        match self.limbs.split_last() {
            None => "0".to_string(),
            Some((top, rest)) => {
                let mut out = format!("{:x}", top);
                for limb in rest.iter().rev() {
                    write!(out, "{:016x}", limb).unwrap();
                }
                out
            }
        }
    }

    fn from_hex_string(hex: &str) -> BigNat {
        let mut result = BigNat::default();
        for c in hex.chars() {
            let digit = c.to_digit(16).expect("invalid hex digit") as u64;
            for _ in 0..4 {
                result.shift_left_one();
            }
            if digit != 0 {
                if result.limbs.is_empty() {
                    result.limbs.push(0);
                }
                result.limbs[0] |= digit;
            }
        }
        result
    }
}

impl Ord for BigNat {
    fn cmp(&self, other: &Self) -> Ordering {
        self.limbs
            .len()
            .cmp(&other.limbs.len())
            .then_with(|| self.limbs.iter().rev().cmp(other.limbs.iter().rev()))
    }
}

impl PartialOrd for BigNat {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// -------------------------------------------------------------------- terms

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Term {
    Var(u32),
    Lam(Box<Term>),
    App(Box<Term>, Box<Term>),
}

pub fn term_size(term: &Term) -> usize {
    match term {
        Term::Var(i) => *i as usize + 1,
        Term::Lam(body) => term_size(body) + 2,
        Term::App(fun, arg) => term_size(fun) + term_size(arg) + 2,
    }
}

pub fn max_de_bruijn_index(term: &Term) -> u32 {
    match term {
        Term::Var(i) => *i,
        Term::Lam(body) => max_de_bruijn_index(body),
        Term::App(fun, arg) => {
            max_de_bruijn_index(fun).max(max_de_bruijn_index(arg))
        }
    }
}

/// Err unless term is a well-formed closed term: every de Bruijn index is >= 1
/// and bound by an enclosing lambda (index <= depth). encode and compress
/// require this: the bijection ranges over closed terms, so only a closed term
/// has a code.
pub fn check_closed(term: &Term, depth: u32) -> Result<(), String> {
    match term {
        Term::Var(i) if *i == 0 => Err("de Bruijn index must be >= 1".to_string()),
        Term::Var(i) if *i > depth => {
            Err("term is not closed (free variable)".to_string())
        }
        Term::Var(_) => Ok(()),
        Term::Lam(body) => check_closed(body, depth + 1),
        Term::App(fun, arg) => {
            check_closed(fun, depth)?;
            check_closed(arg, depth)
        }
    }
}

/// Fold data into a running FNV-1a 64-bit hash, used as a table-file integrity
/// checksum; identical integer arithmetic in the Python and C++ implementations,
/// so a saved table's checksum is the same in every language. Folding
/// incrementally lets save and load stream the file without holding the whole
/// body in memory.
const FNV_OFFSET: u64 = 0xcbf29ce484222325;

fn fnv1a64_update(mut h: u64, data: &[u8]) -> u64 {
    for &byte in data {
        h ^= byte as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Magic for the portable binary table format (.lamtab); the last byte is the
/// format version. Shared byte-for-byte with the Python, C++ and Wolfram ports.
const TABLE_MAGIC: &[u8] = b"LAMTAB\x01";

fn put_u32(out: &mut Vec<u8>, v: u32) {
    out.extend_from_slice(&v.to_le_bytes());
}
fn get_u32(d: &[u8], pos: usize) -> u32 {
    u32::from_le_bytes([d[pos], d[pos + 1], d[pos + 2], d[pos + 3]])
}
fn to_bytes_be(v: &BigNat) -> Vec<u8> {
    let nbytes = (v.bit_length() + 7) / 8;
    let mut out = vec![0u8; nbytes];
    for k in 0..nbytes {
        // k counts bytes from the least significant
        let mut byte = 0u8;
        for j in 0..8 {
            if v.bit(8 * k + j) {
                byte |= 1 << j;
            }
        }
        out[nbytes - 1 - k] = byte;
    }
    out
}
fn from_bytes_be(d: &[u8]) -> BigNat {
    let mut v = BigNat::default();
    for &byte in d {
        for _ in 0..8 {
            v.shift_left_one();
        }
        for b in (0..8).rev() {
            if (byte >> b) & 1 != 0 {
                v.set_bit(b);
            }
        }
    }
    v
}

/// context: 0 = top level / lambda body, 1 = left of app, 2 = right of app
pub fn show_term(term: &Term, context: u8) -> String {
    match term {
        Term::Var(i) => i.to_string(),
        Term::Lam(body) => {
            let rendered = format!("λ{}", show_term(body, 0));
            if context > 0 {
                format!("({})", rendered)
            } else {
                rendered
            }
        }
        Term::App(fun, arg) => {
            let rendered =
                format!("{} {}", show_term(fun, 1), show_term(arg, 2));
            if context == 2 {
                format!("({})", rendered)
            } else {
                rendered
            }
        }
    }
}

pub fn show_bits(bits: &str) -> &str {
    if bits.is_empty() {
        "ε"
    } else {
        bits
    }
}

// -------------------------------------------------------------------- table

/// Counting table T(n, m): number of terms of size n whose free indices are
/// all <= m and whose indices nowhere exceed the cap. Rows are stored under
/// the effective context min(m, n-1, cap). Size extension is append-only;
/// changing the cap reuses every row of size <= min(old, new) + 1 (those
/// are cap-independent) and rebuilds only the rest.
pub struct Table {
    cap: Option<usize>,
    rows: Vec<Vec<BigNat>>,
    cum: Vec<BigNat>, // closed terms of size <= n
    zero: BigNat,
}

impl Table {
    pub fn new(index_cap: Option<usize>) -> Table {
        if let Some(cap) = index_cap {
            assert!(cap >= 1, "index cap must be at least 1");
        }
        Table {
            cap: index_cap,
            rows: vec![Vec::new(), Vec::new()],
            cum: vec![BigNat::default(), BigNat::default()],
            zero: BigNat::default(),
        }
    }

    pub fn built_size(&self) -> usize {
        self.rows.len() - 1
    }

    fn effective_context(&self, n: usize, m: usize) -> usize {
        let mut m = m.min(n - 1);
        if let Some(cap) = self.cap {
            m = m.min(cap);
        }
        m
    }

    fn width(&self, n: usize) -> usize {
        self.effective_context(n, n - 1) + 1
    }

    /// T(n, m) for already-built sizes (n <= built_size); m is any context.
    fn count_built(&self, n: usize, m: usize) -> &BigNat {
        if n < 2 {
            &self.zero
        } else {
            &self.rows[n][self.effective_context(n, m)]
        }
    }

    fn cumulative_built(&self, n: usize) -> &BigNat {
        if n < 2 {
            &self.zero
        } else {
            &self.cum[n]
        }
    }

    /// T(n, m) computed from the rows below n; m is an effective context.
    /// Only ever called with n >= 2 (extend starts at rows.len() >= 2; load
    /// validates built >= 2), so n - 2 never underflows.
    fn row_value(&self, n: usize, m: usize) -> BigNat {
        let mut value = BigNat::from_u64(if m == n - 1 { 1 } else { 0 });
        value.add_assign(self.count_built(n - 2, m + 1));
        if n >= 6 {
            for k in 2..=(n - 4) {
                let block =
                    self.count_built(k, m).mul(self.count_built(n - 2 - k, m));
                value.add_assign(&block);
            }
        }
        value
    }

    /// Append rows up to size_limit; existing rows are never touched.
    pub fn extend(&mut self, size_limit: usize) {
        for n in self.rows.len()..=size_limit {
            let row: Vec<BigNat> =
                (0..self.width(n)).map(|m| self.row_value(n, m)).collect();
            let mut cumulative = self.cum.last().unwrap().clone();
            if let Some(first) = row.first() {
                cumulative.add_assign(first);
            }
            self.rows.push(row);
            self.cum.push(cumulative);
        }
    }

    /// Change the de Bruijn cap, reusing all cap-independent rows: rows of
    /// size n are identical under caps K and K' whenever n-1 <= min(K, K').
    pub fn set_index_cap(&mut self, new_cap: Option<usize>) {
        if let Some(cap) = new_cap {
            assert!(cap >= 1, "index cap must be at least 1");
        }
        if new_cap == self.cap {
            return;
        }
        let lower = match (self.cap, new_cap) {
            (Some(a), Some(b)) => a.min(b),
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => unreachable!(),
        };
        let keep = lower + 1;
        let target = self.built_size();
        self.cap = new_cap;
        if self.built_size() > keep {
            self.rows.truncate(keep + 1);
            self.cum.truncate(keep + 1);
        }
        self.extend(target);
    }

    /// Portable binary format shared with the Python, C++ and Wolfram ports
    /// (identical bytes in every language): the 7-byte magic LAMTAB\x01, a
    /// one-byte cap kind (0 unbounded, 1 finite) and 4-byte little-endian cap
    /// value, the 4-byte little-endian built size, then for each size from 2
    /// upward its row of counts, each a 4-byte little-endian length and that
    /// many big-endian magnitude bytes, then an 8-byte little-endian FNV-1a-64
    /// of every preceding byte. Cumulative counts are derivable, so they are
    /// not stored. The body is streamed, so memory stays flat.
    pub fn save_to_file(&self, path: &str) -> std::io::Result<()> {
        use std::io::Write as _;
        let mut out = std::io::BufWriter::new(std::fs::File::create(path)?);
        let mut checksum = FNV_OFFSET;

        let mut header = Vec::new();
        header.extend_from_slice(TABLE_MAGIC);
        header.push(if self.cap.is_some() { 1 } else { 0 });
        put_u32(&mut header, self.cap.map_or(0, |c| c as u32));
        put_u32(&mut header, self.built_size() as u32);
        out.write_all(&header)?;
        checksum = fnv1a64_update(checksum, &header);

        for n in 2..=self.built_size() {
            for v in &self.rows[n] {
                let mag = to_bytes_be(v);
                let mut cell = Vec::new();
                put_u32(&mut cell, mag.len() as u32);
                cell.extend_from_slice(&mag);
                out.write_all(&cell)?;
                checksum = fnv1a64_update(checksum, &cell);
            }
        }
        out.write_all(&checksum.to_le_bytes())?;
        out.flush()
    }

    /// Read a table written by save_to_file (any implementation). The checksum
    /// is verified, so any corruption or truncation is rejected.
    pub fn load_from_file(path: &str) -> Result<Table, String> {
        let data = fs::read(path).map_err(|e| e.to_string())?;
        let header_size = TABLE_MAGIC.len() + 1 + 4 + 4;
        if data.len() < header_size + 8 || &data[..TABLE_MAGIC.len()] != TABLE_MAGIC {
            return Err("not a lambda-binarization table file".to_string());
        }
        let end = data.len() - 8;
        let mut tail = [0u8; 8];
        tail.copy_from_slice(&data[end..]);
        let stored = u64::from_le_bytes(tail);
        if fnv1a64_update(FNV_OFFSET, &data[..end]) != stored {
            return Err("table file failed checksum (corrupt or truncated)".to_string());
        }
        let cap_kind = data[TABLE_MAGIC.len()];
        let cap_value = get_u32(&data, TABLE_MAGIC.len() + 1);
        let built = get_u32(&data, TABLE_MAGIC.len() + 5) as usize;
        if cap_kind > 1 || built < 1 {
            return Err("table file has a malformed header".to_string());
        }
        let cap = if cap_kind == 0 { None } else { Some(cap_value as usize) };
        let mut table = Table::new(cap);
        let mut pos = header_size;
        for n in 2..=built {
            let mut row = Vec::new();
            for _ in 0..table.width(n) {
                if pos + 4 > end {
                    return Err("table file is truncated".to_string());
                }
                let nbytes = get_u32(&data, pos) as usize;
                pos += 4;
                if pos + nbytes > end {
                    return Err("table file is truncated".to_string());
                }
                row.push(from_bytes_be(&data[pos..pos + nbytes]));
                pos += nbytes;
            }
            let mut cumulative = table.cum.last().unwrap().clone();
            if let Some(first) = row.first() {
                cumulative.add_assign(first);
            }
            table.rows.push(row);
            table.cum.push(cumulative);
        }
        if pos != end {
            return Err("table file size does not match its rows".to_string());
        }
        Ok(table)
    }
}

// -------------------------------------------------------------------- codec

fn var_count(table: &Table, n: usize, m: usize) -> u64 {
    if n < 2 || n - 1 > m {
        return 0;
    }
    if let Some(cap) = table.cap {
        if n - 1 > cap {
            return 0;
        }
    }
    1
}

/// Rank of the term within its size class, in context m. The table must
/// already be built up to the term's size.
fn rank_of(table: &Table, term: &Term, m: usize) -> BigNat {
    match term {
        Term::Var(_) => BigNat::default(),
        Term::Lam(body) => {
            let n = term_size(term);
            let mut rank = BigNat::from_u64(var_count(table, n, m));
            rank.add_assign(&rank_of(table, body, m + 1));
            rank
        }
        Term::App(fun, arg) => {
            let n = term_size(term);
            let mut rank = BigNat::from_u64(var_count(table, n, m));
            rank.add_assign(table.count_built(n - 2, m + 1));
            let left_size = term_size(fun);
            for k in 2..left_size {
                let block =
                    table.count_built(k, m).mul(table.count_built(n - 2 - k, m));
                rank.add_assign(&block);
            }
            let left = rank_of(table, fun, m)
                .mul(table.count_built(term_size(arg), m));
            rank.add_assign(&left);
            rank.add_assign(&rank_of(table, arg, m));
            rank
        }
    }
}

/// The rank-th term of size n in context m; total for 0 <= rank < T(n, m).
/// The divmod split below is exact only because rank < count(k,m)*right, so
/// the quotient stays below count(k,m).
fn unrank(table: &Table, mut rank: BigNat, n: usize, m: usize) -> Term {
    if var_count(table, n, m) == 1 {
        if rank.is_zero() {
            return Term::Var((n - 1) as u32);
        }
        rank = rank.sub(&BigNat::from_u64(1));
    }
    let lam_block = table.count_built(n - 2, m + 1);
    if rank < *lam_block {
        return Term::Lam(Box::new(unrank(table, rank, n - 2, m + 1)));
    }
    rank = rank.sub(lam_block);
    if n >= 6 {
        for k in 2..=(n - 4) {
            let right = table.count_built(n - 2 - k, m);
            let block = table.count_built(k, m).mul(right);
            if rank < block {
                let (left_rank, right_rank) = BigNat::divmod(&rank, right);
                return Term::App(
                    Box::new(unrank(table, left_rank, k, m)),
                    Box::new(unrank(table, right_rank, n - 2 - k, m)),
                );
            }
            rank = rank.sub(&block);
        }
    }
    panic!("rank out of range for size class");
}

/// Closed lambda term -> binary string (inverse of decode). Errors if term is
/// not closed or uses an index above the cap.
pub fn encode(table: &mut Table, term: &Term) -> Result<String, String> {
    check_closed(term, 0)?;
    if let Some(cap) = table.cap {
        if max_de_bruijn_index(term) as usize > cap {
            return Err(format!("term uses indices above the table cap {}", cap));
        }
    }
    let n = term_size(term);
    table.extend(n);
    let mut number = table.cumulative_built(n - 1).clone();
    number.add_assign(&rank_of(table, term, 0));
    number.add_assign(&BigNat::from_u64(1));
    Ok(number.to_binary_string()[1..].to_string())
}

/// Binary string -> closed lambda term (total on {0,1}*).
pub fn decode(table: &mut Table, bits: &str) -> Result<Term, String> {
    if !bits.chars().all(|c| c == '0' || c == '1') {
        return Err("input must consist of 0s and 1s".to_string());
    }
    let number =
        BigNat::from_binary_string(&format!("1{}", bits)).sub(&BigNat::from_u64(1));
    let mut n = 4; // the smallest closed term, λ1, has size 4
    loop {
        table.extend(n);
        if *table.cumulative_built(n) > number {
            break;
        }
        n += 1;
    }
    let rank = number.sub(table.cumulative_built(n - 1));
    Ok(unrank(table, rank, n, 0))
}

// -------------------------------------------------------------- compression
//
// Lambda-specific compression of closed terms - the project's separate
// compression axis. A renormalizing range coder over the term grammar with an
// adaptive model: the coder walks the term in pre-order, the choice set at
// each node is {Var(1..m), Lam, App} with m the lambda depth, and termination
// is structural, so no counting table and no size header are needed. Both
// directions run in time linear in the node count. Byte-compatible with
// python/lambda_compress.py and the C++ port.

const KIND_VAR: usize = 0;
const KIND_LAM: usize = 1;
const KIND_APP: usize = 2;

// 32-bit Subbotin range coder. Registers are u32 and arithmetic wraps modulo
// 2^32 (via the wrapping operators), which makes the byte stream identical
// across the Python, C++ and Rust implementations. Every symbol's total
// frequency stays below RC_BOT; the model guarantees that by rescaling, and
// indices above the bucket count are coded bit by bit.
const RC_TOP: u32 = 1 << 24;
const RC_BOT: u32 = 1 << 16;

/// Narrows [low, low+range) per symbol, emitting settled high bytes.
struct RangeEncoder {
    low: u32,
    range: u32,
    out: Vec<u8>,
}

impl RangeEncoder {
    fn new() -> RangeEncoder {
        RangeEncoder { low: 0, range: 0xFFFF_FFFF, out: Vec::new() }
    }

    fn encode(&mut self, c_low: u32, freq: u32, total: u32) {
        self.range /= total;
        self.low = self.low.wrapping_add(c_low.wrapping_mul(self.range));
        self.range = self.range.wrapping_mul(freq);
        loop {
            let renorm = if (self.low ^ self.low.wrapping_add(self.range)) < RC_TOP {
                true
            } else if self.range < RC_BOT {
                self.range = self.low.wrapping_neg() & (RC_BOT - 1);
                true
            } else {
                false
            };
            if !renorm {
                break;
            }
            self.out.push((self.low >> 24) as u8);
            self.low <<= 8;
            self.range <<= 8;
        }
    }

    fn finish(mut self) -> Vec<u8> {
        for _ in 0..4 {
            self.out.push((self.low >> 24) as u8);
            self.low <<= 8;
        }
        self.out
    }
}

/// Mirrors the encoder; the code register steers symbol choice.
struct RangeDecoder<'a> {
    data: &'a [u8],
    pos: usize,
    low: u32,
    range: u32,
    code: u32,
}

impl<'a> RangeDecoder<'a> {
    fn new(data: &'a [u8]) -> RangeDecoder<'a> {
        let mut decoder = RangeDecoder { data, pos: 0, low: 0, range: 0xFFFF_FFFF, code: 0 };
        for _ in 0..4 {
            decoder.code = (decoder.code << 8) | decoder.next_byte() as u32;
        }
        decoder
    }

    fn next_byte(&mut self) -> u8 {
        // past the stream the decoder reads zero bytes
        if self.pos < self.data.len() {
            let byte = self.data[self.pos];
            self.pos += 1;
            byte
        } else {
            0
        }
    }

    fn target(&mut self, total: u32) -> u32 {
        self.range /= total;
        let value = self.code.wrapping_sub(self.low) / self.range;
        if value >= total {
            total - 1
        } else {
            value
        }
    }

    fn consume(&mut self, c_low: u32, freq: u32) {
        self.low = self.low.wrapping_add(c_low.wrapping_mul(self.range));
        self.range = self.range.wrapping_mul(freq);
        loop {
            let renorm = if (self.low ^ self.low.wrapping_add(self.range)) < RC_TOP {
                true
            } else if self.range < RC_BOT {
                self.range = self.low.wrapping_neg() & (RC_BOT - 1);
                true
            } else {
                false
            };
            if !renorm {
                break;
            }
            self.code = (self.code << 8) | self.next_byte() as u32;
            self.low <<= 8;
            self.range <<= 8;
        }
    }
}

/// Counts shared by compress and decompress, updated identically. Kinds are
/// conditioned on the lambda-depth bucket min(m, 3). Variable indices share an
/// 8-entry bucketed table (bucket min(i, 8)); an index in the tail bucket is
/// followed by its offset, coded bit by bit. Counts are halved when a context
/// total reaches RESCALE_LIMIT, keeping every total below the coder's RC_BOT.
struct AdaptiveModel {
    kinds: [[u32; 3]; 4],
    indices: [u32; 8],
}

const MODEL_INCREMENT: u32 = 16;
const INDEX_BUCKETS: usize = 8;
const RESCALE_LIMIT: u32 = 1 << 14;

impl AdaptiveModel {
    fn new() -> AdaptiveModel {
        AdaptiveModel {
            kinds: [[1; 3]; 4],
            indices: [1; 8],
        }
    }

    fn kind_weights(&self, m: usize) -> [u32; 3] {
        let row = &self.kinds[m.min(3)];
        [if m >= 1 { row[0] } else { 0 }, row[1], row[2]]
    }

    fn saw_kind(&mut self, m: usize, kind: usize) {
        let row = &mut self.kinds[m.min(3)];
        row[kind] += MODEL_INCREMENT;
        if row[0] + row[1] + row[2] >= RESCALE_LIMIT {
            for weight in row.iter_mut() {
                *weight = (*weight >> 1).max(1);
            }
        }
    }

    fn saw_index(&mut self, index: usize) {
        self.indices[(index - 1).min(INDEX_BUCKETS - 1)] += MODEL_INCREMENT;
        if self.indices.iter().sum::<u32>() >= RESCALE_LIMIT {
            for weight in self.indices.iter_mut() {
                *weight = (*weight >> 1).max(1);
            }
        }
    }
}

fn encode_symbol(coder: &mut RangeEncoder, weights: &[u32], symbol: usize) {
    let c_low: u32 = weights[..symbol].iter().sum();
    let total: u32 = weights.iter().sum();
    coder.encode(c_low, weights[symbol], total);
}

fn decode_symbol(coder: &mut RangeDecoder, weights: &[u32]) -> Result<usize, String> {
    let total: u32 = weights.iter().sum();
    let target = coder.target(total);
    let mut c_low = 0u32;
    for (symbol, &weight) in weights.iter().enumerate() {
        if target < c_low + weight {
            coder.consume(c_low, weight);
            return Ok(symbol);
        }
        c_low += weight;
    }
    Err("malformed compressed data (code point out of range)".to_string())
}

fn encode_bit(coder: &mut RangeEncoder, bit: u32) {
    coder.encode(bit, 1, 2);
}

fn decode_bit(coder: &mut RangeDecoder) -> u32 {
    let bit = if coder.target(2) >= 1 { 1 } else { 0 };
    coder.consume(bit, 1);
    bit
}

// Elias-gamma over an equiprobable bit model for the rare de Bruijn index
// above the bucket count; the unary length is capped so malformed input cannot
// spin in the prefix.
const GAMMA_MAX_BITS: usize = 40;

fn encode_gamma(coder: &mut RangeEncoder, value: u32) {
    let v = value + 1;
    let n = (31 - v.leading_zeros()) as usize; // bit length of v, minus one
    for _ in 0..n {
        encode_bit(coder, 0);
    }
    encode_bit(coder, 1);
    for i in (0..n).rev() {
        encode_bit(coder, (v >> i) & 1);
    }
}

fn decode_gamma(coder: &mut RangeDecoder) -> Result<u32, String> {
    let mut n = 0usize;
    while decode_bit(coder) == 0 {
        n += 1;
        if n > GAMMA_MAX_BITS {
            return Err("malformed compressed data (index code too long)".to_string());
        }
    }
    let mut v = 1u32;
    for _ in 0..n {
        v = (v << 1) | decode_bit(coder);
    }
    Ok(v - 1)
}

fn encode_index(coder: &mut RangeEncoder, model: &mut AdaptiveModel, m: usize, index: usize) {
    let alphabet = m.min(INDEX_BUCKETS);
    let bucket = (index - 1).min(INDEX_BUCKETS - 1);
    encode_symbol(coder, &model.indices[..alphabet], bucket);
    if bucket == INDEX_BUCKETS - 1 {
        encode_gamma(coder, (index - INDEX_BUCKETS) as u32);
    }
    model.saw_index(index);
}

fn decode_index(
    coder: &mut RangeDecoder,
    model: &mut AdaptiveModel,
    m: usize,
) -> Result<usize, String> {
    let alphabet = m.min(INDEX_BUCKETS);
    let bucket = decode_symbol(coder, &model.indices[..alphabet])?;
    let index = if bucket < INDEX_BUCKETS - 1 {
        bucket + 1
    } else {
        let candidate = INDEX_BUCKETS + decode_gamma(coder)? as usize;
        if candidate > m {
            return Err("malformed compressed data (index exceeds depth)".to_string());
        }
        candidate
    };
    model.saw_index(index);
    Ok(index)
}

fn walk_encode(term: &Term, m: usize, model: &mut AdaptiveModel, coder: &mut RangeEncoder) {
    let weights = model.kind_weights(m);
    match term {
        Term::Var(index) => {
            encode_symbol(coder, &weights, KIND_VAR);
            model.saw_kind(m, KIND_VAR);
            encode_index(coder, model, m, *index as usize);
        }
        Term::Lam(body) => {
            encode_symbol(coder, &weights, KIND_LAM);
            model.saw_kind(m, KIND_LAM);
            walk_encode(body, m + 1, model, coder);
        }
        Term::App(fun, arg) => {
            encode_symbol(coder, &weights, KIND_APP);
            model.saw_kind(m, KIND_APP);
            walk_encode(fun, m, model, coder);
            walk_encode(arg, m, model, coder);
        }
    }
}

/// Bound on the term's nesting depth (App-spines deepen without raising the
/// lambda context m), kept well under every implementation's stack: a stream
/// that keeps nesting hits this limit and errors while the stack is still safe.
/// Real lambda terms stay far below it.
const MAX_DECODE_DEPTH: usize = 12000;

/// Bound on the node count. Both directions are linear in the node count, so a
/// stream that never closes the tree reaches this ceiling cheaply and errors.
const MAX_DECODE_NODES: i64 = 1 << 20;

fn walk_decode(
    m: usize,
    depth: usize,
    model: &mut AdaptiveModel,
    coder: &mut RangeDecoder,
    budget: &mut i64,
) -> Result<Term, String> {
    *budget -= 1;
    if *budget < 0 || depth > MAX_DECODE_DEPTH {
        return Err("malformed compressed data (does not terminate)".to_string());
    }
    let weights = model.kind_weights(m);
    let kind = decode_symbol(coder, &weights)?;
    model.saw_kind(m, kind);
    if kind == KIND_VAR {
        return Ok(Term::Var(decode_index(coder, model, m)? as u32));
    }
    if kind == KIND_LAM {
        return Ok(Term::Lam(Box::new(walk_decode(m + 1, depth + 1, model, coder, budget)?)));
    }
    let fun = walk_decode(m, depth + 1, model, coder, budget)?;
    let arg = walk_decode(m, depth + 1, model, coder, budget)?;
    Ok(Term::App(Box::new(fun), Box::new(arg)))
}

/// Closed lambda term -> compact bytes: the range coder's byte stream, whose
/// four-byte flush tail makes it self-delimiting against the structural end of
/// the walk. Errors on a non-closed term, like encode.
pub fn compress(term: &Term) -> Result<Vec<u8>, String> {
    check_closed(term, 0)?;
    let mut coder = RangeEncoder::new();
    let mut model = AdaptiveModel::new();
    walk_encode(term, 0, &mut model, &mut coder);
    // The decoder reads zero bytes past the stream, so trailing zero bytes are
    // redundant; drop them (an all-zero stream becomes empty).
    let mut out = coder.finish();
    while out.last() == Some(&0) {
        out.pop();
    }
    Ok(out)
}

/// Any byte string terminates: it yields a term or returns Err. decompress's
/// domain is compress outputs; without an integrity check it accepts any
/// bytes, decoding each to some term or signalling an error.
pub fn decompress(data: &[u8]) -> Result<Term, String> {
    let mut coder = RangeDecoder::new(data);
    let mut model = AdaptiveModel::new();
    let mut budget = MAX_DECODE_NODES;
    walk_decode(0, 0, &mut model, &mut coder, &mut budget)
}

// ---------------------------------------------------------------- self-test

fn bits_for_index(number: u64) -> String {
    BigNat::from_u64(number + 1).to_binary_string()[1..].to_string()
}

/// Parse one term of Tromp's binary lambda calculus, for brute-force checks.
fn blc_parse(bits: &[u8], pos: usize) -> Option<(Term, usize)> {
    if pos >= bits.len() {
        return None;
    }
    if bits[pos] == b'0' {
        if pos + 1 >= bits.len() {
            return None;
        }
        let tag = bits[pos + 1];
        let (first, after_first) = blc_parse(bits, pos + 2)?;
        if tag == b'0' {
            return Some((Term::Lam(Box::new(first)), after_first));
        }
        let (second, after_second) = blc_parse(bits, after_first)?;
        return Some((
            Term::App(Box::new(first), Box::new(second)),
            after_second,
        ));
    }
    let mut end = pos;
    while end < bits.len() && bits[end] == b'1' {
        end += 1;
    }
    if end >= bits.len() {
        return None;
    }
    Some((Term::Var((end - pos) as u32), end + 1))
}

fn max_free(term: &Term, depth: i64) -> i64 {
    match term {
        Term::Var(i) => *i as i64 - depth,
        Term::Lam(body) => max_free(body, depth + 1),
        Term::App(fun, arg) => max_free(fun, depth).max(max_free(arg, depth)),
    }
}

fn brute_force_terms(n: usize) -> Vec<Term> {
    let mut found = Vec::new();
    for value in 0..(1u64 << n) {
        let bits: Vec<u8> = (0..n)
            .map(|i| {
                if (value >> (n - 1 - i)) & 1 != 0 {
                    b'1'
                } else {
                    b'0'
                }
            })
            .collect();
        if let Some((term, consumed)) = blc_parse(&bits, 0) {
            if consumed == n && max_free(&term, 0) <= 0 {
                found.push(term);
            }
        }
    }
    found
}

fn brute_force_count(n: usize, cap: Option<usize>) -> u64 {
    brute_force_terms(n)
        .iter()
        .filter(|term| {
            cap.is_none_or(|cap| max_de_bruijn_index(term) as usize <= cap)
        })
        .count() as u64
}

fn self_test() {
    for cap in [None, Some(1), Some(2), Some(5)] {
        let mut table = Table::new(cap);
        table.extend(14);
        for n in 4..=14 {
            assert_eq!(
                *table.count_built(n, 0),
                BigNat::from_u64(brute_force_count(n, cap)),
                "count vs brute force at size {}",
                n
            );
        }
        for number in 0..3000u64 {
            let bits = bits_for_index(number);
            let term = decode(&mut table, &bits).unwrap();
            assert_eq!(encode(&mut table, &term).unwrap(), bits, "round trip");
        }
    }
    let mut unbounded = Table::new(None);
    let mut capped = Table::new(Some(8));
    unbounded.extend(9);
    let agree = unbounded.cumulative_built(9).clone();
    let mut number = 0u64;
    while BigNat::from_u64(number) < agree {
        let bits = bits_for_index(number);
        assert_eq!(
            decode(&mut unbounded, &bits).unwrap(),
            decode(&mut capped, &bits).unwrap(),
            "cap agreement on small sizes"
        );
        number += 1;
    }
    // set_index_cap in every direction agrees with a freshly built table
    for (from, to) in [(2usize, Some(7usize)), (7, Some(2)), (3, None)] {
        let mut observed = Table::new(Some(from));
        observed.extend(30);
        observed.set_index_cap(to);
        let mut fresh = Table::new(to);
        fresh.extend(30);
        for n in 4..=30 {
            assert_eq!(
                observed.count_built(n, 0),
                fresh.count_built(n, 0),
                "cap change {} at size {}",
                from,
                n
            );
        }
    }
    for cap in [Some(5), None] {
        let mut saved = Table::new(cap);
        saved.extend(40);
        let path = temp_table_path("selftest");
        let path = path.to_str().unwrap();
        saved.save_to_file(path).unwrap();
        let mut loaded = Table::load_from_file(path).unwrap();
        for n in 4..=40 {
            assert_eq!(
                loaded.count_built(n, 0),
                saved.count_built(n, 0),
                "save/load at size {}",
                n
            );
        }
        let bits = bits_for_index(1235);
        let term = decode(&mut loaded, &bits).unwrap();
        assert_eq!(encode(&mut loaded, &term).unwrap(), bits, "round trip loaded");
        // corruption is rejected: flipped byte, truncation, bad magic
        let bytes = fs::read(path).unwrap();
        let mid = bytes.len() / 2;
        let mut flipped = bytes.clone();
        flipped[mid] ^= 1;
        fs::write(path, &flipped).unwrap();
        assert!(Table::load_from_file(path).is_err(), "flipped byte");
        fs::write(path, &bytes[..bytes.len() - 3]).unwrap();
        assert!(Table::load_from_file(path).is_err(), "truncated");
        let mut bad_magic = bytes.clone();
        bad_magic[0] ^= 1;
        fs::write(path, &bad_magic).unwrap();
        assert!(Table::load_from_file(path).is_err(), "bad magic");
        fs::remove_file(path).unwrap();
    }
    {
        let mut cases: Vec<Term> = (4..=12).flat_map(brute_force_terms).collect();
        // Highly compressible terms whose node count far exceeds their byte
        // count, exercising the node ceiling that bounds decompression.
        cases.push(lam_chain(48));
        cases.push(lam_chain(1000));
        cases.push(church(2000));
        cases.extend(compression_vector_terms().into_iter().map(|(_, t)| t));
        for term in &cases {
            assert_eq!(
                &decompress(&compress(term).unwrap()).unwrap(),
                term,
                "compression round trip"
            );
        }
    }
    {
        // error paths: encode/decode/compress reject bad input
        let mut table = Table::new(None);
        table.extend(14);
        for bad in ["2", "01x"] {
            assert!(decode(&mut table, bad).is_err(), "decode rejects bad bits");
        }
        let bad_terms = [
            Term::Var(1),
            Term::App(Box::new(Term::Var(1)), Box::new(Term::Var(1))),
            Term::Lam(Box::new(Term::Var(2))),
            Term::Lam(Box::new(Term::Var(0))),
        ];
        for bad in &bad_terms {
            assert!(encode(&mut table, bad).is_err(), "encode rejects non-closed");
            assert!(compress(bad).is_err(), "compress rejects non-closed");
        }
        // decompress stays well-behaved (Ok or Err, never a panic) on garbage
        let mut garbage = vec![0u8, 0, 4, 0];
        garbage.extend([0xFFu8; 128]);
        for blob in [vec![], vec![0u8], vec![0u8, 0, 0], garbage] {
            let _ = decompress(&blob);
        }
    }
    println!("self-test passed");
}

/// A unique temp path in the OS temp directory (portable across Windows,
/// Linux and macOS), so concurrent runs do not collide on one name.
fn temp_table_path(tag: &str) -> std::path::PathBuf {
    std::env::temp_dir()
        .join(format!("lambda_b2l_{}_{}.tmp", tag, std::process::id()))
}

// Peak resident memory in KiB, read from Linux /proc; returns -1 on other
// platforms (this is a benchmark-only detail, not used by the library).
// Each --bench block runs as its own process so VmHWM is that block's alone.
fn peak_rss_kb() -> i64 {
    fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|status| {
            status.lines().find_map(|line| {
                line.strip_prefix("VmHWM:")?
                    .split_whitespace()
                    .next()?
                    .parse::<i64>()
                    .ok()
            })
        })
        .unwrap_or(-1)
}

fn sample_bits(j: usize, length: usize) -> String {  // length >= 5 (caller checks)
    let mut bits: String = (0..5)
        .rev()
        .map(|b| if (j >> b) & 1 != 0 { '1' } else { '0' })
        .collect();
    bits.push_str(&"0".repeat(length - 5));
    bits
}

fn stats_of(xs: &[f64]) -> (f64, f64, f64, f64) {
    let mean = xs.iter().sum::<f64>() / xs.len() as f64;
    let variance = xs.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>()
        / (xs.len() - 1).max(1) as f64;
    let lo = xs.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    (mean, variance.sqrt(), lo, hi)
}

fn bench_block(cap_arg: &str, length: usize) {
    let cap = if cap_arg == "inf" {
        None
    } else {
        Some(cap_arg.parse::<usize>().unwrap())
    };
    let mut table = Table::new(cap);
    let top = BigNat::from_binary_string(&format!("1{}", "1".repeat(length)))
        .sub(&BigNat::from_u64(1));
    let t0 = std::time::Instant::now();
    let mut n_max = 4;
    loop {
        table.extend(n_max);
        if *table.cumulative_built(n_max) > top {
            break;
        }
        n_max += 1;
    }
    let build_seconds = t0.elapsed().as_secs_f64();

    let mut decode_times = Vec::new();
    let mut encode_times = Vec::new();
    let mut sampled = Vec::new();
    for j in 0..64 {
        let bits = sample_bits(j, length);
        let s0 = std::time::Instant::now();
        let term = decode(&mut table, &bits).unwrap();
        decode_times.push(s0.elapsed().as_secs_f64());
        sampled.push((term, bits));
    }
    for (term, bits) in &sampled {
        let s0 = std::time::Instant::now();
        assert_eq!(&encode(&mut table, term).unwrap(), bits, "bench round trip");
        encode_times.push(s0.elapsed().as_secs_f64());
    }

    let tmp = temp_table_path("bench");
    let tmp = tmp.to_str().unwrap();
    table.save_to_file(tmp).unwrap();
    let disk_bytes = fs::metadata(tmp).unwrap().len();
    fs::remove_file(tmp).unwrap();

    let entries: usize = table.rows.iter().map(Vec::len).sum();
    let table_bits: usize = table
        .rows
        .iter()
        .flat_map(|row| row.iter().map(BigNat::bit_length))
        .sum();
    let dec = stats_of(&decode_times);
    let enc = stats_of(&encode_times);
    println!(
        "rust,{},{},{},{:.6},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{},{},{},{}",
        cap_arg, length, n_max, build_seconds,
        dec.0 * 1e6, dec.1 * 1e6, dec.2 * 1e6, dec.3 * 1e6,
        enc.0 * 1e6, enc.1 * 1e6, enc.2 * 1e6, enc.3 * 1e6,
        entries, table_bits, disk_bytes, peak_rss_kb()
    );
}

fn church(n: usize) -> Term {
    let mut body = Term::Var(1);
    for _ in 0..n {
        body = Term::App(Box::new(Term::Var(2)), Box::new(body));
    }
    Term::Lam(Box::new(Term::Lam(Box::new(body))))
}

/// A nest of n lambdas over Var(1): n+1 nodes whose compressed size is a few
/// bytes, so its node count far exceeds its byte count.
fn lam_chain(n: usize) -> Term {
    let mut body = Term::Var(1);
    for _ in 0..n {
        body = Term::Lam(Box::new(body));
    }
    body
}

/// Fixed term set shared with the Python and C++ implementations; the
/// compressed bytes must match across languages exactly.
fn compression_vector_terms() -> Vec<(String, Term)> {
    let s_comb = Term::Lam(Box::new(Term::Lam(Box::new(Term::Lam(Box::new(
        Term::App(
            Box::new(Term::App(Box::new(Term::Var(3)), Box::new(Term::Var(1)))),
            Box::new(Term::App(Box::new(Term::Var(2)), Box::new(Term::Var(1)))),
        ),
    ))))));
    let y_half = Term::Lam(Box::new(Term::App(
        Box::new(Term::Var(2)),
        Box::new(Term::App(Box::new(Term::Var(1)), Box::new(Term::Var(1)))),
    )));
    let y_comb = Term::Lam(Box::new(Term::App(
        Box::new(y_half.clone()),
        Box::new(y_half),
    )));
    let mut repetitive = s_comb.clone();
    for _ in 0..5 {
        repetitive = Term::App(Box::new(repetitive.clone()), Box::new(repetitive));
    }
    let mut bits = vec![b'0'; 192];
    for b in 0..32 {
        if (987654322u64 >> b) & 1 != 0 {
            bits[191 - b] = b'1';
        }
    }
    let mut table = Table::new(None);
    let uniform = decode(&mut table, &String::from_utf8(bits).unwrap()).unwrap();
    vec![
        ("S".to_string(), s_comb),
        ("Y".to_string(), y_comb),
        ("church10".to_string(), church(10)),
        ("church100".to_string(), church(100)),
        ("rep32S".to_string(), repetitive),
        ("uniform192".to_string(), uniform),
    ]
}

fn print_compress_vectors() {
    for (name, term) in compression_vector_terms() {
        let hex: String = compress(&term)
            .unwrap()
            .iter()
            .map(|byte| format!("{:02x}", byte))
            .collect();
        println!("{}\t{}", name, hex);
    }
}

fn print_vectors() {
    for (cap, limit) in [(None, 500u64), (Some(3), 300u64)] {
        let mut table = Table::new(cap);
        println!(
            "# cap={}",
            cap.map_or("inf".to_string(), |c| c.to_string())
        );
        for number in 0..limit {
            let bits = bits_for_index(number);
            let term = decode(&mut table, &bits).unwrap();
            println!("{}\t{}\t{}", number, show_bits(&bits), show_term(&term, 0));
        }
    }
}

/// Structural FNV-1a 64-bit digest of a term, computed iteratively in
/// pre-order (function subterm before argument). The Python, C++ and Rust
/// fuzz drivers compute it identically, so the same input bytes must yield the
/// same digest in every language; any divergence is a cross-implementation bug.
fn fuzz_digest(root: &Term) -> (u64, u64) {
    let prime: u64 = 0x100000001b3;
    let mut h: u64 = 0xcbf29ce484222325;
    let mut nodes: u64 = 0;
    let mut stack: Vec<&Term> = vec![root];
    while let Some(t) = stack.pop() {
        nodes += 1;
        match t {
            Term::Var(i) => {
                h = (h ^ 0x56).wrapping_mul(prime);
                for b in 0..4 {
                    h = (h ^ (((*i >> (8 * b)) & 0xFF) as u64)).wrapping_mul(prime);
                }
            }
            Term::Lam(body) => {
                h = (h ^ 0x4C).wrapping_mul(prime);
                stack.push(body);
            }
            Term::App(fun, arg) => {
                h = (h ^ 0x41).wrapping_mul(prime);
                stack.push(arg); // argument pushed first
                stack.push(fun); // function popped first
            }
        }
    }
    (h, nodes)
}

/// Decompress every hex-blob line of a file; emit one result line each, either
/// "OK\t<digest>\t<nodes>" or "ERR" for any malformed input. Drives the
/// cross-language differential fuzz against the Python reference.
fn run_fuzz(path: &str) -> Result<(), String> {
    let data = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let mut out = String::new();
    for raw in data.lines() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        let bytes = if line.len() % 2 == 0 {
            (0..line.len())
                .step_by(2)
                .map(|i| u8::from_str_radix(&line[i..i + 2], 16))
                .collect::<Result<Vec<u8>, _>>()
                .ok()
        } else {
            None
        };
        match bytes {
            None => out.push_str("ERR\n"),
            Some(b) => match decompress(&b) {
                Ok(term) => {
                    let (h, nodes) = fuzz_digest(&term);
                    out.push_str(&format!("OK\t{:016x}\t{}\n", h, nodes));
                }
                Err(_) => out.push_str("ERR\n"),
            },
        }
    }
    print!("{}", out);
    Ok(())
}

fn run(args: &[String]) -> Result<(), String> {
    match args.get(1).map(String::as_str) {
        Some("--vectors") => print_vectors(),
        Some("--fuzz") => run_fuzz(args.get(2).ok_or("--fuzz needs a path")?)?,
        Some("--compress-vectors") => print_compress_vectors(),
        Some("--bench") => {
            let length: usize = args[3].parse().map_err(|_| "bad length")?;
            if length < 5 {
                return Err("--bench length must be >= 5".to_string());
            }
            bench_block(&args[2], length);
        }
        Some("--save-table") => {
            let path = &args[2];
            let mut table = Table::new(Some(5));
            table.extend(40);
            table.save_to_file(path).map_err(|e| e.to_string())?;
            println!("saved cap-5 size-40 table to {}", path);
        }
        Some("--load-table") => {
            let loaded = Table::load_from_file(&args[2])?;
            let mut fresh = Table::new(loaded.cap);
            fresh.extend(loaded.built_size());
            for n in 2..=loaded.built_size() {
                assert_eq!(
                    loaded.count_built(n, 0),
                    fresh.count_built(n, 0),
                    "loaded table at size {}",
                    n
                );
            }
            println!(
                "table file OK (cap {}, size {})",
                loaded.cap.map_or("inf".to_string(), |c| c.to_string()),
                loaded.built_size()
            );
        }
        _ => {
            self_test();
            let mut table = Table::new(None);
            println!("first strings of the canonical bijection:");
            for number in 0..8u64 {
                let bits = bits_for_index(number);
                let term = decode(&mut table, &bits).unwrap();
                println!("  {}  ->  {}", show_bits(&bits), show_term(&term, 0));
            }
        }
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if let Err(error) = run(&args) {  // clean message, not a panic/abort
        eprintln!("error: {}", error);
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    /// `cargo test` entry point: the embedded self-test validates counts
    /// against exhaustive enumeration, round trips both directions across
    /// caps, table save/load with corruption detection, and the compression
    /// round trip, panicking via assert on any failure.
    #[test]
    fn self_test_passes() {
        super::self_test();
    }

    /// A cap below 1 admits no terms, so it is rejected at construction.
    #[test]
    #[should_panic(expected = "index cap must be at least 1")]
    fn rejects_cap_below_one() {
        let _ = super::Table::new(Some(0));
    }
}
