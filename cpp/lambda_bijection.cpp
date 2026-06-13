// Binary2Lambda — lambda_bijection.cpp: self-contained single-file C++
// implementation.
//
// Bijection between binary strings and closed untyped lambda terms.
// Every binary string (including the empty string) denotes exactly one
// closed lambda term and vice versa; encode and decode share one
// incrementally-built counting table.
//
// Specification (identical to the Python and Wolfram Language versions):
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
// of size <= K+1. A finite cap shrinks the table from Theta(n^2) to
// Theta(K n) entries.
//
// Tables can be saved to and loaded from disk in a portable text format
// shared with the Python implementation (see Table::saveToFile).
//
// Arbitrary-precision arithmetic is provided by the self-contained BigNat
// class below, keeping this file dependency-free.  It needs only C++17 and a
// 64-bit integer type, so it builds on GCC, Clang and MSVC across Windows,
// Linux and macOS (a 64x64->128 multiply uses __int128 / _umul128 where
// available and a portable 32-bit split otherwise).  For heavy workloads the
// same interface maps onto GMP (mpz_class): BigNat's operations are a strict
// subset of mpz semantics, so swapping in a thin mpz_class alias is
// mechanical and removes the main constant-factor cost.
//
// Non-ASCII output (the λ and ε glyphs) is written as explicit UTF-8 bytes,
// so the program emits UTF-8 regardless of the compiler's source charset.
//
// The file also contains the project's compression layer (namespace
// lambda_compress): a renormalizing range coder over the term grammar with an
// adaptive model, byte-compatible with python/lambda_compress.py. Compression
// is a separate axis from the bijection - it needs compress + decompress, not
// "every bit string is a term".
//
// Build:  g++ -std=c++17 -O2 -o lambda_bijection_cpp lambda_bijection.cpp
//         (clang++ likewise; MSVC: cl /std:c++17 /EHsc /O2 lambda_bijection.cpp)
// Run:    ./lambda_bijection_cpp                     (self-test and demo)
//         ./lambda_bijection_cpp --vectors           (cross-language vectors)
//         ./lambda_bijection_cpp --compress-vectors  (compression vectors)
//         ./lambda_bijection_cpp --save-table <path> (write a sample table)
//         ./lambda_bijection_cpp --load-table <path> (load and verify a table)
//         ./lambda_bijection_cpp --bench <cap> <L>   (one benchmark block, CSV)

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <random>
#include <string>
#include <utility>
#include <vector>
#if defined(_MSC_VER)
#include <intrin.h>  // _umul128 on MSVC (no __int128 there)
#endif

namespace lambda_bijection {

// ----------------------------------------------------------------- numbers

// Arbitrary-precision unsigned integer, little-endian base-2^64 limbs.
// Only the operations the bijection needs: + , - (no underflow), * ,
// divmod, comparison, and binary-string conversion.  The arithmetic uses
// only 64-bit limbs plus a 64x64->128 multiply helper, so it builds on any
// C++17 compiler (GCC, Clang and MSVC); no 128-bit integer type is required.
class BigNat {
 public:
  BigNat() = default;
  explicit BigNat(std::uint64_t value) {
    if (value != 0) limbs_.push_back(value);
  }

  bool isZero() const { return limbs_.empty(); }

  int compare(const BigNat& other) const {
    if (limbs_.size() != other.limbs_.size())
      return limbs_.size() < other.limbs_.size() ? -1 : 1;
    for (std::size_t i = limbs_.size(); i-- > 0;)
      if (limbs_[i] != other.limbs_[i])
        return limbs_[i] < other.limbs_[i] ? -1 : 1;
    return 0;
  }
  bool operator==(const BigNat& o) const { return compare(o) == 0; }
  bool operator<(const BigNat& o) const { return compare(o) < 0; }
  bool operator<=(const BigNat& o) const { return compare(o) <= 0; }

  BigNat& operator+=(const BigNat& other) {
    limbs_.resize(std::max(limbs_.size(), other.limbs_.size()), 0);
    std::uint64_t carry = 0;  // always 0 or 1
    for (std::size_t i = 0; i < limbs_.size(); ++i) {
      const std::uint64_t addend =
          i < other.limbs_.size() ? other.limbs_[i] : 0;
      std::uint64_t sum = limbs_[i] + addend;
      std::uint64_t carryOut = sum < limbs_[i];
      sum += carry;
      carryOut += sum < carry;
      limbs_[i] = sum;
      carry = carryOut;
    }
    if (carry != 0) limbs_.push_back(carry);
    return *this;
  }
  BigNat operator+(const BigNat& other) const {
    BigNat result = *this;
    result += other;
    return result;
  }

  BigNat operator-(const BigNat& other) const {  // requires *this >= other
    assert(other <= *this && "BigNat subtraction would underflow");
    BigNat result = *this;
    std::uint64_t borrow = 0;  // always 0 or 1
    for (std::size_t i = 0; i < result.limbs_.size(); ++i) {
      const std::uint64_t sub = i < other.limbs_.size() ? other.limbs_[i] : 0;
      const std::uint64_t diff = result.limbs_[i] - sub;
      std::uint64_t borrowOut = result.limbs_[i] < sub;
      const std::uint64_t afterBorrow = diff - borrow;
      borrowOut += diff < borrow;
      result.limbs_[i] = afterBorrow;
      borrow = borrowOut;
    }
    result.trim();
    return result;
  }

  BigNat operator*(const BigNat& other) const {
    if (isZero() || other.isZero()) return BigNat();
    std::vector<std::uint64_t> out(limbs_.size() + other.limbs_.size(), 0);
    for (std::size_t i = 0; i < limbs_.size(); ++i) {
      std::uint64_t carry = 0;
      for (std::size_t j = 0; j < other.limbs_.size(); ++j) {
        std::uint64_t high, low;
        mul64(limbs_[i], other.limbs_[j], high, low);
        low += out[i + j];
        high += low < out[i + j];
        low += carry;
        high += low < carry;
        out[i + j] = low;
        carry = high;  // the full product never overflows the high limb
      }
      out[i + other.limbs_.size()] = carry;
    }
    BigNat result;
    result.limbs_ = std::move(out);
    result.trim();
    return result;
  }

  static std::pair<BigNat, BigNat> divmod(const BigNat& num,
                                          const BigNat& den) {
    if (den.isZero()) throw std::logic_error("division by zero");
    BigNat quotient, remainder;
    for (std::size_t i = num.bitLength(); i-- > 0;) {
      remainder.shiftLeftOneBit();
      if (num.bit(i)) remainder.setBit(0);
      if (den <= remainder) {
        remainder = remainder - den;
        quotient.setBit(i);
      }
    }
    return {quotient, remainder};
  }

  std::size_t bitLength() const {
    if (limbs_.empty()) return 0;
    return limbs_.size() * 64 - __builtin_clzll(limbs_.back());
  }
  std::uint64_t toUint64() const {  // requires bitLength() <= 64
    assert(bitLength() <= 64 && "BigNat::toUint64 truncates a value >= 2^64");
    return limbs_.empty() ? 0 : limbs_[0];
  }
  bool bit(std::size_t i) const {
    std::size_t limb = i / 64;
    return limb < limbs_.size() && ((limbs_[limb] >> (i % 64)) & 1u) != 0;
  }
  void setBit(std::size_t i) {
    std::size_t limb = i / 64;
    if (limb >= limbs_.size()) limbs_.resize(limb + 1, 0);
    limbs_[limb] |= std::uint64_t(1) << (i % 64);
  }
  void shiftLeftOneBit() {
    std::uint64_t carry = 0;
    for (auto& limb : limbs_) {
      std::uint64_t next = limb >> 63;
      limb = (limb << 1) | carry;
      carry = next;
    }
    if (carry != 0) limbs_.push_back(carry);
  }

  std::string toBinaryString() const {  // most significant bit first
    if (isZero()) return "0";
    std::string bits;
    for (std::size_t i = bitLength(); i-- > 0;) bits.push_back(bit(i) ? '1' : '0');
    return bits;
  }
  static BigNat fromBinaryString(const std::string& bits) {
    BigNat result;
    for (char c : bits) {
      result.shiftLeftOneBit();
      if (c == '1') result.setBit(0);
    }
    return result;
  }

  std::string toHexString() const {  // lowercase, no leading zeros
    if (isZero()) return "0";
    std::ostringstream out;
    out << std::hex << limbs_.back();
    for (std::size_t i = limbs_.size() - 1; i-- > 0;) {
      out.width(16);
      out.fill('0');
      out << limbs_[i];
    }
    return out.str();
  }
  static BigNat fromHexString(const std::string& hex) {
    BigNat result;
    for (char c : hex) {
      int digit;
      if (c >= '0' && c <= '9') digit = c - '0';
      else if (c >= 'a' && c <= 'f') digit = c - 'a' + 10;
      else if (c >= 'A' && c <= 'F') digit = c - 'A' + 10;
      else throw std::invalid_argument("invalid hex digit in table file");
      for (int b = 0; b < 4; ++b) result.shiftLeftOneBit();
      if (digit != 0) {
        if (result.limbs_.empty()) result.limbs_.push_back(0);
        result.limbs_[0] |= static_cast<std::uint64_t>(digit);
      }
    }
    return result;
  }

 private:
  // 64x64 -> 128 multiply, as (high, low).  Fast paths for GCC/Clang
  // (__int128) and MSVC x64 (_umul128); a portable 32-bit split otherwise.
  static void mul64(std::uint64_t a, std::uint64_t b, std::uint64_t& high,
                    std::uint64_t& low) {
#if defined(__SIZEOF_INT128__)
    const unsigned __int128 product = static_cast<unsigned __int128>(a) * b;
    low = static_cast<std::uint64_t>(product);
    high = static_cast<std::uint64_t>(product >> 64);
#elif defined(_MSC_VER) && defined(_M_X64)
    low = _umul128(a, b, &high);
#else
    const std::uint64_t aLo = a & 0xFFFFFFFFu, aHi = a >> 32;
    const std::uint64_t bLo = b & 0xFFFFFFFFu, bHi = b >> 32;
    const std::uint64_t ll = aLo * bLo, lh = aLo * bHi;
    const std::uint64_t hl = aHi * bLo, hh = aHi * bHi;
    const std::uint64_t mid = (ll >> 32) + (lh & 0xFFFFFFFFu) + (hl & 0xFFFFFFFFu);
    low = (ll & 0xFFFFFFFFu) | (mid << 32);
    high = hh + (lh >> 32) + (hl >> 32) + (mid >> 32);
#endif
  }

  void trim() {
    while (!limbs_.empty() && limbs_.back() == 0) limbs_.pop_back();
  }
  std::vector<std::uint64_t> limbs_;
};

// ------------------------------------------------------------------- terms

struct Term;
using TermPtr = std::shared_ptr<const Term>;

struct Term {
  enum class Kind { Var, Lam, App };
  Kind kind;
  int index = 0;        // Var: de Bruijn index, >= 1
  TermPtr left, right;  // Lam: left = body; App: left = fun, right = arg
};

inline TermPtr var(int index) {
  return std::make_shared<const Term>(Term{Term::Kind::Var, index, {}, {}});
}
inline TermPtr lam(TermPtr body) {
  return std::make_shared<const Term>(
      Term{Term::Kind::Lam, 0, std::move(body), {}});
}
inline TermPtr app(TermPtr fun, TermPtr arg) {
  return std::make_shared<const Term>(
      Term{Term::Kind::App, 0, std::move(fun), std::move(arg)});
}

inline int termSize(const TermPtr& term) {
  switch (term->kind) {
    case Term::Kind::Var: return term->index + 1;
    case Term::Kind::Lam: return termSize(term->left) + 2;
    default: return termSize(term->left) + termSize(term->right) + 2;
  }
}

inline int maxDeBruijnIndex(const TermPtr& term) {
  switch (term->kind) {
    case Term::Kind::Var: return term->index;
    case Term::Kind::Lam: return maxDeBruijnIndex(term->left);
    default:
      return std::max(maxDeBruijnIndex(term->left),
                      maxDeBruijnIndex(term->right));
  }
}

inline bool sameTerm(const TermPtr& a, const TermPtr& b) {
  if (a->kind != b->kind) return false;
  switch (a->kind) {
    case Term::Kind::Var: return a->index == b->index;
    case Term::Kind::Lam: return sameTerm(a->left, b->left);
    default:
      return sameTerm(a->left, b->left) && sameTerm(a->right, b->right);
  }
}

// Throw unless term is a well-formed closed term: every de Bruijn index is
// >= 1 and bound by an enclosing lambda (index <= depth).  encode/compress
// require this: the bijection ranges over closed terms, so only a closed term
// has a code.
inline void checkClosed(const TermPtr& term, int depth = 0) {
  switch (term->kind) {
    case Term::Kind::Var:
      if (term->index < 1)
        throw std::invalid_argument("de Bruijn index must be >= 1");
      if (term->index > depth)
        throw std::invalid_argument("term is not closed (free variable)");
      return;
    case Term::Kind::Lam:
      checkClosed(term->left, depth + 1);
      return;
    default:
      checkClosed(term->left, depth);
      checkClosed(term->right, depth);
  }
}

// FNV-1a 64-bit, used as a table-file integrity checksum; identical integer
// arithmetic in the Python and Rust implementations, so a saved table's
// checksum is the same in every language.  Folding incrementally lets save
// and load stream the file without holding the whole body in memory.
constexpr std::uint64_t kFnvOffset = 0xCBF29CE484222325ULL;

inline std::uint64_t fnv1a64Update(std::uint64_t h, const char* data,
                                   std::size_t len) {
  for (std::size_t i = 0; i < len; ++i) {
    h ^= static_cast<unsigned char>(data[i]);
    h *= 0x100000001B3ULL;
  }
  return h;
}
inline std::uint64_t fnv1a64Update(std::uint64_t h, const std::string& data) {
  return fnv1a64Update(h, data.data(), data.size());
}

// Magic for the portable binary table format (.lamtab); the last byte is the
// format version. Shared byte-for-byte with the Python, Rust and Wolfram ports.
const std::string kTableMagic("LAMTAB\x01", 7);

// Little-endian u32 and big-endian magnitude helpers for that format.
inline void putU32(std::string& out, std::uint32_t v) {
  out.push_back(static_cast<char>(v & 0xFF));
  out.push_back(static_cast<char>((v >> 8) & 0xFF));
  out.push_back(static_cast<char>((v >> 16) & 0xFF));
  out.push_back(static_cast<char>((v >> 24) & 0xFF));
}
inline std::uint32_t getU32(const std::string& d, std::size_t pos) {
  return static_cast<std::uint8_t>(d[pos]) |
         (static_cast<std::uint32_t>(static_cast<std::uint8_t>(d[pos + 1])) << 8) |
         (static_cast<std::uint32_t>(static_cast<std::uint8_t>(d[pos + 2])) << 16) |
         (static_cast<std::uint32_t>(static_cast<std::uint8_t>(d[pos + 3])) << 24);
}
inline std::string toBytesBigEndian(const BigNat& v) {
  const std::size_t nbytes = (v.bitLength() + 7) / 8;
  std::string out(nbytes, '\0');
  for (std::size_t k = 0; k < nbytes; ++k) {  // k counts bytes from the LSB
    unsigned char byte = 0;
    for (int j = 0; j < 8; ++j)
      if (v.bit(8 * k + j)) byte |= static_cast<unsigned char>(1u << j);
    out[nbytes - 1 - k] = static_cast<char>(byte);
  }
  return out;
}
inline BigNat fromBytesBigEndian(const std::string& d, std::size_t pos,
                                 std::size_t n) {
  BigNat v;
  for (std::size_t i = 0; i < n; ++i) {
    for (int b = 0; b < 8; ++b) v.shiftLeftOneBit();
    const unsigned char byte = static_cast<unsigned char>(d[pos + i]);
    for (int b = 7; b >= 0; --b)
      if ((byte >> b) & 1u) v.setBit(static_cast<std::size_t>(b));
  }
  return v;
}

// context: 0 = top level / lambda body, 1 = left of app, 2 = right of app
inline std::string showTerm(const TermPtr& term, int context = 0) {
  switch (term->kind) {
    case Term::Kind::Var:
      return std::to_string(term->index);
    case Term::Kind::Lam: {
      // "\xCE\xBB" is the UTF-8 encoding of λ, written as explicit bytes so
      // it is independent of the compiler's source-charset handling.
      std::string rendered = "\xCE\xBB" + showTerm(term->left, 0);
      return context > 0 ? "(" + rendered + ")" : rendered;
    }
    default: {
      std::string rendered =
          showTerm(term->left, 1) + " " + showTerm(term->right, 2);
      return context == 2 ? "(" + rendered + ")" : rendered;
    }
  }
}

inline std::string showBits(const std::string& bits) {
  return bits.empty() ? "\xCE\xB5" : bits;  // "\xCE\xB5" is UTF-8 for ε
}

// ------------------------------------------------------------------- table

// Counting table T(n, m): number of terms of size n whose free indices are
// all <= m and whose indices nowhere exceed the cap. Rows are stored under
// the effective context min(m, n-1, cap). Size extension is append-only;
// changing the cap reuses every row of size <= min(old, new) + 1 (those are
// cap-independent) and rebuilds only the rest.
class Table {
 public:
  explicit Table(std::optional<int> indexCap = std::nullopt)
      : cap_(indexCap), rows_(2), cum_(2) {
    if (cap_ && *cap_ < 1)
      throw std::invalid_argument("index cap must be at least 1");
  }

  int builtSize() const { return static_cast<int>(rows_.size()) - 1; }
  std::optional<int> indexCap() const { return cap_; }

  long long entryCount() const {
    long long total = 0;
    for (const auto& row : rows_) total += static_cast<long long>(row.size());
    return total;
  }
  long long totalBitLength() const {
    long long total = 0;
    for (const auto& row : rows_)
      for (const auto& value : row)
        total += static_cast<long long>(value.bitLength());
    return total;
  }

  void extend(int sizeLimit) {
    for (int n = static_cast<int>(rows_.size()); n <= sizeLimit; ++n) {
      std::vector<BigNat> row;
      for (int m = 0; m < width(n); ++m) row.push_back(rowValue(n, m));
      BigNat cumulative = cum_.back();
      if (!row.empty()) cumulative += row[0];
      rows_.push_back(std::move(row));
      cum_.push_back(std::move(cumulative));
    }
  }

  void setIndexCap(std::optional<int> newCap) {
    if (newCap && *newCap < 1)
      throw std::invalid_argument("index cap must be at least 1");
    if (newCap == cap_) return;
    const int lower = cap_ && newCap ? std::min(*cap_, *newCap)
                                     : (cap_ ? *cap_ : *newCap);
    const int keep = lower + 1;  // rows of size <= keep are cap-independent
    const int target = builtSize();
    cap_ = newCap;
    if (builtSize() > keep) {
      rows_.resize(keep + 1);
      cum_.resize(keep + 1);
    }
    extend(target);
  }

  // count()/closedCumulative() return a reference into rows_/cum_ and may
  // call extend(), which only appends: references into existing inner rows
  // survive (their buffers are moved, not copied), but a returned cum_
  // reference is invalidated by the next extend(), so consume it before any
  // further call.  rankOf/unrank never trigger extension because encode/decode
  // pre-extend to the term's size first; that is what keeps these safe.
  const BigNat& count(int n, int m) {
    if (n < 2) return zero();
    if (n > builtSize()) extend(n);
    return rows_[n][effectiveContext(n, m)];
  }

  const BigNat& closedCumulative(int n) {  // closed terms of size <= n
    if (n < 2) return zero();
    if (n > builtSize()) extend(n);
    return cum_[n];
  }

  // Portable binary format shared with the Python, Rust and Wolfram ports
  // (identical bytes in every language): the 7-byte magic LAMTAB\x01, a
  // one-byte cap kind (0 unbounded, 1 finite) and 4-byte little-endian cap
  // value, the 4-byte little-endian built size, then for each size from 2
  // upward its row of counts, each a 4-byte little-endian length and that many
  // big-endian magnitude bytes, then an 8-byte little-endian FNV-1a-64 of every
  // preceding byte. Cumulative counts are derivable, so they are not stored.
  // The body is streamed, so memory stays flat.
  void saveToFile(const std::string& path) const {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write " + path);
    std::uint64_t checksum = kFnvOffset;
    auto emit = [&](const std::string& chunk) {
      out.write(chunk.data(), static_cast<std::streamsize>(chunk.size()));
      checksum = fnv1a64Update(checksum, chunk);
    };
    std::string header = kTableMagic;
    header.push_back(cap_ ? char(1) : char(0));
    putU32(header, cap_ ? static_cast<std::uint32_t>(*cap_) : 0u);
    putU32(header, static_cast<std::uint32_t>(builtSize()));
    emit(header);
    for (int n = 2; n <= builtSize(); ++n) {
      for (const BigNat& v : rows_[n]) {
        const std::string mag = toBytesBigEndian(v);
        std::string cell;
        putU32(cell, static_cast<std::uint32_t>(mag.size()));
        cell += mag;
        emit(cell);
      }
    }
    std::string trailer;
    for (int i = 0; i < 8; ++i)
      trailer.push_back(static_cast<char>((checksum >> (8 * i)) & 0xFF));
    out.write(trailer.data(), static_cast<std::streamsize>(trailer.size()));
  }

  static Table loadFromFile(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read " + path);
    const std::string data((std::istreambuf_iterator<char>(in)),
                           std::istreambuf_iterator<char>());
    const std::size_t headerSize = kTableMagic.size() + 1 + 4 + 4;
    if (data.size() < headerSize + 8 ||
        data.compare(0, kTableMagic.size(), kTableMagic) != 0)
      throw std::runtime_error("not a lambda-binarization table file");
    const std::size_t end = data.size() - 8;
    const std::uint64_t checksum = fnv1a64Update(kFnvOffset, data.data(), end);
    std::uint64_t stored = 0;
    for (int i = 0; i < 8; ++i)
      stored |= static_cast<std::uint64_t>(static_cast<std::uint8_t>(data[end + i]))
                << (8 * i);
    if (checksum != stored)
      throw std::runtime_error("table file failed checksum (corrupt)");
    const std::uint8_t capKind = static_cast<std::uint8_t>(data[kTableMagic.size()]);
    const std::uint32_t capValue = getU32(data, kTableMagic.size() + 1);
    const std::uint32_t built = getU32(data, kTableMagic.size() + 5);
    if (capKind > 1 || built < 1)
      throw std::runtime_error("table file has a malformed header");
    Table table(capKind == 0 ? std::optional<int>()
                             : std::optional<int>(static_cast<int>(capValue)));
    std::size_t pos = headerSize;
    for (int n = 2; n <= static_cast<int>(built); ++n) {
      std::vector<BigNat> row;
      for (int j = 0; j < table.width(n); ++j) {
        if (pos + 4 > end) throw std::runtime_error("table file is truncated");
        const std::uint32_t nbytes = getU32(data, pos);
        pos += 4;
        if (pos + nbytes > end) throw std::runtime_error("table file is truncated");
        row.push_back(fromBytesBigEndian(data, pos, nbytes));
        pos += nbytes;
      }
      BigNat cumulative = table.cum_.back();
      if (!row.empty()) cumulative += row[0];
      table.rows_.push_back(std::move(row));
      table.cum_.push_back(std::move(cumulative));
    }
    if (pos != end)
      throw std::runtime_error("table file size does not match its rows");
    return table;
  }

 private:
  static const BigNat& zero() {
    static const BigNat kZero;
    return kZero;
  }

  int effectiveContext(int n, int m) const {
    m = std::min(m, n - 1);
    if (cap_) m = std::min(m, *cap_);
    return m;
  }
  int width(int n) const { return effectiveContext(n, n - 1) + 1; }

  // T(n, m) from already-built smaller rows; m is an effective context.
  BigNat rowValue(int n, int m) const {
    BigNat value(m == n - 1 ? 1 : 0);                 // the variable n-1
    value += builtCount(n - 2, m + 1);                // abstractions
    for (int k = 2; k <= n - 4; ++k)                  // applications
      value += builtCount(k, m) * builtCount(n - 2 - k, m);
    return value;
  }
  const BigNat& builtCount(int n, int m) const {
    if (n < 2) return zero();
    return rows_[n][effectiveContext(n, m)];
  }

  std::optional<int> cap_;
  std::vector<std::vector<BigNat>> rows_;
  std::vector<BigNat> cum_;
};

// ------------------------------------------------------------------- codec

inline int varCount(const Table& table, int n, int m) {
  if (n < 2 || n - 1 > m) return 0;
  if (table.indexCap() && n - 1 > *table.indexCap()) return 0;
  return 1;
}

inline BigNat rankOf(Table& table, const TermPtr& term, int m) {
  if (term->kind == Term::Kind::Var) return BigNat();
  const int n = termSize(term);
  BigNat rank(varCount(table, n, m));
  if (term->kind == Term::Kind::Lam) {
    rank += rankOf(table, term->left, m + 1);
    return rank;
  }
  rank += table.count(n - 2, m + 1);
  const int leftSize = termSize(term->left);
  for (int k = 2; k < leftSize; ++k)
    rank += table.count(k, m) * table.count(n - 2 - k, m);
  rank += rankOf(table, term->left, m) *
          table.count(termSize(term->right), m);
  rank += rankOf(table, term->right, m);
  return rank;
}

// Precondition: 0 <= rank < T(n, m).  The divmod split below is exact only
// because rank < count(k,m)*right, so the quotient stays below count(k,m).
inline TermPtr unrank(Table& table, BigNat rank, int n, int m) {
  if (varCount(table, n, m) == 1) {
    if (rank.isZero()) return var(n - 1);
    rank = rank - BigNat(1);
  }
  const BigNat lamBlock = table.count(n - 2, m + 1);
  if (rank < lamBlock) return lam(unrank(table, std::move(rank), n - 2, m + 1));
  rank = rank - lamBlock;
  for (int k = 2; k <= n - 4; ++k) {
    const BigNat right = table.count(n - 2 - k, m);
    const BigNat block = table.count(k, m) * right;
    if (rank < block) {
      auto [leftRank, rightRank] = BigNat::divmod(rank, right);
      return app(unrank(table, std::move(leftRank), k, m),
                 unrank(table, std::move(rightRank), n - 2 - k, m));
    }
    rank = rank - block;
  }
  throw std::logic_error("rank out of range for size class");
}

// Closed lambda term -> binary string (inverse of decode).  Throws if term
// is not closed or uses an index above the cap.
inline std::string encode(Table& table, const TermPtr& term) {
  checkClosed(term);
  if (table.indexCap() && maxDeBruijnIndex(term) > *table.indexCap())
    throw std::invalid_argument("term uses indices above the table cap");
  const int n = termSize(term);
  table.extend(n);
  BigNat number = table.closedCumulative(n - 1) + rankOf(table, term, 0);
  number += BigNat(1);
  return number.toBinaryString().substr(1);
}

// Binary string -> closed lambda term (total on {0,1}*).
inline TermPtr decode(Table& table, const std::string& bits) {
  for (char c : bits)
    if (c != '0' && c != '1')
      throw std::invalid_argument("input must consist of 0s and 1s");
  const BigNat number = BigNat::fromBinaryString("1" + bits) - BigNat(1);
  int n = 4;  // the smallest closed term, λ1, has size 4
  while (table.closedCumulative(n) <= number) ++n;
  return unrank(table, number - table.closedCumulative(n - 1), n, 0);
}

}  // namespace lambda_bijection

// ------------------------------------------------------------- compression
//
// Lambda-specific compression of closed terms - the project's separate
// compression axis. A renormalizing range coder over the term grammar with an
// adaptive model: the coder walks the term in pre-order, the choice set at
// each node is {Var(1..m), Lam, App} with m the lambda depth, and termination
// is structural, so no counting table and no size header are needed. Both
// directions run in time linear in the node count. Byte-compatible with
// python/lambda_compress.py.
namespace lambda_compress {

using lambda_bijection::app;
using lambda_bijection::lam;
using lambda_bijection::Term;
using lambda_bijection::TermPtr;
using lambda_bijection::var;

constexpr std::size_t kKindVar = 0, kKindLam = 1, kKindApp = 2;

// 32-bit Subbotin range coder. Registers are unsigned 32-bit and arithmetic
// wraps modulo 2^32, which makes the byte stream identical across the Python,
// C++ and Rust implementations. Every symbol's total frequency stays below
// kBot; the model guarantees that by rescaling, and indices above the bucket
// count are coded bit by bit.
constexpr std::uint32_t kTop = 1u << 24;
constexpr std::uint32_t kBot = 1u << 16;

// Narrows [low, low+range) per symbol, emitting settled high bytes.
class RangeEncoder {
 public:
  void encode(std::uint32_t cLow, std::uint32_t freq, std::uint32_t total) {
    range_ /= total;
    low_ += cLow * range_;
    range_ *= freq;
    while ((low_ ^ (low_ + range_)) < kTop ||
           (range_ < kBot && ((range_ = (0u - low_) & (kBot - 1)), true))) {
      out_.push_back(static_cast<std::uint8_t>(low_ >> 24));
      low_ <<= 8;
      range_ <<= 8;
    }
  }
  std::vector<std::uint8_t> finish() {
    for (int i = 0; i < 4; ++i) {
      out_.push_back(static_cast<std::uint8_t>(low_ >> 24));
      low_ <<= 8;
    }
    return std::move(out_);
  }

 private:
  std::uint32_t low_ = 0, range_ = 0xFFFFFFFFu;
  std::vector<std::uint8_t> out_;
};

// Mirrors the encoder; the code register steers symbol choice.
class RangeDecoder {
 public:
  explicit RangeDecoder(const std::vector<std::uint8_t>& data) : data_(data) {
    for (int i = 0; i < 4; ++i) code_ = (code_ << 8) | nextByte();
  }
  std::uint32_t target(std::uint32_t total) {
    range_ /= total;
    const std::uint32_t value = (code_ - low_) / range_;
    return value >= total ? total - 1 : value;
  }
  void consume(std::uint32_t cLow, std::uint32_t freq) {
    low_ += cLow * range_;
    range_ *= freq;
    while ((low_ ^ (low_ + range_)) < kTop ||
           (range_ < kBot && ((range_ = (0u - low_) & (kBot - 1)), true))) {
      code_ = (code_ << 8) | nextByte();
      low_ <<= 8;
      range_ <<= 8;
    }
  }

 private:
  std::uint8_t nextByte() {  // past the stream the decoder reads zero bytes
    return pos_ < data_.size() ? data_[pos_++] : 0;
  }
  const std::vector<std::uint8_t>& data_;
  std::size_t pos_ = 0;
  std::uint32_t low_ = 0, range_ = 0xFFFFFFFFu, code_ = 0;
};

// Counts shared by compress and decompress, updated identically. Kinds are
// conditioned on the lambda-depth bucket min(m, 3). Variable indices share an
// 8-entry bucketed table (bucket min(i, 8)); an index in the tail bucket is
// followed by its offset, coded bit by bit. Counts are halved when a context
// total reaches kRescaleLimit, keeping every total below the coder's kBot.
struct AdaptiveModel {
  std::uint32_t kinds[4][3];
  std::uint32_t indices[8];
  static constexpr std::uint32_t kIncrement = 16;
  static constexpr int kIndexBuckets = 8;
  static constexpr std::uint32_t kRescaleLimit = 1u << 14;

  AdaptiveModel() {
    for (auto& row : kinds)
      for (auto& weight : row) weight = 1;
    for (auto& weight : indices) weight = 1;
  }
  std::array<std::uint32_t, 3> kindWeights(int m) const {
    const auto& row = kinds[std::min(m, 3)];
    return {m >= 1 ? row[0] : 0, row[1], row[2]};
  }
  void sawKind(int m, std::size_t kind) {
    auto& row = kinds[std::min(m, 3)];
    row[kind] += kIncrement;
    if (row[0] + row[1] + row[2] >= kRescaleLimit)
      for (auto& weight : row) weight = std::max<std::uint32_t>(weight >> 1, 1);
  }
  void sawIndex(int index) {
    indices[std::min(index - 1, kIndexBuckets - 1)] += kIncrement;
    std::uint32_t sum = 0;
    for (auto weight : indices) sum += weight;
    if (sum >= kRescaleLimit)
      for (auto& weight : indices)
        weight = std::max<std::uint32_t>(weight >> 1, 1);
  }
};

inline void encodeSymbol(RangeEncoder& coder, const std::uint32_t* weights,
                         std::size_t n, std::size_t symbol) {
  std::uint32_t cLow = 0, total = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (i < symbol) cLow += weights[i];
    total += weights[i];
  }
  coder.encode(cLow, weights[symbol], total);
}

inline std::size_t decodeSymbol(RangeDecoder& coder,
                                const std::uint32_t* weights, std::size_t n) {
  std::uint32_t total = 0;
  for (std::size_t i = 0; i < n; ++i) total += weights[i];
  const std::uint32_t target = coder.target(total);
  std::uint32_t cLow = 0;
  for (std::size_t symbol = 0; symbol < n; ++symbol) {
    if (target < cLow + weights[symbol]) {
      coder.consume(cLow, weights[symbol]);
      return symbol;
    }
    cLow += weights[symbol];
  }
  throw std::runtime_error("malformed compressed data (code point out of range)");
}

inline void encodeBit(RangeEncoder& coder, std::uint32_t bit) {
  coder.encode(bit, 1, 2);
}

inline std::uint32_t decodeBit(RangeDecoder& coder) {
  const std::uint32_t bit = coder.target(2) >= 1 ? 1 : 0;
  coder.consume(bit, 1);
  return bit;
}

// Elias-gamma over an equiprobable bit model for the rare de Bruijn index
// above the bucket count; the unary length is capped so malformed input
// cannot spin in the prefix.
constexpr int kGammaMaxBits = 40;

inline void encodeGamma(RangeEncoder& coder, std::uint32_t value) {
  const std::uint32_t v = value + 1;
  int n = 0;
  while ((v >> (n + 1)) != 0) ++n;  // bit length of v, minus one
  for (int i = 0; i < n; ++i) encodeBit(coder, 0);
  encodeBit(coder, 1);
  for (int i = n - 1; i >= 0; --i) encodeBit(coder, (v >> i) & 1u);
}

inline std::uint32_t decodeGamma(RangeDecoder& coder) {
  int n = 0;
  while (decodeBit(coder) == 0)
    if (++n > kGammaMaxBits)
      throw std::runtime_error("malformed compressed data (index code too long)");
  std::uint32_t v = 1;
  for (int i = 0; i < n; ++i) v = (v << 1) | decodeBit(coder);
  return v - 1;
}

inline void encodeIndex(RangeEncoder& coder, AdaptiveModel& model, int m,
                        int index) {
  const int alphabet = std::min(m, AdaptiveModel::kIndexBuckets);
  const int bucket = std::min(index - 1, AdaptiveModel::kIndexBuckets - 1);
  encodeSymbol(coder, model.indices, static_cast<std::size_t>(alphabet),
               static_cast<std::size_t>(bucket));
  if (bucket == AdaptiveModel::kIndexBuckets - 1)
    encodeGamma(coder, static_cast<std::uint32_t>(
                           index - AdaptiveModel::kIndexBuckets));
  model.sawIndex(index);
}

inline int decodeIndex(RangeDecoder& coder, AdaptiveModel& model, int m) {
  const int alphabet = std::min(m, AdaptiveModel::kIndexBuckets);
  const int bucket = static_cast<int>(
      decodeSymbol(coder, model.indices, static_cast<std::size_t>(alphabet)));
  int index;
  if (bucket < AdaptiveModel::kIndexBuckets - 1) {
    index = bucket + 1;
  } else {
    index = AdaptiveModel::kIndexBuckets + static_cast<int>(decodeGamma(coder));
    if (index > m)
      throw std::runtime_error("malformed compressed data (index exceeds depth)");
  }
  model.sawIndex(index);
  return index;
}

inline void walkEncode(const TermPtr& term, int m, AdaptiveModel& model,
                       RangeEncoder& coder) {
  const auto weights = model.kindWeights(m);
  switch (term->kind) {
    case Term::Kind::Var:
      encodeSymbol(coder, weights.data(), 3, kKindVar);
      model.sawKind(m, kKindVar);
      encodeIndex(coder, model, m, term->index);
      break;
    case Term::Kind::Lam:
      encodeSymbol(coder, weights.data(), 3, kKindLam);
      model.sawKind(m, kKindLam);
      walkEncode(term->left, m + 1, model, coder);
      break;
    default:
      encodeSymbol(coder, weights.data(), 3, kKindApp);
      model.sawKind(m, kKindApp);
      walkEncode(term->left, m, model, coder);
      walkEncode(term->right, m, model, coder);
  }
}

// Bound on the term's nesting depth (App-spines deepen without raising the
// lambda context m), kept well under every implementation's stack: a stream
// that keeps nesting hits this limit and errors while the stack is still safe.
// Real lambda terms stay far below it.
constexpr int kMaxDecodeDepth = 12000;

// Bound on the node count.  Both directions are linear in the node count, so a
// stream that never closes the tree reaches this ceiling cheaply and errors.
constexpr long long kMaxDecodeNodes = 1LL << 20;

inline TermPtr walkDecode(int m, int depth, long long& budget,
                          AdaptiveModel& model, RangeDecoder& coder) {
  if (--budget < 0 || depth > kMaxDecodeDepth)
    throw std::runtime_error("malformed compressed data (does not terminate)");
  const auto weights = model.kindWeights(m);
  const std::size_t kind = decodeSymbol(coder, weights.data(), 3);
  model.sawKind(m, kind);
  if (kind == kKindVar)
    return var(decodeIndex(coder, model, m));
  if (kind == kKindLam)
    return lam(walkDecode(m + 1, depth + 1, budget, model, coder));
  TermPtr fun = walkDecode(m, depth + 1, budget, model, coder);
  return app(std::move(fun), walkDecode(m, depth + 1, budget, model, coder));
}

// Closed lambda term -> compact bytes: the range coder's byte stream, whose
// four-byte flush tail makes it self-delimiting against the structural end of
// the walk.
inline std::vector<std::uint8_t> compress(const TermPtr& term) {
  lambda_bijection::checkClosed(term);
  RangeEncoder coder;
  AdaptiveModel model;
  walkEncode(term, 0, model, coder);
  // The decoder reads zero bytes past the stream, so trailing zero bytes are
  // redundant; drop them (an all-zero stream becomes empty).
  std::vector<std::uint8_t> out = coder.finish();
  while (!out.empty() && out.back() == 0) out.pop_back();
  return out;
}

// Any byte string terminates: it yields a term or throws std::runtime_error.
// decompress's domain is compress outputs; without an integrity check it
// accepts any bytes, decoding each to some term or signalling an error.
inline TermPtr decompress(const std::vector<std::uint8_t>& data) {
  RangeDecoder coder(data);
  AdaptiveModel model;
  long long budget = kMaxDecodeNodes;
  return walkDecode(0, 0, budget, model, coder);
}

}  // namespace lambda_compress

// --------------------------------------------------------------- self-test

namespace {

using namespace lambda_bijection;

std::string bitsForIndex(std::uint64_t number) {
  BigNat value(number + 1);
  return value.toBinaryString().substr(1);
}

// Parse one term of Tromp's binary lambda calculus, for brute-force checks.
std::optional<std::pair<TermPtr, std::size_t>> blcParse(
    const std::string& bits, std::size_t pos) {
  if (pos >= bits.size()) return std::nullopt;
  if (bits[pos] == '0') {
    if (pos + 1 >= bits.size()) return std::nullopt;
    const char tag = bits[pos + 1];
    auto first = blcParse(bits, pos + 2);
    if (!first) return std::nullopt;
    if (tag == '0') return {{lam(first->first), first->second}};
    auto second = blcParse(bits, first->second);
    if (!second) return std::nullopt;
    return {{app(first->first, second->first), second->second}};
  }
  std::size_t end = pos;
  while (end < bits.size() && bits[end] == '1') ++end;
  if (end >= bits.size()) return std::nullopt;
  return {{var(static_cast<int>(end - pos)), end + 1}};
}

int maxFree(const TermPtr& term, int depth = 0) {
  switch (term->kind) {
    case Term::Kind::Var: return term->index - depth;
    case Term::Kind::Lam: return maxFree(term->left, depth + 1);
    default:
      return std::max(maxFree(term->left, depth), maxFree(term->right, depth));
  }
}

std::vector<TermPtr> bruteForceTerms(int n) {
  std::vector<TermPtr> found;
  for (std::uint64_t value = 0; value < (std::uint64_t(1) << n); ++value) {
    std::string bits(n, '0');
    for (int i = 0; i < n; ++i)
      if ((value >> (n - 1 - i)) & 1u) bits[i] = '1';
    auto parsed = blcParse(bits, 0);
    if (!parsed || parsed->second != static_cast<std::size_t>(n)) continue;
    if (maxFree(parsed->first) > 0) continue;
    found.push_back(parsed->first);
  }
  return found;
}

std::uint64_t bruteForceCount(int n, std::optional<int> cap) {
  std::uint64_t found = 0;
  for (const TermPtr& term : bruteForceTerms(n))
    if (!cap || maxDeBruijnIndex(term) <= *cap) ++found;
  return found;
}

TermPtr church(int n) {
  TermPtr body = var(1);
  for (int i = 0; i < n; ++i) body = app(var(2), std::move(body));
  return lam(lam(std::move(body)));
}

// A nest of n lambdas over var(1): n+1 nodes whose compressed size is a few
// bytes, so its node count far exceeds its byte count.
TermPtr lamChain(int n) {
  TermPtr body = var(1);
  for (int i = 0; i < n; ++i) body = lam(std::move(body));
  return body;
}

// Fixed term set shared with the Python and Rust implementations; the
// compressed bytes must match across languages exactly.
std::vector<std::pair<std::string, TermPtr>> compressionVectorTerms() {
  TermPtr sComb = lam(lam(lam(app(app(var(3), var(1)), app(var(2), var(1))))));
  TermPtr yHalf = lam(app(var(2), app(var(1), var(1))));
  TermPtr repetitive = sComb;
  for (int i = 0; i < 5; ++i) repetitive = app(repetitive, repetitive);
  std::string bits(192, '0');
  for (int b = 0; b < 32; ++b)
    if ((987654322u >> b) & 1u) bits[191 - b] = '1';
  Table table;
  TermPtr uniform = decode(table, bits);
  return {{"S", sComb}, {"Y", lam(app(yHalf, yHalf))},
          {"church10", church(10)}, {"church100", church(100)},
          {"rep32S", repetitive}, {"uniform192", uniform}};
}

void require(bool condition, const std::string& what) {
  if (!condition) throw std::logic_error("self-test failed: " + what);
}

template <typename Fn>
bool throwsException(Fn&& fn) {
  try {
    fn();
    return false;
  } catch (const std::exception&) {
    return true;
  }
}

// A unique temp path in the OS temp directory (portable across Windows,
// Linux and macOS), so concurrent runs do not collide on one name.
std::string tempTablePath(const std::string& tag) {
  std::random_device rng;
  std::ostringstream name;
  name << "lambda_b2l_" << tag << "_" << std::hex << rng() << rng() << ".tmp";
  return (std::filesystem::temp_directory_path() / name.str()).string();
}

void selfTest() {
  const std::optional<int> caps[] = {std::nullopt, 1, 2, 5};
  for (auto cap : caps) {
    const std::string capStr = cap ? std::to_string(*cap) : "inf";
    Table table(cap);
    for (int n = 4; n <= 14; ++n)
      require(table.count(n, 0) == BigNat(bruteForceCount(n, cap)),
              "count vs brute force cap=" + capStr + " n=" + std::to_string(n));
    for (std::uint64_t number = 0; number < 3000; ++number) {
      const std::string bits = bitsForIndex(number);
      require(encode(table, decode(table, bits)) == bits,
              "round trip cap=" + capStr + " number=" + std::to_string(number));
    }
  }
  Table unbounded;
  Table capped(8);
  const BigNat agree = unbounded.closedCumulative(9);
  for (std::uint64_t number = 0; BigNat(number) < agree; ++number) {
    const std::string bits = bitsForIndex(number);
    require(sameTerm(decode(unbounded, bits), decode(capped, bits)),
            "cap agreement number=" + std::to_string(number));
  }
  // set_index_cap in every direction agrees with a freshly built table
  for (const auto& [from, to] : std::vector<std::pair<int, std::optional<int>>>{
           {2, 7}, {7, 2}, {3, std::nullopt}}) {
    Table observed(from);
    observed.extend(30);
    observed.setIndexCap(to);
    Table fresh(to);
    fresh.extend(30);
    for (int n = 4; n <= 30; ++n)
      require(observed.count(n, 0) == fresh.count(n, 0),
              "cap change " + std::to_string(from) + " n=" + std::to_string(n));
  }
  for (const std::optional<int> cap : {std::optional<int>(5), std::optional<int>()}) {
    Table saved(cap);
    saved.extend(40);
    const std::string path = tempTablePath("selftest");
    saved.saveToFile(path);
    Table loaded = Table::loadFromFile(path);
    for (int n = 4; n <= 40; ++n)
      require(loaded.count(n, 0) == saved.count(n, 0),
              "save/load n=" + std::to_string(n));
    const std::string bits = bitsForIndex(1235);
    require(encode(loaded, decode(loaded, bits)) == bits,
            "round trip on loaded table");
    // corruption is rejected: flipped byte, truncation, bad magic
    std::string bytes;
    {
      std::ifstream in(path, std::ios::binary);
      bytes.assign((std::istreambuf_iterator<char>(in)),
                   std::istreambuf_iterator<char>());
    }
    auto writeBytes = [&](const std::string& b) {
      std::ofstream out(path, std::ios::binary);
      out.write(b.data(), static_cast<std::streamsize>(b.size()));
    };
    std::string flipped = bytes;
    flipped[flipped.size() / 2] ^= 1;
    writeBytes(flipped);
    require(throwsException([&] { Table::loadFromFile(path); }),
            "load rejects a flipped byte");
    writeBytes(bytes.substr(0, bytes.size() - 3));
    require(throwsException([&] { Table::loadFromFile(path); }),
            "load rejects truncation");
    std::string badMagic = bytes;
    badMagic[0] ^= 1;
    writeBytes(badMagic);
    require(throwsException([&] { Table::loadFromFile(path); }),
            "load rejects bad magic");
    std::remove(path.c_str());
  }
  {
    std::vector<TermPtr> cases;
    for (int n = 4; n <= 12; ++n)
      for (const TermPtr& term : bruteForceTerms(n)) cases.push_back(term);
    // Highly compressible terms whose node count far exceeds their byte count,
    // exercising the node ceiling that bounds decompression.
    cases.push_back(lamChain(48));
    cases.push_back(lamChain(1000));
    cases.push_back(church(2000));
    for (const auto& [name, term] : compressionVectorTerms())
      cases.push_back(term);
    for (const TermPtr& term : cases)
      require(sameTerm(lambda_compress::decompress(
                           lambda_compress::compress(term)),
                       term),
              "compression round trip");
  }
  {  // error paths: encode/decode/compress reject bad input
    Table table;
    table.extend(14);
    for (const std::string& badBits : {std::string("2"), std::string("01x")})
      require(throwsException([&] { decode(table, badBits); }),
              "decode rejects bad bits");
    const std::vector<TermPtr> badTerms = {var(1), app(var(1), var(1)),
                                           lam(var(2)), lam(var(0))};
    for (const TermPtr& bad : badTerms) {
      require(throwsException([&] { encode(table, bad); }),
              "encode rejects a non-closed term");
      require(throwsException([&] { lambda_compress::compress(bad); }),
              "compress rejects a non-closed term");
    }
    // decompress stays well-behaved (returns a term or throws) on garbage
    std::vector<std::vector<std::uint8_t>> blobs = {
        {}, {0x00}, {0x00, 0x00, 0x00}, {0x00, 0x00, 0x04, 0x00}};
    blobs.back().insert(blobs.back().end(), 128, 0xFF);
    for (const auto& blob : blobs)
      try {
        lambda_compress::decompress(blob);
      } catch (const std::exception&) {
      }
  }
  std::cout << "self-test passed\n";
}

// Peak resident memory in KiB, read from Linux /proc; returns -1 on other
// platforms (this is a benchmark-only detail, not used by the library).
// Each --bench block runs as its own process so VmHWM is that block's alone.
long long peakRssKb() {
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line))
    if (line.rfind("VmHWM:", 0) == 0) {
      std::istringstream fields(line.substr(6));
      long long kb = -1;
      fields >> kb;
      return kb;
    }
  return -1;
}

std::string sampleBits(int j, int length) {  // length >= 5 (checked by caller)
  std::string bits;
  for (int b = 4; b >= 0; --b) bits.push_back(((j >> b) & 1) != 0 ? '1' : '0');
  return bits + std::string(length - 5, '0');
}

struct Stats {
  double mean, stdev, min, max;
};

Stats statsOf(const std::vector<double>& xs) {
  double sum = 0, lo = xs[0], hi = xs[0];
  for (double x : xs) {
    sum += x;
    lo = std::min(lo, x);
    hi = std::max(hi, x);
  }
  const double mean = sum / xs.size();
  double variance = 0;
  for (double x : xs) variance += (x - mean) * (x - mean);
  return {mean,
          xs.size() > 1 ? std::sqrt(variance / (xs.size() - 1)) : 0.0,
          lo, hi};
}

void benchBlock(const std::string& capArg, int length) {
  using Clock = std::chrono::steady_clock;
  const std::optional<int> cap =
      capArg == "inf" ? std::optional<int>()
                      : std::optional<int>(std::stoi(capArg));
  Table table(cap);
  const BigNat top =
      BigNat::fromBinaryString("1" + std::string(length, '1')) - BigNat(1);
  const auto t0 = Clock::now();
  int nMax = 4;
  while (table.closedCumulative(nMax) <= top) ++nMax;
  const double buildSeconds =
      std::chrono::duration<double>(Clock::now() - t0).count();

  std::vector<double> decodeTimes, encodeTimes;
  std::vector<std::pair<TermPtr, std::string>> sampled;
  for (int j = 0; j < 64; ++j) {
    const std::string bits = sampleBits(j, length);
    const auto s0 = Clock::now();
    TermPtr term = decode(table, bits);
    decodeTimes.push_back(
        std::chrono::duration<double>(Clock::now() - s0).count());
    sampled.emplace_back(std::move(term), bits);
  }
  for (const auto& [term, bits] : sampled) {
    const auto s0 = Clock::now();
    require(encode(table, term) == bits, "bench round trip");
    encodeTimes.push_back(
        std::chrono::duration<double>(Clock::now() - s0).count());
  }

  const std::string tmp = tempTablePath("bench");
  table.saveToFile(tmp);
  const long long diskBytes =
      static_cast<long long>(std::filesystem::file_size(tmp));
  std::filesystem::remove(tmp);

  const Stats dec = statsOf(decodeTimes), enc = statsOf(encodeTimes);
  std::printf(
      "cpp,%s,%d,%d,%.6f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,"
      "%lld,%lld,%lld,%lld\n",
      capArg.c_str(), length, nMax, buildSeconds, dec.mean * 1e6,
      dec.stdev * 1e6, dec.min * 1e6, dec.max * 1e6, enc.mean * 1e6,
      enc.stdev * 1e6, enc.min * 1e6, enc.max * 1e6, table.entryCount(),
      table.totalBitLength(), diskBytes, peakRssKb());
}

void printVectors() {
  struct Mode { std::optional<int> cap; int limit; };
  for (const Mode& mode : {Mode{std::nullopt, 500}, Mode{3, 300}}) {
    Table table(mode.cap);
    std::cout << "# cap=" << (mode.cap ? std::to_string(*mode.cap) : "inf")
              << "\n";
    for (int number = 0; number < mode.limit; ++number) {
      const std::string bits = bitsForIndex(number);
      std::cout << number << "\t" << showBits(bits) << "\t"
                << showTerm(decode(table, bits)) << "\n";
    }
  }
}

// Structural FNV-1a 64-bit digest of a term, computed iteratively in
// pre-order (function subterm before argument). The Python, C++ and Rust
// fuzz drivers compute it identically, so the same input bytes must yield the
// same digest in every language; any divergence is a cross-implementation bug.
std::pair<std::uint64_t, std::uint64_t> fuzzDigest(const TermPtr& root) {
  std::uint64_t h = 0xcbf29ce484222325ULL, nodes = 0;
  const std::uint64_t prime = 0x100000001b3ULL;
  auto mix = [&](std::uint64_t byte) { h = (h ^ byte) * prime; };
  std::vector<const Term*> stack{root.get()};
  while (!stack.empty()) {
    const Term* t = stack.back();
    stack.pop_back();
    ++nodes;
    switch (t->kind) {
      case Term::Kind::Var: {
        mix(0x56);
        const std::uint32_t i = static_cast<std::uint32_t>(t->index);
        for (int b = 0; b < 4; ++b) mix((i >> (8 * b)) & 0xFFu);
        break;
      }
      case Term::Kind::Lam:
        mix(0x4C);
        stack.push_back(t->left.get());
        break;
      case Term::Kind::App:
        mix(0x41);
        stack.push_back(t->right.get());  // argument pushed first
        stack.push_back(t->left.get());   // function popped first
        break;
    }
  }
  return {h, nodes};
}

// Decompress every hex-blob line of a file; emit one result line each, either
// "OK\t<digest>\t<nodes>" or "ERR" for any malformed input. Drives the
// cross-language differential fuzz against the Python reference.
void runFuzz(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  auto hexVal = [](char c) -> int {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
  };
  std::string line, out;
  while (std::getline(in, line)) {
    while (!line.empty() &&
           (line.back() == '\r' || line.back() == ' ' || line.back() == '\t'))
      line.pop_back();
    if (line.empty()) continue;
    std::vector<std::uint8_t> bytes;
    bool ok = line.size() % 2 == 0;
    for (std::size_t i = 0; ok && i < line.size(); i += 2) {
      const int hi = hexVal(line[i]), lo = hexVal(line[i + 1]);
      if (hi < 0 || lo < 0) { ok = false; break; }
      bytes.push_back(static_cast<std::uint8_t>(hi * 16 + lo));
    }
    if (!ok) { out += "ERR\n"; continue; }
    try {
      const auto digest = fuzzDigest(lambda_compress::decompress(bytes));
      char buf[64];
      std::snprintf(buf, sizeof buf, "OK\t%016llx\t%llu\n",
                    static_cast<unsigned long long>(digest.first),
                    static_cast<unsigned long long>(digest.second));
      out += buf;
    } catch (const std::exception&) {
      out += "ERR\n";
    }
  }
  std::cout << out;
}

int run(int argc, char** argv) {
  const std::string mode = argc > 1 ? argv[1] : "";
  if (mode == "--vectors") {
    printVectors();
    return 0;
  }
  if (mode == "--fuzz" && argc > 2) {
    runFuzz(argv[2]);
    return 0;
  }
  if (mode == "--compress-vectors") {
    for (const auto& [name, term] : compressionVectorTerms()) {
      std::cout << name << "\t";
      for (std::uint8_t byte : lambda_compress::compress(term))
        std::printf("%02x", byte);
      std::cout << "\n";
    }
    return 0;
  }
  if (mode == "--bench" && argc > 3) {
    const int length = std::stoi(argv[3]);
    if (length < 5) throw std::invalid_argument("--bench length must be >= 5");
    benchBlock(argv[2], length);
    return 0;
  }
  if (mode == "--save-table" && argc > 2) {
    Table table(5);
    table.extend(40);
    table.saveToFile(argv[2]);
    std::cout << "saved cap-5 size-40 table to " << argv[2] << "\n";
    return 0;
  }
  if (mode == "--load-table" && argc > 2) {
    Table loaded = Table::loadFromFile(argv[2]);
    Table fresh(loaded.indexCap());
    fresh.extend(loaded.builtSize());
    for (int n = 2; n <= loaded.builtSize(); ++n)
      require(loaded.count(n, 0) == fresh.count(n, 0), "loaded table");
    std::cout << "table file OK (cap "
              << (loaded.indexCap() ? std::to_string(*loaded.indexCap())
                                    : "inf")
              << ", size " << loaded.builtSize() << ")\n";
    return 0;
  }
  selfTest();
  Table table;
  std::cout << "first strings of the canonical bijection:\n";
  for (std::uint64_t number = 0; number < 8; ++number) {
    const std::string bits = bitsForIndex(number);
    std::cout << "  " << showBits(bits) << "  ->  "
              << showTerm(decode(table, bits)) << "\n";
  }
  return 0;
}

}  // namespace

// Define LAMBDA_BINARIZATION_NO_MAIN before including this file to use it as a
// library: the bijection API stays in namespace lambda_bijection and the
// compression API in namespace lambda_compress, with no entry point of its own.
#ifndef LAMBDA_BINARIZATION_NO_MAIN
int main(int argc, char** argv) {
  try {
    return run(argc, argv);
  } catch (const std::exception& error) {  // clean message, not abort/core dump
    std::cerr << "error: " << error.what() << "\n";
    return 1;
  }
}
#endif
