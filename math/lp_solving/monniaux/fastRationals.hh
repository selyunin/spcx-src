#ifndef _FASTRATIONALS_H
#define _FASTRATIONALS_H

#include <iostream>
#include <gmpxx.h>
#include <limits.h>
#include <stdlib.h>
#include <errno.h>
#include <assert.h>

typedef int32_t word;
typedef uint32_t uword;
typedef int64_t lword;
typedef uint64_t ulword;

#ifdef INLINE_FAST_ARITHMETIC
#define FAST_RATIONAL_INLINE inline
#else
#define FAST_RATIONAL_INLINE
#endif

static inline size_t djb2(size_t a, size_t b) {
  return (a << 5) + a + b;
}

class FastInteger {
  bool has_mpz, has_word;
  word num;
  mpz_t mpz;

  void kill_mpz() {
    if (has_mpz) {
      mpz_clear(mpz);
      has_mpz = false;
    }
  }

  void make_mpz() {
    if (!has_mpz) {
      assert(has_word);
      mpz_init(mpz);
      has_mpz=true;
      mpz_set_si(mpz, num);
    }
  }

  void force_make_mpz() const {
    const_cast<FastInteger*>(this)->make_mpz();
  }

  void make_erase_mpz() {
    if (!has_mpz) {
      mpz_init(mpz);
      has_mpz=true;
    }
  }

  void make_word() {
    if (mpz_fits_slong_p(mpz)) {
      num = mpz_get_si(mpz);
      has_word = true;
    } else {
      has_word = false;
    }
  }

  FAST_RATIONAL_INLINE void init_read_string(const char *s);

public:
  FastInteger(const char *s) {
    init_read_string(s);
  }

  FastInteger(const std::string& s) {
    init_read_string(s.c_str());
  }

  FAST_RATIONAL_INLINE FastInteger(const mpz_class& x);
  FAST_RATIONAL_INLINE FastInteger(const mpz_t x);
  FAST_RATIONAL_INLINE FastInteger(const FastInteger& x);
  FAST_RATIONAL_INLINE FastInteger& operator=(const FastInteger& x);

  FastInteger() : has_mpz(false), has_word(true), num(0) {
  }

  ~FastInteger() {
    kill_mpz();
  }

  FastInteger(word x) : has_mpz(false), has_word(true), num(x) {
  }

  friend void FAST_RATIONAL_INLINE addition(FastInteger& dst, const FastInteger& a, const FastInteger& b);
  friend void FAST_RATIONAL_INLINE subtraction(FastInteger& dst, const FastInteger& a, const FastInteger& b);
  friend void FAST_RATIONAL_INLINE multiplication(FastInteger& dst, const FastInteger& a, const FastInteger& b);
  friend void FAST_RATIONAL_INLINE division(FastInteger& dst, const FastInteger& a, const FastInteger& b);
  friend void FAST_RATIONAL_INLINE divexact(FastInteger& dst, const FastInteger& a, const FastInteger& b);
  friend void FAST_RATIONAL_INLINE modulo(FastInteger& dst, const FastInteger& a, const FastInteger& b);
  friend void lcm(FastInteger& dst, const FastInteger& a, const FastInteger& b);  friend void gcd(FastInteger& dst, const FastInteger& a, const FastInteger& b);

  friend void FAST_RATIONAL_INLINE additionAssign(FastInteger& a, const FastInteger& b);
  friend void FAST_RATIONAL_INLINE subtractionAssign(FastInteger& a, const FastInteger& b);
  friend void FAST_RATIONAL_INLINE multiplicationAssign(FastInteger& a, const FastInteger& b);
  friend void FAST_RATIONAL_INLINE divisionAssign(FastInteger& a, const FastInteger& b);
  friend void FAST_RATIONAL_INLINE divexactAssign(FastInteger& a, const FastInteger& b);
  friend void FAST_RATIONAL_INLINE moduloAssign(FastInteger& a, const FastInteger& b);

  void print_details(std::ostream& out) const;
  void print(std::ostream& out) const;

  double get_d() const;

  FAST_RATIONAL_INLINE bool operator==(const FastInteger& b) const;

  FAST_RATIONAL_INLINE FastInteger operator-() const;

private:
  static inline int compare(lword a, lword b) {
    if (a < b) return -1;
    else if (a > b) return 1;
    else return 0;
  }

public:
  FAST_RATIONAL_INLINE int compare(const FastInteger& b) const;
  FAST_RATIONAL_INLINE int sign() const;

  bool operator<(const FastInteger& b) const {
    return compare(b) < 0;
  }

  bool operator>(const FastInteger& b) const {
    return compare(b) > 0;
  }

  bool operator<=(const FastInteger& b) const {
    return compare(b) <= 0;
  }

  bool operator>=(const FastInteger& b) const {
    return compare(b) >= 0;
  }

  bool operator!=(const FastInteger& b) const {
    return !(*this == b);
  }

  FAST_RATIONAL_INLINE unsigned size() const;

  friend class FastRational;

  FastInteger operator+(const FastInteger& b) const {
    FastInteger dest; 
    addition(dest, *this, b);
    return dest;
  }

  FastInteger operator-(const FastInteger& b) const {
    FastInteger dest;
    subtraction(dest, *this, b);
    return dest;
  }

  FastInteger operator*(const FastInteger& b) const {
    FastInteger dest;
    multiplication(dest, *this, b);
    return dest;
  }

  FastInteger operator/(const FastInteger& b) const {
    FastInteger dest;
    division(dest, *this, b);
    return dest;
  }

  FastInteger operator%(const FastInteger& b) const {
    FastInteger dest;
    modulo(dest, *this, b);
    return dest;
  }

  FastInteger& operator+=(const FastInteger& b) {
    additionAssign(*this, b);
    return *this;
  }

  FastInteger& operator-=(const FastInteger& b) {
    subtractionAssign(*this, b);
    return *this;
  }

  FastInteger& operator*=(const FastInteger& b) {
    multiplicationAssign(*this, b);
    return *this;
  }

  FastInteger& operator%=(const FastInteger& b) {
    moduloAssign(*this, b);
    return *this;
  }

  inline FastInteger& operator/=(const FastInteger& b) {
    divisionAssign(*this, b);
    return *this;
  }

  size_t hash() const {
    if (has_word) {
      return num;
    } else {
      return mpz_get_si(mpz);
    }
  }

  bool isZero() const {
    if (has_word) {
      return num==0;
    } else {
      return mpz_sgn(mpz)==0;
    }
  }
};

inline std::ostream& operator<<(std::ostream& out, const FastInteger& r) {
  r.print(out);
  return out;
}

class FastRational {
  bool has_mpq, has_word;
  word num;
  uword den;
  mpq_t mpq;


  void kill_mpq() {
    if (has_mpq) {
      mpq_clear(mpq);
      has_mpq = false;
    }
  }

  void make_mpq() {
    if (!has_mpq) {
      assert(has_word);
      mpq_init(mpq);
      has_mpq=true;
      mpz_set_si(mpq_numref(mpq), num);
      mpz_set_ui(mpq_denref(mpq), den);
    }
  }

  void force_make_mpq() const {
    const_cast<FastRational*>(this)->make_mpq();
  }

  void make_erase_mpq() {
    if (!has_mpq) {
      mpq_init(mpq);
      has_mpq=true;
    }
  }

  void make_word() {
    assert(has_mpq);
    if (mpz_fits_slong_p(mpq_numref(mpq))
	&& mpz_fits_ulong_p(mpq_denref(mpq))) {
      num = mpz_get_si(mpq_numref(mpq));
      den = mpz_get_ui(mpq_denref(mpq));
      has_word = true;
    } else {
      has_word = false;
    }
  }

public:
#if 0
  FastRational(const mpz_class& x) {
    if (mpz_fits_slong_p(x.get_mpz_t())) {
      num = mpz_get_si(x.get_mpz_t());
      den = 1;
      has_word = true;
      has_mpq = false;
    } else {
      has_word = false;
      mpq_init(mpq);
      has_mpq = true;
      mpz_set(mpq_numref(mpq), x.get_mpz_t());
    }
  }
#endif

  FAST_RATIONAL_INLINE FastRational(const FastInteger& x);
  FAST_RATIONAL_INLINE FastRational(const FastRational& x);
  FAST_RATIONAL_INLINE FastRational& operator=(const FastRational& x);

  FastRational() : has_mpq(false), has_word(true), num(0), den(1) {
  }

  ~FastRational() {
    kill_mpq();
  }

  FastRational(word x) : has_mpq(false), has_word(true), num(x), den(1) {
  }

  friend FAST_RATIONAL_INLINE void addition(FastRational& dst, const FastRational& a, const FastRational& b);
  friend FAST_RATIONAL_INLINE void subtraction(FastRational& dst, const FastRational& a, const FastRational& b);
  friend FAST_RATIONAL_INLINE void multiplication(FastRational& dst, const FastRational& a, const FastRational& b);
  friend FAST_RATIONAL_INLINE void division(FastRational& dst, const FastRational& a, const FastRational& b);

  friend FAST_RATIONAL_INLINE void additionAssign(FastRational& a, const FastRational& b);
  friend FAST_RATIONAL_INLINE void subtractionAssign(FastRational& a, const FastRational& b);
  friend FAST_RATIONAL_INLINE void multiplicationAssign(FastRational& a, const FastRational& b);
  friend FAST_RATIONAL_INLINE void divisionAssign(FastRational& a, const FastRational& b);

  void print_details(std::ostream& out) const;
  void print(std::ostream& out) const;

  FAST_RATIONAL_INLINE FastRational(word n, uword d);

  FAST_RATIONAL_INLINE double get_d() const;

  FAST_RATIONAL_INLINE bool operator==(const FastRational& b) const;

  FAST_RATIONAL_INLINE FastRational operator-() const;

private:
  static inline int compare(lword a, lword b) {
    if (a < b) return -1;
    else if (a > b) return 1;
    else return 0;
  }

public:
  FAST_RATIONAL_INLINE int compare(const FastRational& b) const;
  FAST_RATIONAL_INLINE int sign() const;

  bool operator<(const FastRational& b) const {
    return compare(b) < 0;
  }

  bool operator>(const FastRational& b) const {
    return compare(b) > 0;
  }

  bool operator<=(const FastRational& b) const {
    return compare(b) <= 0;
  }

  bool operator>=(const FastRational& b) const {
    return compare(b) >= 0;
  }

  bool operator!=(const FastRational& b) const {
    return !(*this == b);
  }

  FAST_RATIONAL_INLINE unsigned size() const;

  FAST_RATIONAL_INLINE FastInteger get_num() const {
    if (has_word) return FastInteger(num);
    else return FastInteger(mpq_numref(mpq));
  }

  FAST_RATIONAL_INLINE FastInteger get_den() const {
    if (has_word) return FastInteger(den);
    else return FastInteger(mpq_numref(mpq));
  }

  bool isWellFormed() const;

  operator mpq_class() const;

#if 0
  mpz_class get_num() const {
    if (has_word) return mpz_class(num);
    mpz_class x;
    mpz_set(x.get_mpz_t(), mpq_numref(mpq));
    return x;
  }

  mpz_class get_den() const {
    if (has_word) return mpz_class(den);
    mpz_class x;
    mpz_set(x.get_mpz_t(), mpq_denref(mpq));
    return x;
  }

  operator mpq_class() const {
    if (has_word) {
      mpq_class r;
      mpz_set_si(mpq_numref(r.get_mpq_t()), num);
      mpz_set_ui(mpq_denref(r.get_mpq_t()), den);
      return r;
    }
    mpq_class r;
    mpq_set(r.get_mpq_t(), mpq);
    return r;
  }
#endif

  FastRational operator+(const FastRational& b) const {
    FastRational dest;
    assert(isWellFormed());
    assert(b.isWellFormed());
    addition(dest, *this, b);
    assert(dest.isWellFormed());
    return dest;
  }
  
  FastRational operator-(const FastRational& b) const {
    FastRational dest;
    assert(isWellFormed());
    assert(b.isWellFormed());
    subtraction(dest, *this, b);
    assert(dest.isWellFormed());  
    return dest;
  }

  FastRational operator*(const FastRational& b) const {
    FastRational dest;
    assert(isWellFormed());
    assert(b.isWellFormed());
    multiplication(dest, *this, b);
    assert(dest.isWellFormed());  
    return dest;
  }

  FastRational operator/(const FastRational& b) const {
    FastRational dest;
    assert(isWellFormed());
    assert(b.isWellFormed());
    division(dest, *this, b);
    assert(dest.isWellFormed());  
    return dest;
  }

  FastRational& operator+=(const FastRational& b) {
    assert(isWellFormed());
    assert(b.isWellFormed());
    additionAssign(*this, b);
    assert(isWellFormed());  
    return *this;
  }

  FastRational& operator-=(const FastRational& b) {
    assert(isWellFormed());
    assert(b.isWellFormed());
    subtractionAssign(*this, b);
    assert(isWellFormed());  
    return *this;
  }

  FastRational& operator*=(const FastRational& b) {
    assert(isWellFormed());
    assert(b.isWellFormed());
    multiplicationAssign(*this, b);
    assert(isWellFormed());
    return *this;
  }

  FastRational& operator/=(const FastRational& b) {
    assert(isWellFormed());
    assert(b.isWellFormed());
    divisionAssign(*this, b);
    assert(isWellFormed());
    return *this;
  }

  FAST_RATIONAL_INLINE FastRational inverse() const;

  size_t hash() const {
    if (has_word) {
      return djb2(num, den);
    } else {
      return djb2(mpz_get_si(mpq_numref(mpq)), mpz_get_si(mpq_denref(mpq)));
    }
  }

  bool isZero() const {
    if (has_word) {
      return num==0;
    } else {
      return mpq_sgn(mpq)==0;
    }
  }
};

inline std::ostream& operator<<(std::ostream& out, const FastRational& r) {
  r.print(out);
  return out;
}

inline FastRational inverse(const FastRational& a) {
  return a.inverse();
}

inline int sign(const FastRational& x) {
  return x.sign();
}

namespace std {

  template <typename object> struct hash;

  template<> struct hash<FastInteger> {
    size_t operator()(const FastInteger& x) {
      return x.hash();
    }
  };

  template<> struct hash<FastRational> {
    size_t operator()(const FastRational& x) {
      return x.hash();
    }
  };
} 
#ifdef INLINE_FAST_ARITHMETIC
#include "fastRationals.hpp"
#endif
#endif
