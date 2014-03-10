#include <stdlib.h>
#include <assert.h>

FAST_RATIONAL_INLINE void FastInteger::init_read_string(const char *s) {
  errno = 0;
  num = strtol(s, 0, 10);
  if (errno != ERANGE) {
    has_word = true;
    has_mpz = false;
  } else {
    has_word = false;
    mpz_init_set_str(mpz, s, 10);
    has_mpz = true;
  }
}

#if 0
FAST_RATIONAL_INLINE FastInteger::FastInteger(const mpz_class& x) {
  if (mpz_fits_slong_p(x.get_mpz_t())) {
    num = mpz_get_si(x.get_mpz_t());
    has_word = true;
    has_mpz = false;
  } else {
    has_word = false;
    has_mpz = true;
    mpz_init_set(mpz, x.get_mpz_t());
  }
}
#endif

FAST_RATIONAL_INLINE FastInteger::FastInteger(const mpz_t x) {
  if (mpz_fits_slong_p(x)) {
    num = mpz_get_si(x);
    has_word = true;
    has_mpz = false;
  } else {
    has_word = false;
    has_mpz = true;
    mpz_init_set(mpz, x);
  }
}

FAST_RATIONAL_INLINE FastInteger::FastInteger(const FastInteger& x) {
  num = x.num;
  has_word = x.has_word;
  if ( (has_mpz = x.has_mpz) ) {
    mpz_init(mpz);
    mpz_set(mpz, x.mpz);
  }
}

FAST_RATIONAL_INLINE FastInteger& FastInteger::operator=(const FastInteger& x) {
  num = x.num;
  has_word = x.has_word;
  if (x.has_mpz) {
    if (! has_mpz) {
      mpz_init(mpz);
      has_mpz = true;
    }
    mpz_set(mpz, x.mpz);
  } else {
    if (has_mpz) {
      mpz_clear(mpz);
      has_mpz = false;
    }
  }
  return *this;
}

FAST_RATIONAL_INLINE bool FastInteger::operator==(const FastInteger& b) const {
  if (has_word && b.has_word) {
    return num == b.num;
  }
  force_make_mpz();
  b.force_make_mpz();
  return mpz_cmp(mpz, b.mpz)==0;
}

FAST_RATIONAL_INLINE FastInteger FastInteger::operator-() const {
  if (has_word && num > LONG_MIN) {
    return FastInteger(-num);
  } else {
    force_make_mpz();
    FastInteger x;
    x.has_word = false;
    x.has_mpz = true;
    mpz_init(x.mpz);
    mpz_neg(x.mpz, mpz);
    return x;
  }
}

FAST_RATIONAL_INLINE int FastInteger::compare(const FastInteger& b) const {
  if (has_word && b.has_word) {
    return compare(num, b.num);
  }
  force_make_mpz();
  b.force_make_mpz();
  return mpz_cmp(mpz, b.mpz);
}

FAST_RATIONAL_INLINE int FastInteger::sign() const {
  if (has_word) {
    if (num < 0) return -1;
    else if (num > 0) return 1;
    else return 0;
  } else {
    assert(has_mpz);
    return mpz_sgn(mpz);
  }
}

FAST_RATIONAL_INLINE FastRational::FastRational(const FastInteger& x) {
  if (x.has_word) {
    num = x.num;
    den = 1;
    has_word = true;
    has_mpq = false;
  } else {
    assert(x.has_mpz);
    has_word = false;
    has_mpq = true;
    mpq_init(mpq);
    mpz_set(mpq_numref(mpq), x.mpz);
  }
}

FAST_RATIONAL_INLINE FastRational::FastRational(const FastRational& x) {
  num = x.num;
  den = x.den;
  has_mpq = false;
  if (! (has_word = x.has_word)) {
    assert(x.has_mpq);
    has_mpq = true;
    mpq_init(mpq);
    mpq_set(mpq, x.mpq);
  }
}

FAST_RATIONAL_INLINE FastRational& FastRational::operator=(const FastRational& x) {
  num = x.num;
  den = x.den;
  has_word = x.has_word;
  if (!has_word) {
    assert(x.has_mpq);
    if (!has_mpq) {
      mpq_init(mpq);
      has_mpq = true;
    }
    mpq_set(mpq, x.mpq);
  } else {
    if (has_mpq) {
      mpq_clear(mpq);
      has_mpq = false;
    }
  }
  return *this;
}

FAST_RATIONAL_INLINE bool FastRational::operator==(const FastRational& b) const {
  if (has_word && b.has_word) {
    return num == b.num && den == b.den;
  }
  force_make_mpq();
  b.force_make_mpq();
  return mpq_equal(mpq, b.mpq);
}

FAST_RATIONAL_INLINE FastRational FastRational::operator-() const {
  if (has_word && num > LONG_MIN) {
      return FastRational(-num, den);
  } else {
    force_make_mpq();
    FastRational x;
    x.has_word = false;
    x.has_mpq = true;
    mpq_init(x.mpq);
    mpq_neg(x.mpq, mpq);
    return x;
  }
}

FAST_RATIONAL_INLINE int FastRational::compare(const FastRational& b) const {
  if (has_word && b.has_word) {
    if (b.den == den) {
      return compare(num, b.num);
    } else {
      return compare(lword(num)*b.den, lword(b.num)*den);
    }
  }
  force_make_mpq();
  b.force_make_mpq();
  return mpq_cmp(mpq, b.mpq);
}

FAST_RATIONAL_INLINE int FastRational::sign() const {
  if (has_word) {
    if (num < 0) return -1;
    else if (num > 0) return 1;
    else return 0;
  } else {
    assert(has_mpq);
    return mpq_sgn(mpq);
  }
}

inline uword absVal(word x) {
  return x>=0 ? x : -x;
}

inline ulword absVal(lword x) {
  return x>=0 ? x : -x;
}

#ifdef FAST_GCD
uword gcd(uword a, uword b);
ulword gcd(ulword a, ulword b);

#else
#ifdef BINARY_GCD
template<typename integer> integer gcd(integer a, integer b) {
  if(a==0) return b;
  if(b==0) return a;

  if (a > b) {
    a = a%b;
  } else if (a < b) {
    b = b%a;
  } else return a;

  if(a==0) return b;
  if(b==0) return a;

  unsigned shift=0;
  while((a&1)==0 && (b&1)==0) {
    a>>=1;
    b>>=1;
    shift++;
  }
  while((a&1)==0) a>>=1;

  while(b>0) {
    while((b&1)==0) b>>=1;
    if (b>a) {
      b-=a;
    } else {
      integer d=a-b;
      a=b;
      b=d;
    }
  }

  return a<<shift;
}
template<ulword> ulword gcd(ulword a, ulword b);
template<uword> uword gcd(uword a, uword b);

#else
template<typename integer> integer gcd(integer a, integer b) {
  if(a==0) return b;
  if(b==0) return a;

  if (b > a) {
    integer c = a;
    a = b;
    b = c;
  }
  while(true) {
    integer r = a%b;
    if (r==0) return b;
    a = b;
    b = r;
  }
}
template<ulword> ulword gcd(ulword a, ulword b);
template<uword> uword gcd(uword a, uword b);
#endif
#endif

#define CHECK_WORD(var, value) \
  do { \
    lword tmp = value; \
    if (tmp < LONG_MIN/2 || tmp > LONG_MAX/2) { \
      goto overflow; \
    } \
    var = tmp;\
  } while(0)

#ifdef NDEBUG
#define CHECK_POSITIVE(value)
#else
#define CHECK_POSITIVE(value) \
  if (value < 1) abort()
#endif

#define CHECK_UWORD(var, value) \
  do { \
    CHECK_POSITIVE(value); \
    ulword tmp = value; \
    if (tmp > LONG_MAX) { \
      goto overflow; \
    } \
    var = tmp;\
  } while(0)

#define COMPUTE_WORD(var, value) \
  word var; CHECK_WORD(var, value)
#define COMPUTE_UWORD(var, value) \
  uword var; CHECK_UWORD(var, value)

#define INTEGER_OP(name, op, mpzop) \
FAST_RATIONAL_INLINE void name(FastInteger& dst, const FastInteger& a, const FastInteger& b){ \
  if (a.has_word && b.has_word) { \
    word zn; \
    CHECK_WORD(zn, lword(a.num) op lword(b.num)); \
    dst.num = zn; \
    return ; \
  } \
overflow: \
  a.force_make_mpz(); \
  b.force_make_mpz(); \
  dst.make_erase_mpz(); \
  mpz_##mpzop(dst.mpz, a.mpz, b.mpz); \
  dst.make_word(); \
} \
\
FAST_RATIONAL_INLINE void name##Assign(FastInteger& a, const FastInteger& b){ \
  if (a.has_word && b.has_word) { \
    word zn; \
    CHECK_WORD(zn, lword(a.num) op lword(b.num)); \
    a.num = zn; \
    return ; \
  } \
overflow: \
  a.make_mpz(); \
  b.force_make_mpz(); \
  mpz_##mpzop(a.mpz, a.mpz, b.mpz); \
  a.make_word(); \
}

INTEGER_OP(addition, +, add)
INTEGER_OP(subtraction, -, sub)
INTEGER_OP(multiplication, *, mul)
INTEGER_OP(division, /, cdiv_q)
INTEGER_OP(modulo, %, cdiv_r)
INTEGER_OP(divexact, /, divexact)

FAST_RATIONAL_INLINE void lcm(FastInteger& dst, const FastInteger& a, const FastInteger& b){
  if (a.has_word && b.has_word) {
    word zn;
    uword aAbs = absVal(a.num), bAbs = absVal(b.num);
    CHECK_WORD(zn, (ulword(aAbs / gcd(aAbs, bAbs)))*bAbs);
    dst.num = zn;
    dst.kill_mpz();
    return ;
  }
overflow:
  a.force_make_mpz();
  b.force_make_mpz();
  dst.make_erase_mpz();
  mpz_lcm(dst.mpz, a.mpz, b.mpz);
  dst.make_word();
}

FAST_RATIONAL_INLINE void gcd(FastInteger& dst, const FastInteger& a, const FastInteger& b){
  if (a.has_word && b.has_word) {
    dst.num = gcd(absVal(a.num), absVal(b.num));
    dst.kill_mpz();
    return ;
  }

  a.force_make_mpz();
  b.force_make_mpz();
  dst.make_erase_mpz();
  mpz_gcd(dst.mpz, a.mpz, b.mpz);
  dst.make_word();
}

FAST_RATIONAL_INLINE bool FastRational::isWellFormed() const
{
    return (has_word || has_mpq) &&
      (!has_word || (den!=0 && gcd(absVal(num), den)==1))
      && (!has_mpq || mpz_sgn(mpq_denref(mpq))!=0);
}

FAST_RATIONAL_INLINE FastRational::FastRational(word n, uword d) : has_mpq(false), has_word(true) {
  assert(d > 0);
  if (n == 0) {
    num = 0;
    den = 1;
  } else if (n > 0) {
    word common = gcd(uword(n), d);
    num = n/common;
    den = d/common;
  } else {
    word common = gcd(uword(-n), d);
    num = n/common;
    den = d/common;
  }
}

FAST_RATIONAL_INLINE void addition(FastRational& dst, const FastRational& a, const FastRational& b){
  if (a.has_word && b.has_word) {
    if (b.num == 0) {
      dst.num = a.num;
      dst.den = a.den;
    } else if (a.num == 0) {
      dst.num = b.num;
      dst.den = b.den;
    } else if (b.den == 1) {
      CHECK_WORD(dst.num, lword(a.num) + lword(b.num)*a.den);
      dst.den = a.den;
    } else if (a.den == 1) {
      CHECK_WORD(dst.num, lword(b.num) + lword(a.num)*b.den);
      dst.den = b.den;
    } else {
      lword n = lword(a.num)*b.den + lword(b.num)*a.den;
      ulword d = ulword(a.den) * b.den;
      lword common = gcd(absVal(n), d);
      word zn;
      uword zd;
      if (common > 1) {
	CHECK_WORD(zn, n/common);
	CHECK_UWORD(zd, d/common);
      } else {
	CHECK_WORD(zn, n);
	CHECK_UWORD(zd, d);
      }
      dst.num = zn;
      dst.den = zd;
    }
    dst.kill_mpq();
    return;
  }

 overflow:
  a.force_make_mpq();
  b.force_make_mpq();
  dst.make_erase_mpq();
  mpq_add(dst.mpq, a.mpq, b.mpq);
  dst.make_word();
}

FAST_RATIONAL_INLINE void subtraction(FastRational& dst, const FastRational& a, const FastRational& b){
  if (a.has_word && b.has_word) {
    if (b.num == 0) {
      dst.num = a.num;
      dst.den = a.den;
    } else if (a.num == 0) {
      dst.num = -b.num;
      dst.den = b.den;
    } else if (b.den == 1) {
      CHECK_WORD(dst.num, lword(a.num) - lword(b.num)*a.den);
      dst.den = a.den;
    } else if (a.den == 1) {
      CHECK_WORD(dst.num, -lword(b.num) + lword(a.num)*b.den);
      dst.den = b.den;
    } else {
      lword n = lword(a.num)*b.den - lword(b.num)*a.den;
      ulword d = ulword(a.den) * b.den;
      lword common = gcd(absVal(n), d);
      word zn;
      uword zd;
      if (common > 1) {
	CHECK_WORD(zn, n/common);
	CHECK_UWORD(zd, d/common);
      } else {
	CHECK_WORD(zn, n);
	CHECK_UWORD(zd, d);
      }
      dst.num = zn;
      dst.den = zd;
    }
    dst.kill_mpq();
    return;
  }

 overflow:
  a.force_make_mpq();
  b.force_make_mpq();
  dst.make_erase_mpq();
  mpq_sub(dst.mpq, a.mpq, b.mpq);
  dst.make_word();
}

FAST_RATIONAL_INLINE void multiplication(FastRational& dst, const FastRational& a, const FastRational& b){
  if ((a.has_word && a.num==0) || (b.has_word && b.num==0)) {
    dst.num=0;
    dst.den=1;
    dst.has_word = true;
    dst.kill_mpq();
    return;
  }

  if (a.has_word && a.num==1 && a.den==1) {
    dst = b;
    return;
  }
  if (b.has_word && b.num==1 && b.den==1) {
    dst = a;
    return;
  }

  if (a.has_word && b.has_word) {
    word zn;
    uword zd;
#if 1
    word common1 = gcd(absVal(a.num), b.den), common2 = gcd(a.den, absVal(b.num));
    lword k1, k2, k3, k4;
    if (common1 > 1) {
      k1 = lword(a.num)/common1;
      k4 = ulword(b.den)/common1;
    } else {
      k1 = lword(a.num);
      k4 = ulword(b.den);
    }
    if (common2 > 1) {
      k2 = lword(b.num)/common2;
      k3 = ulword(a.den)/common2;
    } else {
      k2 = lword(b.num);
      k3 = ulword(a.den);
    }
    CHECK_WORD(zn, k1 * k2);
    CHECK_UWORD(zd, k3 * k4);
#else
    lword zn1 = lword(a.num)*b.num;
    ulword zd1 = ulword(a.den)*b.den;
    ulword common = gcd(absVal(zn1), zd1);
    if (common > 1) {
      zn1 /= common;
      zd1 /= common;
    }      
    CHECK_WORD(zn, zn1);
    CHECK_UWORD(zd, zd1);
#endif
    dst.num = zn;
    dst.den = zd;
    dst.has_word = true;
    dst.kill_mpq();
    return;
  }

 overflow:
  a.force_make_mpq();
  b.force_make_mpq();
  dst.make_erase_mpq();
  mpq_mul(dst.mpq, a.mpq, b.mpq);
  dst.make_word();
}

FAST_RATIONAL_INLINE void division(FastRational& dst, const FastRational& a, const FastRational& b){
  if (a.has_word && b.has_word) {
    word common1 = gcd(absVal(a.num), absVal(b.num));
    word common2 = gcd(a.den, b.den);
    word zn;
    uword zd;
    CHECK_WORD(zn, (lword(a.num)/common1) * (b.den/common2));
    CHECK_UWORD(zd, ulword(absVal(b.num)/common1) * (a.den/common2));
    dst.num = zn;
    dst.den = zd;
    if (b.num < 0) dst.num=-dst.num;
    dst.has_word = true;
    dst.kill_mpq();
    return;
  }

 overflow:
  a.force_make_mpq();
  b.force_make_mpq();
  dst.make_erase_mpq();
  mpq_div(dst.mpq, a.mpq, b.mpq);
  dst.make_word();
}

FAST_RATIONAL_INLINE double FastRational::get_d() const {
  if (has_word) {
    return double(num)/double(den);
  } else {
    assert(has_mpq);
    return mpq_get_d(mpq);
  }
}


FAST_RATIONAL_INLINE void additionAssign(FastRational& a, const FastRational& b){
  if (b.has_word) {
    if (b.num == 0) return;
    if (a.has_word) {
      if (b.den == 1) {
	CHECK_WORD(a.num, lword(a.num) + lword(b.num)*a.den);
      } else if (a.num == 0) {
	a.num = b.num;
	a.den = b.den;
      } else {
	lword n = lword(a.num)*b.den + lword(b.num)*a.den;
	ulword d = ulword(a.den) * b.den;
	lword common = gcd(absVal(n), d);
	word zn;
	uword zd;
	if (common > 1) {
	  CHECK_WORD(zn, n/common);
	  CHECK_UWORD(zd, d/common);
	} else {
	  CHECK_WORD(zn, n);
	  CHECK_UWORD(zd, d);
	}
	a.num = zn;
	a.den = zd;
      }
      a.kill_mpq();
      return;
    }
  }

 overflow:
  a.make_mpq();
  b.force_make_mpq();
  mpq_add(a.mpq, a.mpq, b.mpq);
  a.make_word();
}

FAST_RATIONAL_INLINE void subtractionAssign(FastRational& a, const FastRational& b){
  if (a.has_word && b.has_word) {
    word common = gcd(a.den, b.den);
    COMPUTE_WORD(n1, lword(a.num) * (b.den / common));
    COMPUTE_WORD(n2, lword(b.num) * (a.den / common));
    lword n = lword(n1) - lword(n2);
    ulword d = ulword(a.den) * (b.den / common);
    common = gcd(absVal(n), d);
    word zn;
    uword zd;
    CHECK_WORD(zn, n/common);
    CHECK_UWORD(zd, d/common);
    a.num = zn;
    a.den = zd;
    a.kill_mpq();
    return;
  }

 overflow:
  a.make_mpq();
  b.force_make_mpq();
  mpq_sub(a.mpq, a.mpq, b.mpq);
  a.make_word();
}

FAST_RATIONAL_INLINE void multiplicationAssign(FastRational& a, const FastRational& b){
  if (a.has_word && b.has_word) {
    word common1 = gcd(absVal(a.num), b.den);
    word common2 = gcd(a.den, absVal(b.num));
    word zn;
    uword zd;
    CHECK_WORD(zn, lword(a.num/common1) * (b.num/common2));
    CHECK_UWORD(zd, ulword(a.den/common2) * (b.den/common1));
    a.num = zn;
    a.den = zd;
    a.kill_mpq();
    return;
  }

 overflow:
  a.make_mpq();
  b.force_make_mpq();
  mpq_mul(a.mpq, a.mpq, b.mpq);
  a.make_word();
}

FAST_RATIONAL_INLINE void divisionAssign(FastRational& a, const FastRational& b){
  if (a.has_word && b.has_word) {
    word common1 = gcd(absVal(a.num), absVal(b.num)),
      common2 = gcd(a.den, b.den);
    word zn;
    uword zd;
    if (b.num < 0) {
      CHECK_WORD(zn, -lword(a.num/common1) * (b.den/common2));
    } else {
      CHECK_WORD(zn, lword(a.num/common1) * (b.den/common2));
    }
    CHECK_UWORD(zd, ulword(absVal(b.num)/common1) * (a.den/common2));
    a.den = zd;
    a.num = zn;
    a.kill_mpq();
    return;
  }

 overflow:
  a.make_mpq();
  b.force_make_mpq();
  mpq_div(a.mpq, a.mpq, b.mpq);
  a.make_word();
}

FAST_RATIONAL_INLINE unsigned FastInteger::size() const {
  if (has_word) return 32;
  return mpz_sizeinbase(mpz, 2);
}

FAST_RATIONAL_INLINE unsigned FastRational::size() const {
  if (has_word) return 64;
  return mpz_sizeinbase(mpq_numref(mpq), 2) +
    mpz_sizeinbase(mpq_denref(mpq), 2);
}

FAST_RATIONAL_INLINE FastRational FastRational::inverse() const {
  FastRational dest;
  if (has_word) {
    assert(num != 0);
    word zn;
    uword zd;
    if (num > 0) {
      CHECK_WORD(zn, den);
      CHECK_UWORD(zd, num);
    } else {
      CHECK_WORD(zn, -lword(den));
      CHECK_UWORD(zd, -lword(num));
    }
    dest.num = zn;
    dest.den = zd;
    return dest;
  }
 overflow:
  dest.has_word = false;
  force_make_mpq();
  dest.make_erase_mpq();
  mpq_inv(dest.mpq, mpq);
  return dest;
}

FAST_RATIONAL_INLINE FastRational::operator mpq_class() const {
  if (has_word) {
    return mpq_class(num)/den;
  } else {
    assert(has_mpq);
    return mpq_class(mpq);
  }
}
