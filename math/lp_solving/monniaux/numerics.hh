/***
Fancy numerical template types:
* <type> with added +infinity
* <type> with added -infinity
* x + epsilon.y with epsilon an infinitesimal

David Monniaux, VERIMAG, 2008
$Id$
 ***/

#define BOOST_VECTORS

#ifndef _NUMERICS_HH
#define _NUMERICS_HH

//#define _ISOC99_SOURCE

#include <gmp.h>
#include <gmpxx.h>

#include <limits>
#include <cassert>
#include <iostream>

#include <math.h>
#include "stdlib.h"

#ifdef FAST_RATIONALS
#include "fastRationals.hh"
#endif

#ifdef SIMPLEX_USE_INTERVAL_DOUBLE
#include <boost/numeric/interval.hpp>

using namespace boost::numeric;
using namespace interval_lib;
#endif


template<class number> struct numeric_type {
};

template <> struct numeric_type<double> {
  static double plusInfinity() {
    return std::numeric_limits<double>::infinity();
  }

  static double minusInfinity() {
    return -std::numeric_limits<double>::infinity();
  }

  static bool isPlusInfinity(double x) {
    return isinf(x)==1;
  }

  static bool isMinusInfinity(double x) {
    return isinf(x)==-1;
  }

  static bool isFinite(double x) {
    return isfinite(x);
  }

  static const bool hasSize = false;
  static unsigned size(double) {
    return 1;
  }

  static const char* name() { return "double"; }
  static double epsilon() { return 1e-13; }
  static double minus_epsilon() { return -1e-13; }
};

template<> struct numeric_type<mpq_class> {
  const static bool hasSize = true;
  static unsigned size(const mpq_class& x) {
    return mpz_sizeinbase(x.get_num().get_mpz_t(), 2)
      + mpz_sizeinbase(x.get_den().get_mpz_t(), 2);
  }
  static const char* name() { return "mpq_class"; }
  static mpq_class epsilon() { return mpq_class(0); }
  static mpq_class minus_epsilon() { return mpq_class(0); }
};

#ifdef FAST_RATIONALS
template<> struct numeric_type<FastRational> {
  const static bool hasSize = true;
  static unsigned size(const FastRational& x) {
    return x.size();
  }
  static const char* name() { return "FastRational"; }
  static FastRational epsilon() { return FastRational(); }
  static FastRational minus_epsilon() { return FastRational(); }
};
#endif
inline bool areReasonablyClose(double x, const double y) {
  return (fabs(x-y) < 1E-6);
}

#ifdef SIMPLEX_USE_INTERVAL_DOUBLE
template<class scalar> bool areReasonablyClose(const interval<scalar>& a, const interval<scalar>& b) {
  return areReasonablyClose(a.lower(), b.upper()) &&
    areReasonablyClose(a.upper(), b.lower());
}
#endif

/**
INFINITIES
**/

template <class finiteNumber> class WithPlusInfinity {
  finiteNumber finiteValue;
  bool isInfinite;

private:
  WithPlusInfinity(const finiteNumber& finiteValue0, bool isInfinite0) :
    finiteValue(finiteValue0), isInfinite(isInfinite0) {
  }

public:
  bool isPlusInfinity() const {
    return isInfinite;
  }

  WithPlusInfinity(const finiteNumber& finiteValue0) :
    finiteValue(finiteValue0), isInfinite(false) {
  }

  static const WithPlusInfinity<finiteNumber> plusInfinity() {
    return WithPlusInfinity(finiteNumber(0), true);
  }

  bool operator<=(const WithPlusInfinity<finiteNumber>& x) const {
    if (x.isInfinite) return true;
    if (isInfinite) return false;
    return finiteValue <= x.finiteValue;
  }

  bool operator<(const WithPlusInfinity<finiteNumber>& x) const {
    if (x.isInfinite) return true;
    if (isInfinite) return false;
    return finiteValue < x.finiteValue;
  }

  bool operator>=(const WithPlusInfinity<finiteNumber>& x) const {
    return x <= *this;
  }

  bool operator>(const WithPlusInfinity<finiteNumber>& x) const {
    return x < *this;
  }

  bool operator==(const finiteNumber& x) const {
    if (isInfinite) return false;
    return finiteValue == x;
  }

  bool operator<=(const finiteNumber& x) const {
    if (isInfinite) return false;
    return finiteValue <= x;
  }

  bool operator<(const finiteNumber& x) const {
    if (isInfinite) return false;
    return finiteValue < x;
  }

  bool operator>=(const finiteNumber& x) const {
    if (isInfinite) return true;
    return finiteValue >= x;
  }

  bool operator>(const finiteNumber& x) const {
    if (isInfinite) return true;
    return finiteValue > x;
  }

  operator finiteNumber() const {
    assert(! isInfinite);
    return finiteValue;
  }

  const finiteNumber& getFiniteValue() const {
    assert(! isInfinite);
    return finiteValue;
  }

  finiteNumber& getFiniteValue() {
    assert(! isInfinite);
    return finiteValue;
  }
};

template <class finiteNumber> inline bool operator==(const finiteNumber& x, const WithPlusInfinity<finiteNumber>& y) {
  return y == x;
}

template <class finiteNumber> inline bool operator<=(const finiteNumber& x, const WithPlusInfinity<finiteNumber>& y) {
  return y >= x;
}

template <class finiteNumber> inline bool operator>=(const finiteNumber& x, const WithPlusInfinity<finiteNumber>& y) {
  return y <= x;
}

template <class finiteNumber> inline bool operator<(const finiteNumber& x, const WithPlusInfinity<finiteNumber>& y) {
  return y > x;
}

template <class finiteNumber> inline bool operator>(const finiteNumber& x, const WithPlusInfinity<finiteNumber>& y) {
  return y < x;
}

template <class finiteNumber> inline std::ostream& operator<<(std::ostream& out, const WithPlusInfinity<finiteNumber>& x) {
  if (x.isPlusInfinity()) out << "+INF";
  else out << (finiteNumber) x;
  return out;
}

template <class finiteNumber> struct numeric_type<WithPlusInfinity<finiteNumber> > {
  const static WithPlusInfinity<finiteNumber> plusInfinity() {
    return WithPlusInfinity<finiteNumber>::plusInfinity();
  }

  static bool isPlusInfinity(const WithPlusInfinity<finiteNumber>& x){
    return x.isPlusInfinity();
  }

  static bool isFinite(const WithPlusInfinity<finiteNumber>& x){
    return !x.isPlusInfinity();
  }

  const static bool hasSize = numeric_type<finiteNumber>::hasSize;
  unsigned static size(WithPlusInfinity<finiteNumber> x) {
    return numeric_type<finiteNumber>::size(x);
  }
};

template <class finiteNumber> class WithMinusInfinity {
  finiteNumber finiteValue;
  bool isInfinite;

private:
  WithMinusInfinity(const finiteNumber& finiteValue0, bool isInfinite0) :
    finiteValue(finiteValue0), isInfinite(isInfinite0) {
  }

public:
  bool isMinusInfinity() const {
    return isInfinite;
  }

  WithMinusInfinity(const finiteNumber& finiteValue0) :
    finiteValue(finiteValue0), isInfinite(false) {
  }

  static const WithMinusInfinity<finiteNumber> minusInfinity() {
    return WithMinusInfinity(finiteNumber(0), true);
  }

  bool operator>=(const WithMinusInfinity<finiteNumber>& x) const {
    if (x.isInfinite) return true;
    if (isInfinite) return false;
    return finiteValue >= x.finiteValue;
  }

  bool operator>(const WithMinusInfinity<finiteNumber>& x) const {
    if (x.isInfinite) return true;
    if (isInfinite) return false;
    return finiteValue > x.finiteValue;
  }

  bool operator<=(const WithMinusInfinity<finiteNumber>& x) const {
    return x >= *this;
  }

  bool operator<(const WithMinusInfinity<finiteNumber>& x) const {
    return x > *this;
  }

  bool operator==(const finiteNumber& x) const {
    if (isInfinite) return false;
    return finiteValue == x;
  }

  bool operator>=(const finiteNumber& x) const {
    if (isInfinite) return false;
    return finiteValue >= x;
  }

  bool operator>(const finiteNumber& x) const {
    if (isInfinite) return false;
    return finiteValue > x;
  }

  bool operator<=(const finiteNumber& x) const {
    if (isInfinite) return true;
    return finiteValue <= x;
  }

  bool operator<(const finiteNumber& x) const {
    if (isInfinite) return true;
    return finiteValue < x;
  }

  operator finiteNumber() const {
    assert(! isInfinite);
    return finiteValue;
  }

  const finiteNumber& getFiniteValue() const {
    assert(! isInfinite);
    return finiteValue;
  }

  finiteNumber& getFiniteValue() {
    assert(! isInfinite);
    return finiteValue;
  }
};

template <class finiteNumber> inline bool operator==(const finiteNumber& x, const WithMinusInfinity<finiteNumber>& y) {
  return y == x;
}

template <class finiteNumber> inline bool operator<=(const finiteNumber& x, const WithMinusInfinity<finiteNumber>& y) {
  return y >= x;
}

template <class finiteNumber> inline bool operator>=(const finiteNumber& x, const WithMinusInfinity<finiteNumber>& y) {
  return y <= x;
}

template <class finiteNumber> inline bool operator<(const finiteNumber& x, const WithMinusInfinity<finiteNumber>& y) {
  return y > x;
}

template <class finiteNumber> inline bool operator>(const finiteNumber& x, const WithMinusInfinity<finiteNumber>& y) {
  return y < x;
}

template <class finiteNumber> inline std::ostream& operator<<(std::ostream& out, const WithMinusInfinity<finiteNumber>& x) {
  if (x.isMinusInfinity()) out << "-INF";
  else out << (finiteNumber) x;
  return out;
}

template <class finiteNumber> struct numeric_type<WithMinusInfinity<finiteNumber> > {
  const static WithMinusInfinity<finiteNumber> minusInfinity() {
    return WithMinusInfinity<finiteNumber>::minusInfinity();
  }

  static bool isMinusInfinity(const WithMinusInfinity<finiteNumber>& x){
    return x.isMinusInfinity();
  }

  static bool isFinite(const WithMinusInfinity<finiteNumber>& x){
    return !x.isMinusInfinity();
  }

  const static bool hasSize = numeric_type<finiteNumber>::hasSize;
  unsigned static size(WithMinusInfinity<finiteNumber> x) {
    return numeric_type<finiteNumber>::size(x);
  }
};

template <class finiteNumber> inline bool operator<=(const WithMinusInfinity<finiteNumber>& x,
					const WithPlusInfinity<finiteNumber>& y) {
  if (x.isMinusInfinity() || y.isPlusInfinity()) return true;
  return (finiteNumber) x <= (finiteNumber) y;
}

template <class finiteNumber> inline bool operator>=(const WithPlusInfinity<finiteNumber>& y, const WithMinusInfinity<finiteNumber>& x) {
  return x <= y;
}

template <class finiteNumber> inline bool operator<(const WithMinusInfinity<finiteNumber>& x,
					const WithPlusInfinity<finiteNumber>& y) {
  if (x.isMinusInfinity() || y.isPlusInfinity()) return true;
  return (finiteNumber) x < (finiteNumber) y;
}

template <class finiteNumber> inline bool operator>(const WithPlusInfinity<finiteNumber>& y, const WithMinusInfinity<finiteNumber>& x) {
  return x < y;
}

template <class finiteNumber> inline bool operator<=(const WithPlusInfinity<finiteNumber>& x,
					const WithMinusInfinity<finiteNumber>& y) {
  if (x.isPlusInfinity() || y.isMinusInfinity()) return false;
  return (finiteNumber) x <= (finiteNumber) y;
}

template <class finiteNumber> inline bool operator>=(const WithMinusInfinity<finiteNumber>& y, const WithPlusInfinity<finiteNumber>& x) {
  return x <= y;
}

template <class finiteNumber> inline bool operator<(const WithPlusInfinity<finiteNumber>& x,
					const WithMinusInfinity<finiteNumber>& y) {
  if (x.isPlusInfinity() || y.isMinusInfinity()) return false;
  return (finiteNumber) x < (finiteNumber) y;
}

template <class finiteNumber> inline bool operator>(const WithMinusInfinity<finiteNumber>& y, const WithPlusInfinity<finiteNumber>& x) {
  return x < y;
}
/**
x + k epsilon
**/

template <class standardNumber> struct WithInfinitesimal{
  standardNumber standardValue, epsilonCoefficient;

  WithInfinitesimal() {
  }

  WithInfinitesimal(const standardNumber& v) :
    standardValue(v), epsilonCoefficient(0) {
  }

  WithInfinitesimal(const standardNumber& v, const standardNumber& eps) :
    standardValue(v), epsilonCoefficient(eps) {
  }

  WithInfinitesimal<standardNumber> operator+(const WithInfinitesimal <standardNumber>& v) const {
    return WithInfinitesimal<standardNumber>(standardValue+v.standardValue, epsilonCoefficient+v.epsilonCoefficient);
  }

  WithInfinitesimal<standardNumber> operator-(const WithInfinitesimal <standardNumber>& v) const {
    return WithInfinitesimal<standardNumber>(standardValue-v.standardValue, epsilonCoefficient-v.epsilonCoefficient);
  }

  WithInfinitesimal<standardNumber> operator*(const standardNumber& s) const {
    return WithInfinitesimal<standardNumber>(standardValue*s, epsilonCoefficient*s);
  }

  WithInfinitesimal<standardNumber> operator/(const standardNumber& s) const {
    return WithInfinitesimal<standardNumber>(standardValue/s, epsilonCoefficient/s);
  }

  WithInfinitesimal<standardNumber>& operator*=(const standardNumber& s) {
    standardValue *= s;
    epsilonCoefficient *= s;
    return *this;
  }

  WithInfinitesimal<standardNumber>& operator+=(const WithInfinitesimal<standardNumber>& s) {
    standardValue += s.standardValue;
    epsilonCoefficient += s.epsilonCoefficient;
    return *this;
  }

  WithInfinitesimal<standardNumber>& operator-=(const WithInfinitesimal<standardNumber>& s) {
    standardValue -= s.standardValue;
    epsilonCoefficient -= s.epsilonCoefficient;
    return *this;
  }

  WithInfinitesimal<standardNumber> operator-() const {
    return WithInfinitesimal<standardNumber>(-standardValue, -epsilonCoefficient);
  }

  const WithInfinitesimal<standardNumber>& getFiniteValue() const {
    return *this;
  }
};

template<class standardNumber>
inline bool operator==(const WithInfinitesimal<standardNumber>& u,
		       const standardNumber& v) {
    return u.standardValue==v && u.epsilonCoefficient==0;
}

template<class standardNumber>
inline bool operator==(const WithInfinitesimal<standardNumber>& u,
		       const WithInfinitesimal<standardNumber>& v) {
    return u.standardValue==v.standardValue && u.epsilonCoefficient==v.epsilonCoefficient;
}

template<class standardNumber>
inline bool operator!=(const WithInfinitesimal<standardNumber>& u,
		       const WithInfinitesimal<standardNumber>& v) {
  return !(u==v);
}

template<class standardNumber>
inline bool operator<=(const WithInfinitesimal<standardNumber>& u,
		       const WithInfinitesimal<standardNumber>& v) {
  if (u.standardValue < v.standardValue) return true;
  if (u.standardValue > v.standardValue) return false;
  return (u.epsilonCoefficient <= v.epsilonCoefficient);
}

template<class standardNumber>
inline bool operator<(const WithInfinitesimal<standardNumber>& u,
		      const WithInfinitesimal<standardNumber>& v) {
  if (u.standardValue < v.standardValue) return true;
  if (u.standardValue > v.standardValue) return false;
  return (u.epsilonCoefficient < v.epsilonCoefficient);
}

template<class standardNumber>
inline bool operator>=(const WithInfinitesimal<standardNumber>& u,
		       const WithInfinitesimal<standardNumber>& v) {
  return (v<=u);
}

template<class standardNumber>
inline bool operator>(const WithInfinitesimal<standardNumber>& u,
		      const WithInfinitesimal<standardNumber>& v) {
  return (v<u);
}

template <class standardNumber> struct numeric_type<WithInfinitesimal<standardNumber> > {
  static WithInfinitesimal<standardNumber> plusInfinity() {
    return WithInfinitesimal<standardNumber>(numeric_type<standardNumber>::plusInfinity());
  }
  static WithInfinitesimal<standardNumber> minusInfinity() {
    return WithInfinitesimal<standardNumber>(numeric_type<standardNumber>::minusInfinity());
  }

  static bool isPlusInfinity(const WithInfinitesimal<standardNumber>& x) {
    return numeric_type<standardNumber>::isPlusInfinity(x.standardNumber);
  }

  static bool isMinusInfinity(const WithInfinitesimal<standardNumber>& x) {
    return numeric_type<standardNumber>::isMinusInfinity(x.standardNumber);
  }

  const static bool hasSize = numeric_type<standardNumber>::hasSize;
  static unsigned size(const numeric_type<WithInfinitesimal<standardNumber> >& x) {
    return x.size();
  }
};

template <class standardNumber> inline WithInfinitesimal<standardNumber> operator*(const standardNumber& s, const WithInfinitesimal<standardNumber>& v) {
  return v*s;
}

template <class standardNumber> inline std::ostream& operator<<(std::ostream& out, WithInfinitesimal<standardNumber> v) {
  if (v.epsilonCoefficient != 0) {
    out << "(+ " << v.standardValue << " (* " << v.epsilonCoefficient << " epsilon))";
  } else {
    out << v.standardValue;
  }
  return out;
}

template <class number> inline bool areReasonablyClose(const WithInfinitesimal<number>& x, const WithInfinitesimal<number>& y) {
  return areReasonablyClose(x.standardValue, y.standardValue) &&
    areReasonablyClose(x.epsilonCoefficient, y.epsilonCoefficient);
}

#ifdef FAST_RATIONALS
typedef FastInteger integer;
typedef FastRational rational;
#else
typedef mpz_class integer;
typedef mpq_class rational;
#endif

inline double inverse(double x) {
  return 1./x;
}

inline mpq_class inverse(const mpq_class& x) {
  return 1/x;
}

inline int sign(const mpq_class& x) {
  return mpq_sgn(x.get_mpq_t());
} 

inline int sign(double x) {
  if (x < 0.) return -1;  else if (x > 0.) return 1;
  else return 0;
}

#if defined(__GNUC__) && defined(VECTOR_OPS)
#include "vecNumerics.hh"
#endif

#endif
