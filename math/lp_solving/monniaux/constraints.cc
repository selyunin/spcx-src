#include "constraints.hh"

std::ostream& operator<<(std::ostream& out, comparator op) {
  switch (op) {
  case LESS_OR_EQUAL: out << "<="; break;
  case GREATER_OR_EQUAL: out << ">="; break;
  }
  return out;
}

#ifndef FAST_RATIONALS
inline void lcm(mpz_class& dst, const mpz_class &a, const mpz_class &b) {
  mpz_lcm(dst.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
}

inline void gcd(mpz_class& dst, const mpz_class &a, const mpz_class &b) {
  mpz_gcd(dst.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
}

inline void divexact(mpz_class& dst, const mpz_class &a, const mpz_class &b) {
  mpz_divexact(dst.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
}
#endif

static integer getCommonDenominator(const LinearForm<rational>& lhs) {
  integer commonDenominator(1);
  for(LinearForm<rational>::const_iterator xi=lhs.begin();
      xi != lhs.end(); xi++) {
    lcm(commonDenominator, commonDenominator, (*xi).get_den());
  }
  assert(commonDenominator != 0);
  return commonDenominator;
}

static integer getCommonFactor(const LinearForm<rational>& lhs) {
  LinearForm<rational>::const_iterator xi=lhs.begin();
  assert((*xi).get_den() == 1);
  for(; xi != lhs.end() && sign(*xi)==0; xi++) ;
  if (xi == lhs.end()) return 1;
  integer commonFactor((*xi).get_num());
  for(; xi!=lhs.end(); xi++) {
    if (commonFactor == 1) break;
    assert((*xi).get_den() == 1);
    if ((*xi).get_num() != 0) {
      gcd(commonFactor, commonFactor, (*xi).get_num());
    }
  }
  assert(commonFactor != 0);
  return commonFactor;
}

LinearConstraint::LinearConstraint(const LinearForm<rational>& lhs0,
				   comparator op0,
				   const rational& rhs0) :
    lhs(lhs0.getVariableNames()), op(op0) {

    // Reduce nonconstant terms to common denominator and scale
    // the equation to only have integer nonconstant coefficients.
  integer commonDenominator = getCommonDenominator(lhs0);
  if (commonDenominator != 1) {
#ifdef BOOST_VECTORS
    for(LinearForm<rational>::const_iterator xi=
	  lhs0.begin(); xi != lhs0.end(); xi++) {
      unsigned i=xi.index();
      const rational& k = *xi; 
#else
    for(unsigned i=0; i<lhs0.size(); i++) {
      const rational& k = lhs0[i];
#endif
      integer factor;
      divexact(factor, commonDenominator, k.get_den());
      lhs[i] = k.get_num() * factor;
    }
    
    rhs = rational(rhs0) * commonDenominator;
  } else {
    integer commonFactor = getCommonFactor(lhs0);
    
    if (commonFactor != 1) {
#ifdef BOOST_VECTORS
      for(LinearForm<rational>::const_iterator xi=lhs0.begin();
	  xi != lhs0.end(); xi++) {
        const rational& k = *xi; 
	const unsigned i = xi.index();
#else
      for(unsigned i=0; i<lhs0.size(); i++) {
        const rational& k = lhs0[i];
#endif
	integer z;
	divexact(z, rational(k).get_num(), commonFactor);
	lhs[i] = z;
      }
      
      rhs = rational(rhs0) / commonFactor;
    } else {
#ifdef BOOST_VECTORS
      for(LinearForm<rational>::const_iterator xi=lhs0.begin();
	  xi != lhs0.end(); xi++) {
        const rational& k = *xi; 
	const unsigned i = xi.index();
#else
      for(unsigned i=0; i<lhs0.size(); i++) {
        const rational& k = lhs0[i];
#endif
	lhs[i] = rational(k).get_num();
      }
      rhs = rational(rhs0);
    }
  }
  
  if (leadingSignIsNegative(lhs)) {
    lhs = -lhs;
    rhs = -rhs;
    op = -op;
  }
}

template class LinearForm<integer>;
template std::ostream& operator<<(std::ostream& out, const LinearForm<integer>& form);

std::ostream& operator<<(std::ostream& out, const LinearConstraint& ct) {
  out << "(" << ct.getComparator() << " " << ct.linearPart() << " " << ct.rhsPart() << ")";
  return out;
}

LinearConstraint LinearConstraint::
  genComparison(const AffineForm<rational> &a,
		comparator op,
		const AffineForm<rational> &b) {
  AffineForm<rational> delta = a-b;
  return LinearConstraint(delta.linearPart(), op, -delta.constantPart());
}

template class std::less<LinearForm<integer> >;

bool std::less<LinearConstraint>::operator() (const LinearConstraint& x,
					      const LinearConstraint& y) {
  if (x.getComparator() < y.getComparator()) return true;
  if (y.getComparator() < x.getComparator()) return false;
  if (std::less<LinearForm<integer> >() (x.linearPart(), y.linearPart()))
    return true;
  if (std::less<LinearForm<integer> >() (y.linearPart(), x.linearPart()))
    return false;
  return x.rhsPart() < y.rhsPart();
}
