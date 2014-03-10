#ifndef _CONSTRAINTS_HH
#define _CONSTRAINTS_HH

#include "numerics.hh"
#include "linearForm.hh"
#ifdef FAST_RATIONALS
#include "fastRationals.hh"
#endif
#include <iostream>

enum comparator { LESS_OR_EQUAL,
		  GREATER_OR_EQUAL };
inline comparator operator-(comparator op) {
  switch (op) {
  case LESS_OR_EQUAL: return GREATER_OR_EQUAL;
  case GREATER_OR_EQUAL: return LESS_OR_EQUAL;
  default: assert(false);
  }
}

std::ostream& operator<<(std::ostream& out, comparator op);

class LinearConstraint {
public:
  
  typedef LinearForm<integer> linearPartType;

private:
  linearPartType lhs;
  comparator op;
  rational rhs;

  static bool leadingSignIsNegative(const LinearForm<integer>& lhs) {
    for(LinearForm<integer>::const_iterator it=lhs.begin();
	it != lhs.end(); it++) {
#ifdef FAST_RATIONALS
      const int s = (*it).sign();
      if (s < 0) return true;
      else if (s > 0) return false;
#else
      if (*it < 0) return true;
      if (*it > 0) return false;
#endif
    }
    return false;
  }

public:
  bool operator==(const LinearConstraint& x) const {
    return getComparator() == x.getComparator() &&
      rhsPart() == x.rhsPart() &&
      linearPart() == x.linearPart(); 
  }

  LinearConstraint(const LinearForm<rational>& lhs0,
		   comparator op0,
		   const rational& rhs0);

  const linearPartType& linearPart() const {
    return lhs;
  }

  const rational& rhsPart() const {
    return rhs;
  }

  comparator getComparator() const {
    return op;
  }

  static LinearConstraint genComparison(const AffineForm<rational> &a,
				  comparator op,
				  const AffineForm<rational> &b);

private:
  static inline size_t djb2(size_t a, size_t b) {
    return (a << 5) + a + b;
  }

public:
  size_t hash() const {
	  // gff
return 0;//    return djb2(djb2(std::hash<LinearForm<integer> >()(linearPart()), getComparator()), std::hash<rational>()(rhsPart()));
  }
};

std::ostream& operator<<(std::ostream& out, const LinearConstraint& ct);

namespace std {
  template<> struct less<LinearConstraint> {
    bool operator() (const LinearConstraint& x,
	             const LinearConstraint& y);
  };
};

#endif
