#include "integerSimplex.hh"

integer floor(const mpq_class& x, bool& isExact) {
  integer q, r;
  mpz_fdiv_qr(q.get_mpz_t(), r.get_mpz_t(), x.get_num().get_mpz_t(), x.get_den().get_mpz_t());
  isExact = (r == 0);
  return q;
}

integer floor(const mpq_class& x) {
  integer q;
  mpz_fdiv_q(q.get_mpz_t(), x.get_num().get_mpz_t(), x.get_den().get_mpz_t());
  return q;
}

integer floorInfinitesimal(const coordinate& x) {
  mpq_class low;
  bool isExact;

  low = floor(x.standardValue, isExact);
  if (!isExact || x.epsilonCoefficient >= 0) return low;
  else return (low - 1);
}

upperBound floorBound(const upperBound& x) {
  if (x.isPlusInfinity()) {
    return x;
  } else {
    return upperBound(coefficient(floorInfinitesimal(x.getFiniteValue())));
  }
}

void IntegerSimplex::assertUpperBound(unsigned i, const upperBound& bound) {
  if (getIntegerFlag(i)) {
    MY_SIMPLEX::assertUpperBound(i, floorBound(bound));
  } else {
    MY_SIMPLEX::assertUpperBound(i, bound);
  }
}

void IntegerSimplex::relaxUpperBound(unsigned i, const upperBound& bound) {
  if (getIntegerFlag(i)) {
    MY_SIMPLEX::relaxUpperBound(i, floorBound(bound));
  } else {
    MY_SIMPLEX::relaxUpperBound(i, bound);
  }
}

integer ceil(const mpq_class& x, bool& isExact) {
  integer q, r;
  mpz_cdiv_qr(q.get_mpz_t(), r.get_mpz_t(), x.get_num().get_mpz_t(), x.get_den().get_mpz_t());
  isExact = (r == 0);
  return q;
}

integer ceil(const mpq_class& x) {
  integer q;
  mpz_cdiv_q(q.get_mpz_t(), x.get_num().get_mpz_t(), x.get_den().get_mpz_t());
  return q;
}

integer ceilInfinitesimal(const coordinate& x) {
  mpq_class high;
  bool isExact;

  high = ceil(x.standardValue, isExact);
  if (!isExact || x.epsilonCoefficient <= 0) return high;
  else return (high + 1);
}

lowerBound ceilBound(const lowerBound& x) {
  if (x.isMinusInfinity()) {
    return x;
  } else {
    return lowerBound(coefficient(ceilInfinitesimal(x.getFiniteValue())));
  }
}

void IntegerSimplex::assertLowerBound(unsigned i, const lowerBound& bound) {
  if (getIntegerFlag(i)) {
    MY_SIMPLEX::assertLowerBound(i, ceilBound(bound));
  } else {
    MY_SIMPLEX::assertLowerBound(i, bound);
  }
}

void IntegerSimplex::relaxLowerBound(unsigned i, const lowerBound& bound) {
  if (getIntegerFlag(i)) {
    MY_SIMPLEX::relaxLowerBound(i, ceilBound(bound));
  } else {
    MY_SIMPLEX::relaxLowerBound(i, bound);
  }
}

bool IntegerSimplex::checkIfBadIndex(unsigned i, integer& lowSplit) {
  if (!isIntegerVariable[i]) return false;
  const coordinate& c = getCurrent(i);

  bool isExact;
  lowSplit = floor(c.standardValue, isExact);
  if (!isExact) return true;
  if (c.epsilonCoefficient >= 0) {
    lowSplit--;
  }
  return false;
}

bool IntegerSimplex::findBadIndex(unsigned& badIndex, integer& lowSplit) {
  for(unsigned i=0; i<isIntegerVariable.size(); i++) {
    if (checkIfBadIndex(i, lowSplit)) {
      badIndex = i;
      return true;
    }
  }
  return false;
}

// CURRENTLY INCORRECT
bool IntegerSimplex::check() {
  // If then rational relaxation has no solution, then the mixed system has no solution.
  if (! MY_SIMPLEX::check()) return false;

  unsigned badIndex;
  integer lowSplit;
  if (! findBadIndex(badIndex, lowSplit)) {
    // no integer variables have noninteger assignments
    return true;
  }

  integer highSplit = lowSplit+1;

#ifdef TRACE
  std::cout << "split: ("
	    << variableName(badIndex) << " <= " << lowSplit << ") \\/ ("
	    << variableName(badIndex) << " >= " << highSplit << ")"
	    << std::endl;
#endif

  {
    upperBound savedUpperBound = getUpperBound(badIndex);
    MY_SIMPLEX::assertUpperBound(badIndex, coordinate(coefficient(lowSplit)));
    bool lowSplitAnswer = check();
    MY_SIMPLEX::relaxUpperBound(badIndex, savedUpperBound);
    //std::cout << "END OF SPLIT" << std::endl << *this << std::endl;

    if (lowSplitAnswer) return true;
  }
  
  {
    lowerBound savedLowerBound = getLowerBound(badIndex);
    //std::cout << "BEFORE ASSERT" << std::endl << *this << std::endl;
    MY_SIMPLEX::assertLowerBound(badIndex, coordinate(coefficient(highSplit)));
    //std::cout << "AFTER ASSERT" << std::endl << *this;
    bool highSplitAnswer = check();
    MY_SIMPLEX::relaxLowerBound(badIndex, savedLowerBound);
    
    return highSplitAnswer;
  }

}

void IntegerSimplex::printStatistics(std::ostream& out) {
  MY_SIMPLEX::printStatistics(out);
}
