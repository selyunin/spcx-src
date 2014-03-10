/***
Incremental simplex algorithm: float/rational layered implementation

Template instantiations and code

David Monniaux, VERIMAG, 2008
$Id$
***/

#include "fastSimplex.hh"
#include "simplex.hpp"
#include "timings.h"

// explicit instantiation
template class LinearForm<coefficientQ>;
template class Simplex<coefficientQ, coordinateQ, upperBoundQ, lowerBoundQ>;

template std::ostream& operator<<(std::ostream& out, const LinearForm<coefficientQ>& form);
template std::ostream& operator<<(std::ostream& out, const Simplex<coefficientQ, coordinateQ, upperBoundQ, lowerBoundQ>& simplex);

#ifndef USE_GLPK
template class LinearForm<coefficientF>;
template class Simplex<coefficientF, coordinateF, upperBoundF, lowerBoundF>;
#endif

template std::ostream& operator<<(std::ostream& out, const LinearForm<coefficientF>& form);
template std::ostream& operator<<(std::ostream& out, const Simplex<coefficientF, coordinateF, upperBoundF, lowerBoundF>& simplex);

// Implementation
void FastSimplex::assignBasic(unsigned index, const LinearForm<coefficientQ>& form)
{
   SimplexQ::assignBasic(index, form);
   simplexF.assignBasic(index, LinearForm<coefficientF>(form));
}

#ifdef USE_GLPK
static inline upperBoundF upperBoundQtoF(const upperBoundQ& bound) {
  if (bound.isPlusInfinity()) {
    return std::numeric_limits<double>::infinity();
  } else {
    const coordinateQ x(bound);
    double y= x.standardValue.get_d();
    if (x.epsilonCoefficient < 0) {
      return nextafter(y, -std::numeric_limits<double>::infinity());
    } else {
      return y;
    }
  }
}

static inline lowerBoundF lowerBoundQtoF(const lowerBoundQ& bound) {
  if (bound.isMinusInfinity()) {
    return -std::numeric_limits<double>::infinity();
  } else {
    const coordinateQ x(bound);
    double y= x.standardValue.get_d();
    if (x.epsilonCoefficient > 0) {
      return nextafter(y, std::numeric_limits<double>::infinity());
    } else {
      return y;
    }
  }
}

#else

static inline upperBoundF upperBoundQtoF(const upperBoundQ& bound) {
  if (bound.isPlusInfinity()) {
    return upperBoundF(std::numeric_limits<double>::infinity());
  } else {
    const coordinateQ x(bound);
    return upperBoundF(x.standardValue.get_d(), x.epsilonCoefficient.get_d());
  }
}

static inline lowerBoundF lowerBoundQtoF(const lowerBoundQ& bound) {
  if (bound.isMinusInfinity()) {
    return lowerBoundF(-std::numeric_limits<double>::infinity());
  } else {
    const coordinateQ x(bound);
    return lowerBoundF(x.standardValue.get_d(), x.epsilonCoefficient.get_d());
  }
}
#endif

bool FastSimplex::assertUpperBound(unsigned i, const upperBound& bound) {
  if (SimplexQ::assertUpperBound(i, bound)) {
    simplexF.setBounds(i, lowerBoundQtoF(getLowerBound(i)),
		       upperBoundQtoF(bound));
    return true;
  } else {
    return false;
  }
}

bool FastSimplex::assertLowerBound(unsigned i, const lowerBound& bound) {
  if (SimplexQ::assertLowerBound(i, bound)) {
    simplexF.setBounds(i, lowerBoundQtoF(bound),
		       upperBoundQtoF(getUpperBound(i)));
    return true;
  } else {
    return false;
 }
}

#if 0
void FastSimplex::relaxUpperBound(unsigned i, const upperBound& bound) {
  simplexF.relaxUpperBound(i, upperBoundQtoF(bound));
  SimplexQ::relaxUpperBound(i, bound);
}

void FastSimplex::relaxLowerBound(unsigned i, const lowerBound& bound) {
  simplexF.relaxLowerBound(i, lowerBoundQtoF(bound));
  SimplexQ::relaxLowerBound(i, bound);
}
#endif

bool FastSimplex::reduceContradiction() {
  TIMING_DECLARE();
  TIMING_BEGINS();
    
  unsigned constraintCount=0;
  bool hasReduced = false;

  for(unsigned i=0; i<numberOfVariables(); i++) {
    if (! getUpperBound(i).isPlusInfinity()) {
      simplexF.setBounds(i, lowerBoundQtoF(getLowerBound(i)), numeric_type<double>::plusInfinity());
      if (simplexF.check()) {
	simplexF.setBounds(i, lowerBoundQtoF(getLowerBound(i)), upperBoundQtoF(getUpperBound(i)));
	constraintCount++;
      } else {
	SimplexQ::relaxUpperBound(i, numeric_type<upperBound>::plusInfinity());
	hasReduced = true;
      }
    }
    
    if (! getLowerBound(i).isMinusInfinity()) {
      simplexF.setBounds(i, numeric_type<double>::minusInfinity(), upperBoundQtoF(getUpperBound(i)));
      if (simplexF.check()) {
	simplexF.setBounds(i, lowerBoundQtoF(getLowerBound(i)), upperBoundQtoF(getUpperBound(i)));
	constraintCount++;
      } else {
	SimplexQ::relaxLowerBound(i, numeric_type<lowerBound>::minusInfinity());
	hasReduced = true;
      }
    }
    }
    
  std::cout << "constraints: " << constraintCount << std::endl;
  TIMING_ENDS(reduceTimeF);
  return hasReduced;
}

bool FastSimplex::fastCheck() {
  bool ret;
  TIMING_DECLARE();
  if (hasImmediateAnswer()) {
    return SimplexQ::check();
  }
  TIMING_BEGINS();
  ret = simplexF.fastCheck();
  TIMING_ENDS(checkTimeF);
  if (ret) return true;

  //reduceContradiction();
  
  TIMING_BEGINS();
  SimplexQ::forcePivot(simplexF.variableTypes(), true);
  TIMING_ENDS(transferTimeFtoQ);
  
  TIMING_BEGINS();
  ret= SimplexQ::check();
  TIMING_ENDS(checkTimeQ);
  
  return ret;
}

bool FastSimplex::check() {
  TIMING_DECLARE();

  TIMING_BEGINS();
  bool ret = checkTrivial();
  TIMING_ENDS(precheck);
  if (!ret) return false;

  if (hasImmediateAnswer()) {
    return SimplexQ::check();
  }
  TIMING_BEGINS();
  simplexF.check();
  TIMING_ENDS(checkTimeF);

  //reduceContradiction();
  
  TIMING_BEGINS();
  SimplexQ::forcePivot(simplexF.variableTypes(), true);
  TIMING_ENDS(transferTimeFtoQ);
  
  TIMING_BEGINS();
  ret= SimplexQ::check();
  TIMING_ENDS(checkTimeQ);
  
  return ret;
}

std::ostream& operator<<(std::ostream& out, const FastSimplex& s) {
  out << "FLOAT SIMPLEX:";
  out << s.getSimplexF();
  out << "EXACT SIMPLEX:";
  out << s.getSimplexQ();
  return out;
}

void FastSimplex::printStatistics(std::ostream& out) const {
  out << "FLOAT SIMPLEX: ";
  getSimplexF().printStatistics(out);
  out << "EXACT SIMPLEX: ";
  getSimplexQ().printStatistics(out);
#ifdef TIMINGS
  out << "Times: precheck=" << precheck << " checkTimeF=" << checkTimeF << " reduceTimeF=" << reduceTimeF << ", transferTimeFtoQ=" << transferTimeFtoQ << ", checkTimeQ=" << checkTimeQ << std::endl;
#endif
}

void FastSimplex::setInternalState(const internalState& state) {
  SimplexQ::setInternalState(state);

  for(unsigned i=0; i<numberOfVariables(); i++) {
    simplexF.setBounds(i, lowerBoundQtoF(getLowerBound(i)), upperBoundQtoF(getUpperBound(i)));
  }
}

void FastSimplex::resetBounds() {
  SimplexQ::resetBounds();

  for(unsigned i=0; i<numberOfVariables(); i++) {
    simplexF.setBounds(i,
		       -std::numeric_limits<double>::infinity(),
		       +std::numeric_limits<double>::infinity());
  }
}
