/***
Incremental simplex algorithm: float/rational layered implementation

Declarations

David Monniaux, VERIMAG, 2008
$Id$
***/

#include "simplex.hh"

#define USE_GLPK
#ifdef USE_GLPK
#include "GLPKSimplex.hh"
#endif

typedef rational coefficientQ;
typedef WithInfinitesimal<coefficientQ> coordinateQ;
typedef WithPlusInfinity<coordinateQ> upperBoundQ;
typedef WithMinusInfinity<coordinateQ> lowerBoundQ;

inline double numericConversion(const coefficientQ& q) {
  return q.get_d();
}

inline bool areReasonablyClose(const coefficientQ& x, const coefficientQ& y) {
  return x==y;
}

typedef Simplex<coefficientQ, coordinateQ, upperBoundQ, lowerBoundQ> SimplexQ;

#ifdef USE_GLPK
typedef double coefficientF, coordinateF, upperBoundF, lowerBoundF;
typedef GLPKSimplex SimplexF;
#else
typedef double coefficientF;
typedef WithInfinitesimal<coefficientF> coordinateF;
typedef coordinateF upperBoundF;
typedef coordinateF lowerBoundF;
typedef Simplex<coefficientF, coordinateF, upperBoundF, lowerBoundF> SimplexF;
#endif

typedef coefficientQ coefficient;
typedef coordinateQ coordinate;
typedef upperBoundQ upperBound;
typedef lowerBoundQ lowerBound;

class FastSimplex : SimplexQ {
  SimplexF simplexF;

#ifdef TIMINGS
  double checkTimeQ, checkTimeF, transferTimeFtoQ, precheck, reduceTimeF;
#endif

public:
  void resetBounds();

  typedef SimplexQ::Row Row;
  typedef SimplexQ::internalState internalState;

  const internalState& getInternalState() const {
    return SimplexQ::getInternalState();
  }

  void setInternalState(const internalState& state);

  const std::string& variableName(unsigned i) const {
    return SimplexQ::variableName(i);
  }

  const SimplexF& getSimplexF() const {
    return simplexF;
  }

  const SimplexQ& getSimplexQ() const {
    return *((SimplexQ*) this);
  }

  FastSimplex(const std::vector<std::string>& variableNames, unsigned trueVariables, bool hasTracker) :
    SimplexQ(variableNames, trueVariables, hasTracker),
    simplexF(variableNames, trueVariables)
  {
#ifdef TIMINGS
    checkTimeQ = checkTimeF = transferTimeFtoQ = precheck = reduceTimeF = 0.0;
#endif
  }

  unsigned numberOfVariables() const {
    return SimplexQ::numberOfVariables();
  }

  void assignBasic(unsigned index, const Row& form);

  void checkWellFormed() const {
    simplexF.checkWellFormed();
    SimplexQ::checkWellFormed();
  }

  bool assertUpperBound(unsigned i, const upperBound& bound);
  bool assertLowerBound(unsigned i, const lowerBound& bound);
#if 0
  void relaxUpperBound(unsigned i, const upperBound& bound);
  void relaxLowerBound(unsigned i, const lowerBound& bound);
#endif

  const Row & getTrackerRow(unsigned i) const {
    return SimplexQ::getTrackerRow(i);
  }

  const Row & getForm(unsigned i) const {
    return SimplexQ::getForm(i);
  }

  bool isSatisfied() const {
    return SimplexQ::isSatisfied();
  }

  bool isBasicVariable(unsigned i) const {
    return SimplexQ::isBasicVariable(i);
  }
  bool isNonBasicVariable(unsigned i) const {
    return SimplexQ::isNonBasicVariable(i);
  }
  unsigned getContradictionRow() const {
    return SimplexQ::getContradictionRow();
  }

  ContradictionKind getContradictionKind() const {
    return SimplexQ::getContradictionKind();
  }

  const upperBound& getUpperBound(unsigned i) const {
    return SimplexQ::getUpperBound(i);
  }

  const lowerBound& getLowerBound(unsigned i) const {
    return SimplexQ::getLowerBound(i);
  }

  const coordinate& getCurrent(unsigned i) const {
    return SimplexQ::getCurrent(i);
  }

#if 0
  void forcePivot(const std::vector<variableType>& coordinateTypes) {
    simplexF.forcePivot(coordinateTypes);
    SimplexQ::forcePivot(coordinateTypes);
  }
#endif

  bool hasImmediateAnswer() const {
    return SimplexQ::hasImmediateAnswer();
  }
  bool check();
  bool fastCheck();
  bool reduceContradiction();

  void printStatistics(std::ostream& out) const;

  unsigned nonBasicOccurences(unsigned nonBasic) const {
    return SimplexQ::nonBasicOccurences(nonBasic);
  }
};

std::ostream& operator<<(std::ostream& out, const FastSimplex& s);
