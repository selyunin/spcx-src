#include <sstream>
#include <algorithm>

/***
Incremental simplex algorithm

Template algorithmic code

David Monniaux, VERIMAG, 2008
$Id$
***/

/**
SIMPLEX
**/

#ifndef NDEBUG
static unsigned countBasicVariables(const std::vector<variableType>& vt) {
  unsigned count=0;
  for(unsigned k=0; k<vt.size(); k++) {
    if (vt[k] == BASIC) count++;
  }
  return count;
}
#endif

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
adjustCurrent(const std::vector<variableType>& newVariableTypes) {
  for(variableSet::iterator xkIter=nonBasicVariables.begin();
      xkIter != nonBasicVariables.end();
      xkIter++) {
    unsigned k=*xkIter;
    
    switch (newVariableTypes[k]) {
    case UPPER:
      updateNonBasic(k, state.upperBounds[k].getFiniteValue());
      break;
    case LOWER:
      updateNonBasic(k, state.lowerBounds[k].getFiniteValue());
      break;
    case ZERO:
      updateNonBasic(k, coordinate(0));
      break;
    case BASIC:
      break;
    }
  }
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
unsigned Simplex <coefficient, coordinate, upperBound, lowerBound> ::
nonBasicOccurences(unsigned nonBasic) const {
  assert(isNonBasicVariable(nonBasic));
  unsigned count = 0;
  for(variableSet::const_iterator xj = basicVariables.begin();
      xj != basicVariables.end(); xj++) {
    unsigned j = *xj;
    if (forms[j].hasNonZero(nonBasic)) {
      count++;
    }
  }
  return count;
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
allNonBasicOccurences(std::vector<std::pair<unsigned, unsigned> >& counts) const {
  for(variableSet::const_iterator xi=basicVariables.begin();
      xi != basicVariables.end(); xi++) {
    unsigned i=*xi;
    const Row& row = forms[i];
    for(typename Row::const_iterator it=row.begin(); it!=row.end(); it++) {
      if (! (*it).isZero()) {
	counts[it.index()].first++;
      }
    }
  }
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
forcePivot(const std::vector<variableType>& newVariableTypes,
	   bool checkUnsat) {
  assert(nrVariables == newVariableTypes.size());
  assert(basicVariables.size() == countBasicVariables(newVariableTypes));

  typedef std::pair<unsigned, unsigned> varP;
  typedef std::vector<varP> varV;
  varV basicVariablesToPivot;
  
  for(variableSet::iterator xj = basicVariables.begin();
      xj != basicVariables.end(); xj++) {
    unsigned j = *xj;
    if (newVariableTypes[j] != BASIC) {
      basicVariablesToPivot.push_back(varP(forms[j].filled(), j));
    }
  }
  std::sort(basicVariablesToPivot.begin(), basicVariablesToPivot.end());

#define SORT_NONBASIC
#ifdef SORT_NONBASIC
  varV nonBasicVariablesToPivot;
#if 0
  std::vector<std::pair<unsigned, unsigned> > counts(nrVariables);
  for(unsigned k=0; k<nrVariables; k++) {
    counts[k].first=0;
    counts[k].second=k;
  }
  allNonBasicOccurences(counts);
#else
  for(variableSet::iterator xk = nonBasicVariables.begin();
      xk != nonBasicVariables.end(); xk++) {
    unsigned k = *xk;
    if (newVariableTypes[k] == BASIC) {
      nonBasicVariablesToPivot.push_back(varP(nonBasicOccurences(k), k));
    }
  }
#endif
  std::sort(nonBasicVariablesToPivot.begin(), nonBasicVariablesToPivot.end());
#endif

  bool hasPivotedOnce;
  do {
    hasPivotedOnce = false;
    for(varV::const_iterator xj = basicVariablesToPivot.begin();
	xj != basicVariablesToPivot.end(); xj++) {
      unsigned j = xj->second, k=0;
      if (!isBasicVariable(j)) continue;

      bool found=false;

#ifdef SORT_NONBASIC
      for(varV::const_iterator xk = nonBasicVariablesToPivot.begin();
	  !found && xk != nonBasicVariablesToPivot.end(); xk++) {
	k = xk->second;
#else
      for(variableSet::const_iterator xk = nonBasicVariables.begin();
	  !found && xk != nonBasicVariables.end(); xk++) {
	k = *xk;
	if (newVariableTypes[k] != BASIC) continue;
#endif
	if (isBasicVariable(k)) continue;
	if (forms[j].hasNonZero(k)) found=true;
      }

      if(found) {
	pivot(j, k, checkUnsat);
	hasPivotedOnce = true;
	if (!isSatisfied()) {
	  //std::cout << "found contradiction when pivoting row " <<j<< std::endl;
	  break;
	}
      }
    }
  } while (hasPivotedOnce);

  adjustCurrent(newVariableTypes);
  recomputeBasicVariables();
  checkWellFormed();
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
std::vector<variableType>
  Simplex <coefficient, coordinate, upperBound, lowerBound> :: variableTypes() const {
  std::vector<variableType> types(nrVariables, ZERO);

  for(unsigned k=0; k<nrVariables; k++) {
    if (isBasicVariable(k)) {
      types[k] = BASIC;
    } else {
      if (state.current[k] == state.upperBounds[k]) {
	types[k] = UPPER;
      } else if (state.current[k] == state.lowerBounds[k]) {
	types[k] = LOWER;
      } else {
	assert (state.current[k] == coefficient(0));
      }
    }
  }
  return types;
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
recomputeBrokenVariableSet() {
  for(unsigned k=0; k<nrVariables; k++) {
    if (isBasicVariable(k) && (state.current[k] > state.upperBounds[k]
	  || state.current[k] < state.lowerBounds[k])) {
      brokenVariables.insert(k);
    } else {
      brokenVariables.erase(k);
    }
  }  
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
recomputeBasicVariables() {
  for(variableSet::const_iterator it=basicVariables.begin();
      it != basicVariables.end(); it++) {
    coordinate x(0);
    const Row& form=getForm(*it);

    for(typename Row::const_iterator jt=form.begin();
	jt != form.end(); jt++) {
      if ((*jt).isZero()) continue;
      x += *jt * state.current[jt.index()];
    }

    state.current[*it] = x;
  }
  recomputeBrokenVariableSet();
}

template <class coefficient,
  class coordinate,
  class upperBound,
  class lowerBound>
  void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
  resetBounds() {
  state.contradictionKind = PROBLEM;
  state.hasSolution = true;
  brokenVariables.clear();
  for(unsigned i=0; i<nrVariables; i++) {
    state.lowerBounds[i] = numeric_type<lowerBound>::minusInfinity();
    state.upperBounds[i] = numeric_type<upperBound>::plusInfinity();
  }
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
pivot(unsigned i, unsigned j, bool checkUnsat) {
  assert(isBasicVariable(i));
  assert(!isBasicVariable(j));
  assert(forms[i].hasNonZero(j));

#ifdef TRACE_PIVOT
  std::cout << numeric_type<coefficient>::name()
	    << " pivoting " << i << " and " << j
	    << " coeff=" << coefficient(forms[i][j])
	    << " low=" << state.lowerBounds[i]
	    << " cur=" << state.current[i]
	    << " high=" << state.upperBounds[i]
	    << std::endl;
#endif

  /* Express x[j] in terms of x[i], instead of the converse */
  const coefficient p = -inverse(coefficient(forms[i][j]));

  forms[j].assignTimes(p, forms[i]);
  forms[j].erase(j);

  forms[j][i]=-p;

  if (tracker != 0) {
    (*tracker)[j].assignTimes(p, (*tracker)[i]);
  }

  basicVariables.insert(j);
  nonBasicVariables.erase(j);
  nonBasicVariables.insert(i);
  basicVariables.erase(i);

  /*
    Main pivoting. This has a n^2 complexity and it's the main
    element in the computational costs.
  */

  /*
    Collecting these indices allows for better parallel balancing
    of the next loop.
   */
  std::vector<unsigned> basicVariableIndices;
  for(variableSet::iterator xuIter=basicVariables.begin();
      xuIter != basicVariables.end();
      xuIter++) {
    unsigned xu = *xuIter;
    if (xu==j) continue;
    if (forms[xu].hasNonZero(j)) {
      basicVariableIndices.push_back(xu);
    }
  }

#ifdef _OPENMP
#pragma omp parallel for if (basicVariableIndices.size()*forms[j].filled() > 1000)
#endif
  for(int z=0; z<(int) basicVariableIndices.size(); z++) {
    unsigned k = basicVariableIndices[z];
    assert (isBasicVariable(k));

    if (k == j) continue;
    const coefficient factor=forms[k][j];
    assert (factor!=0);

    forms[k].erase(j);
    forms[k].addTimes(factor, forms[j]);

    if (tracker != 0) {
      (*tracker)[k].addTimes(factor, (*tracker)[j]);
    }

  }  

  forms[i].clear();

  pivotCounter++;

  if (checkUnsat) {
    for(int z=0; z<(int) basicVariableIndices.size(); z++) {
      unsigned k = basicVariableIndices[z];
      if(checkBoundsOnRow(k)) return;
    }
    checkBoundsOnRow(j);
  }
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
pivotAndUpdate(unsigned i, unsigned j, const coordinate& newValue) {
  assert(i < nrVariables);
  assert(j < nrVariables);
  assert(isBasicVariable(i));
  assert(! isBasicVariable(j));

  assert(newValue >= state.lowerBounds[i]);
  assert(newValue <= state.upperBounds[i]);

  const coefficient pivotCoeff = forms[i][j];
  if (pivotCoeff == 0) {
    abort();
  }

  coordinate theta = (newValue-state.current[i])/pivotCoeff;
  state.current[i] = newValue;
  brokenVariables.erase(i);

  for(variableSet::iterator xkIter=basicVariables.begin();
      xkIter != basicVariables.end();
      xkIter++) {
    unsigned k=*xkIter;
    if (k == i) continue;
    typename Row::const_iterator xjIter = forms[k].find(j);
    if (xjIter == forms[k].end() || xjIter.index()!=j) continue;
    const coefficient& coeff = *xjIter;
    if (coeff.isZero()) continue;
    applyDelta(k, coeff*theta);
  }
  pivot(i, j);

  applyDelta(j, theta);
  checkGoodBasicVariable(j);

  pivotAndUpdateCounter++;
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
applyBasicVariable(unsigned j, const coordinate& newValue) {
  assert(isBasicVariable(j));
  if (newValue < state.lowerBounds[j] || newValue > state.upperBounds[j]) {
    brokenVariables.insert(j);
  } else {
    brokenVariables.erase(j);
  }
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
applyDelta(unsigned j, const coordinate& delta) {
  applyBasicVariable(j, state.current[j] += delta);
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
updateNonBasic(unsigned i, const coordinate& newValue) {
  assert(! isBasicVariable(i));
  //TODO brokenVariables.insert(i);

  coordinate delta = newValue - state.current[i];
  for(variableSet::iterator xjIter=basicVariables.begin();
      xjIter != basicVariables.end();
      xjIter++) {
    unsigned j=*xjIter;
    typename Row::const_iterator xiIter = forms[j].find(i);
    if (xiIter == forms[j].end() || xiIter.index()!=i) continue;
    const coefficient& k = *xiIter;
    if (k.isZero()) continue;
    applyDelta(j, k*delta);
  }
  state.current[i]=newValue;
}

static std::vector<std::string> *alphaVariableNames(unsigned count) {
  std::vector<std::string> *o = new std::vector<std::string>(count);
  for(unsigned i=0; i<count; i++) {
    std::ostringstream oss;
    oss << "alpha";
    oss << i;
    (*o)[i]=oss.str();
  }
  return o;
}

// Constructor with tracker (2)
template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
Simplex <coefficient, coordinate, upperBound, lowerBound> ::
Simplex(const std::vector<std::string>& variableNames0, unsigned, bool hasTracker) :
    nrVariables(variableNames0.size()),
    variableNames(&variableNames0),

    forms(nrVariables, Row (variableNames0)),

    state(nrVariables),

    trackerDimension(variableNames0.size()),
    trackerVariableNames(alphaVariableNames(trackerDimension)),
    tracker(hasTracker ?
	    new std::vector<Row >
	    (nrVariables, Row (*trackerVariableNames))
	    : 0),

#ifndef SIMPLEX_USE_SET_OF_UNSIGNED
    basicVariables(nrVariables),
    nonBasicVariables(nrVariables),
    brokenVariables(nrVariables),
#endif

    pivotAndUpdateCounter(0), pivotCounter(0)
{
  const unsigned nrVariables = variableNames0.size();

  for(unsigned i=0; i<nrVariables; i++) {
    nonBasicVariables.insert(i);
  }

  if (hasTracker) {
    for(unsigned i=0; i<nrVariables; i++) {
      (*tracker)[i][i]=1;
    }
  }
}


#if !defined(NDEBUG) && !defined(SIMPLEX_USE_INTERVAL_DOUBLE)
template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
checkGoodBasicVariable(unsigned i) const {
  assert(isBasicVariable(i));
  assert(forms[i][i]==0);
  coordinate total(0);
  for(unsigned j=0; j<nrVariables; j++) {
    if (isNonBasicVariable(j)) {
      assert(!isBasicVariable(j));
      total += forms[i][j]*state.current[j];
    } else {
      assert(isBasicVariable(j));
      assert(forms[i][j]==0);
    }
  }
#if 0
  assert(areReasonablyClose(state.current[i], total));
#endif
}
#endif

#ifndef NDEBUG
template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
checkWellFormed() const {
  /* All variables should be either basic, or nonbasic */
  for(unsigned i=0; i< nrVariables; i++) {
    assert(isBasicVariable(i) || isNonBasicVariable(i));
    assert(! (isBasicVariable(i) && isNonBasicVariable(i)));
  }

  for(variableSet::const_iterator xiIter=basicVariables.begin();
      xiIter != basicVariables.end();
      xiIter++) {
    unsigned i=*xiIter;
    checkGoodBasicVariable(i);
  }

  if (isSatisfied()) {
    for(variableSet::const_iterator xiIter=basicVariables.begin();
	xiIter != basicVariables.end();
	xiIter++) {
      unsigned i=*xiIter;
      if (!isBrokenVariable(i)) {
	assert(state.current[i] <= state.upperBounds[i]);
	assert(state.current[i] >= state.lowerBounds[i]);
      } else {
	assert(state.current[i] < state.lowerBounds[i] ||
	       state.current[i] > state.upperBounds[i]);
      }
    }
  }

  if (isSatisfied() && hasImmediateAnswer()) {
    /* Check if the answer fits the bounds */
    for(unsigned i=0; i<nrVariables; i++) {
      if(state.current[i] > state.upperBounds[i] ||
         state.current[i] < state.lowerBounds[i]) {
	std::cerr << "fail bound at " << i << std::endl;
	abort();
      }
    }
  } else {
    /* Check if bounds are coherent */
    for(unsigned i=0; i<nrVariables; i++) {
      assert(state.lowerBounds[i] <= state.upperBounds[i]);
    }
  }
}
#endif

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
bool Simplex <coefficient, coordinate, upperBound, lowerBound> ::
assertUpperBound(unsigned i, const upperBound& bound) {
  /* Asserting a wider constraint than the current one does nothing. */
  if (bound >= state.upperBounds[i]) return false;

  /* Asserting an upper bound strictly less than the lower bound
     makes the problem unsatisfiable. */
  if (bound < state.lowerBounds[i]) {
    //std::cout << "broke lower bound:" << bound << "<" << state.lowerBounds[i] << std::endl;
    foundContradiction(i, BETWEEN_BOUNDS_CONTRADICTION);
  } else {
    state.upperBounds[i] = bound;
    if (state.current[i] > bound) {
      if (isBasicVariable(i)) {
	brokenVariables.insert(i);
      } else {
	updateNonBasic(i, bound.getFiniteValue());
      }
    }
  }
  return true;
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
bool Simplex <coefficient, coordinate, upperBound, lowerBound> ::
assertLowerBound(unsigned i, const lowerBound& bound) {
  /* Asserting a wider constraint than the current one does nothing. */
  if (bound <= state.lowerBounds[i]) return false;

  /* Asserting a lower bound strictly greater than the upper bound
     makes the problem unsatisfiable. */
  if (bound > state.upperBounds[i]) {
    //std::cout << "broke upper bound:" << bound << "<" << state.upperBounds[i] << std::endl;
    foundContradiction(i, BETWEEN_BOUNDS_CONTRADICTION);
  } else {
    state.lowerBounds[i] = bound;
    if (state.current[i] < bound) {
      if (isBasicVariable(i)) {
	brokenVariables.insert(i);
      } else {
	updateNonBasic(i, bound.getFiniteValue());
      }
    }
  }
  return true;
}

/**
This is where everything happens.
**/
template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
propagate() {
  typedef std::pair<unsigned, unsigned> varP;
  typedef std::vector<varP> varV;
  varV variablesToPivot;

#if 1
  for(unsigned j=0; j<nrVariables; j++) {
    variablesToPivot.push_back(varP(
      isBasicVariable(j) ? forms[j].filled() : 0, j));
  }
  allNonBasicOccurences(variablesToPivot);
#else
  for(unsigned j=0; j<nrVariables; j++) {
    variablesToPivot.push_back(varP(
      isBasicVariable(j) ? forms[j].filled() : nonBasicOccurences(j), j));
  }
#endif

  std::sort(variablesToPivot.begin(), variablesToPivot.end());

  while (true) {
#ifdef TRACE
    std::cout << "ITERATION" << std::endl << *this << std::endl;
#endif
    checkWellFormed();
    if (!isSatisfied()) break;

    int i=-1;
    for(varV::const_iterator it = variablesToPivot.begin();
	it != variablesToPivot.end(); it++) {
      if (isBrokenVariable(it->second)) {
	i=it->second;
	break;
      }
    }
    if (i < 0) break;

    assert(isBasicVariable(i));
    brokenVariables.erase(i);

    if (state.current[i] < state.lowerBounds[i]) {
      unsigned j;
      bool found=false;
      for(varV::const_iterator xjIter = variablesToPivot.begin();
	  !found && xjIter != variablesToPivot.end();
	  xjIter++) {
	j=xjIter->second;
	if (isBasicVariable(j)) continue;

	typename Row::const_iterator xjIter= forms[i].find(j);
	if (xjIter==forms[i].end() || xjIter.index() != j) continue;
	int sgn = sign(*xjIter);

	if ((sgn > 0 && state.current[j] < state.upperBounds[j])
	    || (sgn < 0 && state.current[j] > state.lowerBounds[j])) {
	  found=true;
	  break;
	}
      }

      if (!found) {
	foundContradiction(i, LOWER_BOUND_CONTRADICTION);
	assert(cannotSatisfyLowerBound(getContradictionRow()));
      } else {
	//if (forms[i][j] == 0) abort();
	pivotAndUpdate(i, j, state.lowerBounds[i].getFiniteValue());
      }
    } else if (state.current[i] > state.upperBounds[i]) {
      unsigned j;
      bool found=false;
      for(varV::const_iterator xjIter = variablesToPivot.begin();
	  !found && xjIter != variablesToPivot.end();
	  xjIter++) {
	j=xjIter->second;
	if (isBasicVariable(j)) continue;

	typename Row::const_iterator xjIter= forms[i].find(j);
	if (xjIter==forms[i].end() || xjIter.index() != j) continue;
	int sgn = sign(*xjIter);

	if ((sgn > 0 && state.current[j] > state.lowerBounds[j])
	    || (sgn < 0 && state.current[j] < state.upperBounds[j])) {
	  found=true;
	  break;
	}
      }

      if (!found) {
	foundContradiction(i, UPPER_BOUND_CONTRADICTION);
	assert(cannotSatisfyUpperBound(getContradictionRow()));
      } else {
	//if (forms[i][j] == 0) abort();
	pivotAndUpdate(i, j, state.upperBounds[i].getFiniteValue());
      }
    } else {
      assert(false && (state.current[i] < state.lowerBounds[i] || state.current[i] > state.upperBounds[i]));
    }
  }
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
std::ostream& operator<<(std::ostream& out, const Simplex <coefficient, coordinate, upperBound, lowerBound> & simplex) {
  for(unsigned i=0; i<simplex.numberOfVariables(); i++) {
    if (simplex.isBasicVariable(i)) {
      out << simplex.variableName(i) << " = " << simplex.getForm(i);
    } else {
      out << simplex.variableName(i) << " is nonbasic";
    }
    if (simplex.hasTracker()) {
      out << " [" << simplex.getTrackerRow(i) << "]";
    }
    out << std::endl;
  }
  
  for(unsigned i=0; i<simplex.numberOfVariables(); i++) {
    out << simplex.getLowerBound(i) << " <= (" << simplex.variableName(i) << " = " << simplex.getCurrent(i) << ") <= " << simplex.getUpperBound(i) << std::endl;
  }

  out << "broken variables:";
  for(variableSet::const_iterator v=simplex.getBrokenVariables().begin();
      v!=simplex.getBrokenVariables().end(); v++) {
    out << " " << simplex.variableName(*v);
  }
  out << std::endl;

  out << "is satisfied: " << std::boolalpha << simplex.isSatisfied() << std::endl;
  if (!simplex.isSatisfied()) {
    out << "contradiction row: " << simplex.getContradictionRow()
	<< " kind: " << simplex.getContradictionKind() << std::endl;
  }
  return out;
}


template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
bool Simplex <coefficient, coordinate, upperBound, lowerBound> ::
checkTrivial() {
  for(variableSet::const_iterator it=basicVariables.begin();
      it != basicVariables.end(); it++) {
    if (checkBoundsOnRow(*it)) return false;
  }
  return true;
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
bool Simplex <coefficient, coordinate, upperBound, lowerBound> :: check()
{
  if (isSatisfied()) {
    propagate();
  }
  return isSatisfied();
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
printStatistics(std::ostream& out) const {
  out << "pivot: " << pivotCounter << ", pivot+update: " << pivotAndUpdateCounter;
  if (numeric_type<coefficient>::hasSize) {
    unsigned size=0;
    for(unsigned i=0; i<forms.size(); i++) {
      for(typename Row::const_iterator xj=forms[i].begin();
	  xj != forms[i].end(); xj++) {
	size = std::max(size, numeric_type<coefficient>::size(*xj));
      }
    }
    out << ", max size: " << size; 
  }
  out << std::endl;
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound>::setInternalState
(const Simplex <coefficient, coordinate, upperBound, lowerBound>::internalState& newState) {
  //std::cout << "before setting internal state" << std::endl;
  checkWellFormed();
  state = newState;
  recomputeBasicVariables();
  //std::cout << "after recomputing basic variables" << std::endl;
  checkWellFormed();
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
assignBasic(unsigned index, const Row& form) {
  if (index >= forms.size()) {
    std::cerr << "assigning too big a linear form" << std::endl;
    abort();
  }

  forms[index] = form;
  basicVariables.insert(index);
  nonBasicVariables.erase(index);
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
bool Simplex <coefficient, coordinate, upperBound, lowerBound> ::
cannotSatisfyUpperBound(unsigned row) const {
  if (numeric_type<upperBound>::isPlusInfinity(state.upperBounds[row])) return false;
  const Row& form = forms[row];
  coordinate sum(0);
  coefficient prod;
  for(typename Row::const_iterator xi = form.begin();
      xi != form.end(); xi++) {
    const unsigned i = xi.index();
    const coefficient& coef = *xi;
    int sgn = sign(coef);
    if (sgn > 0) {
      const lowerBound& bound = state.lowerBounds[i];
      if (numeric_type<lowerBound>::isMinusInfinity(bound)) return false;
      sum += coef * coordinate(bound);
    } else if (sgn < 0) {
      const upperBound& bound = state.upperBounds[i];
      if (numeric_type<upperBound>::isPlusInfinity(bound)) return false;
      sum += coef * coordinate(bound);
    }
  }
  if (sum > state.upperBounds[row]) {
    //std::cout << "sum = " << sum << " ubound = " << state.upperBounds[row] << std::endl;
    //std::cout << "detected beyond upper bound" << std::endl;
    //std::cout.flush();
    return true;
  } else {
    return false;
  }
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
bool Simplex <coefficient, coordinate, upperBound, lowerBound> ::
cannotSatisfyLowerBound(unsigned row) const {
  if (numeric_type<lowerBound>::isMinusInfinity(state.lowerBounds[row])) return false;
  const Row& form = forms[row];
  coordinate sum(0);
  coefficient prod;
  for(typename Row::const_iterator xi = form.begin();
      xi != form.end(); xi++) {
    const unsigned i = xi.index();
    const coefficient& coef = *xi;
    int sgn = sign(coef);
    if (sgn > 0) {
      const upperBound& bound = state.upperBounds[i];
      if (numeric_type<upperBound>::isPlusInfinity(bound)) return false;
      sum += coef * coordinate(bound);
    } else if (sgn < 0) {
      const lowerBound& bound = state.lowerBounds[i];
      if (numeric_type<lowerBound>::isMinusInfinity(bound)) return false;
      sum += coef * coordinate(bound);
    }
  }
  if (sum < state.lowerBounds[row]) {
    //std::cout << "sum = " << sum << " lbound = " << state.lowerBounds[row] << std::endl;
    //std::cout << "detected beyond lower bound" << std::endl;
    //std::cout.flush();
    return true;
  } else {
    return false;
  }
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
bool Simplex <coefficient, coordinate, upperBound, lowerBound> ::
checkBoundsOnRow(unsigned k) {
  //std::cout << "check bounds on row " << k << std::endl;
  if (!isBasicVariable(k)) return false;
  if (cannotSatisfyUpperBound(k)) {
    foundContradiction(k, UPPER_BOUND_CONTRADICTION);
    return true;
  }
  if (cannotSatisfyLowerBound(k)) {
    foundContradiction(k, LOWER_BOUND_CONTRADICTION);
    return true;
  }
  return false;
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
relaxUpperBound(unsigned i, const upperBound& bound) {
  /* Refuse to relax to a stricter constraint */
  assert (bound >= state.upperBounds[i]);

  state.upperBounds[i] = bound;
  if (isBrokenVariable(i) && state.current[i] <= bound && state.current[i] >= state.lowerBounds[i]) {
    brokenVariables.erase(i);
  }
}

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
void Simplex <coefficient, coordinate, upperBound, lowerBound> ::
relaxLowerBound(unsigned i, const lowerBound& bound) {
  /* Refuse to relax to a stricter constraint */
  assert (bound <= state.lowerBounds[i]);

  state.lowerBounds[i] = bound;
  if (isBrokenVariable(i) && state.current[i] >= bound && state.current[i] <= state.upperBounds[i]) {
    brokenVariables.erase(i);
  }
}
