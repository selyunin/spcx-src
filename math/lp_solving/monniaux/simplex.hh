/***
Incremental simplex algorithm

Declarations

David Monniaux, VERIMAG, 2008
$Id$
***/

#ifndef _SIMPLEX_HH
#define _SIMPLEX_HH

#include "numerics.hh"
#include <vector>
#include <set>
#include <iostream>
#include "linearForm.hh"

template <class coefficient> std::ostream& operator<<(std::ostream& out, const LinearForm<coefficient>& form);

/**
SIMPLEX
**/
enum variableType { BASIC, ZERO, LOWER, UPPER };

enum ContradictionKind { PROBLEM, LOWER_BOUND_CONTRADICTION, UPPER_BOUND_CONTRADICTION, BETWEEN_BOUNDS_CONTRADICTION };

#ifdef SIMPLEX_USE_SET_OF_UNSIGNED
typedef std::unordered_set<unsigned> variableSet;
#else
#include "unsignedSet.hh"
typedef unsignedSet variableSet;
#endif

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
class Simplex {
  unsigned nrVariables;
  const std::vector<std::string>* variableNames;

  std::vector<LinearForm<coefficient> > forms;

  class _internalState {
    std::vector<lowerBound> lowerBounds;
    std::vector<upperBound> upperBounds;
    std::vector<coordinate> current;
    friend class Simplex<coefficient, coordinate, upperBound, lowerBound>;

    unsigned contradictionRow;
    ContradictionKind contradictionKind;
    bool hasSolution;

    _internalState(unsigned nrVariables) :
      lowerBounds(nrVariables, numeric_type<lowerBound>::minusInfinity()),
      upperBounds(nrVariables, numeric_type<upperBound>::plusInfinity()),
      current(nrVariables, coefficient(0)),
      hasSolution(true) {
    }
  };

  _internalState state;

  void foundContradiction(unsigned row, ContradictionKind kind) {
    state.hasSolution = false;
    state.contradictionRow = row;
    state.contradictionKind = kind;
  }

public:
  void resetBounds();

  typedef LinearForm<coefficient> Row;

  typedef _internalState internalState;

  unsigned getContradictionRow() const {
    return state.contradictionRow;
  }

  ContradictionKind getContradictionKind() const {
    return state.contradictionKind;
  }

  const internalState& getInternalState() const {
    return state;
  }

  void setInternalState(const internalState& newState);

  unsigned nonBasicOccurences(unsigned nonBasic) const;
  void allNonBasicOccurences(std::vector<std::pair<unsigned, unsigned> >& counts) const;

private:
  unsigned trackerDimension;
  std::vector<std::string> *trackerVariableNames;
  std::vector<Row > *tracker;

  variableSet basicVariables, nonBasicVariables, brokenVariables;

  unsigned pivotAndUpdateCounter, pivotCounter;

public:
  Simplex(const std::vector<std::string>& variableNames0, unsigned trueVariables, bool hasTracker=false);

  ~Simplex() {
    if (tracker!=0) {
      delete tracker;
      delete trackerVariableNames;
    }
  }

  std::vector<variableType> variableTypes() const;

  bool hasBrokenVariables() const {
    return (brokenVariables.begin() != brokenVariables.end());
  }

#ifdef NDEBUG
  void checkWellFormed() const { }
#else
  void checkWellFormed() const;
#endif

  void assignBasic(unsigned index, const Row& form);

  const Row& getForm(unsigned i) const {
    return forms[i];
  }

  unsigned numberOfVariables() const {
    return nrVariables;
  }

  bool isBasicVariable(unsigned i) const {
    return (basicVariables.find(i)!=basicVariables.end());
  }

  bool isNonBasicVariable(unsigned i) const {
    return (nonBasicVariables.find(i)!=nonBasicVariables.end());
  }

  bool isBrokenVariable(unsigned i) const {
    return (brokenVariables.find(i)!=brokenVariables.end());
  }

  const std::string& variableName(unsigned i) const {
    return (*variableNames)[i];
  }
  
  bool isSatisfied() const {
    return state.hasSolution;
  }

  const lowerBound& getLowerBound(unsigned i) const {
    return state.lowerBounds[i];
  }

  const upperBound& getUpperBound(unsigned i) const {
    return state.upperBounds[i];
  }

  const coordinate& getCurrent(unsigned i) const {
    return state.current[i];
  }

  /*
    The tracker answer:
    |alpha[i]| = coefficient for inequality[i]
    if positive, use upper bound inequality
    if negative, use lower bound inequality
   */
  const Row & getTrackerRow(unsigned i) const {
    assert(tracker!=0);
    return (*tracker)[i];
  }

  bool hasTracker() const {
    return tracker != 0;
  }

  bool assertUpperBound(unsigned i, const upperBound& bound);
  bool assertLowerBound(unsigned i, const lowerBound& bound);
  void relaxUpperBound(unsigned i, const upperBound& bound);
  void relaxLowerBound(unsigned i, const lowerBound& bound);

  void adjustCurrent(const std::vector<variableType>& newVariableTypes);
  void forcePivot(const std::vector<variableType>& coordinateTypes,
		  bool checkUnsat=false);
  /*void forcePivotOld(const std::vector<variableType>& coordinateTypes,
    bool checkUnsat=false);*/

  const variableSet& getBrokenVariables() const {
    return brokenVariables;
  }

  bool check();
  bool fastCheck() {
    return check();
  }

  bool checkTrivial();

  bool hasImmediateAnswer() const {
    return !hasBrokenVariables();
  }

  void printStatistics(std::ostream& out) const;

  bool cannotSatisfyUpperBound(unsigned row) const;
  bool cannotSatisfyLowerBound(unsigned row) const;
  bool checkBoundsOnRow(unsigned k);

private:
  void applyDelta(unsigned j, const coordinate& newValue);
  void applyBasicVariable(unsigned j, const coordinate& newValue);
  void updateNonBasic(unsigned i, const coordinate& newValue);
  void pivot(unsigned i, unsigned j, bool checkUnsat = false);
  void pivotAndUpdate(unsigned i, unsigned j, const coordinate& newValue);
#if !defined(NDEBUG) && !defined(SIMPLEX_USE_INTERVAL_DOUBLE)
  void checkGoodBasicVariable(unsigned i) const;
#else
    void checkGoodBasicVariable(unsigned) const {}
#endif
  void recomputeBrokenVariableSet();
  void recomputeBasicVariables();
  void propagate();
};

template <class coefficient,
	  class coordinate,
	  class upperBound,
	  class lowerBound>
std::ostream& operator<<(std::ostream& out,
  const Simplex<coefficient, coordinate, upperBound, lowerBound>& simplex);

#endif
