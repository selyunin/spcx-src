#include "fastSimplex.hh"

class IntegerSimplex : MY_SIMPLEX {
  std::vector<bool> isIntegerVariable;

  inline bool checkIfBadIndex(unsigned index, integer& lowSplit);
  bool findBadIndex(unsigned& badIndex, integer& lowSplit);

public:
  const MY_SIMPLEX& getRationalSimplex() const {
    return *((MY_SIMPLEX*) this);
  }

  IntegerSimplex(const std::vector<std::string>& variableNames) :
    MY_SIMPLEX(variableNames), isIntegerVariable(variableNames.size()) {
  }

  void setIntegerFlag(unsigned i, bool flag) {
    isIntegerVariable[i] = flag;
  }

  bool getIntegerFlag(unsigned i) const {
    return isIntegerVariable[i];
  }

  void assignBasic(unsigned index, const LinearForm<coefficient>& form) {
    MY_SIMPLEX::assignBasic(index, form);
  }

  void checkWellFormed() const {
    MY_SIMPLEX::checkWellFormed();
  }

  void assertUpperBound(unsigned i, const upperBound& bound);
  void assertLowerBound(unsigned i, const lowerBound& bound);
  void relaxUpperBound(unsigned i, const upperBound& bound);
  void relaxLowerBound(unsigned i, const lowerBound& bound);

  const coordinate& getCurrent(unsigned i) const {
    return MY_SIMPLEX::getCurrent(i);
  }

  const upperBound& getUpperBound(unsigned i) const {
    return MY_SIMPLEX::getUpperBound(i);
  }

  const lowerBound& getLowerBound(unsigned i) const {
    return MY_SIMPLEX::getLowerBound(i);
  }

  void forcePivot(const std::vector<variableType>& coordinateTypes) {
    MY_SIMPLEX::forcePivot(coordinateTypes);
  }

  bool check();

  void printStatistics(std::ostream& out);
};

std::ostream& operator<<(std::ostream& out, const IntegerSimplex& s);
