#include <glpk.h>
#include <assert.h>
#include <string>
#include <vector>
#include "simplex.hh"

void print_status(glp_prob *lp);

class GLPKSimplex {
  glp_prob *lp;
  unsigned nrAllVariables;
  unsigned nrTrueVariables;

public:
  GLPKSimplex(const std::vector<std::string>& variableNames,
	      unsigned nrTrueVariables);

  ~GLPKSimplex() {
    glp_delete_prob(lp);
  }

  std::vector<variableType> variableTypes() const;

  void setBounds(unsigned i, double lb, double ub);

  void printStatistics(std::ostream& out) const;
  bool check();
  bool fastCheck();
  bool checkWellFormed() const {
    return true; // TODO
  }
  void assignBasic(unsigned i, const LinearForm<double>& form);
  void printStatus() const {
    print_status(const_cast<glp_prob*>(lp));
  }
};

inline std::ostream& operator<<(std::ostream& out, const GLPKSimplex& simplex) {
  return out;
} // TODO

