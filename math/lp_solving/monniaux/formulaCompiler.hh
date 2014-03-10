#ifndef _FORMULACOMPILER_HH
#define _FORMULACOMPILER_HH

#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include "constraints.hh"
#include "Solver.h"

typedef AffineForm<rational> CompiledExpression;

std::ostream& operator<<(std::ostream& out, Lit lit);

template <typename source, typename indexType> class AutoMapper {
  typedef std::map<source, indexType> mapper;
  mapper sourceToIndex;
  std::vector <source> indexToSource;

protected:
  virtual indexType nextIndex() = 0;

public:
  indexType get(const source& x, bool forcedNew=false) {
    if (forcedNew) {
      indexType index = nextIndex();
      indexToSource.push_back(x);
      return index;
    } else {
      typename mapper::iterator it = sourceToIndex.find(x);
      if (it == sourceToIndex.end()) {
	indexType index = nextIndex();
	sourceToIndex.insert(it, std::pair<source, indexType>(x, index));
	indexToSource.push_back(x);
	return index;
      } else {
	return it->second;
      }
    }
  }

  const source& revGet(indexType x) const {
    return indexToSource[x];
  }

  const std::vector <source>& rev() const {
    return indexToSource;
  }

  source& revGet(indexType x) {
    return indexToSource[x];
  }
  
  AutoMapper() {
  }

  void clear() {
    sourceToIndex.clear();
    indexToSource.clear();
  }

  size_t size() {
    return rev().size();
  }
};

template <typename source, typename indexType> class AutoIncMapper :
  public AutoMapper <source, indexType> {

  indexType index;

protected:
  indexType nextIndex() {
    return index++;
  }

public:
  AutoIncMapper() : index(0) {
  }

  void clear() {
    AutoMapper <source, indexType>::clear();
    index = 0;
  }
};


class Proposition {
public:
  enum Kind { PROPOSITIONAL, LINEAR };
private:
  Kind kind;
  union {
    struct {
      LinearConstraint *linearConstraint;
      unsigned simplexVariable;
    };

    struct {
      bool isBasicProposition;
      std::string *name;
    };
  };

  size_t hashV;

public:
  size_t hash() const {
    return hashV;
  }

  bool operator==(const Proposition& x) const {
    if (kind == PROPOSITIONAL) {
      return x.kind == PROPOSITIONAL && *name == *x.name;
    } else {
      return x.kind == LINEAR && *linearConstraint == *x.linearConstraint;
    }
  }

  void setSimplexVariable(unsigned variable) {
    simplexVariable = variable;
  }

  unsigned getSimplexVariable() const {
    return simplexVariable;
  }

  Kind getKind() const {
    return kind;
  }

  bool getIsBasicProposition() const {
    return isBasicProposition;
  }

  const LinearConstraint& getLinearConstraint() const {
    assert(kind == LINEAR);
    return *linearConstraint;
  }

  const std::string& getPropositionName() const {
    assert(kind == PROPOSITIONAL);
    return *name;
  }

private:
  static inline size_t djb2(size_t a, size_t b) {
    return (a << 5) + a + b;
  }

public:
  Proposition(const LinearConstraint& ct) :
    kind(LINEAR),
    linearConstraint(new LinearConstraint(ct)),
    hashV(djb2(LINEAR, ct.hash())) {
  }

  Proposition(const std::string& name0, bool isBasicProposition0) :
    kind(PROPOSITIONAL),
    isBasicProposition(isBasicProposition0),
    name(new std::string (name0)),
    hashV(djb2(PROPOSITIONAL, std::hash<std::string> () (name0))) {
  }

  Proposition(const Proposition& prop) {
    switch (prop.kind) {
    case LINEAR:
      kind = LINEAR;
      linearConstraint = new LinearConstraint(*prop.linearConstraint);
      break;
    case PROPOSITIONAL:
      kind = PROPOSITIONAL;
      name = new std::string(*prop.name);
      isBasicProposition = prop.isBasicProposition;
      break;
    default:
      assert(false);
    }
  }

  ~Proposition() {
    switch (kind) {
    case LINEAR:
      delete linearConstraint;
      break;
    case PROPOSITIONAL:
      delete name;
      break;
    }
  }

  bool isBasic() const {
    return kind != PROPOSITIONAL || isBasicProposition;
  }
};

namespace std {
  template <> struct less<Proposition> {
    bool operator() (const Proposition& x, const Proposition& y) const {
      if (x.getKind() < y.getKind()) return true;
      if (x.getKind() > y.getKind()) return false;
      switch (x.getKind()) {
      case Proposition::PROPOSITIONAL:
	return x.getPropositionName() < y.getPropositionName();
      case Proposition::LINEAR:
	return std::less<LinearConstraint>() (x.getLinearConstraint(), y.getLinearConstraint());
      }
      abort();
    }
  };

  template<> struct hash<Proposition> {
    size_t operator()(const Proposition& a) const {
      return a.hash();
    }
  };
}

class PropositionMapper : public AutoMapper <Proposition, Var> {
  Solver *solver;

protected:
  Var nextIndex() {
    return solver->newVar();
  }

public:
  Var get(const Proposition& x, bool forcedNew=false) {
    Var var = AutoMapper <Proposition, Var>::get(x, forcedNew);
    solver->setDecisionVar(var, x.getKind()==Proposition::LINEAR
			  || x.getIsBasicProposition());
    return var;
  }

  PropositionMapper(Solver *solver0) : solver(solver0) {
  }
};

std::ostream& operator<<(std::ostream& out, const Proposition& prop);

struct Environment {
  typedef std::map <std::string, Lit> propVarMap_t;
  propVarMap_t propVarMap;
  typedef std::map <std::string, CompiledExpression> numVarMap_t;
  numVarMap_t numVarMap;
};

class FormulaCompiler {
  Solver solver;
  AutoIncMapper<std::string, unsigned> numVars;
  PropositionMapper propVars;
  Lit trueVarLit;

  typedef LinearForm<integer> Form;
  AutoIncMapper<Form, unsigned> forms;

private:
  void collectForms();

public:
  FormulaCompiler();

  const std::vector<Form>& getForms() {
    return forms.rev();
  }

  const Form & getForm(unsigned i) {
    return forms.revGet(i);
  }

  bool solve() {
    vec<Lit> assumptions;
    return solver.solve(assumptions);
  }

  bool solve(vec<Lit>& assumptions);

  lbool modelValue(Lit x) const {
    return solver.modelValue(x);
  }

  lbool modelValue(Var x) const {
    return solver.modelValue(Lit(x));
  }

  Lit trueVar() const {
    return trueVarLit;
  }

  Lit getPropLit(const Proposition& prop);

  unsigned getNumVar(const std::string& name) {
    return numVars.get(name);
  }

  const std::vector<Proposition>& propMeanings() const {
    return propVars.rev();
  }

  const std::vector<std::string>& numVarNames() const {
    return numVars.rev();
  }

  void addClause(vec<Lit>& clause);
  void addBinaryClause(Lit x, Lit y);
  void genAnd(Lit target, std::vector<Lit>& operands);
  void genAnd(Lit target, Lit x, Lit y) {
    std::vector<Lit> lits(2);
    lits[0]=x;
    lits[1]=y;
    genAnd(target, lits);
  }
  void genIff(Lit c, Lit a, Lit b);
  void genIfThenElse(Lit out, Lit condition, Lit e1, Lit e2);
  void genTrivialTheoryPropagation();
};

std::ostream& operator<<(std::ostream& out, lbool v);
#endif
