#include "formulaCompiler.hh"
#include "smtlib_syntax.hh"
#include <stdlib.h>
#include <stack>
#include <sstream>
#include "fastSimplex.hh"
#include "timings.h"

//#define PROPNAMES
//#define PRINT_CLAUSE

std::ostream& operator<<(std::ostream& out, const Proposition& prop) {
  switch (prop.getKind()) {
  case Proposition::LINEAR:
    out << prop.getLinearConstraint();
    break;
  case Proposition::PROPOSITIONAL:
    out << prop.getPropositionName();
    break;
  }
  return out;
}

Lit FormulaCompiler::getPropLit(const Proposition& prop) {
  if (prop.getKind() == Proposition::LINEAR &&
      prop.getLinearConstraint().linearPart().isZero()) {
    bool x;
    if (prop.getLinearConstraint().getComparator()==LESS_OR_EQUAL){
      x = rational(0) <= prop.getLinearConstraint().rhsPart();
    } else {
      x = rational(0) >= prop.getLinearConstraint().rhsPart();
    }
    return x ? trueVarLit : ~trueVarLit;
  }
  return Lit(propVars.get(prop, !prop.isBasicProposition));
}

Lit FLet::compile(FormulaCompiler& compiler, Environment& env) const {
  Lit bindingV = binding->compile(compiler, env);
  Environment::propVarMap_t::iterator it= env.propVarMap.find(boundVariable);
  if (it == env.propVarMap.end()) {
    it = env.propVarMap.insert(it,
      std::pair<std::string, Lit> (boundVariable, bindingV));
    Lit ret = subFormula->compile(compiler, env);
    //std::cout << "flet=" << ret << std::endl;
    env.propVarMap.erase(it);
    return ret;
  } else {
    Lit savedValue = it->second;
    it->second=bindingV;
    Lit ret = subFormula->compile(compiler, env);
    it->second=savedValue;
    return ret;
  }
}

Lit And::compile(FormulaCompiler& compiler, Environment& env) const {
  if (begin() == end()) return compiler.trueVar();

#ifdef PROPNAMES
  std::ostringstream name;
  name << "(and";
#endif
  std::vector<Lit> operands;

  for(const_iterator it=begin(); it!=end(); it++) {
    Lit lit = (*it)->compile(compiler, env);
    if (lit==~compiler.trueVar()) return ~compiler.trueVar();
    operands.push_back(lit);
#ifdef PROPNAMES
    name << " " << lit;
#endif
  }
#ifdef PROPNAMES
  name << ")";
#endif

  Lit target(compiler.getPropLit(Proposition(
#ifdef PROPNAMES
    name.str()
#else
    ""
#endif
    , false)));
  compiler.genAnd(target, operands);
  //std::cout << "and=" << target << std::endl;
  return target;
}

Lit Or::compile(FormulaCompiler& compiler, Environment& env) const {
  if (begin() == end()) return ~compiler.trueVar();

#ifdef PROPNAMES
  std::ostringstream name;
  name << "(or";
#endif
  std::vector<Lit> operands;

  for(const_iterator it=begin(); it!=end(); it++) {
    Lit lit = (*it)->compile(compiler, env);
    if (lit==compiler.trueVar()) return compiler.trueVar();
    operands.push_back(~lit);
#ifdef PROPNAMES
    name << " " << lit;
#endif
  }
#ifdef PROPNAMES
  name << ")";
#endif

  Lit target(compiler.getPropLit(Proposition(
#ifdef PROPNAMES
    name.str()
#else
    ""
#endif
    , false)));
  compiler.genAnd(~target, operands);
  return target;
}

static Lit binaryIffCompile(FormulaCompiler& compiler,
			    Lit x, Lit y) {
#ifdef PROPNAMES
  std::ostringstream name;
  name << "(iff " << x << " " << y << ")";
#endif
  Lit target(compiler.getPropLit(Proposition(
#ifdef PROPNAMES
    name.str()
#else
    ""
#endif
    , false)));
  compiler.genIff(target, x, y);
  return target;
}

Lit Iff::compile(FormulaCompiler& compiler, Environment& env) const {
  const_iterator first=begin();
  if (first==end()) return compiler.trueVar();
  const_iterator second=first;
  second++;
  if (second==end()) return compiler.trueVar();
  const_iterator third=second;
  third++;
  if (third==end()) {
    return binaryIffCompile(compiler,
			    (*first)->compile(compiler, env),
			    (*second)->compile(compiler, env));
  } else {
    std::vector<Lit> outputs;
    Lit firstLit = (*first)->compile(compiler, env);
#ifdef PROPNAMES
    std::ostringstream name;
    name << "(iff " << firstLit;
#endif
    while(second != end()) {
      Lit secondLit = (*second)->compile(compiler, env);
#ifdef PROPNAMES
      name << " " << secondLit;
#endif
      outputs.push_back(binaryIffCompile(compiler, firstLit, secondLit));
      first=second;
      firstLit=secondLit;
      second++;
    }
#ifdef PROPNAMES
    name << ")";
#endif
    Lit target(compiler.getPropLit(Proposition(
#ifdef PROPNAMES
      name.str()
#else
      ""
#endif
      , false)));
    compiler.genAnd(target, outputs);
    return target;
  }
}

Lit PropositionalVariable::compile(FormulaCompiler& compiler, Environment& env) const { // TODO
  Environment::propVarMap_t::const_iterator it = env.propVarMap.find(name);
  if (it != env.propVarMap.end()) {
    return it->second;
  } else {
    return compiler.getPropLit(Proposition(name, true));
  }
}

Lit Not::compile(FormulaCompiler& compiler, Environment& env) const {
  return ~f->compile(compiler, env);
}

static Lit binaryEqualityCompile(FormulaCompiler& compiler,
				 const CompiledExpression& x,
				 const CompiledExpression& y) {
  std::vector<Lit> outputs;
  outputs.push_back(compiler.getPropLit(Proposition(LinearConstraint::genComparison(x, LESS_OR_EQUAL, y))));
  outputs.push_back(compiler.getPropLit(Proposition(LinearConstraint::genComparison(x, GREATER_OR_EQUAL, y))));

#ifdef PROPNAMES
  std::ostringstream name;
  name << "(= " << x << " " << y << ")";
#endif
  Lit target(compiler.getPropLit(Proposition(
#ifdef PROPNAMES
    name.str()
#else
    ""
#endif
    , false)));
  compiler.genAnd(target, outputs);
  return target;
}

Lit Equality::compile(FormulaCompiler& compiler, Environment& env) const {
  const_iterator first=begin();
  if (first==end()) return compiler.trueVar();
  const_iterator second=first;
  second++;
  if (second==end()) return compiler.trueVar();

  std::vector<Lit> outputs;
  CompiledExpression firstExpr = (*first)->compile(compiler, env);
  std::ostringstream name;
#ifdef PROPNAMES
  name << "(= " << firstExpr;
#endif

  while(second != end()) {
    CompiledExpression secondExpr = (*second)->compile(compiler, env);
#ifdef PROPNAMES
    name << " " << secondExpr;
#endif
    outputs.push_back(compiler.getPropLit(Proposition(LinearConstraint::genComparison(firstExpr, LESS_OR_EQUAL, secondExpr))));
    outputs.push_back(compiler.getPropLit(Proposition(LinearConstraint::genComparison(firstExpr, GREATER_OR_EQUAL, secondExpr))));
    
    first=second;
    firstExpr=secondExpr;
    second++;
  }
#ifdef PROPNAMES
  name << ")";
#endif
  Lit target(compiler.getPropLit(Proposition(
#ifdef PROPNAMES
    name.str()
#else
    ""
#endif
    , false)));
  compiler.genAnd(target, outputs);
  return target;
}

Lit PairwiseDistinct::compile(FormulaCompiler& compiler, Environment& env) const {
  std::vector<CompiledExpression> compiledExprs;
  for(const_iterator it=begin(); it!=end(); it++) {
    compiledExprs.push_back((*it)->compile(compiler, env));
  }

  std::vector<CompiledExpression>::const_iterator first=compiledExprs.begin();
  if (first==compiledExprs.end()) return compiler.trueVar();
  {
    std::vector<CompiledExpression>::const_iterator second=first;
    second++;
    if (second==compiledExprs.end()) return compiler.trueVar();
    std::vector<CompiledExpression>::const_iterator third=second;
    third++;
    if (third==compiledExprs.end()) {
      return ~binaryEqualityCompile(compiler, *first, *second);
    }
  }

  std::vector<Lit> outputs;
#ifdef PROPNAMES
  std::ostringstream name;
  name << "(distinct";
#endif
  while (true) {
    std::vector<CompiledExpression>::const_iterator second=first;
    second++;
    if (second==compiledExprs.end()) break;
    for(; second!=compiledExprs.end(); second++) {
      Lit l = ~binaryEqualityCompile(compiler, *first, *second);
      outputs.push_back(l);
#ifdef PROPNAMES
      name << " " << l;
#endif
    }
    first++;
  }

#ifdef PROPNAMES
  name << ")";
#endif
  Lit target(compiler.getPropLit(Proposition(
#ifdef PROPNAMES
    name.str()
#else
    ""
#endif
    , false)));
  compiler.genAnd(target, outputs);
  return target;
}

Lit Comparison::compile(FormulaCompiler& compiler, Environment& env) const {
  CompiledExpression e1C = e1->compile(compiler, env),
    e2C = e2->compile(compiler, env);
  comparator op;
  switch (comparisonKind) {
  case LESS_EQ: case GREATER:
    op = LESS_OR_EQUAL;
    break;
  case GREATER_EQ: case LESS:
    op = GREATER_OR_EQUAL;
    break;
  }
  Lit x(compiler.getPropLit(Proposition(LinearConstraint::genComparison(e1C, op, e2C))));
  return (comparisonKind==LESS || comparisonKind==GREATER) ? ~x : x;
}

Lit ExprLet::compile(FormulaCompiler& compiler, Environment& env) const {
  CompiledExpression bindingV = binding->compile(compiler, env);
  Environment::numVarMap_t::iterator it=env.numVarMap.find(boundVariable);
  if (it == env.numVarMap.end()) {
    it = env.numVarMap.insert(it,
      std::pair<std::string, CompiledExpression> (boundVariable, bindingV));
    Lit ret = subFormula->compile(compiler, env);
    env.numVarMap.erase(it);
    return ret;
  } else {
    CompiledExpression savedValue = it->second;
    it->second=bindingV;
    Lit ret = subFormula->compile(compiler, env);
    it->second=savedValue;
    return ret;
  }
}

Lit FormulaIfThenElse::compile(FormulaCompiler& compiler, Environment& env) const {
  Lit conditionV = condition->compile(compiler, env);
  if (conditionV == compiler.trueVar())
    return e1->compile(compiler, env);
  if (conditionV == ~compiler.trueVar())
    return e2->compile(compiler, env);

  Lit e1V = e1->compile(compiler, env),
    e2V = e2->compile(compiler, env);
#ifdef PROPNAMES
  std::ostringstream name;
  name << "(if_then_else " << conditionV << " " << e1V << " " << e2V << ")";
#endif
  Lit target(compiler.getPropLit(Proposition(
#ifdef PROPNAMES
    name.str()
#else
    ""
#endif
    , false)));
  compiler.genIfThenElse(target, conditionV, e1V, e2V);
  return target;
}

CompiledExpression ExpressionIfThenElse::compile(FormulaCompiler& compiler, Environment& env) const {
  Lit conditionV = condition->compile(compiler, env);
  CompiledExpression e1V = e1->compile(compiler, env),
    e2V = e2->compile(compiler, env);
  CompiledExpression boundVarExpr(compiler.numVarNames());
  boundVarExpr.addMonomial(compiler.getNumVar(variableName()), 1);
  
  Lit cond1V, cond2V;
  {
    std::vector<Lit> cond1;
    cond1.push_back(conditionV);
    cond1.push_back(compiler.getPropLit(Proposition(LinearConstraint::genComparison(boundVarExpr, LESS_OR_EQUAL, e1V))));
    cond1.push_back(compiler.getPropLit(Proposition(LinearConstraint::genComparison(boundVarExpr, GREATER_OR_EQUAL, e1V))));

#ifdef PROPNAMES
    std::ostringstream propName;
    propName << "(and " << conditionV << " (= " << boundVarExpr << " " << e1V << "))";
#endif
    cond1V=compiler.getPropLit(Proposition(
#ifdef PROPNAMES
    propName.str()
#else
    ""
#endif
    , false));
    compiler.genAnd(cond1V, cond1);
  }

  {
    std::vector<Lit> cond2;
    cond2.push_back(~conditionV);
    cond2.push_back(compiler.getPropLit(Proposition(LinearConstraint::genComparison(boundVarExpr, LESS_OR_EQUAL, e2V))));
    cond2.push_back(compiler.getPropLit(Proposition(LinearConstraint::genComparison(boundVarExpr, GREATER_OR_EQUAL, e2V))));

#ifdef PROPNAMES
    std::ostringstream propName;
    propName << "(and " << (~conditionV) << " (= " << boundVarExpr << " " << e2V << "))";
#endif
    cond2V=compiler.getPropLit(Proposition(
#ifdef PROPNAMES
      propName.str()
#else
      ""
#endif
      , false));
    compiler.genAnd(cond2V, cond2);
  }

  compiler.addBinaryClause(cond1V, cond2V);
  return boundVarExpr;
}

CompiledExpression Quotient::compile(FormulaCompiler& compiler, Environment& env) const {
  return e1->compile(compiler, env) / e2->compile(compiler, env);
}

CompiledExpression Product::compile(FormulaCompiler& compiler, Environment& env) const {
  CompiledExpression ret(compiler.numVarNames(), 1);
  for(const_iterator it=begin(); it!=end(); it++) {
    ret *= (*it)->compile(compiler, env);
  }
  return ret;
}

CompiledExpression Sum::compile(FormulaCompiler& compiler, Environment& env) const {
  CompiledExpression ret(compiler.numVarNames());
  for(const_iterator it=begin(); it!=end(); it++) {
    ret += (*it)->compile(compiler, env);
  }
  return ret;
}

CompiledExpression IntegerConstant::compile(FormulaCompiler& compiler, Environment&) const {
  return CompiledExpression(compiler.numVarNames(), value);
}

CompiledExpression NumericVariable::compile(FormulaCompiler& compiler, Environment& env) const {
  Environment::numVarMap_t::const_iterator it = env.numVarMap.find(name);
  if (it != env.numVarMap.end()) {
    return it->second;
  } else {
    CompiledExpression ret(compiler.numVarNames());
    ret.addMonomial(compiler.getNumVar(name), 1);
    return ret;
  }
}

CompiledExpression UnaryNeg::compile(FormulaCompiler& compiler, Environment& env) const {
  return -(f->compile(compiler, env));
}

std::ostream& operator<<(std::ostream& out, Lit lit) {
  if (sign(lit)) {
    out << "!#" << var(lit);
  } else {
    out << "#" << var(lit);
  }
  return out;
}

void print(const vec<Lit> & clause) {
  std::cout << "clause =";
  for(unsigned i=0; i<(unsigned) clause.size(); i++) {
    std::cout << " " << clause[i];
  }
  std::cout << std::endl;  
}

void FormulaCompiler::addClause(vec<Lit>& clause) {
#ifdef PRINT_CLAUSE
  print(clause);
#endif
  solver.addClause(clause);
}

void FormulaCompiler::addBinaryClause(Lit x, Lit y) {
  vec<Lit> clause;
  clause.push(x);
  clause.push(y);
  addClause(clause);
}

void FormulaCompiler::genAnd(Lit target, std::vector<Lit>& operands) {
  for(unsigned i=0; i<operands.size(); i++) {
    addBinaryClause(~target, operands[i]);
  }
  vec<Lit> clause;
  clause.push(target);
  for(unsigned i=0; i<operands.size(); i++) {
    clause.push(~ operands[i]);
  }
  addClause(clause);
}

void FormulaCompiler::genIff(Lit c, Lit a, Lit b) {
  vec<Lit> clause;
  clause.push(~a); clause.push(~b); clause.push(c); 
  addClause(clause); clause.clear();

  clause.push(~a); clause.push(~c); clause.push(b);
  addClause(clause); clause.clear();

  clause.push(~b); clause.push(~c); clause.push(a);
  addClause(clause); clause.clear();

  clause.push(a); clause.push(b); clause.push(c);
  addClause(clause);
}

void FormulaCompiler::genIfThenElse(Lit out, Lit condition, Lit e1, Lit e2) {
  vec<Lit> clause;
  clause.push(~condition); clause.push(~e1); clause.push(~out);
  addClause(clause); clause.clear();

  clause.push(~condition); clause.push(e1); clause.push(out);
  addClause(clause); clause.clear();

  clause.push(condition); clause.push(~e2); clause.push(~out);
  addClause(clause); clause.clear();

  clause.push(condition); clause.push(e2); clause.push(out);
  addClause(clause);
}

FormulaCompiler::FormulaCompiler() :
  propVars(&solver),
  trueVarLit(getPropLit(Proposition("$true", false))) {
  vec<Lit> assumption;
  assumption.push(trueVarLit);
  addClause(assumption);
}

std::ostream& operator<<(std::ostream& out, lbool v) {
  if(v == l_True) {
    out << "true";
  } else if (v == l_False) {
    out << "false";
  } else {
    out << "undef";
  }
  return out;
}

void FormulaCompiler::collectForms() {
  forms.clear();
  for(unsigned i=0; i<numVars.size(); i++) {
    LinearForm<integer> form(numVars.rev());
    form.addMonomial(i, 1);
    forms.get(form);
  }

  for(unsigned i=0; i<propVars.size(); i++) {
    Proposition& prop = propVars.revGet(i);
    if (prop.getKind() != Proposition::LINEAR) continue;
    prop.setSimplexVariable(forms.get(prop.getLinearConstraint().linearPart()));
  }
}

typedef std::map<rational, Var> BoundMapper;
struct ForADirection {
  BoundMapper lowerBounds, upperBounds;
};

//typedef std::unordered_map<HashedLinearForm<integer>, ForADirection> DirectionMapper;
typedef std::map<LinearForm<integer>, ForADirection> DirectionMapper;
 
void FormulaCompiler::genTrivialTheoryPropagation() {
  DirectionMapper directionMapper;
  for(unsigned i=0; i<propVars.size(); i++) {
    const Proposition& prop = propVars.revGet(i);
    if (prop.getKind() != Proposition::LINEAR) continue;

    DirectionMapper::iterator it=directionMapper.find(prop.getLinearConstraint().linearPart());
    if (it == directionMapper.end()) {
      it = directionMapper.insert(std::pair<LinearForm<integer>, ForADirection> (prop.getLinearConstraint().linearPart(), ForADirection())).first;
    }
    BoundMapper& mapper = prop.getLinearConstraint().getComparator() ==
      LESS_OR_EQUAL ? it->second.upperBounds : it->second.lowerBounds;
    bool newInsertion = mapper.insert(std::pair<rational, Var>(prop.getLinearConstraint().rhsPart(), i)).second;
    //TODO look assert(newInsertion);
  }  

  for(DirectionMapper::const_iterator it=directionMapper.begin();
      it != directionMapper.end(); it++) {
#ifdef PRINT_DIRECTION_BOUNDS
    std::cout << it->first << " upper:";
    for(BoundMapper::const_iterator jt=it->second.upperBounds.begin();
	jt != it->second.upperBounds.end(); jt++) {
      std::cout << " " << jt->first;
    }
    std::cout << " lower:";
    for(BoundMapper::const_iterator jt=it->second.lowerBounds.begin();
	jt != it->second.lowerBounds.end(); jt++) {
      std::cout << " " << jt->first;
    }
    std::cout << std::endl;
#endif

    // <= is transitive
    {
      BoundMapper::const_iterator jt=it->second.upperBounds.begin(),
	jte = it->second.upperBounds.end();
      if (jt != jte) {
	BoundMapper::const_iterator jt2=jt;
	jt2++;
	while (jt2 != jte) {
	  addBinaryClause(Lit(jt->second, true), Lit(jt2->second));
	  jt=jt2;
	  jt2++;
	}
      }
    }

    // >= is transitive
    {
      BoundMapper::const_iterator jt=it->second.lowerBounds.begin(),
	jte = it->second.lowerBounds.end();
      if (jt != jte) {
	BoundMapper::const_iterator jt2=jt;
	jt2++;
	while (jt2 != jte) {
	  addBinaryClause(Lit(jt->second), Lit(jt2->second, true));
	  jt=jt2;
	  jt2++;
	}
      }
    }

    // no solution to: x <= a and x >= b when a < b
    for(BoundMapper::const_iterator jt=it->second.upperBounds.begin();
	jt != it->second.upperBounds.end(); jt++) {
      BoundMapper::const_iterator kt=
	it->second.lowerBounds.upper_bound(jt->first);
      if (kt == it->second.lowerBounds.end()) break;
      addBinaryClause(Lit(jt->second, true), Lit(kt->second, true));
    }

    // for all x, x <= a or x >= b when b <= a
    for(BoundMapper::const_iterator jt=it->second.lowerBounds.begin();
	jt != it->second.lowerBounds.end(); jt++) {
      BoundMapper::const_iterator kt=
	it->second.upperBounds.lower_bound(jt->first);
      if (kt == it->second.upperBounds.end()) break;
      addBinaryClause(Lit(jt->second), Lit(kt->second));
    }
  }
}


class lboolAssignment {
  std::vector<lbool> contents;

public:
  lbool& operator[](size_t i) {
    return contents[i];
  }

  const lbool& operator[](size_t i) const {
    return contents[i];
  }

  lboolAssignment(size_t size) : contents(size, l_Undef) {
  }

  void assign(Lit x) {
    Var v = var(x);
    if (toInt(x) < 0) {
      std::cerr << "bad result variable" << std::endl;
      abort();
    }
    lbool oldVal = (*this)[v], newVal = sign(x) ? l_False : l_True;
    if (oldVal == newVal) return;
    if (oldVal != l_Undef) {
      std::cerr << "wrong assignment" << std::endl;
      abort();
    }
    (*this)[v] = newVal;
  }
};

struct LinearTheory : public Theory {
  FormulaCompiler *compiler;

  MY_SIMPLEX simplex;
  
  std::vector<Lit> lowerBoundProps, upperBoundProps;

public:
  LinearTheory(FormulaCompiler *compiler0,
	 std::vector<std::string>& simplexVarNames,
	 unsigned numVars,
	 bool immediateTest);

  void reset();
  void assertLiteral(Lit lit);

  void printStatistics(std::ostream& out) const {
    simplex.printStatistics(out);
  }

  bool check() {
    return simplex.check();
  }

  bool fastCheck() {
    return simplex.fastCheck();
  }

  bool isSatisfied() const {
    return simplex.isSatisfied();
  }

  void fillContradictionClause(vec<Lit>& clause) const;
  void fillContradictionAssignment(lboolAssignment& assignment) const;
};

void LinearTheory::reset() {
  simplex.resetBounds();

  unsigned nrPropVars = compiler->propMeanings().size();
  for(unsigned i=0; i<lowerBoundProps.size(); i++) {
    lowerBoundProps[i] = upperBoundProps[i] = toLit(-1);
  }
}

LinearTheory::LinearTheory(FormulaCompiler *compiler0,
	       std::vector<std::string>& simplexVarNames,
	       unsigned numVarsNr,
	       bool immediateTest) :
  compiler(compiler0),
  simplex(simplexVarNames, numVarsNr, immediateTest),
  lowerBoundProps(simplexVarNames.size()),
  upperBoundProps(simplexVarNames.size()) {
  typedef MY_SIMPLEX::Row Row;
  
  for(unsigned i=numVarsNr; i<compiler0->getForms().size(); i++) {
    const Row& row = compiler0->getForm(i).getForm();
    simplex.assignBasic(i, row);
  }
}

void LinearTheory::assertLiteral(Lit lit) {
  unsigned i = var(lit);
  const Proposition& prop = compiler->propMeanings()[i];
  if (prop.getKind() != Proposition::LINEAR) return;
  if (!sign(lit)) {
    if (prop.getLinearConstraint().getComparator() == LESS_OR_EQUAL) {
      if (simplex.assertUpperBound(prop.getSimplexVariable(), coordinate(prop.getLinearConstraint().rhsPart()))) {
	upperBoundProps[prop.getSimplexVariable()] = Lit(i);
      }
    } else {
      if (simplex.assertLowerBound(prop.getSimplexVariable(), coordinate(prop.getLinearConstraint().rhsPart()))) {
	lowerBoundProps[prop.getSimplexVariable()] = Lit(i);
      }
    }
  } else {
    if (prop.getLinearConstraint().getComparator() == LESS_OR_EQUAL) {
      if (simplex.assertLowerBound(prop.getSimplexVariable(), coordinate(prop.getLinearConstraint().rhsPart(), 1))){
	lowerBoundProps[prop.getSimplexVariable()] = Lit(i, true);
      }
    } else {
      if (simplex.assertUpperBound(prop.getSimplexVariable(), coordinate(prop.getLinearConstraint().rhsPart(), -1))) {
	upperBoundProps[prop.getSimplexVariable()] = Lit(i, true);
      }
    }
  }
}

void LinearTheory::fillContradictionAssignment(lboolAssignment& assignment) const {
  if (simplex.getContradictionKind() == BETWEEN_BOUNDS_CONTRADICTION) {
    assignment.assign(lowerBoundProps[simplex.getContradictionRow()]);
    assignment.assign(upperBoundProps[simplex.getContradictionRow()]);
  } else {
    {
      Lit x = (simplex.getContradictionKind() == UPPER_BOUND_CONTRADICTION)
	? upperBoundProps[simplex.getContradictionRow()]
	: lowerBoundProps[simplex.getContradictionRow()];
      assignment.assign(x);
    }
    
    {
      const LinearForm<rational>& contradiction = simplex.getForm(simplex.getContradictionRow());
      
      for(LinearForm<rational>::const_iterator it=contradiction.begin();
	  it != contradiction.end(); it++) {
	if (*it == 0) continue;
	
	Lit x;
	if (((*it > 0) ==
	     (simplex.getContradictionKind() == UPPER_BOUND_CONTRADICTION))) {
	  x= lowerBoundProps[it.index()];
	} else {
	  x= upperBoundProps[it.index()];
	}
	assignment.assign(x);
      }
    }
    
    {
      const LinearForm<rational>& trk = simplex.getTrackerRow(simplex.getContradictionRow());
      
      for(LinearForm<rational>::const_iterator it=trk.begin();
	  it != trk.end(); it++) {
	const unsigned i=it.index();
	const coefficient& k=*it;
	if (k == 0) continue;
	
	Lit x;
	if (((k > 0) ==
	     (simplex.getContradictionKind()== UPPER_BOUND_CONTRADICTION))){ 
	  x= upperBoundProps[i];
	} else {
	  x= lowerBoundProps[i];
	}
	assignment.assign(x);
      }
    }
   }
}

void LinearTheory::fillContradictionClause(vec<Lit>& clause) const {
  clause.clear();

  unsigned nrPropVars = compiler->propMeanings().size();
  lboolAssignment assignment(nrPropVars);
  fillContradictionAssignment(assignment);

  for(unsigned i=0; i<nrPropVars; i++) {
    if (assignment[i] != l_Undef) {
      clause.push(assignment[i] == l_False ? Lit(i) : Lit(i, true));
    }
  }
}

bool FormulaCompiler::solve(vec<Lit>& assumptions) {
#ifdef TIMINGS
  double gen_theory_time=0., sat_time=0., simplex_setup_time=0., simplex_solve_time=0., simplex_init_time=0.;
#endif

  TIMING_DECLARE();
  TIMING_BEGINS();
  collectForms();
  genTrivialTheoryPropagation();
  TIMING_ENDS(gen_theory_time);

#if 0
  std::cout << "number of forms = " << forms.size() << std::endl;
  for(unsigned i=0; i<forms.size(); i++) {
    std::cout << getForm(i).getForm() << std::endl;
  }
#endif

  std::vector<std::string> simplexVarNames(numVars.rev());

  for(unsigned i=numVars.size(); i<forms.size(); i++) {
#ifdef PROPNAMES
    std::ostringstream os;
    os << "form#" << i;
    simplexVarNames.push_back(os.str());
#else
    simplexVarNames.push_back("");
#endif
  }

  TIMING_BEGINS();
  LinearTheory theory(this, simplexVarNames, numVars.size(), true);
  solver.setTheory(theory);
  TIMING_ENDS(simplex_init_time);

  bool solution=true;
  while (true) {
    std::cout << "SMT loop" << std::endl;

    theory.reset();

    TIMING_BEGINS();
    if (!solver.solve(assumptions)) {
      solution = false;
    }
    TIMING_ENDS(sat_time);
    if (!solution) break;

    std::cout << "starting simplex" << std::endl;

    TIMING_BEGINS();
    for(unsigned i=0; i<propVars.size() && theory.isSatisfied(); i++) {
      switch (solver.modelValue(Lit(i)).toInt()) {
      case 1:
	theory.assertLiteral(Lit(i));
	break;
      case -1:
	theory.assertLiteral(Lit(i, true));
	break;
      }
    }
    TIMING_ENDS(simplex_setup_time);
    
    TIMING_BEGINS();
    bool ret = theory.check();
    TIMING_ENDS(simplex_solve_time);
    
    if (ret) {
      solution = true;
      break;
    }

    vec<Lit> clause;
    theory.fillContradictionClause(clause);
    std::cout << "main feedback clause size: " << clause.size() << std::endl;
    addClause(clause);
  }

#ifdef TIMINGS
  theory.printStatistics(std::cout);

  std::cout << "gen_theory: " << gen_theory_time
	    << " sat: " << sat_time
	    << " simplex_setup: " << simplex_setup_time
	    << " simplex_init: " << simplex_init_time
	    << " simplex_solve: " << simplex_solve_time << std::endl;
#endif
  return solution;
}
