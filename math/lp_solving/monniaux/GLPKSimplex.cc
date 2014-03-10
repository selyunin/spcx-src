#include "GLPKSimplex.hh"
#include <math.h>
#include <stdlib.h>

GLPKSimplex::GLPKSimplex(const std::vector<std::string>& variableNames,
	    unsigned nrTrueVariables0)
 : lp(glp_create_prob()),
   nrAllVariables(variableNames.size()),
   nrTrueVariables(nrTrueVariables0) {
  assert(nrTrueVariables <= nrAllVariables);
  glp_add_cols(lp, nrTrueVariables);
  glp_add_rows(lp, nrAllVariables-nrTrueVariables);

  for(unsigned i=0; i<nrTrueVariables; i++) {
    glp_set_col_name(lp, i+1, variableNames[i].c_str());
    glp_set_col_bnds(lp, i+1, GLP_FR, 0.0, 0.0);
  }

  for(unsigned i=nrTrueVariables; i<nrAllVariables; i++) {
    glp_set_row_name(lp, i-nrTrueVariables+1, variableNames[i].c_str());
  }
}

void GLPKSimplex::setBounds(unsigned i, double lb, double ub) {
  int kind;
  if (isinf(lb)) {
    if (isinf(ub)) {
      kind = GLP_FR;
    } else {
      kind = GLP_UP;
    }
  } else {
    if (isinf(ub)) {
      kind = GLP_LO;
    } else {
      if (lb == ub) {
	kind = GLP_FX;
      } else {
	kind = GLP_DB;
      }
    }
  }
  if (i < nrTrueVariables) {
    glp_set_col_bnds(lp, i+1, kind, lb, ub);
  } else {
    glp_set_row_bnds(lp, i-nrTrueVariables+1, kind, lb, ub);
  }
}

void GLPKSimplex::printStatistics(std::ostream& out) const {
  out << std::endl;
}

void print_status(glp_prob *lp) {
  for(unsigned i=1; i<=(unsigned)glp_get_num_cols(lp); i++) {
    std::cout << "col[" << i << "]: " << glp_get_col_type(lp, i) << std::endl;
  }
  for(unsigned i=1; i<=(unsigned)glp_get_num_rows(lp); i++) {
    std::cout << "row[" << i << "]: " << glp_get_row_type(lp, i) << std::endl;
  }
}

bool GLPKSimplex::check() {
  glp_smcp simplexParameters;
  glp_init_smcp(&simplexParameters);
  simplexParameters.msg_lev = GLP_MSG_OFF;
  simplexParameters.meth = GLP_DUALP;
  glp_simplex(lp, &simplexParameters);

  switch (glp_get_status(lp)) {
  case GLP_OPT:
    //std::cout << "GLPK: sat" << std::endl;
    return true;
  case GLP_NOFEAS:
    //std::cout << "GLPK: unsat" << std::endl;
    return false;
  default:
    //std::cout << "glp_get_status: " << glp_get_status(lp) << std::endl; 
    return true;
  }
}

bool GLPKSimplex::fastCheck() {
  glp_smcp simplexParameters;
  glp_init_smcp(&simplexParameters);
  simplexParameters.msg_lev = GLP_MSG_OFF;
  simplexParameters.meth = GLP_PRIMAL;
  glp_simplex(lp, &simplexParameters);

  switch (glp_get_status(lp)) {
  case GLP_OPT:
    //std::cout << "GLPK: sat" << std::endl;
    return true;
  case GLP_NOFEAS:
    //std::cout << "GLPK: unsat" << std::endl;
    return false;
  default:
    //std::cout << "glp_get_status: " << glp_get_status(lp) << std::endl; 
    return true;
  }
}

static variableType glp_stat_to_variableType(int x) {
  switch(x) {
  case GLP_BS: return BASIC;
  case GLP_NL: case GLP_NS: return LOWER;
  case GLP_NU: return UPPER;
  case GLP_NF: return ZERO;
  default: abort();
  }
}

std::vector<variableType> GLPKSimplex::variableTypes() const {
  std::vector<variableType> ret(nrAllVariables);

  for(unsigned i=0; i<nrTrueVariables; i++) {
    ret[i] = glp_stat_to_variableType(glp_get_col_stat(lp, i+1));
  }

  for(unsigned i=nrTrueVariables; i<nrAllVariables; i++) {
    ret[i] = glp_stat_to_variableType(glp_get_row_stat(lp, i-nrTrueVariables+1));
  }

  return ret;
}

void GLPKSimplex::assignBasic(unsigned var, const LinearForm<double>& form) {
  assert(var >= nrTrueVariables);

  unsigned maxNr = form.size();
  assert(maxNr > 0);

  int indices[maxNr+1];
  double coeffs[maxNr+1];
  unsigned k=1;
  for(LinearForm<double>::const_iterator xi = form.begin();
      xi != form.end(); xi++) {
    assert(k < maxNr);
    indices[k] = xi.index()+1;
    coeffs[k] = *xi;
    k++;
  }

  glp_set_mat_row(lp, var-nrTrueVariables+1, k-1, indices, coeffs);
}
