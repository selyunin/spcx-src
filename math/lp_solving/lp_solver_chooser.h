#ifndef LP_SOLVER_CHOOSER_H_
#define LP_SOLVER_CHOOSER_H_

//#include "lp_solver.h"
#include "math/lp_solving/lp_solver_fm.h"
#include "math/lp_solving/lp_solver_glpk.h"

namespace math {

template <typename scalar_type> class lp_solver;
template <typename scalar_type> class lp_solver_fm;

template <typename scalar_type> class lp_solver_chooser {
public:
	static lp_solver<scalar_type>* get_instance() {
//		static lp_solver_glpk<scalar_type> inst;
		lp_solver<scalar_type>* inst=new lp_solver_fm<scalar_type>();
		return inst;
	}
private:
	lp_solver_chooser(const lp_solver_chooser<scalar_type>&);
	lp_solver_chooser<scalar_type>& operator =(
			const lp_solver_chooser<scalar_type>& other);
};

/** Specialization: Use fm (exact) for Rational. */
template <> class lp_solver_chooser<Rational> {
public:
	static lp_solver<Rational>* get_instance() {
		lp_solver<Rational>* inst=new lp_solver_fm<Rational>();
		return inst;
	}
};

/** Specialization: Use glpk (not exact) for double. */
template <> class lp_solver_chooser<double> {
public:
	static lp_solver<double>* get_instance() {
		lp_solver<double>* inst=new lp_solver_glpk<double>();
		return inst;
	}
};

}

#endif /*LP_SOLVER_CHOOSER_H_*/
