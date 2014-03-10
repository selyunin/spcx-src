/*
 * constr_poly_adaptors.h
 *
 *  Created on: Sep 11, 2009
 *      Author: frehse
 */

#ifndef constr_poly_ADAPTORS_H_
#define constr_poly_ADAPTORS_H_

#include "core/hybrid_automata/adapt_automaton_visitor.h"

namespace hybrid_automata {

class continuous_set_to_constr_poly_adaptor: public adapt_continuous_set_visitor {
public:
	typedef adapt_automaton_visitor::coefficient_type coefficient_type;
	continuous_set_to_constr_poly_adaptor(coefficient_type bool_t,coefficient_type number_t);
	void dispatch(const continuous::support_function_provider* c);
	void dispatch(const continuous::constr_polyhedron<Rational>* c);
	void dispatch(const continuous::constr_polyhedron<double>* c);
	void dispatch(const ppl_polyhedron::continuous_set_PPL_NNC* c);
	void dispatch(const continuous::support_function::sfm_cont_set<global_types::float_type>* c);
	void dispatch(const continuous::spacetime_flowpipe<global_types::float_type>* c);
	void dispatch(const continuous::predicate_continuous_set* c);
	void dispatch(const continuous::continuous_set_simulation<global_types::float_type>* c);
private:
	coefficient_type my_bool_type;
	coefficient_type my_number_type;
};

class dynamics_to_constr_poly_adaptor: public adapt_dynamics_visitor {
public:
	void dispatch(const continuous::ode_affine_dynamics<global_types::rational_type>* c);
	void dispatch(const continuous::ode_affine_dynamics<global_types::float_type>* c);
};

class transform_to_constr_poly_adaptor: public adapt_transform_visitor {
public:
	transform_to_constr_poly_adaptor* clone();
};

class automaton_to_constr_poly_adaptor: public adapt_automaton_visitor {
public:
	automaton_to_constr_poly_adaptor(coefficient_type bool_t,coefficient_type number_t);
	virtual ~automaton_to_constr_poly_adaptor();
	virtual automaton_to_constr_poly_adaptor* clone() const;
	virtual void init();
	virtual void epilogue(hybrid_automaton& h);
};

typedef boost::shared_ptr<automaton_to_constr_poly_adaptor> automaton_to_constr_poly_adaptor_ptr;

}

#endif /* constr_poly_ADAPTORS_H_ */
