/*
 * ppl_adaptors.h
 *
 *  Created on: Sep 11, 2009
 *      Author: frehse
 */

#ifndef PPL_ADAPTORS_H_
#define PPL_ADAPTORS_H_

#include "core/hybrid_automata/adapt_automaton_visitor.h"

namespace hybrid_automata {

class continuous_set_to_PPL_adaptor: public adapt_continuous_set_visitor {
public:
	void dispatch(const continuous::support_function_provider* c);
	void dispatch(const continuous::constr_polyhedron<global_types::rational_type>* c);
	void dispatch(const continuous::constr_polyhedron<global_types::float_type>* c);
	void dispatch(const ppl_polyhedron::continuous_set_PPL_NNC* c);
	void dispatch(const continuous::support_function::sfm_cont_set<global_types::float_type>* c);
	void dispatch(const continuous::spacetime_flowpipe<global_types::float_type>* c);
	void dispatch(const continuous::predicate_continuous_set* c);
	void dispatch(const continuous::continuous_set_simulation<global_types::float_type>* c);
};

class dynamics_to_PPL_adaptor: public adapt_dynamics_visitor {
public:
	void dispatch(const continuous::ode_affine_dynamics<global_types::rational_type>* c);
	void dispatch(const continuous::ode_affine_dynamics<global_types::float_type>* c);
};

class transform_to_PPL_adaptor: public adapt_transform_visitor {
public:
	transform_to_PPL_adaptor* clone();
};

class automaton_to_PPL_adaptor: public adapt_automaton_visitor {
public:
	automaton_to_PPL_adaptor();
	virtual ~automaton_to_PPL_adaptor();
	virtual automaton_to_PPL_adaptor* clone() const;
	virtual void init();
	virtual void epilogue(hybrid_automaton& h);
};

typedef boost::shared_ptr<automaton_to_PPL_adaptor> automaton_to_PPL_adaptor_ptr;

}

#endif /* PPL_ADAPTORS_H_ */
