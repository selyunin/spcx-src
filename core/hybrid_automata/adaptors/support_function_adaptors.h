/*
 * support_function_adaptors.h
 *
 *  Created on: Nov 9, 2009
 *      Author: frehse
 */

#ifndef SUPPORT_FUNCTION_ADAPTORS_H_
#define SUPPORT_FUNCTION_ADAPTORS_H_

#include "core/hybrid_automata/adapt_automaton_visitor.h"

namespace hybrid_automata {

template<typename scalar_type>
class continuous_set_to_supp_f_adaptor: public adapt_continuous_set_visitor {
public:
	void dispatch(const continuous::support_function_provider* c);
	void dispatch(const continuous::constr_polyhedron<Rational>* c);
	void dispatch(const continuous::constr_polyhedron<double>* c);
	void dispatch(const ppl_polyhedron::continuous_set_PPL_NNC* c);
	void dispatch(const continuous::support_function::sfm_cont_set<global_types::float_type>* c);
	void dispatch(const continuous::spacetime_flowpipe<global_types::float_type>* c);
	void dispatch(const continuous::predicate_continuous_set* c);
	void dispatch(const continuous::continuous_set_simulation<global_types::float_type>* c);
};

template<typename scalar_type>
class dynamics_to_supp_f_adaptor: public adapt_dynamics_visitor {
public:
	typedef boost::shared_ptr<dynamics_to_supp_f_adaptor> ptr;

	/** Constructor
	 *
	 * @parameter accept_nonlinear Activates or deactives parsing of nonlinear equations
	 *
	 * Default is true. */
	dynamics_to_supp_f_adaptor(bool accept_nonlinear) : my_accept_nonlinear(accept_nonlinear) {
//		std::cout << "accept nonlinear:" << std::boolalpha << my_accept_nonlinear << std::endl;
	};

	void dispatch(const continuous::ode_affine_dynamics<
			global_types::rational_type>* c);
	void dispatch(const continuous::ode_affine_dynamics<
			global_types::float_type>* c);
	void dispatch(const continuous::relation_dynamics* c);

private:
	bool my_accept_nonlinear;
};

template<typename scalar_type>
class transform_to_supp_f_adaptor: public adapt_transform_visitor {
public:
	transform_to_supp_f_adaptor* clone();
	virtual void dispatch(const continuous::relation_transform* c);
	virtual void dispatch(
			const continuous::constant_bound_time_elapse_transform* c);
	virtual void dispatch(const continuous::reset_affine_transform<
			global_types::rational_type>* c);
	virtual void dispatch(const continuous::reset_affine_transform<
			global_types::float_type>* c);
	virtual void dispatch(const continuous::reset_function_transform* c);
};

class automaton_to_supp_f_adaptor: public adapt_automaton_visitor {
public:
	automaton_to_supp_f_adaptor(coefficient_type bool_t,
			coefficient_type number_t, bool accept_nonlinear);
	virtual ~automaton_to_supp_f_adaptor();
	virtual void init();
	virtual void epilogue(hybrid_automaton& h);
private:
	bool my_accept_nonlinear;
};

typedef boost::shared_ptr<automaton_to_supp_f_adaptor>
		automaton_to_supp_f_adaptor_ptr;

}

#include "support_function_adaptors.hpp"

#endif /* SUPPORT_FUNCTION_ADAPTORS_H_ */
