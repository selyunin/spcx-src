/*
 * compute_transformation.h
 *
 *  Created on: Mar 30, 2010
 *      Author: frehse
 */

#ifndef COMPUTE_TRANSFORMATION_H_
#define COMPUTE_TRANSFORMATION_H_

#include "core/continuous/continuous_set.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transform_declarations.h"

//// includes for resolving typeid
//#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC.h"
//#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
//#include "core/continuous/sfm_poly/sfm_cont_set.h"

namespace continuous {

template<typename transform_type>
continuous_set_ptr compute_transformation(const continuous_set& p,
		const transform_type& t);
template<typename transform_type>
continuous_set_ptr compute_or_assign_transformation(const continuous_set& p,
		const transform_type& t);

/** Transformations without default implementations:
 * (will raise an exception if no implementation is provided by the derived class.)
 *
 * This is a class so that partial specializations can be given for derived
 * classes of continuous_set. */

template<typename T> class transformer {
public:
	static continuous_set::ptr compute(const T& c,
			const constant_bound_time_elapse_transform& t) {
		throw std::runtime_error(
				"missing implementation of constant_bound_time_elapse_transform");
		return continuous_set::ptr();
	}
	;
	static continuous_set::ptr compute_or_assign(T& c,
			const constant_bound_time_elapse_transform& t) {
		return compute(c, t);
	}
	;
	static continuous_set::ptr compute(const T& c,
			const reset_affine_transform<global_types::rational_type>& t) {
		std::string name1; // = typeid(c).name();
		throw std::runtime_error(
				"missing implementation of rat reset_affine_transform on "+name1);
		return continuous_set::ptr();
	}
	;
	static continuous_set::ptr compute_or_assign(T& c,
			const reset_affine_transform<global_types::rational_type>& t) {
		return compute(c, t);
	}
	;
	static continuous_set::ptr compute(const T& c,
			const reset_affine_transform<global_types::float_type>& t) {
		std::string name1; // = typeid(c).name();
		throw std::runtime_error(
				"missing implementation of float reset_affine_transform on "+name1);
		return continuous_set::ptr();
	}
	;
	static continuous_set::ptr compute_or_assign(T& c,
			const reset_affine_transform<global_types::float_type>& t) {
		return compute(c, t);
	}
	;
	static continuous_set::ptr compute(const T& c,
			const reset_function_transform& t) {
		throw std::runtime_error(
				"missing implementation of reset_function_transform");
		return continuous_set::ptr();
	}
	static continuous_set::ptr compute_or_assign(T& c,
			const reset_function_transform& t) {
		return compute(c, t);
	}
	;

};
}

/////////////////////////////////////////////////////////////////////
//// Include template specializations for transforms here
/////////////////////////////////////////////////////////////////////
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_transforms.h"
#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC_transforms.h"

namespace continuous {

/** A class to resolve the type of the
 * continuous set and calls an implementor (transformer
 * class),which can be used for partial specialization.
 */
template<typename transform_type>
class compute_transform_caller: public continuous_set::const_visitor {
public:
	compute_transform_caller(const transform_type& t) :
		my_transf(t), my_res(continuous_set::ptr()) {
	}
	;
	continuous_set::ptr result() {
		return my_res;
	}
	;
	void dispatch(const support_function_provider* c) {
		my_res = transformer<support_function_provider >::compute(*c,
				my_transf);
	}
	;
	void dispatch(const constr_polyhedron<Rational>* c) {
		my_res = transformer<constr_polyhedron<Rational> >::compute(*c,
				my_transf);
	}
	;
	void dispatch(const constr_polyhedron<double>* c) {
		my_res
				= transformer<constr_polyhedron<double> >::compute(*c,
						my_transf);
	}
	;
	void dispatch(const ppl_polyhedron::continuous_set_PPL_NNC* c) {
		my_res = transformer<ppl_polyhedron::continuous_set_PPL_NNC>::compute(
				*c, my_transf);
	}
	;
	void dispatch(const support_function::sfm_cont_set<global_types::float_type>* c) {
		my_res = transformer<support_function::sfm_cont_set<global_types::float_type> >::compute(
				*c, my_transf);
	}
	;
	void dispatch(const spacetime_flowpipe<global_types::float_type>* c) {
		my_res = transformer<spacetime_flowpipe<global_types::float_type> >::compute(
				*c, my_transf);
	}
	;
	void dispatch(const predicate_continuous_set* c) {
		my_res = transformer<predicate_continuous_set>::compute(*c, my_transf);
	}
	;
	void dispatch(const continuous_set_simulation<global_types::float_type>* c){

		my_res = transformer<continuous_set_simulation<global_types::float_type> >::compute(*c, my_transf);
	}
	;
private:
	const transform_type& my_transf;
	continuous_set::ptr my_res;
};

template<typename transform_type>
continuous_set_ptr compute_transformation(const continuous_set& p,
		const transform_type& t) {
	compute_transform_caller<transform_type> C(t);
	p.accept(C);
	return C.result();
}
;

/** A class to resolve the type of the
 * continuous set and calls an implementor (transformer
 * class),which can be used for partial specialization.
 *
 * @todo this violates constness, to be fixed
 * (needs nonconst visitor)
 */
template<typename transform_type>
class compute_or_assign_transform_caller: public continuous_set::const_visitor {
public:
	compute_or_assign_transform_caller(const transform_type& t) :
		my_transf(t), my_res(continuous_set::ptr()) {
	}
	;
	continuous_set::ptr result() {
		return my_res;
	}
	;
	void dispatch(const support_function_provider* c) {
		support_function_provider* d =
				const_cast<support_function_provider*> (c);
		my_res = transformer<support_function_provider >::compute_or_assign(
				*d, my_transf);
	}
	;
	void dispatch(const constr_polyhedron<Rational>* c) {
		constr_polyhedron<Rational>* d =
				const_cast<constr_polyhedron<Rational>*> (c);
		my_res = transformer<constr_polyhedron<Rational> >::compute_or_assign(
				*d, my_transf);
	}
	;
	void dispatch(const constr_polyhedron<double>* c) {
		constr_polyhedron<double>* d =
				const_cast<constr_polyhedron<double>*> (c);
		my_res = transformer<constr_polyhedron<double> >::compute_or_assign(*d,
				my_transf);
	}
	;
	void dispatch(const ppl_polyhedron::continuous_set_PPL_NNC* c) {
		ppl_polyhedron::continuous_set_PPL_NNC* d =
				const_cast<ppl_polyhedron::continuous_set_PPL_NNC*> (c);
		my_res
				= transformer<ppl_polyhedron::continuous_set_PPL_NNC>::compute_or_assign(
						*d, my_transf);
	}
	;
	void dispatch(const support_function::sfm_cont_set<global_types::float_type>* c) {
		support_function::sfm_cont_set<global_types::float_type>* d = const_cast<support_function::sfm_cont_set<
				global_types::float_type>*> (c);
		my_res
				= transformer<
						continuous::support_function::sfm_cont_set<global_types::float_type> >::compute_or_assign(
						*d, my_transf);
	}
	;
	void dispatch(const spacetime_flowpipe<global_types::float_type>* c) {
		spacetime_flowpipe<global_types::float_type>* d = const_cast<spacetime_flowpipe<
				global_types::float_type>*> (c);
		my_res
				= transformer<
						continuous::spacetime_flowpipe<global_types::float_type> >::compute_or_assign(
						*d, my_transf);
	}
	;
	void dispatch(const predicate_continuous_set* c) {
		predicate_continuous_set* d = const_cast<predicate_continuous_set*> (c);
		my_res
				= transformer<continuous::predicate_continuous_set>::compute_or_assign(
						*d, my_transf);
	}
	;
	void dispatch(const continuous_set_simulation<global_types::float_type>* c) {
		 continuous_set_simulation<global_types::float_type>* d = const_cast< continuous_set_simulation<global_types::float_type>*> (c);
			my_res
					= transformer<continuous:: continuous_set_simulation<global_types::float_type> >::compute_or_assign(
							*d, my_transf);
		}
		;
private:
	const transform_type& my_transf;
	continuous_set::ptr my_res;
};

template<typename transform_type>
continuous_set_ptr compute_or_assign_transformation(const continuous_set& p,
		const transform_type& t) {
	compute_or_assign_transform_caller<transform_type> C(t);
	p.accept(C);
	return C.result();
}
;

}

#endif /* COMPUTE_TRANSFORMATION_H_ */
