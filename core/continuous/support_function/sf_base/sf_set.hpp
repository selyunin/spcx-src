/*
 * sf_set.hpp
 *
 *  Created on: Apr 1, 2010
 *      Author: frehse
 */

#ifndef SF_SET_HPP_
#define SF_SET_HPP_

#include "sf_set.h"
#include "sf_unary.h"

#include "math/vdom/vdom_matrix_operators.h"
#include "math/vdom/vdom_vector_operators.h"
#include "math/vdom/affine_map_utility.h"
#include "math/vdom/positional_vdomain.h"

namespace continuous {
namespace support_function {

template<typename scalar_type> sf_set<scalar_type>::sf_set() : my_map(affine_map_ptr()) {
}
;

template<typename scalar_type> sf_set<scalar_type>::sf_set(const affine_map& M) {
	affine_map_ptr mp = affine_map_ptr(new affine_map(M));
	my_map = mp;
}
;

template<typename scalar_type> sf_set<scalar_type>::~sf_set() {
}
;

template<typename scalar_type> typename sf_set<scalar_type>::poly_ptr sf_set<
		scalar_type>::outer_poly(const vector_set& dirs) const {
	poly_ptr res = poly_ptr(new poly());
	for (typename vector_set::const_iterator it = dirs.begin(); it
			!= dirs.end(); ++it) {
		double max_value;
		math::vdom_vector<double> support_vec;
		bool is_empty;
		bool is_bounded;
		math::vdom_vector<scalar_type> l(*it);
		this->compute_support(l.template convert_to<double> (), max_value,
				support_vec, is_empty, is_bounded);
		if (is_empty) {
			return poly_ptr(res->create_empty());
		} else {
			if (is_bounded) {
				lin_constraint con(*it, convert_element<scalar_type> (
						-max_value), LE);
				res->add_constraint(con);
			}
		}
	}
	return res;
}
;

template<typename scalar_type> sf_set<scalar_type> * sf_set<scalar_type>::create_universe() const {
	continuous_set* c = new poly();
	support_function_provider* s = static_cast<support_function_provider*> (c);
	support_function_provider::ptr sp = support_function_provider::ptr(s);
	return new sf_unary<scalar_type>(sp);
}

template<typename scalar_type> sf_set<scalar_type> * sf_set<scalar_type>::create_empty() const {
	continuous_set* c = new poly(poly::empty_poly());
	support_function_provider* s = static_cast<support_function_provider*> (c);
	support_function_provider::ptr sp = support_function_provider::ptr(s);
	return new sf_unary<scalar_type>(sp);
}

template<typename scalar_type> unsigned int sf_set<scalar_type>::get_dim() const {
	return get_variable_ids().size();
}

template<typename scalar_type> void sf_set<scalar_type>::embed_variables(
		const variable_id_set& id_set) {
	// there would be the following implementation, but I'm not sure what
	// it would be good for. embedding is usually followed by intersection
	// which would result in an awkward construction

	// add zero rows to A and b, then turn self into minkowski sum
	// with universe
	throw std::runtime_error("sf_set : missing implementation embed_variables");
}

template<typename scalar_type> void sf_set<scalar_type>::existentially_quantify_variables(
		const variable_id_set& id_set) {
	// remove rows from A and b
	variable_id_set vars = get_variable_ids();
	variable_id_set rem_vars;
	std::set_intersection(vars.begin(), vars.end(), id_set.begin(),
			id_set.end(), std::inserter(rem_vars, rem_vars.begin()));

	std::set<variable> rem_variables = create_variable_set(rem_vars);
	positional_vdomain old_codom = my_map->get_A().codomain();
	positional_vdomain new_codom = my_map->get_A().codomain();
	const positional_vdomain& new_dom = my_map->get_A().domain();
	new_codom.remove_variables(rem_variables);

	if (new_codom.size() > 0) {

		typename affine_map::matrix_type
				newAm(new_codom.size(), new_dom.size());
		typename affine_map::vector_type newbm(new_codom.size());

		// Copy the rows of A and b
		for (unsigned int i = 0; i < new_codom.size(); ++i) {
			unsigned int iold = old_codom.pos(new_codom.get_variable(i));
			for (unsigned int j = 0; i < new_dom.size(); ++j) {
				newAm(i, j) = my_map->get_A()(iold, j);
			}
			newbm[i] = my_map->get_b()[iold];
		}

		typename affine_map::vdom_matrix_type newA(new_codom, new_dom, newAm);
		typename affine_map::vdom_vector_type newb(new_codom, newbm);

		*my_map=affine_map(newA,newb);

	} else {
		throw std::runtime_error(
				"missing implementation existentially_quantify_variables");
	}
}

//template<typename scalar_type> void sf_set<scalar_type>::accept(
//		dispatching::dispatcher<continuous_set_typelist>& d) const {
//	d.dispatch(this);
//}

template<typename scalar_type>
template<typename fun_type> void sf_set<scalar_type>::compute_support_mapped(
		const support_function_provider& a_set,
		const affine_map_const_ptr& a_map,
		const math::vdom_vector<fun_type>& v, fun_type& max_value,
		math::vdom_vector<fun_type>& support_vec, bool& is_empty,
		bool& is_bounded) {
	if (a_map) {
		math::affine_map<fun_type> M = a_map->template convert_to<
				fun_type> ();
		// transform the cost function
		math::vdom_vector<fun_type> lmapped(v * (M.get_A()));
		a_set.compute_support(lmapped, max_value, support_vec, is_empty,
				is_bounded);
		max_value += scalar_product(v, M.get_b());
		// transform the support vector if there is one
		if (support_vec.size() > 0) {
			support_vec = M.map(support_vec);
		}
	} else {
		a_set.compute_support(v, max_value, support_vec, is_empty, is_bounded);
	}
}

//template<typename scalar_type>
//void sf_set<scalar_type>::compute_support_mapped(
//		const support_function_provider& a_set,
//		const affine_map_const_ptr& a_map,
//		const math::vdom_vector<scalar_type>& v, scalar_type& max_value,
//		math::vdom_vector<scalar_type>& support_vec, bool& is_empty,
//		bool& is_bounded) {
//	if (a_map) {
//		// transform the cost function
//		math::vdom_vector<scalar_type> lmapped(v * (a_map->get_A()));
//		a_set.compute_support(lmapped, max_value, support_vec, is_empty,
//				is_bounded);
//		max_value += scalar_product(v, a_map->get_b());
//		// transform the support vector if there is one
//		if (support_vec.size() > 0) {
//			support_vec = a_map->map(support_vec);
//		}
//	} else {
//		a_set.compute_support(v, max_value, support_vec, is_empty, is_bounded);
//	}
//}

template<>
template<>
inline void sf_set<double>::compute_support_mapped(
		const support_function_provider& a_set,
		const affine_map_const_ptr& a_map,
		const math::vdom_vector<double>& v, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) {
	if (a_map) {
		if (!a_map->is_translation()) {
			// transform the cost function
			math::vdom_vector<double> lmapped(v * (a_map->get_A()));
			a_set.compute_support(lmapped, max_value, support_vec, is_empty,
					is_bounded);
			// transform the support vector if there is one
			if (support_vec.size() > 0) {
				support_vec = a_map->map(support_vec);
			}
		} else {
			// A is the identity matrix
			a_set.compute_support(v, max_value, support_vec, is_empty,
					is_bounded);
			if (support_vec.size() > 0) {
				support_vec += a_map->get_b();
			}
		}
		max_value += scalar_product(v, a_map->get_b());
	} else {
		a_set.compute_support(v, max_value, support_vec, is_empty, is_bounded);
	}
}

template<typename scalar_type> void sf_set<scalar_type>::map(
		const affine_map& M) {
	affine_map_ptr mp;
	if (my_map) {
		mp = affine_map_ptr(
				new affine_map(math::concatenate(M, *my_map)));
	} else {
		mp = affine_map_ptr(new affine_map(M));
	}
	my_map = mp;
}

template<typename scalar_type>
typename sf_set<scalar_type>::affine_map_const_ptr sf_set<scalar_type>::get_map() const {
	return my_map;
}


template<typename scalar_type>
typename sf_set<scalar_type>::affine_map_ptr sf_set<scalar_type>::get_map() {
	return my_map;
}

template<typename scalar_type>
void sf_set<scalar_type>::set_map(const affine_map& M) {
	*my_map = M;
}


}
}

#endif /* SF_SET_HPP_ */
