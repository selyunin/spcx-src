/*
 * finite_finite_hyperbox.hpp
 *
 *  Created on: Jan 17, 2011
 *      Author: frehse
 */

#ifndef FINITE_HYPERBOX_HPP_
#define FINITE_HYPERBOX_HPP_

#include "finite_hyperbox.h"

#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"
#include "math/vector_utility.h"
#include "math/vdom/affine_map_utility.h"

namespace continuous {

template<class scalar_type> finite_hyperbox<scalar_type>::finite_hyperbox() {
}

template<class scalar_type> finite_hyperbox<scalar_type>::finite_hyperbox(
		point_type c, point_type g, positional_vdomain dom) :
	index_to_variable_id_map_provider(dom), my_c(c), my_g(g) {
	assert(c.size() == g.size());
	assert(c.size() == dom.size());
	assert(c.size()>=1);

	canonicalize();
}

template<class scalar_type> finite_hyperbox<scalar_type>::finite_hyperbox(
		vdom_vector_type c, vdom_vector_type g) :
	index_to_variable_id_map_provider(c.domain()), my_c(c.get_vector()),
			my_g(g.get_vector()) {
	assert(c.size() == g.size());
	assert(c.size() == c.domain().size());
	assert(c.size()>=1);
	assert(c.domain() == g.domain());

	canonicalize();
}

template<class scalar_type>
finite_hyperbox<scalar_type>::finite_hyperbox(
		const finite_hyperbox<scalar_type>& box) :
	index_to_variable_id_map_provider(box), my_c(box.my_c), my_g(box.my_g) {
	canonicalize();
}

template<class scalar_type>
finite_hyperbox<scalar_type>& finite_hyperbox<scalar_type>::operator=(
		const finite_hyperbox<scalar_type>& box) {
	index_to_variable_id_map_provider::operator=(box);
	my_c = box.my_c;
	my_g = box.my_g;
	return *this;
}

template<class scalar_type> continuous_set::ptr finite_hyperbox<scalar_type>::get_ptr() {
	continuous_set::ptr p = boost::enable_shared_from_this<finite_hyperbox<
			scalar_type> >::shared_from_this();
	return p;
}

template<class scalar_type> continuous_set::const_ptr finite_hyperbox<
		scalar_type>::get_const_ptr() const {
	continuous_set::const_ptr p = boost::enable_shared_from_this<
			finite_hyperbox<scalar_type> >::shared_from_this();
	return p;
}

template<class scalar_type> finite_hyperbox<scalar_type>* finite_hyperbox<
		scalar_type>::create_universe() const {
	throw std::runtime_error(
			"finite_hyperbox<scalar_type>::create_universe() not allowed");
	return new finite_hyperbox<scalar_type> (point_type(), point_type(),
			positional_vdomain());
}

template<class scalar_type> finite_hyperbox<scalar_type>* finite_hyperbox<
		scalar_type>::clone() const {
	return new finite_hyperbox<scalar_type> (*this);
}

template<class scalar_type> finite_hyperbox<scalar_type>* finite_hyperbox<
		scalar_type>::create_empty() const {
	throw std::runtime_error(
			"finite_hyperbox<scalar_type>::create_empty() not allowed");
	point_type new_g = my_g;
	new_g[0] = -scalar_type(1);
	return new finite_hyperbox<scalar_type> (my_c, new_g, domain());
}

template<class scalar_type> int finite_hyperbox<scalar_type>::get_memory() const {
	throw std::runtime_error(
			"finite_hyperbox<scalar_type>::get_memory() missing implementation");
	return 0;
}

template<class scalar_type> unsigned int finite_hyperbox<scalar_type>::get_dim() const {
	return get_index_to_variable_id_map()->dimensions();
}

template<class scalar_type> math::tribool finite_hyperbox<scalar_type>::is_empty() const {
	if (my_g.size() == 0)
		return true;
	return my_g[0] <= -value_type(1);
}

template<class scalar_type> math::tribool finite_hyperbox<scalar_type>::is_universe() const {
	return false;
}

template<class scalar_type> void finite_hyperbox<scalar_type>::canonicalize() {
	for (unsigned int i = 0; i < my_g.size(); ++i) {
		if (math::definitely(math::numeric::is_LT(my_g[i], scalar_type(0)))) {
			set_empty();
			return;
		}
	}
}

template<class scalar_type> void finite_hyperbox<scalar_type>::embed_variables(
		const variable_id_set& id_set) {
	throw std::runtime_error(
			"finite_hyperbox<scalar_type>::embed_variables() not allowed");
}

template<class scalar_type> void finite_hyperbox<scalar_type>::existentially_quantify_variables(
		const variable_id_set& id_set) {
	index_to_variable_id_map_ptr p =
			get_index_to_variable_id_map()->get_map_with_ids_removed(id_set);

	// remap my_l and my_u
	index_to_index_bimap adapt_map = get_index_to_index_mapping(p,
			get_index_to_variable_id_map());

	// it doesn't matter whether things are updated or not, since
	// they keep their (updated or not) value
	my_c = math::map(my_g, adapt_map);
	my_c = math::map(my_g, adapt_map);
	// since existentially_quantify_variables only removes variables,
	// we don't need to define any uninitialized elements

	// fix _index_to_variable_id_map_ptr
	set_index_to_variable_id_map(p);
}

template<class scalar_type>
void finite_hyperbox<scalar_type>::compute_support_impl(
		const math::vdom_vector<scalar_type>& l, value_type& max_value,
		math::vdom_vector<value_type>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	if (this->is_empty()) {
		is_empty = true;
		return;
	} else {
		is_empty = false;
		is_bounded = true;
		const index_to_variable_id_map_ptr& my_map =
				get_index_to_variable_id_map();
		const index_to_variable_id_map_ptr& l_map =
				l.get_index_to_variable_id_map();

		static scalar_type zero = scalar_type(0);
		unsigned int K = my_map->dimensions();

		max_value = zero;
		// initialize support vector with zero, so irrelevant variables take
		// the value zero
		support_vec = math::vdom_vector<value_type>(domain(),
				math::vector<value_type>(K));

		bool remap_variables = (l_map != my_map);
		value_type max_x;

		if (remap_variables) {
			typename math::lin_expression<scalar_type>::size_type N =
					l_map->dimensions();

			// need to zero the support vector
			for (typename math::vdom_vector<value_type>::size_type i = 0; i
								< K; ++i) {
				support_vec[i]=zero;
			}

			size_type my_index;
			bool found;
			for (typename math::lin_expression<scalar_type>::size_type i = 0; i
					< N; ++i) {
				if (l[i] != zero) {
					my_index = my_map->check_for_index(l_map->get_id(i), found);
					if (found) {
						max_x = (l[i] > zero) ? my_c[my_index] + my_g[my_index]
								: my_c[my_index] - my_g[my_index];
					} else {
						is_bounded = false;
						return;
					}
					max_value += l[i] * max_x;
					support_vec[my_index] = max_x;
				} else {
					support_vec[my_index] = my_c[my_index];
				}
			}
		} else {
			for (typename math::lin_expression<scalar_type>::size_type i = 0; i
					< K; ++i) {
				if (l[i] != zero) {
					max_x = (l[i] > zero) ? my_c[i] + my_g[i] : my_c[i]
							- my_g[i];

					max_value += l[i] * max_x;
					support_vec[i] = max_x;
				} else {
					support_vec[i] = my_c[i];
				}
			}
		}
	}
}

template<typename scalar_type> bool finite_hyperbox<scalar_type>::computes_support_vector() const {
	return true;
}
;

template<class scalar_type>
void finite_hyperbox<scalar_type>::compute_support(
		const math::vdom_vector<Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	value_type max_val_st;
	math::vdom_vector<value_type> support_vec_st;
	math::vdom_vector<scalar_type> l_scalar_type = l.template convert_to<
			scalar_type> ();

	compute_support_impl(l_scalar_type, max_val_st, support_vec_st, is_empty,
			is_bounded);
	if (!is_empty && is_bounded) {
		max_value = convert_element<Rational> (max_val_st);
		support_vec = support_vec_st.template convert_to<Rational> ();
	}
}

template<class scalar_type>
void finite_hyperbox<scalar_type>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	value_type max_val_st;
	math::vdom_vector<value_type> support_vec_st;
	math::vdom_vector<scalar_type> l_scalar_type = l.template convert_to<
			scalar_type> ();

	compute_support_impl(l_scalar_type, max_val_st, support_vec_st, is_empty,
			is_bounded);
	if (!is_empty && is_bounded) {
		max_value = convert_element<double> (max_val_st);
		support_vec = support_vec_st.template convert_to<double> ();
	}
}

template<>
inline void finite_hyperbox<double>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {

	compute_support_impl(l, max_value, support_vec, is_empty, is_bounded);

}

template<class scalar_type>
typename math::lin_constraint_system<scalar_type>::const_ptr finite_hyperbox<
		scalar_type>::get_constraints() const {
	typename math::lin_constraint_system<scalar_type>::ptr cons(
			new math::lin_constraint_system<scalar_type>());

	math::lin_constraint<scalar_type> c;
	if (!math::definitely(is_empty())) {
		for (size_type i = 0; i < get_dim(); ++i) {
			math::lin_expression<scalar_type> e(get_index_to_variable_id_map());
			e[i] = scalar_type(1);
			e.set_inh_coeff(-(my_c[i] + my_g[i]));
			c = math::lin_constraint<scalar_type>(e, LE);
			cons->insert(c);

			e[i] = -scalar_type(1);
			e.set_inh_coeff(my_c[i] - my_g[i]);
			c = math::lin_constraint<scalar_type>(e, LE);
			cons->insert(c);
		}
	} else {
		c = math::lin_constraint<scalar_type>::zero_dim_false();
		cons->insert(c);
	}
	return cons;
}

template<class scalar_type>
void finite_hyperbox<scalar_type>::add_constraint(
		const math::lin_constraint<scalar_type> &c, bool check_redundancy) {
	throw std::runtime_error(
			"finite_hyperbox<scalar_type>::add_constraint() missing implementation");
}

template<class scalar_type>
void finite_hyperbox<scalar_type>::remove_redundant_constraints() {
	// do nothing
}

template<class scalar_type>
const typename finite_hyperbox<scalar_type>::point_type& finite_hyperbox<
		scalar_type>::get_c() const {
	return my_c;
}

template<class scalar_type>
const typename finite_hyperbox<scalar_type>::point_type& finite_hyperbox<
		scalar_type>::get_g() const {
	return my_g;
}

template<class scalar_type>
typename finite_hyperbox<scalar_type>::point_type finite_hyperbox<
		scalar_type>::lower() const {
	return my_c-my_g;
}

template<class scalar_type>
typename finite_hyperbox<scalar_type>::point_type finite_hyperbox<
		scalar_type>::upper() const {
	return my_c+my_g;
}

template<class scalar_type> typename finite_hyperbox<scalar_type>::vdom_vector_type finite_hyperbox<
		scalar_type>::get_c_dom() const {
	return vdom_vector_type(my_c, get_index_to_variable_id_map());
}

template<class scalar_type> typename finite_hyperbox<scalar_type>::vdom_vector_type finite_hyperbox<
		scalar_type>::get_g_dom() const {
	return vdom_vector_type(my_g, get_index_to_variable_id_map());
}

template<class scalar_type>
void finite_hyperbox<scalar_type>::reorder(const positional_vdomain& dom) {
	math::vdom_vector<value_type> new_c(dom, my_c);
	math::vdom_vector<value_type> new_g(dom, my_g);
	new_c.reorder(dom);
	new_g.reorder(dom);
	my_c = new_c.get_vector();
	my_g = new_g.get_vector();
	set_domain(dom);
}

template<class scalar_type>
void finite_hyperbox<scalar_type>::set_empty() {
	if (domain().size() > 0) {
		point_type c;
		point_type g;
		my_c = point_type(domain().size(), value_type(0)); // a vector of ones
		my_g = point_type(domain().size(), -value_type(1)); // a vector of zeroes
	} else {
		throw std::runtime_error(
				"finite_hyperbox<scalar_type>::empty_box() not allowed with empty domain");
	}
}

template<class scalar_type>
finite_hyperbox<scalar_type> finite_hyperbox<scalar_type>::empty_box(
		const positional_vdomain& dom) {
	point_type c;
	point_type g;
	if (dom.size() > 0) {
		c = point_type(dom.size(), value_type(0)); // a vector of ones
		g = point_type(dom.size(), -value_type(1)); // a vector of zeroes
		return finite_hyperbox<scalar_type> (c, g, dom);
	} else {
		throw std::runtime_error(
				"finite_hyperbox<scalar_type>::empty_box() not allowed with empty domain");
		return finite_hyperbox<scalar_type> (c, g, dom);
	}
}

template<class scalar_type>
bool finite_hyperbox<scalar_type>::operator==(
		const finite_hyperbox<scalar_type>& h) const {
	if (is_empty())
		return h.is_empty();
	else if (h.domain() == domain())
		return get_c() == h.get_c() && get_g() == h.get_g();
	else {
		return get_c_dom() == h.get_c_dom() && get_g_dom() == h.get_g_dom();
	}
}

template<class scalar_type>
finite_hyperbox<scalar_type> finite_hyperbox<scalar_type>::operator+(
		const finite_hyperbox<scalar_type>& h) const {
	if (domain() == h.domain()) {
		return finite_hyperbox<scalar_type> (get_c() + h.get_c(),
				get_g() + h.get_g(), domain());
	} else {
		return finite_hyperbox<scalar_type> (get_c_dom() + h.get_c_dom(),
				get_g_dom() + h.get_g_dom());
	}
}

template<class scalar_type>
finite_hyperbox<scalar_type>& finite_hyperbox<scalar_type>::operator+=(
		const finite_hyperbox<scalar_type>& h) {
	if (domain() == h.domain()) {
		my_c += h.get_c();
		my_g += h.get_g();
	} else {
		*this = *this + h;
	}
	return *this;
}

template<class scalar_type>
bool finite_hyperbox<scalar_type>::contains(const vdom_vector_type& p) const {
	// p is contained if abs(p-my_c)<=my_g.
	// Recall that my_g>=0.

	// copy so we can remap
	// @todo this is not optimal
	vdom_vector_type q(p);
	q.remap(domain());
	// take the absolute value of x
	point_type x = vec_abs(q.get_vector() - my_c);
	return !math::numeric::is_lex_LT(my_g, x);
}

}

#endif /* FINITE_HYPERBOX_HPP_ */
