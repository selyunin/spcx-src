/*
 * hyperbox.hpp
 *
 *  Created on: Nov 25, 2009
 *      Author: frehse
 */

#ifndef HYPERBOX_HPP_
#define HYPERBOX_HPP_

#include "hyperbox.h"

#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"
#include "math/vector_utility.h"
#include "math/vdom/affine_map_utility.h"

namespace continuous {

template<class scalar_type> hyperbox<scalar_type>::hyperbox() {
}

template<class scalar_type> hyperbox<scalar_type>::hyperbox(point_type l,
		point_type u, index_to_variable_id_map_ptr iimap) : index_to_variable_id_map_provider(iimap),
	my_l(l), my_u(u) {
	assert(l.size() == u.size());
	assert(l.size() == iimap->dimensions());

	my_zero_dim_empty = false; // the only way to create a zero_dim_empty is via empty_box

	//assert(check_consistency());
}


template<class scalar_type> hyperbox<scalar_type>::hyperbox(point_type l,
		point_type u, const positional_vdomain& dom) : index_to_variable_id_map_provider(dom),
	my_l(l), my_u(u) {
	assert(l.size() == u.size());
	assert(l.size() == dom.size());

	my_zero_dim_empty = false; // the only way to create a zero_dim_empty is via empty_box

	//assert(check_consistency());
}

template<class scalar_type> hyperbox<scalar_type>::hyperbox(
		finite_point_type l, finite_point_type u,
		index_to_variable_id_map_ptr iimap) : index_to_variable_id_map_provider(iimap) {
	my_l =
			l.template convert_to<scalar_with_infinity<scalar_type> > ();
	my_u =
			u.template convert_to<scalar_with_infinity<scalar_type> > ();

	my_zero_dim_empty = false; // the only way to create a zero_dim_empty is via empty_box

	//assert(check_consistency());
}

template<class scalar_type> hyperbox<scalar_type>::hyperbox(
		index_to_variable_id_map_ptr iimap) : index_to_variable_id_map_provider(iimap) {
	my_l = point_type(iimap->dimensions(), value_type::NaN());
	my_u = my_l;

	my_zero_dim_empty = false; // the only way to create a zero_dim_empty is via empty_box

	//assert(check_consistency());
}

template<class scalar_type> hyperbox<scalar_type>::hyperbox(
		const positional_vdomain& dom) : index_to_variable_id_map_provider(dom) {
	my_l = point_type(dom.size(), value_type::NaN());
	my_u = my_l;

	my_zero_dim_empty = false; // the only way to create a zero_dim_empty is via empty_box

	//assert(check_consistency());
}

template<class scalar_type>
hyperbox<scalar_type>::hyperbox(const hyperbox<scalar_type>& box) :
	index_to_variable_id_map_provider(box), my_l(box.my_l), my_u(box.my_u),
			my_zero_dim_empty(box.my_zero_dim_empty) {
}

template<class scalar_type>
hyperbox<scalar_type>& hyperbox<scalar_type>::operator=(const hyperbox<
		scalar_type>& box) {
	index_to_variable_id_map_provider::operator=(box);
	my_l = box.my_l;
	my_u = box.my_u;
	my_zero_dim_empty = box.my_zero_dim_empty;
	return *this;
}

template<class scalar_type> hyperbox<scalar_type>::~hyperbox() {
}

template<class scalar_type> continuous_set::ptr hyperbox<scalar_type>::get_ptr() {
	continuous_set::ptr p = boost::enable_shared_from_this<
			hyperbox<scalar_type> >::shared_from_this();
	return p;
}

template<class scalar_type> continuous_set::const_ptr hyperbox<scalar_type>::get_const_ptr() const {
	continuous_set::const_ptr p = boost::enable_shared_from_this<hyperbox<
			scalar_type> >::shared_from_this();
	return p;
}

template<class scalar_type> hyperbox<scalar_type>* hyperbox<scalar_type>::create_universe() const {
	return new hyperbox<scalar_type> ();
}

template<class scalar_type> hyperbox<scalar_type>* hyperbox<scalar_type>::clone() const {
	return new hyperbox<scalar_type> (*this);
}

template<class scalar_type> hyperbox<scalar_type>* hyperbox<scalar_type>::create_empty() const {
	return new hyperbox<scalar_type> (empty_box());
}

template<class scalar_type> int hyperbox<scalar_type>::get_memory() const {
	throw std::runtime_error(
			"hyperbox<scalar_type>::get_memory() missing implementation");
	return 0;
}

template<class scalar_type> unsigned int hyperbox<scalar_type>::get_dim() const {
	return get_index_to_variable_id_map()->dimensions();
}

template<class scalar_type> math::tribool hyperbox<scalar_type>::is_empty() const {
	if (my_zero_dim_empty)
		return true;
	else
		return math::numeric::is_lex_LT(get_u(),get_l());
}

template<class scalar_type> math::tribool hyperbox<scalar_type>::is_universe() const {
	if (is_empty()) return false;
	for (unsigned int i = 0; i < get_dim(); ++i) {
		if (get_l(i).is_finite())
			return false;
		if (get_u(i).is_finite())
			return false;
	}
	return true;
}

template<class scalar_type> void hyperbox<scalar_type>::embed_variables(
		const variable_id_set& id_set) {

	typename point_type::size_type old_dim =
			get_index_to_variable_id_map()->dimensions();
	// add the ids to _index_to_variable_id_map_ptr
	index_to_variable_id_map_ptr p =
			get_index_to_variable_id_map()->get_map_with_ids_added(id_set);
	set_index_to_variable_id_map(p);
	// resize the points (new dimensions are added at the end)
	my_l.resize(p->dimensions());
	my_u.resize(p->dimensions());
	// set bounds to +/- infinity in new dimensions
	for (typename point_type::size_type i = old_dim; i < p->dimensions(); ++i) {
		set_l(i, value_type::neg_infty());
		set_u(i, value_type::pos_infty());
	}

	// fix zero_dim_emptiness flag if necessary
	if (my_zero_dim_empty) {
		my_zero_dim_empty = (get_dim() > 0);
	}
}

template<class scalar_type> void hyperbox<scalar_type>::existentially_quantify_variables(
		const variable_id_set& id_set) {
	if (!my_zero_dim_empty) {
		index_to_variable_id_map_ptr
				p = get_index_to_variable_id_map()->get_map_with_ids_removed(
						id_set);

		// remap my_l and my_u
		index_to_index_bimap adapt_map = get_index_to_index_mapping(p,
				get_index_to_variable_id_map());

		// it doesn't matter whether things are updated or not, since
		// they keep their (updated or not) value
		my_l=math::map(my_l,adapt_map);
		my_u=math::map(my_u,adapt_map);
		// since existentially_quantify_variables only removes variables,
		// we don't need to define any uninitialized elements

		// fix _index_to_variable_id_map_ptr
		set_index_to_variable_id_map(p);
	}
}

template<class scalar_type>
void hyperbox<scalar_type>::compute_support(const math::vdom_vector<
		scalar_type>& l, value_type& max_value, math::vdom_vector<
		value_type>& support_vec) const {
	if (is_empty()) {
		max_value = value_type::NaN();
	} else {
		index_to_variable_id_map_ptr my_map = get_index_to_variable_id_map();
		index_to_variable_id_map_ptr l_map = l.get_index_to_variable_id_map();

		max_value = value_type(0);
		// initialize support vector with zero, so irrelevant variables take
		// the value zero
		support_vec = math::vdom_vector<value_type>(
				math::vector<value_type>(my_map->dimensions(),
						value_type(0)), my_map);

		bool remap_variables = (l_map != my_map);
		value_type max_x;
		for (typename math::lin_expression<scalar_type>::size_type i = 0; i
				< l_map->dimensions(); ++i) {
			if (l[i] != scalar_type(0)) {
				size_type my_index = i;
				bool found = true;
				if (remap_variables) {
					my_index = my_map->check_for_index(l_map->get_id(i), found);
				}
				if (found) {
					if (l[i] > scalar_type(0)) {
						// use upper bound
						max_x = get_u(my_index);
					} else {
						// use lower bound
						max_x = get_l(my_index);
					}
				} else {
					// use + infinity (because we're maximising)
					max_x = value_type::pos_infty();
				}
				max_value += value_type(l[i]) * max_x;
				support_vec[my_index] = max_x;
			}
		}
	}
}

template<typename scalar_type> bool hyperbox<scalar_type>::computes_support_vector() const {
	return true;
}
;

template<class scalar_type>
void hyperbox<scalar_type>::compute_support(const math::vdom_vector<
		Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	value_type max_val_with_infty;
	math::vdom_vector<value_type> support_vec_with_infty;
	math::vdom_vector<scalar_type> l_scalar_type =
			l.template convert_to<scalar_type> ();

	compute_support(l_scalar_type, max_val_with_infty, support_vec_with_infty);
	if (max_val_with_infty.is_NaN()) {
		is_empty = true;
		is_bounded = true;
	} else {
		is_empty = false;
		if (max_val_with_infty.is_finite()) {
			is_bounded = true;
			max_value
					= convert_element<Rational> (max_val_with_infty.get_val());
			support_vec = math::vdom_vector<Rational>(
					get_index_to_variable_id_map());
			for (unsigned int i = 0; i < support_vec_with_infty.size(); ++i) {
				support_vec[i] = convert_element<Rational> (
						support_vec_with_infty[i].get_val());
			}
		} else {
			is_bounded = false;
		}
	}
}

template<class scalar_type>
void hyperbox<scalar_type>::compute_support(const math::vdom_vector<
		double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	value_type max_val_with_infty;
	math::vdom_vector<value_type> support_vec_with_infty;
	math::vdom_vector<scalar_type> l_scalar_type =
			l.template convert_to<scalar_type> ();

	compute_support(l_scalar_type, max_val_with_infty, support_vec_with_infty);
	if (max_val_with_infty.is_NaN()) {
		is_empty = true;
		is_bounded = true;
	} else {
		is_empty = false;
		if (max_val_with_infty.is_finite()) {
			is_bounded = true;
			max_value = convert_element<double> (max_val_with_infty.get_val());
			support_vec = math::vdom_vector<double>(
					get_index_to_variable_id_map());
			for (unsigned int i = 0; i < support_vec_with_infty.size(); ++i) {
				support_vec[i] = convert_element<double> (
						support_vec_with_infty[i].get_val());
			}
		} else {
			is_bounded = false;
		}
	}
}

template<class scalar_type>
typename math::lin_constraint_system<scalar_type>::const_ptr hyperbox<
		scalar_type>::get_constraints() const {
	typename math::lin_constraint_system<scalar_type>::ptr cons(
			new math::lin_constraint_system<scalar_type>());

	math::lin_constraint<scalar_type> c;
	if (!math::definitely(is_empty())) {
		for (size_type i = 0; i < get_dim(); ++i) {
			math::lin_expression<scalar_type> e(
					get_index_to_variable_id_map());
			if (get_u(i).is_finite()) {
				e[i] = scalar_type(1);
				e.set_inh_coeff(-get_u(i).get_val());
				c = math::lin_constraint<scalar_type>(e, LE);
				cons->insert(c);
			}
			if (get_l(i).is_finite()) {
				e[i] = -scalar_type(1);
				e.set_inh_coeff(get_l(i).get_val());
				c = math::lin_constraint<scalar_type>(e, LE);
				cons->insert(c);
			}
		}
	} else {
		c = math::lin_constraint<scalar_type>::zero_dim_false();
		cons->insert(c);
	}
	return cons;
}

template<class scalar_type>
void hyperbox<scalar_type>::add_constraint(const math::lin_constraint<
		scalar_type> &c, bool check_redundancy) {
	throw std::runtime_error(
			"hyperbox<scalar_type>::add_constraint() missing implementation");
}

template<class scalar_type>
void hyperbox<scalar_type>::remove_redundant_constraints() {
	// do nothing
}

template<class scalar_type>
bool hyperbox<scalar_type>::is_l_finite() const {
	for (unsigned int i = 0; i < get_dim(); ++i) {
		if (!get_l(i).is_finite())
			return false;
	}
	return true;
}

template<class scalar_type>
bool hyperbox<scalar_type>::is_u_finite() const {
	for (unsigned int i = 0; i < get_dim(); ++i) {
		if (!get_u(i).is_finite())
			return false;
	}
	return true;
}

template<class scalar_type>
bool hyperbox<scalar_type>::is_finite() const {
	return (is_l_finite() && is_u_finite()) || is_empty();
}

template<class scalar_type>
typename hyperbox<scalar_type>::point_type hyperbox<scalar_type>::compute_center() const {
	if (is_empty())
		throw std::runtime_error("can't compute center of empty hyperbox");
	typename hyperbox<scalar_type>::point_type center;
	try {
		center = (get_l() + get_u()) / scalar_with_infinity<scalar_type> (
				scalar_type(2));
	} catch ( std::exception& e )  {
		std::stringstream s;
		s << *this;
		throw basic_exception("Can't compute center of the following hyperbox: "
				+ s.str(), e);
	}
	return center;
}

template<class scalar_type>
typename hyperbox<scalar_type>::vdom_vector_type hyperbox<scalar_type>::compute_finite_center() const {
	typename hyperbox<scalar_type>::point_type nonfinite_center = compute_center();
	typename hyperbox<scalar_type>::vdom_vector_type center;
	variable_id_set infinite_dims;
	try {
//		center = vdom_vector_type(domain(),
//				nonfinite_center.template convert_to<scalar_type> ());
		// convert variable by variable so we can get a nicer error message
		center = vdom_vector_type(domain());
		for (unsigned int i = 0; i < domain().size(); ++i) {
			if (nonfinite_center[i].is_finite()) {
				center[i] = nonfinite_center[i].get_val();
			} else {
				infinite_dims.insert(domain().get_variable(i).get_id());
			}
		}
	} catch (std::exception& e) {
		std::stringstream s;
		s << *this << " with (non-finite) center at " << nonfinite_center;
		throw basic_exception("Can't compute center of the following hyperbox:\n" + s.str(), e);
	}
	if (!infinite_dims.empty()) {
		std::stringstream sv;
		print_variable_id_set(sv,infinite_dims);
		throw basic_exception(
				"Unbounded values for variable(s) "+sv.str()+".");
	}
	return center;
}

template<class scalar_type>
const typename hyperbox<scalar_type>::point_type& hyperbox<scalar_type>::get_l() const {
	return my_l;
}

template<class scalar_type>
const typename hyperbox<scalar_type>::point_type& hyperbox<scalar_type>::get_u() const {
	return my_u;
}

template<class scalar_type> typename hyperbox<scalar_type>::vdom_vector_type hyperbox<
		scalar_type>::get_finite_l() const {
	return vdom_vector_type(my_l.template convert_to<scalar_type> (),
			get_index_to_variable_id_map());
}

template<class scalar_type> typename hyperbox<scalar_type>::vdom_vector_type hyperbox<
		scalar_type>::get_finite_u() const {
	return vdom_vector_type(my_u.template convert_to<scalar_type> (),
			get_index_to_variable_id_map());
}

template<class scalar_type>
hyperbox<scalar_type> hyperbox<scalar_type>::empty_box(
		index_to_variable_id_map_ptr iimap) {
	point_type l;
	point_type u;
	if (iimap->dimensions() > 0) {
		l = point_type(iimap->dimensions(), value_type(1)); // a vector of ones
		u = point_type(iimap->dimensions(), value_type(0)); // a vector of zeroes
		return hyperbox<scalar_type> (l, u, iimap);
	} else {
		hyperbox<scalar_type> new_h(l, u, iimap);
		new_h.my_zero_dim_empty = true;
		return new_h;
	}
}

template<class scalar_type>
hyperbox<scalar_type> hyperbox<scalar_type>::empty_box(
		const positional_vdomain& dom) {
	return empty_box(dom.get_index_to_variable_id_map());
}

template<class scalar_type>
const typename hyperbox<scalar_type>::value_type& hyperbox<scalar_type>::get_l(
		size_type i) const {
	return my_l[i];
}

template<class scalar_type>
const typename hyperbox<scalar_type>::value_type& hyperbox<scalar_type>::get_u(
		size_type i) const {
	return my_u[i];
}

template<class scalar_type>
void hyperbox<scalar_type>::set_l(size_type i, value_type v) {
	my_l[i] = v;
}

template<class scalar_type>
void hyperbox<scalar_type>::set_u(size_type i, value_type v) {
	my_u[i] = v;
}

template<class scalar_type>
void hyperbox<scalar_type>::set_zero_dim_empty() {
	my_l = point_type();
	my_u = point_type();
	my_zero_dim_empty = true;
}

template<class scalar_type>
bool hyperbox<scalar_type>::check_consistency() {
	for (size_type i = 0; i < get_index_to_variable_id_map()->dimensions(); ++i) {
		if (get_l(i).is_pos_infinity() || get_u(i).is_neg_infinity()) {
			return false;
		}
	}
	return true;
}

template<class scalar_type>
hyperbox<scalar_type> hyperbox<scalar_type>::operator+(
		const hyperbox<scalar_type>& h) const {
	if (domain() == h.domain()) {
		return hyperbox<scalar_type> (get_l() + h.get_l(),
				get_u() + h.get_u(), domain());
	} else {
		using namespace math;
		position_map f1; position_map f2;
		positional_vdomain common_map = compose(domain(), h.domain(), f1, f2);
		vdom_vector<value_type> l1(domain(),get_l());
		vdom_vector<value_type> u1(domain(),get_u());
		vdom_vector<value_type> l2(h.domain(),h.get_l());
		vdom_vector<value_type> u2(h.domain(),h.get_u());

		l1.reorder(common_map);
		l2.reorder(common_map);
		u1.reorder(common_map);
		u2.reorder(common_map);

		vdom_vector<value_type> l = l1+l2;
		vdom_vector<value_type> u = u1+u2;

		return hyperbox<scalar_type> (l.get_vector(),u.get_vector(),common_map);
	}
}

template<class scalar_type>
hyperbox<scalar_type>& hyperbox<scalar_type>::operator+=(
		const hyperbox<scalar_type>& h) {
	if (domain() == h.domain()) {
		my_l += h.my_l;
		my_u += h.my_u;
	} else {
		*this = *this + h;
	}
	return *this;
}

template<class scalar_type>
void hyperbox<scalar_type>::reorder(const positional_vdomain& dom) {
	math::vdom_vector<value_type> new_l(dom,my_l);
	math::vdom_vector<value_type> new_u(dom,my_u);
	new_l.reorder(dom);
	new_u.reorder(dom);
	my_l=new_l.get_vector();
	my_u=new_u.get_vector();
	set_domain(dom);
}
;

}

#endif /* HYPERBOX_HPP_ */
