/*
 * bounding_box.hpp
 *
 *  Created on: Nov 27, 2009
 *      Author: frehse
 */

#ifndef BOUNDING_BOX_HPP_
#define BOUNDING_BOX_HPP_

#include "bounding_box.h"

namespace continuous {

template<typename scalar_type>
bounding_box<scalar_type>::bounding_box(
		const support_function_provider::const_ptr& s) :
	hyperbox<scalar_type> (index_to_variable_id_map::get_map_with_ids(
			s->get_variable_ids())), my_root(s), my_l_is_uptodate(false),
			my_u_is_uptodate(false) {
	if (s->get_dim() == 0 && s->is_empty())
		this->set_zero_dim_empty();
}

template<typename scalar_type>
bounding_box<scalar_type>::~bounding_box() {
}

template<typename scalar_type>
typename bounding_box<scalar_type>::value_type bounding_box<scalar_type>::compute_sup(
		size_type j, bool pos) const {
	scalar_type x;
	math::vdom_vector<scalar_type> v;
	math::vdom_vector<scalar_type> l(
			this->get_index_to_variable_id_map());
	if (pos)
		l[j] = scalar_type(1);
	else
		l[j] = scalar_type(-1);
	bool is_empty, is_bounded;
	my_root->compute_support(l, x, v, is_empty, is_bounded);
	//std::cout << l << ":" << x << " at " << v << ", empty:" << is_empty << ", bounded:" << is_bounded << std::endl;
	if (is_empty) {
		return value_type::NaN();
	} else {
		if (is_bounded) {
			if (pos)
				return value_type(x);
			else
				return value_type(-x);
		} else if (pos)
			return value_type::pos_infty();
		else
			return value_type::neg_infty();
	}
}

template<typename scalar_type>
const typename bounding_box<scalar_type>::point_type& bounding_box<scalar_type>::get_l() const {
	update_l();
	return hyperbox<scalar_type>::get_l();
}

template<typename scalar_type>
const typename bounding_box<scalar_type>::point_type& bounding_box<scalar_type>::get_u() const {
	update_u();
	return hyperbox<scalar_type>::get_u();
}

template<typename scalar_type>
const typename bounding_box<scalar_type>::value_type& bounding_box<scalar_type>::get_l(
		size_type i) const {
	if (!my_l_is_uptodate && hyperbox<scalar_type>::get_l(i).is_NaN()) {
		update_l(i);
	}
	return hyperbox<scalar_type>::get_l(i);
}

template<typename scalar_type>
const typename bounding_box<scalar_type>::value_type& bounding_box<scalar_type>::get_u(
		size_type i) const {
	if (!my_u_is_uptodate && hyperbox<scalar_type>::get_u(i).is_NaN()) {
		update_u(i);
	}
	return hyperbox<scalar_type>::get_u(i);
}

template<typename scalar_type>
void bounding_box<scalar_type>::update_l(size_type i) const {
	bounding_box<scalar_type>* h =
			const_cast<bounding_box<scalar_type>*> (this);
	h->set_l(i, compute_sup(i, false));
	//std::cerr << "updating " << i <<" " << compute_sup(i, false);
}

template<typename scalar_type>
void bounding_box<scalar_type>::update_u(size_type i) const {
	bounding_box<scalar_type>* h =
			const_cast<bounding_box<scalar_type>*> (this);
	h->set_u(i, compute_sup(i, true));
}

template<typename scalar_type>
void bounding_box<scalar_type>::update_l() const {
	if (!my_l_is_uptodate) {
		for (size_type i = 0; i < this->get_dim(); ++i) {
			update_l(i);
		}
		bounding_box<scalar_type>* h =
				const_cast<bounding_box<scalar_type>*> (this);
		h->my_l_is_uptodate = true;
	}
}

template<typename scalar_type>
void bounding_box<scalar_type>::update_u() const {
	if (!my_u_is_uptodate) {
		for (size_type i = 0; i < this->get_dim(); ++i) {
			update_u(i);
		}
		bounding_box<scalar_type>* h =
				const_cast<bounding_box<scalar_type>*> (this);
		h->my_u_is_uptodate = true;
	}
}

template<typename scalar_type>
hyperbox<scalar_type> compute_bounding_box(const support_function_provider& s) {
	typedef typename hyperbox<scalar_type>::value_type value_type;

	bool is_empty = false;
	bool is_bounded = true;

	/* Positive and negative direction */
	variable_id_set vars = s.get_variable_ids();
	positional_vdomain dom(vars);
	typename hyperbox<scalar_type>::point_type l_point(dom.size());
	typename hyperbox<scalar_type>::point_type u_point(dom.size());

	for (int i = 0; i < 2; ++i) {
		for (unsigned int j = 0; j < dom.size(); ++j) {
			typedef typename support_type<scalar_type>::type s_type;
			math::vdom_vector<s_type> l(dom);
			l[j] = s_type(-1 + 2 * i); // first time with -1, second time with +1

			s_type x;
			math::vdom_vector<s_type> v;
			s.compute_support(l, x, v, is_empty, is_bounded);
			//std::cout << l << ":" << x << " at " << v << ", empty:" << is_empty << ", bounded:" << is_bounded << std::endl;
			if (is_empty) {
				// return an empty box
				return hyperbox<scalar_type>::empty_box(dom);
			} else {
				if (is_bounded) {
					if (i == 0) { // direction -1
						l_point[j] = typename hyperbox<scalar_type>::value_type(-x);
					} else { // direction 1
						u_point[j] = typename hyperbox<scalar_type>::value_type(x);
					}
				} else if (i == 0) { // direction -1
					l_point[j] = value_type::neg_infty();
				} else { // direction 1
					u_point[j] = value_type::pos_infty();
				};
			}
		}
	}
	return hyperbox<scalar_type> (l_point, u_point, dom);
}

template<typename scalar_type>
hyperbox<scalar_type> compute_bounding_box(const support_function_provider& s,
		const math::affine_map<scalar_type>& t) {
	typedef typename hyperbox<scalar_type>::value_type value_type;

	bool is_empty = false;
	bool is_bounded = true;

	positional_vdomain dom = t.domain();
	positional_vdomain codom = t.codomain();

	typename hyperbox<scalar_type>::point_type l_point(codom.size());
	typename hyperbox<scalar_type>::point_type u_point(codom.size());

	/* Positive and negative direction */
	math::vdom_vector<scalar_type> l(dom);

	for (unsigned int j = 0; j < codom.size(); ++j) {
		typedef typename support_type<scalar_type>::type s_type;
		math::vector<scalar_type> dir = t.get_A().vector_from_row(j);
		if (!dir.is_zero()) {
			for (int i = 0; i < 2; ++i) {
				//l[j] = s_type(-1 + 2 * i);
				// first time with -1, second time with +1
				if (i == 0) {
					l.set_vector(-dir);
				} else {
					l.set_vector(dir);
				}

				s_type x;
				math::vdom_vector<s_type> v;
				s.compute_support(l, x, v, is_empty, is_bounded);

				if (i == 1) {
					x += t.get_b()[j];
				} else {
					x -= t.get_b()[j];
				}

				//std::cout << l << ":" << x << " at " << v << ", empty:" << is_empty << ", bounded:" << is_bounded << std::endl;
				if (is_empty) {
					// return an empty box
					return hyperbox<scalar_type>::empty_box(codom);
				} else {
					if (is_bounded) {
						if (i == 0) { // direction -1
							l_point[j] = value_type(-x);
						} else { // direction 1
							u_point[j] = value_type(x);
						}
					} else if (i == 0) { // direction -1
						l_point[j] = value_type::neg_infty();
					} else { // direction 1
						u_point[j] = value_type::pos_infty();
					};
				}
			}
		} else {
			// There's only the b component
			l_point[j] = value_type(t.get_b()[j]);
			u_point[j] = value_type(t.get_b()[j]);
		}
	}
	return hyperbox<scalar_type> (l_point, u_point, codom);
}

}

#endif /* BOUNDING_BOX_HPP_ */
