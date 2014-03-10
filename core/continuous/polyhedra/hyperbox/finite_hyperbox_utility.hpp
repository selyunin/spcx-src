/*
 * finite_hyperbox_utility.hpp
 *
 *  Created on: Jan 18, 2011
 *      Author: frehse
 */

#ifndef FINITE_HYPERBOX_UTILITY_HPP_
#define FINITE_HYPERBOX_UTILITY_HPP_

#include "finite_hyperbox_utility.h"

namespace continuous {

template<typename scalar_type> finite_hyperbox<scalar_type> finite_bounding_box(
		const support_function_provider& s) {
	scalar_type xp, xn;
	math::vdom_vector<scalar_type> v;
	bool is_empty = false;
	bool is_bounded = true;

	typename finite_hyperbox<scalar_type>::point_type c_point(s.get_dim());
	typename finite_hyperbox<scalar_type>::point_type g_point(s.get_dim());

	/* Positive and negative direction */
	variable_id_set vars = s.get_variable_ids();
	positional_vdomain dom(vars);
	for (unsigned int j = 0; j < dom.size(); ++j) {
		math::vdom_vector<scalar_type> ln(dom);
		math::vdom_vector<scalar_type> lp(dom);
		ln[j] = -scalar_type(1);
		lp[j] = scalar_type(1);
		s.compute_support(ln, xn, v, is_empty, is_bounded);
		if (!is_bounded) {
			std::stringstream ss;
			ss << "unbounded in direction " << ln << std::endl;
			ss << "in set: " << s << std::endl;
			throw std::runtime_error(
					"finite_bounding_box not allowed with unbounded set:\n"
							+ ss.str());
		};
		s.compute_support(lp, xp, v, is_empty, is_bounded);
		if (!is_bounded) {
			std::stringstream ss;
			ss << "unbounded in direction " << lp;
			ss << "in set: " << s << std::endl;
			throw std::runtime_error(
					"finite_bounding_box not allowed with unbounded set:\n"
							+ ss.str());
		};
		//std::cout << l << ":" << x << " at " << v << ", empty:" << is_empty << ", bounded:" << is_bounded << std::endl;
		if (is_empty) {
			// return an empty box
			return finite_hyperbox<scalar_type>::empty_box(dom);
		} else {
			if (is_bounded) {
				c_point[j] = (xp - xn) / scalar_type(2);
				g_point[j] = xp - c_point[j];
			} else {
				throw std::runtime_error(
						"finite_bounding_box not allowed with unbounded set");
			};
		}
	}
	return finite_hyperbox<scalar_type> (c_point, g_point, dom);
}

template<typename scalar_type> finite_hyperbox<scalar_type> finite_bounding_box(
		const support_function_provider& s,
		const math::affine_map<scalar_type>& t) {
	scalar_type xp, xn;
	math::vdom_vector<scalar_type> v;
	bool is_empty = false;
	bool is_bounded = true;

	positional_vdomain dom = t.domain();
	positional_vdomain codom = t.codomain();

	typename finite_hyperbox<scalar_type>::point_type c_point(codom.size());
	typename finite_hyperbox<scalar_type>::point_type g_point(codom.size());

	/* Positive and negative direction */
	for (unsigned int j = 0; j < codom.size(); ++j) {
		math::vdom_vector<scalar_type> lp(dom, t.get_A().vector_from_row(j));
		math::vdom_vector<scalar_type> ln(-lp);
		s.compute_support(ln, xn, v, is_empty, is_bounded);
		if (!is_bounded) {
			std::stringstream ss;
			ss << "unbounded in direction " << ln << std::endl;
			ss << "in set: " << s << std::endl;
			throw std::runtime_error(
					"finite_bounding_box not allowed with unbounded set:\n"
							+ ss.str());
		};
		s.compute_support(lp, xp, v, is_empty, is_bounded);
		if (!is_bounded) {
			std::stringstream ss;
			ss << "unbounded in direction " << lp;
			ss << "in set: " << s << std::endl;
			throw std::runtime_error(
					"finite_bounding_box not allowed with unbounded set:\n"
							+ ss.str());
		};
		xp += t.get_b()[j];
		xn -= t.get_b()[j];
		//std::cout << l << ":" << x << " at " << v << ", empty:" << is_empty << ", bounded:" << is_bounded << std::endl;
		if (is_empty) {
			// return an empty box
			return finite_hyperbox<scalar_type>::empty_box(codom);
		} else {
			if (is_bounded) {
				c_point[j] = (xp - xn) / scalar_type(2);
				g_point[j] = xp - c_point[j];
			} else {
				throw std::runtime_error(
						"finite_bounding_box not allowed with unbounded set");
			};
		}
	}
	return finite_hyperbox<scalar_type> (c_point, g_point, codom);
}

template<typename scalar_type> finite_hyperbox<scalar_type> finite_symmetric_bounding_box(
		const support_function_provider& s,
		const math::affine_map<scalar_type>& t) {
	scalar_type xp, xn;
	math::vdom_vector<scalar_type> v;
	bool is_empty = false;
	bool is_bounded = true;

	/* Positive and negative direction */
	positional_vdomain dom = t.domain();
	positional_vdomain codom = t.codomain();

	typename finite_hyperbox<scalar_type>::point_type c_point(codom.size(),
			scalar_type(0));
	typename finite_hyperbox<scalar_type>::point_type g_point(codom.size());

	for (unsigned int j = 0; j < codom.size(); ++j) {
		math::vdom_vector<scalar_type> lp(dom, t.get_A().vector_from_row(j));
		math::vdom_vector<scalar_type> ln(-lp);
		s.compute_support(ln, xn, v, is_empty, is_bounded);
		if (!is_bounded) {
			std::stringstream ss;
			ss << "unbounded in direction " << ln << std::endl;
			ss << "in set: " << s << std::endl;
			throw std::runtime_error(
					"finite_symmetric_bounding_box not allowed with unbounded set:\n"
							+ ss.str());
		};
		s.compute_support(lp, xp, v, is_empty, is_bounded);
		if (!is_bounded) {
			std::stringstream ss;
			ss << "unbounded in direction " << lp;
			ss << "in set: " << s << std::endl;
			throw std::runtime_error(
					"finite_symmetric_bounding_box not allowed with unbounded set:\n"
							+ ss.str());
		};
		xp += t.get_b()[j];
		xn -= t.get_b()[j];
		//std::cout << ln << ":" << xn << "..." << lp << ":" << xp << " at " << v << ", empty:" << is_empty << ", bounded:" << is_bounded << std::endl;
		if (is_empty) {
			// return an empty box
			return finite_hyperbox<scalar_type>::empty_box(codom);
		} else {
			// Use absolute value so that the generator is not negative
			// (which bef def it shouldn't be)
			if (xn > xp)
				g_point[j] = abs(xn);
			else
				g_point[j] = abs(xp);
		}
	}
	return finite_hyperbox<scalar_type> (c_point, g_point, codom);
}

template<typename scalar_type> hyperbox<scalar_type> construct_hyperbox(
		const finite_hyperbox<scalar_type>& box) {
	typename hyperbox<scalar_type>::finite_point_type l = box.get_c() - box.get_g();
	typename hyperbox<scalar_type>::finite_point_type u = box.get_c() + box.get_g();
	return hyperbox<scalar_type>(l, u, box.get_index_to_variable_id_map());
}

}

#endif /* FINITE_HYPERBOX_UTILITY_HPP_ */
