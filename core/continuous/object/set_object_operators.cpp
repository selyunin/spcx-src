/*
 * set_object_operators.cpp
 *
 *  Created on: Mar 25, 2010
 *      Author: frehse
 */

#include "set_object_operators.h"
#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"
#include "core/continuous/continuous_set_operators.h"
#include "core/continuous/continuous_set_operator_implementations/compute_transformation.h"
#include "core/continuous/support_function/sf_base/sf_set.h"
#include "core/continuous/support_function/sf_derived/sf_chull.h"
#include "core/continuous/support_function/sf_derived/sf_sum.h"
#include "core/continuous/support_function/template_directions/choose_directions.h"
#include "core/continuous/polyhedra/hyperbox/hyperbox.h"

namespace continuous {
namespace object {

void affine_transform_assign(set_object& X, const rational_matrix& A,
		const rational_vector& b) {
	if (sf_rational::ptr xsf = boost::dynamic_pointer_cast<sf_rational>(
			X.get_impl())) {
		math::affine_map<global_types::rational_type> M(A, b);
		xsf->map(M);
	} else {
		math::affine_map<global_types::rational_type> M(A, b);
		reset_affine_transform<global_types::rational_type> T(M);
		continuous_set_ptr y_ptr = compute_transformation(*X.get_impl(), T);
		X = set_object(y_ptr, X.get_number_type());
	}
}

void affine_transform_assign(set_object& X, const float_matrix& A,
		const float_vector& b) {
	if (sf_float::ptr xsf = boost::dynamic_pointer_cast<sf_float>(
			X.get_impl())) {
		math::affine_map<global_types::float_type> M(A, b);
		xsf->map(M);
	} else {
		math::affine_map<global_types::float_type> M(A, b);
		reset_affine_transform<global_types::float_type> T(M);
		continuous_set_ptr y_ptr = compute_transformation(*X.get_impl(), T);
		X = set_object(y_ptr, X.get_number_type());
	}
}

void existentially_quantify(set_object& X, const std::set<variable>& vars) {
	variable_id_set vis = create_variable_id_set(vars);
	X.get_impl()->existentially_quantify_variables(vis);
}

set_object intersection(const set_object& X, const set_object& Y) {
	continuous_set_ptr y_ptr = compute_intersection(X.get_impl(), Y.get_impl());
	return set_object(y_ptr, X.get_number_type());
}

template<template<typename > class sf_type>
set_object create_sf_binary(const set_object& X, const set_object& Y) {
	set_object Xs = support_function(X);
	set_object Ys = support_function(Y);
	sf_float::const_ptr xsf = boost::dynamic_pointer_cast<const sf_float>(
			Xs.get_impl());
	sf_float::const_ptr ysf = boost::dynamic_pointer_cast<const sf_float>(
			Ys.get_impl());
	sf_rational::const_ptr xsr =
			boost::dynamic_pointer_cast<const sf_rational>(Xs.get_impl());
	sf_rational::const_ptr ysr =
			boost::dynamic_pointer_cast<const sf_rational>(Ys.get_impl());
	if (xsf && ysf) {
		continuous_set_ptr s = continuous_set_ptr(new sf_type<
				global_types::float_type> (xsf, ysf));
		return set_object(s, X.get_number_type());
	} else if (xsr && ysr) {
		throw std::runtime_error("missing impl");
		return set_object();
	} else {
		throw std::runtime_error("object does not provide a support function");
		return set_object();
	}
}

set_object support_minkowski_sum(const set_object& X, const set_object& Y) {
	return create_sf_binary<support_function::sf_sum> (X, Y);
}

set_object support_convex_hull(const set_object& X, const set_object& Y) {
	return create_sf_binary<support_function::sf_chull> (X, Y);
}

template<typename T>
struct create_sf_set_implementor {
	static continuous_set_ptr implement(support_function_provider::const_ptr p) {
		return continuous_set_ptr(new support_function::sf_unary<T>(p));
	}
	;
};

continuous_set_ptr create_sf_set(support_function_provider::const_ptr p,
		global_types::coefficient_type t) {
	return global_types::coefficient_type_caller<continuous_set_ptr,
			support_function_provider::const_ptr, create_sf_set_implementor>::call(
			p, t);
}
;

set_object support_function(const set_object& X) {
	// Test whether it's already an sf_set, in which we do nothing
	if (boost::dynamic_pointer_cast<const support_function::sf_set<
			global_types::float_type> >(X.get_impl())) {
		set_object Y(X);
		return Y;
	} else if (boost::dynamic_pointer_cast<const support_function::sf_set<
			global_types::rational_type> >(X.get_impl())) {
		set_object Y(X);
		return Y;
	} else if (support_function_provider::const_ptr sprov=boost::dynamic_pointer_cast<const support_function_provider>(X.get_impl())) {
		continuous_set_ptr s = create_sf_set(sprov,
				global_types::coefficient_type(X.get_number_type()));
		return set_object(s, X.get_number_type());
	} else {
		throw std::runtime_error("object does not provide a support function");
		return set_object();
	}
}

set_object outer_poly(const set_object& X, const set_object& Y) {
	typedef support_function::sf_set<global_types::float_type> sf_float;
	typedef polyhedron<global_types::float_type> poly_float;
	typedef math::lin_constraint_system<global_types::float_type>
			lin_constraints_float;

	set_object Z = support_function(X);
	if (sf_float::const_ptr sf=boost::dynamic_pointer_cast<const sf_float>(Z.get_impl())) {
		sf_float::vector_set dirs;
		positional_domain dom=X.get_positional_domain();
		if (poly_float::const_ptr dir_poly=boost::dynamic_pointer_cast<const poly_float>(Y.get_impl())) {
			lin_constraints_float::const_ptr cons = dir_poly->get_constraints();
			for (lin_constraints_float::const_iterator it = cons->begin(); it
					!= cons->end(); ++it) {
				sf_float::vector_type v=it->get_normal();
				v.reorder(dom); // map to domain of X
				dirs.insert(v);
			}
			sf_float::poly_ptr p = sf->outer_poly(dirs);
			assert(p);
			set_object res(p, X.get_number_type());
			return res;
		} else {
			throw std::runtime_error("object does not provide directions");
			return set_object();
		}
	} else if (boost::dynamic_pointer_cast<const support_function::sf_set<
			global_types::rational_type> >(Z.get_impl())) {
		throw std::runtime_error("missing impl");
		return set_object();
	} else {
		throw std::runtime_error("object does not provide a support function");
		return set_object();
	}
}

set_object outer_poly(const set_object& X, const float_directions& dirs) {
	set_object Z = support_function(X);
	if (sf_float::const_ptr sf=boost::dynamic_pointer_cast<const sf_float>(Z.get_impl())) {
		sf_float::poly_ptr p = sf->outer_poly(dirs);
		assert(p);
		set_object res(p, X.get_number_type());
		return res;
	} else {
		throw std::runtime_error("object does not provide a support function");
		return set_object();
	}
}

float_directions uni_directions(const positional_domain& dom, unsigned int n) {
	float_directions dirs;
	typedef std::list<math::vector<global_types::float_type> > vdir_type;
	vdir_type vdirs;

	unsigned int dim = dom.size();
	support_function::add_uniform_directions<global_types::float_type>(dim, n,
			vdirs);

	for (vdir_type::const_iterator it = vdirs.begin(); it != vdirs.end(); ++it) {
		math::vdom_vector<global_types::float_type> d(dom, *it);
		dirs.insert(d);
	}

	return dirs;
}

float_directions box_directions(const positional_domain& dom) {
	float_directions dirs;
	typedef std::list<math::vector<global_types::float_type> > vdir_type;
	vdir_type vdirs;

	unsigned int dim = dom.size();
	support_function::add_box_directions(dim, vdirs);

	for (vdir_type::const_iterator it = vdirs.begin(); it != vdirs.end(); ++it) {
		math::vdom_vector<global_types::float_type> d(dom, *it);
		dirs.insert(d);
	}

	return dirs;
}

float_directions oct_directions(const positional_domain& dom) {
	float_directions dirs;
	typedef std::list<math::vector<global_types::float_type> > vdir_type;
	vdir_type vdirs;

	unsigned int dim = dom.size();
	support_function::add_octagonal_directions(dim, vdirs);

	for (vdir_type::const_iterator it = vdirs.begin(); it != vdirs.end(); ++it) {
		math::vdom_vector<global_types::float_type> d(dom, *it);
		dirs.insert(d);
	}

	return dirs;
}

}
}

std::ostream& operator<<(std::ostream& os,
		const continuous::object::float_directions& dirs) {
	os << "{";
	for (continuous::object::float_directions::const_iterator it = dirs.begin(); it
			!= dirs.end(); ++it) {
		if (it != dirs.begin())
			os << ",";
		os << (*it);
	}
	os << "}";

	return os;
}

