/*
 * constr_polyhedron.hpp
 *
 *  Created on: Dec 23, 2009
 *      Author: frehse
 */

#ifndef CONSTR_POLYHEDRON_HPP_
#define CONSTR_POLYHEDRON_HPP_

#include "constr_polyhedron.h"

#include "math/vdom/vdom_matrix_operators.h"

namespace continuous {

template<typename scalar_type>
constr_polyhedron<scalar_type>::~constr_polyhedron() {
}

template<typename scalar_type>
constr_polyhedron<scalar_type>* constr_polyhedron<scalar_type>::clone() const {
	// Create a new object of the same type
	constr_polyhedron<scalar_type>* ip = new my_type();

	// deep copy of member objects
	ip->my_poly = my_poly_ptr(new my_poly_type(*my_poly));
	return ip;
}

template<typename scalar_type>
constr_polyhedron<scalar_type>* constr_polyhedron<scalar_type>::create_universe() const {
	return new my_type();
}

template<typename scalar_type>
constr_polyhedron<scalar_type> constr_polyhedron<scalar_type>::empty_poly() {
	constr_polyhedron<scalar_type> res;
	res.add_constraint(math::lin_constraint<scalar_type>::unsatisfiable_constraint());
	return res;
}

template<typename scalar_type>
constr_polyhedron<scalar_type>* constr_polyhedron<scalar_type>::create_empty() const {
	constr_polyhedron<scalar_type>* p(new constr_polyhedron<scalar_type>(empty_poly()));
	return p;
}

template<typename scalar_type>
int constr_polyhedron<scalar_type>::get_memory() const {
	throw std::runtime_error("feature not implemented yet");
}

template<typename scalar_type>
const variable_id_set& constr_polyhedron<scalar_type>::get_variable_ids() const {
	return my_poly->get_variable_ids();
}

template<typename scalar_type>
dimension_t constr_polyhedron<scalar_type>::get_dim() const {
	/* get the variables of all the constraints and count them. */
	return get_variable_ids().size();
}

template<typename scalar_type>
bool constr_polyhedron<scalar_type>::is_trivially_empty() const {
	typename my_poly_type::const_iterator it = my_poly->begin();
	typename my_poly_type::const_iterator en = my_poly->end();
	for (; it != en; ++it) {
		if (!math::maybe(it->is_satisfiable())) {
			return true;
		}
	}
	return false;
}

template<typename scalar_type>
math::tribool constr_polyhedron<scalar_type>::contains(const point_type& x) const {
	return my_poly->is_satisfied(x);
}

template<typename scalar_type>
math::tribool constr_polyhedron<scalar_type>::is_empty() const {
	if (my_poly->size() == 0)
		return false;
	else {
		if (is_trivially_empty()) {
			// std::cout << "triv";
			return true;
		} else {
			//my_up_to_date=false;
			math::lin_expression<scalar_type> f = my_poly->begin()->get_l();
			typename math::lp_solver<scalar_type>::lp_result res;
			maximize(f, res);
//			std::cout << "constraints: " << my_poly << std::endl;
//			std::cout << "empty "<< global_types::type_identifier<scalar_type>::name << ": " << f << " unsat:" << res.is_unsat << " opt: " << res.max_value << " vec: " << res.support_vec <<  std::endl;
			return res.is_unsat;

			// check whether the result violates any constraints (in the tribool sense)
			//return this->contains(res.support_vec);
		}
	}
}

template<typename scalar_type>
math::tribool constr_polyhedron<scalar_type>::is_universe() const {
	math::tribool result(true);
	for (typename my_poly_type::const_iterator it = my_poly->begin(); it
			!= my_poly->end(); ++it) {
		result = result && it->is_always_satisfied();
	}
	return result;
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::simplify() {
	typename math::lp_solver<scalar_type>::lp_result res;
	math::lin_constraint_system<scalar_type> &the_constrs = *my_poly;

	if (is_empty()) {
		my_poly=my_poly_ptr(new math::lin_constraint_system<scalar_type>());
		my_poly->insert(math::lin_constraint<scalar_type>::unsatisfiable_constraint());
	}
	else {

		typename math::lin_constraint_system<scalar_type>::iterator it =
			the_constrs.begin();

	while (it != the_constrs.end()) {
		math::lin_expression<scalar_type> f = it->get_l();
		scalar_type b=f.get_inh_coeff();
		f.set_inh_coeff(scalar_type(0));
		if (is_LT_or_LE(it->get_sign())) {
			maximize(f, res);
			// if (res.max_value < scalar_type(0))
			if (res.is_bounded && math::numeric::approx_comparator<scalar_type>::is_definitely_strictly_smaller(res.max_value,-b)) {
				/* redundant constraint removed */
//std::cout << "removing " << f << " <= " << res.max_value << " < " << -b << std::endl;
				it = the_constrs.erase(it);
				my_up_to_date = false;
			} else {
				++it;
			}
		} else if (is_GT_or_GE(it->get_sign())) {
			maximize(f * scalar_type(-1), res);
			// if (scalar_type(0) < scalar_type(-1) * res.max_value) {
			if (res.is_bounded && math::numeric::approx_comparator<scalar_type>::is_definitely_strictly_smaller(res.max_value,b)) {
				/* redundant constraint removed */
//std::cout << "removing " << f << " >= " << -res.max_value << " > " << -b << std::endl;
				it = the_constrs.erase(it);
				my_up_to_date = false;
			} else {
				++it;
			}
		} else {
			++it;
		}
	} }
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::collapse_inequalities() {
	my_poly->collapse_inequalities();
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::remove_redundant_constraints() {
	LOGGERSWOC(DEBUG4,__FUNCTION__,"Removing redundant polyhedron constraints");

	using math::definitely;
	using math::maybe;

	math::lin_constraint<scalar_type> con;
	math::lin_expression<scalar_type> f;
	comparison_operator sign;

	math::lin_constraint_system<scalar_type> &the_constrs = *my_poly;
	int size = the_constrs.size();
	typename math::lin_constraint_system<scalar_type>::iterator it =
			the_constrs.begin();
	int ttl = 0;
	while (ttl < size) {
		con = *it;
		sign = con.get_sign();
		//LOGGER_OS(DEBUG6,__FUNCTION__) << "checking " << con;
		it = the_constrs.erase(it);
		//LOGGER_OS(DEBUG6,__FUNCTION__) << "   next  " << *it;
		my_up_to_date = false;
		if (sign == EQ) {
			con.set_sign(LE);
			if (maybe(is_redundant(con))) {
				con.set_sign(GE);
				if (maybe(is_redundant(con))) {
					/* forget the constraint. */
				} else {
					/* LEQ redundant, GEQ not redundant. */
					the_constrs.push_back(con);
				}
			} else {
				con.set_sign(GE);
				if (maybe(is_redundant(con))) {
					/* LEQ not redundant, GEQ redundant */
					con.set_sign(LE);
					the_constrs.push_back(con);
				} else {
					/* neither ineq redundant so keep eq */
					con.set_sign(EQ);
					the_constrs.push_back(con);
				}
			}
		} else {
			/* it's not an equality, only one case to consider. */
			if (!maybe(is_redundant(con))) {
				the_constrs.push_back(con);
			} else {
				/* forget the constraint */
			}
		}
		++ttl;
	}
	my_up_to_date = false;
}

template<typename scalar_type>
math::tribool constr_polyhedron<scalar_type>::is_redundant(const math::lin_constraint<
		scalar_type> &con) const {
	// @todo take care of strict inequalities
	LOGGERSWOC(DEBUG6,__FUNCTION__,"Redundancy testing for single constraint");

	if (!is_EQ(con.get_sign())) {
		typename math::lp_solver<scalar_type>::lp_result res;
		//		typename math::lin_constraint<scalar_type>::sign sign = con.get_sign();

		// compute max w/o inh. coeff so the comparison is not against
		// zero (better relative error)
		math::lin_expression<scalar_type> l;
		l = math::lin_expression<scalar_type> (con.get_normal());
		maximize(l, res);
		if (res.is_unsat) {
			return true;
		}
		if (!res.is_bounded) {
			return false;
		}
		scalar_type x = res.max_value;
		scalar_type y = -con.get_canonic_inh_coeff();

		// the constraint is redundant if x <= y
		//LOGGER_OS(DEBUG6,__FUNCTION__) << con << " test: " << x << " is_LT " << y;
		return math::numeric::is_LE(x, y);
	} else {
		//std::cout << "EQ REDUNDANCY CHECK" << std::endl;
		math::lin_constraint<scalar_type> con2=con;
		con2.set_sign(LE);
		math::tribool red2a = is_redundant(con2);
		if (math::definitely(!red2a)) {
			return false;
		} else {
			con2.set_sign(GE);
			math::tribool red2b = is_redundant(con2);
			return red2a && red2b;
		}
	}
}

template<typename scalar_type> bool constr_polyhedron<scalar_type>::computes_support_vector() const {
	return true;
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::compute_support(const math::vdom_vector<
		Rational> &l, Rational & max_value,
		math::vdom_vector<Rational> &support_vec, bool & is_empty,
		bool & is_bounded) const {
	math::lin_expression<scalar_type> f(l.convert_to<scalar_type> ());
	typename math::lp_solver<scalar_type>::lp_result res;
	maximize(f, res);
	max_value = Rational(res.max_value);
	support_vec = res.support_vec.template convert_to<Rational> ();
	is_empty = res.is_unsat;
	is_bounded = res.is_bounded;
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::compute_support(const math::vdom_vector<
		double>&l, double &max_value, math::vdom_vector<double>&support_vec,
		bool & is_empty, bool & is_bounded) const {
	math::lin_expression<scalar_type> f(l.convert_to<scalar_type> ());
	typename math::lp_solver<scalar_type>::lp_result res;
	maximize(f, res);
	max_value = convert_element<double> (res.max_value);
	support_vec = res.support_vec.template convert_to<double> ();
	is_empty = res.is_unsat;
	is_bounded = res.is_bounded;
}

template<typename scalar_type>
typename math::lin_constraint_system<scalar_type>::const_ptr constr_polyhedron<
		scalar_type>::get_constraints() const {
	return my_poly;
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::embed_variables(
		const variable_id_set & id_set) {
	// don't need to do anything as per
	// code from old constr_polyhedron_implementor.h
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::existentially_quantify_variables(
		const variable_id_set &id_set) {
	math::lin_constraint_system<scalar_type> &the_constrs = *my_poly;
	variable_id_set::const_iterator it;
	for (it = id_set.begin(); it != id_set.end(); ++it) {
		the_constrs = math::fourier_motzkin::eliminate<scalar_type>(
				the_constrs, *it);
	}
	/* Remove id from the domain of the constraints. */
	the_constrs.remove_variables(id_set);
	set_poly(the_constrs);
	my_up_to_date = false;
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::swap(my_type& s) {
	my_poly.swap(s.my_poly);
	//math::lp_solver_user<scalar_type>::swap(static_cast<math::lp_solver_user<scalar_type> >(s));
	math::lp_solver_user<scalar_type>::swap(s);
	std::swap(my_up_to_date, s.my_up_to_date);
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::intersection_assign(const polyhedron<
		scalar_type> &p) {
	add_constraints(*p.get_constraints());
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::add_constraint(const math::lin_constraint<
		scalar_type> &c, bool check_redundancy) {
	if (!check_redundancy || math::maybe(!is_redundant(c))) {
		my_poly->push_back(c);
		my_up_to_date = false;
	}
}

template<typename scalar_type>
variable_id_set constr_polyhedron<scalar_type>::get_primed_variables(
		unsigned int prime_count) const {

	variable_id_set vis, vtmp;
	typename my_poly_type::const_iterator i;
	for (i = my_poly->begin(); i != my_poly->end(); ++i) {
		vtmp = i->get_primed_variables(prime_count);
		vis.insert(vtmp.begin(), vtmp.end());
	}
	return vis;
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::reassign_primedness(unsigned int d,
		unsigned int p) {

	typename my_poly_type::iterator i;
	for (i = my_poly->begin(); i != my_poly->end(); ++i) {
		i->reassign_primedness(d, p);
	}
	my_up_to_date = false;
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::increase_primedness(unsigned int d) {

	typename my_poly_type::iterator i = my_poly->begin();
	for (; i != my_poly->end(); ++i) {
		i->increase_primedness(d);
	}
	my_up_to_date = false;
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::decrease_primedness(unsigned int d) {

	typename my_poly_type::iterator i;
	for (i = my_poly->begin(); i != my_poly->end(); ++i) {
		i->decrease_primedness(d);
	}
	my_up_to_date = false;
}

template<typename scalar_type>
typename constr_polyhedron<scalar_type>::my_poly_type& constr_polyhedron<
		scalar_type>::get_poly() const {
	return *my_poly;
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::set_poly(const my_poly_type &new_poly) {
	my_poly = my_poly_ptr(new my_poly_type(new_poly));
	my_up_to_date = false;
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::set_poly(my_poly_ptr new_poly) {
	my_poly = new_poly;
	my_up_to_date = false;
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::accept(dispatching::dispatcher<
		continuous_set_typelist> &d) const {
	d.dispatch(this);
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::print_double_generators(std::ostream& os) const {
	constr_polyhedron<scalar_type>* p = const_cast<constr_polyhedron<
			scalar_type>*> (this);
	p->remove_redundant_constraints();
	polyhedron<scalar_type>::print_double_generators(os);
}


template<typename scalar_type>
void constr_polyhedron<scalar_type>::print_JVX(std::ostream& os) const {
	constr_polyhedron<scalar_type>* p = const_cast<constr_polyhedron<
			scalar_type>*> (this);
	p->remove_redundant_constraints();
	polyhedron<scalar_type>::print_JVX(os);
}

template<typename scalar_type>
void constr_polyhedron<scalar_type>::maximize(
		const math::lin_expression<scalar_type> &f,
		typename math::lp_solver<scalar_type>::lp_result &res) const {
	assert(my_poly);
	constr_polyhedron<scalar_type>* nonconst_this=const_cast<constr_polyhedron<scalar_type>*>(this);
	math::lp_solver<scalar_type>& solver = nonconst_this->get_lp_solver();
	if (!my_up_to_date || nonconst_this->lp_solver_has_changed()) {
//		std::cout << "setting constraints " << *my_poly << std::endl;
		solver.set_constraints(*my_poly);
		nonconst_this->my_up_to_date = true;
	}
	solver.maximize(f, res);
//	std::cout << "max " << f << " in " << *my_poly << " is " << res.is_unsat << " solver " << nonconst_this->get_lp_solver()->get_solver_name() << std::endl;
}

template<typename scalar_type>
constr_polyhedron<scalar_type>::constr_polyhedron() {
	my_poly = my_poly_ptr(new my_poly_type());
	my_up_to_date = false;
}

template<typename scalar_type>
constr_polyhedron<scalar_type>::constr_polyhedron(math::lp_solver<scalar_type> *s) {
	my_poly = my_poly_ptr(new my_poly_type());
	my_up_to_date = false;
	set_lp_solver(s);
}

template<typename scalar_type>
constr_polyhedron<scalar_type>::constr_polyhedron(my_poly_ptr orig_poly) {
	my_poly = orig_poly;
	my_up_to_date = false;
}

template<typename scalar_type>
constr_polyhedron<scalar_type>::constr_polyhedron(const constr_polyhedron<
		scalar_type> &orig_poly) : math::lp_solver_user<scalar_type>(orig_poly) {
	// make a deep copy of the poly
	my_poly = my_poly_ptr(new my_poly_type(orig_poly.get_poly()));
	my_up_to_date = false;
}

template<typename scalar_type>
constr_polyhedron<scalar_type>&
constr_polyhedron<scalar_type>::operator=(const my_type &orig_poly) {
	math::lp_solver_user<scalar_type>::operator=(orig_poly);
	// make a deep copy of the poly
	my_poly = my_poly_ptr(new my_poly_type(orig_poly.get_poly()));
	my_up_to_date = false;
	return *this;
}


} // namespace continuous

#endif /* CONSTR_POLYHEDRON_HPP_ */
