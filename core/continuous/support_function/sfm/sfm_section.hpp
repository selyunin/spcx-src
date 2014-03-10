/*
 * sfm_section.hpp
 *
 *  Created on: May 3, 2010
 *      Author: frehse
 */

#ifndef SFM_SECTION_HPP_
#define SFM_SECTION_HPP_

#include "sfm_section.h"

#include "math/vdom/vdom_matrix_operators.h"
#include "math/vdom/vdom_vector_operators.h"
#include "math/vdom/affine_map_utility.h"
#include "math/vdom/positional_vdomain.h"

namespace continuous {
namespace support_function {


/** Declaration of static member variables. */
template<class scalar_type> bool sfm_section<scalar_type>::compact_sfms = true;

template<typename scalar_type> sfm_section<scalar_type>::sfm_section(
		const sfm_ptr& s, index_interval intv, bool compact_this) :
	my_set(s), my_intv(intv) {
	if (compact_this && compact_sfms)
		compact();
}
;

template<typename scalar_type> sfm_section<scalar_type>::sfm_section(
		const sfm_ptr& s, index_interval intv, const affine_map& M) :
	sf_set<scalar_type> (M), my_set(s), my_intv(intv) {
	if (compact_sfms)
		compact();
}
;

template<typename scalar_type> sfm_section<scalar_type>::~sfm_section() {
}
;

template<typename scalar_type> sfm_section<scalar_type>* sfm_section<
		scalar_type>::clone() const {
	sfm_ptr new_root(my_set->clone());
	if (this->get_map())
		return new sfm_section<scalar_type> (new_root, my_intv,
				*this->get_map());
	else
		return new sfm_section<scalar_type> (new_root, my_intv);
};

template<typename scalar_type> void sfm_section<scalar_type>::compact() {
	// check if the dynamics are deterministic, and if yes,
	// redefine X0
	LOGGER(DEBUG7,__FUNCTION__,"compacting sfm interval of size "+to_string(my_intv.size()+1));

	// only necessary if the interval doesn't start at the beginning
	if (my_intv.lower().get_val() > 0) {

	// get the dynamics from the sfm
	postc_params<scalar_type> params = my_set->get_postc_prb();
	const postc_params<scalar_type>& old_params = my_set->get_postc_prb();

	support_function_provider::const_ptr U = boost::static_pointer_cast<const support_function_provider>(params.input_set_ptr);
	support_function_provider::const_ptr U_centered;

	positional_vdomain dom = my_set->domain();
	math::vdom_vector<scalar_type> b,u0_shifted;
	math::vdom_vector<scalar_type> old_b = math::vdom_vector<scalar_type>(dom,params.dynamics_b);
	if (!U || is_input_set_point(U,U_centered,b)) {
		if (!U){
			// b = 0
			b = old_b;
		} else {
			b += old_b;
		}

		const math::matrix<scalar_type>& A = params.dynamics_A;

		unsigned int i_first = my_intv.lower().get_val();
		unsigned int i_last = my_intv.upper().get_val();
		unsigned int i_length = i_last - i_first + 1;

		// compute new start time
		scalar_type t_start = 0;
		if (!old_params.delta_vec.empty()) {
			for (unsigned int i = 0; i< i_first; ++i) {
				t_start += params.delta_vec[i];
			}
		} else {
			t_start = scalar_type(i_first) * scalar_type(params.delta);
		}

		// get exponential matrices for t_start
		math::matrix<scalar_type> Phi,Phi1,Phi2;
		math::get_special_matrices(A, t_start, Phi, Phi1, Phi2);
		// Phi = e^(A t_start)
		// new initial set is Phi X + Phi1 b
		u0_shifted = math::vdom_vector<scalar_type>(dom,Phi1 * b.get_vector());
		// define new initial set
		support_function_provider::const_ptr old_initial = boost::static_pointer_cast<const support_function_provider>(params.initial_set_ptr);
		params.initial_set_ptr = support_function_provider::const_ptr(new sf_unary<scalar_type> (
				old_initial, affine_map(math::vdom_matrix<scalar_type>(dom,dom,Phi),u0_shifted)));

		// cut down matrices and vectors of the sfm

		// cut down new time horizon
		params.time_horizon -= t_start;
		// cut down the new delta_vec
		if (!params.delta_vec.empty()) {
			params.delta_vec = std::vector<scalar_type>(old_params.delta_vec.begin()+i_first,old_params.delta_vec.begin()+i_last+1);
		}
		// redefine the new sfm matrix
		typedef typename sfm_cont_set<scalar_type>::matrix_type matrix_type;
		const matrix_type& old_sfm = my_set->get_sfm();
		unsigned int row_size = old_sfm.size1();
		matrix_type new_sfm = old_sfm.project_submatrix(0,row_size,i_first,i_last+1);

		typename sfm_cont_set<scalar_type>::ptr old_set = my_set;

		// create the new sfm
		my_set = sfm_ptr(new sfm_cont_set<scalar_type>(params,
				new_sfm,
				old_set->get_directions(),
				old_set->get_index_to_variable_id_map(),
				old_set->my_nonstate_contraints));

		LOGGER_OS(DEBUG7,__FUNCTION__)
			<< "compacting indices " << i_first << " to " << i_last
					<< std::endl << "old sfm: " << old_set->get_sfm() << std::endl
					<< "new sfm: " << my_set->get_sfm() << std::endl
					<< "old directions: " << my_set->get_directions().size();

		// update the interval
		my_intv = index_interval(0,i_last-i_first);

//		continuous::polyhedron<scalar_type>::set_output_format(continuous::polyhedron<scalar_type>::DOUBLE_GENERATORS);
//		math::vector<scalar_type> test_dir = (old_set->get_directions().begin()->first)+((++old_set->get_directions().begin()))->first;
//		old_set->extend_sfm(test_dir);
//		my_set->extend_sfm(test_dir);
//		std::cout << "test direction " << test_dir << std::endl;



//		std::ofstream myfile;
//		myfile.open ("/tmp/out_tmp_old.tmp");
//		for(unsigned int i=i_first;i<=i_last;i++)
//		{
//			continuous::constr_polyhedron<scalar_type> poly = old_set->get_polytope(i);
//			myfile << poly; myfile << std::endl << std::endl;
//		}
//		myfile.close();
//		myfile.open("/tmp/out_tmp.tmp");
//		for (unsigned int i = 0; i < my_set->get_sfm().size2(); i++) {
//			continuous::constr_polyhedron<scalar_type> poly =
//					my_set->get_polytope(i);
//			myfile << poly;
//			myfile << std::endl << std::endl;
//		}
//		myfile.close();
//		int res = system("graph -TX -C -B -q0.5 /tmp/out_tmp_old.tmp -s -m 2 -q-1 /tmp/out_tmp.tmp");
	}
	}
};


template<typename scalar_type> int sfm_section<scalar_type>::get_memory() const {
	throw std::runtime_error("sfm_section : missing implementation get_memory");
	return 0;
}

template<typename scalar_type> continuous_set_predicate::ptr sfm_section<
		scalar_type>::get_predicate() const {
	throw std::runtime_error(
			"sfm_section : missing implementation get_predicate");
	return continuous_set_predicate::ptr();
}

template<typename scalar_type> math::tribool sfm_section<scalar_type>::is_empty() const {
	return my_intv.is_empty() || my_set->is_empty();
}

template<typename scalar_type>
math::tribool sfm_section<scalar_type>::is_universe() const {
	return my_set->is_universe();
}

template<typename scalar_type>
typename sfm_section<scalar_type>::index_interval sfm_section<scalar_type>::get_interval() const {
	return my_intv;
}

template<typename scalar_type>
typename sfm_section<scalar_type>::sfm_ptr sfm_section<scalar_type>::get_sfm_set() {
	return my_set;
}

template<typename scalar_type> void sfm_section<scalar_type>::print(
		std::ostream& os) const {
	throw std::runtime_error("sfm_section : missing implementation print");
}

template<typename scalar_type> const variable_id_set& sfm_section<scalar_type>::get_variable_ids() const {
	if (!this->get_map() || this->get_map()->is_empty())
		return my_set->get_variable_ids();
	else
		return this->get_map()->codomain().get_variable_ids();
}
template<typename scalar_type> void sfm_section<scalar_type>::reassign_primedness(
		unsigned int, unsigned int) {
	throw std::runtime_error(
			"sfm_section : missing implementation reassign_primedness");
}
template<typename scalar_type> void sfm_section<scalar_type>::increase_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sfm_section : missing implementation increase_primedness");
}
template<typename scalar_type> void sfm_section<scalar_type>::decrease_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sfm_section : missing implementation decrease_primedness");
}

template<typename scalar_type> bool sfm_section<scalar_type>::computes_support_vector() const {
	return my_set->computes_support_vector();
}
;

template<typename scalar_type>
template<typename fun_type> void sfm_section<scalar_type>::compute_support_interval(
		const sfm_ptr& a_set, index_interval intv, const math::vdom_vector<
				fun_type>& l, fun_type& max_value,
		math::vdom_vector<fun_type>& support_vec, bool& is_empty,
		bool& is_bounded) {

	// GF 2012-12-05 removed emptiness check because at this point it is called
	// way too often!
	if (false && a_set->is_empty()) {
		std::cout << "sfm_section compute: empty reached\n";
		//std::cout << "sfm is:\n" << a_set << std::endl;

		is_empty = true;
		is_bounded = true;
	} else {
		// @todo do proper emptyness check
		is_empty = false;
		is_bounded = true;

		std::pair<unsigned int, scalar_type> pr;
		// @todo if l has variables not in the domain of *this, the problem is unbounded
		// The current version will throw.
		if (a_set->domain()==l.domain()) {
			pr = a_set->extend_sfm(l.get_vector());
		} else {
			math::vdom_vector<scalar_type> v = l;
			v.reorder(a_set->domain());
			typename sfm_cont_set<scalar_type>::direction d = v.get_vector();
			pr = a_set->extend_sfm(d);
		}

		unsigned int row_index = pr.first;
		scalar_type f = pr.second;

		const typename sfm_cont_set<scalar_type>::matrix_type& sfm =
				a_set->get_sfm();

		unsigned int idx1 = 0;
		unsigned int idx2 = sfm.size2();

		if (intv.lower().is_finite()) {
			idx1 = intv.lower().get_val();
		}
		if (intv.upper().is_finite()) {
			idx2 = intv.upper().get_val() + 1;
		}

		scalar_type x = sfm(row_index, idx1);
		for (unsigned int j = idx1 + 1; j < idx2; j++) {
			x = std::max(x, sfm(row_index, j));
		}
		max_value = x / f;

//std::cout << "extended sfm with direction " << v << ", obtained max " << max_value << " with scaling factor " << f << std::endl;
	}
}

template<typename scalar_type>
template<typename fun_type> void sfm_section<scalar_type>::compute_support_mapped(
		const sfm_ptr& a_set, index_interval intv,
		const affine_map_const_ptr& a_map,
		const math::vdom_vector<fun_type>& v, fun_type& max_value,
		math::vdom_vector<fun_type>& support_vec, bool& is_empty,
		bool& is_bounded) {
	if (a_map) {
		math::affine_map<fun_type> M = a_map->template convert_to<fun_type> ();
		// transform the cost function
		math::vdom_vector<fun_type> vmapped =
				math::vdom_vector<fun_type>(v * (M.get_A()));
		compute_support_interval(a_set, intv, vmapped, max_value, support_vec,
				is_empty, is_bounded);
		max_value += scalar_product(v, M.get_b());
		// transform the support vector if there is one
		if (support_vec.size() > 0) {
			support_vec = M.map(support_vec);
		}
	} else {
		compute_support_interval(a_set, intv, v, max_value, support_vec,
				is_empty, is_bounded);
	}
}

template<typename scalar_type> void sfm_section<scalar_type>::compute_support(
		const math::vdom_vector<Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	/*
	this->template compute_support_mapped<Rational>(my_set,
			my_intv, this->get_map(), l, max_value, support_vec, is_empty,
			is_bounded);
			*/
	throw std::runtime_error(
			"sfm_section : missing implementation compute_support for Rational");
}

template<typename scalar_type> void sfm_section<scalar_type>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	this->template compute_support_mapped<double> (my_set, my_intv,
			this->get_map(), l, max_value, support_vec, is_empty, is_bounded);
}

} // end of namespace support_function
} // end of namespace continuous

#endif /* SFM_SECTION_HPP_ */
