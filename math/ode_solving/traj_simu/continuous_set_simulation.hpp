#ifndef continuous_set_SIMULATION_HPP_
#define continuous_set_SIMULATION_HPP_

#include "continuous_set_simulation.h"

#include "math/numeric/approx_comparator.h"
#include "math/vdom/vdom_vector_utility.h"
#include "math/vdom/vdom_vector.h"

namespace continuous {

template<typename scalar_type>
 continuous_set_simulation<scalar_type>* continuous_set_simulation<scalar_type>::clone() const {

		return new continuous_set_simulation(*this);
}

/*
template<class scalar_type> continuous_set::ptr continuous_set_simulation<scalar_type>::get_ptr() {
	continuous_set::ptr p = boost::enable_shared_from_this<
			continuous_set_simulation<scalar_type> >::shared_from_this();
	return p;
}

template<class scalar_type> continuous_set::const_ptr continuous_set_simulation<scalar_type>::get_const_ptr() const {
	continuous_set::const_ptr p = boost::enable_shared_from_this<continuous_set_simulation<
			scalar_type> >::shared_from_this();
	return p;
}

*/

template<typename scalar_type>
continuous_set_simulation<scalar_type>* continuous_set_simulation<scalar_type>::create_universe() const{

	throw std::runtime_error("feature not implemented yet");
	return NULL;
}

template<typename scalar_type>
continuous_set_simulation<scalar_type>* continuous_set_simulation<scalar_type>::create_empty() const{

	return new continuous_set_simulation(myhbox);
}

template<typename scalar_type>
int continuous_set_simulation<scalar_type>::get_memory() const{

	throw std::runtime_error("feature not implemented yet");
	return 0;
}

template<typename scalar_type>
dimension_t continuous_set_simulation<scalar_type>::get_dim() const{

	return domain().size();
}

template<typename scalar_type>
math::tribool continuous_set_simulation<scalar_type>::is_empty() const{

	return my_map.size()==0;
}

template<typename scalar_type>
math::tribool continuous_set_simulation<scalar_type>::is_universe() const{
	return false;
}

template<typename scalar_type>
void continuous_set_simulation<scalar_type>::embed_variables(const variable_id_set& id_set){
	throw std::runtime_error("feature not implemented yet");
}

template<typename scalar_type>
void continuous_set_simulation<scalar_type>::existentially_quantify_variables(const variable_id_set& id_set){
	throw std::runtime_error("feature not implemented yet");
}

template<typename scalar_type>
continuous_set_predicate::ptr continuous_set_simulation<scalar_type>::get_predicate() const{
	throw std::runtime_error("feature not implemented yet");
}

template<typename scalar_type>
void continuous_set_simulation<scalar_type>::accept(const_visitor& d) const{
	d.dispatch(this);
}

template<typename scalar_type>
bool continuous_set_simulation<scalar_type>::contains_roots_of(const continuous_set_simulation& d) const{
	//TODO

	for (const_iterator it = d.begin();it!=d.end();++it) {
		if (!root_contains(it->first))
			return false;
	}
	return true;
}

template<typename scalar_type>
bool continuous_set_simulation<scalar_type>::root_contains(const state& v) const{
	// by definition v is contained in a root p if v-p \in myhbox

	for (const_iterator it = begin();it!=end();++it) {
		if (myhbox.contains(v-it->first))
			return true;
	}
	return false;
}

template<typename scalar_type>
continuous_set_ptr continuous_set_simulation<scalar_type>::intersection_with(const continuous_set_simulation<scalar_type> & d) const{
	//TODO
	throw std::runtime_error("feature not implemented yet");
}

template<typename scalar_type>
continuous_set_ptr continuous_set_simulation<scalar_type>::intersection_with(const math::lin_constraint_system<scalar_type> & cons) const{
	typedef typename trajectory::vector_type vector_t;

	// create a new set with the same hbox
	typename continuous_set_simulation<scalar_type>::ptr res(new continuous_set_simulation<scalar_type>(this->get_hbox()));

	for (const_iterator it = begin(); it != end(); ++it) {
		const trajectory& traj = it->second;

		trajectory new_traj(traj.domain());
		state root;
		bool has_points = false;

		// iterate through the states
		for (unsigned int k = 0; k < traj.size(); ++k) {
			state x = traj.get_state(k);
			scalar_type t = traj.get_time(k);

			if (math::maybe(cons.is_satisfied(x))) {
				new_traj.insert(t,x);
				if (!has_points) {
					root = x;
				}
				has_points = true;
			} else {
				if (has_points) {
					// add new_traj to result and start fresh
					res->insert(root,new_traj);
					new_traj = trajectory(traj.domain());
					has_points = false;
				}
			}
		}

		// add the last satisfying part to the output
		if (has_points) {
			// add new_traj to result
			res->insert(root,new_traj);
		}
	}

	return res;
}


template<typename scalar_type>
void continuous_set_simulation<scalar_type>::insert(
		const state& x, const trajectory& traj) {
	my_map.insert(x,traj);
}

template<typename scalar_type>
const typename continuous_set_simulation<scalar_type>::fhyperbox & continuous_set_simulation<scalar_type>::get_hbox() const{
	return myhbox;
}

template<typename scalar_type>
void continuous_set_simulation<scalar_type>::set_hbox(const fhyperbox & b){
	myhbox(b);
}

template<typename scalar_type>
unsigned int continuous_set_simulation<scalar_type>::size() const {
	return my_map.size();
}

template<typename scalar_type>
typename continuous_set_simulation<scalar_type>::iterator continuous_set_simulation<scalar_type>::begin() {
	return my_map.begin();
}

template<typename scalar_type>
typename continuous_set_simulation<scalar_type>::const_iterator continuous_set_simulation<scalar_type>::begin() const {
	return my_map.begin();
}

template<typename scalar_type>
typename continuous_set_simulation<scalar_type>::iterator continuous_set_simulation<scalar_type>::end() {
	return my_map.end();
}

template<typename scalar_type>
typename continuous_set_simulation<scalar_type>::const_iterator continuous_set_simulation<scalar_type>::end() const {
	return my_map.end();
}

template<typename scalar_type> bool continuous_set_simulation<scalar_type>::computes_support_vector() const {
	return true;
}
;

template<class scalar_type>
void continuous_set_simulation<scalar_type>::compute_support_impl(
		const math::vdom_vector<scalar_type>& l, scalar_type& max_value,
		math::vdom_vector<scalar_type>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	if (this->is_empty()) {
		is_empty = true;
		return;
	} else {
		is_empty = false;
		is_bounded = true;

		// we can remap since the coordinates in variable not in
		// the domain of *this are zero (by embedding).
		math::vdom_vector<scalar_type> lmap=l;
		lmap.remap(this->domain());

		max_value = scalar_type(0);
		// initialize support vector with zero, so irrelevant variables take
		// the value zero
		support_vec = math::vdom_vector<scalar_type>(this->domain());

		bool first=true;
		// run over all trajectories, and all points
		for (const_iterator it = begin(); it != end(); ++it) {
			// get the states in the trajectory as matrix
			const typename trajectory::matrix_type& M = it->second.get_states();
			typename trajectory::vector_type vals = M * lmap.get_vector();
			unsigned int max_index;
			scalar_type current_max_value = math::max(vals, max_index);
			if (first || current_max_value > max_value) {
				first = false;
				max_value = current_max_value;
				support_vec = math::vdom_vector<scalar_type>(this->domain(),
						M.vector_from_row(max_index));
			}
		}
	}
}

template<class scalar_type>
void continuous_set_simulation<scalar_type>::compute_support(const math::vdom_vector<
		Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	scalar_type max_val_st;
	math::vdom_vector<scalar_type> support_vec_st;
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
void continuous_set_simulation<scalar_type>::compute_support(const math::vdom_vector<
		double>& l, double& max_value, math::vdom_vector<double>& support_vec,
		bool& is_empty, bool& is_bounded) const {
	scalar_type max_val_st;
	math::vdom_vector<scalar_type> support_vec_st;
	math::vdom_vector<scalar_type> l_scalar_type = l.template convert_to<
			scalar_type> ();

	compute_support_impl(l_scalar_type, max_val_st, support_vec_st, is_empty,
			is_bounded);
	if (!is_empty && is_bounded) {
		max_value = convert_element<double> (max_val_st);
		support_vec = support_vec_st.template convert_to<double> ();
	}
}
}

#endif /* continuous_set_SIMULATION_HPP_ */
