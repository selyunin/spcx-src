/*
 * sfm_lin_cons_intersection.cpp
 *
 *  Created on: May 29, 2010
 *      Author: ray
 */

namespace continuous{

 namespace support_function {

template<typename scalar_type> sfm_lin_con_intersection<scalar_type>::sfm_lin_con_intersection(
		const typename sfm_cont_set<scalar_type>::const_ptr& s,
		const math::numeric::interval<unsigned int>& my_intv,
		const math::lin_constraint<scalar_type>& con,
		const std::string& minbrak_type,
		const double& inters_error ) :
	sf_unary<scalar_type> (s), my_interval(my_intv), my_con(con), my_minbrak_type(minbrak_type), my_intersection_error(inters_error){
}

template<typename scalar_type> sfm_lin_con_intersection<scalar_type>::sfm_lin_con_intersection(
		const typename sfm_cont_set<scalar_type>::const_ptr& s,
		const math::numeric::interval<unsigned int>& my_intv,
		const math::lin_constraint<scalar_type>& con,
		const affine_map& M, const std::string& minbrak_type,
		const double& inters_error) :
	sf_unary<scalar_type> (s, M), my_interval(my_intv), my_con(con), my_minbrak_type(minbrak_type), my_intersection_error(inters_error) {
}

template<typename scalar_type> sfm_lin_con_intersection<scalar_type>::~sfm_lin_con_intersection() {
}

template<typename scalar_type> sfm_lin_con_intersection<scalar_type>* sfm_lin_con_intersection<
		scalar_type>::clone() const {
	support_function_provider::const_ptr new_root(this->my_set->clone());
	/*
	if (this->get_map())
		return new sfm_lin_con_intersection<scalar_type> (new_root, my_con,
				*this->get_map());
	else
		return new sfm_lin_con_intersection<scalar_type> (new_root, my_con); */
	throw std::runtime_error(
			"sfm_lin_con_intersection : missing implementation clone");
}

template<typename scalar_type> int sfm_lin_con_intersection<scalar_type>::get_memory() const {
	throw std::runtime_error(
			"sfm_lin_con_intersection : missing implementation get_memory");
	return 0;
}

template<typename scalar_type> continuous_set_predicate::ptr sfm_lin_con_intersection<
		scalar_type>::get_predicate() const {
	throw std::runtime_error(
			"sfm_lin_con_intersection : missing implementation get_predicate");
	return continuous_set_predicate::ptr();
}

template<typename scalar_type> void sfm_lin_con_intersection<scalar_type>::print(
		std::ostream& os) const {
	throw std::runtime_error(
			"sfm_lin_con_intersection : missing implementation print");
}

template<typename scalar_type> const variable_id_set& sfm_lin_con_intersection<
		scalar_type>::get_variable_ids() const {
	return sf_unary<scalar_type>::get_variable_ids();
}
template<typename scalar_type> void sfm_lin_con_intersection<scalar_type>::reassign_primedness(
		unsigned int, unsigned int) {
	throw std::runtime_error(
			"sfm_lin_con_intersection : missing implementation reassign_primedness");
}
template<typename scalar_type> void sfm_lin_con_intersection<scalar_type>::increase_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sfm_lin_con_intersection : missing implementation increase_primedness");
}
template<typename scalar_type> void sfm_lin_con_intersection<scalar_type>::decrease_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sfm_lin_con_intersection : missing implementation decrease_primedness");
}

template<typename scalar_type> void sfm_lin_con_intersection<scalar_type>::compute_support(
		const math::vdom_vector<Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	throw std::runtime_error(
			"sfm_lin_con_intersection : missing implementation compute_support Rational");
}

template<typename scalar_type> void sfm_lin_con_intersection<scalar_type>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {

		//@todo incomplete implementation, to check for emptyness and boundedness
		is_empty = false;
		is_bounded = true;

		unsigned int samples; // to catch the # samples in the algo.
		sfm_cont_set<double>::const_ptr sfm_ptr = boost::static_pointer_cast<const sfm_cont_set<double>,const support_function_provider>(this->my_set);
		// take the constness away

		typename sfm_cont_set<scalar_type>::ptr nonconst_sfm =
				boost::const_pointer_cast<sfm_cont_set<scalar_type> >(sfm_ptr);
//		math::numeric::interval<double> min_intv = guard_intersection_support<double>(
	//			nonconst_sfm ,this->my_interval, l, this->get_constraint(), my_minbrak_type, my_intersection_error, samples );

		// This calls the optimized, simultaneous, branch and bound solution for the convex hull of the interval
		math::numeric::interval<double> min_intv = guard_intersection_supp_optz<double>(
						nonconst_sfm ,this->my_interval, l, this->get_constraint(), my_minbrak_type, my_intersection_error, samples );
		if(!min_intv.upper().is_finite())
			throw std::runtime_error("sfm_lin_con_intersection:compute_support: Unbounded sf value for the passed direction");

		if(is_bounded && !is_empty)
			max_value = min_intv.upper().get_val();
}

template<typename scalar_type>
std::vector<typename polyhedron<scalar_type>::ptr> sfm_lin_con_intersection<scalar_type>::get_outer_polys(vector_set dirs) const {
	std::vector<typename polyhedron<scalar_type>::ptr > res(my_interval.size()+1);
	//std::cout <<"Interval size:" << my_interval.size() + 1 << std::endl;

	// initialize empty polytopes
	for(unsigned int i = my_interval.lower().get_val(),j=0;i<=my_interval.upper().get_val();i++,j++){
		res[j] = typename polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
	}

	typename sfm_cont_set<scalar_type>::const_ptr sfm_ptr = boost::static_pointer_cast<const sfm_cont_set<scalar_type>,const support_function_provider>(this->my_set);

	// taking the constness away
	typename sfm_cont_set<scalar_type>::ptr nonconst_sfm =
			boost::const_pointer_cast<sfm_cont_set<scalar_type> >(sfm_ptr);

	// get the support function values for each of the sfm section member, guard intersection

	unsigned int samples;

	for (typename vector_set::const_iterator it = dirs.begin(); it
			!= dirs.end(); ++it) {
//		std::cout << "Calling guard intersection_simult..." << std::endl;

		math::vdom_vector<scalar_type> l = *it;
//		std::cout << "Sampling direction:" << l << std::endl;
		std::vector<scalar_type> sf_values = guard_intersection_simult<scalar_type>(
					nonconst_sfm ,this->my_interval, l, this->get_constraint(), my_minbrak_type, my_intersection_error, samples );
		assert(my_interval.size()+ 1 == sf_values.size());

//		std::cout << "sf_values print:" << std::endl;

/*
		for(unsigned int i = 0; i < sf_values.size(); i++)
			std::cout << sf_values[i] << ", ";
		std::cout << std::endl;
*/
		for(unsigned int i = my_interval.lower().get_val(),j=0;i<=my_interval.upper().get_val();i++,j++){
			math::lin_constraint<scalar_type> con(*it, -sf_values[j], LE);
//			std::cout << "added constraint:" << con << std::endl;
			res[j]->add_constraint(con);
		}
	}
	return res;
}

template<typename scalar_type>
std::vector<typename polyhedron<scalar_type>::ptr> sfm_lin_con_intersection<scalar_type>::get_chull_outer_polys(vector_set dirs, size_t split_size) const {
	assert(split_size >= 1);
	size_t intv_size = my_interval.size() + 1;
	size_t polys_size = intv_size/split_size;
	if(polys_size*split_size < intv_size)
		polys_size++;
	std::vector<typename polyhedron<scalar_type>::ptr > res(polys_size);
	std::cout <<"Polys size:" <<  polys_size << std::endl;

	// initialize empty polytopes
	for(unsigned int j=0;j < polys_size;j++){
		res[j] = typename polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type>());
	}

	typename sfm_cont_set<scalar_type>::const_ptr sfm_ptr = boost::static_pointer_cast<const sfm_cont_set<scalar_type>,const support_function_provider>(this->my_set);

	// taking the constness away
	typename sfm_cont_set<scalar_type>::ptr nonconst_sfm =
			boost::const_pointer_cast<sfm_cont_set<scalar_type> >(sfm_ptr);

	// get the support function values for each of the sfm section member, guard intersection

	unsigned int samples;
	unsigned int l, low = my_interval.lower().get_val();
	unsigned int u , up = my_interval.upper().get_val();

	math::numeric::interval<unsigned int> split_interval;

	for(unsigned int j = 0; j<polys_size; j++){
		l = low + j*split_size;
		split_interval.set_lower(l);
		u = l + split_size - 1;
		if(u > up)
			u = up;
		split_interval.set_upper(u);
		// splitting interval done
		typename support_function_provider::ptr sf_prov_set =
				typename support_function_provider::ptr(new sfm_lin_con_intersection(
						nonconst_sfm, split_interval, get_constraint(), my_minbrak_type, my_intersection_error));
		typename polyhedron<scalar_type>::ptr new_poly =
				typename polyhedron<scalar_type>::ptr(
					new constr_polyhedron<scalar_type> (continuous::compute_outer_poly(*sf_prov_set, dirs)));
		res[j] = new_poly;
	}
	return res;
}
template<typename scalar_type> const math::lin_constraint<scalar_type> & sfm_lin_con_intersection<
		scalar_type>::get_constraint() const {
	return my_con;
}

} // end of namespace support_function
} // end of namespace continuous
