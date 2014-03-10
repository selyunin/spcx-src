#include "sfm_cont_set.h"

#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"
#include "core/continuous/continuous_set_operators.h"
#include "core/continuous/continuous_set_operator_implementations/compute_transformation.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "core/post_operators/sfm_post/support_function.h"
#include "math/numeric/container_comp.h"
#include "math/vdom/vdom_vector.h"
#include "math/vdom/lin_constraint_system.h"
#include "math/lp_solving/fourier_motzkin/fm_elimination.h"

namespace continuous {
namespace support_function {

/** Declaration of static member variables. */
template<class scalar_type> typename sfm_cont_set<scalar_type>::output_format
		sfm_cont_set<scalar_type>::my_output_format =
				typename sfm_cont_set<scalar_type>::output_format();

///** Return a shared_ptr to *this.
// */
//template<class scalar_type> continuous_set::ptr sfm_cont_set<scalar_type>::get_ptr() {
//	continuous_set::ptr p = boost::enable_shared_from_this<continuous_set>::shared_from_this();
//	return p;
//}
//
///** Return a shared_ptr to const *this.
// */
//template<class scalar_type> continuous_set::const_ptr sfm_cont_set<scalar_type>::get_const_ptr() const {
//	continuous_set::const_ptr p = boost::enable_shared_from_this<continuous_set>::shared_from_this();
//	return p;
//}

template<class scalar_type> sfm_cont_set<scalar_type>::sfm_cont_set(
		const postc_params<scalar_type>& postc_prb, const math::matrix<scalar_type>& sfm,
		const direction_store& dirs,
		const index_to_variable_id_map_ptr& pnew_map,
		const math::lin_constraint_system<scalar_type>& nonstate_contraints):my_postc_prb(postc_prb),_my_sfm(sfm),my_directions(dirs),index_to_variable_id_map_provider(pnew_map),my_nonstate_contraints(nonstate_contraints) {
	if (sfm.size1() != dirs.size()) {
		//		std::cout << "#sfm rows" << sfm.size1() <<std::endl;
		//		std::cout << "#directions" << dirs.size() <<std::endl;

		throw std::runtime_error(
				" By SFM semantics, #directions in the sfm should be equal to #rows of sfm, which is violated");
	}
}

template<class scalar_type> sfm_cont_set<scalar_type>::sfm_cont_set(
		const postc_params<scalar_type>& postc_prb, const math::matrix<scalar_type>& sfm,
		const direction_store& dirs,
		const index_to_variable_id_map_ptr& pnew_map):my_postc_prb(postc_prb),_my_sfm(sfm),my_directions(dirs),index_to_variable_id_map_provider(pnew_map) {
	if (sfm.size1() != dirs.size()) {
		//		std::cout << "#sfm rows" << sfm.size1() <<std::endl;
		//		std::cout << "#directions" << dirs.size() <<std::endl;

		throw std::runtime_error(
				" By SFM semantics, #directions in the sfm should be equal to #rows of sfm, which is violated");
	}
}

template<class scalar_type> sfm_cont_set<scalar_type>::sfm_cont_set(
		const math::matrix<scalar_type>& sfm, const direction_store& dirs,
		const index_to_variable_id_map_ptr& pnew_map):my_postc_prb(postc_params<scalar_type>()),_my_sfm(sfm),my_directions(dirs),index_to_variable_id_map_provider(pnew_map) {

}

template<class scalar_type> sfm_cont_set<scalar_type>::~sfm_cont_set() {
}

template<class scalar_type> const math::matrix<scalar_type>& sfm_cont_set<
		scalar_type>::get_sfm() const {
	return _my_sfm;
}

template<class scalar_type> const typename sfm_cont_set<
scalar_type>::direction_store& sfm_cont_set<
		scalar_type>::get_directions() const {
	return my_directions;
}

template<class scalar_type> void sfm_cont_set<scalar_type>::set_sfm(
		math::matrix<scalar_type> new_sfm) {
	_my_sfm = new_sfm;
}

template<class scalar_type> void sfm_cont_set<scalar_type>::set_directions(
		const direction_store& new_dirs) {
	my_directions = new_dirs;
}

template<class scalar_type> const postc_params<scalar_type>& sfm_cont_set<
		scalar_type>::get_postc_prb() const {
	return my_postc_prb;
}

template<class scalar_type> sfm_cont_set<scalar_type>* sfm_cont_set<scalar_type>::create_universe() const {
	math::matrix<scalar_type> new_sfm;
	direction_store new_dirs;

	postc_params<scalar_type> empty_postc_prb;
	continuous_set::const_ptr empty_set_ptr = continuous_set::const_ptr(new constr_polyhedron<scalar_type> ());

	empty_postc_prb.initial_set_ptr = empty_set_ptr;


	//root_const_ptr univ_root(new constr_polyhedron<scalar_type> ());
	return new sfm_cont_set(empty_postc_prb, new_sfm, new_dirs,
			index_to_variable_id_map::empty_map()); // An empty index_to_id_map is created
}

template<class scalar_type> void sfm_cont_set<scalar_type>::remove_empty_columns() {
	unsigned int N = _my_sfm.size2();
	if (N > 0) {
		// Find the empty columns
		std::set<unsigned int> empty_columns;
		for (unsigned int i = 0; i < N; ++i) {
			if (is_column_empty(i)) {
				empty_columns.insert(i);
			}
		}
		if (empty_columns.size() < N) {
			// Redefine the matrix
			unsigned int r = my_directions.size();
			math::matrix<scalar_type> m_i(r, N - empty_columns.size());

			unsigned int k = 0;
			for (unsigned int j = 0; j < N; j++) {
				if (empty_columns.find(j) == empty_columns.end()) {
					for (unsigned int i = 0; i < r; i++) {
						m_i(i, k) = _my_sfm(i, j);
					}
					++k;
				}
			}
			_my_sfm = m_i;
		} else {
			typename sfm_cont_set<scalar_type>::ptr empty_set =
					typename sfm_cont_set<scalar_type>::ptr(create_empty());
			*this = *empty_set;
		}
	}
}

template<class scalar_type> void sfm_cont_set<scalar_type>::intersection_with_constraint(
		const math::lin_constraint<scalar_type>& con, bool add_to_ini) {
	if (add_to_ini) {
		typename constr_polyhedron<scalar_type>::ptr con_poly(new constr_polyhedron<scalar_type>());
		con_poly->add_constraint(con);
		continuous_set::ptr con_set(con_poly);
		// @todo only add constraints to the sets that are concerned by it (how to decide?)
		//my_postc_prb.initial_set_ptr = compute_intersection(my_postc_prb.initial_set_ptr,con_set);
		my_postc_prb.invariant_set_ptr = compute_intersection(my_postc_prb.invariant_set_ptr,con_set);
		//my_postc_prb.input_set_ptr = compute_intersection(my_postc_prb.input_set_ptr,con_set);
	}

	// get the constraint in canonic form a.x <= b
	math::vdom_vector<scalar_type> a = con.get_normal();
	scalar_type b = -con.get_canonic_inh_coeff();
	comparison_operator my_sign = con.get_canonic_sign();

	// Check if the domains are compatible
	if (this->domain()!=a.domain()) {
		if (!this->domain().contains_variables(a.domain())) {
			// then this is a constraint that the sfm can't handle
			// so add them to the nonstate constraints
			my_nonstate_contraints.insert(con);
			return;

//			// there are variables in a that are not in *this,
//			// so we need to add them to *this
//			variable_id_set new_vars = a.domain().get_variable_ids();
//			set_difference_assign(new_vars,this->domain().get_variable_ids());
//			embed_variables(new_vars);
		}
		// now all variables are there, it suffices to reorder
		a.reorder(this->domain());
	}
	// intersect with a.x <= b.
	if (false) // GF 2011-12-02 deactivated because it interferes with precise intersection
	{
	std::pair<unsigned int, scalar_type> pr = this->extend_sfm(a.get_vector());
	unsigned int row_index = pr.first;
	scalar_type row_f = pr.second;

	unsigned int N = get_size();
	for (unsigned int j = 0; j < N; j++) {
		if (_my_sfm(row_index, j) > b * row_f) {
			_my_sfm(row_index, j) = b * row_f;
		}
	}
	}

	if (my_sign == EQ) {
		// we already intersect with
		//    a.x <= b.
		// so now intersect with
		//    a.x >= b.
		math::lin_expression<scalar_type> lopp(-a, b);
		math::lin_constraint<scalar_type> con_opp(lopp, LE);
		intersection_with_constraint(con_opp);
	}
}

template<class scalar_type> void sfm_cont_set<scalar_type>::intersection_with_poly(
		const polyhedron<scalar_type>& poly) {
	typename constr_polyhedron<scalar_type>::ptr con_poly(new constr_polyhedron<scalar_type>());
	con_poly->add_constraints(poly.get_constraints());
	continuous_set::ptr con_set(con_poly);
	// @todo only add constraints to the sets that are concerned by it (how to decide?)
	//my_postc_prb.initial_set_ptr = compute_intersection(my_postc_prb.initial_set_ptr,con_set);
	if (my_postc_prb.invariant_set_ptr)
		my_postc_prb.invariant_set_ptr = compute_intersection(my_postc_prb.invariant_set_ptr,con_set);
	else
		my_postc_prb.invariant_set_ptr = con_set;
	// my_postc_prb.input_set_ptr = compute_intersection(my_postc_prb.input_set_ptr,con_set);

	const math::lin_constraint_system<scalar_type>& cons = *poly.get_constraints();
	for (typename math::lin_constraint_system<scalar_type>::const_iterator it=cons.begin();it!=cons.end();++it) {
		intersection_with_constraint(*it,false);
	}
}

template<class scalar_type> typename sfm_cont_set<scalar_type>::ptr sfm_cont_set<scalar_type>::intersection_with_poly_improved(
		const polyhedron<scalar_type>& poly) const {
	typename sfm_cont_set<scalar_type>::ptr s =
			typename sfm_cont_set<scalar_type>::ptr(clone());
	s->intersection_with_poly(poly);
	return s;
//
//	//std::cout << "sfm:Inside improved intersection:" << std::endl;
//	typename math::lin_constraint_system<scalar_type>::const_ptr lin_cons =
//			poly.get_constraints();
//	unsigned int k = lin_cons->size();
//	unsigned int r = directions.size();
//	unsigned int N = _my_sfm.size2();
//	std::list<math::vector<scalar_type> > new_dirs; // empty list created??
//
//	// count equalities
//	unsigned int eq_count = 0;
//
//	for (typename std::list<math::vector<scalar_type> >::const_iterator
//			it = directions.begin(); it != directions.end(); it++) {
//		new_dirs.push_back(*it);
//	}
//	std::vector<int> redundant(k,-1); // vector to say which constraint is redundant and which one is not.
//
//	unsigned int pos = 0;
//
//	for (typename math::lin_constraint_system<scalar_type>::const_iterator it =
//			lin_cons->begin(); it != lin_cons->end(); ++it,pos++) {
//
//		int h = 0,p = 0;
//
//		if(it->is_equality()){ // count the number of equality constraints.
//			eq_count++;
//		}
//
//		for (typename std::list<math::vector<scalar_type> >::const_iterator
//				iter = directions.begin(); iter != directions.end(); iter++, h++) {
//				if(*iter == it->get_canonic_l().get_vector()){
//					redundant[pos] = h;
//					k--;
//					break;
//				}
//		}
//	}
//	math::vector<int> eq_cons(eq_count,-1);
//
//
//	// redundant vector is initialised at this point
//	//if (k == 0){ // All constraints are redundant
//	//	return ptr(new sfm_cont_set<scalar_type>(get_sfm(),get_directions(), get_index_to_variable_id_map()));
//	//}
//	math::matrix<scalar_type> m_i(r + k + eq_count, N); // k now denotes exactly the number of new non-redundant constraints +
//														// additional terms to include the neg constraint(s) for equality(ies).
//	m_i.submatrix_assign(_my_sfm,0,r,0,N);
////		for (unsigned int i = 0; i < r; i++)
////			for (unsigned int j = 0; j < N; j++)
////				m_i(i, j) = _my_sfm(i, j);
//
//	unsigned int i = r,g=0,c=0;
//
//	for(typename math::lin_constraint_system<scalar_type>::const_iterator it =
//			lin_cons->begin(); it != lin_cons->end(); ++it, g++) {
//		// Check if the constraint is redundant
//		// When the constraint is redundant, get the minimum of the corresponding sfm and the constraint inh_coeff
//		if(it->is_equality()){
//			eq_cons[c] = g; // meaning that gth constraint of the poly is an equality.
//			c++;
//		}
//		if(redundant[g]!= -1){
//			for(unsigned int j = 0; j < N; j++) {
//				scalar_type b = scalar_type(-1) * it->get_canonic_inh_coeff();
//				if(m_i(redundant[g], j) > b)
//					m_i(redundant[g],j) = b;
//			}
//		}
//		else { // non - redundant constraint
//			if (get_index_to_variable_id_map()
//					== it->get_l().get_index_to_variable_id_map()) {
//				new_dirs.push_back(it->get_canonic_l().get_vector()); // directions from the con_poly added.
//			} else {
//				math::vdom_vector<scalar_type> rvec = it->get_canonic_l().get_vdom_vec();
//				// Note that reordering will throw is a nonzero coeff of rvec
//				// has no corresponding id.
//				rvec.reorder(get_index_to_variable_id_map());
//				new_dirs.push_back(rvec.get_vector());
//			}
//			for (unsigned int j = 0; j < N; j++) {
//				m_i(i, j) = scalar_type(-1) * it->get_canonic_inh_coeff(); // = 0? because l above contains the inhomogeneous term?
//			}
//			i++; // sfm row size increased
//		}
//	}
//	// Add equality constraints in the list of directions.
//
//	g = 0;
//	i = 0;
//	if(eq_count != 0){
//		for(typename math::lin_constraint_system<scalar_type>::const_iterator it =
//			lin_cons->begin(); it != lin_cons->end(); ++it, g++) {
//			if(eq_cons[i] == g){ // meaning that the current poly constraint is equality constraint
//				if (get_index_to_variable_id_map()
//									== it->get_l().get_index_to_variable_id_map()) {
//					// push back the complement of the constraint
//					new_dirs.push_back(-it->get_canonic_l().get_vector()); // directions from the con_poly added.
//				}
//				else{
//					math::vdom_vector<scalar_type> rvec = it->get_normal();
//					// Note that reordering will throw if a nonzero coeff of rvec
//					// has no corresponding id.
//					rvec.reorder(get_index_to_variable_id_map());
//					new_dirs.push_back(-rvec.get_vector());
//				}
//				for (unsigned int j = 0; j < N; j++) {
//					scalar_type b = it->get_canonic_inh_coeff(); // without the -ve sign.
//						m_i(redundant[g],j) = b; // add the coefficient in all the col entries.
//				}
//				i++;
//			}
//		}
//	}
//
//	postc_params<scalar_type> new_postc_prb = get_postc_prb(); // keep the postc prb as it is
//
//	sfm_cont_set<scalar_type>* new_sfm = new sfm_cont_set<scalar_type> (
//			new_postc_prb, m_i, new_dirs, get_index_to_variable_id_map());
//	new_sfm->remove_empty_columns();
//	return ptr(new_sfm);
}

template<class scalar_type> bool sfm_cont_set<scalar_type>::contains_initially(
		const continuous_set::const_ptr& s) const {
	bool result;
	if (my_postc_prb.initial_set_ptr && s) {
		result = my_postc_prb.initial_set_ptr->contains(s);
	} else
		result = false;
	//	std::cout << "Does " << *this << " contain " << s << "?"  << std::flush;
	//	if (result) std::cout << "yes"; else std::cout << "no";
	//	std::cout << std::endl;
	return result;
}

template<class scalar_type> bool sfm_cont_set<scalar_type>::contains_initial_set(
		const sfm_cont_set& s) const {
	return contains_initially(s.get_postc_prb().initial_set_ptr);
}

template<class scalar_type> polyhedron_collection<scalar_type> sfm_cont_set<scalar_type>::get_outer_polytope_collection(unsigned int i, unsigned int j) const {
	size_t N = get_size();

	typename polyhedron<scalar_type>::ptr element;
	if (N == 0) {
		// return a universe set
		element = typename polyhedron<scalar_type>::ptr(
				new constr_polyhedron<scalar_type> ());
	} else {
		element = typename polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type> (
				get_polytope(i)));
	}
	polyhedron_collection<scalar_type> coll(element);
	++i;
	for (; i <= j; ++i) {
		element = typename polyhedron<scalar_type>::ptr(new constr_polyhedron<scalar_type> (
						get_polytope(i)));
		coll.insert(element, false); // insert w/o redundancy check
	}
	return coll;
}

template<class scalar_type> polyhedron_collection<scalar_type> sfm_cont_set<scalar_type>::get_outer_polytope_collection() const {
	size_t N = get_size();
	if (N>=1) {
		return get_outer_polytope_collection(0,N-1);
	} else {
		return get_outer_polytope_collection(0,0);
	}
}

template<class scalar_type> constr_polyhedron<scalar_type> sfm_cont_set<
		scalar_type>::get_polytope(unsigned int j) const {
	unsigned int k = my_directions.size();
	unsigned int N = _my_sfm.size2();
	if (j > N - 1)
		throw std::runtime_error(
				"Invalid Parameter j - should be less than N\n");

	constr_polyhedron<scalar_type> poly_j;
	unsigned int i = 0;
	for (typename direction_store::const_iterator
			it = my_directions.begin(); it != my_directions.end(), i < k; it++, i++) {
		math::lin_expression<scalar_type> lexp(it->first, -_my_sfm(it->second, j),
				get_index_to_variable_id_map());
		math::lin_constraint<scalar_type> l_cons(lexp, LE);
		poly_j.add_constraint(l_cons);
	}
	//poly_j.add_constraints(my_nonstate_contraints);
	if (my_postc_prb.invariant_set_ptr) {
		typename polyhedron<scalar_type>::const_ptr inv = boost::static_pointer_cast<const polyhedron<scalar_type> >(my_postc_prb.invariant_set_ptr);
		poly_j.add_constraints(inv->get_constraints());
	}
	return poly_j;
}

template<class scalar_type>
bool sfm_cont_set<scalar_type>::is_column_empty(unsigned int j) const {

	return get_polytope(j).is_empty();
}

template<class scalar_type>
math::tribool sfm_cont_set<scalar_type>::is_empty_outer() const {
	size_t N = get_size();

	if (N == 0) {
		// universe set
		return false;
	}

	// loop until nonempty set found or all steps are processed
	unsigned int i = 0;
	math::tribool res_empty = true;
	while (i < N && math::maybe(res_empty)) {
		res_empty = res_empty && get_polytope(i).is_empty();
		++i;
	}
	return res_empty;
}

//template<class scalar_type> unsigned int sfm_cont_set<scalar_type>::get_first_j(direction d, scalar_type b) const {
//
//	// check if the direction d is already in the list of directions
//	unsigned int i = 0, row_index;
//	bool d_exist = false;
//
//	for (direction_store::const_iterator
//			it = my_directions.begin(); it != my_directions.end(); it++, i++) {
//		if(*it == d){
//			row_index = i;
//			d_exist = true;
//			break;
//		}
//	}
//	if(!d_exist) {
//		sfm_cont_set<scalar_type>* nonconst_this = const_cast<sfm_cont_set<scalar_type>*>(this);
//		row_index = nonconst_this->extend_sfm(d).first;
//	}
//	/* Iterate over the col entries in the row_index of the sfm and return the
//	 * first col in which the entry is less than or equals parameter b
//	 */
//	for(unsigned int j=0;j< _my_sfm.size2();j++){
//		if(_my_sfm(row_index, j) >= b)
//			return j;
//	}
//	return  get_size();
//}
template<class scalar_type> std::pair<unsigned int,scalar_type> sfm_cont_set<scalar_type>::extend_sfm(const direction& d) {

	// check if the direction d is already in the list of directions
	unsigned int i = 0, row_index;
	bool d_exist = false;
	scalar_type d_norm = d.infinity_norm();
	scalar_type f = scalar_type(1)/d_norm;

	// sanity check
	if(_my_sfm.size1() == 0 && _my_sfm.size2() == 0){
		throw std::runtime_error("sfm_cont_set:extend_sfm: Cannot extend uninitialized sfm\n");
	}
	// check if d is already in the store
	direction d_normed =  d/d_norm;
	typename direction_store::const_iterator it = my_directions.find(d_normed);
	if (it != my_directions.end()) {
		// direction found, return row index and scaling factor
		return std::make_pair(it->second,f);
	}
	else{
		//std::cout << "now at directions:" << get_directions().size() << std::endl;
		//std::cout << "Extended in the following direction:" << d << std::endl;

		/* compute the new sfm with the added direction.
		 */
		std::list<direction> new_directions;
		new_directions.push_back(d_normed); // List with only the new direction.

		unsigned int N;
		if(get_size()==0)
			N = (unsigned int) ceil(my_postc_prb.time_horizon/my_postc_prb.delta);
		else
			N = get_size();

		scalar_with_infinity<scalar_type> my_scalar;

		typename continuous::support_function_provider::const_ptr sf_provider_init_ptr,sf_provider_input_ptr ;

		typename math::matrix<scalar_type> new_sfm;
		if (sf_provider_init_ptr = boost::dynamic_pointer_cast<
				const continuous::support_function_provider>(
				my_postc_prb.initial_set_ptr)) {
			if (!my_postc_prb.input_set_ptr || (sf_provider_input_ptr
					= boost::dynamic_pointer_cast<
							const continuous::support_function_provider>(
							my_postc_prb.input_set_ptr))) {
				//std::cout << my_postc_prb.dynamics << std::endl;
				//std::cout << new_directions << std::endl;
				//std::cout << get_index_to_variable_id_map() << std::endl;

				// if this object doesn't yet have an evaluator, create one
				// @todo: check if the state of *this changed
				if (!my_evaluator) {
					my_evaluator = boost::shared_ptr<sf_evaluator<scalar_type> > (new sf_evaluator<scalar_type>(
							my_postc_prb.dynamics_A, my_postc_prb.dynamics_b,
							get_index_to_variable_id_map(), sf_provider_init_ptr,
							sf_provider_input_ptr));
				}

				new_sfm = my_evaluator->compute_matrix(new_directions,
						scalar_type(my_postc_prb.delta), N, my_postc_prb.delta_vec);
			} else {
				throw std::runtime_error(
						"sfm_cont_set::extend_sfm::input_set type not supported for support function computation\n");
			}
		} else {
			throw std::runtime_error(
					"sfm_cont_set::extend_sfm::initial_set type not supported for support function computation\n");
		}
		// add this new_sfms's only row to the end of _my_sfm
//		std::cout << "my_sfm.size1=" << _my_sfm.size1() << ", my_sfm.size2=" << _my_sfm.size2() << std::endl;
//		std::cout << "new_sfm.size1=" << new_sfm.size1() << ", new_sfm.size2=" << new_sfm.size2() << std::endl;
		_my_sfm.resize(_my_sfm.size1()+1,_my_sfm.size2());
		unsigned int last_row = _my_sfm.size1()-1;
		_my_sfm.submatrix_assign(new_sfm,last_row,last_row+1,0,_my_sfm.size2());

		// add to the directions of *this
		my_directions.insert_missing(d_normed,last_row);
		return std::make_pair(last_row,f);
	}
}
template<class scalar_type> constr_polyhedron<scalar_type> sfm_cont_set<
		scalar_type>::compute_template_hull(index_interval intv) const {

	unsigned int i,N;
	if (intv.lower().is_finite())
		i = intv.lower().get_val();
	else
		i = 0;
	if (intv.upper().is_finite())
		N = intv.upper().get_val()+1;
	else
		N = _my_sfm.size2();

	constr_polyhedron<scalar_type> poly_j;
	if (N > i) {
		for (typename direction_store::const_iterator
				it = my_directions.begin(); it != my_directions.end(); ++it) {
			index_type row_index = it->second;
			scalar_type max_entry;
			bool max_uninit = true;
			for (unsigned int j = i; j < N; ++j) {
				bool not_empty = !is_column_empty(j);
				if (not_empty) {
					if (max_uninit || max_entry < _my_sfm(row_index, j)) {
						max_entry = _my_sfm(row_index, j);
						max_uninit = false;
					}
				}
			}
			math::lin_expression<scalar_type> lexp(it->first, -max_entry,
					get_index_to_variable_id_map());
			math::lin_constraint<scalar_type> l_cons(lexp, LE);
			poly_j.add_constraint(l_cons);
		}
	}
	return poly_j;
}

//template<class scalar_type> sfm_cont_set<scalar_type> sfm_cont_set<scalar_type>::affine_transform(
//		const math::affine_map<scalar_type>& M) const {
//
//	math::matrix<scalar_type> A=M.get_A().get_matrix();
//	math::vector<scalar_type> v=M.get_b().get_vector();
//
//	// if M is not in the same domain as *this try to reorder
//	math::affine_map<scalar_type> Mreordered;
//	if (this->domain()!=M.domain() || this->domain()!=M.codomain()) {
//		Mreordered=M;
//		try {
//			Mreordered.reorder(this->domain(),this->domain());
//		} catch (std::exception& e) {
//			std::stringstream sM,sdom;
//			logger::copyfmt_to(sM);
//			logger::copyfmt_to(sdom);
//			sM << M;
//			sdom << this->domain();
//			throw basic_exception("This type of affine map is not supported at the moment, it should be over the domain "+sdom.str()+":\n"+sM.str(),e);
//		}
//		A=Mreordered.get_A().get_matrix();
//		v=Mreordered.get_b().get_vector();
//	}
//
//	unsigned int k = _my_sfm.size1();
//	unsigned int N = _my_sfm.size2();
//
//	if (A.size1() != A.size2())
//		throw std::runtime_error(
//				"Cannot perform linear transform for non-square transformation matrix for the moment");
//	if (A.size1() != (directions.begin())->size()) {
//		std::cerr << "transformation matrix : " << A;
//		std::cerr << "transformation vector : " << v;
//		std::cerr << "first direciton vector : " << *directions.begin();
//		throw std::runtime_error(
//				"sfm_cont_set<scalar_type>::affine_transform: dimension of the transformation matrix should be equal to the vector space dimension");
//	}
//
//	math::matrix<scalar_type> trans_sfm(k, N);
//	bool singular;
//	// Compute the inverse of the linear transformation matrix
//	math::matrix<scalar_type> A_inverse(A.inverse(singular));
//
//	if (singular) { // If the transformation matrix is singular, compute a tight polyhedral over-approximation of the exact map.
//
//		bool is_empty, is_bounded;
//		scalar_type max;
//		math::vdom_vector<scalar_type> sup_vec;
//		//		math::matrix<scalar_type> trans_sfm(k, N);
//
//		std::list<math::vector<scalar_type> > trans_dirs;
//
//		unsigned int i = 0;
//		for (typename std::list<math::vector<scalar_type> >::const_iterator
//				it = directions.begin(); it != directions.end(); it++, i++) {
//			math::vector<scalar_type> t =
//					(A.transpose()).multiply_vector(*it);
//			math::vector<scalar_type> p = *it;
//			//		std::cout << "vector t:" << t << std::endl;
//			//		std::cout << "vector p:" << p << std::endl;
//
//			for (unsigned int j = 0; j < N; j++) {
//				math::vdom_vector<scalar_type> lexp(t, get_index_to_variable_id_map());
//				//			std::cout << "linear expression:" << lexp << std::endl;
//				(get_polytope(j)).compute_support(lexp, max, sup_vec, is_empty,
//						is_bounded);
//				//			 std::cout << "max value:" << max << std::endl;
//				if (!is_bounded) {
//					std::cerr << get_polytope(j);
//					throw std::runtime_error(
//							"Attempt to compute support function of unbounded set");
//				}
//				if (is_empty) {
//					throw std::runtime_error(
//							"Attempt to compute support function of empty sfm column");
//				}
//				//			std::cout << "scalar_prod" << v.scalar_product(p) << std::endl;
//				trans_sfm(i, j) = max + v.scalar_product(p);
//			}
//		}
//		// @todo compute the root's outer poly approx.
//
//		postc_params<scalar_type> new_postc_prb = my_postc_prb; // The prb remains the same, only the initial set changes.
//
//		if (my_postc_prb.initial_set_ptr) {
///*
//			reset_affine_transform<scalar_type> new_transf(A, v,
//					get_index_to_variable_id_map());
//			new_postc_prb.initial_set_ptr = my_postc_prb.initial_set_ptr->compute_transformation(new_transf);
//*/
//// Note: We don't want to compute projection on the root (expensive!)
//			new_postc_prb.initial_set_ptr = continuous_set::const_ptr();
//		}
//
//		sfm_cont_set<scalar_type> new_sfm(new_postc_prb, trans_sfm, directions,
//				get_index_to_variable_id_map());
//		return new_sfm;
//
//	}//end if
//
//	/* If the linear transformation matrix is non-singular, then we compute the exact map using the following result:
//	 * Given a polytope Fx <= b, the polytope after affine transformation can be given as F.A^{-1}x <= b + FA^{-1}v,
//	 * where the affine tranformation is given as x --> Ax+v.
//	 */
//
//	else {
//		vector_list trans_directions;
//		vector_type d;
//		scalar_type c_prime;
//		unsigned int i = 0;
//		for (typename std::list<math::vector<scalar_type> >::const_iterator
//				it = directions.begin(); it != directions.end(); it++, i++) {
//			//compute the new directions iteratively.
//			d = (*it) * A_inverse; // computing AF^{-1} row-wise.
//			trans_directions.push_back(d);
//			// the rhs of the inequality, b + FA^{-1}v is computed componentwise. c_prime contains each component.
//
//			c_prime = scalar_product(d, v);
//
//			// this component is added to the sfm's corresponding component row entries to have the effect of b + FA^{-1}v.
//
//			for (unsigned int j = 0; j < N; j++) {
//				trans_sfm(i, j) = _my_sfm(i, j) + c_prime;
//			}
//		}
//		//root_const_ptr new_root = root_const_ptr();
//		postc_params<scalar_type> new_postc_prb;
//
//		if (my_postc_prb.initial_set_ptr) {
//			math::affine_map<scalar_type> M(domain(),A,v);
//			reset_affine_transform<scalar_type> new_transf(M);
//			new_postc_prb.initial_set_ptr = compute_transformation(*my_postc_prb.initial_set_ptr,new_transf);
//		}
//		sfm_cont_set new_sfm(new_postc_prb, trans_sfm, trans_directions,
//				get_index_to_variable_id_map());
//		//std::cout << "transformed sfm: " << new_sfm << std::endl;
//		return new_sfm;
//	}
//}

template<class scalar_type> sfm_cont_set<scalar_type>* sfm_cont_set<scalar_type>::clone() const {
	return new sfm_cont_set<scalar_type> (*this);
}

template<class scalar_type> sfm_cont_set<scalar_type>* sfm_cont_set<scalar_type>::create_empty() const {
	return empty_set();
}

template<class scalar_type> sfm_cont_set<scalar_type>* sfm_cont_set<scalar_type>::empty_set() {
	math::matrix<scalar_type> new_sfm;
	direction_store new_dirs;

	continuous_set::ptr empty_root = continuous_set::ptr(new constr_polyhedron<scalar_type> ());
	empty_root = continuous_set::ptr(empty_root->create_empty());

	postc_params<scalar_type> new_postc_prb;
	new_postc_prb.initial_set_ptr = empty_root;

	return new sfm_cont_set(new_postc_prb, new_sfm, new_dirs,
			index_to_variable_id_map::empty_map()); // An empty index_to_id_map is created
}

template<class scalar_type> int sfm_cont_set<scalar_type>::get_memory() const {
	throw std::runtime_error(
			"sfm_cont_set<scalar_type>::get_memory() missing implementation");
	return 0;
}

template<class scalar_type> unsigned int sfm_cont_set<scalar_type>::get_dim() const {
	//	throw std::runtime_error("missing implementation");
	return get_index_to_variable_id_map()->dimensions();
}

template<class scalar_type>
size_t sfm_cont_set<scalar_type>::get_size() const {
	return _my_sfm.size2();
}

template<class scalar_type> math::tribool sfm_cont_set<scalar_type>::is_empty() const {
	// The set is empty iff the initial states are empty.
	// An null initial_set_ptr signifies an empty set.
	if (!my_postc_prb.initial_set_ptr || math::definitely(my_postc_prb.initial_set_ptr->is_empty()) || get_size()==0) {
		return true;
	} else {
		// check if first outer poly is empty

		// note that here size>0 because of the previous branch
		// return get_polytope(0).is_empty();
		return math::indeterminate();
	}
	return math::indeterminate();
}

template<class scalar_type> math::tribool sfm_cont_set<scalar_type>::is_universe() const {
	// The set is considered universal if the initial set is universal
	// (we're talking bounded time, so that's the only solution)
	if (my_postc_prb.initial_set_ptr && my_postc_prb.initial_set_ptr->is_universe()) {
		return true;
	} else
		return math::indeterminate();
}

template<class scalar_type> void sfm_cont_set<scalar_type>::embed_variables(
		const variable_id_set& id_set) {
	// add the ids to _index_to_variable_id_map_ptr
	std::size_t old_dim = get_index_to_variable_id_map()->dimensions();
	index_to_variable_id_map_ptr p =
			get_index_to_variable_id_map()->get_map_with_ids_added(id_set);
	set_index_to_variable_id_map(p);

	std::size_t new_dim = get_index_to_variable_id_map()->dimensions();

	// adapt the directions with embedding
	direction_store new_directions;
	for (typename direction_store::const_iterator
			it = my_directions.begin(); it != my_directions.end(); it++) {
		math::vector<scalar_type> dir_new(new_dim, scalar_type(0));
		dir_new.subvector_assign(it->first, 0, old_dim);
		new_directions.insert_missing(dir_new,it->second);
	}
	set_directions(new_directions);
	//	throw std::runtime_error("missing implementation");

	// embed the problem def
	if (my_postc_prb.dynamics_A.size1() > 0) {
		math::matrix<scalar_type> A_new(new_dim, new_dim, scalar_type(0));
		A_new.submatrix_assign(my_postc_prb.dynamics_A, 0, old_dim, 0, old_dim);
		my_postc_prb.dynamics_A = A_new;
	}
	if (my_postc_prb.dynamics_b.size() > 0) {
		math::vector<scalar_type> b_new(new_dim, scalar_type(0));
		b_new.subvector_assign(my_postc_prb.dynamics_b, 0, old_dim);
		my_postc_prb.dynamics_b = b_new;
	}
}

template<class scalar_type> void sfm_cont_set<scalar_type>::existentially_quantify_variables(
		const variable_id_set& id_set) {
	if (!id_set.empty()) {
		throw std::runtime_error(
				"sfm_cont_set<scalar_type>::existentially_quantify_variables: missing implementation");
	}
}

/** Accept a visitor. */
template<class scalar_type> void sfm_cont_set<scalar_type>::accept(
		dispatching::dispatcher<continuous_set_typelist>& d) const {
	d.dispatch(this);
}

/** Output as a stream of characters.
 *
 * @todo Not sure of the useful output format.
 */
template<class scalar_type> void sfm_cont_set<scalar_type>::print(
		std::ostream& os) const {

	typename sfm_cont_set<scalar_type>::output_format of = get_output_format();

	if (of == sfm_cont_set::JVX) {
		os
				<< "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>"
				<< std::endl;
		os
				<< "<!DOCTYPE jvx-model SYSTEM \"http://www.javaview.de/rsrc/jvx.dtd\"> "
				<< std::endl;
		os << "<jvx-model>" << std::endl;
		os << "<geometries>" << std::endl;
		// construct a the series of translated sets
		for (unsigned int i = 0; i < _my_sfm.size2(); i++) {
			os << "   <geometry name=\"states" << i << "\">" << std::endl;
			constr_polyhedron<scalar_type> poly = get_polytope(i);
			poly.set_output_format(polyhedron<scalar_type>::JVX);
			os << poly;
			os << "   </geometry>" << std::endl;
		}
		os << "</geometries>" << std::endl;
		os << "</jvx-model>" << std::endl;
	} else if (of == sfm_cont_set::DOUBLE_GENERATORS) {
		if (get_dim() > 2)
			throw std::runtime_error(
					"Cannot Display a SFM represented Polytope set of dimension greater than 2");
		for (unsigned int i = 0; i < _my_sfm.size2(); i++) {
			constr_polyhedron<scalar_type> poly = get_polytope(i);
			poly.set_output_format(polyhedron<scalar_type>::DOUBLE_GENERATORS);
			os << poly;
			os << std::endl << std::endl;
		}
	} else if (of == sfm_cont_set::TEXTUAL) {
		//		os << of.preamble;
		os << "{";
		os << _my_sfm;

		//		os << of.epilogue;
		os << "}";

		//		os << of.preamble;

		//	os << directions;
		os << "{";
		/*Need to print the position vdom vector */
		positional_vdomain my_dom = positional_vdomain(get_index_to_variable_id_map());
		os << "variable to dimension map:" << my_dom << std::endl;

		for (typename direction_store::const_iterator
				it = my_directions.begin(); it != my_directions.end();) {
			os << it->first;
			if (++it != my_directions.end())
				os << ",";
		}
		//		os << of.epilogue;
		os << "}";
		if (my_postc_prb.initial_set_ptr) {
			os << "<initial_set:" << my_postc_prb.initial_set_ptr << ">";
		} else {
			os << " <no root>";
		}
	}
}

template<typename scalar_type> bool sfm_cont_set<scalar_type>::computes_support_vector() const {
	return false;
}
;

template<class scalar_type>
void sfm_cont_set<scalar_type>::compute_support(
		const math::vdom_vector<Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	throw std::runtime_error(
			"sfm_cont_set<scalar_type>::compute_support: missing implementation");
}

template<class scalar_type>
void sfm_cont_set<scalar_type>::compute_support(const math::vdom_vector<
		double>& l, double& max_value, math::vdom_vector<double>& support_vec,
		bool& is_empty, bool& is_bounded) const {
	// @todo what if there are no columns?
	support_vec=math::vdom_vector<double>();

	if (this->is_empty()) {
		is_empty = true;
		is_bounded = true;
	} else {
		// @todo do proper emptyness check
		is_empty = false;
		is_bounded = true;

		// if l has variables not in the domain of *this, the problem is unbounded
		math::vdom_vector<double> v = l;
		if (v.domain() != this->domain()) {
			if (!this->domain().contains_variables(v.domain())) {
				// check whether there are nonzero variables not in the domain of *this
				math::vdom_vector<double> vtemp = v;
				vtemp.remove_variables(this->domain().get_variable_ids());
				if (!vtemp.is_zero()) {
					is_bounded = false;
					return;
				} else {
					// remap, i.e., ignore excess vars since they have coeff zero
					v.remap(this->domain());
				}
			} else {
				// we just have to put the variables in the right places
				v.reorder(this->get_index_to_variable_id_map());
			}
		}
		direction d = v.get_vector().template convert_to<scalar_type>();

		sfm_cont_set<scalar_type>* nonconst_this = const_cast<sfm_cont_set<
				scalar_type>*> (this);
		std::pair<unsigned int,scalar_type> pr = nonconst_this->extend_sfm(d);
		unsigned int row_index = pr.first;
		scalar_type f = pr.second;

		// find the largest value
		scalar_type x_max = _my_sfm(row_index, 0);
		unsigned int j_max = 0;
		for (unsigned int j = 1; j < _my_sfm.size2(); j++) {
			if (_my_sfm(row_index, j)>x_max) {
				x_max = _my_sfm(row_index, j);
				j_max = j;
			}
		}

		max_value = convert_element<double>(x_max/f);

		// if there is an invariant or nonstate constraints, we need to work
		// with the outer poly instead
		if (false) // GF 2011-12-02 deactivated because it interferes with precise intersection
		if (my_postc_prb.invariant_set_ptr || my_nonstate_contraints.size()>0) {
			// take into account the invariant and nonstate constraints, return the
			// support value of the outer poly

			// find the largest value of the sfm until it is smaller than
			// the value of the outer poly
			std::set<unsigned int> checked;
			double poly_max;
			// careful: the poly returns the actual support value, while x_max is scaled with 1/f!
			bool nonefound;
			do {
				get_polytope(j_max).compute_support(l,poly_max,support_vec,is_empty,is_bounded);
				checked.insert(j_max);
//std::cout << "current candidate:" << j_max << " at " << poly_max << " checked: " << checked.size() << " out of " << get_size() << std::endl;
				// if it's empty, restart the search
				nonefound = true;
				// find the second largest x in the sfm excluding those already checked
				for (unsigned int j = 0; j < _my_sfm.size2(); j++) {
					if ((nonefound || _my_sfm(row_index, j) > x_max))
						if (checked.find(j) == checked.end()) {
							nonefound = false;
							x_max = _my_sfm(row_index, j);
							j_max = j;
						}
				}
			} while (!nonefound && (is_empty || math::numeric::is_GT(convert_element<double>(x_max/f),poly_max)));

			max_value = poly_max;
		}
	}
}

}
}
