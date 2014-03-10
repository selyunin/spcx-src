/*
 * simplex.cpp
 *
 *  Created on: Jul 8, 2011
 *      Author: frehse
 */

#include "simplex.h"

//#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"
//#include "math/vector_utility.h"
////#include "math/vdom/affine_map_utility.h"
////#include "math/vdom/affine_map.h"
////#include "math/matrix.h"
////#include "math/vector.h"
#include "math/vdom/vdom_vector.h"
#include "math/vdom/vdom_matrix.h"

namespace continuous {

template<class scalar_type> simplex<scalar_type>::simplex()  {
}

template<class scalar_type> simplex<scalar_type>::simplex(positional_vdomain dom) :
		//	for (it = vertices.begin(); it!=vertices.end(); it++){
		//		vdom_vector<scalar_type> stateV(domain(),*it);
		//		vdom_vector<scalar_type> deriv(compute_deriv(stateV));
		//		for(int j=0;j<dim;j++){
		//			b[j+i] = deriv[j];
		//			for(int k=0;k<dim;k++){
		//				A(i+j)(k+j*dim) = it->first()[k];
		//			}
		//			A(i+j)(dim+j*dim) = 1;
		//		}
		//		i+=dim;
		//	}
														index_to_variable_id_map_provider(dom) {
}

template<class scalar_type> simplex<scalar_type>::simplex(positional_vdomain dom, scalar_type radius) :
																					index_to_variable_id_map_provider(dom) {
	// let origin be the center of the sphere
	int dim = dom.size();
	typedef typename simplex<scalar_type>::point_type point_type;
	point_type x(dom.size());

	// compute constant edge

	// compute first radius and first sphere center
	scalar_type current_radius = radius;
	point_type current_center(dom.size(),0.0);
	// initial sphere center = origin

	// recursive part
	for(int i=0; i<dim;i++) {
		// new point on
		for (int j=0; j< i ; j++) {
			//current_point[i][j]= dimAlreadySet[j];
			x[j] = current_center[j];

		}
		x[i] = current_center[i]+ current_radius ;
		for (int j=i+1; j<dim ; j++) {
			//other coordinate set to 0
			x[j]= 0;
		}

		this->add_vertice(x);

		// now we compute the hyper plan
		if(i==dim-1){
			// we compute the last vertex
			x[dim - 1] = current_center[dim - 1] - current_radius;
			for (int j=0; j<(dim-1) ; j++) {
				x[j] = current_center[j];
			}
			this->add_vertice(x);
		}
		else{
			// compute the new sphere center
			scalar_type dFromCenter = (dim-i+1)*current_radius/(dim-i);
			current_center[i] =  current_center[i]+ current_radius - dFromCenter;

			// compute the new radius
			scalar_type h = current_radius;
			scalar_type o = current_radius/(dim-i);
			current_radius = sqrt(fabs(h*h - o*o));
		}
	}
}

template<class scalar_type> continuous_set::ptr simplex<scalar_type>::get_ptr() {
	continuous_set::ptr p = boost::enable_shared_from_this<simplex<
			scalar_type> >::shared_from_this();
	return p;
}

template<class scalar_type> continuous_set::const_ptr simplex<
scalar_type>::get_const_ptr() const {
	continuous_set::const_ptr p = boost::enable_shared_from_this<
			simplex<scalar_type> >::shared_from_this();
	return p;
}

template<class scalar_type> simplex<scalar_type>* simplex<
scalar_type>::create_universe() const {
	throw std::runtime_error(
			"simplex<scalar_type>::create_universe() not allowed");
	return new simplex<scalar_type> ();
}

template<class scalar_type> simplex<scalar_type>* simplex<
scalar_type>::clone() const {
	return new simplex<scalar_type> (*this);
}

template<class scalar_type> simplex<scalar_type>* simplex<
scalar_type>::create_empty() const {
	throw std::runtime_error(
			"simplex<scalar_type>::create_empty() not allowed");
	return new simplex<scalar_type> ();
}

template<class scalar_type> int simplex<scalar_type>::get_memory() const {
	throw std::runtime_error(
			"simplex<scalar_type>::get_memory() missing implementation");
	return 0;
}

template<class scalar_type> unsigned int simplex<scalar_type>::get_dim() const {
	return get_index_to_variable_id_map()->dimensions();
}

template<class scalar_type> math::tribool simplex<scalar_type>::is_empty() const {
	return false;
}

template<class scalar_type> math::tribool simplex<scalar_type>::is_universe() const {
	return false;
}

template<class scalar_type> void simplex<scalar_type>::embed_variables(
		const variable_id_set& id_set) {
	throw std::runtime_error(
			"simplex<scalar_type>::embed_variables() not allowed");
}

template<class scalar_type> void simplex<scalar_type>::existentially_quantify_variables(
		const variable_id_set& id_set) {
	//	index_to_variable_id_map_ptr p =
	//			get_index_to_variable_id_map()->get_map_with_ids_removed(id_set);
	//
	//	// remap my_l and my_u
	//	index_to_index_bimap adapt_map = get_index_to_index_mapping(p,
	//			get_index_to_variable_id_map());
	//
	//	// it doesn't matter whether things are updated or not, since
	//	// they keep their (updated or not) value
	//	my_c = math::map(my_g, adapt_map);
	//	my_c = math::map(my_g, adapt_map);
	//	// since existentially_quantify_variables only removes variables,
	//	// we don't need to define any uninitialized elements
	//
	//	// fix _index_to_variable_id_map_ptr
	//	set_index_to_variable_id_map(p);

	throw std::runtime_error(
			"simplex<scalar_type>::existentially_quantify_variables() missing implementation");
}

template<typename scalar_type> bool simplex<scalar_type>::computes_support_vector() const {
	return true;
}
;

template<class scalar_type>
void simplex<scalar_type>::compute_support(const math::vdom_vector<
		Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	throw std::runtime_error("simplex<scalar_type>::compute_support(Rational) missing implementation");
}

template<class scalar_type>
void simplex<scalar_type>::compute_support(const math::vdom_vector<
		double>& l, double& max_value, math::vdom_vector<double>& support_vec,
		bool& is_empty, bool& is_bounded) const {
	throw std::runtime_error(
			"simplex<scalar_type>::compute_support(double) missing implementation");
}

template<class scalar_type>
typename math::lin_constraint_system<scalar_type>::const_ptr simplex<
scalar_type>::get_constraints() const {
	typename math::lin_constraint_system<scalar_type>::ptr cons(
			new math::lin_constraint_system<scalar_type>());

	math::lin_constraint<scalar_type> c;
	if (!math::definitely(is_empty())) {
		for (size_type i = 0; i < get_dim(); ++i) {
			// a linear constraint is given by a_1*x_1+...+a_n*x_n+b <= 0
			math::lin_expression<scalar_type> e(get_index_to_variable_id_map());
			e[i] = scalar_type(1); // set a[i] to 1
			// e.set_inh_coeff(...); // define b
			c = math::lin_constraint<scalar_type>(e, LE);
			cons->insert(c);

			// @todo to be completed
		}
	} else {
		c = math::lin_constraint<scalar_type>::zero_dim_false();
		cons->insert(c);
	}
	return cons;
}

template<class scalar_type>
void simplex<scalar_type>::add_constraint(const math::lin_constraint<
		scalar_type> &c, bool check_redundancy) {
	throw std::runtime_error(
			"simplex<scalar_type>::add_constraint() missing implementation");
}

template<class scalar_type>
void simplex<scalar_type>::remove_redundant_constraints() {
	// do nothing
}

template<class scalar_type>
const typename simplex<scalar_type>::point_set_type& simplex<scalar_type>::get_vertices() const {
	return my_vertices;
}

template<class scalar_type>
void simplex<scalar_type>::add_vertice(const point_type &c) {
	assert(c.size()==get_dim());

	my_vertices.insert(c,false); // false is a dummy value that is unused
}

}

